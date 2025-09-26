package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;
    static int segLen = FloatVector.SPECIES_256.length(); //number of vectors

    static int segN;

    public static long StripedLayout( String D, int lenD, String Q, int lenQ ) {
        long startTime = System.nanoTime();

        segN = (lenQ + segLen - 1) / segLen;

        int paddedLen = segN * segLen;
        float[] floatArrayQ = new float[paddedLen];
        int originalLenQ = Q.length();


        for (int i = 0, t = 0, j = 0; i < originalLenQ; i++) {
            floatArrayQ[t * segLen + j] = Q.charAt(i);
            t++;
            if (t >= segN) {
                j++;
                t = 0;
            }
        }

        float matchScore = 20.0f;
        float unmatchScore = 10.0f;
        float gapOpeningPenalty = 10.0f;
        float gapExtensionPenalty = 3.0f;

        float[] vMatch_ = new float[segLen];
        Arrays.fill(vMatch_, matchScore);
        float[] vUnmatch_ = new float[segLen];
        Arrays.fill(vUnmatch_, -unmatchScore);
        float[] vGapOpen_ = new float[segLen];
        Arrays.fill(vGapOpen_, -gapOpeningPenalty);
        float[] vGapExtend_ = new float[segLen];
        Arrays.fill(vGapExtend_, -gapExtensionPenalty);

        FloatVector vGapO = FloatVector.fromArray(SPECIES, vGapOpen_, 0);
        FloatVector vGapE = FloatVector.fromArray(SPECIES, vGapExtend_, 0);
        FloatVector vMatch = FloatVector.fromArray(SPECIES, vMatch_, 0);
        FloatVector vUnmatch = FloatVector.fromArray(SPECIES, vUnmatch_, 0);

        float[] allZeroes = new float[SPECIES.length()];

        FloatVector[][] vProfile = new FloatVector[lenD][segN];
        for ( int i = 0; i < lenD; i++ ) {
            int residueD = D.charAt(i);
            FloatVector vResidueD = FloatVector.broadcast(SPECIES, residueD);
            for ( int j = 0; j < segN; j++ ) {
                int qSegmentStartIndex = j * segLen;
                FloatVector vResidueQ = FloatVector.fromArray(SPECIES, floatArrayQ, qSegmentStartIndex);
                VectorMask<Float> residueComparisonMask = vResidueD.compare(VectorOperators.EQ, vResidueQ);
                VectorMask<Float> withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN);
                vProfile[i][j] = FloatVector.fromArray(SPECIES, allZeroes,0).blend((vUnmatch.blend(vMatch, residueComparisonMask)), withinRange);
            }
        }

        FloatVector[][] vHResult = new FloatVector[lenD][segN];
        FloatVector[] vHStore = new FloatVector[segN];
        FloatVector[] vHLoad = new FloatVector[segN];
        FloatVector[] vE = new FloatVector[segN];
        FloatVector vMax;
        FloatVector vF, vH, vZero;
        VectorMask<Float> withinRange;

        for ( int i = 0; i < segN; i++ ) {
            vHStore[i] = FloatVector.fromArray(SPECIES, allZeroes, 0);
            vHLoad[i] = FloatVector.fromArray(SPECIES, allZeroes, 0);
            vE[i] = FloatVector.fromArray(SPECIES, vGapExtend_,0);
        }
        vZero = FloatVector.fromArray(SPECIES, allZeroes, 0);
        vMax = FloatVector.fromArray(SPECIES, allZeroes, 0);

        // Outer loop to process the database sequence
        // for i:=0...dbLen
        for ( int i = 0; i < lenD; i++ ) {
            // Initialize F value to zeros. Any errors to vH values will be corrected in the Lazy-F loop.
            vF = FloatVector.fromArray(SPECIES, allZeroes,0);

            // Adjust the last H value to be used in the next segment over
            withinRange = SPECIES.indexInRange(1, (lenQ-0 + segN-1)/segN );
            vH = vHStore[segN - 1].unslice(1, vZero, 0, withinRange);

            // Swap the two H buffers - change to store all
            vHLoad = vHStore;
            vHStore = vHResult[i];

            // Inner loop to process the query sequence
            for (int j = 0; j < segN; j++) {
                // Add the scoring profile to vH
                withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN );
                vH = vH.add(vProfile[i][j], withinRange); // Add score from scoring profile

                // Save the vH values greater than the max
                vMax = vMax.max(vH);

                // Adjust vH with any greater vE or vH values
                vH = vH.max(vE[j]);
                vH = vH.max(vF);

                // Save the vH values off
                vHStore[j] = vH.max(vZero); // compare 0 and copy

                // Calculate the new vE and vF based on the gap penalties for this search
                vH = vH.add(vGapO, withinRange);
                vE[j] = vE[j].add(vGapE, withinRange);
                vE[j] = vE[j].max(vH);
                vF = vF.add(vGapE, withinRange);
                vF = vF.max(vH);

                // Load next vH
                vH = vHLoad[j];
            }

            // --- Lazy-F loop ---
            // Shift the vF left so its values can be used to correct the next segment over
            // probably wrong in paper vF := VF >> 1;
            withinRange = SPECIES.indexInRange(1, (lenQ-0 + segN-1)/segN );
            vF = vF.unslice(1, vZero, 0, withinRange);

            // Correct the vH values until there are no elements in vF that could influence the vH values
            int j = 0;
            withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN );
            VectorMask<Float> vCompare = vF.compare(VectorOperators.GT, vHStore[j], withinRange);
            while(vCompare.anyTrue()) {
                vHStore[j] = vHStore[j].max(vF);
                vMax = vMax.max(vHStore[j]); // probably missed in paper
                vH = vHStore[j].add(vGapO, withinRange);
                vE[j] = vE[j].max(vH); // probably missed in paper
                vF = vF.add(vGapE, withinRange);
                vF = vF.max(vH); // probably missed in paper
                if (++j >= segN) {
                    // If we processed the entire segment, we need to carry the vF values to the next segment
                    withinRange = SPECIES.indexInRange(1, (lenQ-0 + segN-1)/segN );
                    vF = vF.unslice(1, vZero, 0, withinRange);
                    j = 0;
                }

                withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN );
                vCompare = vF.compare(VectorOperators.GT, vHStore[j], withinRange);
            }
        }

        long endTime = System.nanoTime();
        System.out.println("Stripped-Layout max " + vMax.reduceLanes(VectorOperators.MAX));
        //System.out.println("FloatVector[] Left Shift time: " + maxTime);
        return (long) endTime - startTime;
    }
}