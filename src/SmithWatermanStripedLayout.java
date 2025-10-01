package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;
    static int segLen = SPECIES.length();

    static int segN;

    public static long StripedLayout( String T, int lenT, String Q, int lenQ ) {
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

        FloatVector[][] vProfile = new FloatVector[lenT][segN];
        for ( int i = 0; i < lenT; i++ ) {
            int residueT = T.charAt(i);
            FloatVector vResidueT = FloatVector.broadcast(SPECIES, residueT);
            for ( int j = 0; j < segN; j++ ) {
                int qSegmentStartIndex = j * segLen;
                FloatVector vResidueQ = FloatVector.fromArray(SPECIES, floatArrayQ, qSegmentStartIndex);
                VectorMask<Float> residueComparisonMask = vResidueT.compare(VectorOperators.EQ, vResidueQ);
                VectorMask<Float> withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN);
                vProfile[i][j] = FloatVector.fromArray(SPECIES, allZeroes,0).blend((vUnmatch.blend(vMatch, residueComparisonMask)), withinRange);
            }
        }

        FloatVector[][] vHResult = new FloatVector[lenT][segN];
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

        for ( int i = 0; i < lenT; i++ ) {
            vF = FloatVector.fromArray(SPECIES, allZeroes,0);

            withinRange = SPECIES.indexInRange(1, (lenQ-0 + segN-1)/segN );
            vH = vHStore[segN - 1].unslice(1, vZero, 0, withinRange);

            vHLoad = vHStore;
            vHStore = vHResult[i];

            for (int j = 0; j < segN; j++) {
                withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN );
                vH = vH.add(vProfile[i][j], withinRange);

                vMax = vMax.max(vH);

                vH = vH.max(vE[j]);
                vH = vH.max(vF);

                vHStore[j] = vH.max(vZero);

                vH = vH.add(vGapO, withinRange);
                vE[j] = vE[j].add(vGapE, withinRange);
                vE[j] = vE[j].max(vH);
                vF = vF.add(vGapE, withinRange);
                vF = vF.max(vH);

                vH = vHLoad[j];
            }

            withinRange = SPECIES.indexInRange(1, (lenQ-0 + segN-1)/segN );
            vF = vF.unslice(1, vZero, 0, withinRange);

            int j = 0;
            withinRange = SPECIES.indexInRange(0, (lenQ-j + segN-1)/segN );
            VectorMask<Float> vCompare = vF.compare(VectorOperators.GT, vHStore[j], withinRange);
            while(vCompare.anyTrue()) {
                vHStore[j] = vHStore[j].max(vF);
                vMax = vMax.max(vHStore[j]);
                vH = vHStore[j].add(vGapO, withinRange);
                vE[j] = vE[j].max(vH);
                vF = vF.add(vGapE, withinRange);
                vF = vF.max(vH);
                if (++j >= segN) {
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
        return (long) endTime - startTime;
    }
}