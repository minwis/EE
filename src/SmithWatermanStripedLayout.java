package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;
    static int segLen = FloatVector.SPECIES_256.length(); //number of vectors

    static int segN;

    public static long StripedLayout( String D, int lenD, String Q, int lenQ ) {
        long startTime = System.currentTimeMillis();

        segN = (int) Math.ceil((double) (lenQ + segLen - 1) / segLen);

        int paddedLen = segN * segLen;
        float blankResidue = ' ';
        float[] floatArrayQ = new float[paddedLen];
        int originalLenQ = Q.length();
        for (int i = 0; i < originalLenQ; i++) {
            floatArrayQ[i] = Q.charAt(i);
        }
        for (int i = originalLenQ; i < paddedLen; i++) {
            floatArrayQ[i] = blankResidue;
        }

        float matchScore = 20.0f;
        float unmatchScore = -10.0f;
        float gapOpeningPenalty = -10.0f;
        float gapExtensionPenalty = -3.0f;

        float[] vMatch_ = new float[segLen];
        Arrays.fill(vMatch_, matchScore);
        float[] vUnmatch_ = new float[segLen];
        Arrays.fill(vUnmatch_, unmatchScore);
        float[] vGapOpen_ = new float[segLen];
        Arrays.fill(vGapOpen_, gapOpeningPenalty);
        float[] vGapExtend_ = new float[segLen];
        Arrays.fill(vGapExtend_, gapExtensionPenalty);

        FloatVector vGapO = FloatVector.fromArray(SPECIES, vGapOpen_, 0);
        FloatVector vGapE = FloatVector.fromArray(SPECIES, vGapExtend_, 0);
        FloatVector vMatch = FloatVector.fromArray(SPECIES, vMatch_, 0);
        FloatVector vUnmatch = FloatVector.fromArray(SPECIES, vUnmatch_, 0);

        FloatVector[][] vProfile = new FloatVector[lenD][segN];
        for ( int i = 0; i < lenD; i++ ) {
            int residueD = D.charAt(i);
            FloatVector vResidueD = FloatVector.broadcast(SPECIES, residueD);
            for ( int j = 0; j < segN; j++ ) {
                int qSegmentStartIndex = j * segLen;
                FloatVector vResidueQ = FloatVector.fromArray(SPECIES, floatArrayQ, qSegmentStartIndex);
                VectorMask<Float> residueComparisonMask = vResidueD.compare(VectorOperators.EQ, vResidueQ);
                vProfile[i][j] = vUnmatch.blend(vMatch, residueComparisonMask);
            }

        }

        float[] allZeroes = new float[SPECIES.length()];

        FloatVector[] vHStore = new FloatVector[segN];
        FloatVector[] vHLoad = new FloatVector[segN];
        FloatVector[] vE = new FloatVector[segN];
        FloatVector vMax;
                FloatVector vF, vH, vESubvGapE, vFSubvGapE;
        FloatVector[] vHTemp;

        for ( int i = 0; i < segN; i++ ) {
            vHStore[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
            vHLoad[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
            vE[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
        }

        vMax = FloatVector.fromArray(SPECIES, allZeroes, 0);
        for ( int i = 0; i < lenD; i++ ) {

            vF = FloatVector.fromArray(SPECIES, allZeroes,0);
            FloatVector[] vHStoreLeftShift = leftShift(vHStore);

            vH = vHStoreLeftShift[0];

            vHTemp = vHLoad;
            vHLoad = vHStore;
            vHStore = vHTemp;

            for (int j = 0; j < segN; j++) {

                vH = vH.add(vProfile[i][j]); // Add score from scoring profile
                vMax = vH.max(vMax); // Track global max

                vH = vH.max(vE[j]); // Compare with E
                vH = vH.max(vF); // Compare with F
                vHStore[j] = vH; // Store updated H

                // Calculate vE and vF for next iteration
                vESubvGapE = vE[j].sub(vGapE);
                vE[j] = vESubvGapE.max((vH.sub(vGapO).max(0)));
                vFSubvGapE = vF.sub(vGapE);
                vF = vFSubvGapE.max((vH.sub(vGapO).max(0)));

                // Load next vH
                vH = vHLoad[j];
            }

            vF = leftShift(vF);
            int j = 0;

            while ( LazyFLoopCondition(vF, vHStore[j].sub(vF)) ) {
                vHStore[j] = vHStore[j].max(vF);
                vF = vF.sub(vGapE);
                if ( ++j >= segN ) {
                    vF = leftShift(vF);
                    j = 0;
                }
                j++;
            }
        }

        long endTime = System.currentTimeMillis();
        System.out.println(vMax.lane(0));
        //System.out.println("FloatVector[] Left Shift time: " + maxTime);
        return (long) endTime - startTime;
    }

    public static long maxTime = 0;

    public static FloatVector[] leftShift(FloatVector[] v1) {
        long start = System.currentTimeMillis();

        int laneCount = SPECIES.length();
        int vecLen = v1.length;
        FloatVector[] result = new FloatVector[vecLen];

        int[] shuffleIndexes = new int[laneCount];

        float[] carryArr = new float[laneCount];
        VectorMask<Float> lastLaneMask = VectorMask.fromLong(SPECIES, 1L << (laneCount - 1));

        for (int i = 0; i < laneCount - 1; i++) {
            shuffleIndexes[i] = i + 1;
        }
        shuffleIndexes[laneCount - 1] = 0;
        VectorShuffle<Float> laneShuffle = VectorShuffle.fromArray(SPECIES, shuffleIndexes, 0);

        // Carry-over value for rotation between vectors
        float carry = v1[0].lane(0); // first element becomes the last

        for (int i = 0; i < vecLen; i++) {
            FloatVector vec = v1[i];
            FloatVector shuffled = vec.rearrange(laneShuffle);

            // Bring carry from previous vector into last lane
            carryArr[laneCount - 1] = carry;
            FloatVector carryVec = FloatVector.fromArray(SPECIES, carryArr, 0);
            shuffled = shuffled.blend(carryVec, lastLaneMask);

            // Update carry for next vector (it’s the old lane 0)
            carry = vec.lane(0);

            result[i] = shuffled;
        }

        // Put original first value into last vector’s last lane (full rotate)
        result[vecLen - 1] = result[vecLen - 1].withLane(laneCount - 1, carry);

        long end = System.currentTimeMillis();
        if (end - start > maxTime) {
            maxTime = end - start;
        }


        return result;

    }

    public static FloatVector leftShift(FloatVector v1) {
        int[] orderArr = new int[segLen];
        orderArr[0] = segLen - 1;
        for ( int i = 1; i < segLen; i++ ) {
            orderArr[i] = i-1;
        }
        VectorShuffle<Float> shuffle = VectorShuffle.fromArray(SPECIES, orderArr, 0);
        return v1.rearrange(shuffle);
    }

    public static boolean LazyFLoopCondition(FloatVector v1, FloatVector v2) {
        VectorMask<Float> compare = v1.compare(VectorOperators.GT, v2);
        return !compare.allTrue();
    }
}