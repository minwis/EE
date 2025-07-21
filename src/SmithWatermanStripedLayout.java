package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;
    static int segLen = IntVector.SPECIES_256.length(); //number of vectors

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

        //lenQ = paddedLen;

        float matchScore = 20.0f;
        float unmatchScore = -10.0f;
        float gapOpeningPenalty = -10.0f;
        float gapExtensionPenalty = -3.0f;

        float[] vMatch_ = new float[segLen];
        Arrays.fill(vMatch_, matchScore);
        float[] vUnmatch_ = new float[segLen];
        Arrays.fill(vUnmatch_, unmatchScore);
        float[] vGapOpen_ = new float[segLen];
        Arrays.fill(vUnmatch_, gapOpeningPenalty);
        float[] vGapExtend_ = new float[segLen];
        Arrays.fill(vGapExtend_, gapExtensionPenalty);

        FloatVector vGapO = FloatVector.fromArray(SPECIES, vGapOpen_, 0);
        FloatVector vGapE = FloatVector.fromArray(SPECIES, vGapExtend_, 0);
        FloatVector vMatch = FloatVector.fromArray(SPECIES, vMatch_, 0);
        FloatVector vUnmatch = FloatVector.fromArray(SPECIES, vUnmatch_, 0);

        FloatVector[][] vProfile = new FloatVector[lenD][segN];

        long debugTimeStamp1 = System.currentTimeMillis();

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

        long debugTimeStamp2 = System.currentTimeMillis();
        System.out.println("Time taken until before major loop: " + (debugTimeStamp2-debugTimeStamp1));


        float[] allZeroes = new float[segN];

        FloatVector[] vHStore = new FloatVector[segN];
        FloatVector[] vHLoad = new FloatVector[segN];
        FloatVector[] vE = new FloatVector[segN];
        for ( int i = 0; i < segN; i++ ) {
            vHStore[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
            vHLoad[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
            vE[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
        }

        FloatVector vMax = FloatVector.fromArray(SPECIES, allZeroes, 0);

        long loopStartTime = System.currentTimeMillis();
        for ( int i = 0; i < lenD; i++ ) {

            FloatVector vF = FloatVector.fromArray(SPECIES, allZeroes,0);
            FloatVector[] vHStoreLeftShift = leftShift(vHStore);

            FloatVector vH = vHStoreLeftShift[0];

            FloatVector[] vHTemp = vHLoad;
            vHLoad = vHStore;
            vHStore = vHTemp;

            for (int j = 0; j < segN; j++) {

                vH = vH.add(vProfile[i][j]);               // Add score from scoring profile
                vMax = vH.max(vMax);                      // Track global max

                vH = vH.max(vE[j]);                        // Compare with E
                vH = vH.max(vF);                           // Compare with F
                vHStore[j] = vH;                           // Store updated H

                // Calculate vE and vF for next iteration
                vE[j] = vE[j].sub(vGapE).max(vH.sub(vGapO));
                vF = vF.sub(vGapE).max(vH.sub(vGapO));

                // Load next vH
                vH = vHLoad[j];
            }

            vF = leftShift(vF); // vF <<= 1
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
        System.out.println("Loop time: " + (endTime-loopStartTime));
        return (long) endTime - startTime;
    }


    public static FloatVector[] leftShift(FloatVector[] v1) {
        long start = System.currentTimeMillis();
        int laneCount = SPECIES.length();
        int totalLen = v1.length * laneCount;

        float[] flat = new float[totalLen];
        for (int i = 0; i * laneCount < v1.length; i++) {
            v1[i].intoArray(flat, i * laneCount);
        }

        float first = flat[0];
        System.arraycopy(flat, 1, flat, 0, totalLen - 1);
        flat[totalLen - 1] = first;

        FloatVector[] result = new FloatVector[v1.length];
        for (int i = 0; i < v1.length; i++) {
            result[i] = FloatVector.fromArray(SPECIES, flat, i * laneCount);
        }
        long end = System.currentTimeMillis();
        System.out.println("FloatVector[] Left Shift time: " + (end-start));
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