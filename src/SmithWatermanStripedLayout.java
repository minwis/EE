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



        //TODO: MAIN LOOP
        for ( int i = 0; i < lenD; i++ ) {

            vF = FloatVector.fromArray(SPECIES, allZeroes,0);
            FloatVector[] vHStoreLeftShift = leftShift(vHStore);

            vH = vHStoreLeftShift[0];

            vHTemp = vHLoad;
            vHLoad = vHStore;
            vHStore = vHTemp;

            for (int j = 0; j < segN; j++) {

                VectorMask<Float> withinRange = SPECIES.indexInRange(j*segLen, lenQ);
                vH = vH.add(vProfile[i][j], withinRange); // Add score from scoring profile


                vMax = vH.max(vMax);

                vH = vH.max(vE[j]);
                vH = vH.max(vF);
                vHStore[j] = vH;


                vESubvGapE = vE[j].sub(vGapE, withinRange);
                vE[j] = vESubvGapE.max((vH.sub(vGapO, withinRange).max(0)));

                vFSubvGapE = vF.sub(vGapE, withinRange);
                vF = vFSubvGapE.max((vH.sub(vGapO).max(0)));

                // Load next vH
                vH = vHLoad[j];
            }

            vF = leftShift(vF);
            int j = 0;

            while ( ++j < segN && LazyFLoopCondition(vF, vHStore[j].sub(vGapO)) ) {
                VectorMask<Float> withinRange = SPECIES.indexInRange(j*segLen, lenQ);
                vHStore[j] = vHStore[j].max(vF);
                vF = vF.sub(vGapE, withinRange);
            }

            vMax = vH.max(vMax);
        }

        long endTime = System.nanoTime();
        System.out.println(vMax.reduceLanes(VectorOperators.MAX));
        //System.out.println("FloatVector[] Left Shift time: " + maxTime);
        return (long) endTime - startTime;
    }

    public static long maxTime = 0;

    public static FloatVector[] leftShift(FloatVector[] v1) {

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