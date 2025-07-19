package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;
    static int segLen = IntVector.SPECIES_256.length(); //number of vectors

    static int segN;

    public static void main(String[] args) throws ServiceException {

        long startTime = System.currentTimeMillis();

        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415";
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int lenD = D.length();
        String querySequenceName = "P49950";
        String Q = uniProtService.getEntry(querySequenceName).getSequence().getValue();
        int lenQ = Q.length();

        segN = (int) Math.ceil((double) (lenQ + segLen - 1) / segLen);

        int paddedLen = segN * segLen;
        for ( int i = 0; i < paddedLen - lenQ; i++ ) {
            Q += " ";
        }
        lenQ = paddedLen;

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

        float[][] tempProfile = new float[lenD][lenQ];
        char[] charArrayQ = Q.toCharArray();
        float[] floatArrayQ = new float[charArrayQ.length];
        for ( int i = 0; i < charArrayQ.length; i++ ) {
            floatArrayQ[i] = (float) charArrayQ[i];
        }

        //initializing match/mismatch profile
        for ( int i = 0; i < lenD; i++ ) {
            int residueD = D.charAt(i);
            FloatVector vResidueD = FloatVector.broadcast(SPECIES, residueD);

            for ( int j = 0; j < lenQ; j+=segLen ) { //
                VectorMask<Float> mask = SPECIES.indexInRange(j, j+segLen);
                FloatVector vResidueQ = FloatVector.fromArray(SPECIES, floatArrayQ, j, mask);
                VectorMask<Float> residueComparisonMask = vResidueD.compare(VectorOperators.EQ, vResidueQ);
                FloatVector profileVector = vUnmatch.blend(vMatch, residueComparisonMask);
                float[] profileVectorArr = profileVector.toArray();
                for( int k = 0; k < segLen; k++ ) {
                    tempProfile[i][k+j] = profileVectorArr[k];
                }
            }

        }

        float[][] profile = new float[lenD+1][lenQ+1];
        for ( int i = 0; i < lenD; i++ ) {
            for ( int j = 0; j < lenQ; j++ ) {
                profile[i+1][j+1] = tempProfile[i][j];
            }

        }

        /*
        E --> for insertion; from left.
        F --> for deletion; from above
        H --> for match/mismatch; upper-left diagonal
        The query is divided into equal length segments, S.
         */

        //System.out.println("No error so far");

        //Row: target sequence. Column: query sequence. H[row][column]
        //TODO: Outer loop and Inner loop, Lazy-F loop
        float[] allZeroes = new float[segN];

        FloatVector[] vHStore = new FloatVector[segN];
        FloatVector[] vHLoad = new FloatVector[segN];
        FloatVector[] vE = new FloatVector[segN];
        for ( int i = 0; i < segLen; i++ ) {
            vHStore[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
            vHLoad[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
            vE[i] = FloatVector.fromArray(SPECIES, allZeroes,0);
        }

        FloatVector vMax = FloatVector.fromArray(SPECIES, allZeroes, 0);

        for ( int i = 1; i <= lenD; i++ ) {

            FloatVector vF = FloatVector.fromArray(SPECIES, allZeroes,0);
            FloatVector[] vHStoreLeftShift = new FloatVector[0];
            if ( i != 0 ) {
                vHStoreLeftShift = leftShift(vHStore);
            }
            FloatVector vH = vHStoreLeftShift[0];

            FloatVector[] vHTemp = vHLoad;
            vHLoad = vHStore;
            vHStore = vHTemp;


            for (int j = 0; j < segN; j++) {
                vH = vH.add(profile[i][j]);                    // Add score from scoring profile
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

            // Lazy F loop to propagate corrections vertically
            vF = leftShift(vF); // vF <<= 1
            int j = 0;
            while (j < segN) {
                FloatVector vHj = vHStore[j];
                FloatVector vHnew = vHj.max(vF);
                vHStore[j] = vHnew;
                vF = vF.sub(vGapE);

                j++;
                if (j == segN && vF.compare(VectorOperators.GT, vHj.sub(vGapO)).anyTrue()) {
                    vF = leftShift(vF);  // vF <<= 1
                    j = 0;
                }
            }

            long endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;
            System.out.println("Elapsed time: " + elapsedTime + " milliseconds");
        }
    }

    public static FloatVector[] leftShift(FloatVector[] v1) {
        //turn v1 to one-dimensional array
        int laneCount = SPECIES.length();
        int totalLen = v1.length * laneCount;

        float[] flat = new float[totalLen];
        for (int i = 0; i < v1.length; i++) {
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
        //turn v1 to one-dimensional array
        float[] flat = new float[v1.length()];
        v1.intoArray(flat, 0);

        float first = flat[0];
        System.arraycopy(flat, 1, flat, 0, v1.length() - 1);
        flat[v1.length() - 1] = first;

        return FloatVector.fromArray(SPECIES, flat, 0);
    }
}