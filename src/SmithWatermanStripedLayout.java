package uk.ac.ebi.uniprot.dataservice.client.examples;

import jdk.incubator.vector.*;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;

    public static void main(String[] args) throws ServiceException {

        //Extracting sequences + initializing target and query sequences and their lengths
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415";
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int lenD = D.length();
        String querySequenceName = "P49950";
        String Q = uniProtService.getEntry(querySequenceName).getSequence().getValue();
        int lenQ = Q.length();

        int segLen = IntVector.SPECIES_256.length(); //number of vectors

        int segN = (int) Math.ceil((lenQ+segLen - 1) / segLen); //length of each vector

        int paddedLen = segN * segLen;
        for ( int i = 0; i <= paddedLen - lenQ; i++ ) {
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

        FloatVector vGapOpen = FloatVector.fromArray(SPECIES, vGapOpen_, 0);
        FloatVector vGapExtend = FloatVector.fromArray(SPECIES, vGapExtend_, 0);
        FloatVector vMatch = FloatVector.fromArray(SPECIES, vMatch_, 0);
        FloatVector vUnmatch = FloatVector.fromArray(SPECIES, vUnmatch_, 0);

        VectorSpecies<Float> target = FloatVector.SPECIES_PREFERRED;
        VectorSpecies<Float> query = FloatVector.SPECIES_PREFERRED;


        float[][] profile = new float[lenD][lenQ];
        char[] charArrayQ = Q.toCharArray();
        float[] floatArrayQ = new float[charArrayQ.length];
        for ( int i = 0; i < charArrayQ.length; i++ ) {
            floatArrayQ[i] = (float) charArrayQ[i];
        }

        //initializing match/mismatch profile
        for ( int i = 0; i < lenD; i++ ) {
            int residueD = D.charAt(i);
            FloatVector vResidueD = FloatVector.broadcast(target, residueD);

            for ( int j = 0; j < lenQ; j+=segLen ) { //
                VectorMask<Float> mask = query.indexInRange(j, j+segLen);
                FloatVector vResidueQ = FloatVector.fromArray(query, floatArrayQ, j, mask);
                VectorMask<Float> residueComparisonMask = vResidueD.compare(VectorOperators.EQ, vResidueQ);
                FloatVector profileVector = vUnmatch.blend(vMatch, residueComparisonMask);
                float[] profileVectorArr = profileVector.toArray();
                for( int k = 0; k < segLen; k++ ) {
                    profile[i][k+j] = profileVectorArr[k];
                }
            }

        }

        /*
        E --> for insertion; from left.
        F --> for deletion; from above
        H --> for match/mismatch; upper-left diagonal
        The query is divded into equal length segments, S.
         */

        System.out.println("No error");
        
    }
}