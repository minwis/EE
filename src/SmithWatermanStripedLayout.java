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

        int segN = IntVector.SPECIES_256.length();
        System.out.println(segN);

        int segLen = (lenQ+segN - 1) / segN;

        int paddedLen = segN * segLen;
        if ( lenQ > paddedLen ) {
            for ( int i = 0; i < lenQ-paddedLen; i++ ) {
                Q += " ";
            }
        }

        float matchScore = 20.0f;
        float unmatchScore = -10.0f;
        float gapOpeningPenalty = -10.0f;
        float gapExtensionPenalty = -3.0f;

        float[] vMatch_ = new float[segN];
        Arrays.fill(vMatch_, matchScore);
        float[] vUnmatch_ = new float[segN];
        Arrays.fill(vUnmatch_, unmatchScore);
        float[] vGapOpen_ = new float[segN];
        Arrays.fill(vUnmatch_, gapOpeningPenalty);
        float[] vGapExtend_ = new float[segN];
        Arrays.fill(vGapExtend_, gapExtensionPenalty);

        FloatVector vGapOpen = FloatVector.fromArray(SPECIES, vGapOpen_, 0);
        FloatVector vGapExtend = FloatVector.fromArray(SPECIES, vGapExtend_, 0);
        FloatVector vMatch = FloatVector.fromArray(SPECIES, vMatch_, 0);
        FloatVector vUnmatch = FloatVector.fromArray(SPECIES, vUnmatch_, 0);

        VectorSpecies<Float> target = FloatVector.SPECIES_PREFERRED;
        VectorSpecies<Float> query = FloatVector.SPECIES_PREFERRED;


        float[][] profile = new float[lenD+1][lenQ+1];
        char[] charArrayQ = Q.toCharArray();
        float[] floatArrayQ = new float[charArrayQ.length];
        for ( int i = 0; i < charArrayQ.length; i++ ) {
            floatArrayQ[i] = (float) charArrayQ[i];
        }

        for ( int i = 0; i < lenD; i++ ) {
            int residueD = D.charAt(i);
            FloatVector vResidueD = FloatVector.broadcast(target, residueD);
            for ( int j = 0; j < lenQ; j+=segLen ) {
                FloatVector vResidueQ = FloatVector.fromArray(query, floatArrayQ, j);
                VectorMask<Float> residueComparisonMask = vResidueD.compare(VectorOperators.EQ, vResidueQ);
                VectorMask<Float> paddingMask = vResidueQ.compare(VectorOperators.EQ, 0);
                FloatVector profileVector = vMatch.blend(vUnmatch, residueComparisonMask);

            }

        }
        //FloatVector[][] vProfile = FloatVector.fromArray(SPECIES, profile,0);


        /*
        E --> for insertion; from left.
        F --> for deletion; from above
        H --> for match/mismatch; upper-left diagonal
         */

        for ( int i = 0; i < lenD; i++ ) {
            IntVector[] vF = new IntVector[segN];


        }
    }
}