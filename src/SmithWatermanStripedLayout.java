package uk.ac.ebi.uniprot.dataservice.client.examples;

import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

import java.util.Arrays;

public class SmithWatermanStripedLayout {

    public static void main(String[] args) throws ServiceException {

        //Extracting sequences + initializing target and query sequences and their lengths
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415";
        String targetSequence = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int targetSeqLen = targetSequence.length();
        String querySequenceName = "P49950";
        StringBuilder querySequence = new StringBuilder(uniProtService.getEntry(querySequenceName).getSequence().getValue());
        int querySeqLen = querySequence.length();

        int lanePerRegister = 4; //number of simd registers. r

        int stripLength= (querySeqLen+lanePerRegister - 1) / lanePerRegister;
        int paddedLength = lanePerRegister * stripLength;
        if ( querySeqLen > paddedLength ) {
            querySequence.append(" ".repeat(querySeqLen - paddedLength));
        }

        int matchScore = 20;
        int mismatchScore = 10;
        int gapOpeningPenalty = 10;
        int gapExtensionPenalty = 3;

        //initialize M profile
        int[][] substitutionScore = new int[targetSeqLen][querySeqLen];
        for ( int i = 0; i < targetSeqLen; i++ ) {
            for ( int j = 0; j < querySeqLen; j++ ) {
                if ( targetSequence.charAt(i) == querySequence.charAt(j) ) {
                    substitutionScore[i][j] = matchScore;
                }
                else {
                    substitutionScore[i][j] = mismatchScore;
                }
            }
        }

        /*
        E --> for insertion; from left.
        F --> for deletion; from above
        H --> for match/mismatch; upper-left diagonal
         */

        int[] vF = new int[stripLength];
        int[] vH = new int[stripLength];
        int[][] H = new int[targetSeqLen][querySeqLen];
        int Hi = 0;
        for ( int j = 1; j <= targetSeqLen; j++ ) {
            vF = new int[stripLength];

            if ( Hi > querySeqLen ) {
                Hi = -1;
            }
            else {
                Hi += lanePerRegister;
            }

            vH = Arrays.copyOfRange(H[j], Hi, Hi+lanePerRegister);


        }
    }

    /*public static int[] leftShift(int[][] matrix, int[] vector, int i) {

    }*/

}