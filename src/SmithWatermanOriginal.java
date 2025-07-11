package uk.ac.ebi.uniprot.dataservice.client.examples;

import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

//CREDIT: https://github.com/JayakrishnaThota/Sequence-Alignment/blob/master/SmithWaterman.java


import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

public class SmithWatermanOriginal {
    public static void main(String[] args) throws ServiceException {
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415";
        String targetSequence = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int targetSeqLen = targetSequence.length();
        String querySequenceName = "P49950";
        String querySequence = uniProtService.getEntry(querySequenceName).getSequence().getValue();
        int querySeqLen = querySequence.length();

        int[][] H = new int[targetSeqLen+1][querySeqLen+1];
        int[][] E = new int[targetSeqLen+1][querySeqLen+1];
        int[][] F = new int[targetSeqLen+1][querySeqLen+1];

        int matchScore = 20;
        int mismatchScore = 10;
        int gapOpeningPenalty = 10;
        int gapExtensionPenalty = 3;

        int Emax;
        int Fmax;
        int Hvalue = 0;
        for ( int i = 1; i <= targetSeqLen; i++ ) {
            for ( int j = 1; j <= querySeqLen; j++ ) {
                if ( targetSequence.charAt(i-1) == querySequence.charAt(j-1) ) {
                    Hvalue = H[i-1][j-1] + matchScore;
                }
                else {
                    Hvalue = H[i-1][j-1] - mismatchScore;
                }

                Fmax = Math.max(E[i][j-1] - gapExtensionPenalty, H[i][j-1] - gapOpeningPenalty);
                F[i][j] = Fmax;

                Emax = Math.max(E[i-1][j] - gapExtensionPenalty, H[i-1][j] - gapOpeningPenalty);
                E[i][j] = Emax;

                H[i][j] = Math.max(Hvalue, Math.max(Math.max(Fmax, Emax),0));
            }
        }
        System.out.println();
    }

}