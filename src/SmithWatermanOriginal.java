package uk.ac.ebi.uniprot.dataservice.client.examples;

//CREDIT: https://github.com/JayakrishnaThota/Sequence-Alignment/blob/master/SmithWaterman.java


public class SmithWatermanOriginal {

    public static long Original(String targetSequence, int targetSeqLen, String querySequence, int querySeqLen) {
        long startTime = System.currentTimeMillis();

        int[][] H = new int[targetSeqLen+1][querySeqLen+1];
        int[][] E = new int[targetSeqLen+1][querySeqLen+1];
        int[][] F = new int[targetSeqLen+1][querySeqLen+1];

        int matchScore = 20;
        int mismatchScore = 10;
        int gapOpeningPenalty = 10;
        int gapExtensionPenalty = 3;

        int Emax = 0;
        int Fmax = 0;
        int Hvalue = 0;
        int[][] path = new int[targetSeqLen + 1][querySeqLen + 1];
        int maxIndexI = 0;
        int maxIndexJ = 0;
        int recordMax = -1;
        for ( int i = 1; i <= targetSeqLen; i++ ) {
            for ( int j = 1; j <= querySeqLen; j++ ) {
                if ( targetSequence.charAt(i-1) == querySequence.charAt(j-1) ) {
                    Hvalue = H[i-1][j-1] + matchScore;
                }
                else {
                    Hvalue = H[i-1][j-1] - mismatchScore;
                }

                Fmax = Math.max(F[i][j-1] - gapExtensionPenalty, H[i][j-1] - gapOpeningPenalty);
                F[i][j] = Fmax;

                Emax = Math.max(E[i-1][j] - gapExtensionPenalty, H[i-1][j] - gapOpeningPenalty);
                E[i][j] = Emax;

                int currentMax = Math.max(Hvalue, Math.max(Math.max(Fmax, Emax),0));
                if ( currentMax > 0 ) {
                    H[i][j] = currentMax;
                }

                if ( currentMax > recordMax ) {
                    //maxIndexI = i;
                    //maxIndexJ = j;
                    recordMax = currentMax;
                }

                if ( H[i][j] == Hvalue ) {
                    path[i][j] = 1;
                }
                else if ( H[i][j] == Fmax ) {
                    path[i][j] = 2;
                }
                else if (H[i][j] == Emax) {
                    path[i][j] = 3;
                }

            }
        }
        /*
        int i = maxIndexI;
        int j = maxIndexJ;
        String targetMaxRegionReversed = "";
        String queryMaxRegionReversed = "";

        while ( path[i][j] != 0 ) {
            if ( path[i][j] == 1 ) {
                targetMaxRegionReversed += targetSequence.charAt(i-1);
                queryMaxRegionReversed += querySequence.charAt(j-1);
                i--;
                j--;
            }
            else if ( path[i][j] == 2 ) {
                targetMaxRegionReversed += targetSequence.charAt(i-1);
                queryMaxRegionReversed += "-";
                i--;
            }
            else {
                targetMaxRegionReversed += "-";
                queryMaxRegionReversed += querySequence.charAt(i);
                j--;
            }
        }

        i = targetMaxRegionReversed.length() - 1;
        String targetMaxRegion = "";
        String queryMaxRegion = "";
        while ( 0 <= i ) {
            targetMaxRegion += targetMaxRegionReversed.charAt(i);
            queryMaxRegion += queryMaxRegionReversed.charAt(i);
            i--;
        }*/


        long endTime = System.currentTimeMillis();
        return endTime - startTime;
    }


}