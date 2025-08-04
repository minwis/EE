package uk.ac.ebi.uniprot.dataservice.client.examples;

import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.model.factories.DefaultUniProtFactory;
import uk.ac.ebi.kraken.parser.UniProtParser;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanOriginal.Original;
import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanStripedLayout.StripedLayout;


public class MAIN {

    // load from file (***)
    public static String getSequenceFromFile(String filename) {
        File f = new File(filename);
        System.out.println("Input Absolute path: " + f.getAbsolutePath());
        System.out.println("File size in bytes: " + f.length());
        try (FileInputStream fis = new FileInputStream(f)) {
            UniProtEntry eD = UniProtParser.parse(fis, DefaultUniProtFactory.getInstance(), true);
            return eD.getSequence().getValue();
        } catch (IOException e) {
            e.printStackTrace();
            return "";
        }
    }

    public static void main(String[] args) throws ServiceException {

        long timeForStripped = 0;
        long timeForOriginal = 0;
/*
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

       String targetProteinName = "P10415"; //A0A1B0GTW7
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        String querySequenceName = "P49950"; //A0A1L8HYT7
        String Q = uniProtService.getEntry(querySequenceName).getSequence().getValue();
*/
        // load from file (***)
        //String D = getSequenceFromFile("P10415.txt");
        //String Q = getSequenceFromFile("P49950.txt");


        String D = "TGTTACGG";
        String Q = "GGTTGACTA";

        int lenD = D.length();
        int lenQ = Q.length();

        System.out.println("Original Time(ns): " + Original(D, lenD, Q, lenQ));
        System.out.println("Stripped-Layout Time(ns): " + StripedLayout(D, lenD, Q, lenQ));

        /*
        String targetProteinName = "Q8WZ42-1";
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int lenD = D.length();

        Query query = UniProtQueryBuilder.proteinName("Titin");
        QueryResult<UniProtEntry> entries = uniProtService.getEntries(query);



        for ( int i = 0; i < 5; i++ )  { //(entries.hasNext()) Vector ==> arr/list ==> scalar

            UniProtEntry entry = entries.next();
            AlternativeProductsComment altProdComment = null;
            for (Comment comment : entry.getComments()) {
                if (comment.getCommentType() == CommentType.ALTERNATIVE_PRODUCTS) {
                    altProdComment = (AlternativeProductsComment) comment;
                    break;
                }
            }

            String Q = uniProtService.getEntry(String.valueOf(entry.getUniProtId())).getSequence().getValue();
            int lenQ = Q.length();

            timeForStripped+=StripedLayout(D, lenD, Q, lenQ);
            timeForOriginal+=Original(D, lenD, Q, lenQ);

            if ( altProdComment != null ) {
                List<AlternativeProductsIsoform> isoforms = altProdComment.getIsoforms();
                if ( !isoforms.isEmpty() ) {
                    for ( int j = 0; j < isoforms.toArray().length; j++ ) {
                        Q = uniProtService.getEntry(isoforms.get(i).getIds().getFirst().getValue()).getSequence().getValue();
                        lenQ = Q.length();

                        timeForStripped+=StripedLayout(D, lenD, Q, lenQ);
                        timeForOriginal+=Original(D, lenD, Q, lenQ);
                    }
                }
            }

        }


        for ( int i = 0; i < 4; i++ ) {
            System.out.println("Stripped Layout: " + timeForStripped + "\nOriginal: " + timeForOriginal );
        }

        */


    }
}
