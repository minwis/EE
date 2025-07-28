package uk.ac.ebi.uniprot.dataservice.client.examples;

import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.AlternativeProductsComment;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.AlternativeProductsIsoform;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.Comment;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.CommentType;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.QueryResult;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;
import uk.ac.ebi.uniprot.dataservice.query.Query;

import java.util.List;

import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanOriginal.Original;
import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanStripedLayout.StripedLayout;

public class test {

    public static void main(String[] args) throws ServiceException {

        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415"; //A0A1B0GTW7
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int lenD = D.length();
        String querySequenceName = "P49950"; //A0A1L8HYT7
        String Q = uniProtService.getEntry(querySequenceName).getSequence().getValue();
        int lenQ = Q.length();

        System.out.println("Stripped Layout: " + StripedLayout(D, lenD, Q, lenQ) + "\nOriginal: " + Original(D, lenD, Q, lenQ));

        /*
        Query query = UniProtQueryBuilder.proteinName("Titin");
        QueryResult<UniProtEntry> entries = uniProtService.getEntries(query);

        String canonicalSequence;
        String isoformSequence;
        for ( int i = 0; i < 1000; i++ )  { //(entries.hasNext())

            UniProtEntry entry = entries.next();
            canonicalSequence = uniProtService.getEntry(entry.getPrimaryUniProtAccession().getValue()).getSequence().getValue();

            AlternativeProductsComment altProdComment = null;
            for (Comment comment : entry.getComments()) {
                if (comment.getCommentType() == CommentType.ALTERNATIVE_PRODUCTS) {
                    altProdComment = (AlternativeProductsComment) comment;
                    break;
                }
            }



            if ( altProdComment != null ) {
                List<AlternativeProductsIsoform> isoforms = altProdComment.getIsoforms();
                System.out.println(isoforms.toArray().length);
                if ( !isoforms.isEmpty() ) {
                    for ( int j = 0; j < isoforms.toArray().length; j++ ) {
                        isoformSequence = uniProtService.getEntry(isoforms.get(i).getIds().getFirst().getValue()).getSequence().getValue();
                    }
                }
            }



        }*/

        }

}
