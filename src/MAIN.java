package uk.ac.ebi.uniprot.dataservice.client.examples;

import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanOriginal.Original;
import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanStripedLayout.StripedLayout;


public class MAIN {

    public static void main(String[] args) throws ServiceException {
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415";
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        int lenD = D.length();
        String querySequenceName = "P49950";
        String Q = uniProtService.getEntry(querySequenceName).getSequence().getValue();
        int lenQ = Q.length();

        System.out.println("Stripped Layout: " + StripedLayout(D, lenD, Q, lenQ));
        System.out.println("Original: " + Original(D, lenD, Q, lenQ));

    }
}
