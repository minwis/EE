package uk.ac.ebi.uniprot.dataservice.client.examples;

import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.model.factories.DefaultUniProtFactory;
import uk.ac.ebi.kraken.parser.UniProtParser;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanOriginal.Original;
import static uk.ac.ebi.uniprot.dataservice.client.examples.SmithWatermanStripedLayout.StripedLayout;


public class MAIN {

    public static void main(String[] args) throws ServiceException {

        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService();

        String targetProteinName = "P10415"; //A0A1B0GTW7
        String D = uniProtService.getEntry(targetProteinName).getSequence().getValue();
        String querySequenceName = "P49950"; //A0A1L8HYT7
        String Q = uniProtService.getEntry(querySequenceName).getSequence().getValue();

        int lenD = D.length();
        int lenQ = Q.length();
        System.out.println("Size: " + lenD * lenQ + ", D : " + lenD + ", Q : " + lenQ);


        System.out.println("Original Time(ns): " + Original(D, lenD, Q, lenQ));
        System.out.println("Stripped-Layout Time(ns): " + StripedLayout(D, lenD, Q, lenQ));


    }
}
