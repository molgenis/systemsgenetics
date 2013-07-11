/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneNetwork.Service;

import java.io.EOFException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.jws.WebService;
import javax.jws.WebMethod;
import javax.jws.WebParam;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

/**
 *
 * @author R.S.N. Fehrmann
 */
@WebService(serviceName = "GeneNetworkService")
public class GeneNetworkService {

    boolean d_serviceReadyToProcessQuery = false;
    // Containing intialization parameters for GeneNetwork Service
    XMLConfiguration config;
    // Datastructure contains symmetricMatrices
    DataStructure d_dataStructure;    
    // Structure containing links between different gene identifiers
    GeneSynonyms d_geneSynonyms;

    @WebMethod(operationName = "parseXMLInitFile")
    public boolean parseXMLInitFile(@WebParam(name = "parseXMLInitFile") String path) {
        try {
            config = new XMLConfiguration(path);
        } catch (ConfigurationException ex) {
            Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
        }

        return true;
    }
    
    @WebMethod(operationName = "readDataIntoMemory")
    public void readDataIntoMemory() {
        this.readGeneSynonymsIntoMemory();
        this.readSymmetricMatricesIntoMemory();
        d_serviceReadyToProcessQuery = true;
    }
    
    @WebMethod(operationName = "readGeneSynonymsIntoMemory")
    public void readGeneSynonymsIntoMemory() {
        try {
            d_geneSynonyms = new GeneSynonyms(config.getString("supportdata.genesynonyms.path"));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @WebMethod(operationName = "readSymmetricMatricesIntoMemory")
    public boolean readSymmetricMatricesIntoMemory() {

        d_dataStructure = new DataStructure();

        ArrayList<String> pathList = (ArrayList<String>) config.getList("symmetricMatrices.symmetricMatrix.path");
        ArrayList<String> identifierList = (ArrayList<String>) config.getList("symmetricMatrices.symmetricMatrix.identifier");

        for (int i = 0; i < pathList.size(); ++i) {
            try {
                IOSymmetricMatrix<Short> ioSymMat = new IOSymmetricMatrix<Short>(new File(pathList.get(i)));

                SymmetricMatrix<Short> symmetricMatrix = new SymmetricMatrix<Short>();

                ioSymMat.readShorts(symmetricMatrix);

                d_dataStructure.addSymmetricMatrix(symmetricMatrix, identifierList.get(i), pathList.get(i));
            } catch (IllegalArgumentException ex) {
                Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
            } catch (EOFException ex) {
                Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(GeneNetworkService.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        return true;
    }
    
    @WebMethod(operationName = "getSizeOfDataStructure")
    public int getSizeOfDataStructure() {
        return d_dataStructure.size();
    }
    
    @WebMethod(operationName = "getPathOfSymmetricMatrix")
    public String getPathOfSymmetricMatrix(int index) {
        return d_dataStructure.getPath(index);
    }
    
    @WebMethod(operationName = "getIdentifierOfSymmetricMatrix")
    public String getIdentifierOfSymmetricMatrix(int index) {
        return d_dataStructure.getIdentifier(index);
    }
    
    @WebMethod(operationName = "testFunction")
    public String testFunction() {
        return "Works!";
    }
    
}
