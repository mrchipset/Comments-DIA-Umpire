/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package Test;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.SearchResultParser.ProtXMLParser;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Test {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {

        Logger logger = Logger.getRootLogger();
        ConsoleAppender ca = new ConsoleAppender();
        ca.setThreshold(Level.INFO);
        ca.setName("ConsoleLogger_Info");
        ca.setLayout(new PatternLayout("%d %-5p [%c{1}] %m%n"));
        ca.activateOptions();

        logger.getLoggerRepository().resetConfiguration();
        logger.addAppender(ca);
//        String Filename = "F:/Data/SWATH_Umpire_Paper/APMS/SWATH/LongSwath_EIF4aJune7_Biorep3/interact-LongSwath_EIF4aJune7_Biorep3_Q1.tandem.pep.xml";
//        
//        LCMSID IDsummary = new LCMSID("");
//        
//        TPPResult tppresult = new TPPResult(0.05f, 0.05f,"rev");
//        tppresult.ReadSearchResult(IDsummary, Filename, Path + Filename.replace("pep.xml", "prot.xml"));
//        System.out.print("Protein No.:" + IDsummary.ProteinList.size() + "; Assigned Peptide No.:" + IDsummary.AssignedPepIonList.size() + "; All peptide No.:" + IDsummary.GetPepIonList().size() + "\n");
//
//        FileWriter writer = new FileWriter(Path + "Hela_pep.csv");
//        for (PepIonID pep : IDsummary.GetPepIonList().values()) {
//            writer.write(pep.ModSequence + "\t" + pep.NeutralPrecursorMz() + "\t" + pep.GetIDRT() + "\n");
//        }
//        writer.close();
        String Filename = "F:\\Data\\Test\\GPM_combined.pep.inter.iprot.prot.decoy.xml";
        String Filename2 = "F:\\Data\\Test\\interact-LongSwath_EIF4aJune7_Biorep1_Combine.iproph.prot.xml";
        String Filename3 = "F:\\Data\\Test\\LongSwath_EIF4aJune7_Biorep1.tandem.Qcombine.prot.xml";
        
        String fasta = "F:\\Data\\Test\\HEK293_RefV57_cRAPandREV_20130129.fasta";
        LCMSID protID = new LCMSID(Filename);
        protID.FastaPath = fasta;
        ProtXMLParser protxmlparser = new ProtXMLParser(protID, Filename, 0f);
        protID.DecoyTag = "DECOY";
        float rf = protID.GetRFactor(0.2f);
        
        protID.RemoveLowLocalPWProtein(0.5f);
        protID.FilterByProteinDecoyFDRUsingMaxIniProb(protID.DecoyTag, 0.01f / rf);
        //protID.FilterByProteinDecoyFDRUsingMaxLocalPW(ProtIDList, tandemPara.DecoyPrefix, tandemPara.ProtFDR / rf);            
        protID.GenerateIndisProtMap();
        protID.LoadSequence();
        logger.info("Protein No.:" + protID.ProteinList.size());

        protID = new LCMSID(Filename2);
        protID.FastaPath = fasta;
        protxmlparser = new ProtXMLParser(protID, Filename2, 0f);
        protID.DecoyTag = "DECOY";
        rf = protID.GetRFactor(0.2f);
        
        protID.RemoveLowLocalPWProtein(0.5f);
        protID.FilterByProteinDecoyFDRUsingMaxIniProb(protID.DecoyTag, 0.01f / rf);
        //protID.FilterByProteinDecoyFDRUsingMaxLocalPW(ProtIDList, tandemPara.DecoyPrefix, tandemPara.ProtFDR / rf);            
        protID.GenerateIndisProtMap();
        protID.LoadSequence();
        logger.info("Protein No.:" + protID.ProteinList.size());
        
        protID = new LCMSID(Filename3);
        protID.FastaPath = fasta;
        protxmlparser = new ProtXMLParser(protID, Filename3, 0f);
        protID.DecoyTag = "DECOY";
        rf = protID.GetRFactor(0.2f);
        
        protID.RemoveLowLocalPWProtein(0.5f);
        protID.FilterByProteinDecoyFDRUsingMaxIniProb(protID.DecoyTag, 0.01f / rf);
        //protID.FilterByProteinDecoyFDRUsingMaxLocalPW(ProtIDList, tandemPara.DecoyPrefix, tandemPara.ProtFDR / rf);            
        protID.GenerateIndisProtMap();
        protID.LoadSequence();
        logger.info("Protein No.:" + protID.ProteinList.size());
    }
}
