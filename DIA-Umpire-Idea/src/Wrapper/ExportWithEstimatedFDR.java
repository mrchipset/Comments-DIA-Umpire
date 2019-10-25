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
package Wrapper;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.SearchResultParser.ProtXMLParser;
import Utility.ConsoleLogger;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class ExportWithEstimatedFDR {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException, ClassNotFoundException, XmlPullParserException, InterruptedException, SQLException {

        System.out.println("=================================================================================================");
        System.out.println("Umpire search result parser(version: v1.0, 2014.11)");
        if (args.length == 0) {
            System.out.println("command : java -jar -Xmx1G Umpire-SearchResultParser.jar [Options] [Combined ProtXML file] [PepXML files...]");
            System.out.println("");
            System.out.println("ProtXML extension: *.prot.xml or *.ProtXML");
            System.out.println("PepXML extension: *.pep.xml or *.PepXML");
            System.out.println("\nOptions");
            System.out.println("\t-fP\tProtein FDR\tex: -fP0.01 (default: 0.01, no filtering: -1)");
            System.out.println("\t-fp\tPeptide FDR\tex: -fp0.05 (default: 0.01, no filtering: -1)");
            System.out.println("\t-d\tDecoy tag prefix\tex: -dDECOY (default: rev_)");
            System.out.println("\t-fa\tFasta file");
            System.out.println("\t-N\tOutput filename");
            System.out.println("\t-pt\tInitial protein probability filtering threshold\tex: -pt0.5 (default: 0.5, no filtering : -1)");
            System.out.println("\t-rf\tR factor threshold, proteins with protein probablity less than the threshold will be used to estimate the R factor \n\t\tex: -rf0.2 (default: 0.2, do not use R factor: -1)");
            return;
        }

        ConsoleLogger.SetConsoleLogger(Level.INFO);
        ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "parser_debug.log");
        
        float protFDR = 0.01f;
        float pepFDR = 0.01f;
        String DecoyTag = "rev_";
        String Fasta = "";
        String Outputname = "";
        float protprob = 0.5f;
        float rfthreshold = 0.2f;
        String ProtXML = "";
        ArrayList<String> PepXML = new ArrayList<>();

        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("-")) {
                if (args[i].startsWith("-fP")) {
                    protFDR = Float.parseFloat(args[i].substring(3));
                    Logger.getRootLogger().info("Protein FDR:" + protFDR);
                }
                if (args[i].startsWith("-fp")) {
                    pepFDR = Float.parseFloat(args[i].substring(3));
                    Logger.getRootLogger().info("Peptide FDR:" + pepFDR);
                }
                if (args[i].startsWith("-d")) {
                    DecoyTag = args[i].substring(2);
                    Logger.getRootLogger().info("Decoy tag: " + DecoyTag);
                }
                if (args[i].startsWith("-fa")) {
                    Fasta = args[i].substring(3);
                    Logger.getRootLogger().info("Fasta file: " + Fasta);
                }
                if (args[i].startsWith("-N")) {
                    Outputname = args[i].substring(2);
                    Logger.getRootLogger().info("Output filename: " +Outputname);
                }
                if (args[i].startsWith("-pt")) {
                    protprob = Float.parseFloat(args[i].substring(3));
                    Logger.getRootLogger().info("Initial protein probablity filtering threshold: " + protprob);
                }
                if (args[i].startsWith("-rf")) {
                    rfthreshold = Float.parseFloat(args[i].substring(3));
                    Logger.getRootLogger().info("R factor threshold: " + rfthreshold);
                }
            }
            if (args[i].endsWith(".pep.xml") || args[i].endsWith(".PepXML")) {
                PepXML.add(args[i]);
            }
            if (args[i].endsWith(".prot.xml") || args[i].endsWith(".ProtXML")) {
                ProtXML = args[i];
            }
        }
        LCMSID lcmsid = new LCMSID(Outputname);
        lcmsid.FastaPath = Fasta;
        for (String pepxml : PepXML) {
            LCMSID pepxmlid = new LCMSID(pepxml);
            PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);
            if (pepFDR != -1f) {
                pepxmlid.FilterByPepDecoyFDR(DecoyTag, pepFDR);
            }
            Logger.getRootLogger().info("peptide No.:" + pepxmlid.GetPepIonList().size() + "; Peptide level threshold: " + pepxmlid.PepProbThreshold);
            for (PepIonID pepID : pepxmlid.GetPepIonList().values()) {
                lcmsid.AddPeptideID(pepID);
            }
        }

        if (!"".equals(ProtXML)) {
            ProtXMLParser protxmlparser = new ProtXMLParser(lcmsid, ProtXML, 0f);
            lcmsid.DecoyTag = DecoyTag;
            if (protprob != -1f) {
                lcmsid.RemoveLowLocalPWProtein(protprob);
            }
            float rf = 1f;
            if (rfthreshold != -1f) {
                rf = lcmsid.GetRFactor(rfthreshold);
            }
            if (protFDR != -1f) {
                lcmsid.FilterByProteinDecoyFDRUsingMaxIniProb(lcmsid.DecoyTag, protFDR / rf);
            }
            if (!"".equals(Fasta)) {
                lcmsid.LoadSequence();
            }
            lcmsid.ReMapProPep();
        }
        lcmsid.CreateInstanceForAllPepIon();
        lcmsid.ExportProtID(null);
        lcmsid.ExportPepID(null);
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; All peptide ions.:" + lcmsid.GetPepIonList().size());
    }

}
