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
import MSUmpire.SearchResultParser.PepXMLParser;
import Utility.ConsoleLogger;
import Utility.PrintThread;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class GenerateLibUsingSpectraST {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, InterruptedException {
        ConsoleLogger.SetConsoleLogger(Level.DEBUG);
        String Path = "F:\\Data\\ETH_SEC_HEK\\DIA_SpecLib\\";
        String combinelib = "CombineLib";
        File folder = new File(Path);
        ArrayList<String> spliblist = new ArrayList<>();

        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.getName().endsWith(".pep.xml")) {
                LCMSID IDsummary = new LCMSID(fileEntry.getName());
                PepXMLParser pepxml = new PepXMLParser(IDsummary, fileEntry.getAbsolutePath(), 0f);
                IDsummary.DecoyTag = "DECOY";
                IDsummary.FDR = 0.01f;
                IDsummary.FindPepProbThresholdByFDR();
                float threshold = IDsummary.PepProbThreshold;

                Process p = null;
                ArrayList<String> cmdlist = new ArrayList<>();
                cmdlist.add("spectrast");
                cmdlist.add("-c_BIN!");
                cmdlist.add("-c_NPK1");
                cmdlist.add("-cf'Protein!~" + IDsummary.DecoyTag + "'");
                cmdlist.add("-cN\"" + fileEntry.getParent() + "/" + FilenameUtils.getBaseName(fileEntry.getName()) + "_filtered\"");
                spliblist.add(fileEntry.getParent() + "/" + FilenameUtils.getBaseName(fileEntry.getName()) + "_filtered" + ".splib");
                cmdlist.add("-cP" + threshold);
                cmdlist.add(fileEntry.getAbsolutePath());
                String[] cmd = new String[cmdlist.size()];
                cmd = cmdlist.toArray(cmd);
                p = Runtime.getRuntime().exec(cmd);
                Logger.getRootLogger().debug("Command: " + Arrays.toString(cmd));
                PrintThread printThread = new PrintThread(p);
                printThread.start();
                p.waitFor();
                if (p.exitValue() != 0) {
                    Logger.getRootLogger().info("spectrast : " + fileEntry.getAbsolutePath() + " failed");
                    return;
                }
            }
        }
        Process p = null;
        ArrayList<String> cmdlist = new ArrayList<>();
        cmdlist.add("spectrast");
        cmdlist.add("-cAC");
        cmdlist.add("-cJU");
        cmdlist.add("-c_DIS!");
        cmdlist.add("-c_QUO0.1");
        cmdlist.add("-cN" + Path + combinelib);
        for (String lib : spliblist) {
            cmdlist.add(lib);
        }
        String[] cmd = new String[cmdlist.size()];
        cmd = cmdlist.toArray(cmd);
        p = Runtime.getRuntime().exec(cmd);
        Logger.getRootLogger().debug("Command: " + Arrays.toString(cmd));
        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("spectrast : combine library failed");
        }
    }
}
