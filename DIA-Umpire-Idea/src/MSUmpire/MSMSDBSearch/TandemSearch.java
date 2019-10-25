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
package MSUmpire.MSMSDBSearch;

import Utility.PrintThread;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TandemSearch extends MSMSDBSearch{

    TandemParam parameter;

    public TandemSearch(TandemParam parameter) {
        this.parameter = parameter;
    }

    public void RunTandem() throws IOException, InterruptedException {

        if (parameter.SpectrumPath.toLowerCase().endsWith(".mgf")) {
            if (!new File(parameter.SpectrumPath.replace(".mgf", ".mzXML")).exists()) {
                String[] msconvertcmd = {parameter.msconvertpath, "--mzXML", "--32", "-z", parameter.SpectrumPath, "-o", FilenameUtils.getFullPath(parameter.SpectrumPath)};
                Process p = Runtime.getRuntime().exec(msconvertcmd);
                Logger.getRootLogger().info("MGF file coversion by msconvert.exe...." + parameter.SpectrumPath );
                Logger.getRootLogger().debug("Command: " + Arrays.toString(msconvertcmd));
                PrintThread printThread=new PrintThread(p);
                printThread.start();
                p.waitFor();
                if (p.exitValue() != 0) {
                    Logger.getRootLogger().info("msconvert : " + parameter.SpectrumPath + " failed");
                    //PrintOutput(p);
                    return;
                }
                //PrintOutput(p);
            }
            parameter.SpectrumPath = parameter.SpectrumPath.replace(".mgf", ".mzXML");
        }
        parameter.GenerateParamFile();

        if (parameter.Overwrite || !(new File(parameter.RawSearchResult)).exists() || (new File(parameter.RawSearchResult)).length() / 1024 < 2) {
            String[] tandemcmd = {parameter.tandempath, parameter.parameterPath};
            Process p = Runtime.getRuntime().exec(tandemcmd);
            Logger.getRootLogger().info("X!Tandem search...." + parameter.RawSearchResult );
            Logger.getRootLogger().debug("Command: " + Arrays.toString(tandemcmd));
            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("tandem : " + parameter.RawSearchResult + " failed");
                //PrintOutput(p);
                return;
            }
        }

        if (parameter.Overwrite || !(new File(parameter.RawSearchResult + ".pep.xml")).exists()) {
            String[] tandemXMLcmd = {parameter.tandem2XML, parameter.RawSearchResult, parameter.PepXMLPath};
            Process p = Runtime.getRuntime().exec(tandemXMLcmd);
            Logger.getRootLogger().info("Converting .tandem to .pep.xml by tandem2XML...." + parameter.PepXMLPath);
            Logger.getRootLogger().debug("Command: " + Arrays.toString(tandemXMLcmd));
            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("tandem2XML : " + parameter.RawSearchResult + " failed");
                //PrintOutput(p);
                return;
            }
        }

        if (parameter.Overwrite || !(new File(parameter.InteractPepXMLPath)).exists()) {
            new File(parameter.InteractPepXMLPath).delete();
            new File(parameter.ProtXMLPath).delete();

            Process p = null;
            ArrayList<String> cmdlist = new ArrayList<>();
            cmdlist.add(parameter.xinteractpath);
            for (String op : parameter.xinteractpara.split(" ")) {
                cmdlist.add(op);
            }
            cmdlist.add(parameter.PepXMLPath);
            cmdlist.add("-N" + parameter.InteractPepXMLPath);
            String[] xinteractcmd = new String[cmdlist.size()];
            xinteractcmd = cmdlist.toArray(xinteractcmd);
            p = Runtime.getRuntime().exec(xinteractcmd);

            Logger.getRootLogger().info("xinteract PeptideProphet/ProteinProphet analysis...." + parameter.InteractPepXMLPath + " Option:" + parameter.xinteractpara );
            Logger.getRootLogger().debug("Command: " + Arrays.toString(xinteractcmd));
            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("xinteract : " + parameter.PepXMLPath + " failed");
                //PrintOutput(p);
                return;
            }
        }
    }

    @Override
    public void DBSearch() {
        try {
            RunTandem();
        } catch (IOException | InterruptedException ex) {
            Logger.getRootLogger().debug(ex.getMessage());
        }
    }

//    @Override
//    public void GenerateParamFile() {
//        this.parameter.GenerateParamFile();
//    }
//
//    @Override
//    public void SetResultFilePath(String mzXMLfile) {
//        parameter.SetResultFilePath(mzXMLfile);
//    }

    @Override
    public DBSearchParam GetParameter() {
       return parameter;
    }

    public MSMSDBSearch Clone() {
        try {
            TandemParam newparameter = (TandemParam) parameter.clone();
            TandemSearch newsearch = new TandemSearch(newparameter);
            return newsearch;
        } catch (CloneNotSupportedException ex) {
            Logger.getRootLogger().debug(ex.getMessage());
        }
        return null;
    }
}
