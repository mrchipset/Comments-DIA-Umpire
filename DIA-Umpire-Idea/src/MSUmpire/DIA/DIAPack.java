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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.FragmentLib.FragmentLibManager;
import MSUmpire.LCMSBaseStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.MSMSDBSearch;
import MSUmpire.MSMSDBSearch.ProteinProphetCombine;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.MSMSDBSearch.iProphet;
import MSUmpire.MySQLTool.ConnectionManager;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.ProtID;
import MSUmpire.PSMDataStructure.SortedPepListMass;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.SearchResultParser.TPPResult;
import MSUmpire.SpectraST.SpectraSTSearch;
import MSUmpire.UmpireSearchDataStructure.PepIonCandidate;
import MSUmpire.UmpireSearchDataStructure.PepIonLib;
import MSUmpire.UmpireSearchDataStructure.SortedPepCandidate;
import MSUmpire.spectrumparser.DIA_Setting;
import MSUmpire.spectrumparser.mzXMLParser;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.sql.Connection;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIAPack {

    private InstrumentParameter parameter;
    private ConnectionManager connectionManager;
    public String Filename;
    private int NoCPUs = 4;
    public LCMSPeakMS1 ms1lcms;
    public ArrayList<LCMSPeakDIAMS2> DIAWindows;
    public DIA_Setting dIA_Setting = new DIA_Setting();
    private mzXMLParser mzXML;
    public LCMSID DDAIDsummary;
    public LCMSID IDsummary;
    //public HashMap<Integer,Integer> MatchedClusters=new  HashMap<>();
    public HashMap<Integer, Integer> ScanClusterMap_Q1;
    public HashMap<Integer, Integer> ScanClusterMap_Q2;
    public HashMap<Integer, String> ScanClusterMap_Q3;
    public ArrayList<String> iProphPepXMLs = new ArrayList<String>();
    public boolean ExportPrecursorPeak = false;
    public boolean ExportFragmentPeak = false;
    public boolean ExportPeakClusterTable = false;
    public HashMap<Integer, Double> FactorialTable;
    TargetMatchScoring TScoring;
    public DIAStatus status=new DIAStatus();
    
    public class DIAStatus{
        public boolean SignalExtraction=false;
        public boolean UntargetedQuant=false;
        public boolean BuildMappedPep=false;
        public boolean TargetedQuant=false;
    }
    
    PepIonLib IonLib;
    public int Q1Scan = 0;
    public int Q2Scan = 0;
    public int Q3Scan = 0;

    public boolean UseMappedIon = false;
    public boolean FilterMappedIonByProb = true;
    public float MappedIonProbThreshold = 0.95f;
    
    public DIAPack(String Filename, int NoCPUs) throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        super();
        this.Filename = Filename;
        this.NoCPUs = NoCPUs;
    }

    public String GetBaseName(){
        return FilenameUtils.getBaseName(Filename);
    }
    
    public void SetNoCPUs(int cpu) {
        this.NoCPUs = cpu;
    }

    public void SetMySQLConnection(ConnectionManager connectionManager) {
        this.connectionManager = connectionManager;
    }

    public void SetParameter(InstrumentParameter parameter) {
        this.parameter = parameter;
    }

    public void SetDataType(SpectralDataType.DataType datatype) {
        this.dIA_Setting.dataType = datatype;
    }

    public void AddVaribleWindow(XYData window) {
        if (dIA_Setting.DIAWindows == null) {
            dIA_Setting.DIAWindows = new TreeMap<>();
        }
        dIA_Setting.DIAWindows.put(window, new ArrayList<Integer>());
    }

    public void AddMS1Window(XYData window) {
        if (dIA_Setting.MS1Windows == null) {
            dIA_Setting.MS1Windows = new TreeMap<>();
        }
        dIA_Setting.MS1Windows.put(window, new ArrayList<Integer>());
    }

    public void SetWindowSize(float size) {
        dIA_Setting.F_DIA_WindowSize = size;
    }

    public mzXMLParser GetMzXML() {

        if (mzXML == null) {
            try {
                mzXML = new mzXMLParser(Filename, parameter, dIA_Setting.dataType, dIA_Setting, NoCPUs);
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                Logger.getRootLogger().error("Read mzXML file:" + Filename + " failed.");
                System.exit(2);
            }
            dIA_Setting = mzXML.dIA_Setting;
            if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_diasetting.ser").exists()) {
                SaveDIASetting();
            }
        }
        return mzXML;
    }

    public void process() throws SQLException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, FileNotFoundException, DataFormatException, Exception {
        BuildDIAWindows();
        MS1PeakDetection();
        DIAMS2PeakDetection();
    }

    public void BuildDIAWindows() throws IOException, DataFormatException, IOException, IOException, IOException, InterruptedException {
        DIAWindows = new ArrayList<>();
        for (Map.Entry entry : GetMzXML().dIA_Setting.DIAWindows.entrySet()) {
            XYData DiaWinMz = (XYData) entry.getKey();
            LCMSPeakDIAMS2 diawindow = new LCMSPeakDIAMS2(Filename, this, parameter, DiaWinMz, GetMzXML(), NoCPUs);
            diawindow.SetMySQLConnection(connectionManager);
            diawindow.datattype = dIA_Setting.dataType;
            diawindow.ExportFragmentPeak = ExportFragmentPeak;
            DIAWindows.add(diawindow);
        }
    }

    private void ReadFactorialTable() throws IOException {
        FactorialTable = new HashMap<>();
        FactorialTable.put(0, 0d);
        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/factorial.txt");
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        String line = "";
        while ((line = reader.readLine()) != null) {
            Integer a = Integer.parseInt(line.split("\t")[0]);
            double b = Float.parseFloat(line.split("\t")[1]);
            FactorialTable.put(a, b);
        }
    }

    public void Qsplit_iProphet(ArrayList<MSMSDBSearch> searches) throws InterruptedException, IOException {
        String mgfname1 = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf";
        String mgfname2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf";
        String mgfname3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf";
        ArrayList<String> iproPepXMLs = new ArrayList<>();
        for (MSMSDBSearch dbsearch : searches) {
            dbsearch.SetResultFilePath(mgfname1);
        }
        String ipropepxml = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".iproph.pep.xml";
        iProphet ipro = new iProphet(ipropepxml, searches);
        ipro.DoiProphetPepXML();
        if (new File(ipropepxml).exists()) {
            iproPepXMLs.add(ipropepxml);
        }

        for (MSMSDBSearch dbsearch : searches) {
            dbsearch.SetResultFilePath(mgfname2);
        }
        ipropepxml = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".iproph.pep.xml";

        ipro = new iProphet(ipropepxml, searches);
        ipro.DoiProphetPepXML();
        if (new File(ipropepxml).exists()) {
            iproPepXMLs.add(ipropepxml);
        }
        for (MSMSDBSearch dbsearch : searches) {
            dbsearch.SetResultFilePath(mgfname3);
        }
        ipropepxml = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".iproph.pep.xml";
        ipro = new iProphet(ipropepxml, searches);
        ipro.DoiProphetPepXML();
        if (new File(ipropepxml).exists()) {
            iproPepXMLs.add(ipropepxml);
        }
        ipro.DoiProphetCombineProteinXML(tparam, iproPepXMLs, GetCombineiProphProtXML());
    }

    public void iProphet(ArrayList<MSMSDBSearch> searches) throws InterruptedException, IOException {
        iProphet ipro = new iProphet(GetIPROPHETPepXML(), searches);
        ipro.DoiProphetByCombinePepXML();
        ipro.DoiProphetProteinXML();
    }

    public void DIAMS2DBSearch(MSMSDBSearch dbsearch) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, FileNotFoundException, ExecutionException, DataFormatException {
        DIAMS2DBSearch(dbsearch, "");
    }

    public void DIAMS2DBSearch(MSMSDBSearch dbsearch, String mgftag) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, FileNotFoundException, ExecutionException, DataFormatException {
        String mgfname1 = FilenameUtils.getFullPath(Filename) + GetQ1Name() + mgftag + ".mgf";
        String mgfname2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + mgftag + ".mgf";
        String mgfname3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + mgftag + ".mgf";
        ArrayList<String> PepXMLs = new ArrayList<>();
        ArrayList<String> InteractPepXMLs = new ArrayList<>();
        dbsearch.SetResultFilePath(mgfname1);
        PepXMLs.add(dbsearch.GetParameter().PepXMLPath);
        InteractPepXMLs.add(dbsearch.GetParameter().InteractPepXMLPath);
        dbsearch.DBSearch();
        dbsearch.SetResultFilePath(mgfname2);
        PepXMLs.add(dbsearch.GetParameter().PepXMLPath);
        InteractPepXMLs.add(dbsearch.GetParameter().InteractPepXMLPath);
        dbsearch.DBSearch();
        dbsearch.SetResultFilePath(mgfname3);
        PepXMLs.add(dbsearch.GetParameter().PepXMLPath);
        InteractPepXMLs.add(dbsearch.GetParameter().InteractPepXMLPath);
        dbsearch.DBSearch();
        dbsearch.SetCombineFileName(Filename, mgftag);
        dbsearch.CombineProteinProphet(InteractPepXMLs);
        //dbsearch.CombinePepXML(PepXMLs);
    }

    public void SWATHMS2SpecLibSearch(String FastaFile, float corr) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, FileNotFoundException, ExecutionException, DataFormatException {
        Logger.getRootLogger().info("Identifying proteins by SpectraST...");
        String mgfname1 = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + (int) Math.floor(corr * 10) + "_T1.mzXML";
        String mgfname2 = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + (int) Math.floor(corr * 10) + "_T2.mzXML";
        //String mgfname3=FilenameUtils.getBaseName(ScanCollectionName) + "_" + (int)Math.floor(corr*10) + "_T3.mgf";
        String mgfname4 = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + (int) Math.floor(corr * 10) + "_T4.mzXML";

        SpectraSTSearch specsearch = new SpectraSTSearch();
        specsearch.RunSpectraST(mgfname1);
        specsearch.RunSpectraST(mgfname2);
        specsearch.RunSpectraST(mgfname4);

    }

    public void ReadScanNoMapping() throws FileNotFoundException, IOException {
        ScanClusterMap_Q1 = new HashMap<>();
        ScanClusterMap_Q2 = new HashMap<>();
        ScanClusterMap_Q3 = new HashMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q1"));
        BufferedReader reader2 = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q2"));
        BufferedReader reader4 = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q3"));

        
        String line = "";
        int StartNo = 0;
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML").exists()) {
            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML"));
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    StartNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
        }

        line = reader.readLine();
        int offset = StartNo - Integer.parseInt(line.split("_")[0]);
        Integer ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
        Integer ClusterIndex = Integer.parseInt(line.split("_")[1]);
        ScanClusterMap_Q1.put(ScanNo, ClusterIndex);

        while ((line = reader.readLine()) != null) {
            ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
            ClusterIndex = Integer.parseInt(line.split("_")[1]);
            ScanClusterMap_Q1.put(ScanNo, ClusterIndex);
        }

        line = "";
        StartNo = 0;
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML").exists()) {
            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML"));            
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    StartNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
        }
        line = reader2.readLine();
        offset = StartNo - Integer.parseInt(line.split("_")[0]);
        ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
        ClusterIndex = Integer.parseInt(line.split("_")[1]);
        ScanClusterMap_Q2.put(ScanNo, ClusterIndex);

        while ((line = reader2.readLine()) != null) {
            ScanNo = Integer.parseInt(line.split("_")[0]) + offset;
            ClusterIndex = Integer.parseInt(line.split("_")[1]);
            ScanClusterMap_Q2.put(ScanNo, ClusterIndex);
        }

        
        line = "";
        StartNo = 0;
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML").exists()) {
            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML"));
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    StartNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
        }
        line = reader4.readLine();        
        if (line.split(";").length == 3) {
            offset = StartNo - Integer.parseInt(line.split(";")[0]);
            ScanNo = Integer.parseInt(line.split(";")[0]) + offset;
            String WindowClusterIndex = line.split(";")[1] + ";" + line.split(";")[2];
            ScanClusterMap_Q3.put(ScanNo, WindowClusterIndex);
        } else {
            String ClusterIndexS = line.split("_")[1];
            ScanNo = Integer.parseInt(line.split("_")[0]);
            ScanClusterMap_Q3.put(ScanNo, ClusterIndexS);
        }
        while ((line = reader4.readLine()) != null) {            
            if (line.split(";").length == 3) {
                ScanNo = Integer.parseInt(line.split(";")[0]) + offset;
                String WindowClusterIndex = line.split(";")[1] + ";" + line.split(";")[2];
                ScanClusterMap_Q3.put(ScanNo, WindowClusterIndex);
            } else {
                ScanNo = Integer.parseInt(line.split("_")[0]);
                String ClusterIndexS = line.split("_")[1];
                ScanClusterMap_Q3.put(ScanNo, ClusterIndexS);
            }
        }
        reader.close();
        reader2.close();
        reader4.close();
    }

    public String GetQ1Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q1";
        return FilenameUtils.getBaseName(Filename) + "_Q1";
    }

    public String GetQ2Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q2";
        return FilenameUtils.getBaseName(Filename) + "_Q2";
    }

    public String GetQ3Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q3";
        return FilenameUtils.getBaseName(Filename) + "_Q3";
    }

    public String GetQ1Pepxml() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(GetQ1Name()) + ".pep.xml");
    }

    public String GetQ2Pepxml() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(GetQ2Name()) + ".pep.xml");
    }

    public String GetQ3Pepxml() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(GetQ3Name()) + ".pep.xml");
    }

    public String GetIPROPHETPepXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + ".iproph.pep.xml";
    }

    public String GetIPROPHETProtXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + ".iproph.prot.xml";
    }

    public String GetCombinePepXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + "_Combine.pep.xml";
    }

    public String GetCombineProtXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + "_Combine.prot.xml";
    }

    public String GetCombineiProphProtXML() {
        return FilenameUtils.getFullPath(Filename) + "interact-" + FilenameUtils.getBaseName(Filename) + "_Combine.iproph.prot.xml";
    }

    public void GenerateClusterScanNomapping() throws IOException {

        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q1").exists()) {
            String mgfname1 = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf";
            String mgfname2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf";
            String mgfname4 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf";

            BufferedReader reader1 = new BufferedReader(new FileReader(mgfname1));
            BufferedReader reader2 = new BufferedReader(new FileReader(mgfname2));
            BufferedReader reader4 = new BufferedReader(new FileReader(mgfname4));

            FileWriter writer = new FileWriter(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q1");
            FileWriter writer2 = new FileWriter(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q2");
            FileWriter writer4 = new FileWriter(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + ".ScanClusterMapping_Q3");

            BufferedReader mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML"));
            String line = "";
            int ScanNo = 0;
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    ScanNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
            while ((line = reader1.readLine()) != null) {
                if (line.startsWith("TITLE=")) {
                    int ClusterIndex = Integer.parseInt(line.split("ClusterIndex:")[1].split(",")[0]);
                    writer.write(ScanNo + "_" + ClusterIndex + "\n");
                    ScanNo++;
                }
            }
            mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML"));
            line = "";
            ScanNo = 0;
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    ScanNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
            while ((line = reader2.readLine()) != null) {
                if (line.startsWith("TITLE=")) {
                    int ClusterIndex = Integer.parseInt(line.split("ClusterIndex:")[1].split(",")[0]);
                    writer2.write(ScanNo + "_" + ClusterIndex + "\n");
                    ScanNo++;
                }
            }

            mzReader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML"));
            line = "";
            ScanNo = 0;
            while ((line = mzReader.readLine()) != null) {
                if (line.contains("<scan num=")) {
                    String substr = line.substring(line.indexOf("<scan num=") + 11);
                    ScanNo = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    break;
                }
            }
            while ((line = reader4.readLine()) != null) {
                if (line.startsWith("TITLE=")) {
                    int ClusterIndex = Integer.parseInt(line.split("ClusterIndex:")[1].split(",")[0]);
                    String DIAwindow = line.split("TITLE=")[1].split(";")[0];
                    writer4.write(ScanNo + ";" + DIAwindow + ";" + ClusterIndex + "\n");
                    ScanNo++;
                }
            }
            writer.close();
            writer2.close();
            writer4.close();
            reader1.close();
            reader2.close();
            mzReader.close();
            reader4.close();
        }
        ReadScanNoMapping();
    }

    public void FilterDeamidateGlycoPep() throws IOException {
        HashSet<Integer> PeakCurveIndexMS1 = new HashSet<>();
        ArrayList<Integer> SharedPeakCurveMS1 = new ArrayList<>();
        HashSet<Integer> PeakCurveIndexMS2 = new HashSet<>();
        ArrayList<Integer> SharedPeakCurveMS2 = new ArrayList<>();
        ArrayList<Integer> IDPeakClusterMS1 = new ArrayList<>();
        ArrayList<Integer> IDPeakClusterMS2 = new ArrayList<>();
        ArrayList<String> AllPepSeq = new ArrayList<>();
        for (String PepSeq : IDsummary.PeptideList.keySet()) {
            if (!AllPepSeq.contains(PepSeq)) {
                AllPepSeq.add(PepSeq);
            }
        }
        for (String PepSeq : IDsummary.MappedPeptideList.keySet()) {
            if (!AllPepSeq.contains(PepSeq)) {
                AllPepSeq.add(PepSeq);
            }
        }
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (!pepIonID.MS1PeakClusters.isEmpty()) {
                pepIonID.GlycoMS1Valid = true;
            }
            if (!pepIonID.MS2UnfragPeakClusters.isEmpty()) {
                pepIonID.GlycoMS2Valid = true;
            }
        }
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            if (!pepIonID.MS1PeakClusters.isEmpty()) {
                pepIonID.GlycoMS1Valid = true;
            }
            if (!pepIonID.MS2UnfragPeakClusters.isEmpty()) {
                pepIonID.GlycoMS2Valid = true;
            }
        }

        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_GlycoSharePeakResult.xls");
        writer.write("MSlevel\tSharedPeakCurve\tPeptideIndex\tSequence\tPepModSeq\tMass\tCharge\tmz\tRT\tNoClusters\tClusterIndex\tIsoScore\tIDIsoScore\tValid\n");

        for (String PepSeq : AllPepSeq) {
            PeakCurveIndexMS1.clear();
            SharedPeakCurveMS1.clear();
            PeakCurveIndexMS2.clear();
            SharedPeakCurveMS2.clear();
            IDPeakClusterMS1.clear();
            IDPeakClusterMS2.clear();

            if (IDsummary.PeptideList.containsKey(PepSeq)) {
                for (PepIonID pepion : IDsummary.PeptideList.get(PepSeq).values()) {
                    AddShareIndexMS1(pepion, PeakCurveIndexMS1, SharedPeakCurveMS1);
                    AddShareIndexMS2(pepion, PeakCurveIndexMS2, SharedPeakCurveMS2);
                    for (PeakCluster cluster : pepion.MS1PeakClusters) {
                        cluster.IDIsoPatternProb = cluster.GetChiSquareProbByTheoIso(pepion.IsotopicDistrubtion(parameter.MinNoPeakCluster));
                        IDPeakClusterMS1.add(cluster.Index);
                    }
                    for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                        cluster.IDIsoPatternProb = cluster.GetChiSquareProbByTheoIso(pepion.IsotopicDistrubtion(parameter.MinNoPeakCluster));
                        IDPeakClusterMS2.add(cluster.Index);
                    }
                }
            }
            if (IDsummary.MappedPeptideList.containsKey(PepSeq)) {
                for (PepIonID pepion : IDsummary.MappedPeptideList.get(PepSeq).values()) {
                    for (PeakCluster cluster : pepion.MS1PeakClusters) {
                        cluster.IDIsoPatternProb = cluster.GetChiSquareProbByTheoIso(pepion.IsotopicDistrubtion(parameter.MinNoPeakCluster));
                        if (IDPeakClusterMS1.contains(cluster.Index)) {
                            pepion.GlycoMS1Valid = false;
                        }
                    }
                    for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                        cluster.IDIsoPatternProb = cluster.GetChiSquareProbByTheoIso(pepion.IsotopicDistrubtion(parameter.MinNoPeakCluster));
                        if (IDPeakClusterMS2.contains(cluster.Index)) {
                            pepion.GlycoMS2Valid = false;
                        }
                    }
                    if (pepion.GlycoMS1Valid) {
                        AddShareIndexMS1(pepion, PeakCurveIndexMS1, SharedPeakCurveMS1);
                    }
                    if (pepion.GlycoMS2Valid) {
                        AddShareIndexMS2(pepion, PeakCurveIndexMS2, SharedPeakCurveMS2);
                    }
                }
            }

            for (Integer sharedindex : SharedPeakCurveMS1) {
                //System.out.println("Share MS1 peak curve index :" +sharedindex);
                SortedPepListMass ShareIDs = new SortedPepListMass();
                if (IDsummary.PeptideList.containsKey(PepSeq)) {
                    for (PepIonID pepion : IDsummary.PeptideList.get(PepSeq).values()) {
                        for (PeakCluster cluster : pepion.MS1PeakClusters) {
                            AddShareID(cluster, sharedindex, pepion, ShareIDs);
                        }
                    }
                }
                if (IDsummary.MappedPeptideList.containsKey(PepSeq)) {
                    for (PepIonID pepion : IDsummary.MappedPeptideList.get(PepSeq).values()) {
                        for (PeakCluster cluster : pepion.MS1PeakClusters) {
                            AddShareID(cluster, sharedindex, pepion, ShareIDs);
                        }
                    }
                }
                boolean CorrectFound = false;
                float CorrectMass = 0f;
                for (PepIonID pepion : ShareIDs) {
                    writer.write("MS1\t" + sharedindex + "\t" + pepion.Index + "\t" + pepion.Sequence + "\t" + pepion.ModSequence + "\t" + pepion.CalcNeutralPepMass() + "\t" + pepion.Charge + "\t" + pepion.ObservedMz + "\t" + pepion.PeakRT + "\t" + pepion.MS1PeakClusters.size() + "\t");
                    PeakCluster targetcluster = null;

                    for (PeakCluster cluster : pepion.MS1PeakClusters) {
                        for (int i = 0; i < cluster.IsoPeakIndex.length; i++) {
                            Integer index = cluster.IsoPeakIndex[i];
                            if (Objects.equals(sharedindex, index)) {
                                writer.write(cluster.Index + "_" + cluster.mz[i] + ";");
                                if (targetcluster == null || targetcluster.IDIsoPatternProb < cluster.IDIsoPatternProb) {
                                    targetcluster = cluster;
                                }
                            }
                        }
                    }
                    if ((CorrectFound && pepion.NeutralPrecursorMz() > CorrectMass) || targetcluster.IDIsoPatternProb < 0.8f) {
                        pepion.GlycoMS1Valid = false;
                    } else {
                        CorrectFound = true;
                        CorrectMass = pepion.NeutralPrecursorMz();
                    }
                    writer.write("\t" + targetcluster.IsoMapProb + "\t" + targetcluster.IDIsoPatternProb + "\t" + pepion.GlycoMS1Valid + "\n");
                }
            }
            for (Integer sharedindex : SharedPeakCurveMS2) {
                //System.out.println("Share MS2 peak curve index :" + sharedindex);
                SortedPepListMass ShareIDs = new SortedPepListMass();
                if (IDsummary.PeptideList.containsKey(PepSeq)) {
                    for (PepIonID pepion : IDsummary.PeptideList.get(PepSeq).values()) {
                        for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                            AddShareID(cluster, sharedindex, pepion, ShareIDs);
                        }
                    }
                }

                if (IDsummary.MappedPeptideList.containsKey(PepSeq)) {
                    for (PepIonID pepion : IDsummary.MappedPeptideList.get(PepSeq).values()) {
                        for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                            AddShareID(cluster, sharedindex, pepion, ShareIDs);
                        }
                    }
                }

                /////MS2
                boolean CorrectFound = false;
                float CorrectMass = 0f;
                for (PepIonID pepion : ShareIDs) {
                    writer.write("MS2\t" + sharedindex + "\t" + pepion.Index + "\t" + pepion.Sequence + "\t" + pepion.ModSequence + "\t" + pepion.CalcNeutralPepMass() + "\t" + pepion.Charge + "\t" + pepion.ObservedMz + "\t" + pepion.PeakRT + "\t" + pepion.MS1PeakClusters.size() + "\t");

                    PeakCluster targetcluster = null;
                    for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                        for (int i = 0; i < cluster.IsoPeakIndex.length; i++) {
                            Integer index = cluster.IsoPeakIndex[i];
                            if (Objects.equals(sharedindex, index)) {
                                writer.write(cluster.Index + "_" + cluster.mz[i] + ";");
                                if (targetcluster == null || targetcluster.IDIsoPatternProb < cluster.IDIsoPatternProb) {
                                    targetcluster = cluster;
                                }
                            }
                        }
                    }
                    if ((CorrectFound && pepion.NeutralPrecursorMz() > CorrectMass) || targetcluster.IDIsoPatternProb < 0.8f) {
                        pepion.GlycoMS2Valid = false;
                    } else {
                        CorrectFound = true;
                        CorrectMass = pepion.NeutralPrecursorMz();
                    }
                    writer.write("\t" + targetcluster.IsoMapProb + "\t" + targetcluster.IDIsoPatternProb + "\t" + pepion.GlycoMS2Valid + "\n");
                }
            }
        }
        writer.close();

        HashMap<String, ArrayList<PepIonID>> DeglycoSeq = new HashMap<>();
        HashMap<Integer, Boolean> NIndex = new HashMap<>();
        for (String PepSeq : AllPepSeq) {
            DeglycoSeq.clear();
            if (IDsummary.PeptideList.containsKey(PepSeq)) {
                for (PepIonID pepion : IDsummary.PeptideList.get(PepSeq).values()) {
                    String deglyco = pepion.ModSequence.replace("[0.9840698(N)]", "") + "_" + pepion.Charge;
                    if (!DeglycoSeq.containsKey(deglyco)) {
                        DeglycoSeq.put(deglyco, new ArrayList<PepIonID>());
                    }
                    DeglycoSeq.get(deglyco).add(pepion);
                }
            }
            if (IDsummary.MappedPeptideList.containsKey(PepSeq)) {
                for (PepIonID pepion : IDsummary.MappedPeptideList.get(PepSeq).values()) {
                    String deglyco = pepion.ModSequence.replace("[0.9840698(N)]", "") + "_" + pepion.Charge;
                    if (!DeglycoSeq.containsKey(deglyco)) {
                        DeglycoSeq.put(deglyco, new ArrayList<PepIonID>());
                    }
                    DeglycoSeq.get(deglyco).add(pepion);
                }
            }

            for (ArrayList<PepIonID> list : DeglycoSeq.values()) {
                //Multiple glyco forms                
                NIndex.clear();
                for (int i = 0; i < list.get(0).Sequence.length(); i++) {
                    if (String.valueOf(list.get(0).Sequence.charAt(i)).equals("N")) {
                        NIndex.put(i + 1, false);
                    }
                }
                for (PepIonID pepIonID : list) {
                    for (Integer site : NIndex.keySet()) {
                        NIndex.put(site, false);
                    }
                    for (FragmentPeak frag : pepIonID.FragmentPeaks) {
                        String type = String.valueOf(frag.IonType.charAt(0));
                        Integer matchindex = Integer.parseInt(frag.IonType.substring(1));
                        if ("b".equals(type)) {
                            for (Integer index : NIndex.keySet()) {
                                if (index <= matchindex) {
                                    NIndex.put(index, true);
                                }
                            }
                        } else if ("y".equals(type)) {
                            for (Integer index : NIndex.keySet()) {
                                if (index <= matchindex) {
                                    NIndex.put(index, true);
                                }
                            }
                        }
                    }
                    pepIonID.NoNsite = NIndex.size();
                    pepIonID.NoNsiteFragObs = 0;
                    for (Boolean obs : NIndex.values()) {
                        if (obs) {
                            pepIonID.NoNsiteFragObs++;
                        }
                    }
                }
            }
        }
        FileWriter writer2 = new FileWriter(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_GlycoAllResult.xls");
        writer2.write("PeptideIndex\tSequence\tPepModSeq\tMass\tCharge\tmz\tRT\tMS1IsoScore\tMS1IDIsoScore\tMS2IsoScore\tMS2IDIsoScore\tGlycoMS1Valid\tGlycoMS2Valid\tNoSite\tNoSiteObs\n");
        for (PepIonID pepion : IDsummary.GetPepIonList().values()) {
            writer2.write("ID_" + pepion.Index + "\t" + pepion.Sequence + "\t" + pepion.ModSequence + "\t" + pepion.CalcNeutralPepMass() + "\t" + pepion.Charge + "\t" + pepion.ObservedMz + "\t" + pepion.PeakRT + "\t");
            PeakCluster targetcluster = null;
            float IsoMapProb = 0f;
            float IDIsoPatternProb = 0f;
            for (PeakCluster cluster : pepion.MS1PeakClusters) {
                if (targetcluster == null || targetcluster.IDIsoPatternProb < cluster.IDIsoPatternProb) {
                    targetcluster = cluster;
                }
            }
            if (targetcluster != null) {
                IsoMapProb = targetcluster.IsoMapProb;
                IDIsoPatternProb = targetcluster.IDIsoPatternProb;
            }
            writer2.write(IsoMapProb + "\t" + IDIsoPatternProb + "\t");
            IsoMapProb = 0f;
            IDIsoPatternProb = 0f;
            for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                if (targetcluster == null || targetcluster.IDIsoPatternProb < cluster.IDIsoPatternProb) {
                    targetcluster = cluster;
                }
            }
            if (targetcluster != null) {
                IsoMapProb = targetcluster.IsoMapProb;
                IDIsoPatternProb = targetcluster.IDIsoPatternProb;
            }
            writer2.write(IsoMapProb + "\t" + IDIsoPatternProb + "\t" + pepion.GlycoMS1Valid + "\t" + pepion.GlycoMS2Valid + "\t" + pepion.NoNsite + "\t" + pepion.NoNsiteFragObs + "\n");
        }
        for (PepIonID pepion : IDsummary.GetMappedPepIonList().values()) {
            writer2.write("ID_" + pepion.Index + "\t" + pepion.Sequence + "\t" + pepion.ModSequence + "\t" + pepion.CalcNeutralPepMass() + "\t" + pepion.Charge + "\t" + pepion.ObservedMz + "\t" + pepion.PeakRT + "\t");
            PeakCluster targetcluster = null;
            float IsoMapProb = 0f;
            float IDIsoPatternProb = 0f;
            for (PeakCluster cluster : pepion.MS1PeakClusters) {
                if (targetcluster == null || targetcluster.IDIsoPatternProb < cluster.IDIsoPatternProb) {
                    targetcluster = cluster;
                }
            }
            if (targetcluster != null) {
                IsoMapProb = targetcluster.IsoMapProb;
                IDIsoPatternProb = targetcluster.IDIsoPatternProb;
            }
            writer2.write(IsoMapProb + "\t" + IDIsoPatternProb + "\t");
            IsoMapProb = 0f;
            IDIsoPatternProb = 0f;
            for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
                if (targetcluster == null || targetcluster.IDIsoPatternProb < cluster.IDIsoPatternProb) {
                    targetcluster = cluster;
                }
            }
            if (targetcluster != null) {
                IsoMapProb = targetcluster.IsoMapProb;
                IDIsoPatternProb = targetcluster.IDIsoPatternProb;
            }
            writer2.write(IsoMapProb + "\t" + IDIsoPatternProb + "\t" + pepion.GlycoMS1Valid + "\t" + pepion.GlycoMS2Valid + "\t" + pepion.NoNsite + "\t" + pepion.NoNsiteFragObs + "\n");
        }

        writer2.close();

        Logger.getRootLogger().info("No. of ID peptide ions:" + IDsummary.GetPepIonList().size() + " No. of mapped ions:" + IDsummary.GetMappedPepIonList().size());
        ArrayList<PepIonID> removlist = new ArrayList<>();

        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (!pepIonID.GlycoMS1Valid && !pepIonID.GlycoMS2Valid) {
                removlist.add(pepIonID);
            }
        }
        for (PepIonID pepIonID : removlist) {
            IDsummary.GetPepIonList().remove(pepIonID.GetKey());
            IDsummary.PeptideList.get(pepIonID.Sequence).remove(pepIonID.GetKey());
        }
        removlist.clear();
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            if (!pepIonID.GlycoMS1Valid && !pepIonID.GlycoMS2Valid) {
                removlist.add(pepIonID);
            }
        }
        for (PepIonID pepIonID : removlist) {
            IDsummary.GetMappedPepIonList().remove(pepIonID.GetKey());
            IDsummary.MappedPeptideList.get(pepIonID.Sequence).remove(pepIonID.GetKey());
        }

        Logger.getRootLogger().info("No. of ID peptide ions:" + IDsummary.GetPepIonList().size() + " No. of mapped ions:" + IDsummary.GetMappedPepIonList().size());

//        removlist=new ArrayList<>();
//                
//        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
//            if (pepIonID.NoNsiteFragObs == 0) {
//                removlist.add(pepIonID);
//            }
//        }
//        for(PepIonID pepIonID : removlist){
//            IDsummary.GetPepIonList().remove(pepIonID.GetKey());
//            IDsummary.PeptideList.get(pepIonID.Sequence).remove(pepIonID.GetKey());
//        }
//        removlist.clear();
//        for(PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()){
//             if (pepIonID.NoNsiteFragObs == 0) {
//                removlist.add(pepIonID);
//            }
//        }
//        for(PepIonID pepIonID : removlist){
//            IDsummary.GetMappedPepIonList().remove(pepIonID.GetKey());
//            IDsummary.MappedPeptideList.get(pepIonID.Sequence).remove(pepIonID.GetKey());
//        }
//        
//        System.out.println("No. of ID peptide ions:"+IDsummary.GetPepIonList().size() +" No. of mapped ions:"+IDsummary.GetMappedPepIonList().size());
//        
        IDsummary.ClearProPeplist();
        IDsummary.AssignProtForPepIon();
        IDsummary.AssignProtForMappedIon();
        IDsummary.GenearteAssignIonList();
    }

    private void AddShareID(PeakCluster cluster, Integer sharedindex, PepIonID pepion, SortedPepListMass ShareIDs) {
        for (Integer index : cluster.IsoPeakIndex) {
            if (Objects.equals(sharedindex, index)) {
                ShareIDs.add(pepion);
            }
        }
    }

    private void AddShareIndexMS1(PepIonID pepion, HashSet<Integer> PeakCurveIndexMS1, ArrayList<Integer> SharedPeakCurveMS1) {
        for (PeakCluster cluster : pepion.MS1PeakClusters) {
            for (Integer index : cluster.IsoPeakIndex) {
                if (index > 0) {
                    if (!PeakCurveIndexMS1.contains(index)) {
                        PeakCurveIndexMS1.add(index);
                    } else {
                        SharedPeakCurveMS1.add(index);
                    }
                }
            }
        }
    }

    private void AddShareIndexMS2(PepIonID pepion, HashSet<Integer> PeakCurveIndexMS2, ArrayList<Integer> SharedPeakCurveMS2) {

        for (PeakCluster cluster : pepion.MS2UnfragPeakClusters) {
            for (Integer index : cluster.IsoPeakIndex) {
                if (index > 0) {
                    if (!PeakCurveIndexMS2.contains(index)) {
                        PeakCurveIndexMS2.add(index);
                    } else {
                        SharedPeakCurveMS2.add(index);
                    }
                }
            }
        }
    }

    public void MapPSM2Cluster() throws IOException, SQLException {
        GenerateClusterScanNomapping();
        for (PSM psm : IDsummary.PSMList.values()) {
            int ClusterIndex = -1;
            if (psm.GetRawNameString().equals(FilenameUtils.getBaseName(GetQ1Name()))) {
                ClusterIndex = ScanClusterMap_Q1.get(psm.ScanNo);
                PeakCluster Cluster = ms1lcms.PeakClusters.get(ClusterIndex - 1);
                Cluster.AssignedPepIon = psm.Sequence;
            } else if (psm.GetRawNameString().equals(FilenameUtils.getBaseName(GetQ2Name()))) {
                ClusterIndex = ScanClusterMap_Q2.get(psm.ScanNo);
                PeakCluster Cluster = ms1lcms.PeakClusters.get(ClusterIndex - 1);
                Cluster.AssignedPepIon = psm.Sequence;
            }
        }
    }

    public void AssignQuant() throws IOException, SQLException {
        AssignQuant(true);
    }

    public void AssignQuant(boolean export) throws IOException, SQLException {
        Logger.getRootLogger().info("Assign peak cluster to identified peptides");
        GenerateClusterScanNomapping();
        //ReadScanNomapping();
        ExecutorService executorPool = null;
        for (PeakCluster cluster : ms1lcms.PeakClusters) {
            cluster.Identified = false;
        }
        //UpdateProcess progress = new UpdateProcess();
        //progress.SetTotal(IDsummary.GetPepIonList().size());
        //Thread thread = new Thread(progress);
        //thread.start();
        //executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            pepIonID.MS1PeakClusters = new ArrayList<>();
            pepIonID.MS2UnfragPeakClusters = new ArrayList<>();
        }
        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
            DIA_window_Quant dia_w = new DIA_window_Quant(GetQ1Name(), GetQ2Name(), GetQ3Name(), ScanClusterMap_Q1, ScanClusterMap_Q2, ScanClusterMap_Q3, ms1lcms, DIAWindow, IDsummary, NoCPUs);
            //executorPool.execute(dia_w);
            dia_w.run();
        }
//        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }

        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //executorPool = Executors.newFixedThreadPool(1);
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            DIAAssignQuantUnit quantunit = new DIAAssignQuantUnit(pepIonID, ms1lcms, parameter);
            executorPool.execute(quantunit);
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }
        //thread = null;
        //progress.ClearMSG();
        //progress = null;

        //ms1lcms.DetermineMS1PeakClusterScore();
        //ms1lcms.ExportPeakCluster();
        if (export) {
            ExportID();
        }
    }

    public void AssignMappedPepQuant(boolean export, FragmentLibManager libManager) throws IOException, SQLException {
        AssignMappedPepQuant(export, libManager, 1.1f);
    }

    public void AssignMappedPepQuant(boolean export, FragmentLibManager libManager, float ReSearchProb) throws IOException, SQLException {
        if (IDsummary.GetMappedPepIonList().isEmpty()) {
            Logger.getRootLogger().error("There is no peptide ion for targeted re-extraction.");
            return;
        }
        GenerateClusterScanNomapping();
//        FragmentLibManager libManager = new FragmentLibManager(LibID,connectionManager);
//        libManager.ReadFromDB();        
        ExecutorService executorPool = null;

        //UpdateProcess progress = new UpdateProcess();
        //progress.SetTotal(IDsummary.GetMappedPepIonList().size());
        //Thread thread = new Thread(progress);
        //thread.start();
        TScoring=new TargetMatchScoring(Filename,libManager.LibID);
                
        Logger.getRootLogger().info("No. of identified peptide ions:" + IDsummary.GetPepIonList().size());
        Logger.getRootLogger().info("No. of mapped peptide ions:" + IDsummary.GetMappedPepIonList().size());
        ArrayList<PepIonID> SearchList = new ArrayList<>();
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            if (libManager.PeptideFragmentLib.containsKey(pepIonID.GetKey()) && libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() >= 3 && Math.max(pepIonID.MS1AlignmentLocalProbability, pepIonID.MS2AlignmentLocalProbability) < ReSearchProb) {
                pepIonID.CreateQuantInstance(parameter.MaxNoPeakCluster);
                pepIonID.MS1PeakClusters = new ArrayList<>();
                pepIonID.MS2UnfragPeakClusters = new ArrayList<>();
                pepIonID.MS1AlignmentLocalProbability = 0f;
                pepIonID.MS1AlignmentProbability = 0f;
                pepIonID.MS2AlignmentLocalProbability = 0f;
                pepIonID.MS2AlignmentProbability = 0f;
                SearchList.add(pepIonID);                
            }
        }
        Logger.getRootLogger().info("No. of searchable peptide ions:" + SearchList.size());

        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
            Logger.getRootLogger().info("Assigning clusters for peak groups in MS2 isolation window:" + FilenameUtils.getBaseName(DIAWindow.ScanCollectionName));

            if (!DIAWindow.ReadPeakCluster() || !DIAWindow.ReadPrecursorFragmentClu2Cur()) {
                Logger.getRootLogger().warn("Reading results for " + DIAWindow.ScanCollectionName + " failed");
                System.exit(2);
            }

            executorPool = Executors.newFixedThreadPool(NoCPUs);
            for (PepIonID pepIonID : SearchList) {
                if (DIAWindow.DIA_MZ_Range.getX() <= pepIonID.NeutralPrecursorMz() && DIAWindow.DIA_MZ_Range.getY() >= pepIonID.NeutralPrecursorMz()) {
                    if (libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() >= 3) {
                        UmpireSpecLibMatch matchunit = new UmpireSpecLibMatch(ms1lcms, DIAWindow, pepIonID, libManager.GetFragmentLib(pepIonID.GetKey()), libManager.GetDecoyFragmentLib(pepIonID.GetKey()), parameter);
                        executorPool.execute(matchunit);
                        TScoring.libTargetMatches.add(matchunit);
                    } else {
                        Logger.getRootLogger().warn("skipping " + pepIonID.GetKey() + ", it has only " + libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() + " matched fragments");
                    }
                }
            }
            
            for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
                if (libManager.PeptideFragmentLib.containsKey(pepIonID.GetKey()) && DIAWindow.DIA_MZ_Range.getX() <= pepIonID.NeutralPrecursorMz() && DIAWindow.DIA_MZ_Range.getY() >= pepIonID.NeutralPrecursorMz()) {
                    if (libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() >= 3) {
                        UmpireSpecLibMatch matchunit = new UmpireSpecLibMatch(ms1lcms, DIAWindow, pepIonID, libManager.GetFragmentLib(pepIonID.GetKey()), libManager.GetDecoyFragmentLib(pepIonID.GetKey()), parameter);
                        matchunit.IdentifiedPeptideIon = true;
                        executorPool.execute(matchunit);
                        TScoring.libIDMatches.add(matchunit);
                    } else {
                        Logger.getRootLogger().warn("skipping " + pepIonID.GetKey() + ", it has only " + libManager.GetFragmentLib(pepIonID.GetKey()).FragmentGroups.size() + " matched fragments");
                    }
                }
            }
            executorPool.shutdown();
            while (!executorPool.isTerminated()) {
            }
            DIAWindow.ClearAllPeaks();
        }

        TScoring.Process();
        TScoring=null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            DIAAssignQuantUnit quantunit = new DIAAssignQuantUnit(pepIonID, ms1lcms, parameter);
            executorPool.execute(quantunit);
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        if (export) {
            ExportID();
            if (connectionManager != null) {
                ExportMappedIDQuant();
            }
        }
    }

    public void ExtractFragmentForMS1PeakCluster() {
        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
            Logger.getRootLogger().info("Extracting fragments for MS1 clusters in MS2 isolation window:" + FilenameUtils.getBaseName(DIAWindow.ScanCollectionName));

            if (!DIAWindow.ReadPrecursorFragmentClu2Cur()) {
                Logger.getRootLogger().error("Reading results for " + DIAWindow.ScanCollectionName + " failed");
                System.exit(2);
            }
            for (PeakCluster cluster : ms1lcms.PeakClusters) {
                DIAWindow.ExtractFragmentForPeakCluser(cluster);
            }
        }
    }

    

    public void ReadIDFromDBAndSkylineQuant(String mzxmlname, String SkylineExport, LCMSID RefProtID) throws SQLException, IOException {

        this.IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(mzxmlname));
        this.IDsummary.FastaPath = RefProtID.FastaPath;
        this.IDsummary.ReadFromSkylineExport(SkylineExport);
        this.IDsummary.CreateProtByRefID(RefProtID);
        this.IDsummary.GenearteAssignIonList();

        connectionManager.CloseConnection();
    }

    private void AssignMS2Cluster() throws SQLException {
        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
            DIAWindow.ReadPeakCluster();
            for (PepIonID peptide : IDsummary.GetPepIonList().values()) {
                if (DIAWindow.DIA_MZ_Range.getX() <= peptide.ObservedMz && DIAWindow.DIA_MZ_Range.getY() >= peptide.ObservedMz) {

                    for (String indexint : peptide.MS2ClusIndex.split(";")) {
                        if (!"".equals(indexint)) {
                            try {
                                PeakCluster cluster = DIAWindow.PeakClusters.get(Integer.parseInt(indexint) - 1);
                                if (Math.abs(cluster.TargetMz() - peptide.ObservedMz) < 0.1f && cluster.Charge == peptide.Charge) {
                                    peptide.MS2UnfragPeakClusters.add(cluster);
                                }
                            } catch (Exception e) {
                                //System.err.println(peptide.ObservedMz+"\t"+DIAWindow.DIA_MZ_Range.getX()+"\t"+DIAWindow.DIA_MZ_Range.getY());
                            }
                        }
                    }
                }
                if (peptide.PeakRT == -1f) {
                    PeakCluster target = null;
                    for (PeakCluster cluster : peptide.MS2UnfragPeakClusters) {
                        if (target == null || cluster.PeakHeight[0] > target.PeakHeight[0]) {
                            target = cluster;
                        }
                    }
                    if (target != null) {
                        peptide.PeakRT = target.PeakHeightRT[0];
                        peptide.ObservedMz = target.mz[0];
                    }
                }
            }
            if (IDsummary.GetMappedPepIonList() != null) {
                for (PepIonID peptide : IDsummary.GetMappedPepIonList().values()) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= peptide.ObservedMz && DIAWindow.DIA_MZ_Range.getY() >= peptide.ObservedMz) {
                        for (String indexint : peptide.MS2ClusIndex.split(";")) {
                            if (!"".equals(indexint)) {
                                try {
                                    PeakCluster cluster = DIAWindow.PeakClusters.get(Integer.parseInt(indexint) - 1);
                                    if (Math.abs(cluster.TargetMz() - peptide.ObservedMz) < 0.1f && cluster.Charge == peptide.Charge) {
                                        peptide.MS2UnfragPeakClusters.add(cluster);
                                    }
                                } catch (Exception e) {
                                    //System.err.println(peptide.ObservedMz + "\t" + DIAWindow.DIA_MZ_Range.getX() + "\t" + DIAWindow.DIA_MZ_Range.getY());
                                }
                            }
                        }
                    }
                    if (peptide.PeakRT == -1f) {
                        PeakCluster target = null;
                        for (PeakCluster cluster : peptide.MS2UnfragPeakClusters) {
                            if (target == null || cluster.PeakHeight[0] > target.PeakHeight[0]) {
                                target = cluster;
                            }
                        }
                        if (target != null) {
                            peptide.PeakRT = target.PeakHeightRT[0];
                            peptide.ObservedMz = target.mz[0];
                        }
                    }
                }
            }
        }
        for (PepIonID peptide : IDsummary.GetMappedPepIonList().values()) {
            if (StringUtils.countMatches(peptide.MS2ClusIndex, ";") != peptide.MS2UnfragPeakClusters.size()) {
                System.err.println("MS2 peak cluster not found " + peptide.GetKey() + ";" + peptide.MS2ClusIndex);
            }
        }
        for (PepIonID peptide : IDsummary.GetPepIonList().values()) {
            if (StringUtils.countMatches(peptide.MS2ClusIndex, ";") != peptide.MS2UnfragPeakClusters.size()) {
                System.err.println("MS2 peak cluster not found " + peptide.GetKey() + ";" + peptide.MS2ClusIndex);
            }
        }
    }

    private void CheckPSM() {
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (pepIonID.GetPSMList().isEmpty()) {
                for (PSM psm : templcmsid.PSMList.values()) {
                    if ((psm.Sequence == null ? pepIonID.Sequence == null : psm.Sequence.equals(pepIonID.Sequence)) && psm.Charge == pepIonID.Charge) {
                        if (Math.abs(pepIonID.NeutralPrecursorMz() - psm.NeutralPrecursorMz()) < 0.1f) {
                            psm.ModSeq = pepIonID.ModSequence;
                            pepIonID.AddPSM(psm);
                        }
                    }
                }
                if (pepIonID.GetPSMList().isEmpty()) {
                    Logger.getRootLogger().error("PSM list is empty:" + pepIonID.GetKey());
                }
            }
            if (Math.abs(pepIonID.ObservedMass() - pepIonID.CalcNeutralPepMass()) > 0.1f) {
                Logger.getRootLogger().error(pepIonID.ObservedMass());
                float mass = pepIonID.CalcNeutralPepMass();
                Logger.getRootLogger().error(mass);
                pepIonID.GetPepFactory().estimateTheoreticMass();
                mass = pepIonID.CalcNeutralPepMass();
                Logger.getRootLogger().error(mass);
            }
        }
    }

    private void CheckOverlap() {
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (IDsummary.GetMappedPepIonList().containsKey(pepIonID.GetKey())) {
                Logger.getRootLogger().error("overlap");
            }
        }
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
            if (IDsummary.GetPepIonList().containsKey(pepIonID.GetKey())) {
                Logger.getRootLogger().error("overlap");
            }
        }
    }

    public boolean ReadSerializedLCMSID() throws Exception {
        return ReadSerializedLCMSID("");
    }
    public boolean ReadSerializedLCMSID(String tag) throws Exception {
        this.IDsummary = LCMSID.ReadLCMSIDSerialization(Filename,tag);
        if (this.IDsummary == null) {
            return false;
        }
        if (ms1lcms != null) {
            ms1lcms.IDsummary = IDsummary;
        }
        this.IDsummary.Filename=Filename;
        return true;
    }

    TandemParam tparam;

    public void SetTandemParam(TandemParam param) {
        tparam = param;
    }

    public void ReadIDFromDB(String mzxmlname, String FastaFile, ConnectionManager connectionManager) throws SQLException, XmlPullParserException, IOException, ClassNotFoundException, InterruptedException {

        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzxmlname) + "...");
        this.IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(mzxmlname));
        this.IDsummary.FastaPath = FastaFile;
        Connection connection = connectionManager.GetConnection();
        this.IDsummary.ReadFromDBPepIon(connection);

        if (UseMappedIon) {
            this.IDsummary.ReadFromDBMappedPepIon(connection, FilterMappedIonByProb, MappedIonProbThreshold);
            this.IDsummary.ReadFromDBMappedPepFragments(connection);
        }

        this.IDsummary.ReadFromDBPSM(connection);
        this.IDsummary.ReadFromDBPepFragments(connection);
        this.IDsummary.ReadFromDBProt(connection);
        this.IDsummary.AssignProtForPepIon();
        if (UseMappedIon) {
            this.IDsummary.AssignProtForMappedIon();
        }
        IDsummary.GenearteAssignIonList();
        IDsummary.GenerateIndisProtMap();
//        if (ms1lcms != null) {
//            ms1lcms.IDsummary = IDsummary;
//            ms1lcms.AssignMS1Cluster();
//        }
//        if (DIAWindows != null) {
//            AssignMS2Cluster();
//        }
//        AssignFragments();
        connectionManager.CloseConnection();
        Logger.getRootLogger().info("Protein No.:" + IDsummary.ProteinList.size() + "; Assigned Peptide No.:" + IDsummary.AssignedPepIonList.size() + "; All peptide No.:" + IDsummary.GetPepIonList().size());
    }

    public void ParseTPPByRefID(DBSearchParam searchPara, LCMSID RefID, boolean UseRefIDProt) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {

        iProphPepXMLs = new ArrayList<>();
        String PepXMLPath1 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ1Name() + ".pep.xml");
        String PepXMLPath2 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ2Name() + ".pep.xml");
        String PepXMLPath3 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ3Name() + ".pep.xml");
        iProphPepXMLs.add(PepXMLPath1);
        iProphPepXMLs.add(PepXMLPath2);
        iProphPepXMLs.add(PepXMLPath3);

        if (!UseRefIDProt) {
            String CombinedProt = FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_Combine.prot.xml";
            CombineProteinProphetByRefID(searchPara, CombinedProt, RefID);
        } else {
            TPPResult tppresult = new TPPResult(searchPara.PepFDR, searchPara.ProtFDR, searchPara.DecoyPrefix);
            IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
            IDsummary.FastaPath = searchPara.FastaPath;
            tppresult.ReadSearchResultByRefIDUseRefProtID(IDsummary, iProphPepXMLs, RefID, UseMappedIon);
        }
        if (ms1lcms != null) {
            this.ms1lcms.IDsummary = IDsummary;
        }
        //CheckPSM();
    }

    public void ParseCombinePepXML(DBSearchParam searchPara, LCMSID ReferenceID) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.FastaPath = searchPara.FastaPath;
        TPPResult tppresult = new TPPResult(searchPara.PepFDR, searchPara.ProtFDR, searchPara.DecoyPrefix);
        tppresult.ReadSearchResultByRefIDProt(IDsummary, GetCombinePepXML(), ReferenceID);
        Logger.getRootLogger().info("No. of peptide ions:" + IDsummary.ProteinList.size() + "; Assigned Peptide No.:" + IDsummary.AssignedPepIonList.size() + "; All peptide No.:" + IDsummary.GetPepIonList().size());
        if (ms1lcms != null) {
            this.ms1lcms.IDsummary = IDsummary;
        }
    }
    LCMSID templcmsid;

    public void ParsePepXMLtoLCMSID(DBSearchParam searchPara) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {

        iProphPepXMLs.clear();
        String PepXMLPath1 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ1Name() + ".pep.xml");
        String PepXMLPath2 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ2Name() + ".pep.xml");
        String PepXMLPath3 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ3Name() + ".pep.xml");
        iProphPepXMLs.add(PepXMLPath1);
        iProphPepXMLs.add(PepXMLPath2);
        iProphPepXMLs.add(PepXMLPath3);
        templcmsid = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        for (String pepxml : iProphPepXMLs) {
            PepXMLParser pepxmlparser = new PepXMLParser(templcmsid, pepxml, 0f);
            templcmsid.FilterByPepDecoyFDR(searchPara.DecoyPrefix, searchPara.PepFDR);
            Logger.getRootLogger().info("No. of peptide ions:" + templcmsid.GetPepIonList().size() + "; Peptide level threshold: " + templcmsid.PepProbThreshold);

        }
    }

    public void ParsePepXML(DBSearchParam searchPara) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {

        String PepXMLPath1 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ1Name() + ".pep.xml");
        String PepXMLPath2 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ2Name() + ".pep.xml");
        String PepXMLPath3 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ3Name() + ".pep.xml");
        iProphPepXMLs.add(PepXMLPath1);
        iProphPepXMLs.add(PepXMLPath2);
        iProphPepXMLs.add(PepXMLPath3);
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.FastaPath = searchPara.FastaPath;
        for (String pepxml : iProphPepXMLs) {
            LCMSID pepxmlid = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
            PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);
            pepxmlid.FilterByPepDecoyFDR(searchPara.DecoyPrefix, searchPara.PepFDR);
            Logger.getRootLogger().info("No. of peptide ions:" + pepxmlid.GetPepIonList().size() + "; Peptide level threshold: " + pepxmlid.PepProbThreshold);
            for (PepIonID pepID : pepxmlid.GetPepIonList().values()) {
                IDsummary.AddPeptideID(pepID);
            }
        }
        IDsummary.ReMapProPep();

        Logger.getRootLogger().info("Total number of peptide ions:" + IDsummary.GetPepIonList().size());
        if (ms1lcms != null) {
            this.ms1lcms.IDsummary = IDsummary;
        }
    }

    public void LableMatchedFragmentFromMS1Cluster(String mgfname, HashMap<Integer, Integer> ScanClusterMap) {
        for (PepIonID peptide : IDsummary.GetPepIonList().values()) {
            PSM psm = null;
            for (PSM psm2 : peptide.GetPSMList()) {
                if (psm == null || psm2.Probability > psm.Probability) {
                    psm = psm2;
                }
            }
            if (psm != null && psm.GetRawNameString().equals(mgfname)) {
                int ClusterIndex = ScanClusterMap.get(psm.ScanNo);
                for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= psm.ObserPrecursorMz() && DIAWindow.DIA_MZ_Range.getY() >= psm.ObserPrecursorMz()) {
                        if (DIAWindow.FragmentsClu2Cur.containsKey(ClusterIndex)) {
                            ArrayList<PrecursorFragmentPairEdge> fragments = DIAWindow.FragmentsClu2Cur.get(ClusterIndex);
                            ArrayList<PrecursorFragmentPairEdge> newlist = new ArrayList<>();
                            ArrayList<Float> CorrArrayList = new ArrayList<>();
                            for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
                                CorrArrayList.add(fragmentClusterUnit.Correlation);
                            }
                            Collections.sort(CorrArrayList);
                            Collections.reverse(CorrArrayList);

                            for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
                                int CorrRank = 0;
                                for (int intidx = 0; intidx < CorrArrayList.size(); intidx++) {
                                    if (CorrArrayList.get(intidx) <= fragmentClusterUnit.Correlation) {
                                        CorrRank = intidx + 1;
                                        break;
                                    }
                                }
                                if (fragmentClusterUnit.ComplementaryFragment || (fragmentClusterUnit.Correlation >= parameter.CorrThreshold && CorrRank <= parameter.FragmentRank && fragmentClusterUnit.FragmentMS1Rank <= parameter.PrecursorRank && fragmentClusterUnit.ApexDelta <= parameter.ApexDelta)) {
                                    newlist.add(fragmentClusterUnit);
                                }
                            }

                            HashSet<Integer> InlcudeIndex = new HashSet<>();
                            for (PrecursorFragmentPairEdge fragmentClusterUnit : newlist) {
                                if (!InlcudeIndex.contains(fragmentClusterUnit.PeakCurveIndexB)) {
                                    InlcudeIndex.add(fragmentClusterUnit.PeakCurveIndexB);
                                    for (FragmentPeak frag : peptide.FragmentPeaks) {
                                        if (InstrumentParameter.CalcPPM(fragmentClusterUnit.FragmentMz, frag.FragMZ) < parameter.MS2PPM) {
                                            if (!DIAWindow.MatchedFragmentMap.containsKey(fragmentClusterUnit.PeakCurveIndexB)) {
                                                ArrayList<PrecursorFragmentPairEdge> ClusterList = new ArrayList<>();
                                                DIAWindow.MatchedFragmentMap.put(fragmentClusterUnit.PeakCurveIndexB, ClusterList);
                                            }
                                            DIAWindow.MatchedFragmentMap.get(fragmentClusterUnit.PeakCurveIndexB).add(fragmentClusterUnit);
////                                        if (InstrumentParameter.CalcPPM(fragmentClusterUnit.FragmentMz, frag.FragMZ) < parameter.MS2PPM) {
////                                            if (DIAWindow.WindowID.equals("1074_1100") && fragmentClusterUnit.PeakCurveIndexB == 45063) {
////                                                System.out.println(frag.IonType);
////                                                System.out.println(frag.FragMZ);
////                                                System.out.println(frag.Charge);
////                                            }
////                                        }
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void LableMatchedFragmentFromUnfragMS1Cluster(String mgfname, HashMap<Integer, String> ScanClusterMap) {
        for (PepIonID peptide : IDsummary.GetPepIonList().values()) {
            PSM psm = null;
            for (PSM psm2 : peptide.GetPSMList()) {
                if (psm == null || psm2.Probability > psm.Probability) {
                    psm = psm2;
                }
            }
            if (psm != null && psm.GetRawNameString().equals(mgfname)) {
                int ClusterIndex = Integer.parseInt(ScanClusterMap.get(psm.ScanNo));
                for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= psm.ObserPrecursorMz() && DIAWindow.DIA_MZ_Range.getY() >= psm.ObserPrecursorMz()) {
                        if (DIAWindow.UnFragIonClu2Cur.containsKey(ClusterIndex)) {
                            for (PrecursorFragmentPairEdge fragmentClusterUnit : DIAWindow.UnFragIonClu2Cur.get(ClusterIndex)) {
                                for (FragmentPeak frag : peptide.FragmentPeaks) {
                                    if (InstrumentParameter.CalcPPM(fragmentClusterUnit.FragmentMz, frag.FragMZ) < parameter.MS2PPM) {
                                        if (!DIAWindow.MatchedFragmentMap.containsKey(fragmentClusterUnit.PeakCurveIndexB)) {
                                            ArrayList<PrecursorFragmentPairEdge> ClusterList = new ArrayList<>();
                                            DIAWindow.MatchedFragmentMap.put(fragmentClusterUnit.PeakCurveIndexB, ClusterList);
                                        }
                                        DIAWindow.MatchedFragmentMap.get(fragmentClusterUnit.PeakCurveIndexB).add(fragmentClusterUnit);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void ReportSharedFrag() throws InterruptedException, IOException, ParserConfigurationException, SAXException, FileNotFoundException, ExecutionException, XmlPullParserException, DataFormatException, ClassNotFoundException {
        Logger.getRootLogger().info("Removing matched fragments for 2nd stage searching.......");
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            DIAwindow.ReadIfProcessed();
        }
        LableMatchedFragmentFromMS1Cluster(GetQ1Name(), ScanClusterMap_Q1);
        LableMatchedFragmentFromMS1Cluster(GetQ2Name(), ScanClusterMap_Q2);
        //LableMatchedFragmentFromUnfragMS1Cluster(GetQ3Name(), ScanClusterMap_Q3);
        //ms1lcms.ReadPeakCluster();
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_SharedFragment.xls");
        writer.write("ID\tWindow\tFragIndex\tPrecursorIndex\tmz\tCorr\tApexDelta\tintensity\tPrecursorMz\tPrecursorCharge\n");
        for (LCMSPeakDIAMS2 DIAWindow : DIAWindows) {
            for (ArrayList<PrecursorFragmentPairEdge> fragmments : DIAWindow.MatchedFragmentMap.values()) {
                for (PrecursorFragmentPairEdge frag : fragmments) {
                    //PeakCluster cluster=ms1lcms.PeakClusters.get(frag.PeakCurveIndexA-1);
                    writer.write(DIAWindow.WindowID + "_" + frag.PeakCurveIndexB + "\t" + DIAWindow.WindowID + "\t" + frag.PeakCurveIndexB + "\t" + frag.PeakCurveIndexA + "\t" + frag.FragmentMz + "\t" + frag.Correlation + "\t" + frag.ApexDelta + "\t" + frag.Intensity + "\n");
                    //writer.write(DIAWindow.WindowID+"_"+ frag.PeakCurveIndexB+"\t"+DIAWindow.WindowID+"\t"+frag.PeakCurveIndexB+"\t"+frag.PeakCurveIndexA+"\t"+frag.FragmentMz+"\t"+frag.Correlation+"\t"+frag.ApexDelta+"\t"+frag.Intensity+"\t"+cluster.TargetMz()+"\t"+cluster.Charge+"\n");
                }
            }
        }
        writer.close();
    }

    public void RemoveMatchedFragment() throws InterruptedException, IOException, ParserConfigurationException, SAXException, FileNotFoundException, ExecutionException, XmlPullParserException, DataFormatException, ClassNotFoundException {
        Logger.getRootLogger().info("Removing matched fragments for 2nd stage searching.......");
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            DIAwindow.ReadIfProcessed();
        }
        LableMatchedFragmentFromMS1Cluster(GetQ1Name(), ScanClusterMap_Q1);
        LableMatchedFragmentFromMS1Cluster(GetQ2Name(), ScanClusterMap_Q2);
        LableMatchedFragmentFromUnfragMS1Cluster(GetQ3Name(), ScanClusterMap_Q3);
        GenerateMGF_2ndStage();
    }

    public void GenerateRawMGF() throws IOException, Exception {
        HashMap<Integer, ArrayList<PseudoMSMSProcessing>> ScanList = new HashMap<>();
        HashMap<String, PseudoMSMSProcessing> UnfragScanList = new HashMap<>();
        parameter.BoostComplementaryIon = false;
        ExecutorService executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            for (PeakCluster ms1cluster : ms1lcms.PeakClusters) {
                if (DIAwindow.DIA_MZ_Range.getX() <= ms1cluster.TargetMz() && DIAwindow.DIA_MZ_Range.getY() >= ms1cluster.TargetMz() && DIAwindow.FragmentsClu2Cur.containsKey(ms1cluster.Index)) {
                    PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms1cluster, DIAwindow.FragmentsClu2Cur.get(ms1cluster.Index), parameter);
                    executorPool.execute(mSMSProcessing);
                    if (!ScanList.containsKey(ms1cluster.Index)) {
                        ScanList.put(ms1cluster.Index, new ArrayList<PseudoMSMSProcessing>());
                    }
                    ScanList.get(ms1cluster.Index).add(mSMSProcessing);
                }
            }
            for (PeakCluster ms1cluster : DIAwindow.PeakClusters) {
                if (DIAwindow.DIA_MZ_Range.getX() <= ms1cluster.TargetMz() && DIAwindow.DIA_MZ_Range.getY() >= ms1cluster.TargetMz() && DIAwindow.UnFragIonClu2Cur.containsKey(ms1cluster.Index)) {
                    PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms1cluster, DIAwindow.UnFragIonClu2Cur.get(ms1cluster.Index), parameter);
                    executorPool.execute(mSMSProcessing);
                    UnfragScanList.put(DIAwindow.WindowID + ";" + ms1cluster.Index, mSMSProcessing);
                }
            }
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }

        ReadScanNoMapping();
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".Raw.mgf";
        FileWriter mgfWriter = new FileWriter(mgffile, false);

        for (int ScanNo = 0; ScanNo < ScanClusterMap_Q1.size(); ScanNo++) {
            int ClusterIndex = ScanClusterMap_Q1.get(ScanNo);
            XYPointCollection Scan = new XYPointCollection();
            PseudoMSMSProcessing mSMSProcessing = null;
            for (PseudoMSMSProcessing MS2Processing : ScanList.get(ClusterIndex)) {
                mSMSProcessing = MS2Processing;
                for (PrecursorFragmentPairEdge fragmentClusterUnit : MS2Processing.fragments) {
                    Scan.AddPointKeepMaxIfValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity);
                }
            }
            StringBuilder mgfString = new StringBuilder();
            mgfString.append("BEGIN IONS\n");
            mgfString.append("PEPMASS=" + mSMSProcessing.ms1cluster.TargetMz() + "\n");
            mgfString.append("CHARGE=" + mSMSProcessing.ms1cluster.Charge + "+\n");
            mgfString.append("RTINSECONDS=" + mSMSProcessing.ms1cluster.PeakHeightRT[0] * 60f + "\n");
            mgfString.append("TITLE=ClusterIndex:" + mSMSProcessing.ms1cluster.Index + "\n");
            for (int i = 0; i < Scan.PointCount(); i++) {
                mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
            }
            mgfString.append("END IONS\n\n");
            mgfWriter.write(mgfString.toString());
        }
        mgfWriter.close();

        ////////////////////////////////////////////////////////////////////////////////
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".Raw.mgf";
        FileWriter mgfWriter2 = new FileWriter(mgffile2, false);

        for (int ScanNo = 0; ScanNo < ScanClusterMap_Q2.size(); ScanNo++) {
            int ClusterIndex = ScanClusterMap_Q2.get(ScanNo);
            XYPointCollection Scan = new XYPointCollection();
            PseudoMSMSProcessing mSMSProcessing = null;
            for (PseudoMSMSProcessing MS2Processing : ScanList.get(ClusterIndex)) {
                mSMSProcessing = MS2Processing;
                for (PrecursorFragmentPairEdge fragmentClusterUnit : MS2Processing.fragments) {
                    Scan.AddPointKeepMaxIfValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity);
                }
            }
            StringBuilder mgfString = new StringBuilder();
            mgfString.append("BEGIN IONS\n");
            mgfString.append("PEPMASS=" + mSMSProcessing.ms1cluster.TargetMz() + "\n");
            mgfString.append("CHARGE=" + mSMSProcessing.ms1cluster.Charge + "+\n");
            mgfString.append("RTINSECONDS=" + mSMSProcessing.ms1cluster.PeakHeightRT[0] * 60f + "\n");
            mgfString.append("TITLE=ClusterIndex:" + mSMSProcessing.ms1cluster.Index + "\n");
            for (int i = 0; i < Scan.PointCount(); i++) {
                mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
            }
            mgfString.append("END IONS\n\n");
            mgfWriter2.write(mgfString.toString());
        }

        mgfWriter2.close();

        ////////////////////////////////
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".Raw.mgf";
        FileWriter mgfWriter3 = new FileWriter(mgffile3, false);
        mzXMLParser Q3mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        Q3mzxml.GetAllScanCollectionMS2Only(true, false);
        for (int ScanNo = 0; ScanNo < ScanClusterMap_Q3.size(); ScanNo++) {
            String key = ScanClusterMap_Q3.get(ScanNo);
            XYPointCollection Scan = new XYPointCollection();
            PseudoMSMSProcessing mSMSProcessing = UnfragScanList.get(key);

            for (PrecursorFragmentPairEdge fragmentClusterUnit : mSMSProcessing.fragments) {
                Scan.AddPointKeepMaxIfValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity);
            }

            StringBuilder mgfString = new StringBuilder();
            mgfString.append("BEGIN IONS\n");
            mgfString.append("PEPMASS=" + mSMSProcessing.ms1cluster.TargetMz() + "\n");
            mgfString.append("CHARGE=" + mSMSProcessing.ms1cluster.Charge + "+\n");
            mgfString.append("RTINSECONDS=" + mSMSProcessing.ms1cluster.PeakHeightRT[0] * 60f + "\n");
            mgfString.append("TITLE=ClusterIndex:" + mSMSProcessing.ms1cluster.Index + "\n");
            for (int i = 0; i < Scan.PointCount(); i++) {
                mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
            }
            mgfString.append("END IONS\n\n");
            mgfWriter3.write(mgfString.toString());
        }
        mgfWriter3.close();
    }

    public void GenerateMGF_2ndStage() throws IOException, InterruptedException {
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            DIAwindow.GenerateMGF(ms1lcms);
        }
        RenameMGF("_RC");
    }

    public void ParseSearchEngineResult(MSMSDBSearch dbsearch) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        ParseSearchEngineResult(dbsearch, "");
    }

    public void ParseSearchEngineResult(MSMSDBSearch dbsearch, String mgftag) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        String mgfname1 = FilenameUtils.getFullPath(Filename) + GetQ1Name() + mgftag + ".mgf";
        String mgfname2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + mgftag + ".mgf";
        String mgfname3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + mgftag + ".mgf";
        dbsearch.SetResultFilePath(mgfname1);
        iProphPepXMLs.add(dbsearch.GetParameter().InteractPepXMLPath);
        dbsearch.SetResultFilePath(mgfname2);
        iProphPepXMLs.add(dbsearch.GetParameter().InteractPepXMLPath);
        dbsearch.SetResultFilePath(mgfname3);
        iProphPepXMLs.add(dbsearch.GetParameter().InteractPepXMLPath);
        TPPResult tppresult = new TPPResult(dbsearch.GetParameter().PepFDR, dbsearch.GetParameter().ProtFDR, dbsearch.GetParameter().DecoyPrefix);
        dbsearch.SetCombineFileName(Filename);
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.FastaPath = dbsearch.GetParameter().FastaPath;
        tppresult.ReadSearchResult(IDsummary, iProphPepXMLs, dbsearch.GetParameter().CombinedProt);
        if (ms1lcms != null) {
            this.ms1lcms.IDsummary = IDsummary;
        }
    }

    public void ParseSearchEngineResultiProphet(DBSearchParam param) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException, Exception {
        iProphPepXMLs.add(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".iproph.pep.xml");
        iProphPepXMLs.add(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".iproph.pep.xml");
        iProphPepXMLs.add(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".iproph.pep.xml");
        TPPResult tppresult = new TPPResult(param.PepFDR, param.ProtFDR, param.DecoyPrefix);
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.FastaPath = param.FastaPath;
        tppresult.ReadSearchResult(IDsummary, iProphPepXMLs, GetCombineiProphProtXML());
//        mzXMLParser mgfname1 = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
//        mzXMLParser mgfname2 =  new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
//        mzXMLParser mgfname3 = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
//        HashMap<String, mzXMLParser> mgfs=new HashMap<>();
//        mgfs.put(GetQ1Name(), mgfname1);
//        mgfs.put(GetQ2Name(), mgfname2);
//        mgfs.put(GetQ3Name(), mgfname3);
//        for(PSM psm : IDsummary.PSMList.values()){
//            if(psm.RetentionTime==-1f){
//                psm.RetentionTime=mgfs.get(psm.GetRawNameString()).GetScanElutionTimeMap().get(psm.ScanNo);
//                psm.NeighborMaxRetentionTime=psm.RetentionTime;
//            }
//        }
        if (ms1lcms != null) {
            this.ms1lcms.IDsummary = IDsummary;
        }
    }

    public void ParseTPP(DBSearchParam searchPara) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {

        String PepXMLPath1 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ1Name() + ".pep.xml");
        String PepXMLPath2 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ2Name() + ".pep.xml");
        String PepXMLPath3 = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(Filename) + "interact-" + GetQ3Name() + ".pep.xml");
        iProphPepXMLs.add(PepXMLPath1);
        iProphPepXMLs.add(PepXMLPath2);
        iProphPepXMLs.add(PepXMLPath3);
        TPPResult tppresult = new TPPResult(searchPara.PepFDR, searchPara.ProtFDR, searchPara.DecoyPrefix);

        String CombinedProt = FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_combine.prot.xml";
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.FastaPath = searchPara.FastaPath;
        tppresult.ReadSearchResult(IDsummary, iProphPepXMLs, CombinedProt);
        if (ms1lcms != null) {
            this.ms1lcms.IDsummary = IDsummary;
        }
    }

    public void CombineProteinProphetByRefID(DBSearchParam searchpara, String CombinedProt, LCMSID RefID) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, ClassNotFoundException {
        ProteinProphetCombine protprophet = new ProteinProphetCombine();
        if (!new File(CombinedProt).exists()) {
            protprophet.ProteinProphetCombineSearch(searchpara, iProphPepXMLs, CombinedProt);
        }
        TPPResult tppresult = new TPPResult(searchpara.PepFDR, searchpara.ProtFDR, searchpara.DecoyPrefix);
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.DecoyTag = "rev_";
        IDsummary.FastaPath = searchpara.FastaPath;
        tppresult.ReadSearchResultByRefID(IDsummary, iProphPepXMLs, CombinedProt, RefID);
    }

    public void EstimateNumberOfCandidatesTargeted(HashMap<String, PepIonID> peplist) throws InterruptedException, ExecutionException, IOException, FileNotFoundException, ParserConfigurationException, SAXException, DataFormatException, ClassNotFoundException, XmlPullParserException, Exception {
        BuildDIAWindows();
        mzXMLParser T1mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        mzXMLParser T2mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        mzXMLParser T4mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);

        SortedPepCandidate IonMzLib = new SortedPepCandidate();

        for (PepIonID pep : peplist.values()) {
            PepIonCandidate candidate = new PepIonCandidate();
            candidate.peptide = pep.GetPepFactory();
            candidate.Charge = pep.Charge;
            candidate.SetMz(pep.NeutralPrecursorMz());
            candidate.RT = pep.GetIDRT();
            IonMzLib.add(candidate);
        }

        int count = 0;
        int total = 0;
        T1mzxml.GetAllScanCollectionMS2Only(true, true);
        ScanCollection Scans = T1mzxml.scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            for (LCMSPeakDIAMS2 swath : DIAWindows) {
                if (swath.DIA_MZ_Range.getX() <= scan.PrecursorMz && swath.DIA_MZ_Range.getY() >= scan.PrecursorMz) {
                    float lowmz = swath.DIA_MZ_Range.getX();
                    float highmz = swath.DIA_MZ_Range.getY();
                    ArrayList<PepIonCandidate> candidates = IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
                    for (PepIonCandidate candidate : candidates) {
                        if (Math.abs(candidate.RT - scan.RetentionTime) < 1f) {
                            count++;
                        }
                    }
                }
            }
        }
        total += Scans.ScanHashMap.size();
        T2mzxml.GetAllScanCollectionMS2Only(true, true);
        Scans = T2mzxml.scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            for (LCMSPeakDIAMS2 swath : DIAWindows) {
                if (swath.DIA_MZ_Range.getX() <= scan.PrecursorMz && swath.DIA_MZ_Range.getY() >= scan.PrecursorMz) {
                    float lowmz = swath.DIA_MZ_Range.getX();
                    float highmz = swath.DIA_MZ_Range.getY();
                    ArrayList<PepIonCandidate> candidates = IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
                    for (PepIonCandidate candidate : candidates) {
                        if (Math.abs(candidate.RT - scan.RetentionTime) < 1f) {
                            count++;
                        }
                    }
                }
            }
        }
        total += Scans.ScanHashMap.size();
        T4mzxml.GetAllScanCollectionMS2Only(true, true);
        Scans = T4mzxml.scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            for (LCMSPeakDIAMS2 swath : DIAWindows) {
                if (swath.DIA_MZ_Range.getX() <= scan.PrecursorMz && swath.DIA_MZ_Range.getY() >= scan.PrecursorMz) {
                    float lowmz = swath.DIA_MZ_Range.getX();
                    float highmz = swath.DIA_MZ_Range.getY();
                    ArrayList<PepIonCandidate> candidates = IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
                    for (PepIonCandidate candidate : candidates) {
                        if (Math.abs(candidate.RT - scan.RetentionTime) < 1f) {
                            count++;
                        }
                    }
                }
            }
        }
        total += Scans.ScanHashMap.size();
        Logger.getRootLogger().info(Filename + ": Average number of candidates (targeted):" + ((double) count / total));
    }

    public void EstimateNumberOfCandidates(String FastaFile) throws InterruptedException, ExecutionException, IOException, FileNotFoundException, ParserConfigurationException, SAXException, DataFormatException, ClassNotFoundException, XmlPullParserException, Exception {

        PepIonLib lib = new PepIonLib(FastaFile);
        lib.GetWholeIonLib();
        Logger.getRootLogger().info("Ion library size:" + lib.IonMzLib.size() + "\n");
        mzXMLParser T1mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        mzXMLParser T2mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        mzXMLParser T4mzxml = new mzXMLParser(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mzXML", parameter, SpectralDataType.DataType.DDA, null, NoCPUs);

        int count = 0;
        int total = 0;
        T1mzxml.GetAllScanCollectionMS2Only(true, true);
        ScanCollection Scans = T1mzxml.scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            float lowmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, parameter.MS1PPM);
            float highmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, -parameter.MS1PPM);
            ArrayList<PepIonCandidate> candidates = lib.IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
            count += candidates.size();
        }
        total += Scans.ScanHashMap.size();
        T2mzxml.GetAllScanCollectionMS2Only(true, true);
        Scans = T2mzxml.scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            float lowmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, parameter.MS1PPM);
            float highmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, -parameter.MS1PPM);
            ArrayList<PepIonCandidate> candidates = lib.IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
            count += candidates.size();
        }
        total += Scans.ScanHashMap.size();
        T4mzxml.GetAllScanCollectionMS2Only(true, true);
        Scans = T4mzxml.scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            float lowmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, parameter.MS1PPM);
            float highmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, -parameter.MS1PPM);
            ArrayList<PepIonCandidate> candidates = lib.IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
            count += candidates.size();
        }
        total += Scans.ScanHashMap.size();

        Logger.getRootLogger().info("\n" + Filename + ": Average number of candidates:" + (double) count / total + "\n");
    }

    private void CalculateOverlap(LCMSID summary1, LCMSID summary2) {

        int proteinoverlap = 0;
        int peptideoverlap = 0;
        for (ProtID protein1 : summary1.ProteinList.values()) {
            for (ProtID protein2 : summary2.ProteinList.values()) {
                if (protein1.getAccNo() == null ? protein2.getAccNo() == null : protein1.getAccNo().equals(protein2.getAccNo())) {
                    proteinoverlap++;
                }
            }
        }

        for (PepIonID peptide1 : summary1.AssignedPepIonList.values()) {
            for (PepIonID peptide2 : summary2.AssignedPepIonList.values()) {
                if (peptide1.GetKey() == null ? peptide2.GetKey() == null : peptide1.GetKey().equals(peptide2.GetKey())) {
                    peptideoverlap++;
                }
            }
        }
        Logger.getRootLogger().info("Overlap between " + FilenameUtils.getBaseName(summary1.mzXMLFileName) + " and " + FilenameUtils.getBaseName(summary2.mzXMLFileName) + " - protein:" + proteinoverlap + " peptide: " + peptideoverlap + "\n");
    }

    private void MergeIDSummary(LCMSID summary) throws ClassNotFoundException, InterruptedException, XmlPullParserException, IOException {
        for (String key : summary.ProteinList.keySet()) {
            ProtID protein = summary.ProteinList.get(key);
            if (!this.IDsummary.ProteinList.containsKey(key)) {
                this.IDsummary.AddProtID(protein);
            } else {
                ProtID existprotein = this.IDsummary.ProteinList.get(key);
                for (String pepkey : protein.PeptideID.keySet()) {
                    PepIonID pepIonID = protein.PeptideID.get(pepkey);
                    if (!existprotein.PeptideID.containsKey(pepkey)) {
                        existprotein.AddPeptideID(pepIonID);
                    } else {
                        PepIonID existpep = existprotein.PeptideID.get(pepkey);
                        for (PSM psm : pepIonID.GetPSMList()) {
                            existpep.AddPSM(psm);
                        }
                    }
                }
            }
        }
        IDsummary.GenearteAssignIonList();
    }

    public void ReadIonLib(DBSearchParam tandemsearch) throws IOException, FileNotFoundException, ClassNotFoundException, InterruptedException, XmlPullParserException {
        IonLib = new PepIonLib(tandemsearch.FastaPath);
    }

    public void UmpireSearch(DBSearchParam searchPara) throws SQLException, IOException, XmlPullParserException, FileNotFoundException, ClassNotFoundException, InterruptedException {
        Logger.getRootLogger().info("Loading all peaks.....");
        ReadIonLib(searchPara);
        ReadFactorialTable();
        new File(FilenameUtils.getFullPath(ms1lcms.ScanCollectionName) + FilenameUtils.getBaseName(ms1lcms.ScanCollectionName) + "Search.txt").delete();

        //ms1lcms.GenerateIsolatedPeakCurve();        
        DIAWindows = new ArrayList<>();
        int count = 1;
        Logger.getRootLogger().info("Search peptides for SWATH:");
        for (XYData swathmz : GetMzXML().dIA_Setting.DIAWindows.keySet()) {
            Logger.getRootLogger().info("(" + count + "/" + GetMzXML().dIA_Setting.DIAWindows.size() + ")...");
            LCMSPeakDIAMS2 swath = new LCMSPeakDIAMS2(Filename, this, parameter, swathmz, GetMzXML(), NoCPUs);
            swath.SetMySQLConnection(connectionManager);
            //swath.ReadPeakCluster();
            swath.ReadPrecursorFragmentClu2Cur();
            //swath.ReadUnFragCurveClusterResultFromDB();
            swath.UmpireSearch(ms1lcms, FactorialTable, IonLib, searchPara);
            //SwathWindows.add(dia);
            count++;
        }
    }

    public void ClearStructure(){
        ms1lcms=null;
        DIAWindows=null;
        dIA_Setting=null;
        mzXML=null;
    }
            
    public void BuildStructure() throws SQLException, FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        //System.out.println("building peak clusters.....");
        LoadDIASetting();
        ms1lcms = new LCMSPeakMS1(Filename, NoCPUs);
        ms1lcms.datattype = dIA_Setting.dataType;
        ms1lcms.SetParameter(parameter);
        ms1lcms.SetMySQLConnection(connectionManager);

        if (IDsummary != null) {
            ms1lcms.IDsummary = IDsummary;
        }

        DIAWindows = new ArrayList<>();
        //System.out.print("Querying SWATH:");
        if (dIA_Setting.DIAWindows == null || dIA_Setting.DIAWindows.isEmpty()) {
            GetMzXML();
        }
        for (XYData swathmz : dIA_Setting.DIAWindows.keySet()) {
            LCMSPeakDIAMS2 dia = new LCMSPeakDIAMS2(Filename, this, parameter, swathmz, GetMzXML(), NoCPUs);
            dia.SetMySQLConnection(connectionManager);
            dia.datattype = dIA_Setting.dataType;
            DIAWindows.add(dia);
        }
    }

    public void ReadFromDBMS1Only(boolean PeakCurveIncluded) throws SQLException {
        Logger.getRootLogger().info("Loading peak data from MySQL..");
        ms1lcms = new LCMSPeakMS1(Filename, NoCPUs);
        ms1lcms.SetParameter(parameter);
        ms1lcms.SetMySQLConnection(connectionManager);
        if (PeakCurveIncluded) {
            ms1lcms.ReadPeakCurveFromDB(false);
        }
        ms1lcms.ReadPeakClusterFromDB();
        DIAWindows = new ArrayList<>();
        int count = 1;
        Logger.getRootLogger().info("Creating DIA windows:");
        for (XYData swathmz : GetMzXML().dIA_Setting.DIAWindows.keySet()) {
            //System.out.print("(" + count + "/" + GetMzXML().SWATHWindows.size() + ")...");
            LCMSPeakDIAMS2 dia = new LCMSPeakDIAMS2(Filename, this, parameter, swathmz, GetMzXML(), NoCPUs);
            dia.SetMySQLConnection(connectionManager);
            DIAWindows.add(dia);
            count++;
        }
    }

    private float GetGrowthFactor(PeakCluster MS1PeakCluster, XYData SWATHwindow) {
        float maxintensity = 0f;
        for (PeakCluster peakCluster : ms1lcms.PeakClusters) {
            if (peakCluster.TargetMz() >= SWATHwindow.getX() && peakCluster.TargetMz() <= SWATHwindow.getY() && Math.abs(peakCluster.PeakHeightRT[0] - MS1PeakCluster.PeakHeightRT[0]) < 0.5f) {

                if (maxintensity < peakCluster.PeakHeight[0]) {
                    maxintensity = peakCluster.PeakHeight[0];
                }
            }
        }
        float growth = maxintensity / MS1PeakCluster.PeakHeight[0];
        if (growth < 1) {
            Logger.getRootLogger().info("growth<1:" + Thread.currentThread().getStackTrace()[2].getLineNumber());
        }
        return growth;
    }

    public boolean MGFgenerated() {
        if (new File(FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf").exists()) {
            return true;
        }
        return false;
    }

    private void MS1PeakDetection() throws SQLException, InterruptedException, ExecutionException, IOException, ParserConfigurationException, SAXException, FileNotFoundException, Exception {
        RemoveMGF();
        ms1lcms = new LCMSPeakMS1(Filename, NoCPUs);
        ms1lcms.datattype = dIA_Setting.dataType;
        ms1lcms.SetParameter(parameter);
        ms1lcms.SetMySQLConnection(connectionManager);
        ms1lcms.SetMS1Windows(dIA_Setting.MS1Windows);
        ms1lcms.CreatePeakFolder();
        ms1lcms.ExportFragmentPeak = ExportPrecursorPeak;
        //ms1lcms.SwathWindows = DIAWindows;
        ms1lcms.SetmzXML(GetMzXML());
        Logger.getRootLogger().info("Processing MS1 peak detection");
        ms1lcms.AssignIDResult(DDAIDsummary);
        ms1lcms.ExportPeakClusterTable = ExportPeakClusterTable;
        ms1lcms.PeakClusterDetection();

        Logger.getRootLogger().info("==================================================================================");
    }

    private void RemoveMGF() {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf";
        File file = new File(mgffile);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.delete();
        }
        mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf.temp";
        mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf.temp";
        mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        file = new File(mgffile);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.delete();
        }
    }

    public void ExportID() throws SQLException, IOException {
        ExportID("");
    }
    public void ExportID(String tag) throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        IDsummary.WriteLCMSIDSerialization(Filename, tag);

        if (connectionManager != null) {
            IDsummary.ExportPepID(connectionManager);
            IDsummary.ExportProtID(connectionManager);
            IDsummary.ExportPepFragmentPeak(connectionManager);
            IDsummary.ExportMappedPepID(connectionManager);
            IDsummary.ExportMappedPepFragmentPeak(connectionManager);
        }
    }

    public void ExportMappedIDQuant() throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(Filename))) {
            return;
        }
        IDsummary.ExportMappedPepID(connectionManager);
        IDsummary.ExportMappedPepFragmentPeak(connectionManager);
    }

    public void DIAMS2PeakDetection() throws SQLException, IOException, InterruptedException, ExecutionException, FileNotFoundException, Exception {
        int count = 1;
        //CreateSWATHTables();
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            Logger.getRootLogger().info("Processing DIA MS2 (mz range):" + DIAwindow.DIA_MZ_Range.getX() + "_" + DIAwindow.DIA_MZ_Range.getY() + "( " + (count++) + "/" + GetMzXML().dIA_Setting.DIAWindows.size() + " )");
            DIAwindow.ExportFragmentPeak = ExportFragmentPeak;
            DIAwindow.ExportPeakClusterTable = ExportPeakClusterTable;
            DIAwindow.PeakDetectionPFGrouping(ms1lcms);
            DIAwindow.ClearAllPeaks();
            Logger.getRootLogger().info("==================================================================================");
        }
        RenameMGF("");
        //}
    }

    private void RenameMGF(String tag) {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf.temp";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf.temp";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        File file = new File(mgffile);
        file = new File(mgffile);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
    }

    public void DiscardMappedIonTable() throws SQLException {
        if (connectionManager != null) {
            Connection connection = connectionManager.GetConnection();
//         this.IDsummary.DiscardMappedPepIonTable(connection);
            try (Statement state = connection.createStatement()) {
                state.execute("DROP TABLE " + FilenameUtils.getBaseName(Filename) + "_MappedPepIonIDs");

            } catch (SQLException ex) {
                Logger.getRootLogger().error("(Discarding table: " + FilenameUtils.getBaseName(Filename) + "_MappedPepIonIDs failed.)\n");
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
            connectionManager.CloseConnection();
        }
    }

    public void GenerateMassCalibrationRTMap() {
        try {
            ms1lcms.GenerateMassCalibrationRTMap();
        } catch (IOException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            DIAwindow.Masscalibrationfunction = ms1lcms.Masscalibrationfunction;
        }
    }

    public void ReplaceProtByRefIDByTheoPep(LCMSID protID) {
        IDsummary.GenerateProteinByRefIDByTheoPep(protID, UseMappedIon);
        //IDsummary.GenearteAssignIonList();        
    }

    public void SaveParams() {
        parameter.WriteParamSerialization(Filename);
    }

    public void SaveDIASetting() {
        dIA_Setting.WriteParamSerialization(Filename);
    }

    public boolean LoadParams() {
        parameter = InstrumentParameter.ReadParametersSerialization(Filename);
        return parameter != null;
    }

    public boolean LoadDIASetting() {
        dIA_Setting = DIA_Setting.ReadDIASettingSerialization(Filename);
        return dIA_Setting != null;
    }

    public InstrumentParameter GetParameter() {
        return parameter;
    }

    public void ParseiProphet(DBSearchParam param) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        IDsummary = new LCMSID(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename));
        IDsummary.FastaPath = param.FastaPath;
        TPPResult tppresult = new TPPResult(param.PepFDR, param.ProtFDR, param.DecoyPrefix);
        tppresult.ReadSearchResult(IDsummary, GetIPROPHETPepXML(), GetIPROPHETProtXML());
    }

}
