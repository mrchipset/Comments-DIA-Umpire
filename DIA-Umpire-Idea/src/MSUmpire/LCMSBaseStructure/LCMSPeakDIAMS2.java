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
package MSUmpire.LCMSBaseStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.DIA.CorrCalcCluster2ClusterUnit;
import MSUmpire.DIA.DIAPack;
import MSUmpire.DIA.IsolatedCurveCorrCalcUnit;
import MSUmpire.DIA.MS1PeakGroupSearch;
import MSUmpire.DIA.PseudoMSMSProcessing;
import MSUmpire.DIA.UmpireMatchUnit;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.PeptidePeakClusterDetection.PDHandlerDIAMS2;
import MSUmpire.UmpireSearchDataStructure.PepIonCandidate;
import MSUmpire.UmpireSearchDataStructure.PepIonLib;
import MSUmpire.spectrumparser.mzXMLParser;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakDIAMS2 extends LCMSPeakBase {

    //public HashMap<String, PeakCluster> UnIDpeakClusters;    
    public XYData DIA_MZ_Range;
    public String WindowID;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> FragmentsClu2Clu;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> FragmentsClu2Cur;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> UnFragIonClu2Clu;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> UnFragIonClu2Cur;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> IsolatedMS1FragmentClustering;
    //public HashMap<Integer, ArrayList<PeakOverlapRegion>> FragmentMSDistance;
    public HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> MatchedFragmentMap = new HashMap<>();
    private DIAPack parentDIA;

    public boolean FragmentGroupByCluster = false;

    public LCMSPeakDIAMS2(String Filename, DIAPack parentDIA, InstrumentParameter parameter, XYData SwathMZ, mzXMLParser mzxml, int NoCPUs) {
        this.DIA_MZ_Range = SwathMZ;
        this.WindowID = (int) Math.floor(SwathMZ.getX()) + "_" + (int) Math.floor(SwathMZ.getY());
        this.mzxml = mzxml;
        this.ScanCollectionName = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + (int) Math.floor(SwathMZ.getX()) + "_" + (int) Math.floor(SwathMZ.getY());
        this.ParentmzXMLName = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename);
        this.parentDIA=parentDIA;
        this.parameter = parameter;
        this.MaxNoPeakCluster = parameter.MaxMS2NoPeakCluster;
        //this.MaxNoPeakCluster = parameter.MaxNoPeakCluster;
        this.MinNoPeakCluster = parameter.MinMS2NoPeakCluster;
        //this.MinNoPeakCluster = parameter.MinNoPeakCluster;
        this.StartCharge = parameter.MS2StartCharge;
        this.EndCharge = parameter.MS2EndCharge;
        //this.StartCharge = parameter.StartCharge;
        //this.EndCharge = parameter.EndCharge;
        this.MiniIntensity = parameter.MinMSMSIntensity;
        this.SNR = parameter.MS2SNThreshold;
        this.NoCPUs = NoCPUs;
    }

    @Override
    public void ClearAllPeaks() {
        //System.out.print("Releasing peak data..(Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)\n");
        PeakClusters = null;
        MZSortedClusters = null;
        ApexRTSortedClusters = null;
        //KeepUnfragClusterOnly();
        PeakCurveListMZ = null;
        PeakCurveListRT = null;
        IsolatedPeakCurves = null;
        FragmentsClu2Cur = null;
        UnFragIonClu2Cur = null;
        System.gc();
        //System.out.print("Peak data is released (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)\n");
    }

    private void KeepUnfragClusterOnly() {
        ArrayList<PeakCluster> removelist = new ArrayList<>();
        for (PeakCluster peakCluster : PeakClusters) {
            if (peakCluster.Charge < parameter.StartCharge || peakCluster.Charge > parameter.EndCharge || peakCluster.TargetMz() < DIA_MZ_Range.getX() && peakCluster.TargetMz() > DIA_MZ_Range.getY()) {
                removelist.add(peakCluster);
            } else {
                peakCluster.IsoPeaksCurves = null;
            }
        }
        for (PeakCluster peakCluster : removelist) {
            PeakClusters.remove(peakCluster);
        }
    }

    public void PeakDetectionPFGrouping(LCMSPeakMS1 ms1lcms) throws SQLException, InterruptedException, ExecutionException, IOException, FileNotFoundException, Exception {
        if (!(Resume && ReadIfProcessed())) {
            PDHandlerDIAMS2 swathdetection = new PDHandlerDIAMS2(this, NoCPUs, ms1lcms, parameter.MS2PPM);
            swathdetection.MSlevel = 2;
            if (datattype != SpectralDataType.DataType.pSMART) {
                swathdetection.DetectPeakCurves(GetScanCollection());
                if (PeakCurveListMZ.isEmpty()) {
                    Logger.getRootLogger().info("No peak detected...................");
                    return;
                }
                if (ExportFragmentPeak) {
                    ExportPeakCurveResult();
                }
                ExportPeakCluster();
                GenerateRTSoretedClusterList(false);
                swathdetection.FragmentGrouping();
            }
            else{
                //////pSMART////////////////
                swathdetection.pSMARTGrouping(GetScanCollection());
            }
        }
        GenerateMGF(ms1lcms);
    }

    public void UmpireSearch(LCMSPeakMS1 ms1lcms, HashMap<Integer, Double> FactorialTable, PepIonLib IonLib, DBSearchParam searchpara) throws IOException, XmlPullParserException {

        IonLib.GetIonLib(DIA_MZ_Range.getX());
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        ArrayList<UmpireMatchUnit> ResultList = new ArrayList<>();
        int count = 0;
        for (PeakCluster ms1cluster : ms1lcms.PeakClusters) {
            if (DIA_MZ_Range.getX() <= ms1cluster.TargetMz() && DIA_MZ_Range.getY() >= ms1cluster.TargetMz()) {
                if (FragmentsClu2Cur.containsKey(ms1cluster.Index)) {
                    ArrayList<PrecursorFragmentPairEdge> fragments = FragmentsClu2Cur.get(ms1cluster.Index);
                    UmpireMatchUnit match = new UmpireMatchUnit(ms1cluster, fragments, IonLib, searchpara, FactorialTable);
                    ResultList.add(match);
                    executorPool.execute(match);
                }
            }
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ms1lcms.ScanCollectionName) + FilenameUtils.getBaseName(ms1lcms.ScanCollectionName) + "Search.txt", true);
        writer.write("ms1index\tms1mz\tIDseq\tmaxint\tSeq\tMod\tDecoy\tScore\tRank\tionmz\tBmatch\tYmatch\tBintensitysum\tYintensitysum\tBcorrsum\tYcorrsum\tBcoverage\tYcoverage\n");

        for (UmpireMatchUnit match : ResultList) {
            for (int i = 0; i < match.Ranking.size(); i++) {
                MS1PeakGroupSearch search = match.Ranking.get(i);
                writer.write(match.ms1cluster.Index + "\t" + match.ms1cluster.TargetMz() + "\t" + match.ms1cluster.AssignedPepIon + "\t" + search.maxint + "\t" + search.candidate.peptide.getSequence() + "\t" + search.candidate.GetModString() + "\t" + search.candidate.DecoyPep() + "\t" + search.CalScore() + "\t" + (i + 1) + "\t" + search.candidate.NeutralPrecursorMz() + "\t" + search.BMatch + "\t" + search.YMatch + "\t" + search.Bintensitysum + "\t" + search.Yintensitysum + "\t" + search.Bcorrsum + "\t" + search.Ycorrsum + "\t" + search.Bconverage + "\t" + search.Yconverage + "\n");
            }
            for (PepIonCandidate candidate : match.candidates) {
                candidate.ReleaseFragments();
            }
        }
        writer.close();
    }

    public String GetQ1Name() {
        return FilenameUtils.getBaseName(ParentmzXMLName) + "_Q1";
    }

    public String GetQ2Name() {
        return FilenameUtils.getBaseName(ParentmzXMLName) + "_Q2";
    }

    public String GetQ3Name() {
        return FilenameUtils.getBaseName(ParentmzXMLName) + "_Q3";
    }

    private void PrepareMGF_MS1Cluster(LCMSPeakMS1 ms1lcms) throws IOException {
        
        ArrayList<PseudoMSMSProcessing> ScanList=new ArrayList<>();
        ExecutorService executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PeakCluster ms1cluster : ms1lcms.PeakClusters) {
            if (DIA_MZ_Range.getX() <= ms1cluster.TargetMz() && DIA_MZ_Range.getY() >= ms1cluster.TargetMz() && FragmentsClu2Cur.containsKey(ms1cluster.Index)) {
                PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms1cluster, FragmentsClu2Cur.get(ms1cluster.Index), parameter);
                executorPool.execute(mSMSProcessing);
                ScanList.add(mSMSProcessing);                
            }
        }        
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }        
        
        String mgffile = FilenameUtils.getFullPath(ParentmzXMLName) + GetQ1Name() + ".mgf.temp";
        String mgffile2 = FilenameUtils.getFullPath(ParentmzXMLName) + GetQ2Name() + ".mgf.temp";
        FileWriter mapwriter = new FileWriter(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q1",true);
        FileWriter mapwriter2 = new FileWriter(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q2",true);
            
        FileWriter mgfWriter = new FileWriter(mgffile, true);
        FileWriter mgfWriter2 = new FileWriter(mgffile2, true);
        
        for (PseudoMSMSProcessing mSMSProcessing : ScanList) {
            if(MatchedFragmentMap.size()>0){
                mSMSProcessing.RemoveMatchedFrag(MatchedFragmentMap);
            }
            
            XYPointCollection Scan = mSMSProcessing.GetScan();

            if (Scan != null && Scan.PointCount() > parameter.MinFrag) {
                StringBuilder mgfString = new StringBuilder();
                if (mSMSProcessing.ms1cluster.IsotopeComplete(3)) {
                    parentDIA.Q1Scan++;
                    mgfString.append("BEGIN IONS\n");
                    mgfString.append("PEPMASS=").append(mSMSProcessing.ms1cluster.TargetMz()).append("\n");
                    mgfString.append("CHARGE=").append(mSMSProcessing.ms1cluster.Charge).append("+\n");
                    mgfString.append("RTINSECONDS=").append(mSMSProcessing.ms1cluster.PeakHeightRT[0] * 60f).append("\n");
                    mgfString.append("TITLE=").append(GetQ1Name()).append(".").append(parentDIA.Q1Scan).append(".").append(parentDIA.Q1Scan).append(".").append(mSMSProcessing.ms1cluster.Charge).append("\n");
                    for (int i = 0; i < Scan.PointCount(); i++) {
                        mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
                    }
                    mgfString.append("END IONS\n\n");                    
                    mapwriter.write(parentDIA.Q1Scan + "_" + mSMSProcessing.ms1cluster.Index + "\n");
                    mgfWriter.write(mgfString.toString());
                } else if (mSMSProcessing.ms1cluster.IsotopeComplete(2)) {
                    parentDIA.Q2Scan++;
                    mgfString.append("BEGIN IONS\n");
                    mgfString.append("PEPMASS=").append(mSMSProcessing.ms1cluster.TargetMz()).append("\n");
                    mgfString.append("CHARGE=").append(mSMSProcessing.ms1cluster.Charge).append("+\n");
                    mgfString.append("RTINSECONDS=").append(mSMSProcessing.ms1cluster.PeakHeightRT[0] * 60f).append("\n");
                    mgfString.append("TITLE=").append(GetQ2Name()).append(".").append(parentDIA.Q2Scan).append(".").append(parentDIA.Q2Scan).append(".").append(mSMSProcessing.ms1cluster.Charge).append("\n");
                    for (int i = 0; i < Scan.PointCount(); i++) {
                        mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
                    }
                    mgfString.append("END IONS\n\n");                    
                    mapwriter2.write(parentDIA.Q2Scan + "_" + mSMSProcessing.ms1cluster.Index + "\n");
                    mgfWriter2.write(mgfString.toString());
                }
            }
        }
        mgfWriter2.close();
        mgfWriter.close();
        mapwriter.close();
        mapwriter2.close();
    }

    private void PrepareMGF_UnfragmentIon() throws IOException {
        String mgffile4 = FilenameUtils.getFullPath(ParentmzXMLName) + GetQ3Name() + ".mgf.temp";
        FileWriter mgfWriter4 = new FileWriter(mgffile4, true);
        FileWriter mapwriter3 = new FileWriter(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + ".ScanClusterMapping_Q3",true);

        ArrayList<PseudoMSMSProcessing> ScanList = new ArrayList<>();
        ExecutorService executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PeakCluster ms1cluster : PeakClusters) {
            if (DIA_MZ_Range.getX() <= ms1cluster.TargetMz() && DIA_MZ_Range.getY() >= ms1cluster.TargetMz() && UnFragIonClu2Cur.containsKey(ms1cluster.Index)) {
                PseudoMSMSProcessing mSMSProcessing = new PseudoMSMSProcessing(ms1cluster, UnFragIonClu2Cur.get(ms1cluster.Index), parameter);
                executorPool.execute(mSMSProcessing);
                ScanList.add(mSMSProcessing);
            }
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }
        
        for (PseudoMSMSProcessing mSMSProcessing : ScanList) {
            if(MatchedFragmentMap.size()>0){
                mSMSProcessing.RemoveMatchedFrag(MatchedFragmentMap);
            }
            XYPointCollection Scan = mSMSProcessing.GetScan();

            if (Scan!=null && Scan.PointCount() > parameter.MinFrag) {
                parentDIA.Q3Scan++;
                StringBuilder mgfString = new StringBuilder();
                mgfString.append("BEGIN IONS\n");
                mgfString.append("PEPMASS=" + mSMSProcessing.ms1cluster.TargetMz() + "\n");
                mgfString.append("CHARGE=" + mSMSProcessing.ms1cluster.Charge + "+\n");
                mgfString.append("RTINSECONDS=" + mSMSProcessing.ms1cluster.PeakHeightRT[0] * 60f + "\n");
                mgfString.append("TITLE=").append(GetQ3Name()).append(".").append(parentDIA.Q3Scan).append(".").append(parentDIA.Q3Scan).append(".").append(mSMSProcessing.ms1cluster.Charge).append("\n");
                //mgfString.append("TITLE=" + WindowID + ";ClusterIndex:" + mSMSProcessing.ms1cluster.Index + "\n");
                //mgfString.append("TITLE=" GetQ3Name() + WindowID + ";ClusterIndex:" + mSMSProcessing.ms1cluster.Index + "\n");
                
                for (int i = 0; i < Scan.PointCount(); i++) {
                    mgfString.append(Scan.Data.get(i).getX()).append(" ").append(Scan.Data.get(i).getY()).append("\n");
                }
                mgfString.append("END IONS\n\n");
                mgfWriter4.write(mgfString.toString());
                mapwriter3.write(parentDIA.Q3Scan + ";" + WindowID+";"+mSMSProcessing.ms1cluster.Index + "\n");
            }
        }
        mgfWriter4.close();
        mapwriter3.close();
    }

    public void GenerateMGF(LCMSPeakMS1 ms1lcms) throws IOException, InterruptedException {
        PrepareMGF_MS1Cluster(ms1lcms);
        PrepareMGF_UnfragmentIon();
    }
  

    public ScanCollection GetScanCollection() throws InterruptedException, ExecutionException, IOException {
        return mzxml.GetScanCollectionDIAMS2(DIA_MZ_Range, true);
    }

    public boolean ReadIfProcessed() {

        return (ReadPrecursorFragmentClu2Cur() & ReadPeakCluster());
//        
//        Connection connection = connectionManager.GetConnection();
//        Statement state = connection.createStatement();
//        ResultSet rsCluster = state.executeQuery("SHOW TABLES LIKE '" + FilenameUtils.getBaseName(ScanCollectionName) + "_SwathPeakCorr'");
//        if (rsCluster.next()) {
//            if (FragmentGroupByCluster) {
//                ReadFragmentClu2CluFromDB();
//            } else {
//                ReadFragmentClu2CurFromDB();
//                ReadUnFragCurveClusterResultFromDB();
//            }
//            return true;
//        }
//        return false;
    }

    private void ReadFragmentClu2CluFromDB() throws SQLException {
        Connection connection = connectionManager.GetConnection();
        Statement state = connection.createStatement();
        FragmentsClu2Clu = new HashMap<>();
        //state.execute("Delete From "+FilenameUtils.getBaseName(ScanCollectionName) + "_SwathPeakCorr where Corr<0.2");
        ResultSet rsCluster = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Clus");
        while (rsCluster.next()) {
            PrecursorFragmentPairEdge framentClusterUnit = new PrecursorFragmentPairEdge();
            framentClusterUnit.PeakCurveIndexA = rsCluster.getInt("PeakClusterA");
            framentClusterUnit.PeakCurveIndexB = rsCluster.getInt("PeakClusterB");
            framentClusterUnit.Correlation = rsCluster.getFloat("Corr");
            framentClusterUnit.FragmentMz = rsCluster.getFloat("FragmentMz");
            framentClusterUnit.Intensity = rsCluster.getFloat("FragInt");
            //framentClusterUnit.FragmentGroupidx=rsCluster.getInt("FragGrpIdx");
            //framentClusterUnit.FragmentConsensusCorr=rsCluster.getFloat("FragCorr");

            if (!FragmentsClu2Clu.containsKey(framentClusterUnit.PeakCurveIndexA)) {
                ArrayList<PrecursorFragmentPairEdge> list = new ArrayList<>();
                FragmentsClu2Clu.put(framentClusterUnit.PeakCurveIndexA, list);
            }
            FragmentsClu2Clu.get(framentClusterUnit.PeakCurveIndexA).add(framentClusterUnit);
        }
        state.close();
        connectionManager.CloseConnection();
    }

    HashMap<Integer, ArrayList<Float>> FragmentMS1Ranking;

    public void BuildFragmentMS1ranking() {
        FragmentMS1Ranking = new HashMap<>();
        for (int clusterindex : FragmentsClu2Cur.keySet()) {
            for (PrecursorFragmentPairEdge framentClusterUnit : FragmentsClu2Cur.get(clusterindex)) {
                if (!FragmentMS1Ranking.containsKey(framentClusterUnit.PeakCurveIndexB)) {
                    ArrayList<Float> scorelist = new ArrayList<>();
                    FragmentMS1Ranking.put(framentClusterUnit.PeakCurveIndexB, scorelist);
                }
                FragmentMS1Ranking.get(framentClusterUnit.PeakCurveIndexB).add(framentClusterUnit.Correlation);
            }
        }
        for (ArrayList<Float> scorelist : FragmentMS1Ranking.values()) {
            Collections.sort(scorelist);
            Collections.reverse(scorelist);
        }
        for (int clusterindex : FragmentsClu2Cur.keySet()) {
            for (PrecursorFragmentPairEdge framentClusterUnit : FragmentsClu2Cur.get(clusterindex)) {
                ArrayList<Float> scorelist = FragmentMS1Ranking.get(framentClusterUnit.PeakCurveIndexB);
                for (int intidx = 0; intidx < scorelist.size(); intidx++) {
                    if (scorelist.get(intidx) <= framentClusterUnit.Correlation) {
                        framentClusterUnit.FragmentMS1Rank = intidx + 1;
                        framentClusterUnit.FragmentMS1RankScore = (float) framentClusterUnit.FragmentMS1Rank / (float) scorelist.size();
                        break;
                    }
                }
            }
        }
    }

    public void ReadPeakFromDB() {
        try {
            ReadPeakClusterFromDB();
            ReadFragmentClu2CurFromDB();
            ReadUnFragCurveClusterResultFromDB();
        } catch (SQLException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public boolean ReadPrecursorFragmentClu2Cur() {
        return ReadCluster2CurveCorrSerialization() & ReadUnfragmentedCluster2CurveCorrSerialization();
    }

    //Read from MySQL: no longer needed
    private void ReadFragmentClu2CurFromDB() throws SQLException {
        //System.out.print("Loading precursor-fragment pairs table........");
        Connection connection = connectionManager.GetConnection();
        Statement state = connection.createStatement();
        FragmentsClu2Cur = new HashMap<>();
        //state.execute("Delete From "+FilenameUtils.getBaseName(ScanCollectionName) + "_SwathPeakCorr where Corr<0.2");
        ResultSet rsCluster = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur");
        while (rsCluster.next()) {
            PrecursorFragmentPairEdge framentClusterUnit = new PrecursorFragmentPairEdge();
            framentClusterUnit.PeakCurveIndexA = rsCluster.getInt("PeakClusterA");
            framentClusterUnit.PeakCurveIndexB = rsCluster.getInt("PeakClusterB");
            framentClusterUnit.Correlation = rsCluster.getFloat("Corr");
            framentClusterUnit.FragmentMz = rsCluster.getFloat("FragmentMz");
            framentClusterUnit.Intensity = rsCluster.getFloat("FragInt");
            framentClusterUnit.ApexDelta = rsCluster.getFloat("ApexDelta");
            framentClusterUnit.RTOverlapP = rsCluster.getFloat("RTOverlapP");
            //framentClusterUnit.FragmentGroupidx=rsCluster.getInt("FragGrpIdx");
            //framentClusterUnit.FragmentConsensusCorr=rsCluster.getFloat("FragCorr");

            if (!FragmentsClu2Cur.containsKey(framentClusterUnit.PeakCurveIndexA)) {
                ArrayList<PrecursorFragmentPairEdge> list = new ArrayList<>();
                FragmentsClu2Cur.put(framentClusterUnit.PeakCurveIndexA, list);
            }
            FragmentsClu2Cur.get(framentClusterUnit.PeakCurveIndexA).add(framentClusterUnit);
        }

        for (ArrayList<PrecursorFragmentPairEdge> list : FragmentsClu2Cur.values()) {
            Collections.sort(list, new Comparator<PrecursorFragmentPairEdge>() {
                @Override
                public int compare(PrecursorFragmentPairEdge o1, PrecursorFragmentPairEdge o2) {
                    return Float.compare(o1.FragmentMz, o2.FragmentMz);
                }
            });
        }
        state.close();
        connectionManager.CloseConnection();
        //System.out.print("done\n");
    }

    private void ReadUnFragCurveClusterResultFromDB() throws SQLException {
        //System.out.print("Loading precursor-fragment pairs (unfragmented ions) table........");
        Connection connection = connectionManager.GetConnection();
        Statement state = connection.createStatement();
        UnFragIonClu2Cur = new HashMap<>();
        //state.execute("Delete From "+FilenameUtils.getBaseName(ScanCollectionName) + "_SwathPeakCorr where Corr<0.2");
        ResultSet rsCluster = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur");
        while (rsCluster.next()) {
            PrecursorFragmentPairEdge framentClusterUnit = new PrecursorFragmentPairEdge();
            framentClusterUnit.PeakCurveIndexA = rsCluster.getInt("PeakClusterA");
            framentClusterUnit.PeakCurveIndexB = rsCluster.getInt("PeakClusterB");
            framentClusterUnit.Correlation = rsCluster.getFloat("Corr");
            framentClusterUnit.FragmentMz = rsCluster.getFloat("FragmentMz");
            framentClusterUnit.Intensity = rsCluster.getFloat("FragInt");
            framentClusterUnit.ApexDelta = rsCluster.getFloat("ApexDelta");
            framentClusterUnit.RTOverlapP = rsCluster.getFloat("RTOverlapP");
            //framentClusterUnit.FragmentGroupidx=rsCluster.getInt("FragGrpIdx");
            //framentClusterUnit.FragmentConsensusCorr=rsCluster.getFloat("FragCorr");

            if (!UnFragIonClu2Cur.containsKey(framentClusterUnit.PeakCurveIndexA)) {
                ArrayList<PrecursorFragmentPairEdge> list = new ArrayList<>();
                UnFragIonClu2Cur.put(framentClusterUnit.PeakCurveIndexA, list);
            }
            UnFragIonClu2Cur.get(framentClusterUnit.PeakCurveIndexA).add(framentClusterUnit);
        }

        for (ArrayList<PrecursorFragmentPairEdge> list : UnFragIonClu2Cur.values()) {
            Collections.sort(list, new Comparator<PrecursorFragmentPairEdge>() {
                @Override
                public int compare(PrecursorFragmentPairEdge o1, PrecursorFragmentPairEdge o2) {
                    return Float.compare(o1.FragmentMz, o2.FragmentMz);
                }
            });
        }
        state.close();
        connectionManager.CloseConnection();
        //System.out.println("done");
    }

    public void WritePrecursorFragmentGrouping() {
        WriteCluster2CurveCorrSerialization();
        WriteUnfragmentedCluster2CurveCorrSerialization();
    }

    private void WriteCluster2CurveCorrSerialization() {
        //JavaSerializationCluster2CurveWrite();
        FSCluster2CurveWrite();
    }

    private void JavaSerializationCluster2CurveWrite() {
        try {
            Logger.getRootLogger().debug("Writing PrecursorFragmentCorr serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_Clus2Cur.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(FragmentsClu2Cur);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private void FSCluster2CurveWrite() {
        try {
            Logger.getRootLogger().debug("Writing PrecursorFragmentCorr serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(FragmentsClu2Cur);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private boolean ReadCluster2CurveCorrSerialization() {
        //return JavaSerializationCluster2CurveRead();   
        if (!FSCluster2CurveRead()) {
            if (JavaSerializationCluster2CurveRead()) {
                FSCluster2CurveWrite();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSCluster2CurveRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_Clus2Cur.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("Reading PrecursorFragmentCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_Clus2Cur.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            FragmentsClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            if (FSCluster2CurveRead_v158()) {
                WriteCluster2CurveCorrSerialization();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSCluster2CurveRead_v158() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_Clus2Cur.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("v158 PrecursorFragmentCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.serFS");
            de.ruedigermoeller.serialization.FSTObjectInput in = new de.ruedigermoeller.serialization.FSTObjectInput(fileIn);
            FragmentsClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().error("v158 still failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private boolean JavaSerializationCluster2CurveRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("Reading PrecursorFragmentCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            FragmentsClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private void WriteUnfragmentedCluster2CurveCorrSerialization() {
        //JavaSerializationCluster2CurveUnfragWrite();
        FSCluster2CurveUnfragWrite();
    }

    private void FSCluster2CurveUnfragWrite() {
        try {
            Logger.getRootLogger().debug("Writing UnfragPrecursorFragCorr serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_UnfClus2Cur.serFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(UnFragIonClu2Cur);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private void JavaSerializationCluster2CurveUnfragWrite() {
        try {
            Logger.getRootLogger().debug("Writing UnfragPrecursorFragCorr serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_UnfClus2Cur.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(UnFragIonClu2Cur);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private boolean ReadUnfragmentedCluster2CurveCorrSerialization() {
        if (!FSCluster2CurveUnfragRead()) {
            if (JavaSerializationCluster2CurveUnfragRead()) {
                FSCluster2CurveUnfragWrite();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSCluster2CurveUnfragRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("Reading UnfragPrecursorFragCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            UnFragIonClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            if (FSCluster2CurveUnfragRead_v158()) {
                WriteUnfragmentedCluster2CurveCorrSerialization();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSCluster2CurveUnfragRead_v158() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("v158 UnfragPrecursorFragCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.serFS");
            de.ruedigermoeller.serialization.FSTObjectInput in = new de.ruedigermoeller.serialization.FSTObjectInput(fileIn);
            UnFragIonClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().error("v158 still failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private boolean JavaSerializationCluster2CurveUnfragRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("Reading UnfragPrecursorFragCorr serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            UnFragIonClu2Cur = (HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>>) in.readObject();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    public void ExportCluster2CurveCorr() throws IOException, SQLException {
        WriteCluster2CurveCorrSerialization();
        if (connectionManager != null) {
            ExportCluster2CurveCorrCSV();
            CreateCluster2CurveCorrTable();
            ExportTableDBBulkLoader("Clus2Cur",false);
        }
    }

    private void ExportCluster2ClusterCorr(ArrayList<CorrCalcCluster2ClusterUnit> ResultArrayList) throws IOException, SQLException {
        if (connectionManager != null) {
            ExportCluster2ClusterCorrCSV(ResultArrayList);
            CreateCluster2ClusterCorrTable();
            ExportTableDBBulkLoader("Clus2Clus",false);
        }
    }

    private void ExportUnfragmentedCluster2Cluster(ArrayList<CorrCalcCluster2ClusterUnit> ResultArrayList) throws IOException, SQLException {
        if (connectionManager != null) {
            ExportUnfragmentedCluster2ClusterCSV(ResultArrayList);
            CreateUnfClus2ClusTable();
            ExportTableDBBulkLoader("UnfClus2Clus",false);
        }
    }

    public void ExportUnfragmentedClusterCurve() throws IOException, SQLException {
        WriteUnfragmentedCluster2CurveCorrSerialization();
        if (connectionManager != null) {
            ExportUnfragmentedClusterCurveCSV();
            CreateUnfClus2CurTable();
            ExportTableDBBulkLoader("UnfClus2Cur",false);
        }
    }

    private void ExportIsolatedClusterCurve(ArrayList<IsolatedCurveCorrCalcUnit> ResultArrayList) throws IOException, SQLException {
        ExportIsolatedClusterCurveCSV(ResultArrayList);
        CreateCorrIsolatedCurveMatrixTable();
        ExportTableDBBulkLoader("IPeakMS2Corr",false);
    }

    public void ExportCluster2CurveCorrCSV() throws IOException {
        Logger.getRootLogger().debug("Writing PrecursorFragmentCorr result to file:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur.csv");
        writer.write("ID,PeakClusterA,PeakClusterB,Corr,FragmentMz,FragInt,RTOverlapP,ApexDelta\n");
        int id = 0;
        for (ArrayList<PrecursorFragmentPairEdge> result : FragmentsClu2Cur.values()) {
            for (PrecursorFragmentPairEdge region : result) {
                writer.write((id++) + "," + region.PeakCurveIndexA + "," + region.PeakCurveIndexB + "," + region.Correlation + "," + region.FragmentMz + "," + region.Intensity + "," + region.RTOverlapP + "," + region.ApexDelta/*+ "," + region.FragmentConsensusCorr + "," + region.FragmentGroupidx*/ + "\n");
            }
        }
        writer.close();
    }

    private void ExportCluster2ClusterCorrCSV(ArrayList<CorrCalcCluster2ClusterUnit> ResultArrayList) throws IOException {
        Logger.getRootLogger().debug("Writing PrecursorFragmentCorr result to file:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Clus.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Clus.csv");
        writer.write("ID,PeakClusterA,PeakClusterB,Corr,FragmentMz,FragInt\n");
        int id = 0;
        for (CorrCalcCluster2ClusterUnit result : ResultArrayList) {
            for (PrecursorFragmentPairEdge region : result.ResultList) {
                writer.write((id++) + "," + region.PeakCurveIndexA + "," + region.PeakCurveIndexB + "," + region.Correlation + "," + region.FragmentMz + "," + region.Intensity /*+ "," + region.FragmentConsensusCorr + "," + region.FragmentGroupidx*/ + "\n");
            }
        }
        writer.close();
    }

    private void ExportUnfragmentedCluster2ClusterCSV(ArrayList<CorrCalcCluster2ClusterUnit> ResultArrayList) throws IOException {
        Logger.getRootLogger().debug("Writing UnfragPeakCorr result to file:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Clus.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Clus.csv");
        writer.write("ID,PeakClusterA,PeakClusterB,Corr,FragmentMz,FragInt\n");
        int id = 0;
        for (CorrCalcCluster2ClusterUnit result : ResultArrayList) {
            for (PrecursorFragmentPairEdge region : result.ResultList) {
                writer.write((id++) + "," + region.PeakCurveIndexA + "," + region.PeakCurveIndexB + "," + region.Correlation + "," + region.FragmentMz + "," + region.Intensity /*+ "," + region.FragmentConsensusCorr + "," + region.FragmentGroupidx*/ + "\n");
            }
        }
        writer.close();
    }

    private void ExportUnfragmentedClusterCurveCSV() throws IOException {
        Logger.getRootLogger().debug("Writing UnfragPeakCorr result to file:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur.csv");
        writer.write("ID,PeakClusterA,PeakClusterB, Corr,FragmentMz,FragInt,RTOverlapP,ApexDelta\n");
        int id = 0;
        for (ArrayList<PrecursorFragmentPairEdge> result : UnFragIonClu2Cur.values()) {
            for (PrecursorFragmentPairEdge region : result) {
                writer.write((id++) + "," + region.PeakCurveIndexA + "," + region.PeakCurveIndexB + "," + region.Correlation + "," + region.FragmentMz + "," + region.Intensity + "," + region.RTOverlapP + "," + region.ApexDelta/*+ "," + region.FragmentConsensusCorr + "," + region.FragmentGroupidx*/ + "\n");
            }
        }
        writer.close();
    }

    private void ExportIsolatedClusterCurveCSV(ArrayList<IsolatedCurveCorrCalcUnit> ResultArrayList) throws IOException {
        Logger.getRootLogger().debug("Writing IPeakMS2Corr result to file:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_IPeakMS2Corr.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_IPeakMS2Corr.csv");
        writer.write("ID,PeakClusterA,PeakClusterB,Corr,FragmentMz,FragInt\n");
        int id = 0;
        for (IsolatedCurveCorrCalcUnit result : ResultArrayList) {
            for (PrecursorFragmentPairEdge region : result.ResultList) {
                writer.write((id++) + "," + region.PeakCurveIndexA + "," + region.PeakCurveIndexB + "," + region.Correlation + "," + region.FragmentMz + "," + region.Intensity /*+ "," + region.FragmentConsensusCorr + "," + region.FragmentGroupidx*/ + "\n");
            }
        }
        writer.close();
    }

    protected void CreateUnfClus2ClusTable() throws SQLException {
        //super.CreateCorrMatrixTable();
        if (connectionManager == null) {
            return;
        }
        Connection connection = connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Clus;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Clus (ID INT NOT NULL AUTO_INCREMENT , PeakClusterA INT NOT NULL, PeakClusterB INT NOT NULL, Corr DOUBLE NOT NULL, FragmentMz DOUBLE NOT NULL, FragInt DOUBLE NOT NULL, PRIMARY KEY (ID));");
        connectionManager.CloseConnection();
    }

    protected void CreateUnfClus2CurTable() throws SQLException {
        //super.CreateCorrMatrixTable();
        if (connectionManager == null) {
            return;
        }
        Connection connection = connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_UnfClus2Cur (ID INT NOT NULL AUTO_INCREMENT , PeakClusterA INT NOT NULL, PeakClusterB INT NOT NULL, Corr DOUBLE NOT NULL, FragmentMz DOUBLE NOT NULL, FragInt DOUBLE NOT NULL,RTOverlapP DOUBLE NOT NULL,ApexDelta DOUBLE NOT NULL, PRIMARY KEY (ID));");
        connectionManager.CloseConnection();
    }

    protected void CreateCluster2CurveCorrTable() throws SQLException {
        //super.CreateCorrMatrixTable();
        if (connectionManager == null) {
            return;
        }
        Connection connection = connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Cur (ID INT NOT NULL AUTO_INCREMENT , PeakClusterA INT NOT NULL, PeakClusterB INT NOT NULL, Corr DOUBLE NOT NULL, FragmentMz DOUBLE NOT NULL, FragInt DOUBLE NOT NULL,RTOverlapP DOUBLE NOT NULL,ApexDelta DOUBLE NOT NULL, PRIMARY KEY (ID));");
        connectionManager.CloseConnection();
    }

    protected void CreateCluster2ClusterCorrTable() throws SQLException {
        //super.CreateCorrMatrixTable();
        if (connectionManager == null) {
            return;
        }
        Connection connection = connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Clus;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_Clus2Clus (ID INT NOT NULL AUTO_INCREMENT , PeakClusterA INT NOT NULL, PeakClusterB INT NOT NULL, Corr DOUBLE NOT NULL, FragmentMz DOUBLE NOT NULL, FragInt DOUBLE NOT NULL, PRIMARY KEY (ID));");
        connectionManager.CloseConnection();
    }

    protected void CreateCorrIsolatedCurveMatrixTable() throws SQLException {
        //super.CreateCorrMatrixTable();
        if (connectionManager == null) {
            return;
        }
        Connection connection = connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_IPeakMS2Corr;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_IPeakMS2Corr (ID INT NOT NULL AUTO_INCREMENT , PeakClusterA INT NOT NULL, PeakClusterB INT NOT NULL, Corr DOUBLE NOT NULL, FragmentMz DOUBLE NOT NULL, FragInt DOUBLE NOT NULL, PRIMARY KEY (ID));");
        connectionManager.CloseConnection();
    }

    protected void CreateCorrCurveMatrixTable() throws SQLException {
        //super.CreateCorrMatrixTable();
        if (connectionManager == null) {
            return;
        }
        Connection connection = connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_SwathPeakCorr;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_SwathPeakCorr (ID INT NOT NULL AUTO_INCREMENT , PeakClusterA INT NOT NULL, PeakClusterB INT NOT NULL, AStart INT NOT NULL, Aend INT NOT NULL, Bstart INT NOT NULL, Bend INT NOT NULL,Corr DOUBLE NOT NULL, FragmentMz DOUBLE NOT NULL, FragInt DOUBLE NOT NULL, PRIMARY KEY (ID));");
        connectionManager.CloseConnection();
    }

    public void ExtractFragmentForPeakCluser(PeakCluster Cluster) {        
        if (FragmentsClu2Cur.containsKey(Cluster.Index)) {
            Cluster.fraglock.writeLock().lock();
            try {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : FragmentsClu2Cur.get(Cluster.Index)) {
                    if (!Cluster.Fragments.contains(fragmentClusterUnit)) {
                        Cluster.Fragments.add(fragmentClusterUnit);
                    }
                }
            } finally {
                Cluster.fraglock.writeLock().unlock();
            }
        }
    }

    public void ExtractFragmentForUnfragPeakCluser(PeakCluster Cluster) {
        if (UnFragIonClu2Cur.containsKey(Cluster.Index)) {
            Cluster.fraglock.writeLock().lock();
            try {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : UnFragIonClu2Cur.get(Cluster.Index)) {
                    if (!Cluster.Fragments.contains(fragmentClusterUnit)) {
                        Cluster.Fragments.add(fragmentClusterUnit);
                    }
                }
            } finally {
                Cluster.fraglock.writeLock().unlock();
            }
        }
    }
}
