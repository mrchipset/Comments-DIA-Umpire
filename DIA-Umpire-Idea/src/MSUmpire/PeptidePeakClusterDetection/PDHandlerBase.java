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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYZData;
import MSUmpire.LCMSBaseStructure.LCMSPeakBase;
import MSUmpire.MySQLTool.ConnectionManager;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedClusterCollectionClassApexRT;
import MSUmpire.PeakDataStructure.SortedClusterCollectionClassMZ;
import Utility.UpdateProcess;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerBase {
    
    protected HashSet<String> IncludedHashMap;
    protected int NoCPUs = 4;
    public float minSNR;
    public TreeMap<Float, XYData>[] IsotopePatternMap;
    public TreeMap<Float, XYData>[] IsotopePatternFragMap;
    protected LCMSPeakBase LCMSPeakBase;
    protected InstrumentParameter parameter;
    protected ConnectionManager connectionManager;
    protected boolean ReleaseScans = true;
    protected float PPM;
    public int MSlevel=1;

    public PDHandlerBase() {
    }

    public void ClearAllPeaks() {
        LCMSPeakBase.ClearAllPeaks();
    }

    public ArrayList<PeakCurve> FindPeakCurveBymzRT(float mz, float RT, float ppm) {
        ArrayList<PeakCurve> ReturnList = new ArrayList<>();

        int startidx = LCMSPeakBase.PeakCurveListMZ.BinarySearchLower(InstrumentParameter.GetMzByPPM(mz, 1, ppm));
        for (int i = startidx; i < LCMSPeakBase.PeakCurveListMZ.size(); i++) {
            PeakCurve peakcurve = LCMSPeakBase.PeakCurveListMZ.get(i);

            if (InstrumentParameter.CalcPPM(peakcurve.TargetMz, mz) < ppm) {
                if (peakcurve.StartRT() <= RT && peakcurve.EndRT() >= RT) {
                    ReturnList.add(peakcurve);
                }
            } else {
                if (peakcurve.TargetMz > mz) {
                    return ReturnList;
                }
            }
        }
        return ReturnList;
    }

    
    public void ClearRawPeaks() {
        for (PeakCurve peakCurve : LCMSPeakBase.PeakCurveListMZ) {
            peakCurve.CalculateMzVar();
            peakCurve.ReleaseRawPeak();
        }
    }

    protected void FindAllPeakCurve(ScanCollection scanCollection) throws SQLException, IOException {

        IncludedHashMap = new HashSet<>();
        Logger.getRootLogger().info("Processing all scans to detect possible peak curves....");

        //Get the ms1 scanNo array
        //LCMSPeakBase.PeakCurveListMZ = new SortedCurveCollectionMZ();
        
        float preRT = 0f;
        for (int idx = 0; idx < scanCollection.GetScanNoArray(MSlevel).size(); idx++) {
            Integer scanNO = scanCollection.GetScanNoArray(MSlevel).get(idx);
            ScanData scanData = scanCollection.GetScan(scanNO);
            if (idx == 0) {
                preRT = scanData.RetentionTime - 0.01f;
            }
            for (int i = 0; i < scanData.PointCount(); i++) {
                XYData peak = scanData.Data.get(i);
                
                //Include the mz-int pair to a hash
                if (!IncludedHashMap.contains(scanNO + "_" + peak.getX())) {//The peak hasn't been included and checked
                    IncludedHashMap.add(scanNO + "_" + peak.getX());

                    float startmz = peak.getX();
                    float startint = peak.getY();
                    for (int j = i + 1; j < scanData.PointCount(); j++) {
                        XYData currentpeak = scanData.Data.get(j);
                        if (!IncludedHashMap.contains(scanNO + "_" + currentpeak.getX())) {
                            if (InstrumentParameter.CalcPPM(currentpeak.getX(), startmz) <= PPM) {
                                IncludedHashMap.add(scanNO + "_" + currentpeak.getX());

                                if (currentpeak.getY() >= startint) {
                                    startmz = currentpeak.getX();
                                    startint = currentpeak.getY();
                                }
                            } else {
                                break;
                            }
                        }
                    }

                    PeakCurve Peakcurve = new PeakCurve(parameter);
                    Peakcurve.AddPeak(new XYZData(preRT, startmz, scanData.background));
                    Peakcurve.AddPeak(new XYZData(scanData.RetentionTime, startmz, startint));
                    Peakcurve.StartScan = scanNO;

                    int missedScan = 0;

                    //Start with the next MS1 scan to group the mz-int pair within the N ppm window
                    for (int idx2 = idx + 1; idx2 < scanCollection.GetScanNoArray(MSlevel).size() && missedScan < parameter.NoMissedScan; idx2++) {
                        Integer scanNO2 = scanCollection.GetScanNoArray(MSlevel).get(idx2);
                        ScanData scanData2 = scanCollection.GetScan(scanNO2);

                        float currentmz = 0f;
                        float currentint = 0f;

                        if (scanData2.PointCount() == 0) {
                            Peakcurve.AddPeak(new XYZData(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background));
                            missedScan++;
                            continue;
                        }

                        int mzidx = scanData2.GetLowerIndexOfX(Peakcurve.TargetMz);
                        for (int pkidx = mzidx; pkidx < scanData2.Data.size(); pkidx++) {
                            XYData currentpeak = scanData2.Data.get(pkidx);
                            if (!IncludedHashMap.contains(scanNO2 + "_" + currentpeak.getX())) {
                                if (InstrumentParameter.CalcPPM(currentpeak.getX(), Peakcurve.TargetMz) > PPM) {
                                    if (currentpeak.getX() > Peakcurve.TargetMz) {
                                        break;
                                    }
                                } else {
                                    //////////The peak is in the ppm window, select the highest peak
                                    IncludedHashMap.add(scanNO2 + "_" + currentpeak.getX());
                                    if (currentint < currentpeak.getY()) {
                                        currentmz = currentpeak.getX();
                                        currentint = currentpeak.getY();
                                    }
                                }
                            }
                        }
                        if (currentmz == 0f) {
                            Peakcurve.AddPeak(new XYZData(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background));
                            missedScan++;
                        } else {
                            missedScan = 0;
                            Peakcurve.AddPeak(new XYZData(scanData2.RetentionTime, currentmz, currentint));
                            Peakcurve.EndScan = scanNO2;
                        }
                    }
//                    if (Peakcurve.TargetMz > 870.5 && Peakcurve.TargetMz < 872.7 && Peakcurve.StartRT() < 48.8 && Peakcurve.EndRT() > 48.8) {
//                        System.out.println("");
//                    }

                    if (Peakcurve.GetRawSNR() > LCMSPeakBase.SNR && Peakcurve.GetPeakList().size() >= parameter.NoMissedScan + parameter.MinPeakPerPeakCurve+2) {
                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
                    } else {
                        Peakcurve = null;
                    }
                }
            }
            preRT = scanData.RetentionTime;
            if (ReleaseScans) {
                scanData.dispose();
            }
        }

        //System.out.print("PSM removed (PeakCurve generation):" + PSMRemoved );         
        IncludedHashMap.clear();
        IncludedHashMap = null;

        System.gc();
        Logger.getRootLogger().info(LCMSPeakBase.UnSortedPeakCurves.size() + " Peak curves found (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
        //writer.close();
    }

    protected void WaveletDetectMax() {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors");
        Logger.getRootLogger().info("Performing CWT to detect peak regions.....");
        
        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //executorPool = Executors.newFixedThreadPool(1);
        //progress.SetTotal(LCMSPeakBase.PeakCurveListMZ.size());
        //progress.SetTotal(LCMSPeakBase.UnSortedPeakCurves.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        ArrayList<WaveletRegionDetection> ResultList = new ArrayList<>();
        for (PeakCurve Peakcurve : LCMSPeakBase.UnSortedPeakCurves) {
            //if (Peakcurve.StartRT() < 32.6 && Peakcurve.EndRT() > 33 && Peakcurve.TargetMz > 322.68 && Peakcurve.TargetMz < 322.7) {          
             //if (Peakcurve.TargetMz > 870.5 && Peakcurve.TargetMz < 872.7 && Peakcurve.StartRT() < 48.8 && Peakcurve.EndRT() > 48.8) {                        
                    
            WaveletRegionDetection unit = new WaveletRegionDetection(Peakcurve, parameter, progress);
            ResultList.add(unit);

            //unit.run();
            executorPool.execute(unit);
            //}
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        executorPool = null;
        for (WaveletRegionDetection result : ResultList) {
            LCMSPeakBase.PeakCurveListMZ.addAll(result.ResultCurves);
            LCMSPeakBase.PeakCurveListRT.addAll(result.ResultCurves);
        }

        LCMSPeakBase.PeakCurveListMZ.Finalize();
        for (int i = 0; i < LCMSPeakBase.PeakCurveListMZ.size(); i++) {
            LCMSPeakBase.PeakCurveListMZ.get(i).Index = i + 1;
        }
        LCMSPeakBase.PeakCurveListRT.Finalize();
        LCMSPeakBase.PeakCurveListMZ.Finalize();
        LCMSPeakBase.UnSortedPeakCurves.clear();
        LCMSPeakBase.UnSortedPeakCurves = null;
        System.gc();
        Logger.getRootLogger().info(LCMSPeakBase.PeakCurveListMZ.size() + " peak curves left (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
    }

    protected void PeakCurveCorrClustering_V2(XYData mzRange) throws SQLException, IOException {
        Logger.getRootLogger().info("Grouping isotopic peak curves........");

        LCMSPeakBase.PeakClusters = new ArrayList<>();
        LCMSPeakBase.MZSortedClusters = new SortedClusterCollectionClassMZ();
        LCMSPeakBase.ApexRTSortedClusters = new SortedClusterCollectionClassApexRT();

        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //executorPool = Executors.newFixedThreadPool(1);
        ArrayList<PeakCurveClusteringCorrV2Unit> ResultList = new ArrayList<>();

        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;

        //progress.SetTotal(LCMSPeakBase.PeakCurveListMZ.size());
        //Thread thread = new Thread(progress);
        //thread.start();

        for (int i = 0; i < LCMSPeakBase.PeakCurveListMZ.size(); i++) {
            PeakCurve Peakcurve = LCMSPeakBase.PeakCurveListMZ.get(i);
            if (Peakcurve.TargetMz >= mzRange.getX() && Peakcurve.TargetMz <= mzRange.getY()) {                
                //if (Peakcurve.TargetMz > 870.5 && Peakcurve.TargetMz < 872.7 && Peakcurve.StartRT() < 48.8 && Peakcurve.EndRT() > 48.8) {
                    PeakCurveClusteringCorrV2Unit unit = new PeakCurveClusteringCorrV2Unit(Peakcurve, LCMSPeakBase.PeakCurveListMZ, LCMSPeakBase.PeakCurveListRT, parameter, IsotopePatternMap, LCMSPeakBase.StartCharge, LCMSPeakBase.EndCharge, LCMSPeakBase.MaxNoPeakCluster, LCMSPeakBase.MinNoPeakCluster, progress);
                    ResultList.add(unit);
                    executorPool.execute(unit);
                //}
            }
        }

        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }

        //thread = null;
        //progress.ClearMSG();
        //progress = null;

        for (PeakCurveClusteringCorrV2Unit unit : ResultList) {
            for (PeakCluster peakCluster : unit.ResultClusters) {
//                 PeakCurve Peakcurve=peakCluster.MonoIsotopePeak;                            
//                                if (Peakcurve.TargetMz > 697.3 && Peakcurve.TargetMz < 698.3 && Peakcurve.StartRT() < 48 && Peakcurve.EndRT() > 47.8) {
//                                    System.out.println("");
//                                }
                if (!parameter.RemoveGroupedPeaks || !peakCluster.MonoIsotopePeak.ChargeGrouped.contains(peakCluster.Charge)) {
                    peakCluster.Index = LCMSPeakBase.PeakClusters.size() + 1;
                    peakCluster.GetConflictCorr();
                    LCMSPeakBase.PeakClusters.add(peakCluster);
                }
                //}               
            }
        }
        
        ////////////////////////////
//        for(PeakCluster cluster : LCMSPeakBase.PeakClusters){
//                                PeakCurve Peakcurve=cluster.MonoIsotopePeak;                            
//                                if (Peakcurve.TargetMz > 697.3 && Peakcurve.TargetMz < 698.3 && Peakcurve.StartRT() < 48 && Peakcurve.EndRT() > 47.8) {
//                                    System.out.println("");
//                                }
//                            }
//        /////////////////////////////
        ResultList.clear();
        ResultList = null;
        System.gc();
        Logger.getRootLogger().info("No of ion clusters:" + LCMSPeakBase.PeakClusters.size() + " (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");

    }
}
