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

import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.LCMSBaseStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_window_Quant implements  Runnable{

    public LCMSPeakDIAMS2 DIAWindow;
    HashMap<Integer, Integer> ScanClusterMap_Q1;
    HashMap<Integer, Integer> ScanClusterMap_Q2;
    HashMap<Integer, String> ScanClusterMap_Q3;
    String Q1Name;
    String Q2Name;    
    String Q3Name;
    LCMSPeakMS1 ms1lcms;
    LCMSID IDsummary;
    int NoThread=1;
    
    public DIA_window_Quant(String Q1Name,String Q2Name,String Q3Name,HashMap<Integer, Integer> ScanClusterMap_Q1, HashMap<Integer, Integer> ScanClusterMap_Q2, HashMap<Integer, String> ScanClusterMap_Q3, LCMSPeakMS1 ms1lcms,LCMSPeakDIAMS2 DIAWindow,LCMSID IDsummary, int NoThreads){
        this.ScanClusterMap_Q1=ScanClusterMap_Q1;
        this.ScanClusterMap_Q2=ScanClusterMap_Q2;
        this.ScanClusterMap_Q3=ScanClusterMap_Q3;
        this.Q1Name=Q1Name;
        this.Q2Name=Q2Name;
        this.Q3Name=Q3Name;
        this.ms1lcms=ms1lcms;
        this.DIAWindow=DIAWindow;
        this.IDsummary=IDsummary;
        this.NoThread=NoThreads;
    }
    @Override
    public void run() {
        
        //System.out.println("Assigning clusters for peak groups in MS2 isolation window:" + FilenameUtils.getBaseName(DIAWindow.ScanCollectionName));
        //Restore serialization files for old projects
        //DIAWindow.ReadPeakFromDB();
        //DIAWindow.WritePeakClusterSerialization();
        //DIAWindow.WritePrecursorFragmentGrouping();
        ////////////////////////////////////
        //Needed to reactivate
       
        if(!DIAWindow.ReadPeakCluster()){
            Logger.getRootLogger().error("Reading Peak cluster result for " + DIAWindow.ScanCollectionName + " failed");
            return;
        }
        ExecutorService executorPool;
        executorPool = Executors.newFixedThreadPool(NoThread);
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (DIAWindow.DIA_MZ_Range.getX() <= pepIonID.ObservedMz && DIAWindow.DIA_MZ_Range.getY() >= pepIonID.ObservedMz) {
                DIAMapClusterUnit mapunit = new DIAMapClusterUnit(pepIonID, Q1Name, Q2Name, Q3Name, ScanClusterMap_Q1, ScanClusterMap_Q2, ScanClusterMap_Q3, ms1lcms, DIAWindow);
                executorPool.execute(mapunit);
            }
        }
        executorPool.shutdown();
        while (!executorPool.isTerminated()) {
        }
        
         if (DIAWindow.datattype != SpectralDataType.DataType.pSMART) {
            if (!DIAWindow.ReadPrecursorFragmentClu2Cur()) {
                Logger.getRootLogger().error("Reading precursor-fragment results for " + DIAWindow.ScanCollectionName + " failed");
                return;
            }
            
             for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
                 if (DIAWindow.DIA_MZ_Range.getX() <= pepIonID.ObservedMz && DIAWindow.DIA_MZ_Range.getY() >= pepIonID.ObservedMz) {
                     for (PeakCluster cluster : pepIonID.MS1PeakClusters) {
                         DIAWindow.ExtractFragmentForPeakCluser(cluster);
                     }
                     for (PeakCluster Assigncluster : pepIonID.MS2UnfragPeakClusters) {
                         PeakCluster cluster = DIAWindow.PeakClusters.get(Assigncluster.Index - 1);
                         if (cluster.TargetMz() == Assigncluster.TargetMz() || cluster.Charge == Assigncluster.Charge) {
                             DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
                         }
                     }
                 }
             }
        }
        DIAWindow.ClearAllPeaks();
    }
    
}
