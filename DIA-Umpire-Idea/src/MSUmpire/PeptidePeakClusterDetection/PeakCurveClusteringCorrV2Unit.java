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
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedCurveCollectionMZ;
import MSUmpire.PeakDataStructure.SortedCurveCollectionApexRT;
import Utility.UpdateProcess;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurveClusteringCorrV2Unit implements Runnable {

    PeakCurve peakA;
    InstrumentParameter parameter;
    private SortedCurveCollectionMZ PeakCurveListMZ;
    private SortedCurveCollectionApexRT PeakCurveSortedListApexRT;
    private TreeMap<Float, XYData>[] IsotopePatternMap;
    public ArrayList<PeakCluster> ResultClusters = new ArrayList<>();

    private int MaxNoOfClusters;
    private int MinNoOfClusters;
    private int StartCharge;
    private int EndCharge;
    UpdateProcess update;

    public PeakCurveClusteringCorrV2Unit(PeakCurve targeCurve, SortedCurveCollectionMZ PeakCurveListMZ, SortedCurveCollectionApexRT PeakCurveSortedListApexRT, InstrumentParameter parameter, TreeMap<Float, XYData>[] IsotopePatternMap, int StartCharge, int EndCharge, int MaxNoClusters, int MinNoClusters, UpdateProcess update) {
        this.peakA = targeCurve;
        this.PeakCurveListMZ = PeakCurveListMZ;
        this.PeakCurveSortedListApexRT = PeakCurveSortedListApexRT;
        this.parameter = parameter;
        this.IsotopePatternMap = IsotopePatternMap;
        this.MaxNoOfClusters = MaxNoClusters;
        this.MinNoOfClusters = MinNoClusters;
        this.StartCharge = StartCharge;
        this.EndCharge = EndCharge;
        this.update = update;
    }

    @Override
    public void run() {

        int startRTidx = PeakCurveSortedListApexRT.BinarySearchLower(peakA.ApexRT - parameter.ApexDelta);
        int endRTidx = PeakCurveSortedListApexRT.BinarySearchHigher(peakA.ApexRT+ parameter.ApexDelta);
        
        HashSet<Integer> comparelist = GetHash(PeakCurveSortedListApexRT,startRTidx,endRTidx, peakA.Index);
        float Arange = peakA.EndRT() - peakA.StartRT();
        for (int charge = EndCharge; charge >= StartCharge; charge--) {
            PeakCluster peakCluster = new PeakCluster(MaxNoOfClusters, charge);
            peakCluster.IsoPeaksCurves[0] = peakA;
            peakCluster.MonoIsotopePeak=peakA;
            XYData[] Ranges = new XYData[MaxNoOfClusters - 1];
            for (int i = 0; i < MaxNoOfClusters - 1; i++) {
                Entry range = IsotopePatternMap[i].ceilingEntry(peakCluster.NeutralMass());
                if (range == null) {
                    range = IsotopePatternMap[i].lastEntry();
                }
                Ranges[i] = (XYData) range.getValue();
            }

            for (int pkidx = 1; pkidx < MaxNoOfClusters; pkidx++) {
                float ppmthreshold = parameter.MS1PPM + (parameter.MS1PPM * pkidx * 0.5f);
                float lowtheomz = InstrumentParameter.GetMzByPPM(peakA.TargetMz + (pkidx * (1f / charge)), charge, ppmthreshold);
                float uptheomz = InstrumentParameter.GetMzByPPM(peakA.TargetMz + (pkidx * (1f / charge)), charge, -ppmthreshold);
                int startmzidx = PeakCurveListMZ.BinarySearchLower(lowtheomz);

                float theomz = peakA.TargetMz + (pkidx * (1f / charge));
                float maxscore = 0f;
                float maxcorr = 0f;
                float maxoverlap = 0f;
                PeakCurve closestPeak = null;

                for (int mzidx = startmzidx; mzidx < PeakCurveListMZ.size(); mzidx++) {
                    PeakCurve peakB = PeakCurveListMZ.get(mzidx);

                    if (peakB.Index <= peakA.Index) {
                        continue;
                    }
                    if (peakB.TargetMz > uptheomz) {
                        break;
                    }

                    if (comparelist.contains(peakB.Index)) {
                        float Brange = peakB.EndRT() - peakB.StartRT();
                        float OverlapP = 0f;
                        if (peakA.StartRT() >= peakB.StartRT() && peakA.StartRT() <= peakB.FinalEndRT() && peakA.FinalEndRT() >= peakB.FinalEndRT()) {
                            OverlapP = (peakB.FinalEndRT() - peakA.StartRT()) / Brange;

                        } else if (peakA.FinalEndRT() >= peakB.StartRT() && peakA.FinalEndRT() <= peakB.FinalEndRT() && peakA.StartRT() <= peakB.StartRT()) {
                            OverlapP = (peakA.FinalEndRT() - peakB.StartRT()) / Brange;

                        } else if (peakA.StartRT() <= peakB.StartRT() && peakA.FinalEndRT() >= peakB.FinalEndRT()) {
                            OverlapP = 1;

                        } else if (peakA.StartRT() >= peakB.StartRT() && peakA.FinalEndRT() <= peakB.FinalEndRT()) {
                            OverlapP = Arange / Brange;
                        }
                        if (OverlapP > 0.3 && peakA.ApexRT >= peakB.StartRT() && peakA.ApexRT <= peakB.FinalEndRT() && peakB.ApexRT >= peakA.StartRT() && peakB.ApexRT <= peakA.FinalEndRT()) {
                            float ppm = InstrumentParameter.CalcPPM(theomz, peakB.TargetMz);
                            if (ppm < ppmthreshold) {
                                float corr = 0f;
                                try {
                                    corr = PeakCurveCorrCalc.CalPeakCorr(peakA, peakB, parameter.NoPeakPerMin);
                                } catch (IOException ex) {
                                    Logger.getLogger(PeakCurveClusteringCorrV2Unit.class.getName()).log(Level.SEVERE, null, ex);
                                }
                                if (Float.isNaN(corr)) {
                                    corr = 0f;
                                    //System.out.print("Corr=NAN\n");
                                }

                                if (corr > 0f) {
                                    float PeakIntA = peakA.ApexInt;
                                    float PeakIntB = peakB.ApexInt;

                                    float intscore = 0f;
                                    float IntRatio = PeakIntB / PeakIntA;

                                    if (IntRatio > Ranges[pkidx - 1].getY() && IntRatio <= Ranges[pkidx - 1].getX()) {
                                        intscore = 1f;
                                    } else {
                                        if (Math.abs(IntRatio - Ranges[pkidx - 1].getY()) > Math.abs(IntRatio - Ranges[pkidx - 1].getX())) {
                                            intscore = 1 - Math.abs(IntRatio - Ranges[pkidx - 1].getX());
                                        } else {
                                            intscore = 1 - Math.abs(IntRatio - Ranges[pkidx - 1].getY());
                                        }
                                    }
                                    if (intscore < 0f) {
                                        intscore = 0f;
                                    }
                                    float score = ((ppmthreshold - ppm) / ppmthreshold) + corr + intscore;

                                    if (maxscore < score) {
                                        maxscore = score;
                                        closestPeak = peakB;
                                        maxcorr = corr;
                                        maxoverlap = OverlapP;
                                    }
                                }
                            }
                        }
                    }
                }

                if (closestPeak != null) {
                    peakCluster.Corrs[pkidx - 1] = maxcorr;
                    peakCluster.IsoPeaksCurves[pkidx] = closestPeak;
                    peakCluster.OverlapRT[pkidx-1]=maxoverlap;
                    peakCluster.GetSNR(pkidx-1);

                    if (pkidx == 1) {
                        peakCluster.OverlapP = maxoverlap;
                    }
                }
                if (closestPeak == null) {
                    break;
                }
            }
            if (peakCluster.IsotopeComplete(MinNoOfClusters)) {
                peakCluster.CalcPeakArea_V2();
                peakCluster.UpdateIsoMapProb(IsotopePatternMap);
                peakCluster.UpdateIsoMapError(IsotopePatternMap);
                peakCluster.AssignConfilictCorr();
                peakCluster.LeftInt = peakA.GetSmoothedList().Data.get(0).getY();
                peakCluster.RightInt = peakA.GetSmoothedList().Data.get(peakA.GetSmoothedList().PointCount() - 1).getY();
                ResultClusters.add(peakCluster);
                if (peakCluster.Corrs[0] > 0.5f && peakCluster.OverlapP > 0.8f) {
                    for (int i = 1; i < peakCluster.IsoPeaksCurves.length; i++) {
                        PeakCurve peak = peakCluster.IsoPeaksCurves[i];
                        if (peak != null && peakCluster.Corrs[i-1] > 0.5f && peakCluster.OverlapRT[i-1] > 0.8f) {
                            peak.ChargeGrouped.add(charge);
                        }
                    }
                    break;
                }
            }
        }

        comparelist.clear();
        comparelist = null;
        if (update != null) {
            update.Update();
        }
        //System.out.print("....done\n");
    }

    private HashSet<Integer> GetHash(SortedCurveCollectionApexRT list, int startidx, int endidx, int IndexThresholdI) {
        HashSet<Integer> ReturnList = new HashSet<>();
        for (int i = startidx; i <= endidx; i++) {
            if (list.get(i).Index > IndexThresholdI) {
                ReturnList.add(list.get(i).Index);
            }
        }
        return ReturnList;
    }
}
