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
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class QuantitationModule {

    private LCMSPeakMS1 lCMSQuant;
    private InstrumentParameter parameter;

    public QuantitationModule(LCMSPeakMS1 lCMSQuant, InstrumentParameter parameter) {
        this.lCMSQuant = lCMSQuant;
        this.parameter = parameter;
    }

    public void Quantify() throws IOException {
        QuantifyPredictedIon();
    }

    private void ExportResult() throws IOException {
        FileWriter writer = new FileWriter("IDQuant.csv");
        for (PepIonID pepIonID : lCMSQuant.IDsummary.GetPepIonList().values()) {
            if (pepIonID.GetRTSD() < 1 && pepIonID.MS1PeakClusters.size() > 0) {
                writer.write(pepIonID.ModSequence + "," + pepIonID.Charge + "," + pepIonID.PeakArea[0] + "," + pepIonID.PeakHeight[0] + "," + pepIonID.MS1PeakClusters.get(0).Index + "\n");
            }
        }
        writer.close();
    }

    private void QuantifyPredictedIon() {
        for (PepIonID pepIonID : lCMSQuant.IDsummary.GetMappedPepIonList().values()) {
            pepIonID.CreateQuantInstance(lCMSQuant.MaxNoPeakCluster);
            float minRT = Float.MAX_VALUE;
            float maxRT = Float.MIN_VALUE;
            for (float rt : pepIonID.PredictRT) {
                if (rt > maxRT) {
                    maxRT = rt;
                }
                if (rt < minRT) {
                    minRT = rt;
                }
            }
            float mz = pepIonID.NeutralPrecursorMz();
            int startindex = lCMSQuant.MZSortedClusters.BinarySearchLower(InstrumentParameter.GetMzByPPM(mz, pepIonID.Charge, parameter.MS1PPM));

            PeakCluster targetCluster = null;
            float maxscore = Float.MIN_VALUE;

            for (int i = startindex; i < lCMSQuant.MZSortedClusters.size(); i++) {
                PeakCluster candidateCluster = lCMSQuant.MZSortedClusters.get(i);
                if (InstrumentParameter.CalcPPM(mz, candidateCluster.TargetMz()) > parameter.MS1PPM) {
                    if (candidateCluster.TargetMz() > mz) {
                        break;
                    }
                } else {
                    float chisq = candidateCluster.GetChiSquareProbByTheoIso(pepIonID.IsotopicDistrubtion(lCMSQuant.MinNoPeakCluster));
                    float score = chisq + candidateCluster.Corrs[0];
                    if (score > maxscore) {
                        targetCluster = candidateCluster;
                        maxscore = score;
                    }
                }
            }
            if (targetCluster != null) {
                pepIonID.PeakArea = targetCluster.PeakArea;
                pepIonID.PeakHeight = targetCluster.PeakHeight;
                pepIonID.MS1PeakClusters.add(targetCluster);
            }
        }
    }
}
