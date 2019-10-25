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

import java.io.Serializable;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MatchSubscore implements Serializable{
    private static final long serialVersionUID = -6347236539892953039L;    
    public subscore[] Subscores;    
    public int NoEnableSubscores;
    public double[] SubSCoeff;
    public String[] SubSName;
        
    public class subscore implements Serializable{
        private static final long serialVersionUID = 2801998695766081693L;
        public String name;
        public boolean enable;
        public double initialweight;

        public subscore(String name, boolean enable,double iniweight) {
            this.name = name;
            this.enable = enable;
            this.initialweight=iniweight;
        }
    }
    public MatchSubscore() {
        Subscores = new subscore[13];
        Subscores[0] = new subscore("SpecDotProduct", true, 6d);
        Subscores[1] = new subscore("SpecCorrelation", false, 6d);
        Subscores[2] = new subscore("SpecContrastAngle", false, 6d);
        Subscores[3] = new subscore("CorrScore", true, 1d);
        Subscores[4] = new subscore("IntScore", true, 0.1d);
        Subscores[5] = new subscore("PPMScore", true, 4d);
        Subscores[6] = new subscore("ApexDeltaScore", true, 4d);
        Subscores[7] = new subscore("RTOverlapScore", false, 1d);
        Subscores[8] = new subscore("NoMatchB", true, 1d);
        Subscores[9] = new subscore("NoMatchY", true, 1d);
        Subscores[10] = new subscore("MS1Corr", false, 1.5d);
        Subscores[11] = new subscore("RTDiff", false, 1d);
        Subscores[12] = new subscore("PrecursorPPM", true, 1d);
    }
    
    public void InitializeCoeff() {
        NoEnableSubscores = 0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                NoEnableSubscores++;
            }
        }
        SubSCoeff=new double[NoEnableSubscores];
        SubSName=new String[NoEnableSubscores];
        int idx=0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                SubSCoeff[idx] = Subscores[i].initialweight;
                SubSName[idx]=Subscores[i].name;
                idx++;
            }
        }
    }
    public void AssignUmpireScore(PeakGroupScore peakgroup){
        float[] Subs=GetSubScoreArray(peakgroup);
        peakgroup.UmpireScore = 0f;
        int idx=0;
        for (int i = 0; i < Subs.length; i++) {
            if (Subscores[i].enable) {
                peakgroup.UmpireScore += (float) (SubSCoeff[idx++] * Subs[i]);
            }
        }
    }
    
    private float [] GetSubScoreArray(PeakGroupScore peakgroup){
        float[] Subs = new float[13];
        int idx=0;
        Subs[idx++] = peakgroup.SpecDotProduct;
        Subs[idx++] = peakgroup.SpecCorrelation;
        Subs[idx++] = peakgroup.ContrastAngle;
        Subs[idx++] = peakgroup.CorrScore;
        Subs[idx++] = peakgroup.IntScore;
        Subs[idx++] = peakgroup.PPMScore;
        Subs[idx++] = peakgroup.ApexDeltaScore;
        Subs[idx++] = peakgroup.RTOverlapScore;
        Subs[idx++] = peakgroup.NoMatchB;
        Subs[idx++] = peakgroup.NoMatchY;
        Subs[idx++] = peakgroup.MS1Corr;
        Subs[idx++] = peakgroup.RTDiff;
        Subs[idx++] = peakgroup.PrecursorPPM;
        return Subs;
    }
    
    public double[] GetEnableSubScoreArray(PeakGroupScore peakgroup){
        float[] Subs=GetSubScoreArray(peakgroup);
        double[] enablesubsore=new double[NoEnableSubscores];
        int idx=0;
        for (int i = 0; i < Subs.length; i++) {
            if (Subscores[i].enable) {
                enablesubsore[idx++] = Subs[i];
            }
        }    
        return enablesubsore;
    }
}
