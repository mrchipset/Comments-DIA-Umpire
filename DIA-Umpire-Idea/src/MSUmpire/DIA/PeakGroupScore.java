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

import MSUmpire.PeakDataStructure.PeakCluster;
import java.io.Serializable;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakGroupScore implements Serializable{
    private static final long serialVersionUID = 164894161684L;

    public PeakCluster cluster;
    public int MSlevel = 0;

    public float SpecDotProduct=0f;
    public float SpecCorrelation=0f;
    public float ContrastAngle=0f;
    public float CorrScore=0f;
    public float IntScore=0f;
    public float PPMScore=0f;
    public float ApexDeltaScore=0f;
    public float RTOverlapScore=0f;
    public int NoMatchB = 0;
    public int NoMatchY = 0;
    public float MS1Corr=0f;
    public float RTDiff=0f;
    public float PrecursorPPM=0f;
    
    public float MS1Score=0f;    
    public int NoFragmentLib = 0;
    public float MixtureModelProb;
    public float MixtureModelLocalProb;
    public float UmpireScore;
    
    public PeakGroupScore(PeakCluster cluster) {
        this.cluster = cluster;
    }
    
}
