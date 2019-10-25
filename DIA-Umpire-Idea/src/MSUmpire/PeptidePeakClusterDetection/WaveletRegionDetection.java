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
import MSUmpire.PeakDataStructure.PeakCurve;
import Utility.UpdateProcess;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class WaveletRegionDetection implements Runnable {

    PeakCurve curve;
    boolean export;
    UpdateProcess update;
    public ArrayList<PeakCurve> ResultCurves;
    InstrumentParameter parameter;

    public WaveletRegionDetection(PeakCurve curve, InstrumentParameter para, UpdateProcess update) {
        this.curve = curve;
        this.parameter = para;
        this.update = update;
    }

    @Override
    public void run() {
        curve.DoBspline();
        //curve.DoInterpolation();                
        curve.DetectPeakRegion();
        ResultCurves = curve.SeparatePeakByRegion(parameter.SNThreshold);
        if (update != null) {
            update.Update();
        }
    }
}
