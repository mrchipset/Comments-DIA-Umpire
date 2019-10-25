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

import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.LCMSBaseStructure.LCMSPeakBase;
import java.io.*;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerMS1 extends PDHandlerBase {

    public PDHandlerMS1(LCMSPeakBase lcmspeak, int NoCPUs, float PPM) throws SQLException {
        this.NoCPUs = NoCPUs;
        this.PPM = PPM;
        this.LCMSPeakBase = lcmspeak;
        this.parameter = lcmspeak.parameter;
        this.connectionManager = lcmspeak.connectionManager;
    }

    public void DetectPeakCurves(ArrayList<ScanCollection> scanCollections) throws InterruptedException, ExecutionException, IOException, SQLException {        
        ReadPepIsoMS1PatternMap();
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        for (ScanCollection scanCollection : scanCollections) {
            FindAllPeakCurve(scanCollection);
        }
        WaveletDetectMax();
        ClearRawPeaks();
        PeakCurveCorrClustering_V2(new XYData(Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY));
    }
    
    private void ReadPepIsoMS1PatternMap() throws FileNotFoundException, IOException {

        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/IsotopicPatternRange.csv");
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        IsotopePatternMap = new TreeMap[LCMSPeakBase.MaxNoPeakCluster - 1];
        for (int i = 0; i < IsotopePatternMap.length; i++) {
            IsotopePatternMap[i] = new TreeMap<>();
        }
        String line = "";
        while ((line = reader.readLine()) != null) {
            float MW = Float.parseFloat(line.split(",")[0]);

            for (int i = 0; i < LCMSPeakBase.MaxNoPeakCluster - 1; i++) {
                float Mean = Float.parseFloat(line.split(",")[1 + (i * 2)]);
                float SD = Float.parseFloat(line.split(",")[2 + (i * 2)]);

                if (!Float.isNaN(Mean)) {
                    //IsoMapSecond.put(MW,new Normal(MeanSecond, SDSecond));
                    //IsoMapThird.put(MW,new Normal(MeanThird, SDThird));                    
                    IsotopePatternMap[i].put(MW, new XYData(Mean + 3.3f * SD, Mean - 3.3f * SD));
                }
            }
        }
        reader.close();
    }

}
