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
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LinearInterpolation {

    public XYPointCollection Run(XYPointCollection data, int PtNum) {

        XYData[] Smoothdata = new XYData[PtNum];
        float intv = (data.Data.get(data.PointCount() - 1).getX() - data.Data.get(0).getX()) / (float) PtNum;
        float rt = data.Data.get(0).getX();
        for (int i = 0; i < PtNum; i++) {
            Smoothdata[i] = new XYData(intv * i + rt, -1f);
        }
        int index = 0;
        for (XYData point : data.Data) {
            XYData closet = Smoothdata[index];
            boolean found = false;
            for (int i = index; i < PtNum - 1; i++) {
                if (Smoothdata[i].getX() <= point.getX() && Smoothdata[i + 1].getX() > point.getX()) {
                    Smoothdata[i].setY(point.getY());
                    index = i;
                    found = true;
                    break;
                }
            }
            if (!found) {
                Smoothdata[PtNum - 1].setY(point.getY());
                index = PtNum - 1;
            }
        }

        boolean gapfound = false;
        int startidx = 0;
        int endidx = 0;
        float startintensity = Smoothdata[0].getY();
        float endintensity = Smoothdata[0].getY();

        for (int i = 1; i < PtNum; i++) {
            if (gapfound && Smoothdata[i].getY() != -1) {
                endidx = i;
                endintensity = Smoothdata[i].getY();
                Smoothdata[(startidx + endidx) / 2].setY((startintensity + endintensity) / 2);
                i = startidx;
                gapfound = false;
            }
            if (!gapfound && Smoothdata[i].getY() == -1) {
                startidx = i - 1;
                startintensity = Smoothdata[i - 1].getY();
                gapfound = true;
            }
        }
        XYPointCollection returndata = new XYPointCollection();
        for (XYData point : Smoothdata) {
            returndata.AddPoint(point);
        }
        return returndata;
    }

}
