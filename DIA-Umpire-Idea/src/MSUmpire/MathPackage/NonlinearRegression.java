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
package MSUmpire.MathPackage;

import MSUmpire.BaseDataStructure.XYPointCollection;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class NonlinearRegression extends Regression {

    private float EstimationPt;

    public NonlinearRegression(float EstimationPt) {
        this.EstimationPt = EstimationPt;
    }

    @Override
    public float GetR2() {
        ComputeSST();
        ComputeSSR();
        equation.R2value = (SST - SSR) / SST;
        return equation.R2value;
    }

    @Override
    public void SetData(XYPointCollection pointset) {
        this.pointset = pointset;
        equation = new Equation();
        FindEquation();
        GeneratePredictYList();
        GeneratePredictXList();
    }

    private void ComputeSSR() {
        SSR = 0;
        for (int i = 0; i < pointset.PointCount(); i++) {
            SSR += (pointset.Data.get(i).getY() - (GetPredictYByTimelist(pointset.Data.get(i).getX()))) * (pointset.Data.get(i).getY() - (GetPredictYByTimelist(pointset.Data.get(i).getX())));
        }
    }
    public XYPointCollection PredictYList;
    public XYPointCollection PredictXList;

    public float GetPredictYByTimelist(float xvalue) {
        return PredictYList.GetPoinByXCloset(xvalue).getY();
    }

    public float GetPredictXByTimelist(float yvalue) {
        return PredictXList.GetPoinByXCloset(yvalue).getY();
    }

    private float PredictYByLocalRegion(float Xvalue) {
        float error = 0;
        int count = 0;

        for (int i = 0; i < pointset.PointCount(); i++) {
            if (pointset.Data.get(i).getX() > Xvalue - EstimationPt && pointset.Data.get(i).getX() < Xvalue + EstimationPt) {
                error += pointset.Data.get(i).getY() - (float) (GetY(pointset.Data.get(i).getX()));
                count++;
            }
        }
        if (count == 0) {
            for (int i = 0; i < pointset.PointCount(); i++) {
                if (pointset.Data.get(i).getX() > Xvalue - 2 * EstimationPt && pointset.Data.get(i).getX() < Xvalue + 2 * EstimationPt) {
                    error += pointset.Data.get(i).getY() - GetY(pointset.Data.get(i).getX());
                    count++;
                }
            }
        }
        if (count > 0) {
            error /= count;
        }

        return GetY(Xvalue) + error;
    }

    private float PredictXByLocalRegion(float Yvalue) {
        float error = 0;
        int count = 0;

        for (int i = 0; i < pointset.PointCount(); i++) {
            if (pointset.Data.get(i).getY() > Yvalue - EstimationPt && pointset.Data.get(i).getY() < Yvalue + EstimationPt) {
                error += pointset.Data.get(i).getX() - (float) (GetX(pointset.Data.get(i).getY()));
                count++;
            }
        }
        if (count == 0) {
            for (int i = 0; i < pointset.PointCount(); i++) {
                if (pointset.Data.get(i).getY() > Yvalue - 2 * EstimationPt && pointset.Data.get(i).getY() < Yvalue + 2 * EstimationPt) {
                    error += pointset.Data.get(i).getX() - GetX(pointset.Data.get(i).getY());
                    count++;
                }
            }
        }
        if (count > 0) {
            error /= count;
        }
        return GetX(Yvalue) + error;
    }

    public void GeneratePredictYList() {
        PredictYList = new XYPointCollection();

        float gap = 0.2f;
        for (int i = 0; i < (int) (max_x + 2 - min_x) * 5; i++) {
            float x = min_x - 1 + gap * i;
            PredictYList.AddPoint(x, PredictYByLocalRegion(x));
        }
    }

    public void GeneratePredictXList() {
        PredictXList = new XYPointCollection();
        float gap = 0.2f;
        for (int i = 0; i < (int) (max_y + 2 - min_y) * 5; i++) {
            float y = min_y - 1 + gap * i;
            PredictXList.AddPoint(y, PredictXByLocalRegion(y));
        }
    }
}
