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
public class PearsonCorr {

    public float CalcCorrNeighborBin(XYPointCollection CollectionA, XYPointCollection CollectionB) {
        Regression regression = new Regression();

        int num = (int) ((Math.min(CollectionA.Data.get(CollectionA.PointCount() - 1).getX(), CollectionB.Data.get(CollectionB.PointCount() - 1).getX()) - Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX())) * 100);
        float timeinterval = 1 / 100f;

        float[] arrayA = new float[num];
        float[] arrayB = new float[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        for (int i = 0; i < num - 1; i++) {
            float low = start + i * timeinterval;
            float up = start + (i + 1) * timeinterval;

            for (int j = 0; j < CollectionA.PointCount(); j++) {
                if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                    float intenlow = CollectionA.Data.get(j).getY() * (1 - (CollectionA.Data.get(j).getX() - low) / timeinterval);
                    float intenup = CollectionA.Data.get(j).getY() * (1 - (up - CollectionA.Data.get(j).getX()) / timeinterval);
                    if (intenlow > arrayA[i]) {
                        arrayA[i] = intenlow;
                    }
                    if (intenup > arrayA[i + 1]) {
                        arrayA[i + 1] = intenup;
                    }
                } else if (CollectionA.Data.get(j).getX() > up) {
                    break;
                }
            }

            for (int j = 0; j < CollectionB.PointCount(); j++) {
                if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                    float intenlow = CollectionB.Data.get(j).getY() * (1 - (CollectionB.Data.get(j).getX() - low) / timeinterval);
                    float intenup = CollectionB.Data.get(j).getY() * (1 - (up - CollectionB.Data.get(j).getX()) / timeinterval);
                    if (intenlow > arrayB[i]) {
                        arrayB[i] = intenlow;
                    }
                    if (intenup > arrayB[i + 1]) {
                        arrayB[i + 1] = intenup;
                    }
                } else if (CollectionB.Data.get(j).getX() > up) {
                    break;
                }
            }
        }

        XYPointCollection pointset = new XYPointCollection();
        for (int i = 0; i < num; i++) {
            if (arrayA[i] > 0 && arrayB[i] > 0) {
                pointset.AddPoint(arrayA[i], arrayB[i]);
            }
        }

        float R2 = 0f;

        if (pointset.PointCount() > 5) {
            regression.SetData(pointset);
            if (regression.equation.Mvalue > 0) {
                R2 = regression.GetR2();
            }
        }
        pointset.dispose();
        regression.dispose();
        regression = null;
        arrayA = null;
        arrayB = null;
        return R2;
    }

    public float CalcCorr(XYPointCollection CollectionA, XYPointCollection CollectionB, int NoPointPerInterval) {
        Regression regression = new Regression();

        int num = Math.max(CollectionA.PointCount(), CollectionB.PointCount()) / 2;
        float timeinterval = 2f / (float) NoPointPerInterval;
        if (num < 1) {
            return 0f;
        }

//        double[] arrayA = new double[num];
//        double[] arrayB = new double[num];
        float[] arrayA = new float[num];
        float[] arrayB = new float[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        int i = 0;
        float low = start;
        float up = start + timeinterval;

        for (int j = 0; j < CollectionA.PointCount(); j++) {
            while (CollectionA.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                if (CollectionA.Data.get(j).getY() > arrayA[i]) {
                    arrayA[i] = CollectionA.Data.get(j).getY();
                }
            }
        }
        i = 0;
        low = start;
        up = start + timeinterval;
        for (int j = 0; j < CollectionB.PointCount(); j++) {
            while (CollectionB.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                if (CollectionB.Data.get(j).getY() > arrayB[i]) {
                    arrayB[i] = CollectionB.Data.get(j).getY();
                }
            }
        }

        for (int idx = 1; idx < num - 1; idx++) {
            if (arrayA[idx] == 0f) {
                arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
            }
            if (arrayB[idx] == 0f) {
                arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
            }
        }

        XYPointCollection pointset = new XYPointCollection();
        for (int idx = 0; idx < num; idx++) {
            if (arrayA[idx] > 0 && arrayB[idx] > 0) {
                pointset.AddPoint(arrayA[idx], arrayB[idx]);
            }
        }

        float R2 = 0f;
        if (pointset.PointCount() > 5) {
            regression.SetData(pointset);
            if (regression.equation.Mvalue > 0) {
                R2 = regression.GetR2();
            }
        }
        //PearsonsCorrelation correlation=new PearsonsCorrelation();
        //double R2_2=correlation.correlation(arrayA, arrayB);

        pointset.dispose();
        regression.dispose();
        regression = null;
        arrayA = null;
        arrayB = null;
        return R2;
    }
}
