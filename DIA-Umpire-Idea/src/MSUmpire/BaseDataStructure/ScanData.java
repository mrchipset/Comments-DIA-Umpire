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
package MSUmpire.BaseDataStructure;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ScanData extends XYPointCollection{
    public int Num;
    public int MsLevel;
    public float RetentionTime;
    public float StartMz;
    public float EndMz;
    public float BasePeakMz;
    public float BasePeakIntensity;
    private float _totIonCurrent = 0f;
    public float PrecursorMz;
    public int PrecursorCharge;
    public String ActivationMethod;
    public float PrecursorIntensity;
    public String Scantype;
    public int precision;
    public String compressionType;
    public boolean centroided;
    public int precursorScanNum;
    public int PeaksCountString;
    public float background = 0f;
    public String MGFTitle;
    public XYPointCollection TopPeakScan;
    public float windowWideness;
    public String scanType;

    public void Centroiding(int Resolution, float MinMZ) {
        CentroidingbyLocalMaximum(Resolution, MinMZ);
        centroided = true;
    }

    public void GenerateTopPeakScanData(int toppeaks) {
        SortedXYCollectionClass Intsorted = new SortedXYCollectionClass();
        for (int i = 0; i < Data.size; i++) {
            Intsorted.add(new XYData(Data.get(i).getY(), Data.get(i).getX()));
        }
        TopPeakScan = new XYPointCollection();
        for (int i = Intsorted.size() - 1; TopPeakScan.PointCount() < toppeaks && i >= 0; i--) {
            XYData peak = (XYData) Intsorted.get(i);
            TopPeakScan.AddPoint(new XYData(peak.getY(), peak.getX()));
        }
    }

    public void Normalization() {
        if (MaxY != 0) {
            for (int i = 0; i < PointCount(); i++) {
                XYData pt = Data.get(i);
                pt.setY(pt.getY() / MaxY);
            }
        }
    }
   

    public void RemoveSignalBelowBG() {
        SortedXYCollectionClass newData = new SortedXYCollectionClass();
        for (int i = 0; i < Data.size(); i++) {
            if (Data.get(i).getY() > background) {
                newData.add(Data.get(i));
            }
        }
        Data.clear();
        Data = newData;
        Data.Finalize();
    }

    public float TotIonCurrent() {
        if (_totIonCurrent == 0f) {
            for (int i = 0; i < PointCount(); i++) {
                _totIonCurrent += Data.get(i).getY();
            }
        }
        return _totIonCurrent;
    }

    public void SetTotIonCurrent(float totioncurrent) {
        _totIonCurrent = totioncurrent;
    }

    public ScanData CloneScanData() {
        ScanData newscanData = new ScanData();
        for (XYData pt : Data) {
            newscanData.AddPoint(pt);
        }
        newscanData.Num = Num;
        newscanData.MsLevel = MsLevel;
        newscanData.RetentionTime = RetentionTime;
        newscanData.StartMz = StartMz;
        newscanData.EndMz = EndMz;
        newscanData.BasePeakMz = BasePeakMz;
        newscanData.BasePeakIntensity = BasePeakIntensity;
        newscanData.SetTotIonCurrent(_totIonCurrent);
        newscanData.PrecursorMz = PrecursorMz;
        newscanData.PrecursorCharge = PrecursorCharge;
        newscanData.ActivationMethod = ActivationMethod;
        newscanData.PrecursorIntensity = PrecursorIntensity;
        newscanData.Scantype = Scantype;
        newscanData.precision = precision;
        newscanData.compressionType = compressionType;
        newscanData.centroided = centroided;
        return newscanData;
    }

    public ScanData GetNewSubScanBymzRange(float startmz, float endmz) {
        ScanData newScanData = CloneScanData();
        newScanData.Data = GetSubSetByXRange(startmz, endmz).Data;
        return newScanData;
    }
}
