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
package MSUmpire.MySQLTool;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.PeakDataStructure.PeakCurve;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakClusterDBReader {

    public PeakCluster GetAnnotatedPeakClusterByIndex(ConnectionManager connectionManager, int ClusterIndex, String Filename, InstrumentParameter parameter) throws SQLException {
        Statement state = connectionManager.GetConnection().createStatement();
        ResultSet rs = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(Filename) + "_AnnotatedPeakCluster WHERE Cluster_Index='" + ClusterIndex + "'");

        rs.next();
        int FirstInt = rs.getInt("PeakIdx1");
        int SecondInt = rs.getInt("PeakIdx2");
        int ThirdInt = rs.getInt("PeakIdx3");
        int FourthInt = rs.getInt("PeakIdx4");
        int Charge = rs.getInt("Charge");
        float Corr1 = rs.getFloat("Corr2");
        float Corr2 = rs.getFloat("Corr3");
        float Corr3 = rs.getFloat("Corr4");
        PeakCluster cluster = new PeakCluster(4, Charge);
        cluster.Corrs[0] = Corr1;
        cluster.Corrs[1] = Corr2;
        cluster.Corrs[2] = Corr3;
        PeakCurveDBReader curveDBReader = new PeakCurveDBReader(connectionManager);
        PeakCurve FirstCurve = curveDBReader.Query(Filename, FirstInt, parameter, true);
        PeakCurve SecondCurve = curveDBReader.Query(Filename, SecondInt, parameter, true);
        PeakCurve ThirdCurve = curveDBReader.Query(Filename, ThirdInt, parameter, true);
        if (Corr3 != 0) {
            PeakCurve FourthCurve = curveDBReader.Query(Filename, FourthInt, parameter, true);
            cluster.IsoPeaksCurves[3] = FourthCurve;
        }
        cluster.IsoPeaksCurves[0] = FirstCurve;
        cluster.MonoIsotopePeak = FirstCurve;
        cluster.IsoPeaksCurves[1] = SecondCurve;
        cluster.IsoPeaksCurves[2] = ThirdCurve;

        state.close();
        state = null;
        connectionManager.CloseConnection();
        return cluster;
    }

    public void GetGroupedFragments(PeakCluster cluster, InstrumentParameter parameter, ConnectionManager connectionManager, String SWATHMS2Filename) throws SQLException {
        Statement state = connectionManager.GetConnection().createStatement();
        ResultSet rs = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(SWATHMS2Filename) + "_Clus2Cur WHERE PeakClusterA='" + cluster.Index + "'");

        if (cluster.Fragments == null) {
            cluster.Fragments = new ArrayList<>();
        }
        PeakCurveDBReader curveDBReader = new PeakCurveDBReader(connectionManager);
        while (rs.next()) {
            PrecursorFragmentPairEdge framentClusterUnit = new PrecursorFragmentPairEdge();
            framentClusterUnit.PeakCurveIndexA = rs.getInt("PeakClusterA");
            framentClusterUnit.PeakCurveIndexB = rs.getInt("PeakClusterB");
            framentClusterUnit.Correlation = rs.getFloat("Corr");
            framentClusterUnit.FragmentMz = rs.getFloat("FragmentMz");
            framentClusterUnit.Intensity = rs.getFloat("FragInt");
            framentClusterUnit.ApexDelta = rs.getFloat("ApexDelta");
            framentClusterUnit.RTOverlapP = rs.getFloat("RTOverlapP");

            //groupedFragment.peakCurve = curveDBReader.Query(SWATHMS2Filename, framentClusterUnit.PeakCurveIndexB, parameter, false);
            cluster.Fragments.add(framentClusterUnit);
        }
        connectionManager.CloseConnection();
    }

    public PeakCluster GetPeakClusterByIndex(ConnectionManager connectionManager, int ClusterIndex, String Filename, InstrumentParameter parameter) throws SQLException {
        Statement state = connectionManager.GetConnection().createStatement();
        ResultSet rs = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(Filename) + "_PeakCluster WHERE Cluster_Index='" + ClusterIndex + "'");

        rs.next();
        int FirstInt = rs.getInt("PeakIdx1");
        int SecondInt = rs.getInt("PeakIdx2");
        int ThirdInt = rs.getInt("PeakIdx3");
        //int FourthInt = rs.getInt("PeakIdx4");
        int Charge = rs.getInt("Charge");
        float Corr1 = rs.getFloat("Corr2");
        float Corr2 = rs.getFloat("Corr3");
        //float Corr3 = rs.getFloat("Corr4");
        PeakCluster cluster = new PeakCluster(3, Charge);
        cluster.Corrs[0] = Corr1;
        cluster.Corrs[1] = Corr2;
        //cluster.Corrs[2] = Corr3;

        cluster.PeakHeightRT[0] = rs.getFloat("PeakHeightRT1");
        cluster.PeakHeightRT[1] = rs.getFloat("PeakHeightRT2");
        cluster.PeakHeightRT[2] = rs.getFloat("PeakHeightRT3");

        PeakCurveDBReader curveDBReader = new PeakCurveDBReader(connectionManager);
        PeakCurve FirstCurve = curveDBReader.Query(Filename, FirstInt, parameter, true);
        PeakCurve SecondCurve = curveDBReader.Query(Filename, SecondInt, parameter, true);
        if (ThirdInt != 0) {
            PeakCurve ThirdCurve = curveDBReader.Query(Filename, ThirdInt, parameter, true);
            cluster.IsoPeaksCurves[2] = ThirdCurve;
        }
//        if (Corr3 != 0) {
//            PeakCurve FourthCurve = curveDBReader.Query(Filename, FourthInt, parameter);
//            cluster.IsoPeaksCurves[3] = FourthCurve;
//        }
        cluster.IsoPeaksCurves[0] = FirstCurve;
        cluster.MonoIsotopePeak = FirstCurve;
        cluster.IsoPeaksCurves[1] = SecondCurve;

        state.close();
        state = null;
        return cluster;
    }
}
