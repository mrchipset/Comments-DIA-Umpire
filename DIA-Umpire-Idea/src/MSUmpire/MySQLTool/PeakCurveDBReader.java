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
import MSUmpire.PeakDataStructure.PeakCurve;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurveDBReader {

    private ConnectionManager connectionManager = null;

    public PeakCurveDBReader(ConnectionManager connectionManager) {
        this.connectionManager = connectionManager;
    }

    public PeakCurve QueryMS2(String Filename, int PeakCurveIndex, InstrumentParameter parameter, float startMZ, float endMZ, boolean CloseConnection) throws SQLException {
        Statement state = connectionManager.GetConnection().createStatement();
        ResultSet rsCurve = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(Filename) + "_MS2PeakCurve WHERE Curve_index='" + PeakCurveIndex + "' AND StartMZ=" + startMZ + " AND EndMZ=" + endMZ);

        rsCurve.next();
        PeakCurve peakCurve = new PeakCurve(parameter);
        peakCurve.Index = PeakCurveIndex;
        peakCurve.TargetMz = rsCurve.getFloat("mz");
        peakCurve.ReadPeakResultMySQL(rsCurve);
//            String PSMString=rsCurve.getString("PSMs");
//                if (PSMString != "") {
//                    String[] ScanNOs = PSMString.split("_");
//                    try {
//                        for (int i = 0; i < ScanNOs.length; i++) {
//                            peakCurve.PSMs.add(summary.PSMList.get(Integer.parseInt(ScanNOs[i])));
//                        }
//                    } catch (Exception e) {
//                        System.out.print("Get Corr error\n");
//                    }
//                }
        state.close();
        if (CloseConnection) {
            connectionManager.CloseConnection();
        }
        return peakCurve;
    }

    public PeakCurve Query(String Filename, int PeakCurveIndex, InstrumentParameter parameter, boolean CloseConnection) throws SQLException {
        Statement state = connectionManager.GetConnection().createStatement();
        ResultSet rsCurve = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(Filename) + "_PeakCurve WHERE Curve_index='" + PeakCurveIndex + "'");

        rsCurve.next();
        PeakCurve peakCurve = new PeakCurve(parameter);
        peakCurve.Index = PeakCurveIndex;
        peakCurve.TargetMz = rsCurve.getFloat("mz");
        peakCurve.ReadPeakResultMySQL(rsCurve);
//            String PSMString=rsCurve.getString("PSMs");
//                if (PSMString != "") {
//                    String[] ScanNOs = PSMString.split("_");
//                    try {
//                        for (int i = 0; i < ScanNOs.length; i++) {
//                            peakCurve.PSMs.add(summary.PSMList.get(Integer.parseInt(ScanNOs[i])));
//                        }
//                    } catch (Exception e) {
//                        System.out.print("Get Corr error\n");
//                    }
//                }
        state.close();
        if (CloseConnection) {
            connectionManager.CloseConnection();
        }
        return peakCurve;
    }
}
