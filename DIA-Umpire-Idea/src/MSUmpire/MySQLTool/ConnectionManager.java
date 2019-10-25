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

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ConnectionManager {

    private Connection connection;

    public ConnectionManager(String URL, int Port, String DBName, String UserName, String PWD) {

        this.URLString = URL;
        this.Port = Port;
        this.DBname = DBName;
        this.UserName = UserName;
        this.PWD = PWD;
        try {
            Connection tempConnection = DriverManager.getConnection("jdbc:mysql://" + URLString + ":" + Port, UserName, PWD);
            tempConnection.createStatement().execute("CREATE DATABASE IF NOT EXISTS " + DBname);
            tempConnection.close();
        } catch (SQLException ex) {
            Logger.getLogger(ConnectionManager.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public Connection GetConnection() throws SQLException {
        if (connection == null || connection.isClosed()) {
            connection = DriverManager.getConnection("jdbc:mysql://" + URLString + ":" + Port + "/" + DBname, UserName, PWD);
        }
        return connection;
    }
    private String URLString;
    private int Port;
    private String UserName;
    private String PWD;
    private String DBname;

    public void CloseConnection() throws SQLException {
        if (connection != null) {
            connection.close();
            connection = null;
        }
    }
}
