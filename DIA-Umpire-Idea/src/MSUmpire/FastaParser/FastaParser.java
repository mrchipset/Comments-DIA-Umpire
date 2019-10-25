/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of computational medicine and bioinformatics, 
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
package MSUmpire.FastaParser;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaParser implements Serializable{
    private static final long serialVersionUID = 19398249L;

    public HashMap<String, String[]> ProteinList;
    
    public FastaParser(String filename){
        ProteinList=new HashMap<>();
        try {
            Parse(filename);
        } catch (IOException ex) {
            Logger.getLogger(FastaParser.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void Parse(String filename) throws FileNotFoundException, IOException {

        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String ACC = "";
        String des="";
        StringBuilder Seq = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (!"".equals(ACC) && !"".equals(Seq.toString())) {
                    ProteinList.put(ACC, new String[]{Seq.toString(),des});
                    Seq=new StringBuilder();
                }
                ACC = line.trim().split(" ")[0].substring(1);
                des=line.replace(">"+ACC+" ", "");
            } else {
                if (!"".equals(line.trim())) {
                    Seq.append(line);
                }
            }
        }
        if (!"".equals(ACC) && !"".equals(Seq.toString())) {
            ProteinList.put(ACC, new String[]{Seq.toString(),des});
        }
        reader.close();
    }
    
}
