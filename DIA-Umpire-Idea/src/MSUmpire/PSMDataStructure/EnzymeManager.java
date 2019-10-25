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
package MSUmpire.PSMDataStructure;

import com.compomics.util.experiment.biology.Enzyme;
import com.compomics.util.experiment.biology.EnzymeFactory;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class EnzymeManager {

    private static EnzymeManager enzymeManager;
    EnzymeFactory enzymeFactory = null;

    public EnzymeManager() throws XmlPullParserException, IOException {
        if (enzymeFactory == null) {
            enzymeFactory = EnzymeFactory.getInstance();
            InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/enzymes.xml");
            String tmpfile = "enzymes.xml";
            InputStreamToFile convert = new InputStreamToFile();

            File enzymeFile = convert.GetFile(is, tmpfile);
            enzymeFactory.importEnzymes(enzymeFile);
        }
    }

    public static EnzymeManager GetInstance() throws XmlPullParserException, IOException {
        if (enzymeManager == null) {
            enzymeManager = new EnzymeManager();
        }
        return enzymeManager;
    }

    public Enzyme GetTrypsin() {
        return enzymeFactory.getEnzyme("Trypsin");
    }
}
