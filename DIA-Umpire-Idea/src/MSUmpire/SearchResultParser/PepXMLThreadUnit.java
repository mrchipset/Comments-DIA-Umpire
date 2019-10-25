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
package MSUmpire.SearchResultParser;

import MSUmpire.PSMDataStructure.LCMSID;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepXMLThreadUnit implements Runnable {

    public LCMSID lCMSPepID;
    String mzxml;
    String pepxml;
    float pepscore;

    public PepXMLThreadUnit(String mzxml, String pepxml, float pepscore) {
        this.mzxml = mzxml;
        this.pepxml = pepxml;
        this.pepscore = pepscore;
    }

    @Override
    public void run() {
        lCMSPepID = new LCMSID(mzxml);
        try {
            PepXMLParser parser = new PepXMLParser(lCMSPepID, pepxml, pepscore);
        } catch (Exception ex) {
           Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
}
