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
import com.vseravno.solna.SolnaParser;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ProtXMLParser {

    public LCMSID SingleLCMSID;
    public String FileName;
    public float threshold = 0f;

    public ProtXMLParser(LCMSID singleLCMSID, String FileName, float threshold) throws ParserConfigurationException, SAXException, IOException, ClassNotFoundException, XmlPullParserException, InterruptedException {
        this.SingleLCMSID = singleLCMSID;
        this.FileName = FileName;
        this.threshold = threshold;
        Logger.getRootLogger().info("Parsing protXML: " + FileName + ".....");
        ParseSAX();
        SingleLCMSID.DetermineAssignIonListByProtPepSeq();
        SingleLCMSID.FixProteinWithDecoyHead();
        SingleLCMSID.SetGroupProbForNonDecoyGroupHead();
        //System.out.print("done\n");
    }


   private void ParseSAX() throws ParserConfigurationException, SAXException, IOException, XmlPullParserException {
        File fXmlFile = new File(FileName);
        if (!fXmlFile.exists()) {
            Logger.getRootLogger().info("File :" + FileName + " doesn't exist\n");
            return;
        }
        FileInputStream inputStream = new FileInputStream(FileName);
        SolnaParser parser = new SolnaParser();
        ProtXMLParseHandler handler = new ProtXMLParseHandler(SingleLCMSID,threshold);
        parser.addHandler("/protein_summary/protein_group", handler);
        parser.parse(inputStream);        
    }    
}
