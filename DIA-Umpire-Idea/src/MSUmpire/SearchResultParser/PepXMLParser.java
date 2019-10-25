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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepXMLParser {

    private LCMSID singleLCMSID;
    public String FileName;
    public float threshold = 0f;
    public float StartRT = 0f;
    public float EndRT = 9999f;
    public boolean FilteredID = false;

    public PepXMLParser(LCMSID singleLCMSID, String FileName, float threshold, float StartRT, float EndRT) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException {
        this.singleLCMSID = singleLCMSID;
        this.FileName = FileName;
        this.threshold = threshold;
        this.StartRT = StartRT;
        this.EndRT = EndRT;
        Logger.getRootLogger().info("Parsing pepXML: " + FileName + "....");
        ParseSAX();
    }

    public PepXMLParser(LCMSID singleLCMSID, String FileName, float threshold) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException {
        this.singleLCMSID = singleLCMSID;
        this.FileName = FileName;
        this.threshold = threshold;
        Logger.getRootLogger().info("Parsing pepXML: " + FileName + "....");
        try {
            ParseSAX();
        } catch (Exception e) {
            Logger.getRootLogger().info("Parsing pepXML: " + FileName + " failed. Trying to fix the file...");
            insert_msms_run_summary(new File(FileName));
            ParseSAX();
        }

        //System.out.print("done\n");
    }

    private void ParseSAX() throws ParserConfigurationException, SAXException, IOException, XmlPullParserException {
        File fXmlFile = new File(FileName);
        if (!fXmlFile.exists()) {
            Logger.getRootLogger().error("File :" + FileName + " doesn't exist");
            return;
        }
        FileInputStream inputStream = new FileInputStream(FileName);
        SolnaParser parser = new SolnaParser();
        PepXMLParseHandler handler = new PepXMLParseHandler(singleLCMSID, StartRT, EndRT, threshold);
        parser.addHandler("/msms_pipeline_analysis/msms_run_summary/spectrum_query", handler);
        parser.addHandler("/msms_pipeline_analysis/msms_run_summary/search_summary", handler);
        parser.parse(inputStream);
    }

    public void insert_msms_run_summary(File inFile) throws FileNotFoundException, IOException {
        // temp file
        File outFile = new File(FileName.replace("pep.xml", "_fixed.pep.xml"));
        FileName = FileName.replace("pep.xml", "_fixed.pep.xml");

        // input
        FileInputStream fis = new FileInputStream(inFile);
        BufferedReader in = new BufferedReader(new InputStreamReader(fis));

        // output         
        FileOutputStream fos = new FileOutputStream(outFile);
        PrintWriter out = new PrintWriter(fos);

        String thisLine = "";
        String lastLine = "";

        while ((thisLine = in.readLine()) != null) {
            if (!"".equals(lastLine)) {
                out.println(lastLine);
            }
            lastLine = thisLine;
        }
        if (!"".equals(lastLine)) {
            out.println("</msms_run_summary>");
            out.println(lastLine);
        }
        out.flush();
        out.close();
        in.close();
    }
}
