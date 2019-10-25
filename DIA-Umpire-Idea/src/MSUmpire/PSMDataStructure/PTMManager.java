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

import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.io.SerializationUtils;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PTMManager {

    private PTMFactory ptmFactory;
    private static PTMManager pTMManager;
    private String tempptmfile="ptmFactory-3.28.24.cus";
    private String tmpfilefolder=System.getProperty("user.home") + "/.umpire";

    public static PTMManager GetInstance() throws XmlPullParserException, IOException {
        if (pTMManager == null) {
            pTMManager = new PTMManager();
        }
        return pTMManager;
    }
       
    private PTMManager() throws XmlPullParserException, IOException { 
        
        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/mods.xml");
        String tmpfile = "mods.xml";
        InputStreamToFile convert = new InputStreamToFile();

        File ptmFile = convert.GetFile(is, tmpfile);
        PTMFactory.setSerializationFolder(tmpfilefolder);
        
        PTMFactory.getInstance().clearFactory();
        PTMFactory.getInstance().importModifications(ptmFile, false, true);
        ptmFactory = PTMFactory.getInstance();   
        if (!ptmFactory.getDefaultModifications().isEmpty()) {
            SaveTempFile();
        } else {
            Logger.getRootLogger().error("Modification map file is empty");
        }
//        for (int i = 0; i < ptmFactory.getPTMs().size(); i++) {   
//            String name=ptmFactory.getPTMs().get(i);
//            PTM ptm = ptmFactory.getPTM(name);
//            System.out.println(ptm.getName()+"_"+ptm.getMass());
//        }             
    }

    private void SaveTempFile() throws IOException {
        File factoryFile = new File(tmpfilefolder, tempptmfile);
        if (!factoryFile.getParentFile().exists()) {
            factoryFile.getParentFile().mkdir();
        }
        SerializationUtils.writeObject(ptmFactory, factoryFile);
    }
    
    
    
    public void ImportUserMod(String file) throws XmlPullParserException, IOException{
        File usermod=new File(file);        
        if (usermod.exists()) {
            ptmFactory.importModifications(usermod, true,false);
            if (!ptmFactory.getDefaultModifications().isEmpty()) {
                SaveTempFile();
            }
            else{
                Logger.getRootLogger().error("Modification map file is empty");;
            }
        }        
    }
    
    
    public static ArrayList<ModificationMatch> TranslateModificationString(String ModificationString) {        
        ArrayList<ModificationMatch> modlist = new ArrayList<>();
        if (ModificationString != null && !"".equals(ModificationString)) {
            String[] Mods = ModificationString.split(";");
            for (String mod : Mods) {
                String ptmstring = mod.substring(0, mod.indexOf("("));
                int site = Integer.parseInt(mod.substring(mod.indexOf("(") + 1, mod.indexOf(")")));
                modlist.add(new ModificationMatch(ptmstring, true, site));
            }
        }
        return modlist;
    }
    
    public PTM GetPTM(String AA, float massdiff) {

        double smallmassdiff = Double.MAX_VALUE;
        PTM smallestdiffptm = null;
        for (int i = 0; i < ptmFactory.getPTMs().size(); i++) {
            String name = ptmFactory.getPTMs().get(i);
            PTM ptm = ptmFactory.getPTM(name);
            boolean sitecorrect = false;
            if (("C-term".equals(AA) && name.contains("c-term")) || ("N-term".equals(AA) && name.contains("n-term"))) {
                sitecorrect = true;
            }
            for (Character residue : ptm.getPattern().getAminoAcidsAtTarget()) {
                if (String.valueOf(residue).equals(AA)) {
                    sitecorrect = true;
                }
            }
            if (sitecorrect) {
                double diff = Math.abs(ptm.getMass() - massdiff);
                if (diff < 1f) {
                    if (diff < smallmassdiff) {
                        smallmassdiff = diff;
                        smallestdiffptm = ptm;
                    }
                }
            }
        }

        return smallestdiffptm;
    }
    
    public PTM GetPTMByName(String modname) {

        for (int i = 0; i < ptmFactory.getPTMs().size(); i++) {
            String name = ptmFactory.getPTMs().get(i);
            PTM ptm = ptmFactory.getPTM(name);
            if (ptm.getShortName() == null ? modname == null : ptm.getShortName().equals(modname)) {
                return ptm;
            }
        }
        return null;
    }        
}
