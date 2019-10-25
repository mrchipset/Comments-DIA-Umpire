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
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ModStringConvert {

    //n[43]AAAAAAGAGPEM[147]VR

    public static void ConvertTPPModString(String modString, ArrayList<ModificationMatch> Modifications) throws XmlPullParserException, IOException {
        String Sequence = modString.replace("n[", "[").replace("c[", "").replaceAll("[\\[0-9\\]]", "");
        String ConvertString = Sequence;
        while (modString.contains("[")) {
            String site = String.valueOf(modString.charAt(modString.indexOf("[") - 1));
            int idx = -1;
            float massdiff = 0f;
            if (site.equals("n")) {
                site = "N-term";
                idx = 1;
                massdiff = (float) (Float.parseFloat(modString.substring(modString.indexOf("[")+1, modString.indexOf("]"))) - ElementaryIon.proton.getTheoreticMass());
                modString = modString.replaceFirst("[n\\[0-9\\]]", "");
            } else if (site.equals("c")) {
                site = "C-term";
                idx = Sequence.length();
                massdiff = (float) (Float.parseFloat(modString.substring(modString.indexOf("["), modString.indexOf("]"))) - ElementaryIon.proton.getTheoreticMass());
                modString = modString.replaceFirst("[c\\[0-9\\]]", "");
            } else {
                idx = modString.indexOf("[") - 1;
                site = String.valueOf(modString.charAt(idx - 1));
                AminoAcid aa = AminoAcid.getAminoAcid(site.charAt(0));
                massdiff = (float) (Float.parseFloat(modString.substring(modString.indexOf("["), modString.indexOf("]"))) - aa.monoisotopicMass);
                modString = modString.replaceFirst("[\\[0-9\\]]", "");
            }

            PTM ptm = PTMManager.GetInstance().GetPTM(site, massdiff);
            Modifications.add(new ModificationMatch(ptm.getName(), true, idx));
            massdiff = (float) ptm.getMass();
            String tmp = new DecimalFormat("#.###").format(massdiff);
            ConvertString = AddModIntoSeq(ConvertString, "[" + tmp + "(" + site + ")]", idx - 1);
        }
    }
    
    
    

    public static String AddModIntoSeq(String seq, String modstring, int index) {
        boolean inmod = false;

        if (index == -1) {
            return modstring + seq;
        }
        int countidx = 0;
        for (int i = 0; i < seq.length(); i++) {
            if (seq.charAt(i) == '[') {
                inmod = true;
            } else if (seq.charAt(i) == ']') {
                inmod = false;
            } else if (!inmod) {
                if (countidx == index) {
                    return seq.substring(0, i) + modstring + seq.substring(i);
                }
                countidx++;
            }
        }
        return seq;
    }
}
