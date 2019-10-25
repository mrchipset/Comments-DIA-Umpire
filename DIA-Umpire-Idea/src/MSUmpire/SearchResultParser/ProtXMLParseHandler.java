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
import MSUmpire.PSMDataStructure.ModStringConvert;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.ProtID;
import com.vseravno.solna.SolnaHandler;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ProtXMLParseHandler implements SolnaHandler<Element> {

    public ProtXMLParseHandler(LCMSID SingleLCMSID, float threshold) {        
        this.SingleLCMSID=SingleLCMSID;
        this.threshold=threshold;
    }
    float threshold;
    LCMSID SingleLCMSID;

    @Override
    public void handle(Element node) throws Exception {
        switch(node.getNodeName()){
            case "protein_group": {
                ParseProteinGroupNode(node);
                break;
            }
        }
    }
    
    private void ParseProteinGroupNode(Element protgroupnode) throws XmlPullParserException, IOException, ClassNotFoundException, InterruptedException {

        int groupindex = Integer.parseInt(protgroupnode.getAttributes().getNamedItem("group_number").getNodeValue());
        float groupprob = Float.parseFloat(protgroupnode.getAttributes().getNamedItem("probability").getNodeValue());
        int alphabet = 1;
        
        for (int j = 0; j < protgroupnode.getChildNodes().getLength(); j++) {
            Node protnode = protgroupnode.getChildNodes().item(j);
            if ("protein".equals(protnode.getNodeName())) {
                ProtID proid = new ProtID();
                proid.setAccNo(protnode.getAttributes().getNamedItem("protein_name").getNodeValue());
                proid.IndisProteins.add(proid.getAccNo());
                proid.UniProtID = proid.getAccNo();
                proid.Probability = Float.parseFloat(protnode.getAttributes().getNamedItem("probability").getNodeValue());
                proid.GroupProb=groupprob;              
                
                if (proid.Probability < threshold) {
                    continue;
                }
                proid.ProteinGroup = String.valueOf(groupindex) + "_" + (alphabet++);

                for (int k = 0; k < protnode.getChildNodes().getLength(); k++) {
                    Node child = protnode.getChildNodes().item(k);

                    switch (child.getNodeName()) {
                        case ("indistinguishable_protein"): {
                            proid.IndisProteins.add(child.getAttributes().getNamedItem("protein_name").getNodeValue());
                            proid.IndisProtDes.add(child.getChildNodes().item(1).getAttributes().getNamedItem("protein_description").getNodeValue().replace(",", "_"));
                            break;
                        }
                        case ("annotation"): {
                            proid.SetDescription(child.getAttributes().getNamedItem("protein_description").getNodeValue().replace(",", "_"));
                            proid.IndisProtDes.add(child.getAttributes().getNamedItem("protein_description").getNodeValue().replace(",", "_"));
                            break;
                        }
                        case ("parameter"): {
                            if ("prot_length".equals(child.getAttributes().getNamedItem("name").getNodeValue())) {
                                proid.ProteinLength = Integer.parseInt(child.getAttributes().getNamedItem("value").getNodeValue());
                            }
                            break;
                        }
                        case ("peptide"): {
                            //iProphet
                            if (child.getAttributes().getNamedItem("charge") == null || "0".equals(child.getAttributes().getNamedItem("charge").getNodeValue())) {
                                for (int l = 0; l < child.getChildNodes().getLength(); l++) {
                                    Node chidpep = child.getChildNodes().item(l);
                                    if ("indistinguishable_peptide".equals(chidpep.getNodeName())) {
                                        PepIonID pepIonID = new PepIonID();
                                        pepIonID.Sequence = child.getAttributes().getNamedItem("peptide_sequence").getNodeValue();
                                        pepIonID.MaxProbability = Float.parseFloat(child.getAttributes().getNamedItem("initial_probability").getNodeValue());
                                        pepIonID.Is_NonDegenerate = "Y".equals(child.getAttributes().getNamedItem("is_nondegenerate_evidence").getNodeValue());
                                        pepIonID.Weight = Float.parseFloat(child.getAttributes().getNamedItem("weight").getNodeValue());
                                        pepIonID.GroupWeight = Float.parseFloat(child.getAttributes().getNamedItem("group_weight").getNodeValue());
                                        pepIonID.Charge = Integer.parseInt(chidpep.getAttributes().getNamedItem("charge").getNodeValue());
                                        String modseq = pepIonID.Sequence;
                                        for (int n = 0; n < chidpep.getChildNodes().getLength(); n++) {
                                            Node mod = chidpep.getChildNodes().item(n);

                                            if ("modification_info".equals(mod.getNodeName())) {
                                                if (mod.getAttributes().getNamedItem("mod_nterm_mass") != null) {
                                                    float n_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_nterm_mass").getNodeValue());
                                                    String n_term = new DecimalFormat("#.###").format(n_termmass);
                                                    modseq = ModStringConvert.AddModIntoSeq(modseq, "[" + n_term + "(N-term)]", -1);
                                                }
                                                if (mod.getAttributes().getNamedItem("mod_cterm_mass") != null) {
                                                    float c_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_cterm_mass").getNodeValue());
                                                    String c_term = new DecimalFormat("#.###").format(c_termmass);
                                                    modseq = ModStringConvert.AddModIntoSeq(modseq, "[" + c_term + "(C-term)]", pepIonID.Sequence.length() - 1);
                                                }

                                                for (int m = 0; m < mod.getChildNodes().getLength(); m++) {
                                                    Node aa = mod.getChildNodes().item(m);
                                                    if ("mod_aminoacid_mass".equals(aa.getNodeName())) {
                                                        float tmpmass = Float.parseFloat(aa.getAttributes().getNamedItem("mass").getNodeValue());
                                                        String tmp = new DecimalFormat("#.###").format(tmpmass);
                                                        int idx = Integer.parseInt(aa.getAttributes().getNamedItem("position").getNodeValue());
                                                        modseq = ModStringConvert.AddModIntoSeq(modseq, "[" + tmp + "(" + pepIonID.Sequence.charAt(idx - 1) + ")]", idx - 1);
                                                    }
                                                }
                                            }
                                        }
                                        pepIonID.ModSequence = modseq;
                                        if (proid.MaxIniProb < pepIonID.MaxProbability) {
                                            proid.MaxIniProb = pepIonID.MaxProbability;
                                        }
                                        SingleLCMSID.ProtXMLPepIonList.put(pepIonID.GetKey(), pepIonID);
                                        proid.ProtPeptideID.put(pepIonID.GetKey(), pepIonID);
                                        if (!proid.ProtPepSeq.contains(pepIonID.Sequence)) {
                                            proid.ProtPepSeq.add(pepIonID.Sequence);
                                        }
                                        break;
                                    }
                                }
                            } else {
                                PepIonID pepIonID = new PepIonID();
                                pepIonID.Sequence = child.getAttributes().getNamedItem("peptide_sequence").getNodeValue();
                                pepIonID.MaxProbability = Float.parseFloat(child.getAttributes().getNamedItem("initial_probability").getNodeValue());
                                pepIonID.Is_NonDegenerate = "Y".equals(child.getAttributes().getNamedItem("is_nondegenerate_evidence").getNodeValue());
                                pepIonID.Weight = Float.parseFloat(child.getAttributes().getNamedItem("weight").getNodeValue());
                                pepIonID.GroupWeight = Float.parseFloat(child.getAttributes().getNamedItem("group_weight").getNodeValue());
                                pepIonID.Charge = Integer.parseInt(child.getAttributes().getNamedItem("charge").getNodeValue());

                                //pepIonID.NeutralPepMass = Float.parseFloat(child.getAttributes().getNamedItem("calc_neutral_pep_mass").getNodeValue());
                                String modseq = pepIonID.Sequence;

                                for (int l = 0; l < child.getChildNodes().getLength(); l++) {
                                    Node mod = child.getChildNodes().item(l);

                                    if ("modification_info".equals(mod.getNodeName())) {
                                        if (mod.getAttributes().getNamedItem("mod_nterm_mass") != null) {
                                            float n_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_nterm_mass").getNodeValue());
                                            String n_term = new DecimalFormat("#.###").format(n_termmass);
                                            modseq = ModStringConvert.AddModIntoSeq(modseq, "[" + n_term + "(N-term)]", -1);
                                        }
                                        if (mod.getAttributes().getNamedItem("mod_cterm_mass") != null) {
                                            float c_termmass = Float.parseFloat(mod.getAttributes().getNamedItem("mod_cterm_mass").getNodeValue());
                                            String c_term = new DecimalFormat("#.###").format(c_termmass);
                                            modseq = ModStringConvert.AddModIntoSeq(modseq, "[" + c_term + "(C-term)]", pepIonID.Sequence.length() - 1);
                                        }

                                        for (int m = 0; m < mod.getChildNodes().getLength(); m++) {
                                            Node aa = mod.getChildNodes().item(m);
                                            if ("mod_aminoacid_mass".equals(aa.getNodeName())) {
                                                float tmpmass = Float.parseFloat(aa.getAttributes().getNamedItem("mass").getNodeValue());
                                                String tmp = new DecimalFormat("#.###").format(tmpmass);
                                                int idx = Integer.parseInt(aa.getAttributes().getNamedItem("position").getNodeValue());

                                                modseq = ModStringConvert.AddModIntoSeq(modseq, "[" + tmp + "(" + pepIonID.Sequence.charAt(idx - 1) + ")]", idx - 1);
                                            }
                                        }
                                    }
                                }
                                pepIonID.ModSequence = modseq;
                                if (proid.MaxIniProb < pepIonID.MaxProbability) {
                                    proid.MaxIniProb = pepIonID.MaxProbability;
                                }
                                SingleLCMSID.ProtXMLPepIonList.put(pepIonID.GetKey(), pepIonID);
                                proid.ProtPeptideID.put(pepIonID.GetKey(), pepIonID);
                                if (!proid.ProtPepSeq.contains(pepIonID.Sequence)) {
                                    proid.ProtPepSeq.add(pepIonID.Sequence);
                                }
                                break;
                            }
                        }
                    }
                }
                if (proid.ProtPeptideID.size() > 0) {
                    SingleLCMSID.AddProtID(proid);
                    if (!SingleLCMSID.ProteinGroups.containsKey(groupindex)) {
                        SingleLCMSID.ProteinGroups.put(groupindex, new ArrayList<ProtID>());
                    }
                    SingleLCMSID.ProteinGroups.get(groupindex).add(proid);
                }
            }
        }        
    }
    
}
