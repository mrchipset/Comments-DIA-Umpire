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
import MSUmpire.PSMDataStructure.ModificationInfo;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PTMManager;
import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.vseravno.solna.SolnaHandler;
import java.io.IOException;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepXMLParseHandler implements SolnaHandler<Element> {

    public PepXMLParseHandler(LCMSID singleLCMSID, float StartRT, float EndRT, float threshold) {
        this.StartRT = StartRT;
        this.EndRT = EndRT;
        this.singleLCMSID = singleLCMSID;
        this.threshold = threshold;
    }
    float StartRT;
    float EndRT;
    float threshold;
    LCMSID singleLCMSID;

    @Override
    public void handle(Element node) throws Exception {
        switch (node.getNodeName()) {
            case "spectrum_query": {
                ParseSpectrumNode(node);
                break;
            }
            case "search_summary": {
                ParseSearchSummary(node);
                break;
            }
        }
    }

    private void ParseSearchSummary(Element node) throws XmlPullParserException, XmlPullParserException, IOException {
        if (node.getAttributes().getNamedItem("search_engine") != null) {
            singleLCMSID.SearchEngine = node.getAttributes().getNamedItem("search_engine").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msDetector") != null) {
            singleLCMSID.msDetector = node.getAttributes().getNamedItem("msDetector").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msIonization") != null) {
            singleLCMSID.msIonization = node.getAttributes().getNamedItem("msIonization").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msManufacturer") != null) {
            singleLCMSID.msManufacturer = node.getAttributes().getNamedItem("msManufacturer").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msMassAnalyzer") != null) {
            singleLCMSID.msMassAnalyzer = node.getAttributes().getNamedItem("msMassAnalyzer").getNodeValue();
        }
        if (node.getAttributes().getNamedItem("msModel") != null) {
            singleLCMSID.msModel = node.getAttributes().getNamedItem("msModel").getNodeValue();
        }
        for (int k = 0; k < node.getChildNodes().getLength(); k++) {
            if ("search_database".equals(node.getChildNodes().item(k).getNodeName())) {
                singleLCMSID.DataBase = node.getChildNodes().item(k).getAttributes().getNamedItem("local_path").getNodeValue();
            }
            if ("aminoacid_modification".equals(node.getChildNodes().item(k).getNodeName())) {
                if (node.getChildNodes().item(k).getAttributes().getNamedItem("aminoacid") != null) {
                    String site = node.getChildNodes().item(k).getAttributes().getNamedItem("aminoacid").getNodeValue();
                    float mass = (float) Math.round(Float.parseFloat(node.getChildNodes().item(k).getAttributes().getNamedItem("mass").getNodeValue()) * 1000) / 1000;
                    float massdiff2 = Float.parseFloat(node.getChildNodes().item(k).getAttributes().getNamedItem("massdiff").getNodeValue());
                    AminoAcid aa = AminoAcid.getAminoAcid(site.charAt(0));
                    float massdiff = mass - (float) aa.monoisotopicMass;

                    if (massdiff != 0f && Math.abs(massdiff - massdiff2) < 0.1f) {
                        PTM ptm = PTMManager.GetInstance().GetPTM(site, massdiff);
                        if (ptm == null) {
                            Logger.getRootLogger().warn("Warning! modification in pepxml : amino acid " + site + "(mass diff:" + massdiff + ") doesn't exist in the library.\n");
                        } else {
                            singleLCMSID.AddModification(ptm, site);
                        }
                    } else {
                        if (Math.abs(massdiff2 + 17.0265) < 0.01 && Math.abs(mass - 143.0041f) < 0.001 && "C".equals(site)) {
                            PTM ptm = PTMManager.GetInstance().GetPTM(site, massdiff);
                            if (ptm == null) {
                                Logger.getRootLogger().warn("Warning! modification in pepxml : amino acid " + site + "(mass diff:" + massdiff + ") doesn't exist in the library.\n");
                            } else {
                                singleLCMSID.AddModification(ptm, site);
                            }
                        } else {
                            Logger.getRootLogger().warn("Warning! modification in pepxml : amino acid " + site + "(mass: " + mass + ", massdiff:" + massdiff2 + ") ignored.\n");
                        }
                    }
                }
            }
            if ("terminal_modification".equals(node.getChildNodes().item(k).getNodeName())) {

                if (node.getChildNodes().item(k).getAttributes().getNamedItem("terminus") != null) {
                    String site = "";
                    if ("c".equals(node.getChildNodes().item(k).getAttributes().getNamedItem("terminus").getNodeValue().toLowerCase())) {
                        site = "C-term";
                    }
                    if ("n".equals(node.getChildNodes().item(k).getAttributes().getNamedItem("terminus").getNodeValue().toLowerCase())) {
                        site = "N-term";
                    }
                    float massdiff = Float.parseFloat(node.getChildNodes().item(k).getAttributes().getNamedItem("massdiff").getNodeValue());
                    PTM ptm = PTMManager.GetInstance().GetPTM(site, massdiff);
                    if (ptm == null) {
                        Logger.getRootLogger().warn("Warning! term-modification:" + site + "(" + massdiff + ") doesn't exist in the library.\n");
                    } else {
                        singleLCMSID.AddModification(ptm, site);
                    }
                }
            }
        }
    }

    private void ParseSpectrumNode(Element spectrum) throws XmlPullParserException, IOException {
        PSM psm = new PSM();
        psm.SpecNumber = spectrum.getAttributes().getNamedItem("spectrum").getNodeValue();

        psm.ObserPrecursorMass = Float.parseFloat(spectrum.getAttributes().getNamedItem("precursor_neutral_mass").getNodeValue());
        psm.Charge = Integer.parseInt(spectrum.getAttributes().getNamedItem("assumed_charge").getNodeValue());
        psm.ScanNo = Integer.parseInt(spectrum.getAttributes().getNamedItem("start_scan").getNodeValue());
        if (spectrum.getAttributes().getNamedItem("retention_time_sec") != null) {
            psm.RetentionTime = Float.parseFloat(spectrum.getAttributes().getNamedItem("retention_time_sec").getNodeValue()) / 60f;
        }
        psm.NeighborMaxRetentionTime = psm.RetentionTime;

        psm.RawDataName = psm.SpecNumber.substring(0, psm.SpecNumber.indexOf("."));

        for (int k = 0; k < spectrum.getChildNodes().getLength(); k++) {
            Node resultNode = spectrum.getChildNodes().item(k);
            if ("search_result".equals(resultNode.getNodeName())) {
                for (int l = 0; l < resultNode.getChildNodes().getLength(); l++) {
                    Node hitNode = resultNode.getChildNodes().item(l);
                    if ("search_hit".equals(hitNode.getNodeName()) && "1".equals(hitNode.getAttributes().getNamedItem("hit_rank").getNodeValue())) {
                        psm.NeutralPepMass = Float.parseFloat(hitNode.getAttributes().getNamedItem("calc_neutral_pep_mass").getNodeValue());
                        float error = Float.parseFloat(hitNode.getAttributes().getNamedItem("massdiff").getNodeValue());
                        float error0 = Math.abs(error);
                        float error1 = Math.abs(error - 1f);
                        float error2 = Math.abs(error - 2f);
                        float error3 = Math.abs(error - 3f);

                        if (error0 < error1 && error0 < error2 && error0 < error3) {
                            psm.MassError = error;
                        } else if (error1 < error0 && error1 < error2 && error1 < error3) {
                            psm.MassError = error - 1f;
                            psm.ObserPrecursorMass -= 1f;
                        } else if (error2 < error0 && error2 < error1 && error2 < error3) {
                            psm.MassError = error - 2f;
                            psm.ObserPrecursorMass -= 2f;
                        } else if (error3 < error0 && error3 < error1 && error3 < error2) {
                            psm.MassError = error - 3f;
                            psm.ObserPrecursorMass -= 3f;
                        }

                        psm.PreAA = hitNode.getAttributes().getNamedItem("peptide_prev_aa").getNodeValue();
                        psm.NextAA = hitNode.getAttributes().getNamedItem("peptide_next_aa").getNodeValue();
                        psm.MissedCleavage = Integer.parseInt(hitNode.getAttributes().getNamedItem("num_missed_cleavages").getNodeValue());
                        psm.Sequence = hitNode.getAttributes().getNamedItem("peptide").getNodeValue();
                        psm.ModSeq = psm.Sequence;
                        psm.TPPModSeq = psm.Sequence;

                        String ProtACC = hitNode.getAttributes().getNamedItem("protein").getNodeValue();
                        if (!"".equals(ProtACC)) {
                            psm.AddParentProtein(ProtACC);
                        }
                        String altproACC = "";
                        boolean iprophet = false;
                        for (int m = 0; m < hitNode.getChildNodes().getLength(); m++) {
                            Node hitModNode = hitNode.getChildNodes().item(m);

                            switch (hitModNode.getNodeName()) {
                                case ("modification_info"): {
                                    GetModificationInfo(psm, hitModNode);
                                    break;
                                }
                                case ("analysis_result"): {
                                    switch (hitModNode.getAttributes().getNamedItem("analysis").getNodeValue()) {
                                        case "peptideprophet": {
                                            if (!iprophet && hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability") != null) {
                                                psm.Probability = Float.parseFloat(hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability").getNodeValue());
                                            }
                                            break;
                                        }
                                        case "interprophet": {
                                            iprophet = true;
                                            if (hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability") != null) {
                                                psm.Probability = Float.parseFloat(hitModNode.getChildNodes().item(1).getAttributes().getNamedItem("probability").getNodeValue());
                                            }
                                            break;
                                        }
                                    }
                                    break;
                                }
                                case ("search_score"): {
                                    switch (hitModNode.getAttributes().item(0).getNodeValue()) {
                                        case ("hyperscore"): {
                                            psm.hyperscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("nextscore"): {
                                            psm.nextscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("bscore"): {
                                            psm.bscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("yscore"): {
                                            psm.yscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("zscore"): {
                                            psm.zscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("ascore"): {
                                            psm.ascore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("xscore"): {
                                            psm.xscore = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                        case ("expect"): {
                                            psm.expect = Float.parseFloat(hitModNode.getAttributes().item(1).getNodeValue());
                                            break;
                                        }
                                    }
                                    break;
                                }
                                case ("alternative_protein"): {
                                    altproACC = hitModNode.getAttributes().getNamedItem("protein").getNodeValue();
                                    if (!"".equals(altproACC)) {
                                        psm.AddParentProtein(altproACC);
                                    }
                                    break;
                                }
                            }
                        }
                        if (psm.Probability > threshold) {
                            singleLCMSID.AddPSM(psm);
                        }
                    }
                }
            }
        }
    }

    private void GetModificationInfo(PSM psmid, Node node) throws XmlPullParserException, XmlPullParserException, XmlPullParserException, XmlPullParserException, XmlPullParserException, IOException {
        String PepSeq = psmid.Sequence;
        String modseq = psmid.Sequence;
        String TPPmodseq = node.getAttributes().getNamedItem("modified_peptide").getNodeValue();

        if (node.getAttributes().getNamedItem("mod_nterm_mass") != null) {
            float mass = Float.parseFloat(node.getAttributes().getNamedItem("mod_nterm_mass").getNodeValue());
            ModificationInfo matchmod = null;
            float massdiff = Float.MAX_VALUE;
            for (ModificationInfo mod : singleLCMSID.ModificationList.values()) {
                if (mod.site.equals("N-term")) {
                    float diff = Math.abs(mod.mass - mass);
                    if (diff < massdiff) {
                        massdiff = diff;
                        matchmod = mod;
                    }
                }
            }
            if (matchmod != null) {
                psmid.Modifications.add(new ModificationMatch(matchmod.modification.getName(), true, 1));
                modseq = ModStringConvert.AddModIntoSeq(modseq, matchmod.GetKey(), -1);
            } else {
                Logger.getRootLogger().warn("Modification [" + mass + " @ nterm] for spectrum: " + psmid.SpecNumber + " not found in the library:");
            }
        }
        if (node.getAttributes().getNamedItem("mod_cterm_mass") != null) {
            float mass = Float.parseFloat(node.getAttributes().getNamedItem("mod_cterm_mass").getNodeValue());

            ModificationInfo matchmod = null;
            float massdiff = Float.MAX_VALUE;
            for (ModificationInfo mod : singleLCMSID.ModificationList.values()) {
                if (mod.site.equals("C-term")) {
                    float diff = Math.abs(mod.mass - mass);
                    if (diff < massdiff) {
                        massdiff = diff;
                        matchmod = mod;
                    }
                }
            }
            if (matchmod != null) {
                psmid.Modifications.add(new ModificationMatch(matchmod.modification.getName(), true, psmid.Sequence.length()));
                modseq = ModStringConvert.AddModIntoSeq(modseq, matchmod.GetKey(), psmid.Sequence.length() - 1);
            } else {
                Logger.getRootLogger().warn("Modification [" + mass + " @ cterm] for spectrum: " + psmid.SpecNumber + " not found in the library:");
            }
        }

        for (int i = 0; i < node.getChildNodes().getLength(); i++) {
            if ("mod_aminoacid_mass".equals(node.getChildNodes().item(i).getNodeName())) {
                int idx = Integer.parseInt(node.getChildNodes().item(i).getAttributes().getNamedItem("position").getNodeValue());
                String site = String.valueOf(PepSeq.charAt(idx - 1));
                float mass = Float.parseFloat(node.getChildNodes().item(i).getAttributes().getNamedItem("mass").getNodeValue());

                ModificationInfo matchmod = null;
                float massdiff = Float.MAX_VALUE;
                for (ModificationInfo mod : singleLCMSID.ModificationList.values()) {
                    if (mod.site.equals(site)) {
                        float diff = Math.abs(mod.mass - mass);
                        if (diff < massdiff) {
                            massdiff = diff;
                            matchmod = mod;
                        }
                    }
                }
                if (matchmod != null) {
                    psmid.Modifications.add(new ModificationMatch(matchmod.modification.getName(), true, idx));
                    modseq = ModStringConvert.AddModIntoSeq(modseq, matchmod.GetKey(), idx - 1);
                } else {
                    Logger.getRootLogger().warn("Modification [" + mass + " @ " + site + "] for spectrum: " + psmid.SpecNumber + " not found in the library:");
                }
            }
        }
        psmid.ModSeq = modseq;
        psmid.TPPModSeq = TPPmodseq;
        //UpdateFromLuciphor(psmid, modseq);
    }

    private void UpdateFromLuciphor(PSM psmid, String modseq) throws NumberFormatException {
        if (singleLCMSID.LuciphorResult != null && singleLCMSID.LuciphorResult.containsKey(psmid.SpecNumber)) {
            String line = singleLCMSID.LuciphorResult.get(psmid.SpecNumber);
            if (Integer.parseInt(line.split("\t")[5]) < Integer.parseInt(line.split("\t")[6])) {
                for (int i = 0; i < 2; i++) {
                    boolean isdecoy = line.split("\t")[12 + i].equals("1");
                    if (!isdecoy) {
                        String lumodseq = line.split("\t")[2 + i];
                        String resultseq = modseq.replace("[79.96637(S)]", "").replace("[79.96637(Y)]", "").replace("[79.96633(T)]", "");

                        while (lumodseq.contains("[167]")) {
                            int aaindex = StringUtils.countMatches(lumodseq.substring(0, lumodseq.indexOf("[167]")).replaceAll("\\[(.*?)\\]", "$1"), "S");
                            int cont = 0;
                            for (int idx = 0; idx < resultseq.length(); idx++) {
                                if ("S".equals(String.valueOf(resultseq.charAt(idx)))) {
                                    cont++;
                                    if (cont == aaindex) {
                                        resultseq = resultseq.substring(0, idx) + "[79.96637(S)]" + resultseq.substring(idx);
                                    }
                                }
                            }
                            lumodseq = lumodseq.substring(0, lumodseq.indexOf("[167]")) + lumodseq.substring(lumodseq.indexOf("[167]") + 5);
                        }
                        while (lumodseq.contains("[181]")) {
                            int aaindex = StringUtils.countMatches(lumodseq.substring(0, lumodseq.indexOf("[181]")).replaceAll("\\[(.*?)\\]", "$1"), "T");
                            int cont = 0;
                            for (int idx = 0; idx < resultseq.length(); idx++) {
                                if ("T".equals(String.valueOf(resultseq.charAt(idx)))) {
                                    cont++;
                                    if (cont == aaindex) {
                                        resultseq = resultseq.substring(0, idx) + "[79.96633(T)]" + resultseq.substring(idx);
                                    }
                                }
                            }
                            lumodseq = lumodseq.substring(0, lumodseq.indexOf("[181]")) + lumodseq.substring(lumodseq.indexOf("[181]") + 5);
                        }
                        while (lumodseq.contains("[243]")) {
                            int aaindex = StringUtils.countMatches(lumodseq.substring(0, lumodseq.indexOf("[243]")).replaceAll("\\[(.*?)\\]", "$1"), "Y");
                            int cont = 0;
                            for (int idx = 0; idx < resultseq.length(); idx++) {
                                if ("Y".equals(String.valueOf(resultseq.charAt(idx)))) {
                                    cont++;
                                    if (cont == aaindex) {
                                        resultseq = resultseq.substring(0, idx) + "[79.96637(Y)]" + resultseq.substring(idx);
                                    }
                                }
                            }
                            lumodseq = lumodseq.substring(0, lumodseq.indexOf("[243]")) + lumodseq.substring(lumodseq.indexOf("[243]") + 5);
                        }
                        psmid.LuciphorLFLR = Float.parseFloat(line.split("\t")[7].replace("NA", "1"));
                        psmid.LuciphorFLR = Float.parseFloat(line.split("\t")[8].replace("NA", "1"));
                        psmid.LuciphorScore = Float.parseFloat(line.split("\t")[10 + i].replace("NA", "0"));
                        psmid.ModSeq = resultseq;
                        psmid.TPPModSeq = lumodseq;
                        return;
                    }
                }
            }
        }
    }

}
