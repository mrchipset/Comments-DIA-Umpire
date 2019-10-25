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
package MSUmpire.PostProcessing;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.FragmentSelection;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.ProtID;
import Utility.DateTimeTag;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ExportTable {

    HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments = new HashMap<>();
    HashMap<String, Float> ProteinFragMap = new HashMap<>();
    HashSet<String> IdentifiedPepMap = new HashSet<>();
    ArrayList<LCMSID> FileList;
    LCMSID CombineProtID;
    String WorkFolder;
    FragmentSelection fragselection;
       
    public ExportTable(String folder, ArrayList<LCMSID> FileList, HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments, LCMSID CombineProtID, FragmentSelection fragselection){
        this.WorkFolder=folder;
        this.FileList=FileList;
        this.CombineProtID=CombineProtID;
        this.IDSummaryFragments=IDSummaryFragments;
        this.fragselection=fragselection;
    }
    
    public void Export(int TopNPep, int TopNFrag, float Freq) throws IOException {
        Logger.getRootLogger().info("Writing reports...............");
            FileWriter proWriter = new FileWriter(WorkFolder + "ProtSummary_" + DateTimeTag.GetTag() + ".xls");
            FileWriter pepWriter = new FileWriter(WorkFolder + "PeptideSummary_" + DateTimeTag.GetTag() + ".xls");
            FileWriter NOWriter = new FileWriter(WorkFolder + "IDNoSummary_" + DateTimeTag.GetTag() + ".xls");
            FileWriter fragWriter = new FileWriter(WorkFolder + "FragSummary_" + DateTimeTag.GetTag() + ".xls");
            //Fragment Summary//////////////////////////////////////////
            for (LCMSID IDsummary : FileList) {
                HashMap<String, FragmentPeak> FragMap = IDSummaryFragments.get(IDsummary.mzXMLFileName);
                for (String key : CombineProtID.ProteinList.keySet()) {
                    if (IDsummary.ProteinList.containsKey(key)) {
                        ProtID protein = IDsummary.ProteinList.get(key);
                        for (PepIonID pep : protein.PeptideID.values()) {
                            for (FragmentPeak frag : pep.FragmentPeaks) {
                                if (!ProteinFragMap.containsKey(key + ";" + pep.GetKey() + ";" + frag.IonType)) {
                                    ProteinFragMap.put(key + ";" + pep.GetKey() + ";" + frag.IonType, frag.FragMZ);
                                }
                                frag.Prob1 = pep.MaxProbability;
                                frag.Prob2 = Math.max(pep.MS1AlignmentLocalProbability, pep.MS2AlignmentLocalProbability);
                                frag.RT = pep.PeakRT;
                                FragMap.put(key + ";" + pep.GetKey() + ";" + frag.IonType, frag);
                            }
                        }
                    }
                }
            }

            fragWriter.write("Fragment Key\tProtein\tPeptide\tFragment\tFragMz\t");
            for (LCMSID IDSummary : FileList) {
                String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
                fragWriter.write(file + "_RT\t" + file + "_Spec_Centric_Prob\t" + file + "_Pep_Centric_Prob\t" + file + "_Intensity\t" + file + "_Corr\t" + file + "_PPM\t");
            }
            fragWriter.write("\n");
            for (String key : ProteinFragMap.keySet()) {
                fragWriter.write(key + "\t" + key.split(";")[0] + "\t" + key.split(";")[1] + "\t" + key.split(";")[2] + "\t" + ProteinFragMap.get(key) + "\t");
                for (LCMSID IDSummary : FileList) {
                    if (IDSummaryFragments.get(IDSummary.mzXMLFileName).containsKey(key)) {
                        FragmentPeak fragmentPeak = IDSummaryFragments.get(IDSummary.mzXMLFileName).get(key);
                        fragWriter.write(fragmentPeak.RT + "\t" + fragmentPeak.Prob1 + "\t" + fragmentPeak.Prob2 + "\t" + fragmentPeak.intensity + "\t" + fragmentPeak.corr + "\t" + fragmentPeak.ppm + "\t");
                    } else {
                        fragWriter.write("\t\t\t\t\t\t");
                    }
                }
                fragWriter.write("\n");
            }
            fragWriter.close();

            //ProteinSummary/////////////
            NOWriter.write("File\tNo. Proteins\tNo. peptide ions (Spec-centric)\tNo. peptide ions (Pep-centric)\tNo. proein assoc. ions\n");
            proWriter.write("Protein Key\t");
            for (LCMSID IDSummary : FileList) {
                NOWriter.write(FilenameUtils.getBaseName(IDSummary.mzXMLFileName) + "\t" + IDSummary.ProteinList.size() + "\t" + IDSummary.GetPepIonList().size() + "\t" + IDSummary.GetMappedPepIonList().size() + "\t" + IDSummary.AssignedPepIonList.size() + "\n");
                String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
                proWriter.write(file + "_Prob\t" + file + "_Peptides\t" + file + "_PSMs\t" + file + "_MS1_iBAQ\t" + file + "_Top" + TopNPep + "pep/Top" + TopNFrag + "fra ,Freq>" + Freq + "\t");
            }
            proWriter.write("\n");

            NOWriter.close();
            for (String key : CombineProtID.ProteinList.keySet()) {
                proWriter.write(key + "\t");
                for (LCMSID IDsummary : FileList) {
                    if (IDsummary.ProteinList.containsKey(key)) {
                        ProtID protein = IDsummary.ProteinList.get(key);
                        proWriter.write(protein.Probability + "\t" + protein.PeptideID.size() + "\t" + protein.GetSpectralCount() + "\t" + protein.GetAbundanceByMS1_IBAQ() + "\t" + protein.GetPepAbundanceByTopCorrFragAcrossSample(fragselection.TopPeps.get(protein.getAccNo()), fragselection.TopFrags) + "\t");
                    } else {
                        proWriter.write("\t\t\t\t\t");
                    }
                }
                proWriter.write("\n");
            }
            proWriter.close();

            ////PepSummary///////////////////////////////////
            for (LCMSID IDsummary : FileList) {
                HashMap<String, FragmentPeak> FragMap = new HashMap<>();
                IDSummaryFragments.put(IDsummary.mzXMLFileName, FragMap);

                for (String key : IDsummary.GetPepIonList().keySet()) {
                    if (!IdentifiedPepMap.contains(key)) {
                        IdentifiedPepMap.add(key);
                    }
                }
                for (String key : IDsummary.GetMappedPepIonList().keySet()) {
                    if (!IdentifiedPepMap.contains(key)) {
                        IdentifiedPepMap.add(key);
                    }
                }
            }

            pepWriter.write("Peptide Key\tSequence\tModSeq\tProteins\tmz\tCharge\tMaxProb\t");

            for (LCMSID IDSummary : FileList) {
                String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
                pepWriter.write(file + "_Spec_Centric_Prob\t" + file + "_Pep_Centric_Prob\t" + file + "_PSMs\t" + file + "_RT\t" + file + "_MS1\t" + file + "_Top" + TopNFrag + "fra\t");
            }
            pepWriter.write("\n");

            for (String key : IdentifiedPepMap) {
                pepWriter.write(key + "\t");
                float maxprob = 0f;
                boolean output = false;
                for (LCMSID IDSummary : FileList) {
                    if (IDSummary.GetPepIonList().containsKey(key)) {
                        PepIonID peptide = IDSummary.GetPepIonList().get(key);
                        if (!output) {
                            pepWriter.write(peptide.Sequence + "\t" + peptide.ModSequence + "\t" + peptide.ParentProteins() + "\t" + peptide.NeutralPrecursorMz() + "\t" + peptide.Charge + "\t");
                            output = true;
                        }
                        if (peptide.MaxProbability > maxprob) {
                            maxprob = peptide.MaxProbability;
                        }
                    }
                }
                pepWriter.write(maxprob + "\t");
                for (LCMSID IDSummary : FileList) {
                    if (IDSummary.GetPepIonList().containsKey(key)) {
                        PepIonID peptide = IDSummary.GetPepIonList().get(key);
                        pepWriter.write(peptide.MaxProbability + "\t-1\t" + peptide.GetSpectralCount() + "\t" + peptide.PeakRT + "\t" + peptide.PeakHeight[0] + "\t" + peptide.GetPepAbundanceByTopCorrFragAcrossSample(fragselection.TopFrags.get(peptide.GetKey())) + "\t");
                    } else if (IDSummary.GetMappedPepIonList().containsKey(key)) {
                        PepIonID peptide = IDSummary.GetMappedPepIonList().get(key);
                        pepWriter.write("-1\t" + Math.max(peptide.MS1AlignmentLocalProbability, peptide.MS2AlignmentLocalProbability) + "\t" + peptide.GetSpectralCount() + "\t" + peptide.PeakRT + "\t" + peptide.PeakHeight[0] + "\t" + peptide.GetPepAbundanceByTopFragments(6) + "\t");
                    } else {
                        pepWriter.write("\t\t\t\t\t\t");
                    }
                }
                pepWriter.write("\n");
            }
            pepWriter.close();

        Logger.getRootLogger().info("Job done");
        Logger.getRootLogger().info("=================================================================================================");
    }

}
