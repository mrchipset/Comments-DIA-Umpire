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

import MSUmpire.PeakDataStructure.PeakCluster;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.IonFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.general.IsotopicDistribution;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.MolecularFormula;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepIonID implements Serializable {
    private static final long serialVersionUID = 33203837556783L;

    public int Index;
    public String ModSequence;
    public float Weight;
    public boolean Is_NonDegenerate;
    public int Charge;
    public float GroupWeight;
    public transient float FilteringWeight;
    public String Sequence;
    public float MaxProbability = -1f;
    public boolean GlycoMS1Valid=false;
    public boolean GlycoMS2Valid=false;
    public String MS1ClusIndex = "";
    public String MS2ClusIndex = "";
    public int NoNsite;
    public int NoNsiteFragObs=-1;
    public int IsDecoy=-1;//added 0828, needs to be deleted for older serialization
    
    private ArrayList<PSM> PSMList;
    public ArrayList<String> ParentProtID_PepXML = new ArrayList<>();
    public ArrayList<ProtID> ParentProt_ProtXML = new ArrayList<>();
    public ArrayList<PeakCluster> MS1PeakClusters = new ArrayList<>();
    public ArrayList<PeakCluster> MS2UnfragPeakClusters = new ArrayList<>();
    public float PeakArea[];
    public float PeakHeight[];
    public float PeakClusterScore;
    private float RT = -1f;
    public float PeakRT = -1f;
    public ArrayList<Float> PredictRT = new ArrayList<>();
    //public int PeakClusterIndex = -1;
    private float RTSD = -1f;
    private ArrayList<Ion> Fragments;
    private SortedFragmentSet FragmentMZSet;
    public ArrayList<ModificationMatch> Modifications = new ArrayList<>();
    public ArrayList<FragmentPeak> FragmentPeaks = new ArrayList<>();    
    private Peptide peptide = null;
    private float mz = -1f;
    public float ObservedMz;
    public String TPPModSeq;
    public float MS1AlignmentProbability = -1f;
    public float MS2AlignmentProbability = -1f;
    public float MS1AlignmentLocalProbability= -1f;
    public float MS2AlignmentLocalProbability= -1f;

    public Peptide GetPepFactory() {
        if (peptide == null) {
            ReFreshPepFactory();
        }
        return peptide;
    }
    
    public void ClearPepFragFactory(){
        peptide=null;
        FragmentMZSet=null;
        Fragments=null;
    }

    public float GetMaxLuciphorScore() {
        float score = 0f;
        for (PSM psm : PSMList) {
            if (psm.LuciphorScore > score) {
                score = psm.LuciphorScore;
            }
        }
        return score;
    }

    public float GetMinLuciphorLFLR() {
        float flr = 1f;
        for (PSM psm : PSMList) {
            if (psm.LuciphorLFLR < flr) {
                flr = psm.LuciphorLFLR;
            }
        }
        return flr;
    }

    public float GetMinLuciphorFLR() {
        float flr = 1f;
        for (PSM psm : PSMList) {
            if (psm.LuciphorFLR < flr) {
                flr = psm.LuciphorFLR;
            }
        }
        return flr;
    }

    public void ReFreshPepFactory() {
        peptide = new Peptide(Sequence, Modifications);
        mz = -1f;
    }

    public String GetMS1ClusIndex() {
        if ("".equals(MS1ClusIndex)) {
            for (PeakCluster cluster : MS1PeakClusters) {
                MS1ClusIndex += cluster.Index + ";";
            }
        }
        return MS1ClusIndex;
    }

    public String GetMS2ClusIndex() {
        if ("".equals(MS2ClusIndex)) {
            for (PeakCluster cluster : MS2UnfragPeakClusters) {
                MS2ClusIndex += cluster.Index + ";";
            }
        }
        return MS2ClusIndex;
    }

    public PepIonID ClonePepIonID() {
        PepIonID newpeptide = new PepIonID();
        newpeptide.Sequence = Sequence;
        newpeptide.Charge = Charge;
        newpeptide.Is_NonDegenerate = Is_NonDegenerate;
        newpeptide.ModSequence = ModSequence;
        newpeptide.Modifications = (ArrayList<ModificationMatch>)Modifications.clone();
        if (PeakArea != null) {
            newpeptide.PeakArea = new float[PeakArea.length];
            newpeptide.PeakHeight = new float[PeakHeight.length];
        }
//        if (PSMList != null) {
//            newpeptide.PSMList = (ArrayList<PSM>) PSMList.clone();
//        }
        return newpeptide;
    }

    
    public ArrayList<Ion> GetFragments() {
        if (Fragments == null) {
            SetFragments();
        }
        return Fragments;
    }
    
    public float GetPepAbundanceByTopCorrFragAcrossSample(ArrayList<String> PepFrag) {
        float totalabundance = 0f;
        int count=0;
        if (PepFrag!=null && FragmentPeaks != null) {
            for (FragmentPeak frag : FragmentPeaks) {
                if (PepFrag.contains(frag.GetFragKey())) {
                    totalabundance += frag.intensity;
                    count++;
                }
            }
        }
        if(count<2){
            totalabundance=0;
        }
        return totalabundance;
    }
    
    public float GetPepAbundanceByTopCorrFragments(int topN) {        
        ArrayList<FragmentPeak> includelist = new ArrayList<>();
        for (int i = 0; i < topN; i++) {
            FragmentPeak bestfrag=null;
            for (FragmentPeak fragment : FragmentPeaks) {
                if (!includelist.contains(fragment) && (bestfrag == null || fragment.intensity > bestfrag.intensity)) {
                    bestfrag = fragment;
                }
            }
            if (bestfrag != null) {
                includelist.add(bestfrag);
            }
        }
        float totalabundance = 0f;
        for (FragmentPeak frag : includelist ) {
            totalabundance+=frag.intensity;
        }
        return totalabundance;
    }
     
    public float GetPepAbundanceByTopFragments(int topN) {        
        ArrayList<FragmentPeak> includelist = new ArrayList<>();
        for (int i = 0; i < topN; i++) {
            FragmentPeak bestfrag=null;
            for (FragmentPeak fragment : FragmentPeaks) {
                if (!includelist.contains(fragment) && (bestfrag == null || fragment.intensity > bestfrag.intensity)) {
                    bestfrag = fragment;
                }
            }
            if (bestfrag != null) {
                includelist.add(bestfrag);
            }
        }
        float totalabundance = 0f;
        for (FragmentPeak frag : includelist ) {
            totalabundance+=frag.intensity;
        }
        return totalabundance;
    }

    public float GetPepAbundanceByAllFragments() {
        float totalabundance = 0f;
        for (FragmentPeak fragment : FragmentPeaks) {
            totalabundance += fragment.intensity;
        }
        return totalabundance;
    }

    public PSM GetBestExpectValuePSM() {
        PSM top = GetPSMList().get(0);
        for (int psmidx = 1; psmidx < GetPSMList().size(); psmidx++) {
            if (GetPSMList().get(psmidx).expect < top.expect) {
                top = GetPSMList().get(psmidx);
            }
        }
        return top;
    }

    public PSM GetBestPSM() {
        PSM bestpsm = null;
        for (PSM psm : GetPSMList()) {
            if (bestpsm == null || psm.Probability > bestpsm.Probability) {
                bestpsm = psm;
            } else if (psm.Probability == bestpsm.Probability && psm.expect < bestpsm.expect) {
                bestpsm = psm;
            }
        }
        return bestpsm;
    }

    public String GetKey() {
        return ModSequence + "_" + Charge;
    }
    
    public boolean IsDecoy(String decoytag) {
        if (IsDecoy == -1) {
            IsDecoy = 1;
            for (String pro : ParentProtID_PepXML) {
                if (!pro.startsWith(decoytag)) {
                    IsDecoy = 0;
                    break;
                }
            }
        }
        return IsDecoy==1;
    }

    public float LookUpFragmentMZ(PeptideFragmentIon.IonType fragmentIonType, int num) {        
        for (Ion frag : GetFragments()) {
            if (frag.getType() == fragmentIonType && "".equals(frag.getNeutralLossesAsString()) && ((PeptideFragmentIon) frag).getNumber() == num) {
                return (float) (frag.getTheoreticMass() + ElementaryIon.proton.getTheoreticMass());
            }
        }
        return 0f;
    }

    public SortedFragmentSet GetFragmentMZ() {

        if (FragmentMZSet == null) {
            FragmentMZSet = new SortedFragmentSet();
            double protonMass = ElementaryIon.proton.getTheoreticMass();
            for (Ion frag : GetFragments()) {
                float targetmz = (float) frag.getTheoreticMz(1);
                PeptideFragment fragment = new PeptideFragment();
                fragment.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                fragment.Charge = 1;
                fragment.FragMZ = targetmz;
                //fragment.Number = ((PeptideFragmentIon) frag).getNumber();
                FragmentMZSet.add(fragment);
                float targetmz2 = (float) frag.getTheoreticMz(2);
                PeptideFragment fragment2 = new PeptideFragment();
                fragment2.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                fragment2.Charge = 2;
                fragment2.FragMZ = targetmz2;
                //fragment2.Number = ((PeptideFragmentIon) frag).getNumber();
                FragmentMZSet.add(fragment2);
            }
        }
        return FragmentMZSet;
    }
    
    private void SetFragments() {
        peptide = GetPepFactory();
        HashMap<Integer, HashMap<Integer, ArrayList<Ion>>> allfragment = IonFactory.getInstance().getFragmentIons(peptide);
        Fragments = new ArrayList<>();
        for (Ion frag : allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.B_ION)) {
            Fragments.add(frag);
        }
        for (Ion frag : allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.Y_ION)) {
            Fragments.add(frag);
        }
    }

    public void CalcRT() {
        if(PSMList.isEmpty()){
            return ;
        }
        RT = 0f;
        for (PSM psm : PSMList) {
            RT += psm.NeighborMaxRetentionTime;
        }
        RT /= PSMList.size();

        RTSD = 0f;
        for (PSM psm : PSMList) {
            RTSD += (psm.NeighborMaxRetentionTime - RT) * (psm.NeighborMaxRetentionTime - RT);
        }

        RTSD /= PSMList.size();
        RTSD = (float) Math.sqrt(RTSD);
    }

    public void SetRT(float value) {
        RT = value;
    }

    public float GetRTSD() {
        if (RTSD == -1) {
            CalcRT();
        }
        return RTSD;
    }

    public float GetRT() {
        if (PeakRT != -1) {
            return PeakRT;
        }
        return GetIDRT();
    }
    
    
    public float GetAvgPredictRT(){
        float avg=0f;
        for(float rt : PredictRT){
            avg+=rt;
        }
        return avg/=PredictRT.size();
    }

    
    public String PredictRTString(){
        String rtout="";
        for(float rt : PredictRT){
            rtout+=rt+";";
        }
        return rtout;
    }

    public float GetIDRT() {
        if (RT == -1) {
            CalcRT();
        }
        return RT;
    }

    public void CreateQuantInstance(int NoPeakClusters) {
        PeakArea = new float[NoPeakClusters];
        PeakHeight = new float[NoPeakClusters];
    }

    public PepIonID() {
        PSMList = new ArrayList<>();
    }

    public void AddPSM(PSM psm) {
        if (psm.Probability > MaxProbability) {
            MaxProbability = psm.Probability;
            ObservedMz = psm.ObserPrecursorMz();
        }
        PSMList.add(psm);
        psm.pepIonID = this;
        for (String protein : psm.ParentProtIDs) {
            if (!ParentProtID_PepXML.contains(protein)) {
                ParentProtID_PepXML.add(protein);
            }
        }

        if (Modifications.isEmpty()) {
            AddModifications(psm);
        }
        if (Modifications.size() != psm.Modifications.size()) {
            Logger.getRootLogger().error("Modification size doesn't match.........:"+psm.SpecNumber);
        }
    }

    public void AddModifications(PSM psm) {
        for (ModificationMatch mod : psm.Modifications) {
            Modifications.add(mod);
        }
    }

    
    //private transient int spc=-1;
    
    public int GetSpectralCount(){
//        if(spc==-1){
//            spc=PSMList.size();
//        }
        return PSMList.size();
    }
    
    public ArrayList<PSM> GetPSMList() {
        return PSMList;
    }

    public void SetInfobyPSM(PSM psm) {
        this.Charge = psm.Charge;
        this.ModSequence = psm.ModSeq;
        this.TPPModSeq = psm.TPPModSeq;
        this.Sequence = psm.Sequence;
        this.MaxProbability = psm.Probability;
        this.ObservedMz = psm.ObserPrecursorMz();
    }

    public void SetMz(float mz) {
        this.mz = mz;
    }

    public float ObservedMass() {
        return Charge * (ObservedMz - 1.00727f);
    }

    public float NeutralPrecursorMz() {
        if (mz == -1f) {
            mz = (CalcNeutralPepMass() + Charge * (float) ElementaryIon.proton.getTheoreticMass()) / Charge;
        }
        return mz;
    }

    public float CalcNeutralPepMass() {
        return GetPepFactory().getMass().floatValue();
    }

    public String GetModificationString() {
        String ModificationString = "";
        for (ModificationMatch mod : Modifications) {
            ModificationString += mod.getTheoreticPtm() + "(" + mod.getModificationSite() + ");";
        }
        return ModificationString;
    }

    public MolecularFormula GetMolecularFormula() {
        MolecularFormula formula = new MolecularFormula(GetAASequenceImpl());
        //formula.addMolecularFormula(GetModMolecularFormula());
        return formula;
    }
    AASequenceImpl AAimple = null;

    public AASequenceImpl GetAASequenceImpl() {
        if (AAimple == null) {
            AAimple = new AASequenceImpl(Sequence);
        }
        return AAimple;
    }
    IsotopicDistribution calc = null;

    public float[] IsotopicDistrubtion(int NoOfIsoPeaks) {
        IsotopicDistribution calc = GetIsotopicDistribution();
        Double[] isopeak = calc.getPercMax();
        float firstPattern = (float) (double) isopeak[0];

        float[] TheoIso = new float[NoOfIsoPeaks];
        for (int i = 0; i < NoOfIsoPeaks; i++) {
            TheoIso[i] = (float) (double) (isopeak[i] / firstPattern);
        }
        return TheoIso;
    }

    private IsotopicDistribution GetIsotopicDistribution() {
        if (calc == null) {
            calc = new IsotopicDistribution(GetMolecularFormula());
            //calc.calculate();
        }
        //AAimple.addModification(new ModificationImplementation(Modifications , Modifications, new HashMap(), 0));
        return calc;
    }

    public String ParentProteins() {
        String proteins="";
        for(ProtID prot : ParentProt_ProtXML){
            proteins+=prot.getAccNo()+";";
        }
        return proteins;
    }

    public float GetMS1() {
        if (PeakHeight != null) {
            return PeakHeight[0];
        }
        return 0f;
    }

    public void RemoveRedundantFrag() {
        HashMap<String, FragmentPeak> UniqueFrag = new HashMap<>();
        for (FragmentPeak frag : FragmentPeaks) {
            if (!UniqueFrag.containsKey(frag.GetFragKey())) {
                UniqueFrag.put(frag.GetFragKey(), frag);
            } else if (UniqueFrag.get(frag.GetFragKey()).intensity * UniqueFrag.get(frag.GetFragKey()).corr < frag.intensity * frag.corr) {
                UniqueFrag.put(frag.GetFragKey(), frag);
            }
        }
        FragmentPeaks.clear();
        for (FragmentPeak frag : UniqueFrag.values()) {
            FragmentPeaks.add(frag);
        }
    }

}
