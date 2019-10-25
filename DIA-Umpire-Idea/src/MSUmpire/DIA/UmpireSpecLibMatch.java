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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.LCMSBaseStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.FragmentPeakGroup;
import MSUmpire.PSMDataStructure.PepFragmentLib;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.SpectralProcessingModule.ScoreFunction;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class UmpireSpecLibMatch implements Runnable, Serializable{
    private static final long serialVersionUID = 3164978164L;

    transient PepFragmentLib fragmentLib;
    transient PepFragmentLib decoyfragmentLib;
    transient InstrumentParameter parameter;
    transient LCMSPeakMS1 ms1lcms;
    transient LCMSPeakDIAMS2 DIAWindow;
    PepIonID pepIonID;
    public PeakGroupScore BestMS1Hit;
    public PeakGroupScore BestMS2Hit;
    public PeakGroupScore BestHit;
    public PeakGroupScore BestMS1DecoyHit;
    public PeakGroupScore BestMS2DecoyHit;
    public ArrayList<PeakGroupScore> TargetHits;
    public ArrayList<PeakGroupScore> DecoyHits;
    public boolean IdentifiedPeptideIon = false;

    public UmpireSpecLibMatch(LCMSPeakMS1 ms1lcms, LCMSPeakDIAMS2 DIAWindow, PepIonID pepIonID, PepFragmentLib fragmentLib, PepFragmentLib decoyfragmentLib, InstrumentParameter parameter) {
        this.ms1lcms = ms1lcms;
        this.DIAWindow = DIAWindow;
        this.fragmentLib = fragmentLib;
        this.decoyfragmentLib = decoyfragmentLib;
        this.parameter = parameter;
        this.pepIonID = pepIonID;
    }

    public class MatchFragment {

        public FragmentPeakGroup libfrag;
        public PrecursorFragmentPairEdge peakfrag;

        public MatchFragment(FragmentPeakGroup libfrag, PrecursorFragmentPairEdge peakfrag) {
            this.libfrag = libfrag;
            this.peakfrag = peakfrag;
        }
    }

    private void CalMatchScore(PeakCluster cluster, PepFragmentLib fragmentLib, PeakGroupScore peakscore) {                
        XYPointCollection pointset = new XYPointCollection();
        ArrayList<MatchFragment> matchFragments = new ArrayList<>();        
        for (FragmentPeakGroup frag : fragmentLib.FragmentGroups.values()) {
            PrecursorFragmentPairEdge bestfragment = null;
            for (PrecursorFragmentPairEdge fragmentClusterUnit : cluster.Fragments) {                
                if (InstrumentParameter.CalcPPM(frag.FragMZ, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }
            matchFragments.add(new MatchFragment(frag, bestfragment));
        }

        //float maxint = 0f;
//        for (MatchFragment matchFragment : matchFragments) {
//            if (matchFragment.peakfrag != null && matchFragment.peakfrag.Intensity > maxint) {
//                maxint = matchFragment.peakfrag.Intensity;
//            }
//        }
        for (MatchFragment matchFragment : matchFragments) {
            if (matchFragment.peakfrag != null) {
                pointset.AddPoint(matchFragment.libfrag.GetAvgInt(), matchFragment.peakfrag.Intensity);
                peakscore.IntScore += matchFragment.peakfrag.Intensity;         
                peakscore.CorrScore += matchFragment.peakfrag.Correlation;
                peakscore.ApexDeltaScore += 1 - matchFragment.peakfrag.ApexDelta;
                peakscore.PPMScore += 1 - InstrumentParameter.CalcPPM(matchFragment.peakfrag.FragmentMz, matchFragment.libfrag.FragMZ) / parameter.MS2PPM;
                peakscore.RTOverlapScore += matchFragment.peakfrag.RTOverlapP;
                if (matchFragment.libfrag.IonType.startsWith("b")) {
                    peakscore.NoMatchB++;
                }
                if (matchFragment.libfrag.IonType.startsWith("y")) {
                    peakscore.NoMatchY++;
                }
            } else {
                pointset.AddPoint(matchFragment.libfrag.GetAvgInt(), 0.000000000001f);
            }
        }
        int totmamatch=peakscore.NoMatchB+peakscore.NoMatchY;
        peakscore.NoFragmentLib = matchFragments.size();
        peakscore.IntScore = (float) Math.log(peakscore.IntScore/totmamatch);
        peakscore.CorrScore /= totmamatch;
        peakscore.ApexDeltaScore /= totmamatch;
        peakscore.PPMScore /= totmamatch;
        peakscore.RTOverlapScore /= totmamatch;
                        
        XYPointCollection norP=ScoreFunction.SpectralNormalizationForPairCollection(pointset);
        peakscore.SpecCorrelation=ScoreFunction.CalcSpecCorrForPairPointCollection(norP);
        peakscore.SpecDotProduct=ScoreFunction.CalcDotProductForPairPointCollection(norP);
        peakscore.ContrastAngle=ScoreFunction.CalcSpecContrastAngleForPairPointCollection(norP);        
    }

    @Override
    public void run() {
        if (fragmentLib != null) {
            if (IdentifiedPeptideIon) {
                IDPeak();
                if (decoyfragmentLib != null) {
                    Decoy();
                }
                else{
                    Logger.getRootLogger().error("decoy spectrum is null : "+pepIonID.GetKey());
                }
            } else {
                Target();
                if (decoyfragmentLib != null) {
                    Decoy();
                }
                else{
                    Logger.getRootLogger().error("decoy spectrum is null : "+pepIonID.GetKey());
                }
            }
        }
        else{
            Logger.getRootLogger().warn("lib spectrum is null : "+pepIonID.GetKey());
        }
    }

    public void IDPeak() {
        ArrayList<PeakCluster> clusterList = pepIonID.MS1PeakClusters;
        TargetHits = new ArrayList<>();
        for (PeakCluster Assigncluster : clusterList) {
            PeakCluster cluster=ms1lcms.PeakClusters.get(Assigncluster.Index-1);
            DIAWindow.ExtractFragmentForPeakCluser(cluster);            
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 1;
            peakGroupScore.MS1Score = cluster.MS1Score;
            peakGroupScore.MS1Corr=cluster.Corrs[0];            
             if(!pepIonID.PredictRT.isEmpty()){
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());                        
            }
            else {
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
            }
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());            
            CalMatchScore(cluster, fragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                TargetHits.add(peakGroupScore);
            }
        }
        clusterList = pepIonID.MS2UnfragPeakClusters;
        for (PeakCluster Assigncluster : clusterList) {
            if (DIAWindow.PeakClusters.size() >= Assigncluster.Index) {
                PeakCluster cluster = DIAWindow.PeakClusters.get(Assigncluster.Index - 1);
                if (cluster.TargetMz() == Assigncluster.TargetMz() || cluster.Charge == Assigncluster.Charge) {
                    DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
                    PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
                    peakGroupScore.MSlevel = 2;
                    peakGroupScore.MS1Score = cluster.MS1Score;
                    peakGroupScore.MS1Corr = cluster.Corrs[0];
                    if (!pepIonID.PredictRT.isEmpty()) {
                        peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());
                    } else {
                        peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
                    }
                    peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());
                    CalMatchScore(cluster, fragmentLib, peakGroupScore);
                    if (peakGroupScore.NoMatchB + peakGroupScore.NoMatchY > 0) {
                        TargetHits.add(peakGroupScore);
                    }
                }
            }
        }
    }

    public void Target() {        
        ArrayList<PeakCluster> clusterList = ms1lcms.FindAllPeakClustersForMappedPep(pepIonID);
        TargetHits = new ArrayList<>();
        for (PeakCluster cluster : clusterList) {
            DIAWindow.ExtractFragmentForPeakCluser(cluster);
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 1;
            peakGroupScore.MS1Score = cluster.MS1Score;
            peakGroupScore.MS1Corr=cluster.Corrs[0];
            peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());            
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());            
            CalMatchScore(cluster, fragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                TargetHits.add(peakGroupScore);
            }
        }
        clusterList = DIAWindow.FindAllPeakClustersForMappedPep(pepIonID);
        for (PeakCluster cluster : clusterList) {
            DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 2;
            peakGroupScore.MS1Score = cluster.MS1Score;
            peakGroupScore.MS1Corr=cluster.Corrs[0];
            peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());            
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());
            CalMatchScore(cluster, fragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                TargetHits.add(peakGroupScore);
            }
        }
    }

    public void Decoy() {
        ArrayList<PeakCluster> clusterList = ms1lcms.FindAllPeakClustersForMappedPep(pepIonID);
        DecoyHits = new ArrayList<>();
        for (PeakCluster cluster : clusterList) {
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 1;
            peakGroupScore.MS1Score = cluster.MS1Score;
            peakGroupScore.MS1Corr=cluster.Corrs[0];            
            if(!pepIonID.PredictRT.isEmpty()){
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());                        
            }
            else {
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
            }
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());
            CalMatchScore(cluster, decoyfragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY > 0) {
                DecoyHits.add(peakGroupScore);
            }
        }
        clusterList = DIAWindow.FindAllPeakClustersForMappedPep(pepIonID);
        for (PeakCluster cluster : clusterList) {
            PeakGroupScore peakGroupScore = new PeakGroupScore(cluster);
            peakGroupScore.MSlevel = 2;
            peakGroupScore.MS1Score = cluster.MS1Score;
            peakGroupScore.MS1Corr=cluster.Corrs[0];
             if(!pepIonID.PredictRT.isEmpty()){
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetAvgPredictRT());                        
            }
            else {
                peakGroupScore.RTDiff = Math.abs(cluster.PeakHeightRT[0] - pepIonID.GetIDRT());
            }           
            peakGroupScore.PrecursorPPM = InstrumentParameter.CalcPPM(cluster.TargetMz(), pepIonID.NeutralPrecursorMz());            
            CalMatchScore(cluster, decoyfragmentLib, peakGroupScore);
            if (peakGroupScore.NoMatchB+peakGroupScore.NoMatchY> 0) {
                DecoyHits.add(peakGroupScore);
            }
        }
    }

    public void AssignProbToPepIon() {
        pepIonID.MS1AlignmentLocalProbability = 0f;
        pepIonID.MS1AlignmentProbability = 0f;
        pepIonID.MS2AlignmentLocalProbability = 0f;
        pepIonID.MS2AlignmentProbability = 0f;
        if (BestMS1Hit != null) {
            pepIonID.MS1PeakClusters.add(BestMS1Hit.cluster);
            pepIonID.MS1AlignmentProbability = BestMS1Hit.MixtureModelProb;
            pepIonID.MS1AlignmentLocalProbability = BestMS1Hit.MixtureModelLocalProb;
        }
        if (BestMS2Hit != null) {
            pepIonID.MS2UnfragPeakClusters.add(BestMS2Hit.cluster);
            pepIonID.MS2AlignmentProbability = BestMS2Hit.MixtureModelProb;
            pepIonID.MS2AlignmentLocalProbability = BestMS2Hit.MixtureModelLocalProb;
        }
    }

    public void Ranking() {
        Collections.sort(TargetHits, new Comparator<PeakGroupScore>() {
            @Override
            public int compare(PeakGroupScore o1, PeakGroupScore o2) {
                return -Float.compare(o1.UmpireScore, o2.UmpireScore);
            }
        });

        Collections.sort(DecoyHits, new Comparator<PeakGroupScore>() {
            @Override
            public int compare(PeakGroupScore o1, PeakGroupScore o2) {
                return -Float.compare(o1.UmpireScore, o2.UmpireScore);
            }
        });

        if (!TargetHits.isEmpty()) {
            BestHit = TargetHits.get(0);
        }
        for (int i = 0; i < TargetHits.size(); i++) {
            if (BestMS1Hit != null && BestMS2Hit != null) {
                break;
            }
            if (TargetHits.get(i).MSlevel == 1 && BestMS1Hit == null) {
                BestMS1Hit = TargetHits.get(i);
            }
            if (TargetHits.get(i).MSlevel == 2 && BestMS2Hit == null) {
                BestMS2Hit = TargetHits.get(i);
            }
        }

        for (int i = 0; i < DecoyHits.size(); i++) {
            if (BestMS1DecoyHit != null && BestMS2DecoyHit != null) {
                break;
            }
            if (DecoyHits.get(i).MSlevel == 1 && BestMS1DecoyHit == null) {
                BestMS1DecoyHit = DecoyHits.get(i);
            }
            if (DecoyHits.get(i).MSlevel == 2 && BestMS2DecoyHit == null) {
                BestMS2DecoyHit = DecoyHits.get(i);
            }
        }
    }
}
