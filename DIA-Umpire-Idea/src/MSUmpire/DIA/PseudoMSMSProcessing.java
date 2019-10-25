/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational medicine and bioinformatics, 
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
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PseudoMSMSProcessing implements Runnable{

    InstrumentParameter parameter;
    ArrayList<PrecursorFragmentPairEdge> fragments;
    public PeakCluster ms1cluster;
    float growth=1f;
    
    public PseudoMSMSProcessing(PeakCluster ms1cluster, ArrayList<PrecursorFragmentPairEdge> fragments, InstrumentParameter parameter) {        
        this.parameter = parameter;
        this.fragments=fragments;
        this.ms1cluster=ms1cluster;
    }
        
     public void TopNPeakPerMZRangeSelection(int topN, float MassWindow) {
        
         XYData massrange=new XYData(ms1cluster.NeutralMass()-MassWindow, ms1cluster.NeutralMass());
         PrecursorFragmentPairEdge[] TopList = new PrecursorFragmentPairEdge[topN];
        ArrayList<PrecursorFragmentPairEdge> includelist = new ArrayList<>();
//        for (int i = 0; i < topN; i++) {
//            for (FragmentPeak fragment : FragmentPeaks) {
//                if (!includelist.contains(fragment) && (TopList[i] == null || fragment.intensity > TopList[i].intensity)) {
//                    TopList[i] = fragment;
//                }
//            }
//            if (TopList[i] != null) {
//                includelist.add(TopList[i]);
//            }
//        }
//        for (int i = fragments.size()-1; i >0; i--) {
//            if (InstrumentParameter.CalcPPM(fragments.get(i).FragmentMz, currentmaxfragment.FragmentMz) > parameter.MS2PPM) {
//                fragmentmarked[currentmaxindex] = false;
//                currentmaxindex = i;
//                currentmaxfragment = fragments.get(i);
//            } else if (fragments.get(i).Intensity > currentmaxfragment.Intensity) {
//                currentmaxindex = i;
//                currentmaxfragment = fragments.get(i);
//            }
//        }
//        fragmentmarked[currentmaxindex] = false;
//        for (int i = 0; i < fragments.size(); i++) {
//            if (!fragmentmarked[i]) {
//                fragmentmarked[i] = true;
//                PrecursorFragmentPairEdge startfrag = fragments.get(i);
//
//                boolean groupped = false;
//                for (int charge = 2; charge >= 1; charge--) {
//                    float lastint = startfrag.Intensity;
//                    boolean found = false;
//                    for (int pkidx = 1; pkidx < 5; pkidx++) {
//                        float targetmz = startfrag.FragmentMz + (float) pkidx / charge;
//                        for (int j = i + 1; j < fragments.size(); j++) {
//                            if (!fragmentmarked[j]) {
//                                PrecursorFragmentPairEdge targetfrag = fragments.get(j);
//                                if (InstrumentParameter.CalcPPM(targetfrag.FragmentMz, targetmz) < parameter.MS2PPM * (pkidx * 0.5 + 1)) {
//                                    if (targetfrag.Intensity < lastint) {
//                                        fragmentmarked[j] = true;
//                                        lastint = targetfrag.Intensity;
//                                        found = true;
//                                        break;
//                                    }
//                                } else if (targetfrag.FragmentMz > targetmz) {
//                                    break;
//                                }
//                            }
//                        }
//                        if (!found) {
//                            break;
//                        }
////                    }
////                    if (found) {
////                        groupped = true;
////                        //convert to charge 1 m/z
////                        startfrag.FragmentMz = startfrag.FragmentMz * charge - (charge - 1) * (float) ElementaryIon.proton.getTheoreticMass();
////                        newfragments.add(startfrag);
////                    }
////                }
////                if (!groupped) {
////                    newfragments.add(startfrag);
////                }
////            }
//        }
//        fragments.clear();
//        for (PrecursorFragmentPairEdge fragment : newfragments) {
//            fragments.add(fragment);
//        }
    }
    
    public void DeisotopingForPeakClusterFragment() {
        ArrayList<PrecursorFragmentPairEdge> newfragments = new ArrayList<>();
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.TRUE);
        PrecursorFragmentPairEdge currentmaxfragment = fragments.get(0);
        int currentmaxindex = 0;
        for (int i = 1; i < fragments.size(); i++) {
            if (InstrumentParameter.CalcPPM(fragments.get(i).FragmentMz, currentmaxfragment.FragmentMz) > parameter.MS2PPM) {
                fragmentmarked[currentmaxindex] = false;
                currentmaxindex = i;
                currentmaxfragment = fragments.get(i);
            } else if (fragments.get(i).Intensity > currentmaxfragment.Intensity) {
                currentmaxindex = i;
                currentmaxfragment = fragments.get(i);
            }
        }
        fragmentmarked[currentmaxindex] = false;
        for (int i = 0; i < fragments.size(); i++) {
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                PrecursorFragmentPairEdge startfrag = fragments.get(i);

                boolean groupped = false;
                for (int charge = 2; charge >= 1; charge--) {
                    float lastint = startfrag.Intensity;
                    boolean found = false;
                    for (int pkidx = 1; pkidx < 5; pkidx++) {
                        float targetmz = startfrag.FragmentMz + (float) pkidx / charge;
                        for (int j = i + 1; j < fragments.size(); j++) {
                            if (!fragmentmarked[j]) {
                                PrecursorFragmentPairEdge targetfrag = fragments.get(j);
                                if (InstrumentParameter.CalcPPM(targetfrag.FragmentMz, targetmz) < parameter.MS2PPM * (pkidx * 0.5 + 1)) {
                                    if (targetfrag.Intensity < lastint) {
                                        fragmentmarked[j] = true;
                                        lastint = targetfrag.Intensity;
                                        found = true;
                                        break;
                                    }
                                } else if (targetfrag.FragmentMz > targetmz) {
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            break;
                        }
                    }
                    if (found) {
                        groupped = true;
                        //convert to charge 1 m/z
                        startfrag.FragmentMz = startfrag.FragmentMz * charge - (charge - 1) * (float) ElementaryIon.proton.getTheoreticMass();
                        if (startfrag.FragmentMz <= ms1cluster.NeutralMass()) {
                            newfragments.add(startfrag);
                        }
                    }
                }
                if (!groupped) {
                    newfragments.add(startfrag);
                }
            }
        }
        
        fragments=newfragments;
        SortFragmentByMZ();
    }
    
    
    public void RemoveMatchedFrag(HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> MatchedFragmentMap) {
        ArrayList<PrecursorFragmentPairEdge> newlist = new ArrayList<>();
        for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
            if (!MatchedFragmentMap.containsKey(fragmentClusterUnit.PeakCurveIndexB)) {
                newlist.add(fragmentClusterUnit);
            }
        }
        fragments = newlist;
    }
    
    
    public void FilterByCriteria() {
        ArrayList<PrecursorFragmentPairEdge> newlist = new ArrayList<>();
        ArrayList<Float> CorrArrayList = new ArrayList<>();
        for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
            CorrArrayList.add(fragmentClusterUnit.Correlation);
        }
        Collections.sort(CorrArrayList);
        Collections.reverse(CorrArrayList);

        for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
            int CorrRank = 0;
            for (int intidx = 0; intidx < CorrArrayList.size(); intidx++) {
                if (CorrArrayList.get(intidx) <= fragmentClusterUnit.Correlation) {
                    CorrRank = intidx + 1;
                    break;
                }
            }
            if (fragmentClusterUnit.ComplementaryFragment || (fragmentClusterUnit.Correlation >= parameter.CorrThreshold && CorrRank <= parameter.FragmentRank && fragmentClusterUnit.FragmentMS1Rank <= parameter.PrecursorRank && fragmentClusterUnit.ApexDelta <= parameter.ApexDelta)) {
                newlist.add(fragmentClusterUnit);
            }
        }
        fragments=newlist;
    }
    
    public void BoostComplementaryIon() {
        float totalmass = (float) (ms1cluster.TargetMz() * ms1cluster.Charge - ms1cluster.Charge * ElementaryIon.proton.getTheoreticMass());
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.FALSE);
        for (int i = 0; i < fragments.size(); i++) {
            PrecursorFragmentPairEdge fragmentClusterUnit = fragments.get(i);
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                ArrayList<PrecursorFragmentPairEdge> GroupedFragments = new ArrayList<>();
                GroupedFragments.add(fragmentClusterUnit);
                float complefrag1 = (float) (totalmass - fragmentClusterUnit.FragmentMz + 2f * ElementaryIon.proton.getTheoreticMass());
                if (complefrag1 >= fragmentClusterUnit.FragmentMz) {
                    for (int j = i + 1; j < fragments.size(); j++) {
                        if (!fragmentmarked[j]) {
                            PrecursorFragmentPairEdge fragmentClusterUnit2 = fragments.get(j);
                            if (InstrumentParameter.CalcPPM(complefrag1, fragmentClusterUnit2.FragmentMz) < parameter.MS2PPM) {
                                GroupedFragments.add(fragmentClusterUnit2);
                                fragmentmarked[j] = true;
                            } else if (fragmentClusterUnit2.FragmentMz > complefrag1) {
                                break;
                            }
                        }
                    }
                }

                if (GroupedFragments.size() > 1) {
                    PrecursorFragmentPairEdge bestfragment = GroupedFragments.get(0);
                    for (PrecursorFragmentPairEdge fragment : GroupedFragments) {
                        if (fragment.Intensity > bestfragment.Intensity) {
                            bestfragment = fragment;
                        }
                    }
                    for (PrecursorFragmentPairEdge fragment : GroupedFragments) {
                        fragment.ComplementaryFragment = true;
                        fragment.Intensity = bestfragment.Intensity * growth;
                        fragment.Correlation = bestfragment.Correlation;
                        fragment.ApexDelta = bestfragment.ApexDelta;
                        fragment.RTOverlapP = bestfragment.RTOverlapP;
                        fragment.FragmentMS1Rank = bestfragment.FragmentMS1Rank;
                        fragment.FragmentMS1RankScore = bestfragment.FragmentMS1RankScore;
                    }
                }
            }
        }
    }
    
    public void IdentifyComplementaryIon(float totalmass) {
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.FALSE);
        for (int i = 0; i < fragments.size(); i++) {
            PrecursorFragmentPairEdge fragmentClusterUnit = fragments.get(i);
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                ArrayList<PrecursorFragmentPairEdge> GroupedFragments = new ArrayList<>();
                GroupedFragments.add(fragmentClusterUnit);
                float complefrag1 = (float) (totalmass - fragmentClusterUnit.FragmentMz + 2f * ElementaryIon.proton.getTheoreticMass());
                if (complefrag1 >= fragmentClusterUnit.FragmentMz) {
                    for (int j = i + 1; j < fragments.size(); j++) {
                        if (!fragmentmarked[j]) {
                            PrecursorFragmentPairEdge fragmentClusterUnit2 = fragments.get(j);
                            if (InstrumentParameter.CalcPPM(complefrag1, fragmentClusterUnit2.FragmentMz) < parameter.MS2PPM) {
                                GroupedFragments.add(fragmentClusterUnit2);
                                fragmentmarked[j] = true;
                            } else if (fragmentClusterUnit2.FragmentMz > complefrag1) {
                                break;
                            }
                        }
                    }
                }

                for (PrecursorFragmentPairEdge fragment : GroupedFragments) {
                    fragment.ComplementaryFragment = true;
                }
            }
        }
    }
    
    public XYPointCollection GetScan() {
        XYPointCollection Scan = new XYPointCollection();
        for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
            Scan.AddPointKeepMaxIfCloseValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity * fragmentClusterUnit.Correlation * fragmentClusterUnit.Correlation,parameter.MS2PPM);
        }
        return Scan;
    }
    
    private void SortFragmentByMZ() {
        Collections.sort(fragments, new Comparator<PrecursorFragmentPairEdge>() {
            @Override
            public int compare(PrecursorFragmentPairEdge o1, PrecursorFragmentPairEdge o2) {
                return Float.compare(o1.FragmentMz, o2.FragmentMz);
            }
        });
    }
    
    @Override
    public void run() {
         if (fragments.size() > 1) {             
            SortFragmentByMZ();
             if (parameter.BoostComplementaryIon) {
                 DeisotopingForPeakClusterFragment();
                 BoostComplementaryIon();
             }
            FilterByCriteria();
        }
    }

}
