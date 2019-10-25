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
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYData;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * This class implements the Continuous Wavelet Transform (CWT), Mexican Hat,
 * over raw datapoints of a certain spectrum. After get the spectrum in the
 * wavelet's time domain, we use the local maxima to detect possible peaks in
 * the original raw datapoints.
 */
public class WaveletMassDetector implements Serializable{

    /**
     * Parameters of the wavelet, NPOINTS is the number of wavelet values to use
     * The WAVELET_ESL & WAVELET_ESL indicates the Effective Support boundaries
     */
    private double NPOINTS;
    private int WAVELET_ESL = -5;
    private int WAVELET_ESR = 5;
    static boolean waveletDebug = false;
    private InstrumentParameter parameter;
    public ArrayList<XYData> DataPoint;
    double waveletWindow = 0.3;
    private double[] MEXHAT;
    double NPOINTS_half;

    public WaveletMassDetector(InstrumentParameter parameter, ArrayList<XYData> DataPoint, int NoPoints) {
        this.parameter = parameter;
        this.DataPoint = DataPoint;
        this.NPOINTS = NoPoints;

        double wstep = ((WAVELET_ESR - WAVELET_ESL) / NPOINTS);
        MEXHAT = new double[(int) NPOINTS];

        double waveletIndex = WAVELET_ESL;
        for (int j = 0; j < NPOINTS; j++) {
            // Pre calculate the values of the wavelet
            MEXHAT[j] = cwtMEXHATreal(waveletIndex, waveletWindow, 0.0);
            waveletIndex += wstep;
        }

        NPOINTS_half = NPOINTS / 2;
        d = (int) NPOINTS / (WAVELET_ESR - WAVELET_ESL);
    }
    int d;
    //public ArrayList<XYData>[] waveletCWT;
    public ArrayList<XYData>[] PeakRidge;

    public void Run() {

        //"Intensities less than this value are interpreted as noise",                
        //"Scale level",
        //"Number of wavelet'scale (coeficients) to use in m/z peak detection"
        //"Wavelet window size (%)",
        //"Size in % of wavelet window to apply in m/z peak detection");        
        int maxscale = (int) (Math.max(Math.min((DataPoint.get(DataPoint.size() - 1).getX() - DataPoint.get(0).getX()), parameter.MaxCurveRTRange), 0.5f) * parameter.NoPeakPerMin / (WAVELET_ESR + WAVELET_ESR));

        //waveletCWT = new ArrayList[15];
        PeakRidge = new ArrayList[maxscale];
        //XYData maxint = new XYData(0f, 0f);
        for (int scaleLevel = 0; scaleLevel < maxscale; scaleLevel++) {
            ArrayList<XYData> wavelet = performCWT(scaleLevel * 2 + 5);
            PeakRidge[scaleLevel] = new ArrayList<>();
            //waveletCWT[scaleLevel] = wavelet;
            XYData lastpt = wavelet.get(0);
            XYData localmax = null;
            XYData startpt = wavelet.get(0);

            boolean increasing = false;
            boolean decreasing = false;
            XYData localmaxint = null;

            for (int cwtidx = 1; cwtidx < wavelet.size(); cwtidx++) {
                XYData CurrentPoint = wavelet.get(cwtidx);
                if (CurrentPoint.getY() > lastpt.getY()) {//the peak is increasing
                    if (decreasing) {//first increasing point, last point was a possible local minimum
                        //check if the peak was symetric
                        if (localmax != null && (lastpt.getY() <= startpt.getY() || Math.abs(lastpt.getY() - startpt.getY()) / localmax.getY() < parameter.SymThreshold)) {
                            PeakRidge[scaleLevel].add(localmax);
                            localmax = CurrentPoint;
                            startpt = lastpt;
                        }
                    }
                    increasing = true;
                    decreasing = false;
                } else if (CurrentPoint.getY() < lastpt.getY()) {//peak decreasing
                    if (increasing) {//first point decreasing, last point was a possible local maximum
                        if (localmax == null || localmax.getY() < lastpt.getY()) {
                            localmax = lastpt;
                        }
                    }
                    decreasing = true;
                    increasing = false;
                }
                lastpt = CurrentPoint;
                if (localmaxint == null || CurrentPoint.getY() > localmaxint.getY()) {
                    localmaxint = CurrentPoint;
                }
                if (cwtidx == wavelet.size() - 1 && decreasing) {
                    if (localmax != null && (CurrentPoint.getY() <= startpt.getY() || Math.abs(CurrentPoint.getY() - startpt.getY()) / localmax.getY() < parameter.SymThreshold)) {
                        PeakRidge[scaleLevel].add(localmax);
                    }
                }
            }

            if (!waveletDebug) {
                wavelet.clear();
                //wavelet = null;
            }
        }
    }

    /**
     * Perform the CWT over raw data points in the selected scale level
     *
     *
     */
    private ArrayList<XYData> performCWT(int scaleLevel) {
        int length = DataPoint.size();
        ArrayList<XYData> cwtDataPoints = new ArrayList<XYData>();

        int a_esl = scaleLevel * WAVELET_ESL;
        int a_esr = scaleLevel * WAVELET_ESR;
        double sqrtScaleLevel = Math.sqrt(scaleLevel);
        for (int dx = 0; dx < length; dx++) {
            /*
             * Compute wavelet boundaries
             */
            int t1 = a_esl + dx;
            if (t1 < 0) {
                t1 = 0;
            }
            int t2 = a_esr + dx;
            if (t2 >= length) {
                t2 = (length - 1);
            }

            /*
             * Perform convolution
             */
            float intensity = 0f;
            for (int i = t1; i <= t2; i++) {
                int ind = (int) (NPOINTS_half) + ((int) d * (i - dx) / scaleLevel);
                if (ind < 0) {
                    ind = 0;
                }
                if (ind >= NPOINTS) {
                    ind = (int) NPOINTS - 1;
                }
//                if(i<0 || ind<0){
//                    System.out.print("");
//                }
                intensity += DataPoint.get(i).getY() * MEXHAT[ind];
            }
            intensity /= sqrtScaleLevel;
            // Eliminate the negative part of the wavelet map
            if (intensity < 0) {
                intensity = 0;
            }
            cwtDataPoints.add(new XYData(DataPoint.get(dx).getX(), intensity));
        }
        return cwtDataPoints;
    }

    /**
     * This function calculates the wavelets's coefficients in Time domain
     *
     * @param double x Step of the wavelet
     * @param double a Window Width of the wavelet
     * @param double b Offset from the center of the peak
     */
    private double cwtMEXHATreal(double x, double window, double b) {
        /*
         * c = 2 / ( sqrt(3) * pi^(1/4) )
         */
        double c = 0.8673250705840776;
        double TINY = 1E-200;
        double x2;

        if (window == 0.0) {
            window = TINY;
        }
        //x-b=t
        //window=delta
        x = (x - b) / window;
        x2 = x * x;
        return c * (1.0 - x2) * Math.exp(-x2 / 2);
    }
    /**
     * This function searches for maximums from wavelet data points
     */
}
