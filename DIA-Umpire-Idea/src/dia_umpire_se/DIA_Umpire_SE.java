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
package dia_umpire_se;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.DIA.DIAPack;
import Utility.ConsoleLogger;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_Umpire_SE {

    /**
     * @param args the command line arguments DIA_Umpire parameterfile
     */

    public static void main(String[] args) throws InterruptedException, FileNotFoundException, ExecutionException, IOException, ParserConfigurationException, DataFormatException, SAXException, SQLException, Exception {
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire singal extraction analysis (version: v1.284, 2015.01)");
        /*
        if (args.length != 2) {
            System.out.println("command format error, it should be like: java –jar –Xmx8G DIA_Umpire_SE.jar mzMXL_file diaumpire.se_params");
            return;
        }
        */

         try {
           ConsoleLogger.SetConsoleLogger(Level.INFO);
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_debug.log");
        } catch (Exception e) {
        }
        

//        String parameterfile=args[1];
//        String mzXMLPath=args[0];
        String parameterfile = "E:\\DIA2018\\Tutorial3_DIAUmpire\\parameter_files\\diaumpire_se_Thermo_params.txt";
        String mzXMLPath = "E:\\DIA2018\\Tutorial3_DIAUmpire\\collinsb_X1803_171.mzXML";
        Logger.getRootLogger().info("Parameter file:" + parameterfile);        
        BufferedReader reader = new BufferedReader(new FileReader(parameterfile));
        
        String line = "";
        InstrumentParameter para = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600);
        para.DetermineBGByID = false;
        para.EstimateBG=true;
        int NoCPUs = 2;

        SpectralDataType.DataType dataType = SpectralDataType.DataType.DIA_F_Window;
        String WindowType = "";
        int WindowSize = 25;

        ArrayList<XYData> WindowList = new ArrayList<>();

        boolean ExportPrecursorPeak = false;
        boolean ExportFragmentPeak = false;
        
        //<editor-fold defaultstate="collapsed" desc="Read parameter file">
        while ((line = reader.readLine()) != null) {
            if (!"".equals(line) && !line.startsWith("#")) {
                //System.out.println(line);
                if (line.equals("==window setting begin")) {
                    while (!(line = reader.readLine()).equals("==window setting end")) {
                        if (!"".equals(line)) {
                            WindowList.add(new XYData(Float.parseFloat(line.split("\t")[0]), Float.parseFloat(line.split("\t")[1])));
                        }
                    }
                    continue;
                }
                if (line.split("=").length < 2) {
                    continue;
                }
                String type = line.split("=")[0].trim();
                if(type.startsWith("para.")){
                    type=type.replace("para.","SE.");
                }
                String value = line.split("=")[1].trim();
                switch (type) {
                    case "Thread": {
                        NoCPUs = Integer.parseInt(value);
                        break;
                    }
                    case "ExportPrecursorPeak": {
                        ExportPrecursorPeak = Boolean.parseBoolean(value);
                        break;
                    }
                    case "ExportFragmentPeak": {
                        ExportFragmentPeak = Boolean.parseBoolean(value);
                        break;
                    }
                    
                    //<editor-fold defaultstate="collapsed" desc="instrument parameters">
                    case "RPmax": {
                        para.PrecursorRank = Integer.parseInt(value);
                        break;
                    }
                    case "RFmax": {
                        para.FragmentRank = Integer.parseInt(value);
                        break;
                    }
                    case "CorrThreshold": {
                        para.CorrThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "DeltaApex": {
                        para.ApexDelta = Float.parseFloat(value);
                        break;
                    }                    
                    case "BoostComplementaryIon":{
                        para.BoostComplementaryIon=Boolean.parseBoolean(value);
                        break;
                    }
                    case "AdjustFragIntensity":{
                        para.AdjustFragIntensity=Boolean.parseBoolean(value);
                        break;
                    }                    
                    case "SE.MS1PPM": {
                        para.MS1PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2PPM": {
                        para.MS2PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.SN": {
                        para.SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2SN": {
                        para.MS2SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSIntensity": {
                        para.MinMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSMSIntensity": {
                        para.MinMSMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinRTRange": {
                        para.MinRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MaxNoPeakCluster": {
                        para.MaxNoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinNoPeakCluster": {
                        para.MinNoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinMS2NoPeakCluster": {
                        para.MinMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MaxCurveRTRange": {
                        para.MaxCurveRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.Resolution": {
                        para.Resolution = Integer.parseInt(value);
                        break;
                    }
                    case "SE.RTtol": {
                        para.RTtol = Float.parseFloat(value);
                        break;
                    }
                    case "SE.NoPeakPerMin": {
                        para.NoPeakPerMin = Integer.parseInt(value);
                        break;
                    }
                    case "SE.StartCharge": {
                        para.StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.EndCharge": {
                        para.EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2StartCharge": {
                        para.MS2StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2EndCharge": {
                        para.MS2EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.NoMissedScan": {
                        para.NoMissedScan = Integer.parseInt(value);
                        break;
                    }
                    case "SE.Denoise": {
                        para.Denoise = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.EstimateBG": {
                        para.EstimateBG = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaks": {
                        para.RemoveGroupedPeaks = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.MinFrag":{
                        para.MinFrag = Integer.parseInt(value);
                        break;
                    }
//</editor-fold>
                    case "WindowType": {
                        WindowType = value;
                        switch (WindowType) {
                            case "SWATH": {
                                dataType = SpectralDataType.DataType.DIA_F_Window;
                                break;
                            }
                            case "V_SWATH": {
                                dataType = SpectralDataType.DataType.DIA_V_Window;
                                break;
                            }
                            case "MSX": {
                                dataType = SpectralDataType.DataType.MSX;
                                break;
                            }
                            case "MSE": {
                                dataType = SpectralDataType.DataType.MSe;
                                break;
                            }
                        }
                        break;
                    }
                    case "WindowSize": {
                        WindowSize = Integer.parseInt(value);
                        break;
                    }
                }
            }
        }
//</editor-fold>

        try {            
            File mzxml = new File(mzXMLPath);
            if (mzxml.exists()) {
                long time = System.currentTimeMillis();
                Logger.getRootLogger().info("=================================================================================================");
                Logger.getRootLogger().info("Processing " + mzXMLPath + "....");
                DIAPack DiaFile = new DIAPack(mzXMLPath, NoCPUs);
                DiaFile.SetDataType(dataType);
                DiaFile.SetParameter(para);
                if (dataType == SpectralDataType.DataType.DIA_F_Window) {
                    DiaFile.SetWindowSize(WindowSize);
                } else if (dataType == SpectralDataType.DataType.DIA_V_Window) {
                    for (XYData window : WindowList) {
                        DiaFile.AddVaribleWindow(window);
                    }
                }
                DiaFile.SaveDIASetting();
                DiaFile.SaveParams();     
                DiaFile.ExportPrecursorPeak = ExportPrecursorPeak;
                DiaFile.ExportFragmentPeak = ExportFragmentPeak;
                Logger.getRootLogger().info("Module A: Signal extraction");
                DiaFile.process();
                time = System.currentTimeMillis() - time;
                Logger.getRootLogger().info(mzXMLPath + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            }
            Logger.getRootLogger().info("Job complete");
            Logger.getRootLogger().info("=================================================================================================");

        } catch (Exception e) {
            System.out.print(e.getMessage());
            throw e;
        }
    }

}
