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
package MSUmpire.spectrumparser;

import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIA_Setting implements Serializable{
    private static final long serialVersionUID = 646984181894L;
    public TreeMap<XYData, ArrayList<Integer>> DIAWindows = new TreeMap<XYData, ArrayList<Integer>>();
    public TreeMap<XYData, ArrayList<Integer>> MS1Windows = new TreeMap<XYData, ArrayList<Integer>>();
    public float F_DIA_WindowSize = 25;
    public SpectralDataType.DataType dataType;

    private synchronized void writeObject(java.io.ObjectOutputStream stream) throws java.io.IOException {        
        stream.defaultWriteObject();
        stream.writeInt(DIAWindows.size());
        for(XYData window : DIAWindows.keySet()){
            stream.writeObject(window);
            stream.writeObject(DIAWindows.get(window));
        }        
        stream.writeFloat(F_DIA_WindowSize);
        stream.writeObject(dataType);
    }
    private void readObject(java.io.ObjectInputStream in) throws ClassNotFoundException, IOException {
        in.defaultReadObject();
        DIAWindows=new TreeMap<>();        
        int size=in.readInt();
        for (int i = 0; i < size; i++) {            
            XYData window = (XYData) in.readObject();
            ArrayList<Integer> scans=(ArrayList<Integer>)in.readObject();
            DIAWindows.put(window, scans);
        }
        F_DIA_WindowSize=in.readInt();
        dataType=(SpectralDataType.DataType)in.readObject();
    }
    
    
    public void WriteParamSerialization(String mzXMLFileName) {
        try {
            Logger.getRootLogger().info("Writing DIA setting to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_diasetting.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_diasetting.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(this);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public static DIA_Setting ReadDIASettingSerialization(String filepath) {
        if (!new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_diasetting.ser").exists()) {
            return null;
        }
        try {
            Logger.getRootLogger().debug("Reading DIA setting from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_diasetting.ser...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_diasetting.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            DIA_Setting setting = (DIA_Setting) in.readObject();
            in.close();
            fileIn.close();
            return setting;

        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;        
        }
    }

}
