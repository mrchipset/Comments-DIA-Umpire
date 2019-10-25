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
package MSUmpire.FragmentLib;

import java.util.Random;

/**
 * This class contains a random generator for protein sequences.
 * 
 * @author Bram Minnaert
 */

public class RandomSequenceGeneratorWoKPR {
    
    /**
     * All possible characters
     */
    private static final char[] CHARS = {
    	'A', 'N', 'D', 'C', 'Q', 'E',
		'G', 'H', 'I', 'L', 'M', 'F',
		'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X'};    
    
    /**
     * Number of possible characters
     */
    private static final int NUMBER_OF_CHARS = CHARS.length;
    
    /**
     * Random generator
     */
    private static Random random = new Random();
    
    /**
     * Returns random sequence
     * @param length Size of the sequence
     * @return Random sequence
     */
    public static String generate(int length) {
        StringBuffer buffer = new StringBuffer();
        char randomChar;
        int  randomInt;
        for (int i = 0; i < length ; i++) {
            randomInt = random.nextInt(NUMBER_OF_CHARS);
            randomChar = CHARS[randomInt];
            buffer.append(randomChar);
        }
        return buffer.toString();
    }

    /**
     * Displays 10 random protein sequences with length 50.
     * @param args no args
     */
    public static void main(String[] args) {
        for (int i = 0; i < 10; i++) {
        	System.out.println("S" + i + " = " + generate(50));
        }
    }
}