package org.apidb.ggtools.array;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.lang.Math;
import java.util.Comparator;
import java.util.Hashtable;
import java.io.*;
import java.lang.Exception;
import java.lang.NullPointerException;
import java.util.*;
import java.lang.reflect.Array;
import java.util.Date;

/*
  args[0] = input filename
  args[1] = output filename peaks
  args[2] = output filename smoothed
  args[3] = mult of SD cutoff
  args[4] = num consecutive probes
  args[5] = peak max gap (don't use neighboring probes this far apart or more in the same feature)
  args[6] = smoothing max gap (don't smooth over things separated by this much or more)
*/

public class ChIP_Chip_Peak_Finder {
    public static void main(String[] args) {
	try {
	    theApp = new ChIP_Chip_Peak_Finder();
	    boolean inspect = false;
	    File dataFile = null;
	    if(Array.getLength(args) == 1) {
		dataFile = new File(args[0]);
		inspect = true;
	    }
	    if((Array.getLength(args) < 7 && Array.getLength(args) > 1) || Array.getLength(args) == 0) {
		theApp.printUsage();
		System.exit(0);
	    }
	    boolean output_smoothed = true;
	    long genomelength = 0;
	    String PFN = "";
	    String SFN = "";
	    File outputFile = null;
	    float mult_of_SD_cutoff = 0;
	    int num_consecutive_probes = 0;
	    int feature_maxgap = 0;
	    int smoother_maxgap = 0;
	    boolean gbrowse = false;
	    if(Array.getLength(args) >= 7) {
		if(args[2].equals(false)) {
		    output_smoothed = false;
		} else {
		    SFN = args[2];
		}
		dataFile = new File(args[0]);
		PFN = args[1];
		outputFile = new File(PFN);
		mult_of_SD_cutoff = Float.parseFloat(args[3]);
		try {
		    num_consecutive_probes = Integer.parseInt(args[4]);
		} catch (Exception e) {
		    System.err.println("Error: Argument four must be an integer\n");
		    System.exit(1);
		}
		try {
		    feature_maxgap = Integer.parseInt(args[5]);
		} catch (Exception e) {
		    System.err.println("Error: Argument five must be an integer\n");
		    System.exit(1);
		}
		try {
		    smoother_maxgap = Integer.parseInt(args[6]);
		} catch (Exception e) {
		    System.err.println("Error: Argument six must be an integer\n");
		    System.exit(1);
		}
		int numargs = Array.getLength(args);
		for(int i=7; i<numargs; i++) {
		    if(args[i].equals("-gbrowse"))
			gbrowse = true;
		    if(args[i].equals("-inspect"))
			inspect = true;
		}
	    }

	    BufferedWriter out = null;
	    if(!inspect) {
		out = new BufferedWriter(new FileWriter(outputFile));
	    }
	    System.err.println("");
	    System.err.println("Assessing data file...");
	    Object[] fileSpecs = theApp.assessFile(dataFile);
	    int numChromosomes = (Integer) (fileSpecs[0]);
	    int maxProbesOneSpan = (Integer) (fileSpecs[1]);
	    int numProbes = (Integer) (fileSpecs[2]);
	    Hashtable<String, Integer> contigs = (Hashtable) (fileSpecs[3]);
	    String[] names = (String[]) (fileSpecs[4]);
	    boolean header = (Boolean) (fileSpecs[5]);
	    if(!inspect) {
		out.write("# input file = " + args[0] + "\n");
		out.write("# multiple of SD cutoff = " + mult_of_SD_cutoff + "\n");
		out.write("# feature max gap = " + feature_maxgap + "\n");
		out.write("# smoothing max gap = " + smoother_maxgap + "\n");
		out.write("# num consecutive probes cutoff = " + num_consecutive_probes + "\n");
		out.write("# num chromosomes = " + numChromosomes + "\n");
		String mm = theApp.format((long)maxProbesOneSpan);
		out.write("# max probes on one chromosome in the data file = " + mm + "\n");
		mm = theApp.format((long)numProbes);
		out.write("# num probes in data file = " + mm + "\n");
		out.flush();
	    }
	    else {
		System.out.println("input file = " + args[0]);
		System.out.println("multiple of SD cutoff = " + mult_of_SD_cutoff);
		System.out.println("feature max gap = " + feature_maxgap);
		System.out.println("smoothing max gap = " + smoother_maxgap);
		System.out.println("num consecutive probes cutoff = " + num_consecutive_probes);
		System.out.println("num chromosomes = " + numChromosomes);
		String mm = theApp.format((long)maxProbesOneSpan);
		System.out.println("max probes on one chromosome in the data file = " + mm);
		mm = theApp.format((long)numProbes);
		System.out.println("num probes in data file = " + mm);
	    }
	    System.err.println("");
	    System.err.println("Reading the data file...");
	    Object[] readfileOutput = theApp.readFile(dataFile, numChromosomes, maxProbesOneSpan, numProbes, contigs, header);
	    int[][] locationArray = (int[][])readfileOutput[0];
	    float[][] logratioArray = (float[][])readfileOutput[1];
	    int[] readCnt = (int[])readfileOutput[2];
	    genomelength = (Long)readfileOutput[3];
	    long[] chrMin = ((long[])readfileOutput[4]);
	    long[] chrMax = ((long[])readfileOutput[5]);
	    String mm = theApp.format((long)maxProbesOneSpan);
	    if(!inspect) {
		out.write("# max probes on one chromosome = " + mm + "\n");
	    }
	    else {
		System.out.println("max probes on one chromosome = " + mm);
	    }
	    mm = theApp.format((long)genomelength);
	    if(!inspect) {
		out.write("# genomelength = " + mm + "\n");
	    }
	    else {
		System.out.println("genomelength = " + mm);
	    }
	    float[] distances_between_probes_histogram = new float[51];
	    int MAX = 0;
	    for(int chr=0; chr<numChromosomes; chr++) {
		for(int i=0; i<readCnt[chr]-1; i++) {
		    if(locationArray[chr][i+1] - locationArray[chr][i] > 5000) {
			distances_between_probes_histogram[50]++;
		    }
		    else {
			int bin = Math.round((locationArray[chr][i+1] - locationArray[chr][i]) / 100);
			if(locationArray[chr][i+1] - locationArray[chr][i] > MAX) {
			    MAX = (locationArray[chr][i+1] - locationArray[chr][i]);
			}
			distances_between_probes_histogram[bin]++;
		    }
		}
	    }
	    String[] labels = new String[51];
	    labels[50] = "> 5,000";
	    for(int i=0; i<50; i++) {
		int s = i * 100;
		int e = s + 100;
		labels[i] = s + "-" + e;
	    }
	    if(!inspect) {
		out.write("#\n# breakdown of distances between probes:\n");
		for(int i=0; i<51; i++) {
		    if(distances_between_probes_histogram[i] > 0)
			out.write("# " + labels[i] + ": " + distances_between_probes_histogram[i] + "\n");
		}
		out.write("#");
	    }
	    else {
		System.out.println("\nbreakdown of distances between probes:");
		for(int i=0; i<51; i++) {
		    if(distances_between_probes_histogram[i] > 0)
			System.out.println(labels[i] + ": " + distances_between_probes_histogram[i]);
		}
		System.exit(0);
	    }

	    /*
	    String graph = theApp.TextBarGraph(distances_between_probes_histogram, labels, 30);
	    System.out.println("\n" + graph);
	    */

	    System.err.println("\nSmoothing...");
	    for(int chr=0; chr<numChromosomes; chr++) {
		logratioArray[chr] = theApp.Smooth(logratioArray[chr], locationArray[chr], readCnt[chr], smoother_maxgap);
	    }

	    if(output_smoothed) {
		File outputFile3 = new File(SFN);
		BufferedWriter out3 = new BufferedWriter(new FileWriter(outputFile3));
		if(gbrowse) {
		    out3.write("[ChIP CHIP]\nglyph=xyplot\ngraph_type=points\npoint_symbol=disc\nfgcolor=blackbgcolor=red\nheight=100\nmin_score=-2\nmax_score=3\nlabel=1\nkey=ChIP CHIP\n");
		}
		for(int chr=0; chr<numChromosomes; chr++) {
		    if(gbrowse) {
			out3.write("\nreference=" + names[chr] + "\n");
		    }
		    for(int j=0; j<readCnt[chr]; j++) {
			int start = locationArray[chr][j];
			int end = locationArray[chr][j]+1;
			if(gbrowse) {
			    out3.write("\"ChIP CHIP\"\t\"Log Ratio\"\t" + start + "-" + end + "\tscore=" + logratioArray[chr][j] + "\n");
			}
			else {
			    out3.write(names[chr] + "\t" + start + "-" + end + "\t" + logratioArray[chr][j] + "\n");
			}
		    }
		}
		out3.close();
	    }

	    double SD = theApp.computeSD(logratioArray, readCnt, numChromosomes);
	    double mean = theApp.computeMEAN(logratioArray, readCnt, numChromosomes);
	    double X = SD * mult_of_SD_cutoff + mean;
	    out.write("# SD = " + SD + "\n");
	    out.write("# mean = " + mean + "\n");
	    out.write("# Absolute Cutoff = " + X + "\n");
	    out.write("# PEAKS START HERE\n");
	    if(!gbrowse) {
		out.write("# chr\tstart\tend\tnum.probes\tscore\n");
	    }
	    out.flush();
	    int start = 0;
	    int end = 0;
	    int cnt_consec_probes_exceeding_cutoff = 0;
	    String outstr = "";
	    boolean flag = false;
	    int cnt=0;
	    for(int chr=0; chr<numChromosomes; chr++) {
		cnt_consec_probes_exceeding_cutoff = 0;
		double peakscore = 0;
		outstr = "";
		flag = false;
		for(int j=0; j<readCnt[chr]; j++) {
		    if(logratioArray[chr][j] >= mult_of_SD_cutoff * SD + mean && j < readCnt[chr]-1 && (j==0 || locationArray[chr][j] - locationArray[chr][j-1] <= feature_maxgap)) {
			cnt_consec_probes_exceeding_cutoff++;
			peakscore = peakscore + (logratioArray[chr][j] - mean) / SD;
			if(cnt_consec_probes_exceeding_cutoff == 1) {
			    start = locationArray[chr][j];
			}
			end = locationArray[chr][j];
		    }
		    else {
			if(logratioArray[chr][j] >= mult_of_SD_cutoff * SD + mean && j == readCnt[chr]-1) {
			    end = locationArray[chr][j];
			    cnt_consec_probes_exceeding_cutoff++;
			    peakscore = peakscore + (logratioArray[chr][j] - mean) / SD;
			}
			if(cnt_consec_probes_exceeding_cutoff >= num_consecutive_probes) {
			    String s = mult_of_SD_cutoff + "_" + num_consecutive_probes;
			    s = s.replace(".","_");
			    if(gbrowse) {
				outstr = outstr + names[chr] + "\t" + s + "\t" + "peaks" + "\t" + start + "\t" + end + "\t" + cnt_consec_probes_exceeding_cutoff + "\t" + "." + "\t" + "." + "\t" + cnt + "\n";
			    }
			    else {
				peakscore = (double) Math.round(peakscore * 10000) / (double) 10000;
				outstr = outstr + names[chr] + "\t" + start + "\t" + end + "\t" + cnt_consec_probes_exceeding_cutoff + "\t" + peakscore + "\n";
			    }
			    cnt++;
			    flag = true;
			}
			cnt_consec_probes_exceeding_cutoff = 0;
			peakscore = 0;
		    }
		}
		if(flag) {
		    out.write(outstr);
		}
	    }	    
	    out.flush();
	} catch (Exception e) { 
	    theApp.printUsage();
	    e.printStackTrace(System.err);
            System.exit(0);
	}
    }

    public double computeSD(float[][] logratioArray, int[] readCnt, int numChromosomes) {
	double mean = 0;
	int cnt = 0;
	for(int chr=0; chr<numChromosomes; chr++) {
	    for(int i=0; i<readCnt[chr]; i++) {
		mean = mean + logratioArray[chr][i];
		cnt++;
	    }
	}
	mean = mean / (double) cnt;
	double sd = 0;
	for(int chr=0; chr<numChromosomes; chr++) {
	    for(int i=0; i<readCnt[chr]; i++) {
		sd = sd + (mean - logratioArray[chr][i]) * (mean - logratioArray[chr][i]);
	    }
	}
	sd = sd / (cnt - 1);
	sd = Math.sqrt(sd);

	return sd;
    }

    public double computeMEAN(float[][] logratioArray, int[] readCnt, int numChromosomes) {
	double mean = 0;
	int cnt = 0;
	for(int chr=0; chr<numChromosomes; chr++) {
	    for(int i=0; i<readCnt[chr]; i++) {
		mean = mean + logratioArray[chr][i];
		cnt++;
	    }
	}
	mean = mean / (double) cnt;

	return mean;
    }

    private Object[] readFile(File file, int numChromosomes, int maxProbesOneSpan, int numProbes, Hashtable<String, Integer> contigs, boolean header) {
        FileInputStream inFile = null;
	int[][] locationArray = new int[numChromosomes][maxProbesOneSpan];
	float[][] logratioArray = new float[numChromosomes][maxProbesOneSpan];
	System.out.println("numChromosomes = " + numChromosomes);
	System.out.println("maxProbesOneSpan = " + maxProbesOneSpan);
	int maxReadLength = 0;
	int[] readCnt = new int[numChromosomes];
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The file " + file + " does not seem to exist...");
            System.exit(1);
        }
	long[] chrMin = new long[numChromosomes];
	long[] chrMax = new long[numChromosomes];
	for(int i=0; i<numChromosomes; i++) {
	    chrMin[i] = 100000000;
	    chrMax[i] = 0;
	}
	for(int i=0; i<numChromosomes; i++) {
	    readCnt[i] = 0;
	}
	int flag = 0;
        try {
            inFile = new FileInputStream(file);
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
            String line = "";
	    if(header)
		line = br.readLine();
            while((line = br.readLine()) != null) {
		String[] temp = line.split("\\s+", 0);
		locationArray[contigs.get(temp[0])][readCnt[contigs.get(temp[0])]] = Integer.parseInt(temp[1]);
		logratioArray[contigs.get(temp[0])][readCnt[contigs.get(temp[0])]] = Float.parseFloat(temp[2]);
		if(locationArray[contigs.get(temp[0])][readCnt[contigs.get(temp[0])]] < chrMin[contigs.get(temp[0])])
		    chrMin[contigs.get(temp[0])] = locationArray[contigs.get(temp[0])][readCnt[contigs.get(temp[0])]];
		if(locationArray[contigs.get(temp[0])][readCnt[contigs.get(temp[0])]] > chrMax[contigs.get(temp[0])]) 
		    chrMax[contigs.get(temp[0])] = locationArray[contigs.get(temp[0])][readCnt[contigs.get(temp[0])]];
		readCnt[contigs.get(temp[0])]++;
	    }
            inFile.close();
        } catch(IOException e) {
            e.printStackTrace(System.err);
            System.exit(1);
        }
	for(int chr=0; chr<numChromosomes; chr++) {
	    float[][] temparray = new float[readCnt[chr]][2];
	    for(int i=0; i<readCnt[chr]; i++) {
		temparray[i][0] = locationArray[chr][i];
		temparray[i][1] = logratioArray[chr][i];
	    }
	    Arrays.sort(temparray, theApp.FirstIndexComparator);
	    for(int i=0; i<readCnt[chr]; i++) {
		locationArray[chr][i] = (int)temparray[i][0];
		logratioArray[chr][i] = temparray[i][1];
	    }
	}
	long genomelength = 0;
	for(int i=0; i<numChromosomes; i++) {
	    if(chrMax[i] >= chrMin[i])
		genomelength = genomelength + chrMax[i] - chrMin[i] + 1;
	}

	Object[] RETURN = new Object[6];
	RETURN[0] = locationArray;
	RETURN[1] = logratioArray;
	RETURN[2] = readCnt;
	RETURN[3] = genomelength;
	RETURN[4] = chrMin;
	RETURN[5] = chrMax;
	return RETURN;
    }

    private Object[] assessFile(File file) {
        FileInputStream inFile = null;
        int x = 0;
	int chrNum = 0;
	int numChromosomes = 0;
	Object[] answer = new Object[6];
        /*
	int[] probeCnt = new int[500];
        */
	int maxProbes = 0;
	int numProbes = 0;
	String line = "";
	String firstline = "";
	boolean header = false;
        /*
	for(int i=0; i<500; i++) {
	    probeCnt[i] = 0;
	}
        */
        try {
            inFile = new FileInputStream(file);
        } catch(FileNotFoundException e) {
	    System.err.println("The file " + file + " does not seem to exist...");
            System.exit(1);
        }
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
	    firstline = br.readLine();
	    if(firstline == null) {
		System.err.println("This file is empty...");
		System.exit(1);
	    }
	    String[] temp = firstline.split("\\s+", 0);	
	    int n = Array.getLength(temp);
	    if(n < 3) {
		System.err.println("\nError: the file must have three tab delimited columns, check line one.\n\n");
		System.exit(0);
	    }
	    try {
		Float.parseFloat(temp[1]);
		Float.parseFloat(temp[2]);
	    } catch(Exception e5) {
		header = true;
	    }
            inFile.close();
        } catch(IOException e2) {
	    e2.printStackTrace(System.err);
            System.exit(1);
        }
	int linecounter = 0;
	Hashtable<String, Integer> contigs = new Hashtable<String, Integer>();

        try {
	    inFile = new FileInputStream(file);
	    BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
	    if(header) {
		line = br.readLine();
		linecounter++;
	    }
	    int cntr = 0;
            while((line = br.readLine()) != null) {
		String[] temp = line.split("\\s+", 0);
		if(!contigs.containsKey(temp[0])) {
		    contigs.put(temp[0], new Integer(cntr));
		    cntr++;
		}
	    }
	} catch(IOException e3) {
	    e3.printStackTrace(System.err);
	    System.exit(1);
	}
	Enumeration e = contigs.keys();
	int cnt=0;
	while( e.hasMoreElements() ){
	    e.nextElement();
	    cnt++;
	}
	numChromosomes = cnt;
	String[] names = new String[cnt];
	Enumeration e2 = contigs.keys();
	while( e2.hasMoreElements() ){
	    String tmp = (String) e2.nextElement();
	    names[contigs.get(tmp)] = tmp;
	}

	int[] probeCnt = new int[contigs.size()];

	try {
            inFile = new FileInputStream(file);
            BufferedReader br = new BufferedReader(new InputStreamReader(inFile));
	    if(header) {
		line = br.readLine();
		linecounter++;
	    }
	    int flag = 0;
            while((line = br.readLine()) != null) {
		linecounter++;
		flag = 0;
		String[] temp = line.split("\\s+", 0);
		numProbes++;
		probeCnt[contigs.get(temp[0])]++;
		if(probeCnt[contigs.get(temp[0])] > maxProbes) {
		    maxProbes = probeCnt[contigs.get(temp[0])];
		}
	    }
	    inFile.close();
	} catch(IOException e4) {
	    e4.printStackTrace(System.err);
	    System.exit(1);
	}

	answer[0] = numChromosomes;
	answer[1] = maxProbes;
	answer[2] = numProbes;
	answer[3] = contigs;
	answer[4] = names;
	answer[5] = header;
	return answer;
    }

    public String format (long num) {
        String s = Long.toString(num);
        String answer = "";
        String[] temp = s.split("", 0);

        int len = Array.getLength(temp);
        int count = 0;
        for(int i=len - 1; i>0; i--) {
            answer = temp[i] + answer;
            count++;
            if(count == 3 && i>1) {
                answer = "," + answer;
                count = 0;
            }
        }
        return answer;
    }

    public float[] Smooth (float[] input_vector, int[] location_vector, int readCnt, int max_gap) {

	int w = max_gap;
	int vector_length = readCnt;
	float[] ls = new float[0];
	float b = 0;
	float m = 0;

	if(vector_length<6) {
	    return input_vector;
	}
	float[] temp1 = new float[3];
	float[] temp2 = new float[3];
	float[] smoothed_vector = new float[vector_length];


	//      for the cases notes below:
	//	let d(j,k) be distance between probes j and k
	//	let w = smoother_maxgap

	//  the case 0
	smoothed_vector[0] = input_vector[0];
	if(location_vector[1] - location_vector[0] < w && location_vector[2] - location_vector[1] < w) {
	    temp1 = new float[3];
	    temp2 = new float[3];
	    int shift = location_vector[0] - 1;
	    temp1[0]=input_vector[0];
	    temp1[1]=input_vector[1];
	    temp1[2]=input_vector[2];
	    temp2[0]=location_vector[0] - shift;
	    temp2[1]=location_vector[1] - shift;
	    temp2[2]=location_vector[2] - shift;
	    
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[0] = m * (location_vector[0] - shift) + b;
	}
	else {
	    smoothed_vector[0] = input_vector[0];
	}
	if(((Float) smoothed_vector[0]).isNaN() || ((Float) smoothed_vector[0]).isInfinite()) {
	    smoothed_vector[0] = input_vector[0];
	}
	//  the case n (n=vector length)
	smoothed_vector[vector_length-1] = input_vector[vector_length-1];
	if(location_vector[vector_length-1] - location_vector[vector_length-2] < w && location_vector[vector_length-2] - location_vector[vector_length-3] < w) {
	    temp1 = new float[3];
	    temp2 = new float[3];
	    int shift = location_vector[vector_length-3] - 1;
	    temp1[0]=input_vector[vector_length-3];
	    temp1[1]=input_vector[vector_length-2];
	    temp1[2]=input_vector[vector_length-1];
	    temp2[0]=location_vector[vector_length-3] - shift;
	    temp2[1]=location_vector[vector_length-2] - shift;
	    temp2[2]=location_vector[vector_length-1] - shift;
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[vector_length-1] = m*(location_vector[vector_length-1] - shift) + b;
	}
	if(((Float) smoothed_vector[vector_length-1]).isNaN() || ((Float) smoothed_vector[vector_length-1]).isInfinite()) {
	    smoothed_vector[vector_length-1] = input_vector[vector_length-1];
	}

	//  the case 1

	/*
	Cases:
	   1) d(0,1) <= w and d(1,2) <= w and d(2,3) <= w
	       - smooth using 0,1,2,3
	   2) d(0,1) > w and d(1,2) <= w and d(2,3) <= w
	       - smooth using 1,2,3
	   3) d(0,1) <=w and d(1,2) <= w and d(2,3) > w
	       - smooth using 0,1,2
	   4) all other cases don't smooth
	*/
	smoothed_vector[1] = input_vector[1];
	if(location_vector[1] - location_vector[0] <= w && location_vector[2] - location_vector[1] <= w && location_vector[3] - location_vector[2] <= w) {
	    temp1 = new float[4];
	    temp2 = new float[4];
	    int shift = location_vector[0] - 1;
	    temp1[0]=input_vector[0];
	    temp1[1]=input_vector[1];
	    temp1[2]=input_vector[2];
	    temp1[3]=input_vector[3];
	    temp2[0]=location_vector[0] - shift;
	    temp2[1]=location_vector[1] - shift;
	    temp2[2]=location_vector[2] - shift;
	    temp2[3]=location_vector[3] - shift;

	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[1] = m*(location_vector[1] - shift) + b;
	}
	if(location_vector[1] - location_vector[0] > w && location_vector[2] - location_vector[1] <= w && location_vector[3] - location_vector[2] <= w) {
	    temp1 = new float[3];
	    temp2 = new float[3];
	    int shift = location_vector[1] - 1;
	    temp1[0]=input_vector[1];
	    temp1[1]=input_vector[2];
	    temp1[2]=input_vector[3];
	    temp2[0]=location_vector[1] - shift;
	    temp2[1]=location_vector[2] - shift;
	    temp2[2]=location_vector[3] - shift;
	    
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[1] = m*(location_vector[1] - shift) + b;
	}
	if(location_vector[1] - location_vector[0] <= w && location_vector[2] - location_vector[1] <= w && location_vector[3] - location_vector[2] > w) {
	    temp1 = new float[3];
	    temp2 = new float[3];

	    temp1[0]=input_vector[0];
	    temp1[1]=input_vector[1];
	    int shift = location_vector[0] - 1;
	    temp1[2]=input_vector[2];
	    temp2[0]=location_vector[0] - shift;
	    temp2[1]=location_vector[1] - shift;
	    temp2[2]=location_vector[2] - shift;
	    
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[1] = m*(location_vector[1] - shift) + b;
	}
	if(((Float) smoothed_vector[1]).isNaN() || ((Float) smoothed_vector[1]).isInfinite()) {
	    smoothed_vector[1] = input_vector[1];
	}

	//  the case n-1 (n=vector length)
	/*
	Cases:
	   1) d(n-3,n-2) <= w and d(n-2,n-1) <= w and d(n-1,n) <= w
	       - smooth using n-3, n-2, n-1, n
	   2) d(n-3,n-2) > w and d(n-2,n-1) <= w and d(n-1,n) <= w
	       - smooth using n-2, n-1, n
	   3) d(n-3,n-2) <=w and d(n-2,n-1) <= w and d(n-1,n) > w
	       - smooth using n-3, n-2, n-1
	   4) all other cases don't smooth
	*/

	smoothed_vector[vector_length-2] = input_vector[vector_length-2];
	if(location_vector[vector_length-3] - location_vector[vector_length-4] <= w && location_vector[vector_length-2] - location_vector[vector_length-3] <= w && location_vector[vector_length-1] - location_vector[vector_length-2] <= w) {
	    temp1 = new float[4];
	    temp2 = new float[4];
	    int shift = location_vector[vector_length-4] - 1;
	    temp1[0]=input_vector[vector_length-4];
	    temp1[1]=input_vector[vector_length-3];
	    temp1[2]=input_vector[vector_length-2];
	    temp1[3]=input_vector[vector_length-1];
	    temp2[0]=location_vector[vector_length-4] - shift;
	    temp2[1]=location_vector[vector_length-3] - shift;
	    temp2[2]=location_vector[vector_length-2] - shift;
	    temp2[3]=location_vector[vector_length-1] - shift;
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[vector_length-2]=m*(location_vector[vector_length-2] - shift) + b;
	}
	if(location_vector[vector_length-3] - location_vector[vector_length-4] > w && location_vector[vector_length-2] - location_vector[vector_length-3] <= w && location_vector[vector_length-1] - location_vector[vector_length-2] <= w) {
	    temp1 = new float[3];
	    temp2 = new float[3];
	    int shift = location_vector[vector_length-3] - 1;
	    temp1[0]=input_vector[vector_length-3];
	    temp1[1]=input_vector[vector_length-2];
	    temp1[2]=input_vector[vector_length-1];
	    temp2[0]=location_vector[vector_length-3] - shift;
	    temp2[1]=location_vector[vector_length-2] - shift;
	    temp2[2]=location_vector[vector_length-1] - shift;
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[vector_length-2]=m*(location_vector[vector_length-2] - shift) + b;
	}
	if(location_vector[vector_length-3] - location_vector[vector_length-4] <= w && location_vector[vector_length-2] - location_vector[vector_length-3] <= w && location_vector[vector_length-1] - location_vector[vector_length-2] > w) {
	    temp1 = new float[3];
	    temp2 = new float[3];
	    int shift = location_vector[vector_length-4] - 1;
	    temp1[0]=input_vector[vector_length-4];
	    temp1[1]=input_vector[vector_length-3];
	    temp1[2]=input_vector[vector_length-2];
	    temp2[0]=location_vector[vector_length-4] - shift;
	    temp2[1]=location_vector[vector_length-3] - shift;
	    temp2[2]=location_vector[vector_length-2] - shift;
	    ls = LS(temp2, temp1);
	    b = ls[0];
	    m = ls[1];
	    smoothed_vector[vector_length-2]=m*(location_vector[vector_length-2] - shift) + b;
	}
	if(((Float) smoothed_vector[vector_length-2]).isNaN() || ((Float) smoothed_vector[vector_length-2]).isInfinite()) {
	    smoothed_vector[vector_length-2] = input_vector[vector_length-2];
	}

	// the rest of the cases

	/*

	Cases:
	    1) max distance between probes <= w
	           - smooth using i-2, i-1, i, i+1, i+2
	    2) d(i-1,i) > w
	       a) d(i,i+1) > w or d(i+1,i+2) > w
	           - do not smooth point i
	       b) if d(i,i+1) <=w and d(i+1, i+2) <= w
	           - smooth using i, i+1, i+2
	    3) d(i,i+1) > w
	       a) d(i-2,i-1) > w or d(i-1,i) > w
	           - do not smooth point i
	       b) d(i-2,i-1) <= w and d(i-1,i) <= w
	           - smooth using i-2, i-1, i
	    4) d(i-2,i-1) <= w and d(i-1,i) <= w and d(i,i+1) <= w and d(i+1,i+2) > w
	           - smooth using i-2, i-1, i, i+1
	    5) d(i-2,i-1) > w and d(i-1,i) <= w and d(i,i+1) <= w and d(i+1,i+2) <= w
	           - smooth using i-1, i, i+1, i+2


	*/

	for(int i=2; i<vector_length-2; i++) {
	    smoothed_vector[i] = input_vector[i];
	    if(location_vector[i-1] - location_vector[i-2] <= w && location_vector[i] - location_vector[i-1] <= w && location_vector[i+1] - location_vector[i] <= w && location_vector[i+2] - location_vector[i+1] <=w) {
		temp1 = new float[5];
		temp2 = new float[5];
		int shift = location_vector[i-2] - 1;
		temp1[0]=input_vector[i-2];
		temp1[1]=input_vector[i-1];
		temp1[2]=input_vector[i];
		temp1[3]=input_vector[i+1];
		temp1[4]=input_vector[i+2];
		temp2[0]=location_vector[i-2] - shift;
		temp2[1]=location_vector[i-1] - shift;
		temp2[2]=location_vector[i] - shift;
		temp2[3]=location_vector[i+1] - shift;
		temp2[4]=location_vector[i+2] - shift;
		ls = LS(temp2, temp1);
		b = ls[0];
		m = ls[1];
		smoothed_vector[i]=m*(location_vector[i] - shift) + b;
	    }
	    if(location_vector[i] - location_vector[i-1] > w) {
		if(location_vector[i+1] - location_vector[i] > w || location_vector[i+2] - location_vector[i+1] > w) {
		    smoothed_vector[i] = input_vector[i];
		}
		else {
		    temp1 = new float[3];
		    temp2 = new float[3];
		    int shift = location_vector[i] - 1;
		    temp1[0]=input_vector[i];
		    temp1[1]=input_vector[i+1];
		    temp1[2]=input_vector[i+2];
		    temp2[0]=location_vector[i] - shift;
		    temp2[1]=location_vector[i+1] - shift;
		    temp2[2]=location_vector[i+2] - shift;
		    ls = LS(temp2, temp1);
		    b = ls[0];
		    m = ls[1];
		    smoothed_vector[i]=m*(location_vector[i] - shift) + b;
		}
	    }
	    if(location_vector[i+1] - location_vector[i] > w) {
		if(location_vector[i-1] - location_vector[i-2] > w || location_vector[i-1] - location_vector[i] > w) {
		    smoothed_vector[i] = input_vector[i];
		}
		else {
		    temp1 = new float[3];
		    temp2 = new float[3];
		    int shift = location_vector[i-2] - 1;
		    temp1[0]=input_vector[i-2];
		    temp1[1]=input_vector[i-1];
		    temp1[2]=input_vector[i];
		    temp2[0]=location_vector[i-2] - shift;
		    temp2[1]=location_vector[i-1] - shift;
		    temp2[2]=location_vector[i] - shift;
		    ls = LS(temp2, temp1);
		    b = ls[0];
		    m = ls[1];
		    smoothed_vector[i]=m*(location_vector[i] - shift) + b;
		}
	    }
	    if(location_vector[i-1] - location_vector[i-2] <= w && location_vector[i] - location_vector[i-1] <= w && location_vector[i+1] - location_vector[i] <= w && location_vector[i+2] - location_vector[i+1] > w) {
		    temp1 = new float[4];
		    temp2 = new float[4];
		    int shift = location_vector[i-2] - 1;
		    temp1[0]=input_vector[i-2];
		    temp1[1]=input_vector[i-1];
		    temp1[2]=input_vector[i];
		    temp1[3]=input_vector[i+1];
		    temp2[0]=location_vector[i-2] - shift;
		    temp2[1]=location_vector[i-1] - shift;
		    temp2[2]=location_vector[i] - shift;
		    temp2[3]=location_vector[i+1] - shift;
		    ls = LS(temp2, temp1);
		    b = ls[0];
		    m = ls[1];
		    smoothed_vector[i]=m*(location_vector[i] - shift) + b;
	    }
	    if(location_vector[i-1] - location_vector[i-2] > w && location_vector[i] - location_vector[i-1] <= w && location_vector[i+1] - location_vector[i] <= w && location_vector[i+2] - location_vector[i+1] <= w) {
		    temp1 = new float[4];
		    temp2 = new float[4];
		    int shift = location_vector[i-1] - 1;
		    temp1[0]=input_vector[i-1];
		    temp1[1]=input_vector[i];
		    temp1[2]=input_vector[i+1];
		    temp1[3]=input_vector[i+2];
		    temp2[0]=location_vector[i-1] - shift;
		    temp2[1]=location_vector[i] - shift;
		    temp2[2]=location_vector[i+1] - shift;
		    temp2[3]=location_vector[i+2] - shift;
		    ls = LS(temp2, temp1);
		    b = ls[0];
		    m = ls[1];
		    smoothed_vector[i]=m*(location_vector[i] - shift) + b;
	    }
	    if(((Float) smoothed_vector[i]).isNaN() || ((Float) smoothed_vector[i]).isInfinite()) {
		smoothed_vector[i] = input_vector[i];
	    }
	}

	return smoothed_vector;
    }    
    
    public float[] LS(float[] X_vect, float[] Y_vect) {
	
	int X_vect_length=Array.getLength(X_vect);
	int Y_vect_length=Array.getLength(Y_vect);
	
	if(X_vect_length != Y_vect_length) {
	    System.err.println("ERROR: tried to run the LS subroutine with vectors of differing lengths\n");
	}
	
	float y=0;
	float x=0;
	float yy=0;
	float xx=0;
	float xy=0;
	
	for(int i=0; i<X_vect_length;i++) {
	    x=x+X_vect[i];
	    y=y+Y_vect[i];
	    xx=xx+(X_vect[i])*(X_vect[i]);
	    yy=yy+(Y_vect[i])*(Y_vect[i]);
	    xy=xy+(X_vect[i])*(Y_vect[i]);
	}
	float[] ls = new float[2];
	
	ls[0]=(y*xx-x*xy)/(X_vect_length*xx-x*x); // y-intercept
	ls[1]=(X_vect_length*xy-x*y)/(X_vect_length*xx-x*x);  // slope
	
	return ls;
    }

    public String TextBarGraph(float[] vector, String[] labels, int yResolution) {
	int num_bars = Array.getLength(vector);
	float max = vector[0];
	for(int i=0; i<num_bars; i++) {
	    if(vector[i] > max)
		max = vector[i];
	}
	String graph = "";
	graph = graph + max + "\n";
	float step = max / yResolution;
	for(int y=yResolution; y>=0; y--) {
	    graph = graph + "|";
	    for(int x=0; x<num_bars; x++) {
		int bin = Math.round(vector[x] / step);
		if(bin > y) {
		    graph = graph + "*";
		}
		else {
		    graph = graph + " ";
		}
	    }
	    graph = graph + "\n";
	}
	graph = graph + " ";
	for(int x=0; x<num_bars; x++) {
	    graph = graph + "-";
	}
	graph = graph + "\n ";
	for(int x=0; x<num_bars; x++) {
	    int xx = x+1;
	    if(xx < 10) 
		graph = graph + xx;
	    if(xx >= 10 && xx < 100) {
		int f = Math.round(xx / 10);
		graph = graph + f;
	    }
	    if(xx >= 100 && xx < 1000) {
		int f = Math.round(xx / 100);
		graph = graph + f;
	    }
	    if(xx >= 1000 && xx < 10000) {
		int f = Math.round(xx / 1000);
		graph = graph + f;
	    }
	}
	if(num_bars > 10) {
	    graph = graph + "\n ";
	    for(int x=0; x<num_bars; x++) {
		int xx = x+1;
		if(xx < 10) 
		    graph = graph + " ";
		if(xx >= 10 && xx < 100) {
		    int f = Math.round(xx % 10);
		    graph = graph + f;
		}
		if(xx >= 100 && xx < 1000) {
		    int f = Math.round(xx % 100);
		    graph = graph + f;
		}
		if(xx >= 1000 && xx < 10000) {
		    int f = Math.round(xx % 1000);
		    graph = graph + f;
		}
	    }
	}
	if(num_bars > 100) {
	    graph = graph + "\n ";
	    for(int x=0; x<num_bars; x++) {
		int xx = x+1;
		if(xx < 100) 
		    graph = graph + "  ";
		if(xx >= 100 && xx < 1000) {
		    int f = Math.round(xx % 10);
		    graph = graph + f;
		}
		if(xx >= 1000 && xx < 10000) {
		    int f = Math.round(xx % 100);
		    graph = graph + f;
		}
	    }
	}
	if(num_bars > 1000) {
	    graph = graph + "\n ";
	    for(int x=0; x<num_bars; x++) {
		int xx = x+1;
		if(xx < 1000) 
		    graph = graph + "  ";
		if(xx >= 1000 && xx < 10000) {
		    int f = Math.round(xx % 10);
		    graph = graph + f;
		}
	    }
	}
	graph = graph + "\n\n";
	for(int x=0; x<num_bars; x++) {
	    int xx = x+1;
	    graph = graph + xx + ": " + labels[x] + "\n";
	}
	graph = graph + "\n";
	return graph;
    }

    public void printUsage() {
	    System.err.println("\n========================================USAGE============================================\n\njava -jar ChIP_Chip_Peak_Finder.jar <input file> <peaks output filename> <smoothed data output filename> <SD mult cutoff> <num consec probes> <feature max gap> <smoother max gap> [options]\n\n\n* Input file should have three tab delimited columns: 'chr', 'location', 'log ratio'\n - chr can be any string (without tabs).\n - location is an integer representing the position of the probe on the chromosome.\n - There can be a header, which also must have three tab delimited columns.\n    (Header is auto-detected as long as it doesn't look like data.)\n\nPeaks are called if <num consec probes> neighboring probes all exceed\nM + <SD mult cutoff> * SD, where M is the mean log ratio (over the entire\narray) and SD is the standard deviation of the log ratios.\n\n<feature max gap> is an integer, neighboring probes separated by more\n                  than this amount are not used in the same feature\n\n<smoother max gap> is an interger, neighboring probes separated by more\n                   than this amount are not used to smooth each other\n\nSmoothed data will be written to <smoothed data output filename>, to skip\nwriting smoothed data put 'false' for this argument.\n\n[options] (options come last):\n -gbrowse : Output files will be in gbrowse uploadable format\n -inspect : Output basic array info to stdout, including breakdown of probe distances.\n            (Alternatively: run with just one argument: the data file name.)\n\nNOTE: if you get an 'out of memory' error, try running with the option -Xmx1000m\n(put this right after the word 'java')\n==========================================================================================\n\n");
    }

    public static Comparator FirstIndexComparator = new Comparator() {
	    public int compare(Object array1, Object array2) {
		int r=0;
		if(((float[])array1)[0] - ((float[])array2)[0] < 0)
		    r = -1;
		if(((float[])array1)[0] - ((float[])array2)[0] > 0)
		    r = 1;
		if(((float[])array1)[0] - ((float[])array2)[0] == 0)
		    r = 0;
		return r;
	    }
    };
    
    public static ChIP_Chip_Peak_Finder theApp;
}
