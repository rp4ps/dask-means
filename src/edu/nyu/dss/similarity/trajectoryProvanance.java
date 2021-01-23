package edu.nyu.dss.similarity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.junit.experimental.theories.Theories;
import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kpaths.Util;

public class trajectoryProvanance {
	
	/*
	 * create a set of dataset with map matching, sampling, simplication (ratio), enrichment (ratio), outlier removal (far points), 
	 * segmentation (ratio)
	 * Input: trajectory dataset, output: a new file with multiple operations, and combinations (multiple).
	 */
	static void cleanDataset(Map<Integer, double [][]> trajectoryDataset, Map<Integer, double [][]> mappedtrajectoryDataset, int datasetThreshold, String folder) {
		int i=1;
		String fileString="";
		int sampingRatio = 5;
		for(int trajectoryid: trajectoryDataset.keySet()) {
			int fileno = (int)i/datasetThreshold;
			String fielnameString = folder+fileno+".txt";
			String mappingfielnameString = folder+ fileno+"mm.txt";
			String samplefielnameString = folder+fileno+"sample.txt";
			String simplynameString = folder+fileno+"remove.txt";
			String enrichnameString = folder+fileno+"enrich.txt";
			
			fileString="";
			double [][]trajectory = trajectoryDataset.get(trajectoryid);
			for(int j=0; j<trajectory.length; j++)
				fileString += "-"+trajectory[j][0]+","+trajectory[j][1]+",";
			Util.write(fielnameString, fileString+"\n");
			Util.write(mappingfielnameString, mapmatching(mappedtrajectoryDataset.get(trajectoryid))+"\n");
			Util.write(simplynameString, simplifying(trajectory, 2)+"\n");
			Util.write(enrichnameString, enrichment(trajectory, 2)+"\n");
			if((int)i%datasetThreshold%sampingRatio==0)
				Util.write(samplefielnameString, fileString+"\n");
			// conduct the cleaning and storing into multiple files
			// store label in the title
			i++;
		}
	}
	
	
	static String mapmatching(double [][]trajectory) {
		String fileString = "";
		for(int j=0; j<trajectory.length; j++)
			fileString += trajectory[j][0]+","+trajectory[j][1]+",";
		return fileString;
	}
	
	static String simplifying(double [][]trajectory, int ratio) {
		String fileString = "";
		for(int j=0; j<trajectory.length; j++) {
			if(j%ratio==0)
				fileString += "-"+trajectory[j][0]+","+trajectory[j][1]+",";
		}
		return fileString;
	}
	
	static String enrichment(double [][]trajectory, int ratio) {
		String fileString = "";
		for(int j=0; j<trajectory.length; j++) {
			fileString += "-"+trajectory[j][0]+","+trajectory[j][1]+",";
			if(j==trajectory.length-1)
				continue;
			double rangex = trajectory[j][0] - trajectory[j+1][0];
			double rangey = trajectory[j][1] - trajectory[j+1][1];
			for(int i=1; i<ratio; i++) {
				double x = trajectory[j][0]+rangex/ratio*i;
				double y = trajectory[j][1]+rangey/ratio*i;
				fileString += "-"+x+","+y+",";
			}
		}
		return fileString;
	}
	
	/**/
	void outlier() {
		//run outlier removal
	}
	
	void segmentation() {
		
	}
	
	void clean_raw() {
		
	}
	
	void convertVertex(String vertexfile) {
		
	}
	
	static void convertEdge(String edgefile, String edgeTrajectory, String writeedge, String writevertex) throws FileNotFoundException, IOException {
		File file = new File(edgefile);
		Map<Integer, String> edgefileMap = new HashMap<Integer, String>();
		Map<Integer, String> vertexfileMap = new HashMap<Integer, String>();
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split(";");
				String[] edgeidStrings = splitString[1].split(",");
				String[] edgeidStrings1 = splitString[2].split(",");
				String edgeString = "";
				System.out.println(splitString[1]);
				for(int i =0; i<edgeidStrings.length-1; i++) {
					edgeString += edgeidStrings1[i]+","+edgeidStrings[i]+",";
					System.out.println(edgeidStrings1[i]+","+edgeidStrings[i]);
				}
				String vertexString = edgeidStrings1[0]+","+edgeidStrings[0];
				edgefileMap.put(Integer.valueOf(splitString[0]), edgeString);
				vertexfileMap.put(Integer.valueOf(splitString[0]), vertexString);
			}
		}
		file = new File(edgeTrajectory);
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split("\t");
				String[] edgeidStrings = splitString[1].split(",");
				String edgeString = "";
				String vertexString = "";
				for(int i = 0; i<edgeidStrings.length; i++) {
					int edgeid = Integer.valueOf(edgeidStrings[i]);
					edgeString += edgefileMap.get(edgeid);
					vertexString += vertexfileMap.get(edgeid)+",";
				}
				Util.write(writeedge, edgeString+"\n");
				Util.write(writevertex, vertexString+"\n");
			}
		}
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		File dataFile = new File(args[0]);
		File mappdataFile = new File(args[1]);
		
		String edgeFile = args[3];
		String edgedatasetFile = args[4];
		String writeedgedatasetFile = args[5];
		String writenodedatasetFile = args[6];
	//	convertEdge(edgeFile, edgedatasetFile, writeedgedatasetFile, writenodedatasetFile);
		mappdataFile = new File(writenodedatasetFile);
		Map<Integer, double[][]> trajectory = Framework.readPorto(dataFile);
		Map<Integer, double[][]> mapped = Framework.readPorto(mappdataFile);
		cleanDataset(trajectory, mapped, 1000, args[2]);
	}
}
