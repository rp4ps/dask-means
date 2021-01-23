package au.edu.rmit.trajectory.clustering.kmeans;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kpaths.Util;

public class fairkmeans {
	static void readDataset(String rawdata, String poi, String newfile) throws IOException {
		File file = new File(poi);
		Map<String, String> poiMap = new HashMap<String, String>();
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split("\t");
				String countryString = splitString[4];
				if(countryString.equals("US"))// just present the results in US, we can also shrink the range
					poiMap.put(splitString[0], splitString[1]+","+splitString[2]);
			}
		}
		file = new File(rawdata);
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split("\t");
				if(poiMap.containsKey(splitString[1])) {
					Util.write(newfile, splitString[0]+","+poiMap.get(splitString[1])+"\n");
				}
			}
		}
	}
	
	/*
	 * preprocess the data, generate the normal k-means, group all the people's point,
	 */
	static void FairProcessData(double [][]originData, int userID[], Map<Integer, Integer> userNumber, int dim, String file) {
		Map<Integer, double[]> summapMap = new HashMap<Integer, double[]>();
		for(int i=0; i<originData.length; i++) {
			double []point = originData[i];
			int id = userID[i];
			double sum[] = new double[dim];
			if(summapMap.containsKey(id)) {
				sum = summapMap.get(id);
			}
			for(int j=0; j<dim; j++) {
				sum[j] += point[j];
			}
			summapMap.put(id, sum);
		}
		
		for(int user: summapMap.keySet()) {
			double sum[] = summapMap.get(user);
			String content = "";
			for(int j=0; j<dim; j++) {
				content += Double.toString(sum[j]/userNumber.get(user))+",";
			}
			Util.write(file, content+"\n");
		}
	}
	
	
	/*
	 * we study the quality of generated clusters, distance vs frequency, show whether they are related in the normal k-means.
	 * 
	 * we design several metrics to prove it is good.
	 * 
	 * whether points in each cluster are euqality accessible to the centroids.
	 */
	void effectivenessStudy(ArrayList<cluster> CENTERSEuc, double[][] dataOriginal) {
		// we study the number of people in each group
		
		// we need to give multiple metric to reflect we give more weight to the under-represented person.
		
		// over-representation
		
		// the top-n closers whether are the frequent people.
		
		// whether every one can find the site within a given threshold.
	
	}
	
	
	/*
	 * plot the distance-to-centroid, frequency, to observe,
	 * 
	 * we need to put every cluster's into one, we evaluate
	 */
	void corelations() {
		// we show the normal k-means' and we can observe that 
		// compute the distance
	}
	
	
	/*
	 * find the top-k for each centroid, and observe their frequency shown in the list are balanced.
	 */
	void topkCentorid() {
		
	}
	
	/*
	 * conduct the range query, and see the the range distributions on various.
	 */
	void topkPoint() {
		
	}
							
	// we verify our efficiency mainly, and compare with a baseline through preprocessing to observe the difference.
			
	// We can foresee that normal method will lead to unfair planning, all go to dense areas, check-in, taxi trips dataset all have this bias.

	
	public static void main(String[] args) throws IOException {
		readDataset("/Users/sw160/Documents/datasets/dataset_WWW2019/raw_Checkins_anonymized.txt", "/Users/sw160/Documents/datasets/dataset_WWW2019/raw_POIs.txt", 
				"/Users/sw160/Documents/datasets/dataset_WWW2019/clustering.txt");
	}
}
