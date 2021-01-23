package edu.nyu.dss.similarity;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.channels.NonWritableChannelException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.stream.Stream;

import org.apache.commons.lang3.tuple.Pair;
import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kmeans.indexAlgorithm;
import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;
import java_cup.internal_error;
import org.apache.commons.io.FileUtils;


public class Framework {
	static String edgeString, nodeString, distanceString, fetureString, indexString="";
	static String aString = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
	static Map<Integer, Integer> argoDataMap = new HashMap<Integer, Integer>();
	
	/*
	 * to do: 1) indexing road network, and can augment the edge information for each point.
	 * 2) 
	 */
	
	
	/*
	 * information to compare all the algorithms
	 */
	int numberofComparison = 17;//all the comparison
	double time[]= new double[numberofComparison];// the overall running time
	static int fileNo = 1;
	
	/*
	 * read file and store all the points, check whether its csv, and extract locations only.
	 */
	void readFileLocationsAuto() {

	}
	
	/*
	 * the file name indicates the dimension
	 */
	void reafFileByFileNames(String filename) {
		
	}
	
	/*
	 * the file name indicates the dimension, we divide it into tables based on the poi types.
	 * we have multiple cities poi dataset, and store them in one.
	 * we can do a case study based on this.
	 */
	void readPOIdataset(String filename) {
		
	}
	
	/*
	 * read the shapenet dataset
	 */
	static void readShapeNet(File folder, Map<Integer, double[][]> datasetMap, int limit) throws IOException {
		File[] fileNames = folder.listFiles();
        for(File file : fileNames){
            if(file.isDirectory()){
            	//read the files inside
            	readShapeNet(file, datasetMap, limit);
            }else{
                if(file.getName().contains(".pts")) {
					System.out.println(file.getName());
					long lineNumber = 0;
					try (Stream<String> lines = Files.lines(file.toPath())) {
						lineNumber = lines.count();
					}
					double[][] a = new double[(int) lineNumber][];
					int i=0;
					try (BufferedReader br = new BufferedReader(new FileReader(file))) {
						String strLine;
						while ((strLine = br.readLine()) != null) {
							// System.out.println(strLine);
							String[] splitString = strLine.split(" ");
							a[i] = new double[3];
							a[i][0] = Double.valueOf(splitString[0]);
							a[i][1] = Double.valueOf(splitString[1]);
							a[i][2] = Double.valueOf(splitString[2]);
							i++;
						}
					}
					datasetMap.put(fileNo, a);
					fileNo++;
					if(fileNo>limit)
						break;
                }    
            }
        }
	}
	
	/*
	 * read the trajectory database
	 */
	static Map<Integer, double[][]> readPorto(File file) throws IOException {
		Map<Integer, double[][]> datasetMap = new HashMap<Integer, double[][]>();
		System.out.println("read file " + file.getCanonicalPath());
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			int line = 1; 
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split(",");
				double [][]a = new double[splitString.length/2][];
				for (int i = 0; i < splitString.length; i++) {
					if ((i + 1) % 2 == 0) {
						a[i/2] = new double[2];
						for(int j=i-1; j<i+1;j++) {
							a[i/2][j-i+1] = Double.valueOf(splitString[j]);
						}
					}
				}
				datasetMap.put(line, a);
				line++;
			}
		}
		return datasetMap;
	}

	/*
	 * read a folder and extract the corresponding column of each file inside
	 */
	static Map<Integer, double[][]> readFolder(File folder, int limit) {
		File[] fileNames = folder.listFiles();
		int fileNo = 1;
		Map<Integer, double[][]> datasetMap = new HashMap<Integer, double[][]>();
        for(File file : fileNames){
            if(file.isDirectory()){
            	readFolder(file, limit);
            }else{
                try {
            		String a = file.getName();
            		argoDataMap.put(fileNo, Integer.valueOf(a.substring(0, a.length()-4)));
                    readContent(file, fileNo++, datasetMap);
                } catch (IOException e) {
                    e.printStackTrace();
                }
       
            }
            if(fileNo> limit)
            	break;
        }
        return datasetMap;
	}

	public static void readContent(File file, int fileNo,Map<Integer, double[][]> datasetMap) throws IOException {
		System.out.println("read file " + file.getCanonicalPath());
		long lineNumber = 0;
		try (Stream<String> lines = Files.lines(file.toPath())) {
			lineNumber = lines.count();
		}
		double[][] a = new double[(int) lineNumber-1][];
		int i=0;
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				// System.out.println(strLine);
				String[] splitString = strLine.split(",");
				String aString = splitString[3];
				if (aString.matches("-?\\d+(\\.\\d+)?")) {
					a[i] = new double[2];
					a[i][0] = Double.valueOf(splitString[3]);
					a[i][1] = Double.valueOf(splitString[4]);
					i++;
				}
			}
		}
		datasetMap.put(fileNo, a);
	}
	
	/*
	 * read from the mnist dataset and comparing with the dataset
	 */
	public static Map<Integer, double[][]> readMnist(File file) throws IOException {
		Map<Integer, double[][]> datasetMap = new HashMap<Integer, double[][]>();
		System.out.println("read file " + file.getCanonicalPath());
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			int line = 1; 
			while ((strLine = br.readLine()) != null) {
				// System.out.println(strLine);
				String[] splitString = strLine.split(" ");
			//	System.out.println(strLine);
				double [][]a = new double[28][];
				for (int i = 0; i < splitString.length; i++) {
					if ((i + 1) % 28 == 0) {
						a[i/28] = new double[28];
						for(int j=i-27; j<i+1;j++) {
						//	System.out.println(j-i+27);
							a[i/28][j-i+27] = Double.valueOf(splitString[j]);
						}
					}
				}
				datasetMap.put(line, a);
				line++;
			}
		}
		return datasetMap;
	}
	
	/*
	 * read from the mnist dataset and comparing with the dataset
	 */
	public static Map<Integer, double[][]> readMnist2D(File file) throws IOException {
		Map<Integer, double[][]> datasetMap = new HashMap<Integer, double[][]>();
		System.out.println("read file " + file.getCanonicalPath());
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			int line = 1; 
			while ((strLine = br.readLine()) != null) {
				// System.out.println(strLine);
				String[] splitString = strLine.split(" ");
			//	System.out.println(strLine);
				double [][]a = new double[28][];
				for (int i = 0; i < splitString.length; i++) {
					if ((i + 1) % 28 == 0) {
						a[i/28] = new double[28];
						for(int j=i-27; j<i+1;j++) {
						//	System.out.println(j-i+27);
							a[i/28][j-i+27] = Double.valueOf(splitString[j]);
						}
					}
				}
				datasetMap.put(line, a);
				line++;
			}
		}
		return datasetMap;
	}
	
	
	
	public void testTriangle() throws IOException {
		String foldername = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
		File folder = new File(foldername);
		Map<Integer, double[][]> dataMap = readFolder(folder, 1000);
		folder = new File("/Users/sw160/Downloads/torch-clus/dataset/minist_60000_784");
		dataMap = readMnist(folder);
		for(int a:dataMap.keySet()) {
			if(a>dataMap.size()-3)
				break;
		//	System.out.println(dataMap.get(a)[0].length);
			double distance = Hausdorff.HausdorffExact(dataMap.get(1), dataMap.get(a), dataMap.get(a)[0].length, false);
		//	System.out.println(distance);
			double distance1 = Hausdorff.HausdorffExact(dataMap.get(a), dataMap.get(a+1), dataMap.get(a)[0].length,false);
			double distance2 = Hausdorff.HausdorffExact(dataMap.get(a+1), dataMap.get(a+2), dataMap.get(a)[0].length,false);
			double distance3 = Hausdorff.HausdorffExact(dataMap.get(a), dataMap.get(a+2), dataMap.get(a)[0].length,false);
			if(distance1+distance2<distance3)
				System.out.println("no triangle inequality");//
		}
	}
	
	public static void testTriangleFrechet() throws IOException {
		String foldername = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
		File folder = new File(foldername);
		Map<Integer, double[][]> dataMap = readFolder(folder, 1000);
		folder = new File("/Users/sw160/Downloads/torch-clus/dataset/minist_60000_784");
		dataMap = readMnist(folder);
		for(int a:dataMap.keySet()) {
			if(a>dataMap.size()-3)
				break;
		//	System.out.println(dataMap.get(a)[0].length);
			double distance = Hausdorff.Frechet(dataMap.get(1), dataMap.get(a), dataMap.get(a)[0].length);
		//	System.out.println(distance);
			double distance1 = Hausdorff.Frechet(dataMap.get(a), dataMap.get(a+1), dataMap.get(a)[0].length);
			double distance2 = Hausdorff.Frechet(dataMap.get(a+1), dataMap.get(a+2), dataMap.get(a)[0].length);
			double distance3 = Hausdorff.Frechet(dataMap.get(a), dataMap.get(a+2), dataMap.get(a)[0].length);
			if(distance1+distance2<distance3 || Math.abs(distance1-distance2)>distance3)
				System.out.println("no triangle inequality");//
		}
	}
	
	public static void generateGTforPrediciton(String distanceString, Map<Integer, double[][]> dataMapPorto, Map<Integer, indexNode> indexMap, int dimension, int limit) {
		Util.write(distanceString, "datasetID,dataset2ID,distance,normalizedDistance\n");
		for (int b : dataMapPorto.keySet()) {// we can build index to accelerate
			if (b > limit)
				break;
			for (int a : dataMapPorto.keySet()) {// we can build index to accelerate
				if (a > limit)
					break;
		//		double distance = Hausdorff.HausdorffExact(dataMapPorto.get(b), dataMapPorto.get(a),
		//				dataMapPorto.get(a)[0].length, false);
				
				double distance = AdvancedHausdorff.IncrementalDistanceDirected(dataMapPorto.get(b), dataMapPorto.get(a), dimension, indexMap.get(b), indexMap.get(a), 0, 0, 0.05, null, null);

				// we should change it to our fast one, instead of this slow
				// building index for each dataset
				double pivot_distance = Util.EuclideanDis(indexMap.get(b).getPivot(), indexMap.get(a).getPivot(),
						dimension);
				double max_distance = pivot_distance + indexMap.get(b).getRadius() + indexMap.get(a).getRadius();
				Util.write(distanceString,
						Integer.toString(b) + "," + Integer.toString(a) + "," + Double.toString(distance) + ","
								+ Double.toString((max_distance - distance) / max_distance) + "\n");
			}
		}
	}

	public static void SerializedMappingTable(String file) {
		try {
			FileOutputStream fos = new FileOutputStream(file);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(argoDataMap);
			oos.close();
			fos.close();
			System.out.printf("Serialized HashMap data is saved in hashmap.ser");
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	/*
	 * put the files into a folder
	 */
	public static void UnionFile(int datasetid, String destination) throws FileNotFoundException, IOException {
		File file = new File(aString+"/"+datasetid+".csv");
		File dest = new File(destination+"/"+datasetid+".csv");
		try {
		    FileUtils.copyFile(file, dest);
		} catch (IOException e) {
		    e.printStackTrace();
		}
	}
	
	/*
	 * this reads a single dataset, instead of storing all the dataset in the main memory
	 */
	public static double[][] readSingleFile(int datasetid) throws FileNotFoundException, IOException {
		File file = new File(aString+"/"+datasetid+".csv");
		long lineNumber = 0;
		try (Stream<String> lines = Files.lines(file.toPath())) {
			lineNumber = lines.count();
		}
		double[][] a = new double[(int) lineNumber-1][];
		int i=0;
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				// System.out.println(strLine);
				String[] splitString = strLine.split(",");
				String aString = splitString[3];
				if (aString.matches("-?\\d+(\\.\\d+)?")) {
					a[i] = new double[2];
					a[i][0] = Double.valueOf(splitString[3]);
					a[i][1] = Double.valueOf(splitString[4]);
					i++;
				}
			}
		}
		return a;
	}
	
	/*
	 * normalize the dataset accroding to the latitude and longitude, so we can have a fair Hausdorff and join
	 */
	public void normalization() {
		// get the range first, then normalize by a range
	}
	
	public static void main(String[] args) throws IOException {		
		
		aString = args[0];
		int dimension = 2;
		int capacity = 10;
		int limit = 10010;
		File folder = new File(aString);
		Map<Integer, double[][]> dataMapPorto = new HashMap<Integer, double[][]>();
		if(aString.contains("argo")) {
			dataMapPorto = readFolder(folder, limit);
			edgeString = "./index/dss/argo_edge.txt";
			nodeString = "./index/dss/argo_node.txt";
			capacity = 30;
			limit = 10000;
			distanceString = "./index/dss/argo_haus_distance.txt";
			fetureString = "./index/dss/argo_haus_features.txt";
			indexString = "./index/dss/index/argo/";
			SerializedMappingTable(indexString+"mapping.ser");
		}
		else if(aString.contains("mnist")){
			dataMapPorto = readMnist(folder);
			dimension = 28;
			edgeString = "./index/dss/mnist_edge.txt";
			nodeString = "./index/dss/mnist_node.txt";
			distanceString = "./index/dss/mnist_haus_distance.txt";
			fetureString = "./index/dss/mnist_haus_features.txt";
		}else if(aString.contains("beijing")){
			dataMapPorto = readPorto(folder);
			edgeString = "./index/dss/beijing_edge.txt";
			nodeString = "./index/dss/beijing_node.txt";
			distanceString = "./index/dss/beijing_haus_distance.txt";
			fetureString = "./index/dss/beijing_haus_features.txt";
			indexString = "./index/dss/index/beijing/";
		}else if(aString.contains("footprint")){
			dataMapPorto = readPorto(folder);
			edgeString = "./index/dss/nyc_edge.txt";
			nodeString = "./index/dss/nyc_node.txt";
			distanceString = "./index/dss/nyc_haus_distance.txt";
			fetureString = "./index/dss/nyc_haus_features.txt";
		}else if(aString.contains("Shapenet")){
			dimension = 3;
			limit = 100;
			readShapeNet(folder, dataMapPorto, limit+100);
			edgeString = "./index/dss/shapenet_edge.txt";
			nodeString = "./index/dss/shapenet_node.txt";
			distanceString = "./index/dss/shapenet_haus_distance.txt";
			fetureString = "./index/dss/shapenet_haus_features.txt";
			indexString = "./index/dss/index/shapenet/";
			System.out.println(fileNo);
			
		}else {
			dataMapPorto = readPorto(folder);
			edgeString = "./index/dss/porto_edge.txt";
			nodeString = "./index/dss/porto_node.txt";
		//	limit = 1000;
			distanceString = "./index/dss/porto_haus_distance.txt";
			fetureString = "./index/dss/porto_haus_features.txt";
			indexString = "./index/dss/index/porto/";
		}
	
		
		Map<Integer, indexNode> indexMap = new HashMap<Integer, indexNode>();
		Map<Integer, double[]> featureMap = new HashMap<Integer, double[]>();
		ArrayList<indexNode> indexNodes = new ArrayList<indexNode>();
		Util.write(edgeString,"datasetID,startNode,endNode,distance\n");
		Util.write(nodeString,"datasetID,node,#CoverPoints,Depth,Radius,Leaf,PivotPoint\n");
		indexAlgorithm<Object> indexDSS = new indexAlgorithm<>();
		long startTime1a = System.nanoTime();//exact search: 0, 1, 0
		for(int a:dataMapPorto.keySet()) {//we can build index to accelerate
			if(a>limit+10)
				break;
			indexNode rootBall = indexDSS.buildBalltree2(dataMapPorto.get(a), dimension, capacity, null, null); 
		//	indexDSS.storeFeatures(rootBall, 1, 1, edgeString, nodeString, a);
			indexDSS.setGloabalid();
			indexDSS.storeIndex(rootBall, 1, indexString+String.valueOf(a)+".txt", 0);
			rootBall.setroot(a);//set an id to identify which dataset it belongs to
			indexMap.put(a, rootBall);
			if(a<=limit)
				indexNodes.add(rootBall);
			// get the features of index tree, 
			double features[] = datasetFeatures.getFeature(fetureString, rootBall, rootBall.getTotalCoveredPoints(), capacity, a);
			featureMap.put(a, features);
		}
		long endtimea = System.nanoTime();
		
		/**/
	//	indexAlgorithm<Object> indexDSS = new indexAlgorithm<>();
		for(int a:dataMapPorto.keySet()) {//we can build index to accelerate
			if(a>limit)
				break;
			long startTime1 = System.nanoTime();
			double distance = Hausdorff.HausdorffExact(dataMapPorto.get(1), dataMapPorto.get(a), dataMapPorto.get(a)[0].length, false);
			long endtime = System.nanoTime();	
			System.out.println((endtime-startTime1)/1000000000.0);
			System.out.println(distance);
			
			startTime1 = System.nanoTime();
			distance = Hausdorff.earlyBreaking(dataMapPorto.get(1), dataMapPorto.get(a), dimension, false);
			endtime = System.nanoTime();	
			System.out.println((endtime-startTime1)/1000000000.0);
			System.out.println(distance);
			
			startTime1 = System.nanoTime();
			double pivot_distance = Util.EuclideanDis(indexMap.get(1).getPivot(), indexMap.get(a).getPivot(), dimension);
			double ub = pivot_distance+indexMap.get(1).getRadius() + indexMap.get(a).getRadius();
			distance = Hausdorff.HausdorffWithIndexDFS(dataMapPorto.get(1), dataMapPorto.get(a), dimension, indexMap.get(1), indexMap.get(a), 0, ub, indexDSS);
			endtime = System.nanoTime();	
			System.out.println((endtime-startTime1)/1000000000.0);
			System.out.println(distance);
			
			
			startTime1 = System.nanoTime();
		//	distance = Hausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(a), dimension, indexMap.get(1), indexMap.get(a), 1);
			distance = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(a), dimension, indexMap.get(1), indexMap.get(a), 0);
			endtime = System.nanoTime();
			System.out.println((endtime-startTime1)/1000000000.0);
			System.out.println(distance);
			
			startTime1 = System.nanoTime();//exact search: 0, 1, 0
			Pair<Double, PriorityQueue<queueMain>> resultPair = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(a), dimension, indexMap.get(1), indexMap.get(a), 0, 1, 0, false, 0, false,null, null);
			endtime = System.nanoTime();
			System.out.println((endtime-startTime1)/1000000000.0+ ","+ AdvancedHausdorff.disCompTime+","+AdvancedHausdorff.EuclideanCounter);
			System.out.println(resultPair.getLeft());
			
			startTime1 = System.nanoTime();//exact search: 0, 1, 0
			distance = AdvancedHausdorff.IncrementalDistanceDirected(dataMapPorto.get(1), dataMapPorto.get(a), dimension, indexMap.get(1), indexMap.get(a), 0, 1, 0, null, null);
			endtime = System.nanoTime();
			System.out.println((endtime-startTime1)/1000000000.0+ ","+ AdvancedHausdorff.disCompTime+","+AdvancedHausdorff.EuclideanCounter);
			System.out.println(distance);//+","+resultPair1.getLeft()+","+resultPair2.getLeft()
			System.out.println();
		}
		
		AdvancedHausdorff.setBoundChoice(0);
		//conduct top-k dataset search, the baseline
		long startTime1 = System.nanoTime();
		Search.hausdorffEarylyAbandonging(dataMapPorto, limit+1, indexMap, dimension, limit, false, null, null);//no top-k early breaking
		long endtime = System.nanoTime();
		System.out.println("top-1 search costs: "+(endtime-startTime1)/1000000000.0);
		
		startTime1 = System.nanoTime();
		Search.hausdorffEarylyAbandonging(dataMapPorto, limit+1, indexMap, dimension, limit, true, null, null);
		endtime = System.nanoTime();
		System.out.println("top-1 search costs: "+(endtime-startTime1)/1000000000.0);
		
		// rank all the candidate by lower bounds (Basic and Approximaate), then one by one
		// store the queue into heap with id and hausdorff distance,
		startTime1 = System.nanoTime();
		Map<Integer, Pair<Double, PriorityQueue<queueMain>>> resultPair = new HashMap<Integer, Pair<Double,PriorityQueue<queueMain>>>();
		int datasetID = Search.HausdorffEarylyAbandongingRanking(dataMapPorto, limit+1, indexMap, dimension, limit, resultPair, null, null);
		endtime = System.nanoTime();
		System.out.println("top-1 search costs: "+(endtime-startTime1)/1000000000.0);
		System.out.println();
		
		 
		/* test join */
		startTime1 = System.nanoTime();
		Join.scanning(dataMapPorto.get(limit+1), dataMapPorto.get(datasetID), dimension);
		endtime = System.nanoTime();
		System.out.println("join with top-1 dataset with scanning baseline costs: "+(endtime-startTime1)/1000000000.0);
		
		startTime1 = System.nanoTime();
		Join.joinTableBaseline(dataMapPorto.get(limit+1), dataMapPorto.get(datasetID), 
				indexMap.get(limit+1), indexMap.get(datasetID), dimension, indexDSS);
		endtime = System.nanoTime();
		System.out.println("join with top-1 dataset with nn search baseline costs: "+(endtime-startTime1)/1000000000.0);
		
		startTime1 = System.nanoTime();
		for(int datasetID1:resultPair.keySet()) {
			Pair<Double, PriorityQueue<queueMain>> aPair = resultPair.get(datasetID1);
			Join.IncrementalJoin(dataMapPorto.get(limit+1), dataMapPorto.get(datasetID1), dimension, indexMap.get(limit+1), 
					indexMap.get(datasetID1), 0, 0, 0.01, false, 0, false, aPair.getLeft(), aPair.getRight(), null, null, "haus");
		}
		endtime = System.nanoTime();
		System.out.println("join with top-1 dataset on reused queues costs: "+(endtime-startTime1)/1000000000.0);
		
		
		System.out.println();
		/* test self-join */
		startTime1 = System.nanoTime();//exact search: 0, 1, 0
		AdvancedHausdorff.setParameter(true, false);
		Pair<Double, PriorityQueue<queueMain>> resultPaira = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(1), dimension, indexMap.get(1), indexMap.get(1), 0, 1, 0, false, 0, false, null, null);
		endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0+ ","+ AdvancedHausdorff.disCompTime+","+AdvancedHausdorff.EuclideanCounter);
		System.out.println(resultPaira.getLeft());
		System.out.println();
		
		
		/* test the restored index, it works well, while it is a little slower because of finding the map */
		startTime1 = System.nanoTime();//exact search: 0, 1, 0
		Pair<Double, PriorityQueue<queueMain>> resultPairq = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(3), dimension, indexMap.get(1), indexMap.get(3), 0, 1, 0, false, 0, false,null, null);
		endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0+ ","+ AdvancedHausdorff.disCompTime+","+AdvancedHausdorff.EuclideanCounter);
		System.out.println(resultPairq.getLeft());
		
		startTime1 = System.nanoTime();//exact search: 0, 1, 0
		Map<Integer, Map<Integer, indexNode>> datasetIndex  = indexDSS.restoreIndex(indexString, dimension, dataMapPorto);
		endtime = System.nanoTime();
		System.out.println("building index costs: "+(endtimea-startTime1a)/1000000000.0);
		System.out.println("reconstruction costs: "+(endtime-startTime1)/1000000000.0);
		
		startTime1 = System.nanoTime();//exact search: 0, 1, 0
		Pair<Double, PriorityQueue<queueMain>> resultPair11 = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(3), dimension, datasetIndex.get(1).get(1), datasetIndex.get(3).get(1), 0, 1, 0, false, 0, false, datasetIndex.get(1), datasetIndex.get(3));
		endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0+ ","+ AdvancedHausdorff.disCompTime+","+AdvancedHausdorff.EuclideanCounter);
		System.out.println(resultPair11.getLeft());
		
		
//		indexNodes = new ArrayList<indexNode>();
//		for(int dataid: datasetIndex.keySet()) {
//			indexNodes.add(datasetIndex.get(dataid).get(1));
//		}
		
		/*
		 * test the index of whole data lake
		 */
		System.out.println("indexing datasets");
		indexNode datasetRoot = null;
		Map<Integer, indexNode> datalakeIndex = null;
		indexDSS.setGloabalid();
		File tempFile = new File(indexString+"datalake.txt");
		if(!tempFile.exists()) {
			datasetRoot = indexDSS.indexDatasetKD(indexNodes, dimension, capacity);
			indexDSS.storeDatalakeIndex(datasetRoot, 1, indexString+"datalake.txt", 0);//store the 
			//	int rootNum = indexDSS.getRootCount(datasetRoot);
			//	System.out.println("number of datasets: "+ rootNum);
		}else
			datalakeIndex = indexDSS.restoreDatalakeIndex(indexString+"datalake.txt", dimension);

//		indexDSS.updateDatasetNodeId(datasetRoot, 1);
		 
	//	datalakeIndex = null;
		
		// update index
	//	indexNode rootBall = indexDSS.buildBalltree2(dataMapPorto.get(limit+1), dimension, capacity); 
	//	rootBall.setroot(limit+1);
	//	datasetRoot = indexDSS.insertNewDataset(rootBall, datasetRoot, capacity, dimension);
	//	rootNum = indexDSS.getRootCount(datasetRoot);
	//	System.out.println("number of datasets: "+ rootNum);
	//	dataMapPorto = null;
		int k = 10;
		startTime1 = System.nanoTime();
		// test top-k search, which also works for restored index.
		int queryid = 10;
		HashMap<Integer, Double> result = Search.pruneByIndex(dataMapPorto, datasetRoot, indexMap.get(queryid), queryid, 
				dimension, indexMap, datasetIndex, datasetIndex.get(queryid), datalakeIndex, argoDataMap, k, indexString);
		endtime = System.nanoTime();
		System.out.println("top-1 search costs: "+(endtime-startTime1)/1000000000.0);
		System.out.println(result+", "+argoDataMap.get(result.entrySet().iterator().next().getKey()));
		
		//run framework to generate the index first, then go to the effectiveness study to test offline index.
	}
}
