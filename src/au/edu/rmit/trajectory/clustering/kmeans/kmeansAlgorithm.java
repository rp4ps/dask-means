package au.edu.rmit.trajectory.clustering.kmeans;

import java.io.BufferedReader;
import java.io.File;
import java.text.DecimalFormat;
import java.io.FileDescriptor;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.rmit.trajectory.clustering.kpaths.KPathsOptimization;
import java.util.PriorityQueue;
import au.edu.rmit.trajectory.clustering.kpaths.Util;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;


/*
 * ball-tree, M-tree, and Hierarchical k-means tree can used be extended to answer k-means
 */
@SuppressWarnings("restriction")
public class kmeansAlgorithm<T> extends KPathsOptimization<T>{
	protected ArrayList<cluster> CENTERSEuc; // it stores the k clusters
	double [][]dataMatrix;//|D|*dimension, we use array
	
	double [][]centroidsData; //k*dimension, store the centorid, this can be changed as cache for optimization,
	
	double [][]allBounds;//|D|*(2+k) or |D|*(2+b) or |D|*2, store bounds for every point.
	
	double []dataNorm;//|D|, l2 norm
	double [][]dataDivvector;//|D|*b, blockvector
	double []centroidNorm;//k, l2 norm
	double [][]centoridDivvector; //k*b
	int userID[];// for fair clustering
	Map<Integer, Integer> userNumber;// for fair clustering
	boolean readOwner = false;
	
	double []distanceToFather; //|D|, for index-based methods
	protected int [][]allCentroids = null; //#test*k, storage a set of k center id for test various methods
	double []heapBound; // k, 
	double previoustime = 0;
	int numberofComparison = 18;//all the comparison
	
	String datafilename;
	int dimension = 0;//the dimension of the Euclidean dataset
	int dimension_start = 0;// the dimension start in a file
	int dimension_end = 0;// the dimension end in file
	String split = null;	/*use to split the colume of dataset*/
	
	indexAlgorithm<Object> indexkmeans;
	indexNode root;// the root node
	indexNode rootCentroids;// the root node of centroid index
	int capacity=30;
	
	int prunenode=0;
	long numComputeEuc = 0;
	long boundCompare = 0;
	long dataReach = 0;
	long nodeReach = 0;
	long boundUpdate = 0;
	double timedistanceiter = 0;
	
	double time[]= new double[numberofComparison];// the overall running time
	double assigntime[] = new double[numberofComparison];
	double refinetime[] = new double[numberofComparison];
	long computations[] = new long[numberofComparison];// the number of computations
	
	long boundCompares[] = new long[numberofComparison];//used for bound accessing, also for bound computation.
	long dataAccess[] = new long[numberofComparison];
	long boundUpdates[] = new long[numberofComparison];
	double memoryUsage[] = new double[numberofComparison];
	double distanceTime[] = new double[numberofComparison];
	int counter = 0; // algorithm indicator

	boolean nonkmeansTree = true;
	
	boolean lloyd = false;
	boolean usingIndex = false;
	boolean indexPAMI=false;// travesal from root or not?
	boolean Phillips = false;//Acceleration of k-means and related clustering algorithms, a simple case of elkan
	boolean elkan = false;
	boolean Hamerly = true;//Making k-means even faster, using the drift as bound, SDM'10
	boolean Drake12 = false;//this is for storing b bounds
	boolean annulus = false;
	boolean wsdm14 = false; //
	boolean heap = false;
	boolean Yinyang = false;// ICML'15: Yinyang K-means: A drop-in replacement of the classic K-means with consistent speedup
	boolean Exponion = false; // newling 16
	boolean blockVector= false;// icml 16
	boolean reGroup = false; // regrouping idea of ICAISC 2017
	boolean sdm16 = false;
	boolean pami20 = false;
	
	boolean pckmeans = false;
	boolean pckmeansbound = false;
	boolean pckmeanspointbound = true; // using only bound for points for lightweight k-means
	boolean pckmeanspointboundRecursive = false; //using recursive or queue to traverse the tree
	boolean pckmeansUsinginterbound = true;
	
	boolean runIndexMethodOnly = false;
	boolean runBalltreeOnly = false;
	

	int numberBlockVectors = 2;// this is the default number from the paper ICML 16
	int binSortcenter = 0;
	int iterationTimes;
	int MAXITE = 10;// can be set as bigger number
	int maxIteration = MAXITE;
	double maxDrift = 0;
	
	double []maxdis;//k, for sdm'16 to store the maximum radius for each cluster.
	double [][]tighterdrift;//k*k
	
	/* lightweight k-means using dual-tree method, not worked well, we then only used the assigned*/
	double maxdrift[];// store the maximum drift excluding the one access
	short assigned[]; // store the assigned iterations
	short newassigned[]; // store the assigned iterations
	double []quantilizedLowerBound=null;// create new only when leaf node is not pruned, or not used, optional
	double []quantilizedUpperBound=null;// used for storing bounds on nn
	short []pointcounterPruned=null;
	double unitQuantization;
	
	/* other configuration */
	boolean usingUpperBound = false; // compute the distance to current cluster directly, which can prune more than using upper bound.
	boolean kmeansplusplus = false;// true; //use the kmean++ initialization
	boolean shownSum = false;//print the sse for checking whether the algorithm is right
	int kThresholdPick = 1000; // to judge whether to compute the real distance with assigned cluster, tronto dataset should use have
	
	boolean pami20_flags[] = new boolean[k]; //indicate whether each cluster is stable
	ArrayList<Boolean> PAMI20ClusterNeedScanningArrayList = new ArrayList<Boolean>();
	double []radiusPAMI20;
	double indexingTime = 0;
	double interBoundTime = 0;
	boolean early_terminating;
	
	indexAlgorithm<indexNode> algorithm = new indexAlgorithm<indexNode>();
	/*
	 * set sign to test,
	 */
	void setSign(int sign[]) {
		lloyd = sign[0]==1?true:false;
		usingIndex = sign[1]==1?true:false;
		elkan = sign[2]==1?true:false;
		Hamerly = sign[3]==1?true:false;	
		Drake12 = sign[4]==1?true:false;
		annulus = sign[5]==1?true:false;
		wsdm14 = sign[6]==1?true:false;
		heap = sign[7]==1?true:false;
		Yinyang = sign[8]==1?true:false;
		Exponion = sign[9]==1?true:false;
		blockVector = sign[10]==1?true:false;
		reGroup = sign[11]==1?true:false;
		sdm16 = sign[12]==1?true:false;
		if(lloyd)
			assBoundSign = false;
		else {
			assBoundSign = true;
		}
		if(Yinyang)
			elkan = true;
		if(Hamerly && elkan)//they cannot be true in the same time
			Hamerly = false;
		if(dimension%2==1 || k<5) //we found it has to be an even for the blockvector algorithm, k cannot be 
			blockVector = false;
	}
	
	public kmeansAlgorithm(String []datapath) {
		super(datapath);
		datafilename = datapath[4];
		dimension_start = Integer.valueOf(datapath[5]);// the dimensions
		dimension_end = Integer.valueOf(datapath[6]);
		maxIteration = MAXITE;
		dimension = dimension_end - dimension_start+1;
		if(datapath.length>7)
			split = datapath[7];
	}
	
	public void setScale(int scale) {
		trajectoryNumber = scale;
	}
	
	public void setCapacity(int scale) {
		capacity = scale;
	}
	
	public void setDimension(int dim) {
		dimension = dim;
		dimension_end =dimension_start+dimension-1;
	}
	
	public int getDimension() {
		return dimension;
	}
	
	public void loadDataEuc(String path, int number) throws FileNotFoundException, IOException {
		dataMatrix = new double[trajectoryNumber][];//store the data in this matrix.
		int pointid = 1;
		if(path.contains("fair")) {
			readOwner = true;
			userID = new int[trajectoryNumber];
		}
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(path)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = null;				
				if(split!=null) {
					if (!split.equals("a")) {
						abc = strr.split(split);
					} else if (split.equals("b")) {
						abc = strr.split(" ");
					} else {
						if(strr.contains(","))
							abc = strr.split(",");// when both connected 
						else if(strr.contains(";"))
							abc = strr.split(";");
						else if(strr.contains(" "))
							abc = strr.split(" ");
					}
				}else{
					if(strr.contains(","))
						abc = strr.split(",");// when both connected 
					else if(strr.contains(";"))
						abc = strr.split(";");
					else if(strr.contains(" "))
						abc = strr.split(" ");
				}
				if(abc.length-1 < dimension_end)// the dimension is not right
					continue;
				dataMatrix[pointid-1] = new double[dimension];
				if(readOwner) {// get the user id of point, for fair clustering
					userID[pointid-1] = Integer.valueOf(abc[0]);// read 
				}
				double []point = new double[dimension];
				boolean nan=false;
				for(int i=dimension_start; i<=dimension_end; i++) {					
					if(!abc[i].equals("?")) {
						dataMatrix[pointid-1][i-dimension_start] = Double.valueOf(abc[i]);
						point[i-dimension_start] = Double.valueOf(abc[i]);
					}
					else {
						dataMatrix[pointid-1][i-dimension_start] = 0;
						point[i-dimension_start] = 0;
					}
				}				
				if(nan==true)
					continue;				
				pointid++;
				if(pointid>number)
					break;
			}
			in.close();
		}		
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	//	storeNYCLibSVM();
	//
		if(readOwner) {//for fair k-means
			getUsersCheckInPointNumber("./fairKmeans/userNumber"+trajectoryNumber+".txt");
			String filenString = "./dataset/processedDataset"+trajectoryNumber+".txt";
			File tempFile = new File(filenString);
			if (!tempFile.exists())
				fairkmeans.FairProcessData(dataMatrix, userID, userNumber, dimension, filenString);
		}
	}
	
	
	void storeNYCLibSVM() {
		for(int i=0; i<trajectoryNumber; i++) {
			if(dataMatrix[i][0]>-80 && dataMatrix[i][0]<-70 && dataMatrix[i][1]<50 &&dataMatrix[i][1]>30) {
				String contentString = i+" "+"1:"+Double.toString(dataMatrix[i][0])+" 2:"+Double.toString(dataMatrix[i][1])+"\n";
				Util.write("/Users/sw160/Desktop/NYC_TAXI_Location_2013_Oct_libsvm_million.txt", contentString);
			}
			if(i%10000==0) {
				System.out.println(i);
			}
		}
	}
	/*
	 * read the posting frequency of each point
	 */
	void getUsersCheckInPointNumber(String numberfile) throws FileNotFoundException, IOException {
		File tempFile = new File(numberfile);
		userNumber = new HashMap<Integer, Integer>();
		if (!tempFile.exists()) {
			for (int point : userID) {
				int number = 0;
				if (userNumber.containsKey(point))
					number = userNumber.get(point);
				number++;
				userNumber.put(point, number);
			}
			for (int point : userNumber.keySet()) {
				Util.write(numberfile, point + "," + userNumber.get(point) + "\n");
			}
		}else {
			try (BufferedReader br = new BufferedReader(new FileReader(numberfile))) {
				String strLine;
				while ((strLine = br.readLine()) != null) {
					String[] splitString = strLine.split(",");
					userNumber.put(Integer.valueOf(splitString[0]), Integer.valueOf(splitString[1]));
				}
			}
		}
	}
	
	/*
	 * access the data using point id
	 */
	double []accessPointById(int id){
		dataReach++;
		return dataMatrix[id-1];
	}
	
	/*
	 * compute the distance between any two centers for bound computation
	 */
	public void computeInterCentoridEuc(int k, ArrayList<cluster> Center, double [][]clustData) {
		for(int i=0; i<k; i++) {
			innerCentoridDis[i] = new double[k];
			double []a = clustData[i];		
			double min = Double.MAX_VALUE;
			for(int j=0; j<k; j++) {				
				if(i!=j) {
					double []b = clustData[j];
					long startTime1 = System.nanoTime();
					double distance = Util.EuclideanDis(a, b, a.length);
					long endtime = System.nanoTime();	
					timedistanceiter += (endtime-startTime1)/1000000000.0;
					
					numComputeEuc++;
					innerCentoridDis[i][j] = distance;
					if(distance<min) {
						min = distance;
					}
				}
			}
			interMinimumCentoridDis[i] = min;
		}
		for (int i = 0; i < k; i++) {//the distance in each group
			innerCentoridDisGroup[i] = new double[group.size()];
			for (int groupid : group.keySet()) {
				ArrayList<Integer> arrayList = group.get(groupid);
				double min = Double.MAX_VALUE;
				for (int centerid : arrayList) {
					if (innerCentoridDis[i][centerid] < min) {
						min = innerCentoridDis[i][centerid];
					}
				}
				innerCentoridDisGroup[i][groupid] = min;
			}
		}
	}
	
	// we optimize the above algorithm using knn and join, 
	public void computeInterCentoridEuckNN(int k, ArrayList<cluster> Center, double [][]clustData) {
	//	AdvancedHausdorff.IncrementalDistance(point1xys, point2xys, k, X, Y, splitOption, fastMode, error, reverse, directDis, topkEarlyBreaking, nodelist, nodelist1);
		indexAlgorithm<indexNode> algorithm = new indexAlgorithm<indexNode>();
		double[] minDistnearestID = new double[4];
		for(int i=0; i<k; i++) {
			for(int j=0; j<4;j++)
				minDistnearestID[j] = Double.MAX_VALUE;
			if(iterationTimes>0)
				minDistnearestID[0] = interMinimumCentoridDis[i] + group_drift[i] + maxDrift;// the bound for centroid index pruning
			double []a = clustData[i];
			algorithm.TwoNearestNeighborSearchBall(a, rootCentroids, dimension, clustData, minDistnearestID);
			interMinimumCentoridDis[i] = minDistnearestID[0];
		}
	}
	
	
	/*
	 * update the lower bound of trajectory toward group i
	 */
	protected void updateSingleLowerBound(int traid, 
			int group_i, double newbound, int groupNumber) {
		if(allBounds[traid-1] == null) {
			allBounds[traid-1] = new double[groupNumber+2];
			allBounds[traid-1][0] = Double.MAX_VALUE;
			for(int i=1; i<groupNumber+2; i++) {
				allBounds[traid-1][i] = 0;//initialize the bound as zero
			}
			if(elkan)
				allBounds[traid-1][group_i+2] = newbound;
		}else {
			if(elkan && allBounds[traid-1][group_i+2] > newbound)
				allBounds[traid-1][group_i+2] = newbound;
		}
	}
	
	/*
	 * initialize the centroid
	 */
	public void initializeClustersRandom(int k){
		CENTERSEuc = new ArrayList<>();
		int unit= trajectoryNumber/k;
		System.out.println(unit+ " "+ k+" "+trajectoryNumber);
		if(centroids==null) {//if there is no given centroids
			Random rand = new Random();
			for(int t=0; t<k; t++) {
				cluster cl = null;
				int  n = rand.nextInt(trajectoryNumber)+1;
				n = t+1;// using same centroids every time.
				double[] cluster = accessPointById(n);
				if(usingIndex) {// if there is an index
					if(t==k-1) {
						if(pckmeanspointboundRecursive)
							cl = new cluster(cluster, dimension);
						else
							cl = new cluster(cluster, root, dimension);//assign the root node to the last cluster
					}
					else {
						cl = new cluster(cluster, dimension);
					}
				}else {					
					int end = (t+1)*unit;
					if(t == k-1)
						end = trajectoryNumber;
					cl = new cluster(cluster, n, t*unit, end, dimension, dataMatrix);//no index
				}
				CENTERSEuc.add(cl);
			}
		}else {
			for(int t=0; t<k; t++) {
				int  n = centroids[t];
				double[] cluster = accessPointById(n);
				cluster cl = null;
				if(usingIndex) {// if there is an index
					if(t==0) {
						if(pckmeanspointboundRecursive)
							cl = new cluster(cluster, dimension);
						else 
							cl = new cluster(cluster, root, dimension);//assign the root node to the first cluster
					}
					else {
						cl = new cluster(cluster, dimension);
					}
				}else {					
					int end =(t+1)*unit;
					if(t == k-1)
						end = trajectoryNumber;
					cl = new cluster(cluster, n, t*unit, end, dimension, dataMatrix);//no index
				}
				CENTERSEuc.add(cl);
			}
		}
		System.out.println("Centroid is initialized");
	}
	
	/*
	 * initialize the k clusters by randomly choosing from existing points
	 * this is for building the indexing
	 */
	public void initializeClustersRandomForIndex(int k, Set<Integer> data) {
		CENTERSEuc = new ArrayList<>();
		Random rand = new Random();
		ArrayList<Integer> arrayList = new ArrayList<>(data);
		int unit = data.size()/k;
		ArrayList<Integer> tempor = new ArrayList<>();
		int num=0;
		System.out.println("data is");
		for(int t=0; t<data.size(); t++) {
			int idx = arrayList.get(t);
			tempor.add(idx);
			if((t+1)%unit == 0 || t==data.size()-1) {
				if(CENTERSEuc.size() == k-1 && t!=data.size()-1) {
					continue;
				}
				int n = rand.nextInt(unit);
				if( t == data.size()-1)
					n = 0;//the last group, we assign the value directly
				idx = tempor.get(n);
				num += tempor.size();
				double[] clusterdata = accessPointById(idx);// we may loose some data here.
				System.out.println("data is"+idx);
				cluster cl = new cluster(clusterdata, tempor, dimension, dataMatrix);//we need to assign the data gradually
				CENTERSEuc.add(cl);
				tempor = new ArrayList<>();
			}
		}
		System.out.println("total initialized points: "+num);
	}
	
	/**
	 * using the algorithms of k-means++ published in soda 2007,
	 * 
	 * it will cost much time
	 */
	private int[] kmeanplusplus(int dimen, int k, int number) {
		double[][] centroidsdata = new double[k][dimen];
		double[] distToClosestCentroid = new double[number];
		double[] weightedDistribution = new double[number]; // cumulative sum of squared distances
		int centoridid[] = new int[k];
		Random gen = new Random();
		int choose = 0;
		for (int c = 0; c < k; c++) {
			if (c == 0)
				choose = gen.nextInt(number);
			else {
				for (int p = 0; p < number; p++) {
					double tempDistance = Util.EuclideanDis(dataMatrix[p], centroidsdata[c - 1], dimen);
					if (c == 1)
						distToClosestCentroid[p] = tempDistance;
					else {
						if (tempDistance < distToClosestCentroid[p])
							distToClosestCentroid[p] = tempDistance;
					}
					if (p == 0)
						weightedDistribution[0] = distToClosestCentroid[0];
					else
						weightedDistribution[p] = weightedDistribution[p - 1] + distToClosestCentroid[p];
				}
				double rand = gen.nextDouble();
				for (int j = number - 1; j > 0; j--) {
					if (rand > weightedDistribution[j - 1] / weightedDistribution[number - 1]) {
						choose = j; // one bigger than the one above
						break;
					} else
						choose = 0;
				}
			}
			for (int i = 0; i < dimen; i++)
				centroidsdata[c][i] = dataMatrix[choose][i];
			centoridid[c] = choose+1;
		}
		return centoridid;
	}
	
	/*
	 * we select the centroids from the leaf nodes, 
	 */
	
	/*
	 * use the k-means to build the index HKM-tree, when a group has a radius less than the 
	 * threshold use the capacity as k to run kmeans
	 */
	public indexNode runIndexbuildQueuePoint(double radius, int capacity, int fanout) throws IOException {
		root = new indexNode(dimension);// this will store the 
		k = fanout;
		GroupedTrajectory = new HashSet<>();
		System.out.println("Building Hiarachical k-means tree...");
		String LOG_DIR = "./index/index.log";
		PrintStream fileOut = new PrintStream(LOG_DIR);
		System.setOut(fileOut);	
		String groupFilename = datafile+"_"+Integer.toString(trajectoryNumber)+"_"+Double.toString(radius)+"_"+Integer.toString(capacity)+"_index";
		String pivotname = groupFilename+".all";
		pivotGroup = new HashMap<>();
		Queue<Pair<Set<Integer>, indexNode>> queue = new LinkedList<>();//queue with
		Set<Integer> KeySet = new HashSet<>();
		for(int i=1; i<=trajectoryNumber; i++)
			KeySet.add(i);
		Pair<Set<Integer>, indexNode> aPair = new ImmutablePair<>(KeySet, root);// build the index using a queue
		queue.add(aPair);
		usingIndex = false;
		assBoundSign = true;// we use the bound here to accelerate the indexing building.
		indexPAMI = false; // whether traverse again
		int firstIteration = 0;
		long startTime1 = System.nanoTime();
		while(!queue.isEmpty()) {
			aPair = queue.poll();
			Set<Integer> candidates = aPair.getKey();
			indexNode fatherNode = aPair.getValue();// the nodes 
			CENTERSEuc = new ArrayList<cluster>();	
			interMinimumCentoridDis = new double[k];
			innerCentoridDis = new double[k][];								
			System.out.println("#data points: "+candidates.size());// the size of the dataset
			initializeClustersRandomForIndex(k, candidates);// initialize the center	
			runkmeans(fanout, candidates);// run k-means to divide into k groups
			String content = "";
			for(int i = 0; i<fanout; i++) {
				cluster node = CENTERSEuc.get(i);
				int nodeCapacity = node.getcoveredPoints().size();
			//	content += Integer.toString(node.getTrajectoryID()) + ","; // we should write the name
			}
			int num = 0;
			for(int i = 0; i<fanout; i++) {
				indexNode childNode = new indexNode(dimension);
				cluster node = CENTERSEuc.get(i);
				candidates = node.getcoveredPoints();
				int nodeCapacity = candidates.size();
				num += nodeCapacity;
				if(nodeCapacity==0)// no trajectory
					continue;
			//	System.out.print(nodeCapacity+";");
				double[] nodeSum = new double[dimension];
				double nodeRadius = node.getRadius(candidates, dataMatrix, nodeSum);//compute the distance
			//	System.out.println(nodeSum[0]);
			//	double[] nodeSum = node.getSumTotal();
				childNode.setRadius(nodeRadius);
				childNode.setSum(nodeSum);
				double[] newpivot = new double[dimension];
				for(int j=0; j<dimension; j++)
					newpivot[j] = nodeSum[j]/nodeCapacity;
				childNode.setPivot(node.getcentroid());
				childNode.setTotalCoveredPoints(nodeCapacity);
				content = nodeRadius + ":";
				for (int idx : candidates) {
					content += Integer.toString(idx) + ",";
				}
				if(nodeCapacity <= capacity || nodeRadius <= radius) {// stop splitting and form the leaf node
				/*	if(firstIteration == 0 && !content.equals("0:0,"))
						Util.write(pivotname, content+"\n");//write all the contents into the pivot table
					if (nodeCapacity >= capacity/2 && nodeRadius <= radius) {
						Util.write(groupFilename, content + "\n");// write the group into file
					}else if(nodeCapacity>0) {					
						GroupedTrajectory.addAll(candidates);// add to another file
					}*/
					childNode.addPoint(candidates);	// this is the leaf node				
				}else {
					aPair = new ImmutablePair<Set<Integer>, indexNode>(candidates, childNode);
					queue.add(aPair);// conduct the iteration again
				}
				fatherNode.addNodes(childNode);//add the nodes.
			}
		//	System.out.println("after clustering: "+num);// the size of the dataset
		}
		long endtime = System.nanoTime();	
		System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
	//	System.out.println("the height is: "+getHeight(root)+", the #points is: "+getcount(root));
		System.out.println("HKM Indexing time");
		System.out.print((endtime-startTime1)/1000000000.0+ "\n");
		return root;//return the root node of the tree.
	}
	
	/*
	 * get the minimum bound as the global bound
	 */
	public double getMinimumLowerbound(double [] bounds, int groupNumber, int grouplocate) {
		double lowerboud = Double.MAX_VALUE;		
		for(int group_j=0; group_j<groupNumber; group_j++) {//get the minimum lower bound of all group
			if(group_j == grouplocate)
				continue;
			double lowerboud_temp = bounds[group_j+2] - group_drift[group_j];
			if(lowerboud > lowerboud_temp)//choose the minimum one
				lowerboud = lowerboud_temp;
		}
		return lowerboud;
	}
	
	/*
	 * store into array, or into cache.
	 */
	public void storeCentroid() {
		for (int j = 0; j < k; j++) {// combine the inverted index for pruning
			centroidsData[j] = new double[dimension];
			double []clustra = CENTERSEuc.get(j).getcentroid();
			for(int d=0; d<dimension; d++) {
				centroidsData[j][d] = clustra[d];
			}
		}
	}
	
	/*
	 * regroup idea, to run regroup again
	 */
	void regroup(int groupNumber) throws IOException {
		int roundRegroup = 5;//regroup (ICAISC'17), every 5 iterations, we regroup
		if(assBoundSign && reGroup && (iterationTimes+1)%roundRegroup == 0) {
			for(int i=0; i<trajectoryNumber; i++) {// reset all lower bounds to be zero, only remain two global bounds
				if(allBounds[i] !=null)
					for(int j=0; j<groupNumber; j++)
						allBounds[i][j+2] = 0;
			}
			if(root!=null)
				setBoundsforRegrouping(root, groupNumber); //reset the bound in the node for next test
			groupInitialClusters(groupNumber, k);
		}
	}
	
	/*
	 * use 1nn and 2nn to assign when centroid is big, we can argumented with a bound.
	 */
	double[] assignPCKmenas(double pivot[], indexNode centroidRoot, double centerMatrix[][], boolean isNode) {
		indexAlgorithm<indexNode> algorithm = new indexAlgorithm<indexNode>();
		double[] minDistnearestID;
		if(isNode) {
			minDistnearestID = new double[4];
			for(int i=0; i<4; i++)
				minDistnearestID[i] = Double.MAX_VALUE; // initilaized it as the bound
			algorithm.TwoNearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
		}else {
			minDistnearestID = new double[2];
			for(int i=0; i<2; i++)
				minDistnearestID[i] = Double.MAX_VALUE;
			algorithm.NearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
		}
		return minDistnearestID;
	}
	
	/*
	 * use 1nn and 2nn to assign when centroid is big, we augmented with multiple bounds to accelerate nn over the centroid index
	 */
	double[] assignPCKmenasBound(double pivot[], indexNode centroidRoot, double centerMatrix[][], boolean isNode, indexNode node, int idx, double radius) {
		double[] minDistnearestID;
		if(isNode) {
			minDistnearestID = new double[4];
			int assignedCluster = node.getAssignedCluster();//get the assigned cluster
			if(assignedCluster != 0) {//compute the distance to assigned cluster
				minDistnearestID[3] = (double)assignedCluster;
				if(iterationTimes<=0 || k <= kThresholdPick || assignedCluster<=0) {
					minDistnearestID[2] = node.getUpperBoundPick();// + group_drift[assignedCluster-1];// the drift bound
				}else
					minDistnearestID[2] = Util.EuclideanDis(pivot, centerMatrix[assignedCluster-1], dimension);// this will cost more time
				minDistnearestID[0] = minDistnearestID[2] + interMinimumCentoridDis[assignedCluster-1];// for the bound on second nearest neighbor
				if(minDistnearestID[2] + radius < interMinimumCentoridDis[assignedCluster-1]/2) {// pruned directly
					return minDistnearestID;
 				}
			}else {
				for(int i=0; i<4; i++)
					minDistnearestID[i] = Double.MAX_VALUE; // initilaized it as the bound
			}
			double temp = node.getUpperBoundPick();// + maxDrift;// the bound based on previous iterations
			if(temp < minDistnearestID[0])
				minDistnearestID[0] = temp;// get the bound from father node
			algorithm.TwoNearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
			if(minDistnearestID[0] - minDistnearestID[2] >= 2 * radius) {
				if((int)minDistnearestID[3]!=assignedCluster) {
					assigned = node.setAssignPICK((short)minDistnearestID[3], assigned);// all assigned to the nearest cluster
				}
			}else{
				node.unAssignPICK();//put as unassigned
			}
		}else {
			minDistnearestID = new double[2];
			if(assigned!=null && assigned[idx-1] != 0) {//check whether can assign directly, and get a tighter bound if not
				int assignedCluster = assigned[idx-1];
				minDistnearestID[1] = (double)assignedCluster;
				if(iterationTimes<=0 || k <= kThresholdPick || assignedCluster<=0)
					minDistnearestID[0] = quantilizedUpperBound[idx-1];// + group_drift[assignedCluster-1]; //the stored drift bound
				else 
					minDistnearestID[0] = Util.EuclideanDis(pivot, centerMatrix[assignedCluster-1], dimension);
				if(minDistnearestID[0] < interMinimumCentoridDis[assignedCluster-1]/2) {
					return minDistnearestID;
				}
				if(quantilizedUpperBound[idx-1] < minDistnearestID[0])
					minDistnearestID[0] = quantilizedUpperBound[idx-1];
			}else {
				for(int i=0; i<2; i++)
					minDistnearestID[i] = Double.MAX_VALUE;
				if(quantilizedUpperBound!=null && iterationTimes==0 && quantilizedUpperBound[idx-1] < minDistnearestID[0])
					minDistnearestID[0] = quantilizedUpperBound[idx-1];
			}
			algorithm.NearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
			if(iterationTimes==0 && minDistnearestID[1]>=k && assigned!=null) {
				minDistnearestID[0] = Double.MAX_VALUE;
				algorithm.NearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
			}
			assigned[idx-1] = (short)minDistnearestID[1];//assigned cluster
		}
		return minDistnearestID;
	}
	
	/*
	 * use 1nn and 2nn to assign when centroid is big, we use a lightweight recursive approach to accelerate
	 */
	double[] assignPCKmenasBoundRecursive(double pivot[], indexNode centroidRoot, double centerMatrix[][], boolean isNode, 
			indexNode node, int idx, double radius, double knnbound, Map<Integer, ArrayList<Integer>> idxNeedsIn, 
			Map<Integer, ArrayList<Integer>> idxNeedsOut, Map<Integer, ArrayList<indexNode>> nodeNeedsIn) {
		double[] minDistnearestID;
		if(isNode) {
			minDistnearestID = new double[4];
			int assignedCluster = node.getAssignedCluster();//get the assigned cluster
			if(assignedCluster != 0) {//compute the distance to assigned cluster
				minDistnearestID[3] = (double)assignedCluster;
				if(iterationTimes<=0 || k <= kThresholdPick || assignedCluster<=0) {// here we use the drift when k is relatively small
					minDistnearestID[2] = knnbound; // + group_drift[assignedCluster-1];// the drift bound
				}else {
					numComputeEuc++;
					minDistnearestID[2] = Util.EuclideanDis(pivot, centerMatrix[assignedCluster-1], dimension);// this will cost more time
				}
				if(pckmeansUsinginterbound && minDistnearestID[2] + radius < interMinimumCentoridDis[assignedCluster-1]/2) {// pruned directly
					prunenode++;
					CENTERSEuc.get((int)minDistnearestID[3]-1).addNode(node);
					return minDistnearestID;
 				}
				minDistnearestID[0] = minDistnearestID[2] + interMinimumCentoridDis[assignedCluster-1];// for the bound on second nearest neighbor
			}else {
				for(int i=0; i<4; i++)
					minDistnearestID[i] = Double.MAX_VALUE; // initilaized it as the bound
			}
			if(knnbound < minDistnearestID[0])
				minDistnearestID[0] = knnbound;// get the bound from father node
			algorithm.TwoNearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
			if(minDistnearestID[0] - minDistnearestID[2] >= 2 * radius) {// can be assigned
				if((int)minDistnearestID[3] != assignedCluster) {
					assigned = node.setAssignPICK((short)minDistnearestID[3], assigned);// all assigned to the nearest cluster
				}
				prunenode++;
				CENTERSEuc.get((int)minDistnearestID[3]-1).addNode(node);
			}else{
				node.unAssignPICK();//put as unassigned
				lightweight(minDistnearestID, radius, pivot, centroidRoot, centerMatrix, isNode, node, idx, idxNeedsIn, idxNeedsOut, nodeNeedsIn);// keep exploring child nodes
			}
		}else {
			minDistnearestID = new double[2];
			int assignedCluster = 0;
			if(assigned!=null && assigned[idx-1] != 0) {//check whether can assign directly, and get a tighter bound if not
				assignedCluster = assigned[idx-1];
				minDistnearestID[1] = (double)assignedCluster;
				if(iterationTimes<=0 || k <= kThresholdPick || assignedCluster<=0)// here we use the drift when k is relatively small
					minDistnearestID[0] = knnbound;// + group_drift[assignedCluster-1]; //the stored drift bound
				else {
					numComputeEuc++;
					minDistnearestID[0] = Util.EuclideanDis(pivot, centerMatrix[assignedCluster-1], dimension);
				}
				if(pckmeansUsinginterbound && minDistnearestID[0] < interMinimumCentoridDis[assignedCluster-1]/2) {
					CENTERSEuc.get((int)minDistnearestID[1]-1).addPointToCluster(idx, pivot);
					return minDistnearestID;
				}
				if(knnbound < minDistnearestID[0])
					minDistnearestID[0] = knnbound;
			}else {
				for(int i=0; i<2; i++)
					minDistnearestID[i] = Double.MAX_VALUE;
				if( iterationTimes==0 && knnbound < minDistnearestID[0])
					minDistnearestID[0] = knnbound;
			}
			algorithm.NearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
			if(iterationTimes==0 && minDistnearestID[1]>=k && assigned!=null) {
				minDistnearestID[0] = Double.MAX_VALUE;
				algorithm.NearestNeighborSearchBall(pivot, centroidRoot, dimension, centerMatrix, minDistnearestID);
			}
			assigned[idx-1] = (short)minDistnearestID[1];//assigned cluster
			CENTERSEuc.get((int)minDistnearestID[1]-1).addPointToCluster(idx, pivot);
		}
		return minDistnearestID;
	}
		
	//using recursive way to traversal
	void lightweight(double[] minDistnearestID, double radius, double pivot[], indexNode centroidRoot, double centerMatrix[][], boolean isNode, indexNode node, int idx, Map<Integer, ArrayList<Integer>> idxNeedsIn, 
			Map<Integer, ArrayList<Integer>> idxNeedsOut, Map<Integer, ArrayList<indexNode>> nodeNeedsIn) {
		if(node.isLeaf()) {
			double knnbound = minDistnearestID[2] + radius;
			for(int pointid :node.getpointIdList()) {
				minDistnearestID = assignPCKmenasBoundRecursive(dataMatrix[pointid-1], rootCentroids, centroidsData, false, null, pointid, 0, knnbound, idxNeedsIn, idxNeedsOut, nodeNeedsIn);
			}
		}else {
			double knnbound= minDistnearestID[0] + radius;
			for(indexNode childnode :node.getNodelist()) {
				minDistnearestID = assignPCKmenasBoundRecursive(childnode.getPivot(), rootCentroids, centroidsData, true, childnode, 0, childnode.getRadius(), knnbound, idxNeedsIn, idxNeedsOut, nodeNeedsIn);
			}
		}
	}

	// assignment, all the optimizations are in one functions, Elkan, Hamerly, Yinyang, Newlying's bound will be used here
	// if a node cannot be safely assigned, we need to remove it from current cluster and add the children into the queue
	public void assignmentBounds(int k, int groupNumber) throws IOException {
		numeMovedTrajectories= 0;
		long start = System.nanoTime();
		Map<Integer, ArrayList<Integer>> idxNeedsIn = new HashMap<>();//it stores all the idxs of trajectories that move in
		Map<Integer, ArrayList<Integer>> idxNeedsOut = new HashMap<>();
		Map<Integer, ArrayList<indexNode>> nodeNeedsIn = new HashMap<>();
		storeCentroid();//read from every cluster
		// pami20 compute radius here, later
		if(pami20){// and not stable
			long startTime1 = System.nanoTime();
			for(int ii=0; ii<k; ii++) {
				if(pami20_flags[ii] == false) {
					radiusPAMI20[ii] = CENTERSEuc.get(ii).pami20_calculate_radius(dataMatrix, pami20_flags);
					numComputeEuc += CENTERSEuc.get(ii).getcoveredPoints().size();
					dataReach += CENTERSEuc.get(ii).getcoveredPoints().size();
				}
			}
			long endtime = System.nanoTime();
			timedistanceiter += (endtime-startTime1)/1000000000.0;
		}
		if(iterationTimes>0 && pami20) {
			double [][]previousCentroidDistance = new double[k][k];
			for(int ii=0; ii<k; ii++) {
				for(int jj=0; jj<k; jj++)
					previousCentroidDistance[ii][jj] = innerCentoridDis[ii][jj];
			}
			pami20_compute_centroid_distance(k, centroidsData, previousCentroidDistance, center_drift, radiusPAMI20);
		}else {
			long startInter = System.nanoTime();
			if(pckmeans)// use knn to accelerate the inner bound for every point
				computeInterCentoridEuckNN(k, CENTERSEuc, centroidsData);// we use knn to accelerate
			else if(assBoundSign || pami20) // lloyd's algorithm will not create
				computeInterCentoridEuc(k, CENTERSEuc, centroidsData);//compute the inter centroid bound martix
			long endInter = System.nanoTime();
			interBoundTime += (endInter-startInter)/1000000000.0;
		}
		regroup(groupNumber);
		if(blockVector || annulus || sdm16)
			computeNormDivisionCentroid(numberBlockVectors, 2);//compute the norm of centroids.
		int boundNumber = groupNumber;//allBound's dimension: boundNumber+2
		if(Hamerly)
			boundNumber = 0; //no bound for each center
		maxdis = new double[k];
		long end = System.nanoTime();
		for (int group_i = 0; group_i < groupNumber && !pckmeanspointboundRecursive; group_i++) {//check each group
			ArrayList<Integer> centers = group.get(group_i);//get the belonging 
			for (int centerID : centers){ //check each center in the group
				Set<indexNode> nodeList = CENTERSEuc.get(centerID).getcoveredNodes();		
				Set<Integer> pointlist = CENTERSEuc.get(centerID).getcoveredPoints();	// the list of points
				if(pointlist.isEmpty() && nodeList.isEmpty()) //there is no point or node in this cluster
					continue;
				//pami20, if cluster is stable, continue
				if(iterationTimes>0 && pami20 && PAMI20ClusterNeedScanningArrayList.get(centerID))
					continue;
				Queue<Object> queue = new LinkedList<>(pointlist);// create a queue to store all the candidates node or point			
				for(indexNode aIndexNode: nodeList)
					queue.add(aIndexNode);//add all the nodes
				while(!queue.isEmpty()) {
					int idx=0;
					Object aObject = queue.poll();
					double [] tra;
					double [] bounds = null;
					double radius = 0;
					
					indexNode node = null;
					if(aObject instanceof Integer) {// point
						idx = (int) aObject;
						tra = accessPointById(idx);
						if(allBounds[idx-1] != null)
							bounds = allBounds[idx-1];
					}else {// node
						nodeReach++;
						node = (indexNode)aObject;
						tra = node.getPivot();
						radius = node.getRadius();
						bounds = node.getBounds();// get the bounds for index node which has considered the father distance 	
					}
					boolean allPrune = false;//initialized as unpruned
					int newCenterId = centerID;//initialize as the original center							
					if(assBoundSign) {//check the bound one by one, use the global, group, and local method	&& aObject instanceof Integer
						double min_dist = Double.MAX_VALUE;//compute the distance with new center					
						if(bounds == null || indexPAMI || usingUpperBound==false) {//indexPAMI means start from root in every iteration.							
							long startTime1 = System.nanoTime();
							min_dist = Util.EuclideanDis(tra, centroidsData[centerID], dimension);
							long endtime = System.nanoTime();	
							timedistanceiter += (endtime-startTime1)/1000000000.0;
							numComputeEuc++;
							boundUpdate++;
							if(aObject instanceof Integer) 	//update the lower bound of point				
								updateSingleLowerBound(idx, group_i, min_dist, boundNumber);
							else// of node
								node.updateSingleLowerBound(idx, group_i, min_dist, boundNumber);
						}else if(center_drift.containsKey(centerID))
							min_dist = bounds[0] + center_drift.get(centerID); // the upper bound																																			
						double newupperbound = min_dist + 2*radius;// tighten the upper bound
						double lowerbound=0;
						if(bounds != null) {
							if(elkan)
								lowerbound = getMinimumLowerbound(bounds, boundNumber, group_i);	// bound from drift
							else if(Hamerly)
								lowerbound = bounds[1] - maxDrift;// use single bound
						}						
						boundCompare++;
						lowerbound = Math.max(lowerbound, interMinimumCentoridDis[centerID]/2.0);//global bounds	
						double second_min_dist = lowerbound;
						if(lowerbound < newupperbound){	//cannot not pass the global filtering							
							if(bounds!=null && indexPAMI==false && usingUpperBound==true) {
								long startTime1 = System.nanoTime();
								min_dist = Util.EuclideanDis(tra, centroidsData[centerID], dimension);
								long endtime = System.nanoTime();	
								timedistanceiter += (endtime-startTime1)/1000000000.0;
								
								numComputeEuc++;
								boundUpdate++;
								if(aObject instanceof Integer) 	//update the lower bound of point				
									updateSingleLowerBound(idx, group_i, min_dist, boundNumber);
								else// of node
									node.updateSingleLowerBound(idx, group_i, min_dist, boundNumber);
								newupperbound = min_dist + 2*radius;// update the upper bound as distance	
							}
							second_min_dist = Double.MAX_VALUE;
							double radiusExp = 2 * newupperbound + interMinimumCentoridDis[centerID];
							double radiusAnnu = newupperbound;// we do not use the second nearest							
							int[] candidate = null, candidate1= null, candidate2=null;
							if(Exponion && aObject instanceof Integer)
								candidate1 = ExponionBound(centerID, radiusExp);//the ICML'16 exp, used to filter the centroids other than 2 nearest.
							if(annulus && aObject instanceof Integer)
								candidate2 = AnnularBound(idx-1, radiusAnnu);
							if(candidate1 == null && candidate2!=null)
								candidate = candidate2;
							else if (candidate1 != null && candidate2==null)
								candidate = candidate1;
							else if (candidate1 != null && candidate2!=null)//if both, we can get an intersection to filter 
								candidate = Util.getIntersection(candidate1, candidate2, k);
							
							for(int group_j=0; group_j < groupNumber; group_j++) {						
								double localbound = innerCentoridDisGroup[centerID][group_j]/2.0;
								if(bounds!=null) {
									double boundupdate = group_drift[group_j];
									if(sdm16 && aObject instanceof Integer && tighterdrift!=null) {
									//	System.out.println(tighterdrift.length);
										boundupdate = Math.min(boundupdate, tighterdrift[centerID][group_j]);
									}
									if(elkan)
										localbound = Math.max((bounds[group_j+2] - boundupdate), innerCentoridDisGroup[centerID][group_j]/2.0);
								}
								double newlowerbound = localbound;//use the last one
								boundCompare++;
								if( localbound < newupperbound) {//cannot pass the group filtering of bound 							
									ArrayList<Integer> centerCandidates = group.get(group_j);
									for(int center_j: centerCandidates) {// goto the local filtering on center in a group, by checking the candidate list and bounds												
										if(center_j == centerID)// has computed
											continue;
										double localbound_center = innerCentoridDis[centerID][center_j]/2.0;
										boundCompare++;										
										if(blockVector && aObject instanceof Integer) {// use the block vector (ICML'16) to tighten the bound
											localbound_center = Math.max(localbound_center, BlockVectorBound(idx-1, center_j, 2));
										}
										if(localbound_center < newupperbound) {//cannot pass the inner centroid bound prunning
											boundCompare++;
											if(candidate != null && candidate[center_j]==0)//pruned by Exp 
												continue;
											long startTime1 = System.nanoTime();
											double dist = Util.EuclideanDis(tra, centroidsData[center_j], dimension);//if it is a point
											long endtime = System.nanoTime();	
											timedistanceiter += (endtime-startTime1)/1000000000.0;// add a timer to count time spent on distance computation
											numComputeEuc++;
											
											if (min_dist > dist) {
												second_min_dist = min_dist;//pass it to the second one
												min_dist = dist; // maintain the one with min distance
												newCenterId = center_j;											
											}else if(dist < second_min_dist && dist != min_dist){
												second_min_dist = dist;
											}
											if(newlowerbound > dist) {// update the bound with minimum distance in the group
												newlowerbound = dist;
											}
										}else {
											numFillocal++;// local pruning											
										}
									}
								}else {
									numFilGroup += group.get(group_j).size();//pruned 
								}
								if(elkan) {// we can also update the global bound as second_min_dist
									boundUpdate++;
									if(node == null)	//update the lower bound of point
										updateSingleLowerBound(idx, group_j, newlowerbound, boundNumber);
									else
										node.updateSingleLowerBound(idx, group_j, newlowerbound, boundNumber);
								}							
							}							
							if(node != null && second_min_dist - min_dist >= 2*radius) {//pruned a node using real distance
								allPrune = true;
							}
							if(sdm16 && maxdis[newCenterId] < min_dist)
								maxdis[newCenterId] = min_dist;
						}else {//global filtering: all prune
							allPrune = true;
							if(node!=null)
								numFilGlobal += k*node.getTotalCoveredPoints();//k centroids are all pruned, we also have a group of points
							else
								numFilGlobal += k;
							if(sdm16 && maxdis[newCenterId] < newupperbound)
								maxdis[newCenterId] = newupperbound;
						}
						boundUpdate++;
						if(node == null) {//update the upper bound.
							allBounds[idx-1][0] = min_dist;// the upper bound
							allBounds[idx-1][1] = second_min_dist;// the local lower bound
						}else {
							node.setUpperbound(min_dist);
							node.setLowerbound(second_min_dist);
						}						
					}else if(iterationTimes>0 && pami20){//assign by pami20
						newCenterId = pami20_assign(tra, centerID, radiusPAMI20);
					}else
					{//brute force if do not use bounds
						if(pckmeans) {// search the two nearest neighbor
							boolean isnode = aObject instanceof indexNode;
							double[] minDistnearestID;
							if(!pckmeansbound)
								minDistnearestID = assignPCKmenas(tra, rootCentroids, centroidsData, isnode);// no bound and optimizations
							else
								minDistnearestID = assignPCKmenasBound(tra, rootCentroids, centroidsData, isnode,node, idx, radius);
							if (isnode) {// how to prune
								newCenterId = (int)minDistnearestID[3] - 1;//the closest one
								if(minDistnearestID[0] - minDistnearestID[2] >= 2 * radius)
									allPrune = true;
								else
									if(node.isLeaf())
										node.setUpperBoundPick(minDistnearestID[2]);//nearest neighbor
									else
										node.setUpperBoundPick(minDistnearestID[0]);//set the node bound here, the second nearest bound
							}else {
								newCenterId = (int)minDistnearestID[1] - 1;//the closest
							}
						}else{
							double min_dist = Double.MAX_VALUE;
							double second_min_dist = 0;
							for (int j = 0; j < k; j++) {
								long startTime1 = System.nanoTime();
								double dist = Util.EuclideanDis(tra, centroidsData[j], dimension);
								long endtime = System.nanoTime();
								timedistanceiter += (endtime - startTime1) / 1000000000.0;
								numComputeEuc++;
								if (min_dist > dist) {
									second_min_dist = min_dist;// pass it to the second one
									min_dist = dist; // maintain the one with min distance, and second min distance
									newCenterId = j;
								} else if (dist < second_min_dist && dist != min_dist) {
									second_min_dist = dist;
								}
							}
							if (second_min_dist - min_dist >= 2 * radius) {// how to prune
								allPrune = true;
							}
						}
					}
					// insert the pami20 assigmnment here
					if(aObject instanceof Integer) {// the point moves to other center, this should be counted into the time of refinement.
						CENTERSEuc.get(centerID).addPointToCluster(idx, tra);//the first iteration
						if(newCenterId != centerID) {
							numeMovedTrajectories++;
							long startTime1 = System.nanoTime();							
							/*		ArrayList<Integer> idxlist;
							if(idxNeedsIn.containsKey(newCenterId))
								idxlist = idxNeedsIn.get(newCenterId);
							else
								idxlist = new ArrayList<Integer>();
							idxlist.add(idx);
							idxNeedsIn.put(newCenterId, idxlist);// temporal store as we cannot add them the trajectory list which will be scanned later, batch remove later
							if(idxNeedsOut.containsKey(centerID))
								idxlist = idxNeedsOut.get(centerID);
							else
								idxlist = new ArrayList<Integer>();
							idxlist.add(idx);
							idxNeedsOut.put(centerID, idxlist);// temporal store, batch remove later*/
							//		System.out.println(newCenterId);
							idxNeedsIn.put(newCenterId, null);
							idxNeedsOut.put(centerID, null);
							CENTERSEuc.get(newCenterId).addSum(tra);
							CENTERSEuc.get(centerID).minusSum(tra);		
							CENTERSEuc.get(newCenterId).addPointToCluster(idx);
							CENTERSEuc.get(centerID).removePointToCluster(idx);
							long endtime = System.nanoTime();
							assigntime[counter] -= (endtime-startTime1)/1000000000.0;
							refinetime[counter] += (endtime-startTime1)/1000000000.0;
						}
					}
					if(aObject instanceof indexNode) {// this is a node
						if (allPrune == false) {//this node cannot be pruned
							if(!node.getNodelist().isEmpty())
								for (indexNode childnode : node.getNodelist()) {
									if(assBoundSign && node.getBounds()!=null) {
										boundUpdate++;
										childnode.setBounds(node.getBounds(), boundNumber);
									}
									if(pckmeans)
										childnode.updateBoundPick(node.getUpperBoundPick()+node.getRadius());
									queue.add(childnode);// push all the child node or point with farther's bounds to the queue,
								}
							else
								for (int childpoint : node.getpointIdList()) {//update the bounds
									if(assBoundSign && node.getBounds()!=null) {
										boundUpdate++;
										allBounds[childpoint-1] = new double[boundNumber+2];
										allBounds[childpoint-1][0] = node.getBounds()[0] + distanceToFather[childpoint-1];
										for(int i=-1; i<boundNumber; i++) {
											allBounds[childpoint-1][i+2] = node.getBounds()[i+2] - distanceToFather[childpoint-1];
										}																			
									}
									if(pckmeans && assigned!=null)
										quantilizedUpperBound[childpoint-1] = node.getUpperBoundPick() + node.getRadius();//distanceToFather[childpoint-1];// this can be reduced to the radius.
									queue.add(childpoint);// push all the child node or point with father's bounds to the queue,
								}
							CENTERSEuc.get(centerID).removeNode(node);//this should be revised, store them temporally.
						}else {
							prunenode++;
							if(newCenterId != centerID) {//the nodes need to leave current node
								CENTERSEuc.get(centerID).removeNode(node);								
								numeMovedTrajectories += node.getTotalCoveredPoints();		
							}
							CENTERSEuc.get(newCenterId).addNode(node);
						}
					}
				}
			}
		}
		if(pckmeanspointboundRecursive)//call our recursive method here
			assignPCKmenasBoundRecursive(root.getPivot(), rootCentroids, centroidsData, true, root, 0, root.getRadius(), Double.MAX_VALUE, idxNeedsIn, idxNeedsOut, nodeNeedsIn);
	//	updatePointsToCluster(idxNeedsIn, idxNeedsOut);
		if(pami20) {
			pami20_update_flag(idxNeedsIn, idxNeedsOut, radiusPAMI20);
		}
		System.out.println("Initialization time: "+(end-start)/1000000000.0+", #computations: "+numComputeEuc + " pruned nodes: "+ prunenode);	
	}
	
	
	
	/*
	 * pami20: A Fast Adaptive k-means with No Bounds
	 * 1) only depends on centroid distance, rank centroids for each centorid
	 * 2) compute the centroids, calculate the distance with local centroid to get the radius first
	 * 3) we integrate the algorithm into our assigment
	 */
	void pami20_compute_centroid_distance(int k, double [][]clustData,
			double [][]previousCentroidDistance, Map<Integer, Double> center_drift, double []radius) {
		for(int i=0; i<k; i++) {
			innerCentoridDis[i] = new double[k];
			double []a = clustData[i];		
			double min = Double.MAX_VALUE;
			for(int j=0; j<k; j++) {				
				if(i!=j) {
					double distance = 0;
					double []b = clustData[j];
					if(previousCentroidDistance[i][j] >= 2*radius[i]+center_drift.get(i) + center_drift.get(j)) {
						boundCompare++;
						distance = previousCentroidDistance[i][j]- center_drift.get(i) - center_drift.get(j);
					}
					else {
						long startTime1 = System.nanoTime();
						distance = Util.EuclideanDis(a, b, a.length);
						long endtime = System.nanoTime();	
						timedistanceiter += (endtime-startTime1)/1000000000.0;					
						numComputeEuc++;
					}
					innerCentoridDis[i][j] = distance;
					if(distance<min) {
						min = distance;
					}
				}
			}
			interMinimumCentoridDis[i] = min;
		}
	}
	
	/*
	 * indicate stable cluster
	 */
	ArrayList<Boolean> pami20_update_flag(Map<Integer, ArrayList<Integer>> idxNeedsIn, Map<Integer, ArrayList<Integer>> idxNeedsOut, double []radius) {
		//check whether points move and centroid changes
		for(int i=0; i<k; i++) {
			boundUpdate++;
			if(!idxNeedsIn.containsKey(i) && !idxNeedsOut.containsKey(i))
				pami20_flags[i] = true;//no point move in or out, stable
			else {
				pami20_flags[i] = false;
			}
		}
		//whether the center needs to scan every points
		PAMI20ClusterNeedScanningArrayList = new ArrayList<Boolean>();
		for(int i=0; i<k; i++) {
			boolean needScanning = true;
			for (int j=0; j<k; j++) {
				if(i != j && innerCentoridDis[i][j] < 2*radius[i])
					if(pami20_flags[j] == false) {
						needScanning = false;
						break;
					}
			}
			boundUpdate++;
			PAMI20ClusterNeedScanningArrayList.add(needScanning);
		}
		return PAMI20ClusterNeedScanningArrayList;//utilize this for each cluster before reading every point inside
	}
	
	int pami20_assign(double []tra, int centerid, double []radius) {
		int newCenterId = centerid;
		long startTime1 = System.nanoTime();
		double localDist = Util.EuclideanDis(tra, centroidsData[centerid], dimension);
		long endtime = System.nanoTime();
		numComputeEuc++;
		timedistanceiter += (endtime-startTime1)/1000000000.0;
		boundCompare+=1;
		if(localDist < 0.5*interMinimumCentoridDis[centerid])
			return centerid;
		double min_dist = localDist;		
		
		for (int j=0; j<k; j++) {
			boundCompare+=2;
			if(centerid == j || innerCentoridDis[centerid][j] > 2*radius[centerid] 
					|| localDist<innerCentoridDis[centerid][j]/2)
				continue;
			startTime1 = System.nanoTime();
			double dist = Util.EuclideanDis(tra, centroidsData[j], dimension);
			endtime = System.nanoTime();	
			timedistanceiter += (endtime-startTime1)/1000000000.0;
			numComputeEuc++;							
			if (min_dist > dist) {
				min_dist = dist; // maintain the one with min distance, and second min distance
				newCenterId = j;
			}
		}
		return newCenterId;
	}
	

	/*
	 * update the points in each cluster
	 */
	public void updatePointsToCluster(Map<Integer, ArrayList<Integer>> idxNeedsIn, Map<Integer, ArrayList<Integer>> idxNeedsOut) {
		for(int idx: idxNeedsIn.keySet()) {
			ArrayList<Integer> idxs = idxNeedsIn.get(idx);
			cluster newCluster = CENTERSEuc.get(idx);
			newCluster.mergePointToCluster(idxs);
		}
		for(int idx: idxNeedsOut.keySet()) {
			ArrayList<Integer> idxs = idxNeedsOut.get(idx);
			cluster newCluster = CENTERSEuc.get(idx);
			newCluster.removePointToCluster(idxs);
		}
	}
	
	/*
	 * update the nodes in each cluster
	 */
	public void updateNodesToCluster(Map<Integer, ArrayList<indexNode>> nodeNeedsIn) {
		for(int idx: nodeNeedsIn.keySet()) {
			ArrayList<indexNode> idxs = nodeNeedsIn.get(idx);
			cluster newCluster = CENTERSEuc.get(idx);
			newCluster.mergeNodesToCluster(idxs);
		}
	}
	
	/*
	 * get the possible two nearest neighbors, ICML16 Newling
	 */
	public int[] ExponionBound(int i,  double radius) {
		int centroidList[] = new int[k];
		for(int i1=0; i1<k; i1++) {
			if(i1==i)
				continue;
			if(innerCentoridDis[i][i1] <= radius)//the distance from the other centroids to current centroid
				centroidList[i1] = 1;
		}
		return centroidList;
	}
	
	/*
	 * rank all centroids based on the norm distance to the origin, return the possible centroids, Drake'13
	 */
	public int[] AnnularBound(int dataID, double radius) {
		int centroidList[] = new int[k];
		for(int i1=0; i1<k; i1++) {
			if(Math.abs(centroidNorm[i1] - dataNorm[dataID])<=radius)
				centroidList[i1] = 1;
		}
		return centroidList;
	}
	
	public Map<Integer, ArrayList<Integer>> GroupsBuilt(){
		groupsBuilt = new HashMap<>();
		int count = 0;
		for(int i=0; i<k; i++) {			
			ArrayList<Integer> arrayList = new ArrayList<Integer>();
			for(int ids: CENTERSEuc.get(i).getcoveredPoints()) {
				arrayList.add(ids-1);
				count++;
			}			
			groupsBuilt.put(i, arrayList);
		}
		System.out.println(count);
		return groupsBuilt;
	}
	
	/*
	 * divide k clusters into t groups when k not equals to t, ICML'15, this can also be called when needs regrouping
	 */
	protected void groupInitialClusters(int t, int k) throws IOException {
		group = new HashMap<>();
		centerGroup = new HashMap<>();
		if(t==k) {// when k is small, we do not divide into too many groups
			for(int i = 0;i<k; i++) {
				ArrayList<Integer> a = new ArrayList<>();
				a.add(i);
				group.put(i, a);
				centerGroup.put(i, i);
			}
		}else {
			String LOG_DIR = "./seeds/Groupcentroid.txt";
			PrintStream fileOut = new PrintStream(LOG_DIR);
			System.setOut(fileOut);	
			for(int i = 0; i<k; i++) {
				double[] aa = CENTERSEuc.get(i).getcentroid();
				int counter1 = 0;
				for(double id: aa) {
					System.out.print(id);
					if(counter1++ < aa.length-1)
						System.out.print(",");
				}
				System.out.println();
			}
			String LOG_DIR1 = "./seeds/Groupcentroid.log";
			PrintStream fileOut1 = new PrintStream(LOG_DIR1);
			System.setOut(fileOut1);
			String args[] = new String[7];
			args[0] = LOG_DIR;
			args[1] = Integer.toString(t);// groups
			args[2] = Integer.toString(k);// number of data points
			args[5] = "0";
			args[6] = Integer.toString(dimension_end-dimension_start);
			kmeansAlgorithm<?> run2 = new kmeansAlgorithm(args);
			run2.staticKmeans(false, false, false);
			group = run2.GroupsBuilt();			
			for(int i = 0; i<t; i++) {
				ArrayList<Integer> centers = group.get(i);
				for(int center: centers) {
					centerGroup.put(center, i);
				}
			}
			System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
		}
		System.out.println("The centroids are grouped");
	}
	
	/*
	 * compute the sum distance after refinement
	 */
	public double getSumDis() {
		double sum = 0;
		Set<Integer> allPoints = new HashSet<>();
		for(int i=0; i<k; i++) {
			if(readOwner)// update the fair algorithm
				sum += CENTERSEuc.get(i).computeFairSum(dataMatrix, allPoints, userID, userNumber);
			else
				sum += CENTERSEuc.get(i).computeSum(dataMatrix, allPoints);
		}
	//	System.out.println("aaaa"+allPoints.size());
		return sum;
	}
	
	public int getPunedNodes() {
		int sum = 0;
		for(int i=0; i<k; i++) {
			sum += CENTERSEuc.get(i).computeNumberofNode();
		}
		return sum;
	}
	
	
	//candidate set stores the keyset of the dataset
	public int runkmeans(int k, Set<Integer> candidateset) throws IOException {
		int groupNumber = k;
		if(k>10 && Yinyang) {//used for the grouping
			groupNumber = k/10;
			if(groupNumber>10)
				groupNumber = 10;
		}
		numComputeEuc = 0;
		boundCompare = 0;
		nodeReach = 0;
		dataReach = 0;
		boundUpdate = 0;
		timedistanceiter = 0;
		groupInitialClusters(groupNumber, k); // Step 1: divide k centroid into t groups
		interMinimumCentoridDis = new double[k];
		innerCentoridDis = new double[k][];
		innerCentoridDisGroup = new double[k][];
		centroidsData = new double[k][];
		centoridData = new ArrayList<>();
		allBounds = new double[trajectoryNumber][];
		center_drift = new HashMap<Integer, Double>();
		group_drift = new double[groupNumber];
		heapBound = new double[k];
		pami20_flags = new boolean[k];
		radiusPAMI20 = new double[k];
		int sortID[][] = null;
		if(blockVector || annulus)
			computeNormDivision(numberBlockVectors, 2);//compute the block vector norm
		if(Drake12)
			sortID = new int[trajectoryNumber][];
		for (int i = 0; i < groupNumber; i++) {
			group_drift[i] = Double.MAX_VALUE;// initialize as max in the begining
		}
		iterationTimes = 0;
		double []iterationtime = new double[TRY_TIMES];
		double []centroidIndexingTime = new double[TRY_TIMES];
		double []iterationdis = new double[TRY_TIMES];
		binSortcenter = (int)((float)k/4);//this can be changed to a number, just like yinyang using 10
		binSortcenter = 10;// just set for pick-means testing, comment it if not
		for(; iterationTimes < TRY_TIMES; iterationTimes++){
			long startTime1 = System.nanoTime();
			if(pckmeans)// we do not maintain any bound across the iterations
				root.setUpperBoundPick(Double.MAX_VALUE);
			if(lloyd || Hamerly || Yinyang || usingIndex || elkan || annulus || Exponion || blockVector || reGroup || sdm16) {
				assignmentBounds(k, groupNumber); // more choice on here
			}else if(heap) {
				assBoundSign = false;
				assignHeapBound();
			}else if(wsdm14) {
				assignWsdm14();
				assBoundSign = false;
			}else if(Drake12) {
				assignSortCenter(sortID);
				assBoundSign = false;
			}else if(pckmeans) {//dual-tree
				recursivePICKmeans(root, rootCentroids, dataMatrix, centroidsData, Double.MAX_VALUE, 0);
			}
			
			long endtime = System.nanoTime();
			assigntime[counter] += (endtime-startTime1)/1000000000.0;
			long startTime = System.nanoTime();
			for(int i=0; i<groupNumber; i++) {
	        	group_drift[i] = 0;// initialize as min
	        }
			double drfitSum = 0;
			for(int i=0; i<k; i++) {
				double drfit = 0;
				if(assBoundSign) {//we use the incremental for bounds only.
					drfit = CENTERSEuc.get(i).extractNewCentroidByMeansIncremental();
				}else if(readOwner){//extract the center based on proportionaly data.
					drfit = CENTERSEuc.get(i).fairkMeansRefinement(userID, userNumber, dataMatrix);
				}else{
					dataReach += CENTERSEuc.get(i).getcoveredPoints().size();
					drfit = CENTERSEuc.get(i).extractNewCentroidByMeans(dataMatrix);
				}
				
				drfitSum += drfit;
				center_drift.put(i, drfit);				
				int groupid = centerGroup.get(i);
				if(group_drift[groupid] < drfit) //update the group drift as maximum
					group_drift[groupid] = drfit;
				if(maxDrift < drfit)
					maxDrift = drfit;
			}
			if(assBoundSign && sdm16 && !lloyd) {
				double[][] centroidsDataold = new double[k][dimension];
				for(int i=0; i<k; i++)
					for(int j=0; j<dimension; j++)
						centroidsDataold[i][j] = centroidsData[i][j];
				storeCentroid();
				tighterDrift(centroidsDataold, centroidsData, maxdis);
			}
			endtime = System.nanoTime();
			
			if(heap) {
				heapUpdateDrift();
			}
			if(Drake12)
				sortCenterUpdateBound(sortID);
			
			if(pckmeans) {// need to refine the index of centroids, very efficient
				long startTime11 = System.nanoTime();
				rootCentroids =  indexkmeans.buildBalltree2(centroidsData, dimension, capacity, userID, userNumber); // capacity
				long endTime11 = System.nanoTime();
				centroidIndexingTime[iterationTimes] = (endTime11-startTime11)/1000000000.0;	
			}
			
			refinetime[counter] +=(endtime-startTime)/1000000000.0;		
			System.out.print("\niteration "+(iterationTimes+1)+", time cost: ");
			System.out.printf("%.5f", (endtime-startTime1)/1000000000.0);
			System.out.println("s");
			iterationtime[iterationTimes] = (endtime-startTime1)/1000000000.0;
			if(drfitSum == 0 || iterationTimes > maxIteration) {//used to terminate as it does not need.
				runrecord.setIterationtimes(iterationTimes+1);
				maxIteration = iterationTimes-1;
				break;//convergence
			}
			if(usingIndex && !indexPAMI && iterationTimes<2) {//scan from the rest
			//	indexPAMI = changeAccessMode();
			}
			if(shownSum) {
				iterationdis[iterationTimes] = getSumDis();
				System.out.println("the sum distance after refinement is: "+iterationdis[iterationTimes]);// this is for verifying the function	
			}
			if(indexPAMI) {//we will empty every cluster, traverse the tree from root again
				boolean first = true;
				for(cluster clus:CENTERSEuc) {
					if(first) {
						clus.reset(root); // reset this as root
						first = false;
					}else {
						clus.reset(null);
					}
				}
			}
			if(pckmeanspointboundRecursive)
				for(cluster clus:CENTERSEuc)
					clus.clearSetPICK();
		}
		System.out.println("#moved points: "+numeMovedTrajectories);
		System.out.println("#computation: "+numComputeEuc);
		System.out.println("#Pruned nodes: "+getPunedNodes());
		
	//	for(int i=0; i<TRY_TIMES; i++)//recoding the running time of each iteration
	//		System.out.println(iterationtime[i]);
	//	System.out.println("Centroid indexing time:");
	//	for(int i=0; i<TRY_TIMES; i++)//recoding the running time of each iteration
	//		System.out.println(centroidIndexingTime[i]);
		System.out.println();
		computations[counter] += numComputeEuc;
		boundCompares[counter] += boundCompare;
		dataAccess[counter] += dataReach + nodeReach; // node's pivot is accessed, also counted as data access
		boundUpdates[counter] += boundUpdate;
		distanceTime[counter] += timedistanceiter;// the time spend on distance computation
		memoryUsage[counter] += getAllMemory(groupNumber, 2);
		System.out.println("Used memory: "+memoryUsage[counter]);
		return iterationTimes;
	}
	
	/*
	 * this is for cascading trick
	 */
	Boolean changeAccessMode() {
		Boolean retravesal = false;
		if(iterationTimes == 0) {//first iteration
			previoustime = time[counter];
		}else {//second iteraiton
			if( time[counter] - previoustime > previoustime) {
				retravesal = true;
			}
		}
		return retravesal;
	}
	
	/*
	 * compute the norm of whole dataset. b is the number of vector, p is the p-norm (2 by default)
	 */
	public void computeNormDivision(int b, int p) {
		dataDivvector = new double[trajectoryNumber][];
		dataNorm = new double[trajectoryNumber];
		int unit = (int)Math.ceil((double)dimension/b);
		if(dimension%2==1)
			unit = unit+1;
		System.out.print("the unit is: "+unit);
		for(int i=0; i<trajectoryNumber; i++) {
			dataNorm[i]=0;
			int m=0;
			dataDivvector[i] =  new double[b];
			for(int j = 0; j<dimension; j++) {
				dataNorm[i] += Math.pow(dataMatrix[i][j], 2);
				if (blockVector) {
					dataDivvector[i][m] += Math.pow(dataMatrix[i][j], p);
					if (j % unit == (unit-1)) {
						dataDivvector[i][m] = Math.pow(dataDivvector[i][m], 1/(float)p);// sqrt p
						m++;
					}
				}
			}
			dataNorm[i] = Math.sqrt(dataNorm[i]);
		}
	}
	
	
	
	/*
	 * compute the norm of whole dataset.b is the number of vector, p is the p-norm (2 by default)
	 */
	public void computeNormDivisionCentroid(int b, int q) {
		centoridDivvector = new double[k][];
		centroidNorm = new double[k];
		int unit = (int)Math.ceil((double)dimension/b);
		if(dimension%2==1)
			unit = unit+1;
		for(int i=0; i<k; i++) {
			int m=0;
			centroidNorm[i] = 0;
			centoridDivvector[i] =  new double[b];
			for(int j = 0; j<dimension; j++) {
				centroidNorm[i] += Math.pow(centroidsData[i][j], 2);
				if (blockVector) {
					centoridDivvector[i][m] += Math.pow(centroidsData[i][j], q);
					if (j % unit == (unit-1)) {
						centoridDivvector[i][m] = Math.pow(centoridDivvector[i][m], 1/(float)q);
						m++;
					}
				}
			}
			centroidNorm[i] = Math.sqrt(centroidNorm[i]);
		}
	}
	
	/*
	 * the product of two vectors, b is the dimension of block vector, 3,4,5 by default.
	 */
	double vectorProduct(int b, double[] vectorA, double[] vectorB) {
		double sum=0;
		for(int i=0; i<b; i++)
			sum += vectorA[i]*vectorB[i];
		return sum;
	}
	
	/*
	 * ICML'16 Speeding up k-means by approximating Euclidean distances via block vectors
	 */
	double BlockVectorBound(int i, int j, int b) {
		return Math.sqrt(Math.pow(dataNorm[i],2) + Math.pow(centroidNorm[j],2) - 
				2*vectorProduct(b, dataDivvector[i], centoridDivvector[j]));
	}
	
	/*
	 * Drake'12, use b bounds instead of k bounds, only maintain the nearest b<k centers.
	 */
	void assignSortCenter(int sortID[][]) {
		Map<Integer, ArrayList<Integer>> cluster_integer = new HashMap<Integer, ArrayList<Integer>>();
		storeCentroid();
		if(allBounds[0] == null) {
			for(int i=0; i<trajectoryNumber; i++) {// initialize all the bounds
				boundUpdate++;
				allBounds[i] = new double[binSortcenter+3];
				allBounds[i][0] = Double.MAX_VALUE;//the upper bound
				for(int j=0; j<binSortcenter+1; j++) {
					allBounds[i][j+2] = 0;
				}
			}
		}
		int m = binSortcenter;// parameter b
		for (cluster aCluster : CENTERSEuc) {			
			for (int trajid : aCluster.getcoveredPoints()) {
				int idx = trajid - 1;
				double min_dist = Double.MAX_VALUE;
				int newCenterId = 0;
				int tunePara = 0;
				boundCompare++;
				if(allBounds[idx][2] < allBounds[idx][0]) {//needs to compute every trajectory
					double []tra = dataMatrix[idx];
					dataReach++;
					for (int a=0; a<k; a++) {	
						long startTime1 = System.nanoTime();
						double dist = Util.EuclideanDis(tra, centroidsData[a], dimension);
						long endtime = System.nanoTime();	
						timedistanceiter += (endtime-startTime1)/1000000000.0;
						
						numComputeEuc++;							
						if (min_dist > dist) {
							min_dist = dist; // maintain the one with min distance, and second min distance
							newCenterId = a;
						}
						if(sortID[idx] == null)
							sortID[idx] = new int[binSortcenter+1];
						boundUpdate++;
						storekBound(sortID[idx], allBounds[idx], a, dist);// update	the b lower bound, 2nd to b sortcenter
					}
					allBounds[idx][0] = min_dist;// set the upper bound
				}else {
					int increment = 0;
					dataReach++;
					double []tra = dataMatrix[idx];
					for (; increment<binSortcenter+1; increment++) {
						if(increment<binSortcenter && allBounds[idx][2+increment] > allBounds[idx][0]) {
							boundCompare++;
							tunePara++;
							continue;//pruned
						}
						int centerID = sortID[idx][increment];
						long startTime1 = System.nanoTime();
						double dist = Util.EuclideanDis(tra, centroidsData[centerID], dimension);
						long endtime = System.nanoTime();	
						timedistanceiter += (endtime-startTime1)/1000000000.0;
						numComputeEuc++;
						if (min_dist > dist) {
							min_dist = dist; // maintain the one with min distance
							newCenterId = centerID;
						}
						boundUpdate++;
						allBounds[idx][increment+2] = dist;
					}
					allBounds[idx][0] = min_dist;// update the upper bound
					storekBound(sortID[idx], allBounds[idx], binSortcenter+1);//resort the bound, no new value is added
				}
				ArrayList<Integer> arrayList = new ArrayList<Integer>();
				if(cluster_integer.containsKey(newCenterId)){
					arrayList = cluster_integer.get(newCenterId);
				}
				arrayList.add(trajid);
				cluster_integer.put(newCenterId, arrayList);// insert pointID, newBound, closest
			//	m = Math.max(m, tunePara);
			}
		}
//		binSortcenter = Math.max(k/8, m);
		for(int i= 0; i<k; i++) {
			CENTERSEuc.get(i).setcoveredPoints(cluster_integer.get(i));		//update the point list in each cluster.
		}
	}
	
	// sort the bound decreasingly, the sortID also needs to be updated
	void storekBound(int sortID[], double bounds[], int length) {
		double temp;
		int tempID;
		for (int i = 1; i < length; i++) {
			for (int j = i; j > 0; j--) {
				if (bounds[j + 2] > bounds[j + 1]) {
					temp = bounds[j + 2];
					tempID = sortID[j];
					bounds[j + 2] = bounds[j + 1];
					sortID[j] = sortID[j-1];
					bounds[j + 1] = temp;
					sortID[j-1] = tempID;
				}
			}
		}
	}
	
	// update the sortID and allBounds.
	void storekBound(int sortID[], double bounds[], int id, double bound) {
		if(id <= binSortcenter) {
			sortID[id] = id;
			bounds[id+2] = bound;
			if(id == binSortcenter)			
				storekBound(sortID, bounds, binSortcenter+1);//descending
		}else {
			int i = 0;
			for (i = 0; i < binSortcenter+1; i++) {
				if (bounds[i + 2] < bound)
					break;
			}			
			if(i>0) {
				for (int j = 0; j < i-1; j++) {
					bounds[j + 2] = bounds[j + 3];
					sortID[j] = sortID[j+1];
				}						
				i--;
				bounds[i + 2] = bound;
				sortID[i] = id;
			}
		}
	}
	
	/*
	 * update the bound after refined based on the drift
	 */
	void sortCenterUpdateBound(int sortID[][]){
		int centerJ = 0;
		double maxDrift = 0;
		for(int center_j=0; center_j<k; center_j++) {
			if(maxDrift < center_drift.get(center_j))
				maxDrift = center_drift.get(center_j);
		}
		for (cluster aCluster : CENTERSEuc) {			
			for (int traid : aCluster.coveredPoints) {
				int idx = traid-1;
				allBounds[idx][0] = allBounds[idx][0] + center_drift.get(centerJ);
				allBounds[idx][2] = allBounds[idx][2] - maxDrift;// the b farthest, minus
				boundUpdate++;
				for(int j = binSortcenter-1; j>=1; j--) {
					double drift = center_drift.get(sortID[idx][j]); 
					allBounds[idx][j+2] = Math.min(allBounds[idx][j+2] - drift, allBounds[idx][j+1]);
				}				
			}			
			centerJ++;
		}
	}	
	
	/*
	 * wsdm 14. conduct the knn for every cluster, and the rest assign to each cluster using baseline.
	 */
	void assignWsdm14() {
		storeCentroid();
		computeInterCentoridEuc(k, CENTERSEuc, centroidsData);//compute the inter centroid bound martix
		HashSet<Integer> assigned = new HashSet<>();// store the assigned 
		for(int i=0; i<k; i++) {
			double centorid[] = CENTERSEuc.get(i).getcentroid();
			double radius = interMinimumCentoridDis[i]/2;// the radius to search
			indexkmeans.setdistanceCompute(0);
			indexkmeans.setNodeAccess(0);
			indexkmeans.setdataAccess(0);
			ArrayList<Integer> result = indexkmeans.SimilaritySearchBall(radius, centorid, root, dimension, dataMatrix);
		//	System.out.println("aaaaa "+indexkmeans.getdistanceCompute());
			numComputeEuc += indexkmeans.getdistanceCompute();
			nodeReach += indexkmeans.getNodeAccess();
			dataReach += indexkmeans.getdataAccess();
			CENTERSEuc.get(i).setcoveredPoints(result);
			assigned.addAll(result);
		}
		System.out.println("#points assigned by similarity search:"+assigned.size());// 
		for(int i=0; i<trajectoryNumber; i++) {
			if(!assigned.contains(i+1)) 
			{
				double min_dist = Double.MAX_VALUE;
				int newCenterId = 0;
				dataReach++;
				double []tra = dataMatrix[i];
				for (int j=0; j<k; j++) {
					long startTime1 = System.nanoTime();
					double dist = Util.EuclideanDis(tra, centroidsData[j], dimension);
					long endtime = System.nanoTime();	
					timedistanceiter += (endtime-startTime1)/1000000000.0;
					
					numComputeEuc++;							
					if (min_dist > dist) {
						min_dist = dist; 
						newCenterId = j;
					}
				}
				CENTERSEuc.get(newCenterId).addPointToCluster(i+1);// the rest will be assigned by computing distance.
			}
		}
	}
	
	/*
	 * 2015 Accelerating Lloyd’s Algorithm for k-Means Clustering, using heaps.
	 */
	void assignHeapBound() {
		storeCentroid();
		numeMovedTrajectories= 0;
		int cluster_id = 0;
		Map<Integer, PriorityQueue<kmeansHeap>> cluster_heap = new HashMap<Integer, PriorityQueue<kmeansHeap>>();//temporal store
		Map<Integer, ArrayList<Integer>> cluster_points = new HashMap<Integer, ArrayList<Integer>>();//temporal store
		for (cluster aCluster : CENTERSEuc) {			
			if (aCluster.emptyheap()) {// put all the points in the cluster into the heap if it is empty
				aCluster.initializeHeap();				
				for (int idx : aCluster.getcoveredPoints()) {
					boundUpdate++;
					kmeansHeap aHeap = new kmeansHeap(idx, -1);
					aCluster.push(aHeap);
				}
			}
			while (!aCluster.emptyheap()) {
				kmeansHeap check = aCluster.poll();
				int pointID = check.getidx();
				double boundGap = check.getBound();
				int closest = -1;
				double newBound = 0;
				boundCompare++;
				if (boundGap > heapBound[cluster_id]) {// stay in current center
					closest = cluster_id;
					newBound = boundGap - center_drift.get(cluster_id);
				} else {
					double tra[] = accessPointById(pointID);
					double min_dis = Double.MAX_VALUE, second_min_dist = 0;					
					for (int center_j = 0; center_j < k; center_j++) {// compute the distance every point
						long startTime1 = System.nanoTime();
						double dist = Util.EuclideanDis(tra, centroidsData[center_j], dimension);
						long endtime = System.nanoTime();	
						timedistanceiter += (endtime-startTime1)/1000000000.0;
						
						
						numComputeEuc++;
						if (min_dis > dist) {// find the closet and second-closed centroid
							second_min_dist = min_dis;
							closest = center_j;
							min_dis = dist;
						} else if (dist < second_min_dist && dist != min_dis) {
							second_min_dist = dist;
						}
					}
					newBound = second_min_dist - min_dis + heapBound[closest];
				}
				boundUpdate++;
				PriorityQueue<kmeansHeap> aHeaps = new PriorityQueue<kmeansHeap>();
				if (cluster_heap.containsKey(closest)) {
					aHeaps = cluster_heap.get(closest);
				}
				aHeaps.add(new kmeansHeap(pointID, newBound));
				cluster_heap.put(closest, aHeaps);// insert pointID, newBound, closest
				
				ArrayList<Integer> aList = new ArrayList<Integer>();
				if (cluster_points.containsKey(closest)) {
					aList = cluster_points.get(closest);
				}
				aList.add(pointID);
				cluster_points.put(closest, aList);// add the point to the closest cluster,
			}
			cluster_id++;
		}
		for(int i= 0; i<k; i++) {
			CENTERSEuc.get(i).setHeap(cluster_heap.get(i));
			CENTERSEuc.get(i).setcoveredPoints(cluster_points.get(i));
		}
	}
	
	
	void heapUpdateDrift() {// update the heap bound after the refinement
		double maxDrift = 0;
		for(int center_j=0; center_j<k; center_j++) {
			if(maxDrift < center_drift.get(center_j))
				maxDrift = center_drift.get(center_j);
		}
		for(int center_j=0; center_j<k; center_j++) {
			boundUpdate++;
			heapBound[center_j] = heapBound[center_j] + center_drift.get(center_j) + maxDrift;
		}
	}
	
	/*
	 * SDM'16, we use it for Elkan according to the algorithm with the tightbound
	 * update the bound which is tighter than the drift, we implement algorithm 2
	 * it can be seen as a tight drift, which can be applied to all the points assigned to c_i before
	 * maxdis is the maximum distance from center to any point in the cluster
	 */
	void tighterDrift(double [][]oldcenter, double [][]newcenter, double maxdis[]) {		
		tighterdrift = new double[k][k];
		for(int i=0; i<k; i++) {
			for(int j=0; j<k; j++) {
				double []gapij = new double[dimension];
				double []gapNewij = new double[dimension];
				double product = 0;
				for(int p=0; p<dimension; p++) {
					gapij[p] = oldcenter[i][p] - oldcenter[j][p];
					gapNewij[p] = newcenter[j][p] - oldcenter[j][p];
					product += gapij[p] * gapNewij[p];
				}
				double tvalue = product/Math.pow(center_drift.get(j),2);
				
				double []pci = new double[dimension];
				for(int p=0; p<dimension; p++) {
					pci[p] += oldcenter[j][p] + gapNewij[p]*tvalue;
				}
				long startTime1 = System.nanoTime();
				double dist = Util.EuclideanDis(pci, oldcenter[i], dimension);
				long endtime = System.nanoTime();	
				timedistanceiter += (endtime-startTime1)/1000000000.0;
				
				numComputeEuc++;
				double cix = dist*2/center_drift.get(j);
				double ciy = 1-2*tvalue;
				double r = maxdis[i]*2/center_drift.get(j);
				double cinorm = Math.pow(centroidNorm[i],2);
				if(cix<=r)
					tighterdrift[i][j] = Math.max(0, Math.min(2, 2*(r-ciy)))*center_drift.get(j)/2;
				else {
					if(ciy>r)
						ciy = ciy - 1;
					tighterdrift[i][j] = 2*(cix*r-ciy*Math.sqrt(cinorm-r*r))/cinorm*center_drift.get(j)/2;
				}
			}
		}
	}
	
	
	/*
	 * get the count of nodes in the tree.
	 */
	public int getNodesCount(indexNode root) {
		if(root == null)
			return 0;
		if(root.isLeaf()) {	
			return 1;
		}else {
			Set<indexNode> listnode = root.getNodelist();
			int max = listnode.size();
			for(indexNode aIndexNode: listnode) {
				max += getNodesCount(aIndexNode);
			}
			return max;
		}
	}
	
	/*
	 * get the count of nodes in the tree.
	 */
	public long getIndexSize(indexNode root) {
		if(root == null)
			return 0;
		if(root.isLeaf()) {	
			return root.getMemory(dimension);
		}else {
			Set<indexNode> listnode = root.getNodelist();
			long max = 0;
			for(indexNode aIndexNode: listnode) {
				max += aIndexNode.getMemory(dimension);
				max += getIndexSize(aIndexNode);
			}
			return max;
		}
	}
	
	/*
	 * get the count of nodes in the lightweight tree.
	 */
	public long getIndexSizePiCK(indexNode root, boolean centorid) {
		if(root == null)
			return 0;
		if(root.isLeaf()) {	
			return root.getMemoryPick(dimension, centorid);
		}else {
			Set<indexNode> listnode = root.getNodelist();
			long max = 0;
			for(indexNode aIndexNode: listnode) {
				max += aIndexNode.getMemoryPick(dimension, centorid);
				max += getIndexSizePiCK(aIndexNode, centorid);
			}
			return max;
		}
	}
	
	/*
	 * initialize the bounds of node
	 */
	public void setBoundsEmpty(indexNode root) {
		root.setBoundsEmpty();
		Set<indexNode> listnode = root.getNodelist();
		for (indexNode aIndexNode : listnode) {
			setBoundsEmpty(aIndexNode);
		}
	}
	
	/*
	 * initialize the bounds of node for regrouping
	 */
	public void setBoundsforRegrouping(indexNode root, int groupNumber) {
		root.setLowerBoundForRegroup(groupNumber);
		Set<indexNode> listnode = root.getNodelist();
		for (indexNode aIndexNode : listnode) {
			setBoundsforRegrouping(aIndexNode, groupNumber);
		}
	}
	
	/*
	 * compute the distance from child (node or point) to father node
	 */
	public void computeFartherToChild(indexNode root) {
		double[] pivot = root.getPivot();
		if(root.isLeaf()) {
			Set<Integer> listpoints = root.getpointIdList();
			for(int pointid: listpoints) {
				double[] childPivot = accessPointById(pointid);
				double distance = Util.EuclideanDis(pivot, childPivot, dimension);
				distanceToFather[pointid-1] = distance; //create a map to store this value
			}
		}else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				double[] childPivot = aIndexNode.getPivot();
				double distance = Util.EuclideanDis(pivot, childPivot, dimension);
				aIndexNode.setdistanceToFarther(distance);
				computeFartherToChild(aIndexNode);
			}
		}
	}
	
	public void staticKmeans(boolean index, boolean bound, boolean tkde02) throws IOException {
		usingIndex = index;
		interBoundTime = 0;
		if(dataMatrix==null)// avoid importing the data every time
			loadDataEuc(datafile, trajectoryNumber);
		if(usingIndex && root == null) {
			root = runIndexbuildQueuePoint(0.01, 20, 10);//radius, capacity, and fanout, we 
		}
		
		if(pckmeans) {//we build centorids for big k.
			assigned = root.setAssignPICK((short)0, assigned);//set all as unassigned
			assigned = new short[trajectoryNumber];
			double centoridMatirx[][]= new double[k][dimension];
			int counter = 0;
			for(int cid: centroids) {
				centoridMatirx[counter++] = dataMatrix[cid-1];
			}
			rootCentroids =  indexkmeans.buildBalltree2(centoridMatirx, dimension, capacity, userID, userNumber); // capacity
			
			if( pckmeanspointbound) {// create to store bound to accelerate knn search
				quantilizedUpperBound = new double[trajectoryNumber];
				for(int i=0; i<trajectoryNumber; i++) {
					quantilizedUpperBound[i] = Double.MAX_VALUE;
				}
			}
		}
		
		initializeClustersRandom(k); 	 //randomly choose k, k-means++, or use index to accelerate.
		printCluster("./logs/pickmeans/"+datafilename+"InitializedCentroids.txt");
		
		assBoundSign = bound;
		indexPAMI = tkde02;
		Set<Integer> KeySet = new HashSet<>();
		for(int i=1; i<=trajectoryNumber; i++)
			KeySet.add(i);
		runkmeans(k, KeySet);
		time[counter] = assigntime[counter] + refinetime[counter];
		if(maxIteration==-1) {//for lloyd algorithm
			time[counter] *= (MAXITE+2);
			assigntime[counter] *= (MAXITE+2);
			refinetime[counter] *= (MAXITE+2);
		}
			
		System.out.println("the overall time is "+time[counter]+"s: "+ assigntime[counter]+" "+refinetime[counter]+"\n");
		System.out.println((time[counter]+indexingTime)+","+ (assigntime[counter]+indexingTime) +","+(indexingTime+interBoundTime)+","+indexingTime+"\n");//
		double datasetspace = trajectoryNumber*dimension*8/(1024.0*1024.0);
		System.out.println((datasetspace+memoryUsage[counter])+","+datasetspace+","+predictedMemory(capacity)+"\n");
		//print the memory cost of dual-tree, with extra bound and bigger tree with capacity 2
		System.out.println("dual-tree cost: "+predictedMemoryDualTree(2)+"\n");
		printCluster("./logs/pickmeans/"+datafilename+"ClusteringResults.txt");
		
		counter++;
		if(root!=null)
			setBoundsEmpty(root); //reset the bound in the node for next test
		cleanMemory();
		runrecord.clear();
		clearStatistics();
	}
	
	/*
	 * we predict the memory cost here, 
	 */
	double predictedMemory(int capacity) {
		double unit = 1024*1024;
		double leafNodeNum = trajectoryNumber/(double)capacity*2; // we assume that node is half filled
		double memorySize = leafNodeNum*(capacity*4+dimension*8+8);
		memorySize += leafNodeNum*(1-Math.pow(2, 1-Math.log10(leafNodeNum)/Math.log10(2)))*(dimension*8+8+4);
		
		if(pckmeans) {//using index to prune
			double centroidLeafNodeNum = k/(double)capacity*2;
			memorySize += centroidLeafNodeNum*(capacity*4+dimension*8+8);
			memorySize += centroidLeafNodeNum*(1-Math.pow(2, 1-Math.log10(centroidLeafNodeNum)/Math.log10(2)))*(dimension*8+8);
		}
		
		memorySize += trajectoryNumber*2 + trajectoryNumber*dimension*8;
		
		if(!pckmeansUsinginterbound)// not using inter bound
			memorySize -= trajectoryNumber*2;
		
		return memorySize/unit;
	}

	/*
	 * we predict the memory cost here,
	 */
	double predictedMemoryDualTree(int capacity) {
		double unit = 1024*1024;
		double leafNodeNum = trajectoryNumber/(double)capacity; // we assume that node is half filled
		double memorySize = leafNodeNum*(capacity*4+dimension*8+8+2*8+4);
		memorySize += leafNodeNum*(1-Math.pow(2, 1-Math.log10(leafNodeNum)/Math.log10(2)))*(dimension*8+8+4+2*8);
		
		if(pckmeans) {//using index to prune
			double centroidLeafNodeNum = k/(double)capacity;
			memorySize += centroidLeafNodeNum*(capacity*4+dimension*8+8);
			memorySize += centroidLeafNodeNum*(1-Math.pow(2, 1-Math.log10(centroidLeafNodeNum)/Math.log10(2)))*(dimension*8+8);
		}
		
		memorySize += trajectoryNumber*2*8;// + trajectoryNumber*dimension*8;//upper bound and lower bound
		
		return memorySize/unit;
	}
	
	public void skipMethods() {//this is for wsdm, sortcenter, 
	//	time[counter] = Double.MAX_VALUE;
		time[counter] = 0;
		assigntime[counter] = 0;
		refinetime[counter] = 0;
		interBoundTime = 0;
		counter++;
	}
	
	void InitializeAllUnifiedCentroid(int maximumk, int group) {
		allCentroids = new int[group][];
		Random rand = new Random();
		for(int i=0; i<group; i++) {
			if(kmeansplusplus) {
				long startTime1 = System.nanoTime();
				allCentroids[i] = kmeanplusplus(dimension, maximumk, trajectoryNumber);//use the SODA 07
				long endtime = System.nanoTime();	
				System.out.println("kmeans++ costs: "+(endtime-startTime1)/1000000000.0);
			}else {
				allCentroids[i] = new int[maximumk]; //it is initialized to use a same centroid set.	
				for(int j=0; j<maximumk; j++) {
					allCentroids[i][j] = rand.nextInt(trajectoryNumber)+1;
				}
			}
		}
	}
	
	/*
	 * write multiple metrics into log files
	 */
	void writelogs(int testTime, String indexname, String tab) throws FileNotFoundException {
	//	String LOG_DIR = "./logs/vldb_logs1_pami20/"+datafilename+"_"+trajectoryNumber+"_"+dimension+"_"+k+"_"+indexname+"_"+capacity+".log";
		String LOG_DIR = "./logs/pickmeans/"+datafilename+"_"+trajectoryNumber+"_"+dimension+"_"+k+"_"+indexname+"_"+capacity+".log";
		PrintStream fileOut = new PrintStream(LOG_DIR);
		System.setOut(fileOut);			
		for(int i=0; i<numberofComparison; i++) {			// show the time speedup over Lloyd algorithm
			System.out.printf("%.2f",time[i]/testTime);
			System.out.print(tab);
		}
		System.out.println();
		for(int i=0; i<numberofComparison; i++) {			// show the time speedup over Lloyd algorithm
			if(time[i]>0.0)
				System.out.printf("%.2f",1/(time[i]/time[0]));
			else
				System.out.printf("0");
			System.out.print(tab);	
		}
		System.out.println();
		System.out.println();
		
		for(int i=0; i<numberofComparison; i++) {			// show the time speedup over Lloyd algorithm
			System.out.printf("%.2f",assigntime[i]/testTime);
			System.out.print(tab);
		}
		System.out.println();
		for(int i=0; i<numberofComparison; i++) {			// show the time speedup over Lloyd algorithm
			if(assigntime[i]>0.0)
				System.out.printf("%.2f",1/(assigntime[i]/assigntime[0]));
			else
				System.out.printf("0");
			System.out.print(tab);	
		}
		System.out.println();
		System.out.println();
		
		for(int i=0; i<numberofComparison; i++) {			// show the time speedup over Lloyd algorithm
			System.out.printf("%.2f",refinetime[i]/testTime);
			System.out.print(tab);
		}
		System.out.println();
		for(int i=0; i<numberofComparison; i++) {			// show the time speedup over Lloyd algorithm
			if(refinetime[i]>0.0)
				System.out.printf("%.2f",1/(refinetime[i]/refinetime[0]));
			else
				System.out.printf("0");
			System.out.print(tab);	
		}
		System.out.println();
		System.out.println();
		
		for(int i=0; i<numberofComparison; i++){		// show the #computation speedup over Lloyd algorithm				
			System.out.print(computations[i]/testTime);
			System.out.print(tab);	
		}
		System.out.println();
		for(int i=0; i<numberofComparison; i++) {			// show the #computation speedup over Lloyd algorithm
			System.out.printf("%.2f",1-computations[i]/(double)computations[0]);
			System.out.print(tab);
		}
		
		System.out.println();
		System.out.println();
		for(int i=0; i<numberofComparison; i++){		// show the #computation speedup over Lloyd algorithm				
			System.out.print(boundCompares[i]/testTime);
			System.out.print(tab);	
		}
		
		System.out.println();
		System.out.println();
		for(int i=0; i<numberofComparison; i++){	// show the #computation speedup over Lloyd algorithm				
			System.out.print(dataAccess[i]/testTime);
			System.out.print(tab);
		}
		
		System.out.println();
		System.out.println();
		for(int i=0; i<numberofComparison; i++){		// show the #computation speedup over Lloyd algorithm				
			System.out.print(boundUpdates[i]/testTime);
			System.out.print(tab);
		}
		
		System.out.println();
		System.out.println();
		for(int i=0; i<numberofComparison; i++){		// show the #computation speedup over Lloyd algorithm				
			System.out.printf("%.2f", memoryUsage[i]/testTime);
			System.out.print(tab);
		}
		
		System.out.println();
		System.out.println();
		for(int i=0; i<numberofComparison; i++){		// show the #computation speedup over Lloyd algorithm				
			System.out.printf("%.2f", distanceTime[i]/testTime);
			System.out.print(tab);
		}
		
		
		for(int i=2; i<numberofComparison; i++) {//only maintain the Lloyd's and Sequential as they will not be affected by the type of index
			computations[i] = 0;
			time[i] = 0;
			assigntime[i] = 0;
			refinetime[i] = 0;
			dataAccess[i] = 0;
			boundCompares[i] = 0;
			boundUpdates[i] = 0;
			memoryUsage[i] = 0;
			distanceTime[i] = 0;
		}
	}
	
	void cleanMemory() {
		allBounds = null;
		centoridData = null;
		centoridDivvector = null;
		centroidNorm = null;
		heapBound = null;
		maxdis = null;
		tighterdrift = null;
		interMinimumCentoridDis = null;
		innerCentoridDis = null;
		innerCentoridDisGroup = null;
		centroidsData = null;
		trajectoryBounds = null;
		center_drift =null;
		group_drift = null;
	}
	
	/*
	 * compute the space usage of bounds
	 */
	double getAllMemory(int groupNumber, int b) {
		double all = 0;
		double unit = 1024.0*1024.0;
		if(assBoundSign) {
			all += k*8/unit;//center drift
		}
		if(usingIndex || wsdm14) {//index size
			all += getIndexSize(root)/unit;
			if(pckmeanspointboundRecursive) {
				all = getIndexSizePiCK(root, false)/unit;
				if(pckmeans)
					all += getIndexSizePiCK(rootCentroids, true)/unit;
				if(pckmeansUsinginterbound)
					all += trajectoryNumber*2/unit;
			}
		}
		if(usingIndex && assBoundSign)
			all += distanceToFather.length*8/unit;
		if(elkan) {
			all += k*8/unit;
			all += k*k*8/2/unit;
			if(Yinyang) {
				all += (long)trajectoryNumber*(groupNumber+2)*8/unit;
				all += groupNumber*8/unit;	//group drift
				all += groupNumber*groupNumber*8/2/unit;
				if(blockVector || annulus) {
					all += k*8/unit; // the centorid norm
				}
				if(blockVector) {
					all += (long)trajectoryNumber*8/unit; // norm of data
					all += (long)trajectoryNumber*b*8/unit; // block vector of data
					all += k*b*8/unit; // block vector of centroid
				}
			}else {
				all += (long)trajectoryNumber*(k+2)*8/unit; //elkan's full bound
			}
			if(sdm16) {
				all += (long)trajectoryNumber*8/unit; // the size of norm;
				all += k*k*8/unit; //tighter drift
			}
		}
		if(Hamerly)
			all += (long)trajectoryNumber*2*8/unit; // two bound size
		if(Drake12) {
			all += (binSortcenter+1)*(long)trajectoryNumber*8/unit; //sort
			all += binSortcenter*(long)trajectoryNumber*4/unit; //sortID
		}
		if(heap) {
			all += k*8/unit;//center's upper bound
			all += (long)trajectoryNumber*8/unit; //the single bound
		}
		if(pami20) {
			all += k*8/unit;//center drift
			all += k*k*8/unit;//inner centroid distance
			all += k*8/unit;
			all += k*8/unit;
			all += k/unit;
		}
		return all;
	}
	
	public void setrunIndexMethodOnly() {
		runIndexMethodOnly = true;
	}
	
	public void setTestBallTree() {
		runBalltreeOnly = true;
	}
	
	
	/*
	 * for analytic use, no need to cover in the unik.
	 */
	public void printCluster(String filename) {
		for(cluster clus: CENTERSEuc) {
			double[] aDoubleWs = clus.getcentroid();
			String contentString = "";
			for(int dim = 0; dim<dimension; dim++)
				contentString += aDoubleWs[dim]+",";
			Util.write(filename, contentString+"\n");
		}
	}
	
	/*
	 * for analytic use, no need to cover in the unik.
	 */
	public void printCluster() {
		for(cluster clus: CENTERSEuc) {
			double[] aDoubleWs = clus.getcentroid();
			for(int dim = 0; dim<dimension; dim++)
				System.out.print(aDoubleWs[dim]+", ");
			System.out.println();
		}
		saveClusterAsFile();
	}
	
	
	/*
	 * output each cluster as a file
	 */
	public void saveClusterAsFile() {
		int counter = 1;
		int point_id = 1;
		Map<Integer, String> id_pointsMap = new HashMap<Integer, String>();
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(datafile)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the point by their id.
				String str = in.nextLine();
				id_pointsMap.put(point_id++, str);
			}
			in.close();
		}		
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}		
		
		for(cluster clus: CENTERSEuc) {
			double[] aDoubleWs = clus.getcentroid();
			Set<Integer> pointSet = clus.getcoveredPoints();
			String filenameString = "./logs/"+datafilename+"/"+String.valueOf(k)+"_"+String.valueOf(counter)+".txt";
			Util.write(filenameString, pointSet.size()+","+aDoubleWs[0]+","+aDoubleWs[1]+"\n");
			for(int pointid: pointSet) {
				Util.write(filenameString, id_pointsMap.get(pointid)+"\n");
			}
			counter++;
		}
	}
	
	
	/*
	 * we test different index, and run the pure index
	 */
	public void experiments_index(int []setK, int testTime) throws IOException, KeySizeException, KeyDuplicateException {
		String LOG_DIR = "./logs/vldb_logs1/"+datafilename+"_"+trajectoryNumber+"_"+dimension+"_"+capacity+"_index.log";
		System.out.println("aaaaa");
	//	plotData.runPlot("");
		loadDataEuc(datafile, trajectoryNumber);	// load the data and create index
		indexkmeans = new indexAlgorithm<>();
		indexNode rootHKT=null, rootMtree=null, rootBall=null, rootCover=null, rootkd = null;
		
		String filenameString = "";
		long startTime1 = System.nanoTime();
		rootHKT = runIndexbuildQueuePoint(0, capacity, 10);//load the dataset and build one index for all testing methods
		long endtime = System.nanoTime();	
		System.out.println("HKT Indexing time");
		System.out.print((endtime-startTime1)/1000000000.0+ ",");
		Util.rewrite(LOG_DIR, "HKT: "+(endtime-startTime1)/1000000000.0+ "\n");
		
		startTime1 = System.nanoTime();
		rootCover = indexkmeans.buildCoverTree(dimension, dataMatrix);
		endtime = System.nanoTime();	
		System.out.println("cover-tree Indexing time");
		System.out.print((endtime-startTime1)/1000000000.0+ ",");
		Util.rewrite(LOG_DIR, "Cover: "+(endtime-startTime1)/1000000000.0+ "\n");
		System.out.println("Cover-tree Indexing space:"+getIndexSize(rootCover)/(1024.0*1024.0));
		
		startTime1 = System.nanoTime();
		rootBall = indexkmeans.buildBalltree2(dataMatrix, dimension, capacity, null, null); // capacity
		endtime = System.nanoTime();	
		System.out.println("ball-tree Indexing time");
		System.out.print((endtime-startTime1)/1000000000.0+ ",");
		Util.rewrite(LOG_DIR, "Ball: "+(endtime-startTime1)/1000000000.0+ "\n");
		System.out.println("ball-tree Indexing space:"+getIndexSize(rootBall)/(1024.0*1024.0));
		
	//	Util.write(LOG_DIR, datafilename+","+dimension+","+trajectoryNumber+", Ball-tree: "
	//			+getIndexSize(rootBall)/(1024.0*1024.0)+", Cover-tree: "+getIndexSize(rootCover)/(1024.0*1024.0)+"\n");
		
		startTime1 = System.nanoTime();
		indexkmeans.buildKDtree(dimension, dataMatrix);//
		endtime = System.nanoTime();	
		System.out.println("kd-tree Indexing time");
		System.out.print((endtime-startTime1)/1000000000.0+ ",");
		Util.rewrite(LOG_DIR, "kd: "+(endtime-startTime1)/1000000000.0+ "\n");
		
		startTime1 = System.nanoTime();
	//	rootMtree = indexkmeans.buildMtree(dataMatrix, dimension, capacity); // m-tree too slow.
		endtime = System.nanoTime();	
		System.out.println("m-tree Indexing time");
		System.out.print((endtime-startTime1)/1000000000.0+ ",");
		Util.rewrite(LOG_DIR, "M-tree: "+(endtime-startTime1)/1000000000.0+ "\n\n");
	//	setK = null;//only test construction
		
		for(int kvalue: setK) {//test various k
			LOG_DIR = "./logs/vldb_logs1/"+datafilename+"_"+trajectoryNumber+"_"+dimension+"_"+k+"_"+capacity+"_clustering.log";

			k = kvalue;
			ArrayList<Pair<indexNode, String>> roots =  new ArrayList<>();
			roots.add(new ImmutablePair<>(rootHKT, "HKT"));
			roots.add(new ImmutablePair<>(rootBall, "BallMetric"));
		//	roots.add(new ImmutablePair<>(rootMtree, "Mtree"));	//too slow
			roots.add(new ImmutablePair<>(rootCover, "CoverTree"));
			InitializeAllUnifiedCentroid(kvalue, testTime);//maximum k, time time
			String content = "";
			for(Pair<indexNode, String> newroot: roots) { 
				root = newroot.getLeft();
				if(root==null)
					continue;
				distanceToFather = new double[trajectoryNumber];//store the point distance
				computeFartherToChild(root);
				String indexname = newroot.getRight();
				if(!indexname.equals("HKT"))
					nonkmeansTree = false;		
				
				System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
				for (int testtime = 0; testtime < testTime; testtime++) {				
					counter = 0;
					centroids = new int[k];
					for(int i=0; i<k; i++)
						centroids[i] = allCentroids[testtime][i];							
					int signBall[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
					setSign(signBall);
					startTime1 = System.nanoTime();
					staticKmeans(true, false, true);
					endtime = System.nanoTime();	
					System.out.println(kvalue+","+dimension+","+trajectoryNumber+","+indexname+" time");
					System.out.print((endtime-startTime1)/1000000000.0+ ",\n");
					content += indexname + ": " +(endtime-startTime1)/1000000000.0+ "\n";
				}
			}
			Util.rewrite(LOG_DIR, content+"\n");
		}
	}
	
	
	/*
	 * we test different baselines, parameters, 
	 */
	public void experiments(int []setK, int testTime) throws IOException, KeySizeException, KeyDuplicateException {
	//	plotData.runPlot("");
		loadDataEuc(datafile, trajectoryNumber);	// load the data and create index
		indexkmeans = new indexAlgorithm<>();
		indexNode rootHKT=null, rootMtree=null, rootBall=null, rootCover=null, rootkd = null;
		if(runBalltreeOnly)
			rootHKT = runIndexbuildQueuePoint(0, capacity, 10);//load the dataset and build one index for all testing methods
	//	String LOG_DIR = "./logs/vldb_logs1_pami20/"+datafilename+"_"+trajectoryNumber+"_"+dimension+"_"+capacity+"_index.log";
		String LOG_DIR = "./logs/pickmeans/"+datafilename+"_"+trajectoryNumber+"_"+dimension+"_"+capacity+"_index.log";
		PrintStream fileOut = new PrintStream(LOG_DIR);
	//	System.setOut(fileOut);	
		long startTime1 = System.nanoTime();
		rootBall = indexkmeans.buildBalltree2(dataMatrix, dimension, capacity, userID, userNumber); // capacity
		long endtime = System.nanoTime();
		indexingTime = (endtime-startTime1)/1000000000.0;
		ArrayList<Double> raidus = new ArrayList<Double>();
		ArrayList<Double> fatherdis = new ArrayList<Double>();
		ArrayList<Double> numPoints = new ArrayList<Double>();
		ArrayList<Double> nodeDepth = new ArrayList<Double>();
		
			
	//	indexkmeans.buildKDtree(dimension, dataMatrix);//
	//	rootMtree = indexkmeans.buildMtree(dataMatrix, dimension, capacity); // m-tree too slow.
	//	rootCover = indexkmeans.buildCoverTree(dimension, dataMatrix);
	//	System.out.println("aaa: "+rootCover.getSum()[0]);
	//	System.out.println("aaa: "+rootBall.getSum()[0]);
	//	System.out.println("aaa: "+indexkmeans.getcount(rootCover));
	//	String groundTruthString = "dataset,k,dimension,trajectoryNumber,height,leafnode,radiusMean,radiusSD,fastest\n";//"grund-truth\n";
		String groundTruthString ="";
		System.out.println(getNodesCount(rootHKT)+" "+ getNodesCount(rootBall) + " "+getNodesCount(rootMtree));
		for(int i=0; i<numberofComparison; i++) {//only maintain the Lloyd's and Sequential as they will not be affected by the type of index
			computations[i] = 0;
			time[i] = 0;
			assigntime[i] = 0;
			refinetime[i] = 0;
			dataAccess[i] = 0;
			boundCompares[i] = 0;
			memoryUsage[i] = 0;
			distanceTime[i] = 0;
		}
		System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
		long starttime = System.nanoTime();
		for(int kvalue: setK) {//test various k
			k = kvalue;
			ArrayList<Pair<indexNode, String>> roots =  new ArrayList<>();
			roots.add(new ImmutablePair<>(rootHKT, "HKT"));
			roots.add(new ImmutablePair<>(rootBall, "BallMetric"));
			roots.add(new ImmutablePair<>(rootMtree, "Mtree"));	
			roots.add(new ImmutablePair<>(rootCover, "CoverTree"));
			InitializeAllUnifiedCentroid(kvalue, testTime);//maximum k, time time
			boolean LloydandSeq = runIndexMethodOnly;
			for(Pair<indexNode, String> newroot: roots) { 
				root = newroot.getLeft();
				if(root==null)
					continue;
				distanceToFather = new double[trajectoryNumber];//store the point distance
				computeFartherToChild(root);
				double leafnode = indexkmeans.getLeafRadius(root, raidus, rootBall.getRadius(),fatherdis,numPoints, nodeDepth, 0)/(trajectoryNumber/(double)capacity);

				String indexname = newroot.getRight();
				if(!indexname.equals("HKT"))
					nonkmeansTree = false;			
				System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
				for (int testtime = 0; testtime < testTime; testtime++) {				
					counter = 0;
					centroids = new int[k];
					for(int i=0; i<k; i++)
						centroids[i] = allCentroids[testtime][i];

					//lloyd, usingIndex, elkan, Hamerly, Drake12, annulus, wsdm14, heap, Yinyang, Exponion, blockVector, reGroup, sdm16
					int[] best_config = new int[13];
					double min_time = 0;
					if(!LloydandSeq) {
					//	min_time = testExisting(best_config);
					//	min_time = testPAMI20(best_config);
					//	min_time = testExistingSelective(best_config);
					//	min_time = testUndiscover(best_config, min_time);// for an future direction
					//	min_time = testUniKSelective(best_config); //test representative
						test_PICKmeans();
					}					
					counter = 13;//skip the first two as the processing is the same, index might be different as the capacity may be different					
				//	testIndex(best_config, min_time);
				//	testPickmeansVariousSetting();
					//run the algorithm with the optimal configuration, as another task to do another classification if it beats both, this
					maxIteration = MAXITE;
					// our own method, using cascading method, all different 
					System.out.println("k:"+k+",capacity:"+capacity+"data scale:"+trajectoryNumber+", test time:"+testtime+"============================= ends here!");//write into logs for comparison						
				}
			//	rankings(raidus, fatherdis, numPoints, nodeDepth, rootBall, leafnode, LloydandSeq, testTime, indexname);
				LloydandSeq = true;
				if(testTime>0)
					writelogs(testTime, indexname, "\t");
			}
		}
	//	Util.write("groundtruth2.csv", groundTruthString);
		long end = System.nanoTime();
	//	System.out.println("evaluation time: "+(end-starttime)/1000000000.0);
	//	Util.write("./logs/vldb_logs1_pami20/"+"time.txt", (end-starttime)/1000000000.0+","+ datafilename+"\n");
		Util.write("./logs/pickmeans/"+"time.txt", (end-starttime)/1000000000.0+","+ datafilename+"\n");
		counter = 0;
	}
	
	void rankings(ArrayList<Double> raidus, ArrayList<Double> fatherdis, ArrayList<Double> numPoints, ArrayList<Double> nodeDepth,
			indexNode rootBall, double leafnode, boolean LloydandSeq, int testTime, String indexname ) throws FileNotFoundException {
		String aString = generateRank(time, numberofComparison);
		DecimalFormat dec = new DecimalFormat("#0.0000");
		double mean = Util.calculateMean(raidus);//
		double depth_mean = Util.calculateMean(nodeDepth);//
		double depth_sd = Util.calculateSD(nodeDepth, depth_mean);//
		double coverPointsMean = Util.calculateMean(numPoints);//
		double coverPointsSD = Util.calculateSD(numPoints, coverPointsMean);//
		double disFatherMean = Util.calculateMean(fatherdis);//
		double disFatherSD = Util.calculateSD(fatherdis, disFatherMean);//
		
		int expectedHeight = (int) (Math.log(trajectoryNumber/(float)capacity) / Math.log(2.0));
		if(rootBall.getRadius()>0) {
			String aaString = datafilename+","+k+","+dimension+","+trajectoryNumber+","+indexkmeans.getHeight(rootBall)/(float)expectedHeight+
				","+dec.format(leafnode)+","+mean+","+Util.calculateSD(raidus, mean)+","+aString+"\n";
		//	groundTruthString += aaString;
		//	Util.write("groundtruth1.csv", aaString);
			//dataset,k,dimension,trajectoryNumber,height,nodecount,leafnode,depthmean,depthsd,radiusmean,radiussd,pointsmean,pointssd,disfathermean,disfathersd
			String content = datafilename+","+k+","+dimension+","+trajectoryNumber+","+indexkmeans.getHeight(rootBall)/(float)expectedHeight+","+getNodesCount(rootBall)/(trajectoryNumber/(double)capacity)+
					","+dec.format(leafnode)+","+dec.format(depth_mean/expectedHeight)+","+dec.format(depth_sd/expectedHeight)+","
					+mean+","+Util.calculateSD(raidus, mean)+","+
					dec.format(coverPointsMean/capacity)+","+dec.format(coverPointsSD/capacity)+","+
				dec.format(disFatherMean)+","+dec.format(disFatherSD);
		//	Util.write("./logs/vldb_logs1_pami20/"+"groundtruth_rank.csv", content+","+aString+"\n");
			Util.write("./logs/pickmeans/"+"groundtruth_rank.csv", content+","+aString+"\n");
		//	Util.write("groundtruth_align.csv", content+"\n");
		}
	}
	
	/*
	 * test selective sequential settings that perform fast 
	 */
	double testExistingSelective(int[] best_sign) throws IOException {
		double min = Double.MAX_VALUE;
		int sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("Lloyd");			//2			
		setSign(sign);
		staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		min = getOptimalConfig(min, sign, best_sign);
		skipMethods();
		
		System.out.println("Sequential-Hamerly");
		int signharmly[] = {0,0,0,1,0,0,0,0,0,0,0,0,0};
		setSign(signharmly);
		staticKmeans(false, true, false); 
		min = getOptimalConfig(min, signharmly, best_sign);
																		
		System.out.println("Drake12-SortCenter");
		int signdrake[] = {0,0,0,0,1,0,0,0,0,0,0,0,0};
		setSign(signdrake);
		staticKmeans(false, true, false);//test the drake12
	//	skipMethods();
		min = getOptimalConfig(min, signdrake, best_sign);
		
		skipMethods();
		skipMethods();
								
		System.out.println("Heap-15");// heap
		int signheap[] = {0,0,0,0,0,0,0,1,0,0,0,0,0};
		setSign(signheap);
		staticKmeans(false, true, false);//test the heap	
		min = getOptimalConfig(min, signheap, best_sign);
		
		System.out.println("Sequential-Yinyang");//yinyang twice
		int signYinyang[] = {0,0,1,0,0,0,0,0,1,0,0,0,0};
		setSign(signYinyang);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signYinyang, best_sign);
		
		skipMethods();
		skipMethods();
		
		System.out.println("Sequential-regroup");//parameter, 
		int signregroup[] = {0,0,1,0,0,0,0,0,1,0,0,1,0};
		setSign(signregroup);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signregroup, best_sign);
		
		skipMethods();					
		skipMethods();
		return min;
	}
	
	/*
	 * test four settings: lloyds, yinyang, index and unik, for table 6
	 */
	double testUniKSelective(int[] best_sign) throws IOException {
		double min = Double.MAX_VALUE;
		int sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("Lloyd");			//2			
		setSign(sign);
		staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		min = getOptimalConfig(min, sign, best_sign);
		skipMethods();
		
		skipMethods();
		
		skipMethods();
		
		skipMethods();
		skipMethods();
								
		skipMethods();
		
		System.out.println("Sequential-Yinyang");//yinyang twice
		int signYinyang[] = {0,0,1,0,0,0,0,0,1,0,0,0,0};
		setSign(signYinyang);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signYinyang, best_sign);
		
		skipMethods();
		skipMethods();
		
		skipMethods();
		
		skipMethods();					
		skipMethods();
		return min;
	}
	
	
	/*
	 * test existing settings of sequential methods.
	 */
	double testExisting(int[] best_sign) throws IOException {
		double min = Double.MAX_VALUE;
		int sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("Lloyd");			//2			
		setSign(sign);
		staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		min = getOptimalConfig(min, sign, best_sign);
		
		System.out.println("Sequential-elkan");//3
		int signelkan[] = {0,0,1,0,0,0,0,0,0,0,0,0,0};
		setSign(signelkan);
		staticKmeans(false, true, false); 
		min = getOptimalConfig(min, signelkan, best_sign);
		
		System.out.println("Sequential-Hamerly");//4
		int signharmly[] = {0,0,0,1,0,0,0,0,0,0,0,0,0};
		setSign(signharmly);
		staticKmeans(false, true, false); 
		min = getOptimalConfig(min, signharmly, best_sign);
																		
		System.out.println("Drake12-SortCenter");//5
		int signdrake[] = {0,0,0,0,1,0,0,0,0,0,0,0,0};
		setSign(signdrake);
		staticKmeans(false, true, false);//test the drake12
	//	skipMethods();
		min = getOptimalConfig(min, signdrake, best_sign);
		
		System.out.println("Sequential-annulus");//6
		int signannulus[] = {0,0,1,0,0,1,0,0,0,0,0,0,0};
		setSign(signannulus);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signannulus, best_sign);
		
		System.out.println("WSDM14");//7
		int signwsdm[] = {0,0,0,0,0,0,1,0,0,0,0,0,0};
		setSign(signwsdm);
		staticKmeans(false, true, false);//test the wsdm14
	//	skipMethods();
		min = getOptimalConfig(min, signwsdm, best_sign);
								
		System.out.println("Heap-15");// heap 8
		int signheap[] = {0,0,0,0,0,0,0,1,0,0,0,0,0};
		setSign(signheap);
		staticKmeans(false, true, false);//test the heap	
		min = getOptimalConfig(min, signheap, best_sign);
		
		System.out.println("Sequential-Yinyang");//yinyang twice 9
		int signYinyang[] = {0,0,1,0,0,0,0,0,1,0,0,0,0};
		setSign(signYinyang);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signYinyang, best_sign);
		
		System.out.println("Sequential-exponion");//2016 10
		int signexponion[] = {0,0,1,0,0,0,0,0,1,1,0,0,0};
		setSign(signexponion);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signexponion, best_sign);
		
		System.out.println("Sequential-block vector");//2016 11
		int signblockvector[] = {0,0,1,0,0,0,0,0,1,0,1,0,0};
		setSign(signblockvector);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signblockvector, best_sign);
		
		System.out.println("Sequential-regroup");//parameter, 12
		int signregroup[] = {0,0,1,0,0,0,0,0,1,0,0,1,0};
		setSign(signregroup);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signregroup, best_sign);
		
		System.out.println("Sequential-sdm16");// we use the tighter drift 13
		int signsdm16[] = {0,0,1,0,0,0,0,0,0,0,0,0,1};
		setSign(signsdm16);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signsdm16, best_sign);
								
		System.out.println("Sequential-full optimizations");// 14 test all the function integrated in one, except sortCenter, wsdm, heap
		int signFull[] = {0,0,1,0,0,1,0,0,1,1,1,1,1};
		setSign(signFull);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, signFull, best_sign);
		
		// add our new algorithm: pck-means, use index to accelerate

		return min;
	}
	
	double testPAMI20(int[] best_sign) throws IOException {
		// we only run lloyd algorithm and pami20
		double min = Double.MAX_VALUE;
		int sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		pami20 = false;
		System.out.println("Lloyd");			//2			
		setSign(sign);
		staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		min = getOptimalConfig(min, sign, best_sign);
		
		System.out.println("PAMI20");			//2	
		pami20 = true;
		staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		min = getOptimalConfig(min, sign, best_sign);
		
		pami20 = false;
		
		skipMethods();
		
		skipMethods();
		
		skipMethods();
		skipMethods();
								
		skipMethods();
		
		skipMethods();
		
		skipMethods();
		skipMethods();
		
		skipMethods();
		
		skipMethods();					
		skipMethods();
		return min;
	}
	
	
	/*
	 * test three choice: pure, single traversal, multiple-traversal
	 */
	void testIndex(int[] best_config, double min) throws IOException {
	//	int signBall[] = {0,1,1,0,0,1,0,0,1,0,1,1,0};		
		System.out.println("Index-PAMI vs."+ min);
	//	for(int i=0; i<best_config.length; i++)
	//		System.out.println(best_config[i]);
		int signBall[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall);
		staticKmeans(true, false, true);//15 the index method that traverse the index every time, but without any bound computation.				
	//	if(time[counter-1] < 2*min) {// how to reduce the costs
			setSign(best_config);
			if(best_config[4]==1)
				System.out.println("cccc");
			prunenode = 0; 
			System.out.println("Index-single");// cascading with a fixed configration, not sure whether, it can beats most, which verify that our method
			staticKmeans(true, true, false);// the unified framework we propose, traversing the tree in single time, and use the basic bound
		
			prunenode = 0;// this is cascading method, usd for 
			System.out.println("Index-multiple");//if the pruning in the first iteration is good, we can set indexPAMI02 as true, traverse the tree from start.
			staticKmeans(true, true, true);// the unified framework we propose, traversing the tree in every iteration
	//	}else {
	//		skipMethods();
	//		skipMethods();
	//	}
	}
	
	/* 
	 * we put all the memory-efficient methods, for our algorithm pick-means submitted to ICML 2021
	 */
	void test_PICKmeans() throws IOException{
		pckmeans = false;
		pckmeansbound = false;
		pckmeanspointboundRecursive = false;
		maxIteration = -1;//just run single time of iteration, to avoid to be too slow.
		int sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("====Lloyd");			//2			
		setSign(sign);
		skipMethods();
	//	staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		maxIteration = MAXITE;
		
		System.out.println("Sequential-elkan");//3
		int signelkan[] = {0,0,1,0,0,0,0,0,0,0,0,0,0};
		setSign(signelkan);
		skipMethods();
	//	staticKmeans(false, true, false); 
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		
		System.out.println("====Sequential-Hamerly");// hamerly
		int signharmly[] = {0,0,0,1,0,0,0,0,0,0,0,0,0};
		setSign(signharmly);
		skipMethods();
	//	staticKmeans(false, true, false);
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
																		
		System.out.println("====Drake12-SortCenter");//5 drake
		int signdrake[] = {0,0,0,0,1,0,0,0,0,0,0,0,0};
		setSign(signdrake);
		skipMethods();
		binSortcenter = (int)((float)k/4);
	//	staticKmeans(false, true, false);//test the drake12
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		
		System.out.println("====Heap-15");// heap 8
		int signheap[] = {0,0,0,0,0,0,0,1,0,0,0,0,0};
		setSign(signheap);
		skipMethods();
	//	staticKmeans(false, true, false);//test the heap
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		
		pckmeanspointboundRecursive = false;
		System.out.println("====Sequential-Yinyang");//yinyang
		int signYinyang[] = {0,0,1,0,0,0,0,0,1,0,0,0,0};
		setSign(signYinyang);
		skipMethods();
	//	staticKmeans(false, true, false);
		memoryUsage[counter] = getAllMemory(k/10, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		
		maxIteration = MAXITE;
		int pami_sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("====PAMI20");//pami 20 method
		pami20 = true;
		setSign(pami_sign);
		skipMethods();
	//	staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		pami20 = false;
		
		
		maxIteration = MAXITE;
		
		System.out.println("====Index-PAMI: scan every centroid");
		int signBall[] = { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		setSign(signBall);
		maxIteration = -1;
		skipMethods();
	//	staticKmeans(true, false, true);
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		maxIteration = MAXITE;
		
		
		System.out.println("====pick-means with no bound on knn and using queue traversal");// without using our accelerating tricks
		pckmeans = true;
		int signBall1[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall1);
		maxIteration = -1;
		skipMethods();
	//	staticKmeans(true, false, true);
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		maxIteration = MAXITE;
		
		System.out.println("====pick-means with bound on knn and using queue traversal");// without using our accelerating tricks
		pckmeans = true;
		pckmeansbound = true;
	//	int signBall1[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall1);
		skipMethods();
	//	staticKmeans(true, false, true);
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		
		System.out.println("====pick-means with no interbound but with recursive traversal");// the most efficent one
		pckmeans = true;
		pckmeansbound = true;
		pckmeanspointboundRecursive = true;
		pckmeansUsinginterbound = false;//not using the inter bound
		int signBall2[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall2);
		skipMethods();
	//	staticKmeans(true, false, true);
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		pckmeans = false;
		
		System.out.println("====pick-means with recursive traversal");// the most efficient one we propose
		maxIteration = 100;
		pckmeans = true;
		pckmeansbound = true;
		pckmeanspointboundRecursive = true;
		pckmeansUsinginterbound = true;
		setSign(signBall2);
		staticKmeans(true, false, true);
		memoryUsage[counter] = getAllMemory(0, 2);
		System.out.println("memory cost:"+memoryUsage[counter]+",");
		pckmeans = false;
	}
	
	
	void testPickmeansVariousSetting() throws IOException {
		int pami_sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("PAMI20");			//2	
		pami20 = true;
		setSign(pami_sign);
	//	staticKmeans(false, false, false);// the lloyd's algorithm, the standard method
		pami20 = false;
		
		//test lloyd algorithm
		double min = Double.MAX_VALUE;
		int sign[] = {1,0,0,0,0,0,0,0,0,0,0,0,0};
		System.out.println("Lloyd");			//2			
		setSign(sign);
	//	staticKmeans(false, false, false);// the lloyd's algorithm, the standard method*/
		
		pckmeans = false;
		int signBall11[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall11);
	//	staticKmeans(true, false, true);
		
		//****new method here ********/
		pckmeans = true;
		int signBall[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall);
	//	staticKmeans(true, false, true);
		pckmeans = false;
		
		//test the pckmeans by setting pckmeans as true, with nn method
		pckmeans = true;// use knn without ub
		pckmeansbound = true;// use knn with ub and queue
		pckmeanspointboundRecursive = true; // use ub and recursive traversal
		int signBall1[] = {0,1,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBall1);
		staticKmeans(true, false, true);
		pckmeans = false;
		
		//test the pckmeans by setting pckmeans as true, with dual tree method
	/*	pckmeans = true;
		int signBalln[] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
		setSign(signBalln);
		staticKmeans(false, false, false);
		pckmeans = false;*/
		// index data points only
	//	staticKmeans(true, false, true);
		//clean the index and assigned arraylist
	}
	
	/*
	 * get the optimal configuration
	 */
	double getOptimalConfig(double min, int[] sign, int[] opt) {
		if (time[counter-1] < min && sign[4] == 0 && sign[7] == 0 && sign[6]==0) {// index cannot use with 4, 6, 7 together
			min = time[counter-1];
			for(int i=0; i<sign.length; i++) {
				opt[i] = sign[i];
			}
		}
		return min;
	}

	
	/*
	 * rank all the method by their running time
	 */
	String generateRank(double time[], int number) {
		HashMap<Integer, Double> orderTime = new HashMap<Integer, Double>();
		for(int ai=0; ai<number; ai++) {//output the rank
			orderTime.put(ai, time[ai]);
		}
		LinkedHashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
		orderTime.entrySet().stream().sorted(Map.Entry.comparingByValue())
		    .forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue()));
		String aString = "";
		for(int a: sortedMap.keySet())
			aString += a+",";
		return aString;
	}
	
	/*
	 * discover the untested combination, future work
	 */
	double testUndiscover(int[] best_sign, double min) throws IOException {
		System.out.println("Sequential-full optimizations without annu");
		int newsign1[] = {0,0,1,0,0,0,0,0,1,1,1,1,1};// without annu
		setSign(newsign1);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, newsign1, best_sign);
			
		System.out.println("Sequential-full optimizations without regroup");
		int newsign2[] = {0,0,1,0,0,1,0,0,1,1,1,0,1};//without regroup
		setSign(newsign2);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, newsign2, best_sign);
		
		System.out.println("Sequential-full optimizations without block");
		int newsign3[] = {0,0,1,0,0,1,0,0,1,1,0,1,1};//without block
		setSign(newsign3);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, newsign3, best_sign);
			
		System.out.println("Sequential-full optimizations without expon");
		int newsign4[] = {0,0,1,0,0,1,0,0,1,0,1,1,1};//without exp
		setSign(newsign4);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, newsign4, best_sign);
		
		System.out.println("Sequential-full optimizations without yinyang");
		int newsign5[] = {0,0,1,0,0,1,0,0,0,1,1,1,1};//without yinyang
		setSign(newsign5);
		staticKmeans(false, true, false);
		min = getOptimalConfig(min, newsign5, best_sign);
		return min;
	}
	
	/*** lightweight k-means for large-scale point cloud ***/
	
	/*
	 * we compute the inner bound using self-join operation, to avoid compute k*k distance, this is a lot when k is large
	 */
	void innerBound(indexNode centroidNode, double [][]centroidMatrix) {
		// it share the same framework with dual-tree, we can call our join method to do this.
	}

	
	// for all the bound, we estimate the quantilized bound, which can save much space if using one byte.
	void initializationPICKmeans() {
		assigned = new short[trajectoryNumber];
		newassigned = new short[trajectoryNumber];
		quantilizedLowerBound = new double[trajectoryNumber];
		quantilizedUpperBound = new double[trajectoryNumber];
		for(int i=0; i<trajectoryNumber; i++)
			quantilizedUpperBound[i] = Double.MAX_VALUE;
		pointcounterPruned = new short[trajectoryNumber];
		maxdrift = new double[k];
		unitQuantization = 0.01;// this needs auto configuration
		for(int i=0; i<k; i++) {
			CENTERSEuc.get(i).clearSet();
		}
	}
	
	// traverse the datapoint index, and get the assginment
	public void recursivePICKmeans(indexNode dataRoot, indexNode centroidNode, double [][]dataMatrix, 
			double [][]centroidMatrix, double ubfather, int centroid) {
		storeCentroid();
		System.out.println("covered points: "+ centroidNode.getTotalCoveredPoints());
		if(iterationTimes==0)
			initializationPICKmeans();
		else {
			for(int assignedID = 0; assignedID<k; assignedID++)
				for(int i=0; i<k; i++)// get the maximum except the assignedID.
					if(i!= assignedID && maxdrift[assignedID] < group_drift[i])
						maxdrift[assignedID] = group_drift[i];
			updatebound();// update lb and ub by the drift
			for(int i=0; i<trajectoryNumber; i++)
				pointcounterPruned[i] = 0;
			dataRoot.cleanIndex(k, group_drift, maxdrift, unitQuantization);// update the assigned for nodes by calling the function
		}
		int counterr = 0;
		for(int i=0; i<trajectoryNumber; i++) {
			assigned[i] = newassigned[i];
			if(newassigned[i]==0)
				counterr++;
		}
		System.out.println(iterationTimes+"aaaa"+counterr);
		dualTreeKmeans(dataRoot, centroidNode, dataMatrix, centroidMatrix, unitQuantization, ubfather, centroid);
		// we use a new one to create
		updateCentroid(newassigned, dataMatrix);
		for(int i=0; i<trajectoryNumber; i++) {
			assigned[i] = newassigned[i];
		}
	}

	
	/*
	 * update the bound of all the nodes
	 */
	private void updatebound() {
		for(int i=0; i<trajectoryNumber; i++) {
		//	quantilizedUpperBound[i] = (byte)((double)quantilizedUpperBound[i]*unitQuantization + group_drift[assigned[i]-1]/unitQuantization+1);
		//	quantilizedLowerBound[i] = (byte)((double)quantilizedUpperBound[i]*unitQuantization - maxdrift[assigned[i]-1]/unitQuantization-1);
			quantilizedUpperBound[i] += group_drift[assigned[i]-1];
			quantilizedLowerBound[i] -= maxdrift[assigned[i]-1];
		}
	}
	
	public short[] setAssign(indexNode datarNode, short nearest) {
		// we can also update the upper bound from lb
		datarNode.setNearestCluster(nearest);
	//	System.out.println(nearest);
		if(!datarNode.isLeaf()) {
			for(indexNode child: datarNode.getNodelist()) {
				newassigned = setAssign(child, nearest);
			}
		}else {
			for(int pointid: datarNode.getpointIdList())
				newassigned[pointid-1] = nearest;
		}
		return newassigned;
	}
	
	public short[] setAllPrune(indexNode datarNode, short prundNum) {
		datarNode.setPrunedCounterSingle(prundNum);
		if(!datarNode.isLeaf()) {
			for(indexNode child: datarNode.getNodelist()) {
				pointcounterPruned = setAllPrune(child, prundNum);
			}
		}else {
			for(int pointid: datarNode.getpointIdList())
				pointcounterPruned[pointid-1]= prundNum;
		}
		return pointcounterPruned;
	}

	/*
	 * A Dual-Tree Algorithm for Fast k-means Clustering With Large k, SDM 2017, we use multiple quantilized bounds
	 * when k is small, we put all the points into the centroid node.
	 */
	void dualTreeKmeans(indexNode dataRoot, indexNode centroidNode, double [][]dataMatrix, 
			double [][]centroidMatrix, double unit, double ubfather, int centroid) {
		int scannedCentrorid = dataRoot.getPrunedCounter();
		if(scannedCentrorid == k) {
			System.out.println("aaaaaaaaaa");
			return;
		}
		short assignedID = dataRoot.getAssignedCluster(); // check which cluster it assigned
		double ub = dataRoot.getUpperBound(unit); // it is maximum as initialized or if not pruned
		double lb = dataRoot.getLowerBound(unit);
		lb=0;// needs investigation, and the ub from father
		if(assignedID > 0 && ub < Math.max(lb, interMinimumCentoridDis[assignedID-1]/2.0)){// maintain in the node, this is for second and following iterations
			newassigned = setAssign(dataRoot, assignedID);
			pointcounterPruned = setAllPrune(dataRoot,(short)k);// set the getPrunedCounter as k
			return;
		}
		double centroidData[] = null;
		double radius = 0;
		short pruneNum = 1;
		if(centroidNode!=null) {
			centroidData = centroidNode.getPivot();
			radius = centroidNode.getRadius();
			pruneNum = (short)centroidNode.getTotalCoveredPoints();
		}else {
			centroidData = centroidMatrix[centroid-1];
		}		 
		double pivotDis = Util.EuclideanDis(dataRoot.getPivot(), centroidData, dimension);
		double newub = pivotDis + dataRoot.getRadius() + radius; // we will change to half ball theory, which is tighter
		if(newub < ub) {//update the bound as we found a tighter bound,
			if(centroidNode==null) {// the centroid is a point
				newassigned = setAssign(dataRoot, (short)centroid); //assign the node
			}
			ub = newub;//Math.min(newub,ubfather); // cannot rank, but priority queue can rank, see which one is good, we investigate both
			quantilizedUpperBound = dataRoot.setUpperBound(ub, unit, quantilizedUpperBound);
		}
		double nodelb = pivotDis - dataRoot.getRadius() - radius;// how to update the bound
		if(nodelb > ub) {// if the lower bound is bigger, we do not need to go deeper of this pair,
			dataRoot.setPrunedCounter(pruneNum, pointcounterPruned);// and all its children, pruned
			short nearest = dataRoot.getnearestID();
			if(scannedCentrorid + pruneNum == k && nearest > 0) {//if it equals k, we assign the node directly, we add the node to assign id
				newassigned = setAssign(dataRoot, nearest);
				pointcounterPruned = setAllPrune(dataRoot, (short)k);
				if(nearest != assignedID) {// assign the node to the sum if not include in the current node, 
					if(assignedID>0)
						CENTERSEuc.get(assignedID-1).removeNode(dataRoot);
					CENTERSEuc.get(nearest-1).addNode(dataRoot);
				}
				System.out.println("bbbbb");
				return;
			}
		} else {// split the centroid nodes and data nodes, and access its child node
			if (dataRoot.isLeaf() && centroidNode == null) {// we will not go deeper as we store our data here
				for (int i : dataRoot.getpointIdList()) {
					short previousNearest = assigned[i-1];
					double currentDis = quantilizedUpperBound[i-1]; // quantilization is not precise
					double second_min_dist = quantilizedLowerBound[i-1];// we store the lb
					if(pointcounterPruned[i-1] == k || (previousNearest>0 && currentDis < Math.max(0, interMinimumCentoridDis[previousNearest-1]/2.0))){// check lb and ub to see whether we can assign directly
						newassigned[i-1] = (short)previousNearest;
						pointcounterPruned[i-1] = (short)k;
						System.out.println("dddd");
						continue;
					}
					double dataPoint[] = dataMatrix[i - 1];
					double dis = Util.EuclideanDis(dataPoint, centroidData, dimension);// pair wise computation and pruning by checking the point bound, and update the node id
					if (dis < currentDis) {
						second_min_dist = currentDis;
						currentDis = dis;
						newassigned[i-1] = (short) centroid;
						setUpperBound(currentDis, unit, i-1);
					} else if (dis < second_min_dist && dis != currentDis) {// store the second nearest, for next iteration pruning and node pruning
						second_min_dist = dis;
					}
					setLowerBound(second_min_dist, unit, i-1);// store the second nearest
					pointcounterPruned[i-1]++;
					if (pointcounterPruned[i-1] >= k) {// all the centroids have been scanned for this point
						if (newassigned[i-1] != previousNearest) {
							if(previousNearest > 0) {
								CENTERSEuc.get(previousNearest-1).deleteSinglePointToCluster(i-1);
								CENTERSEuc.get(previousNearest-1).minusSum(dataPoint);// for incremental use
							}
							CENTERSEuc.get(newassigned[i-1]-1).addPointToCluster(i-1);
							CENTERSEuc.get(newassigned[i-1]-1).addSum(dataPoint);
						}
						System.out.println("cccccc");
					}
				}
			} else {
				if (centroidNode != null) {
					if (centroidNode.isLeaf())
						for (int childpoint : centroidNode.getpointIdList())
							dualTreeKmeans(dataRoot, null, dataMatrix, centroidMatrix, unit, pivotDis+centroidNode.getRadius(), childpoint);
					else if (dataRoot.isLeaf())
						for (indexNode childNode : centroidNode.getNodelist())
							dualTreeKmeans(dataRoot, childNode, dataMatrix, centroidMatrix, unit, pivotDis+centroidNode.getRadius(), 0);
					else
						for (indexNode childNode : dataRoot.getNodelist())
							for (indexNode childNode1 : centroidNode.getNodelist())
								dualTreeKmeans(childNode, childNode1, dataMatrix, centroidMatrix, unit, pivotDis+childNode.getDisFather(), 0);//sort by the distance
				} else
					for (indexNode childNode : dataRoot.getNodelist())
						dualTreeKmeans(childNode, null, dataMatrix, centroidMatrix, unit, pivotDis + childNode.getDisFather(), centroid);
			}
		}
	}
	
	private void setUpperBound(double currentDis, double unit, int i) {
	//	quantilizedUpperBound[i] = (byte)(currentDis/unit+1);
		quantilizedUpperBound[i]= currentDis;
	}

	private void setLowerBound(double second_min_dist, double unit, int i) {
	//	quantilizedLowerBound[i] = (byte)(second_min_dist/unit-1);
		quantilizedLowerBound[i] = second_min_dist;
	}
	
	
	// integrate the lb and father distance bound, and
	double updateCentroid(short assingned[], double[][] datamapDouble) {
		double meanCentroid[][] = new double[k][];
		int number[] = new int[k];
		for(int i=0; i<k; i++) {
			meanCentroid[i] = new double[dimension];
		}
		for(int i=0; i<trajectoryNumber; i++) {
			int clusterID = assingned[i];	
			number[clusterID-1]++;
			for(int j=0; j<dimension; j++) {
				meanCentroid[clusterID-1][j] += datamapDouble[i][j];
			}
		}
		int counterUnassigned = 0;
		for(int j=0; j<k; j++) 
			if(number[j]==0)
				counterUnassigned++;
		
		System.out.println("there are multiple clusters with no points:"+counterUnassigned);
		double drfitSum = 0;
		for(int ii=0; ii<k; ii++) {
			for (int i = 0; i < dimension; i++) {
				meanCentroid[ii][i] /= number[ii];//compute the new
			}
			double drfit = Util.EuclideanDis(meanCentroid[ii], centroidsData[ii], dimension);
			for (int i = 0; i < dimension; i++) {
			//	centroidsData[ii][i] = meanCentroid[ii][i];
			//	System.out.println(meanCentroid[ii][i]);
			}
			drfitSum += drfit;
			center_drift.put(ii, drfit);				
			if(maxDrift < drfit)
				maxDrift = drfit;
			group_drift[ii] = drfit;
		}
		return drfitSum;
	}	
}
