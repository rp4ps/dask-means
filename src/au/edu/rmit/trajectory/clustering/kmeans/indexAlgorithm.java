package au.edu.rmit.trajectory.clustering.kmeans;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kpaths.Util;
import edu.wlu.cs.levy.cg.KDTree;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;
import es.saulvargas.balltrees.BallTreeMatrix;
import scala.reflect.internal.Trees.Return;
import skyline0623.balltree.BallTree;
import skyline0623.balltree.Hypersphere;
import skyline0623.balltree.Point;
import skyline0623.balltree.Process;
import covertree.*;
import java.util.*;

public class indexAlgorithm<E> {

	int distanceCompute = 0;
	int NodeAccess = 0;
	int dataAccess = 0;
	int globalNodeid = 1;
	
	
	public indexAlgorithm() {
		// TODO Auto-generated constructor stub
	}
	
	public void buildKDtree(int dims, double[][] itemMatrix) throws KeySizeException, KeyDuplicateException {
		KDTree<Long> kt = new KDTree<Long>(dims);
		long idx = 1;
		for(double[] point: itemMatrix) {
			kt.insert(point, idx++);//fast construction: point, value
		//	System.out.println("inset kd"+idx);
		}		
		indexNode rootKmeans = new indexNode(dims);
	//	kt.traverseConvert(rootKmeans, kt.getroot(), dims);	// traversing the index is hard for kd-tree
	}
	
	public void setGloabalid() {
		globalNodeid = 1;
	}
	/*
	 * get the count of the tree.
	 */
	public double[] updateSum(indexNode root, int dimension) {
		if(root.isLeaf()) {	
			return root.getSum();
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
		//	System.out.println(listnode.size());
			double []sum = new double[dimension];
			for(indexNode aIndexNode: listnode) {
				double []sumd = updateSum(aIndexNode, dimension);
				for(int i=0; i<dimension; i++)
					sum[i] += sumd[i];
			}
			root.setSum(sum);
			return sum;
		}
	}
	
	/*
	 * update the sum vector of index with fair constraint
	 */
	public double[] updateSumFair(indexNode root, int dimension, int userID[], Map<Integer, Integer> userNumber, double[][] datamapDouble) {
		if(root.isLeaf()) {	
			Set<Integer> listnode = root.getpointIdList();
			double[] sum = new double[dimension];
			for (int aIndexNode : listnode) {
				double[] sumd = datamapDouble[aIndexNode - 1];//this needs to be fixed
			//	System.out.println(aIndexNode + "," +userID.length);
				for (int i = 0; i < dimension; i++) {
					sum[i] += sumd[i] * 1.0/userNumber.get(userID[aIndexNode - 1]);
				}
			}
			root.setSum(sum);
			return sum;
		}else {
			Set<indexNode> listnode = root.getNodelist();
			double []sum = new double[dimension];
			for(indexNode aIndexNode: listnode) {
				double []sumd = updateSumFair(aIndexNode, dimension, userID, userNumber, datamapDouble);
				for(int i=0; i<dimension; i++)
					sum[i] += sumd[i];
			}
			root.setSum(sum);
			return sum;
		}
	}
	/*
	 * update the number of points for fair clustering
	 */
	public double updateCoveredPointsFair(indexNode root, int dimension, int userID[], Map<Integer, Integer> userNumber, double[][] datamapDouble) {
		if(root.isLeaf()) {	
			Set<Integer> listnode = root.getpointIdList();
			double coveredpints = 0;
			for (int aIndexNode: listnode) {
				coveredpints += 1.0/userNumber.get(userID[aIndexNode - 1]);
			}
			root.setTotalCoveredPointsFair(coveredpints);
			return coveredpints;
		}else {
			Set<indexNode> listnode = root.getNodelist();
			double coveredpints = 0;
			for(indexNode aIndexNode: listnode) {
				coveredpints += updateCoveredPointsFair(aIndexNode, dimension, userID, userNumber, datamapDouble);
			}
			root.setTotalCoveredPointsFair(coveredpints);
			return coveredpints;
		}
	}
	
	/*
	 * get the count of the tree.
	 */
	public void updateNodeId(indexNode root, int nodeid) {
		root.setNodeid(nodeid++);
		if(!root.isLeaf()){
			Set<indexNode> nodelist = root.getNodelist();
			for(indexNode aIndexNode: nodelist) {
				updateNodeId(aIndexNode, nodeid++);
			}
		}
	}
	
	/*
	 * get the count of the tree.
	 */
	public void updateDatasetNodeId(indexNode root, int nodeid) {
		root.setNodeid(nodeid++);
		if(!root.isrootLeaf()){
			Set<indexNode> nodelist = root.getNodelist();
			for(indexNode aIndexNode: nodelist) {
				updateDatasetNodeId(aIndexNode, nodeid++);
			}
		}
	}
	
	public indexNode buildMtree(double[][] itemMatrix, int dimension, int capacity) {// too slow
		System.out.println("Building M-tree...");
		long startTime1 = System.nanoTime();
		PointMtree mindex = new PointMtree(capacity);		//capacity
		int idx = 1;
		for(double[] point: itemMatrix) {
			mindex.buildMtree(point, idx++);//create the M-tree
		//	System.out.println("inset M-tree"+idx);
		}
		indexNode rootKmeans = new indexNode(dimension);
		mindex.traverseConvert(rootKmeans, mindex.getRoot(), dimension);		// traversing the index
		updateSum(rootKmeans, dimension);
		long endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0);
	//	System.out.println("the count of M-tree is " + rootKmeans.getTotalCoveredPoints());
		return rootKmeans;
	}
	
	
	public indexNode buildCoverTree(int dims, double[][] datamapEuc) {
		//build the tree based on 
		CoverTree<Integer> cTree = new CoverTree<Integer>();
		int idx=0;
		for(double[] point:  datamapEuc) {
			cTree.insert(idx, point);//fast construction: point, value
		//	System.out.println("inset cover tree"+idx);
			idx++;
		}
		indexNode rootKmeans = new indexNode(dims);
		cTree.traverseConvert(rootKmeans, dims);
		updateSum(rootKmeans, dims);
		return rootKmeans;
		
	}
	
	public indexNode buildBalltree(Map<Integer, double[]> datamapEuc, int dimension, int capacity) {// too slow	
	//	System.out.println("Building Ball-tree...");
		long startTime1 = System.nanoTime();	
		for(int idx: datamapEuc.keySet()) {
			double[] point = datamapEuc.get(idx);
			Process.DIMENSION = dimension;
			Process.INSTANCES.add(new Point(point));
			Process.MAX_INSTANCE_NUM_NOT_SPLIT = capacity;
		//	System.out.println("inset Ball-tree"+idx);
		}
		Hypersphere BALL_TREE = BallTree.buildAnInstance(null);		
		indexNode rootKmeans = new indexNode(dimension);
		BALL_TREE.traverseConvert(rootKmeans, dimension);
	//	computeFartherToChild(rootKmeans);
		long endtime = System.nanoTime();
		System.out.print((endtime-startTime1)/1000000000.0+",");
		System.out.println("the count of Ball-tree is " + 
				rootKmeans.getTotalCoveredPoints()+", the radius is "+rootKmeans.getRadius());
	//	System.out.println("the count of Ball-tree is " + rootKmeans.getTotalCoveredPoints());
		return rootKmeans;
	}
	
	public indexNode buildBalltree2(double[][] itemMatrix, int dimension, int capacity, int userID[], Map<Integer, Integer> userNumber) {// too slow	
		System.out.println("Building Ball-tree using Matrix...");
		long startTime1 = System.nanoTime();
		int deepth = (int) (Math.log(itemMatrix.length)/Math.log(2));//the deepth is computed based on binary tree
		indexNode rootKmeans = BallTreeMatrix.create(itemMatrix, capacity, deepth);//we should not set the deepth too deep
		if(userNumber==null)
			updateSum(rootKmeans, dimension);
		else {
			updateSumFair(rootKmeans, dimension, userID, userNumber, itemMatrix);// update sum fair for useage.
			double a = updateCoveredPointsFair(rootKmeans, dimension, userID, userNumber, itemMatrix);
		}
		updateNodeId(rootKmeans, 1);
		calculateMaxBoundBox(rootKmeans, dimension, itemMatrix);
		calculateMinBoundBox(rootKmeans, dimension, itemMatrix);
		long endtime = System.nanoTime();
		System.out.println("index time cost: "+(endtime-startTime1)/1000000000.0);
		System.out.println("the count of Ball-tree using Matrix is " + 
			rootKmeans.getTotalCoveredPoints()+", the radius is "+rootKmeans.getRadius());
		return rootKmeans;
	}
	
	
	/*
	 * get all the points under a node
	 */
	public ArrayList<Integer> getPointsIDlist(indexNode root) {
		ArrayList<Integer> result = new ArrayList<>();
		if(root == null)
			return null;
		if(root.isLeaf()) {	
			result.addAll(root.getpointIdList());
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				result.addAll(getPointsIDlist(aIndexNode));
			}
		}
		return result;
	}
	
	public void setdistanceCompute(int distanceCompute) {
		this.distanceCompute = distanceCompute;
	}
	
	public int getdistanceCompute() {
		return distanceCompute;
	}
	
	public void setNodeAccess(int nodeAccess) {
		this.NodeAccess = nodeAccess;
	}
	
	public int getNodeAccess() {
		return NodeAccess;
	}
	
	
	public void setdataAccess(int dataAccess) {
		this.dataAccess = dataAccess;
	}
	
	public int getdataAccess() {
		return dataAccess;
	}
	
	/*
	 * search all the points within a distance radius, known as similarity search, using M-tree may be better.
	 */
	public ArrayList<Integer> SimilaritySearchBall(double radius, double point[], indexNode root, int dimension, 
			double[][] itemMatrix) { 
		ArrayList<Integer> result = new ArrayList<>();
		if (root.isLeaf()) {
			for (int id : root.getpointIdList()) {
				dataAccess++;
				double distance = Util.EuclideanDis(itemMatrix[id - 1], point, dimension);
				distanceCompute++;
				if (distance < radius)
					result.add(id);
			}
		}else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				NodeAccess++;
				double distance = Util.EuclideanDis(aIndexNode.getPivot(), point, dimension);
				distanceCompute++;
				if(radius - distance >= aIndexNode.getRadius()) {
				//	System.out.println("aaaaaaaaaaaaaaaaaa");
					result.addAll(getPointsIDlist(aIndexNode));
				}else {// cannot be pruned
					result.addAll(SimilaritySearchBall(radius, point, aIndexNode, dimension, itemMatrix));
				}
			}
		}
		return result;
	}
	
	/*
	 * search the nearest neighbor, used for point
	 */
	public void NearestNeighborSearchBall(double point[], indexNode root, int dimension, 
			double[][] itemMatrix, double []minDistnearestID) { 
		if (root.isLeaf()) {
			for (int id : root.getpointIdList()) {
				double distance = Util.EuclideanDis(itemMatrix[id - 1], point, dimension);
				distanceCompute++;
				if (distance <= minDistnearestID[0]) {
					minDistnearestID[0] = distance;
					minDistnearestID[1] = (double)id;
				}
			}
		}else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				double distance = Util.EuclideanDis(aIndexNode.getPivot(), point, dimension);
				distanceCompute++;
				if(minDistnearestID[0] >= distance - aIndexNode.getRadius()) {
					NearestNeighborSearchBall(point, aIndexNode, dimension, itemMatrix, minDistnearestID);
				}
			}
		}
	}
	
	// search two nearest neighbors, minDistnearestID stores two points
	public void TwoNearestNeighborSearchBall(double point[], indexNode root, int dimension, 
			double[][] itemMatrix, double []minDistnearestID) {
		if (root.isLeaf()) {
			for (int id : root.getpointIdList()) {
				double distance = Util.EuclideanDis(itemMatrix[id - 1], point, dimension);
				distanceCompute++;
				if (distance <= minDistnearestID[0]) {
					minDistnearestID[0] = distance;
					minDistnearestID[1] = (double)id;
					if(minDistnearestID[0] < minDistnearestID[2]) {
						double temp = minDistnearestID[2];
						minDistnearestID[2] = minDistnearestID[0];
						minDistnearestID[0] = temp;
						double temp1 = minDistnearestID[3];
						minDistnearestID[3] = minDistnearestID[1];
						minDistnearestID[1] = temp1;
					}
				}
			}
		}else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				double distance = Util.EuclideanDis(aIndexNode.getPivot(), point, dimension);
				distanceCompute++;
				if(minDistnearestID[0] > (distance - aIndexNode.getRadius())) {
					TwoNearestNeighborSearchBall(point, aIndexNode, dimension, itemMatrix, minDistnearestID);
				}
			}
		}
	}
	
	/*
	 * get the count of points in the tree.
	 */
	public int getcount(indexNode root) {
		if(root == null)
			return 0;
		if(root.isLeaf()) {	
			return root.getpointIdList().size();
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				max += getcount(aIndexNode);
			}
			return max;
		}
	}
	
	/*
	 * get the count of points in the tree.
	 */
	public int getRootCount(indexNode root) {
		if(root == null)
			return 0;
		if(root.isrootLeaf()) {
			return 1;
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				max += getRootCount(aIndexNode);
			}
			return max;
		}
	}
	
	/*
	 * get the highest weight of the tree
	 */
	public int getHeight(indexNode root) {
		if(root.isLeaf())
			return 0;
		else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				if(getHeight(aIndexNode)>max) {
					max = getHeight(aIndexNode);
				}
			}
			return max+1;
		}
	}
	
	
	/*
	 * get the radius of the leaf node tree and leaf number, divided by the maximum radius, used for dataset search
	 */
	public int getLeafRadius(indexNode root, ArrayList<Double> radius, double maximum, ArrayList<Double> distanceToFather, ArrayList<Double> numPoint, ArrayList<Double> nodedept, int depth) {
		if(root.isLeaf()) {
			radius.add(root.getRadius()/maximum);
			distanceToFather.add(root.getDisFather()/maximum);
			numPoint.add((double)(root.getTotalCoveredPoints()));
			nodedept.add((double)depth);
			return 1;
		}else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				max += getLeafRadius(aIndexNode, radius, maximum,distanceToFather,numPoint, nodedept, depth+1);
			}
			return max;
		}
	}
	
	
	/*
	 * store the tree as a directed graph, and each node information 
	 */
	public void storeFeatures(indexNode root, int nodeid, int depth, String edge, String node, int datasetID) {
		// assign an id to the node
		String contentString = Integer.toString(datasetID)+","+Integer.toString(nodeid)+","+Integer.toString(root.getTotalCoveredPoints())+","+Integer.toString(depth)+","+Double.toString(root.getRadius());
		if(root.isLeaf()) {
			//store the leaf node id
			contentString += ",0";// leaf node
		}else {
			Set<indexNode> listnode = root.getNodelist();
			contentString += ",1";//internal node
			for(indexNode aIndexNode: listnode) {
				int assignedNodeID = ++globalNodeid;
				double distance = Util.EuclideanDis(root.getPivot(), aIndexNode.getPivot(), root.getPivot().length);
				aIndexNode.setdistanceToFarther(distance);
			//	System.out.println(aIndexNode.getDisFather());
				// write the edge information into files
				Util.write(edge, Integer.toString(datasetID)+","+Integer.toString(nodeid)+","+Integer.toString(assignedNodeID)+","+Double.toString(aIndexNode.getDisFather())+"\n");
				storeFeatures(aIndexNode, assignedNodeID, depth+1, edge, node, datasetID);
			}
		}
		double[] pivot = root.getPivot();
		for(double a:pivot)
			contentString += ","+Double.toString(a);
		Util.write(node, contentString+"\n");
	}
	
	
	/*
	 * store the tree as a directed graph, and each node information,
	 * data lake index is also supported
	 */
	public void storeIndex(indexNode root, int nodeid, String node, double fartherdis) {
		// assign an id to the node
		String contentString = Integer.toString(nodeid)+","+Integer.toString(root.getTotalCoveredPoints())+","
		+Double.toString(fartherdis)+","+Double.toString(root.getRadius());
		if(root.isLeaf()) {
			contentString += ",0:";// leaf node
			for(int pointid: root.getpointIdList())
				contentString += Integer.toString(pointid)+",";
		}else {
			Set<indexNode> listnode = root.getNodelist();
			contentString += ",1:";//internal node
			for(indexNode aIndexNode: listnode) {
				int assignedNodeID = ++globalNodeid;
				contentString += Integer.toString(assignedNodeID)+",";
				double distance = Util.EuclideanDis(root.getPivot(), aIndexNode.getPivot(), root.getPivot().length);
				aIndexNode.setdistanceToFarther(distance);
				storeIndex(aIndexNode, assignedNodeID, node, distance);
			}
		}
		contentString += ":";// store the pivot
		double[] pivot = root.getPivot();
		for(double a:pivot)
			contentString += Double.toString(a)+",";
		contentString += ":";// store the max//add the max and min bounding box
		for(double max: root.getMBRmax())
			contentString += Double.toString(max)+",";
		for(double min: root.getMBRmin())
			contentString += Double.toString(min)+",";
		Util.write(node, contentString+"\n");
	}
	
	/*
	 * restore the index structure from files, each dataset will have one index
	 */
	public Map<Integer, Map<Integer, indexNode>> restoreIndex(String foldername, int dimension, Map<Integer, double[][]> dataMapPorto) {
		// restore the index based on the index files, store every node in a hash and later set up
		File folder = new File(foldername);
		File[] fileNames = folder.listFiles();
		int numberIndex = fileNames.length-2;
		Map<Integer, Map<Integer, indexNode>> roots =  new HashMap<Integer, Map<Integer, indexNode>>();
        for(int ii=1; ii<numberIndex; ii++){
                try {
                	try (BufferedReader br = new BufferedReader(new FileReader(foldername+Integer.valueOf(ii)+".txt"))) {
            			String strLine;
            			Map<Integer, indexNode> nodeidIndexnodesMap = new HashMap<Integer, indexNode>();
            			while ((strLine = br.readLine()) != null) {
            				String[] splitString = strLine.split(":");
            				String[] basicInfo = splitString[0].split(",");
            				String[] liString = splitString[1].split(",");
            				String[] centerString = splitString[2].split(",");
            				String[] mbr = splitString[3].split(",");
            				indexNode aIndexNode = new indexNode(dimension);
            				aIndexNode.addEveryThing(centerString, basicInfo, liString, mbr, dimension);//
            				int nodeid = Integer.valueOf(basicInfo[0]);
            				nodeidIndexnodesMap.put(nodeid, aIndexNode);
            			}
            			roots.put(ii, nodeidIndexnodesMap);//insert the root nodes,
            		}
                } catch (IOException e) {
                    e.printStackTrace();
                }
        }
        return roots;
	}
	
	/*
	 * restore the index structure from files, each dataset will have one index, instead of reading all index into memeory like above
	 */
	public Map<Integer, indexNode> restoreSingleIndex(String foldername, int ii, int dimension) {
		// restore the index based on the index files, store every node in a hash and later set up
		Map<Integer, indexNode> nodeidIndexnodesMap = new HashMap<Integer, indexNode>();        
		try {
			try (BufferedReader br = new BufferedReader(new FileReader(foldername + Integer.valueOf(ii) + ".txt"))) {
				String strLine;
				while ((strLine = br.readLine()) != null) {
					String[] splitString = strLine.split(":");
					String[] basicInfo = splitString[0].split(",");
					String[] liString = splitString[1].split(",");
					String[] centerString = splitString[2].split(",");
					String[] mbr = splitString[3].split(",");
					indexNode aIndexNode = new indexNode(dimension);
					//add the bounding box
					aIndexNode.addEveryThing(centerString, basicInfo, liString, mbr, dimension);//
					int nodeid = Integer.valueOf(basicInfo[0]);
					nodeidIndexnodesMap.put(nodeid, aIndexNode);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
        return nodeidIndexnodesMap;
	}
	
	/*
	 * store the data lake index
	 */
	public void storeDatalakeIndex(indexNode root, int nodeid, String node, double fartherdis) {
		// assign an id to the node
		String contentString = Integer.toString(nodeid)+","+Integer.toString(root.getTotalCoveredPoints())+","
		+Double.toString(fartherdis)+","+Double.toString(root.getRadius());
		if(root.isrootLeaf()) {
			contentString += ",0:";// leaf node
			contentString += root.getDatasetID()+",";
		}else if (root.getDatasetID()<0){
			Set<indexNode> listnode = root.getNodelist();
			if(root.getDatasetID()==-1)
				contentString += ",1:";//internal node
			else
				contentString += ",2:";//internal node
			for(indexNode aIndexNode: listnode) {
				int assignedNodeID = ++globalNodeid;
				contentString += Integer.toString(assignedNodeID)+",";
				double distance = Util.EuclideanDis(root.getPivot(), aIndexNode.getPivot(), root.getPivot().length);
				aIndexNode.setdistanceToFarther(distance);
				storeDatalakeIndex(aIndexNode, assignedNodeID, node, distance);
			}
		}
		contentString += ":";// store the pivot
		double[] pivot = root.getPivot();
		for(double a:pivot)
			contentString += Double.toString(a)+",";
		contentString += ":";// store the max//add the max and min bounding box
		for(double max: root.getMBRmax())
			contentString += Double.toString(max)+",";
		for(double min: root.getMBRmin())
			contentString += Double.toString(min)+",";
		Util.write(node, contentString+"\n");
	}
	
	
	/*
	 * restore the index structure from files, we read the boundary of each dataset and create index
	 */
	public Map<Integer, indexNode> restoreDatalakeIndex(String foldername, int dimension) {
		// restore the index based on the index files, store every node in a hash and later set up
		Map<Integer, indexNode> nodeidIndexnodesMap = new HashMap<Integer, indexNode>();
		try {
			try (BufferedReader br = new BufferedReader(new FileReader(foldername))) {
				String strLine;
				while ((strLine = br.readLine()) != null) {
					String[] splitString = strLine.split(":");
					String[] basicInfo = splitString[0].split(",");
					String[] liString = splitString[1].split(",");
					String[] centerString = splitString[2].split(",");
					String[] mbr = splitString[3].split(",");
					indexNode aIndexNode = new indexNode(dimension);
					aIndexNode.addEveryThingDatalake(centerString, basicInfo, liString, mbr, dimension);//
					int nodeid = Integer.valueOf(basicInfo[0]);
					nodeidIndexnodesMap.put(nodeid, aIndexNode);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return nodeidIndexnodesMap;
	}
	
	/*
	 * indexing the dataset using a kd-tree way
	 */
	public indexNode indexDatasetKD(ArrayList<indexNode> roots, int dimension, int leafSize) {
		double minBox[] = new double[dimension];
		double maxBox[] = new double[dimension];
		for(int dim=0; dim<dimension; dim++) {
			minBox[dim] = Double.MAX_VALUE;
			maxBox[dim] = Double.MIN_VALUE;
		}
		// get the bounding box of all the nodes first,
		for(indexNode aIndexNode : roots) {
			double pivot[] = aIndexNode.getPivot();
			double radius = aIndexNode.getRadius();
			for(int dim=0; dim<dimension; dim++) {
				if(pivot[dim]+radius > maxBox[dim])
					maxBox[dim] = pivot[dim] + radius;
				if(pivot[dim]-radius < minBox[dim])
					minBox[dim] = pivot[dim] - radius;
			}
		}
		// get the pivot point, and the range in multiple dimension, and the radius
		double pivot[] = new double[dimension];
		int d=0;
		double maxrange = Double.MIN_VALUE;
		double radius = 0;
		for(int dim=0; dim<dimension; dim++) {
			pivot[dim] = (maxBox[dim] + minBox[dim])/2;
			radius += Math.pow((maxBox[dim] - minBox[dim])/2, 2); // the radius should consider d dimension
			if((maxBox[dim] - minBox[dim]) > maxrange) {
				maxrange = maxBox[dim] - minBox[dim];
				d = dim;
			}
		}
		//create a new leaf node and return
		indexNode a = new indexNode(dimension);
		a.setRadius(Math.sqrt(radius));// the radius need to be bigger.
		a.setPivot(pivot);
		a.setMBRmax(maxBox);
		a.setMBRmin(minBox);
		if(roots.size()<=leafSize) {
			for(indexNode root: roots) {
				a.addNodes(root);
			}
			a.setroot(-1);//label it as leaf node
			return a;
		}else {
			// dividing the space by the broadest dimension d
			ArrayList<indexNode> rootleft = new ArrayList<indexNode>();
			ArrayList<indexNode> rootright = new ArrayList<indexNode>();
			for(indexNode aIndexNode: roots) {
			//	System.out.println(aIndexNode.getPivot()[d]+","+pivot[d]);
				if(aIndexNode.getPivot()[d] < pivot[d]) {
					rootleft.add(aIndexNode);
				}else {
					rootright.add(aIndexNode);
				}
			}
			if(rootleft.isEmpty() || rootright.isEmpty()) {
				for(indexNode root: roots) {
					a.addNodes(root);
				}
			}else {
				a.addNodes(indexDatasetKD(rootleft, dimension, leafSize));
				a.addNodes(indexDatasetKD(rootright, dimension, leafSize));
			}
			a.setroot(-2);//lable it as internal node
			return a;
		}
	}
	
	// insert a new node into the index: finding the node that covers the ball most, 
	// and insert to the leaf node, split if exceeding the maximum number of points,
	public indexNode insertNewDataset(indexNode newNode, indexNode root, int capacity, int dimension) {
		double minBox[] = new double[dimension];
		double maxBox[] = new double[dimension];
		double pivot[] = root.getPivot();
		double pivotn[] = newNode.getPivot();
		double newpivot[] = new double[dimension];
		double maxrange = Double.MIN_VALUE;
		int d=0;
		for(int dim=0; dim<dimension; dim++) {
			minBox[dim] = Math.min(pivot[dim]-root.getRadius(), pivotn[dim]-newNode.getRadius());
			maxBox[dim] = Math.min(pivot[dim]+root.getRadius(), pivotn[dim]+newNode.getRadius());
			newpivot[dim] = (maxBox[dim] + minBox[dim])/2;
			if((maxBox[dim] - minBox[dim]) > maxrange) {
				maxrange = maxBox[dim] - minBox[dim];
				d = dim;
			}
		}
		root.setPivot(newpivot);//used the above techniques to get the new pivot and radius
		root.setRadius(maxrange/2);
		if(root.getDatasetID()==-1) {//check whether it covers root nodes directly, the lowest level of the dataset tree
			root.addNodes(newNode);// refine the node by adding this node, 
			if(root.getNodelist().size()>capacity) {//check whether we need to split, 
				indexNode nodeleftIndexNode = new indexNode(dimension);
				indexNode noderightIndexNode = new indexNode(dimension);
				int counter = 0;
				for(indexNode aIndexNode: root.getNodelist()) {
					if(aIndexNode.getPivot()[d] < pivot[d]) {// split based on the pivot point used above
						nodeleftIndexNode.addNodes(aIndexNode);
						counter++;
					}else {
						noderightIndexNode.addNodes(aIndexNode);
					}
				}
				if(counter>0 && counter<root.getNodelist().size()) {
					nodeleftIndexNode.setroot(-1);
					noderightIndexNode.setroot(-1);
					root.clearNodes();//clear the nodes
					root.addNodes(nodeleftIndexNode);
					root.addNodes(noderightIndexNode);
					root.setroot(-2);
				}
			}
		}else if(root.getDatasetID()==-2){
			indexNode tempIndexNode = null;
			double mindis = Double.MAX_VALUE;
			for(indexNode aIndexNode: root.getNodelist()) {// scan every node, and see which one is good, compute the 
				double distance = Util.EuclideanDis(aIndexNode.getPivot(), newNode.getPivot(), aIndexNode.getPivot().length);
				double coverDis = distance + aIndexNode.getRadius() + newNode.getRadius();
				if(coverDis < mindis) {
					mindis = coverDis;
					tempIndexNode = aIndexNode;
				}
			}
			root.removeNode(tempIndexNode);// we need to update this node, remove then add
			tempIndexNode = insertNewDataset(newNode, tempIndexNode, capacity, dimension);
			root.addNodes(tempIndexNode);
		}
		return root;
	}
	
	/*
	 * calculating the bounding box of each node, for bounding box 
	 */
	public double[] calculateMaxBoundBox(indexNode root, int dimension, double[][] dataMatrix) {
		root.mbrmax = new double[dimension];
		if(root.isLeaf()) {
			for(int pointid: root.getpointIdList()) {
				for(int i=0; i<dimension; i++) {
					if(dataMatrix[pointid-1][i]>root.mbrmax[i])
						root.mbrmax[i] = dataMatrix[pointid-1][i];
				}
			}
		}else {
			for(indexNode child: root.getNodelist()) {
				double []childmax = calculateMaxBoundBox(child, dimension, dataMatrix);
				for(int i=0; i<dimension; i++) {
					if(childmax[i]>root.mbrmax[i])
						root.mbrmax[i] = childmax[i];
				}
			}
		}
		return root.mbrmax;
	}
	
	/*
	 * calculating the bounding box of each node, for bounding box 
	 */
	public double[] calculateMinBoundBox(indexNode root, int dimension, double[][] dataMatrix) {
		root.mbrmin = new double[dimension];
		for(int i = 0; i<dimension; i++)
			root.mbrmin[i] = Double.MAX_VALUE;
		if(root.isLeaf()) {
			for(int pointid: root.getpointIdList()) {
				for(int i=0; i<dimension; i++) {
					if(dataMatrix[pointid-1][i]<root.mbrmin[i])
						root.mbrmin[i] = dataMatrix[pointid-1][i];
				}
			}
		}else {
			for(indexNode child: root.getNodelist()) {
				double []childmin = calculateMinBoundBox(child, dimension, dataMatrix);
				for(int i=0; i<dimension; i++) {
					if(childmin[i]<root.mbrmin[i])
						root.mbrmin[i] = childmin[i];
				}
			}
		}
		return root.mbrmin;
	}
}
