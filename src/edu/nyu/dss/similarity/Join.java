package edu.nyu.dss.similarity;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.TreeMap;
import au.edu.rmit.trajectory.clustering.kmeans.indexAlgorithm;
import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;


public class Join {
	
	// join two tables by finding nearest neighbor, can also be used in join
	public static void joinTableBaseline(double[][] query, double[][] dataset, indexNode idxQuery, indexNode idxDataset, int dimension, indexAlgorithm<Object> indexDSS) {
		int i = 1;
		Map<Integer, Integer> joinNearestMap = new TreeMap<Integer, Integer>();
		for(double []queryPoint: query) {
			//search the nearest neighbor;
			double []minDistnearestID = {Double.MAX_VALUE,0};
			indexDSS.NearestNeighborSearchBall(queryPoint, idxDataset, dimension, dataset, minDistnearestID);
			joinNearestMap.put(i++, (int)minDistnearestID[1]);
		}
		System.out.println(joinNearestMap.get(1));
	}
	
	
	/*
	 * 
	 * using the cached distance to prune, ball tree method, knn and also plus a threshold
	 * cached more information, queue and index nodes to accelerate the join
	 * implement a basic method that used ball-tree for knn, but has known the threshold, just show it is promising
	 * 
	 * The advantages are we do not need to scan the tree again, and save much resources using one pipeline, while normal work
	 * needs different methods. We compare with various work and prove the performance
	 */
	void usingIndex(PriorityQueue<queueMain> aHeaps, double Hausdorff) {
		// search one by one, stop if the lb is bigger than the haus dis
		// continue the two-level queue, compress the second queue's lower bound that is bigger than the haus dis
		// then we can cache multiple queues if space allowed, what is the space of queue?
		// once the candidate dataset is selected, user can also specify the threshold that is smaller than Hausdorff for the nearest neighbor
		PriorityQueue<queueMain> aHeapsNew = new PriorityQueue<queueMain>();
		for(queueMain qMain: aHeaps) {
			PriorityQueue<queueSecond> qSecond = qMain.getQueue();
			double ub = qMain.getbound();
			PriorityQueue<queueSecond> qSecondNew = new PriorityQueue<queueSecond>();
			for(queueSecond qnodeS: qSecond) {
				if(qnodeS.getbound()<Hausdorff) {
					qSecondNew.add(qnodeS);
				}
			}
		}
		//get a pure queue
	}
	
	/*
	 * this join will reuse the queue and it is terminated until the queue is empty, future direction
	 * 1. compress the queue by using the nearest neighbor
	 * 2. assign in batch to the nearest point, by changing the splitting option, prefer spliting the dataset nodes
	 * using the internal centroid distance
	 * 3. when can we assign a node to point directly, use our vldb paper idea, i.e., knowing each point's nearest neighbor in advance
	 * 4. the main queue is not needed, unless the hausdorff threshold is updated, we prune some pairs based on the nearest neighbor distance
	 * 
	 * 
	 * this can be used in k-means clustering, when k is large, we can also use this idea 
	 * to prune by adding more bounds between farther and child nodes, and centroid's movement in various iterations
	 * 1 using the internal centroid distance, do not need to
	 * 2 updating cost is high, avoid rebuilding the index again and again.
	 * can be a pami paper, more theory, and combine with "Metric Space Similarity Joins"
	 */
	public static ArrayList<Double> IncrementalJoin(double [][]point1xys, double [][]point2xys, int dimension, indexNode X, indexNode Y, int splitOption, int fastMode, 
    		double error, boolean reverse, double directDis, boolean topkEarlyBreaking, double hausdorff, PriorityQueue<queueMain> aHeaps, 
    		Map<Integer, indexNode> nodelist, Map<Integer, indexNode> nodelist1, String folder) {
    	if(splitOption==0)
    		fastMode = 0;
    	File myObj = new File("./logs/spadas/"+folder+"/distance.txt"); 
    	myObj.delete();
    	ArrayList<Double> distancearray = new ArrayList<Double>();
   // 	System.out.println("the size of queue is "+aHeaps.size());
    	Map<Integer, Integer> joinNearestMap = new TreeMap<Integer, Integer>();
    	while(!aHeaps.isEmpty()) {
    		queueMain mainq = aHeaps.poll();
    		PriorityQueue<queueSecond> secondHeaps = mainq.getQueue();
    		queueSecond secondq = secondHeaps.peek();
    		indexNode anode = mainq.getIndexNode();
    		double ub = mainq.getbound();
    		indexNode bnode;
    		if(secondq!=null)
    			bnode = secondq.getNode();
    		else
    			continue;
    		boolean asplit = AdvancedHausdorff.stopSplitCondition(anode, splitOption);
    		boolean bsplit = AdvancedHausdorff.stopSplitCondition(bnode, splitOption);
    		if( asplit == true &&  bsplit == true) {//if two nodes are null, i.e., they are points // we can change it to radius threshold    			
    		//	System.out.println(secondq.getPointId());
    			joinNearestMap.put(mainq.getpointID(), secondq.getPointId());
    			Util.write("./logs/spadas/"+folder+"/distance.txt", ub+"\n");
    			distancearray.add(ub);
    			//write the distance into file for later effectiveness study
    			if(hausdorff>ub)//this will not happen if users specify one parameter.
    				hausdorff = ub; // updating the new Hausdorff, smaller distance
    		}else{
    			if(secondq.getbound()>=hausdorff)// filter by lower bound to shrink the queue
    				continue;
    			// another trick is to utilize the nearest point in the dataset for each point, 
    			// and empty the queue but only maintain one, if yes, every point and node needs to store its nearest neighbor
    			boolean splitPr = AdvancedHausdorff.splitPriority(anode, bnode, 1, asplit, bsplit);//
    			if(splitPr) {
    				AdvancedHausdorff.traverseX(anode, point1xys, point2xys, secondHeaps, dimension, aHeaps, fastMode, topkEarlyBreaking, directDis, true, hausdorff, nodelist);
    			}else {
    				secondHeaps.poll();
    				ub = AdvancedHausdorff.traverseY(bnode, point1xys, point2xys, secondHeaps, dimension, aHeaps, ub, anode, mainq.getpointID(),fastMode, topkEarlyBreaking, directDis, true, hausdorff, nodelist1);
    			}
    		}
    	}
   // 	System.out.println(joinNearestMap.get(1));
    //	System.out.println(joinNearestMap.get(2));
    	myObj = new File("./logs/spadas/"+folder+"/join.txt"); 
    	myObj.delete();
    	combineJoinableDataset(point1xys, point2xys, joinNearestMap, dimension, "./logs/spadas/"+folder+"/join.txt");
    	return distancearray;
    }
	
	
	static void combineJoinableDataset(double [][]point1xys, double [][]point2xys, Map<Integer, Integer> joinNearestMap, int dimension, String fileName) {
		for(int ii = 1; ii<=point1xys.length; ii++) {
		//	System.out.println(ii);
			if (joinNearestMap.containsKey(ii)) {
				int b = joinNearestMap.get(ii);
				String conString = "";
				for (int jj = 0; jj < dimension; jj++) {
					conString += Double.toString(point1xys[ii - 1][jj]) + ",";
				}
				for (int jj = 0; jj < dimension; jj++) {
					conString += Double.toString(point2xys[b - 1][jj]) + ",";
				}
				conString = conString.substring(0, conString.length() - 1) + "\n";
				Util.write(fileName, conString);
			}else {
			//	System.out.println("aaaaa");
			}
		}
	}
	/*
	 * scan every point, can be used for join and outlier detection
	 */
	public static void scanning(double [][]query, double [][]datasetB, int dimension) {
		int i = 1;
		Map<Integer, Integer> joinNearestMap = new TreeMap<Integer, Integer>();
		for(double []queryPoint: query) {
			int j=1;
			double min_dis = Double.MAX_VALUE;
			int nearest_point = 0;
			for(double []dataPoint: datasetB) {//search the nearest neighbor;
				double distance = Util.EuclideanDis(queryPoint, dataPoint, dimension);
				if(distance<min_dis) {
					min_dis = distance;
					nearest_point = j;
				}
				j++;
			}
			joinNearestMap.put(i++, nearest_point);
		}
		System.out.println(joinNearestMap.get(1));
	}
}
