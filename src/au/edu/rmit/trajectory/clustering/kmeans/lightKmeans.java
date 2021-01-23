package au.edu.rmit.trajectory.clustering.kmeans;


import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kpaths.Util;

public class lightKmeans<T> extends kmeansAlgorithm<T>{

	double maxdrift[];
	short assigned[]; // store the assigned iterations
	short newassigned[]; // store the assigned iterations
	byte []quantilizedLowerBound=null;// create new only when leaf node is not pruned, or not used, optional
	byte []quantilizedUpperBound=null;// 
	short []pointcounterPruned=null;
	double unitQuantization;
	
	public lightKmeans(String[] datapath) {
		super(datapath);
		// TODO Auto-generated constructor stub
	}
	
	/*
	 * build index and only use partial centroid distance to accelerate k-means, saving space, 
	 * applications: data summarization when dataset is too large, and k is also large, such as sampling, 
	 * e.g., point cloud simplification in autonomous vehicles, ipad pro and iphone with lidar scanner
	 */
	
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
		quantilizedLowerBound = new byte[trajectoryNumber];
		quantilizedUpperBound = new byte[trajectoryNumber];
		pointcounterPruned = new short[trajectoryNumber];
		maxdrift = new double[k];
		unitQuantization = 0.01;
	}
	
	public void recursivePICKmeans(indexNode dataRoot, indexNode centroidNode, double [][]dataMatrix, 
			double [][]centroidMatrix, double ubfather, int centroid, double group_drift[], double interMinimumCentoridDis[]) {
		// traverse the datapoint index, and get the assginment
		if(iterationTimes==1)
			initializationPICKmeans();
		else {
			updatebound();// update lb and ub by the drift
		}
		for(int assignedID= 0; assignedID<k; assignedID++)
			for(int i=0; i<k; i++)//get the maximum except the assignedID.
				if(i!= assignedID && maxdrift[assignedID]<group_drift[i-1])
					maxdrift[assignedID] = group_drift[i-1];
		
		dualTreeKmeans(dataRoot, centroidNode, dataMatrix, centroidMatrix, unitQuantization, ubfather, centroid, interMinimumCentoridDis);
		cleanIndex(dataRoot);//update the assigned for nodes by calling the function
		assigned = newassigned;// all the points have 
	}

	
	/*
	 * update the bound of all the nodes
	 */
	private void updatebound() {
		// TODO Auto-generated method stub
		
	}

	/*
	 * A Dual-Tree Algorithm for Fast k-means Clustering With Large k, SDM 2017, we use multiple quantilized bounds
	 * 
	 * when k is small, we put all the points into the centroidnode.
	 */
	void dualTreeKmeans(indexNode dataRoot, indexNode centroidNode, double [][]dataMatrix, 
			double [][]centroidMatrix, double unit, double ubfather, int centroid, double interMinimumCentoridDis[]) {
		int scannedCentrorid = dataRoot.getPrunedCounter();
		if(scannedCentrorid == k)
			return;
		short assignedID = dataRoot.getAssignedCluster(); // check which cluster it assigned
		double ub = dataRoot.getUpperBound(unit); // it is maximum as initialized or if not pruned
		double lb = dataRoot.getLowerBound(unit);
		if(assignedID > 0 && ub < Math.max(lb, interMinimumCentoridDis[assignedID]/2.0)){// maintain in the node, this is for second and following iterations
			dataRoot.setNearestCluster(assignedID);//and it's child node, everyone
			dataRoot.setAllPrune((short)k, pointcounterPruned);// set the getPrunedCounter as k
			return;
		}
		double centroidData[] = null;
		double radius = 0;
		short pruneNum = 1;
		if(centroidNode!=null) {
			centroidData = centroidNode.getPivot();
			radius = centroidNode.getRadius();
			pruneNum = (short)centroidNode.getTotalCoveredPoints();
		}else
			centroidData = centroidMatrix[centroid-1];
		double pivotDis = Util.EuclideanDis(dataRoot.getPivot(), centroidData, dimension);
		double newub = pivotDis + dataRoot.getRadius() + radius; // we will change to half ball theory, which is tighter
		if(newub < ub) {//update the bound as we found a tighter bound,
			if(centroidNode==null) {
				dataRoot.setNearestCluster((short)centroid);
			//	dataRoot.setLowerbound(ub);//??
			}
			ub = newub;//Math.min(newub,ubfather); // cannot rank, but priority queue can rank, see which one is good, we investigate both
			dataRoot.setUpperbound(ub);
		}
		double nodelb = pivotDis - dataRoot.getRadius() - radius;// how to update the bound
		if(nodelb > ub) {// if the lower bound is bigger, we do not need to go deeper of this centroid node,
			dataRoot.setPrunedCounter(pruneNum, pointcounterPruned);// and all its children, pruned
			if(scannedCentrorid + pruneNum == k) {//if it equals k, we assign the node directly, we add the node to assign id
				short nearest = dataRoot.getnearestID();
				if(nearest != assignedID) {// assign the node to the sum if not include in the current node, 
					CENTERSEuc.get(assignedID).removeNode(dataRoot);
					CENTERSEuc.get(nearest).addNode(dataRoot);
				}
				return;
			}
		} else {// split the centroid nodes and datanodes, and access its child node
			if (dataRoot.isLeaf() && centroidNode == null) {// we will not go deeper as we store our data here
				for (int i : dataRoot.getpointIdList()) {
					short previousNearest = assigned[i-1];
					double currentDis = unit*quantilizedUpperBound[i-1];
					double second_min_dist = unit*quantilizedLowerBound[i-1];// we store the lb
					if(pointcounterPruned[i-1] == k || currentDis < Math.max(second_min_dist, interMinimumCentoridDis[newassigned[i-1]]/2.0)){// check lb and ub to see whether we can assign directly
						newassigned[i-1] = (short)previousNearest;
						pointcounterPruned[i-1] = (short)k;
						continue;
					}
					double dataPoint[] = dataMatrix[i - 1];
					double dis = Util.EuclideanDis(dataPoint, centroidData, dimension);// pair wise computation and pruning by checking the point bound, and update the node id
					if (dis < currentDis) {
						second_min_dist = currentDis;
						currentDis = dis;
						newassigned[i-1] = (short) centroid;
						setUpperBound(currentDis, unit, quantilizedUpperBound, i-1);
					} else if (dis < second_min_dist && dis != currentDis) {// store the second nearest, for next iteration pruning and node pruning
						second_min_dist = dis;
					}
					setLowerBound(second_min_dist, unit, quantilizedLowerBound, i-1);// store the second nearest
					if (pointcounterPruned[i-1] + pruneNum == k) {// all the centroids have been scanned for this point
						if (newassigned[i-1] != previousNearest) {
							CENTERSEuc.get(previousNearest).deleteSinglePointToCluster(i-1);
							CENTERSEuc.get(newassigned[i-1]).addPointToCluster(i-1);
							CENTERSEuc.get(newassigned[i-1]).addSum(dataPoint);
							CENTERSEuc.get(previousNearest).minusSum(dataPoint);		
						}
					} else {
						pointcounterPruned[i-1]++;
					}
				}
			} else {
				if (centroidNode != null) {
					if (centroidNode.isLeaf())
						for (int childpoint : centroidNode.getpointIdList())
							dualTreeKmeans(dataRoot, null, dataMatrix, centroidMatrix, unit, pivotDis+centroidNode.getRadius(), childpoint, interMinimumCentoridDis);
					else
						if(radius < dataRoot.getRadius() && !dataRoot.isLeaf())// which one is good, check which node to split by comparing the radius
							for (indexNode childNode : dataRoot.getNodelist())
								dualTreeKmeans(childNode, centroidNode, dataMatrix, centroidMatrix, unit, pivotDis+childNode.getDisFather(), centroid, interMinimumCentoridDis);//sort by the distance
						else
							for (indexNode childNode : centroidNode.getNodelist())
								dualTreeKmeans(dataRoot, childNode, dataMatrix, centroidMatrix, unit, pivotDis+childNode.getDisFather(), centroid, interMinimumCentoridDis);
				} else
					for (indexNode childNode : dataRoot.getNodelist())
						dualTreeKmeans(childNode, null, dataMatrix, centroidMatrix, unit, pivotDis + childNode.getDisFather(), centroid, interMinimumCentoridDis);
			}
		}
	}
	
	/*
	 * this function we search by 2nn and nn in varius iterations, with bound on the range
	 * we implement our pruning based on two nearest neighbor search and bound-based similarity search
	 */
	void dualTreeKmeansNN(indexNode dataRoot, indexNode centroidNode, double [][]dataMatrix, 
			double [][]centroidMatrix, double unit, double ubfather, int centroid) {
		// use our 
	}
	
	private void setUpperBound(double currentDis, double unit, byte[] quantilizedUpperBound2, int i) {
		quantilizedUpperBound2[i] = (byte)(currentDis/unit+1);
	}

	private void setLowerBound(double second_min_dist, double unit, byte[] quantilizedLowerBound2, int i) {
		quantilizedLowerBound2[i] = (byte)(second_min_dist/unit-1);
	}

	// when to split the centroid index, whether we should split all, update the bound of each node
	void cleanIndex(indexNode dataRoot) {
		// update the ub and lb by the drift based on their assignment
		int assignedID = dataRoot.getAssignedCluster();
		double drift = 0, otherDrift = 0;
		if(assignedID < 0) {
			drift = group_drift[assignedID-1];
			otherDrift = maxdrift[assignedID-1];
		}
		
		// update if the assigned point is not 0, update the child's bound, and its assigned
		
		// label the node if not assign by checking pruned nodes,
	}
	
	/*
	 * 	 we will infer the best size for node, if centroids are not that many, we put them in one leaf node
	 */
	void parameterization(int memoryCapacity) {
		// we set the parameter of index leaf node and unit for quantization, and whether to use the lb for maximum pruning, and 
		// unit
		// index capacity
		// storing the lb or not, as it costs much to 
	}
	
	//to do make it run with inner bound first, than integrate the lb and father distance bound, and test with method
}
