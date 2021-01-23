package edu.nyu.dss.similarity;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;

import au.edu.rmit.trajectory.clustering.kmeans.indexAlgorithm;
import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;

/*
 * measuring the distance between two set of points.
 * 
 * bound computation
 */
public class Hausdorff {
	
	static double radiusThreshold;
	static double coveredPointsThreshold;// for evaluation
	
	static boolean smallMainQueue = false;
	static int tightBound = 0;
	static boolean cacheDis = true;
	static int splitCompare = 1; //which can be evaluated later
	static int EuclideanCounter = 0;
	static double disCompTime = 0.0;
	static boolean apporixmate = false;
	
	
	
	static void setBoundChoice(int boundChoice) {
		tightBound = boundChoice;
	}
	/*
	 * compute the distance between two starting points of id.
	 * 
	 * ??? we can also have early breaking version here as we do not need to compute the full distance
	 */
	static double EuclideanDis(double []point1xy, double []point2xy, int dimension) {
		double sum = 0;
		for(int i=0; i<dimension; i++) {
			sum +=  Math.pow((point1xy[i]-point2xy[i]), 2);
		}
		return Math.sqrt(sum);
	}
	
	 /*
     * the exact Hausdorff distance, appliable for small dataset, make it for general use, hau determines whether to use two-level
     * 1) directed
     * 2) outlier
     * 3) two-tier hausdorff distance
     */
    public static double HausdorffExact(double [][]point1xys, double [][]point2xys, int dimension, boolean Hau) {
        double[][] dist_matrix;
        dist_matrix = new double[point1xys.length][point2xys.length];
        double result = 0.0D;
        ArrayList<Double> minDistances1 = new ArrayList<Double>();
        ArrayList<Double> minDistances2 = new ArrayList<Double>();
        int i;
        for (i = 0; i < dist_matrix.length; ++i) {
            for (int j = 0; j < dist_matrix[0].length; ++j) {
				if (dimension > 1 && Hau) {//use Hausdorff for two vectors that may do not have same length
					double[][] transfrom1 = new double[point1xys[i].length][];
					double[][] transfrom2 = new double[point2xys[i].length][];
					for (int ii = 0; ii < point1xys[i].length; ii++) {
						transfrom1[ii] = new double[1];
						transfrom1[ii][0] = point1xys[i][ii];
					}
					for (int ii = 0; ii < point2xys[i].length; ii++) {
						transfrom2[ii] = new double[1];
						transfrom2[ii][0] = point2xys[i][ii];
					}
					dist_matrix[i][j] = HausdorffExact(transfrom1, transfrom2, 1, Hau);
				}else {//use Euclidean for vectors that have the same length.
					dist_matrix[i][j] = EuclideanDis(point1xys[i], point2xys[j], dimension);
				}
            }
        }

        int j;
        double min;
        for (i = 0; i < dist_matrix.length; ++i) {
            min = Double.MAX_VALUE;
            for (j = 0; j < dist_matrix[0].length; ++j) {
                if (dist_matrix[i][j] <= min) {
                    min = dist_matrix[i][j];
                }
            }
            minDistances1.add(min);
        }

        for (i = 0; i < dist_matrix[0].length; ++i) {
            min = Double.MAX_VALUE;
            for (j = 0; j < dist_matrix.length; ++j) {
                if (dist_matrix[j][i] <= min) {
                    min = dist_matrix[j][i];
                }
            }
            minDistances2.add(min);
        }
        Collections.sort(minDistances1);
        Collections.sort(minDistances2);
        double value1 = minDistances1.get(minDistances1.size()-1);// we use the last k for outliers
        double value2 = minDistances2.get(minDistances2.size()-1);
        result = Math.max(value1, value2);
        return result;
    }
    
    /*
     * the discrete frechet distance, check whether it obeys triangle inequality
     */
	public static double Frechet(double[][] t1, double[][]  t2, int dimension) {
		double[][] ca = new double[t2.length][t1.length];
		for (int i = 0; i < t2.length; ++i) {
			for (int j = 0; j < t1.length; ++j) {
				ca[i][j] = -1.0D;
			}
		}
		return c(t2.length - 1, t1.length - 1, ca, t1, t2, dimension);
	}

	private static double c(int i, int j, double[][] ca, double[][] t1, double[][] t2, int dimension) {
		if (ca[i][j] > -1.0D)
			return ca[i][j];
		if (i == 0 && j == 0) {
			ca[i][j] = EuclideanDis(t1[0], t2[0], dimension);
		} else if (j == 0) {
			ca[i][j] = Math.max(c(i - 1, 0, ca, t1, t2, dimension), EuclideanDis(t2[i], t1[0], dimension));
		} else if (i == 0) {
			ca[i][j] = Math.max(c(0, j - 1, ca, t1, t2, dimension), EuclideanDis(t2[0], t1[j], dimension));
		} else {
			ca[i][j] = Math.max(Math.min(Math.min(c(i - 1, j, ca, t1, t2, dimension), c(i - 1, j - 1, ca, t1, t2, dimension)), 
							c(i, j - 1, ca, t1, t2, dimension)),
					EuclideanDis(t2[i], t1[j], dimension));
		}
		return ca[i][j];
	}
    
    /*
     * we compute the initial point by randomly selecting a point.
     * //do the sampling to get a bigger lb for pruning, it is updated when a point shows in the end.
     */
    static double computeInitialLowerBound(double [][]point1xys, double [][]point2xys, int dimension, indexNode X, indexNode Y, int option) {
    	double lb = 0;
    	if(option==1) {// random selected points can be a choice
    		Random rand = new Random();
            // Generate random integers in range 0 to 999 
            int rand_int1 = rand.nextInt(point1xys.length); 
            double []point = point1xys[rand_int1];
            indexAlgorithm<Object> indexDSS = new indexAlgorithm<Object>();
            double []minDistnearestID = {Double.MAX_VALUE,0};
			indexDSS.NearestNeighborSearchBall(point, Y, dimension, point2xys, minDistnearestID);
			lb = minDistnearestID[0];
    	}else if(option==2) {// greedly select points that are far to another datasets
    		
    	}
    	return lb;
    }
    
    /*
     * based on radius, leaf (capacity), or depth, or covered points
     */
    static boolean stopSplitCondition(indexNode a, int option) {
    	if(option==0) {
    		return a.isLeaf();
    	}else if(option==1) {
    		if(a.getRadius()>radiusThreshold)
    			return false;
    		else {
				return true;
			}
    	}else if(option==2) {
    		if(a.getTotalCoveredPoints()>coveredPointsThreshold)
    			return false;
    		else {
				return true;
			}
    	}
    	return false;
    }
    
    /*
     * split the first node, or the second node
     */
    static boolean splitPriority(indexNode a, indexNode b, int option) {
    	/*
    	if(a.isLeaf())
    		return false;
    	else if(b.isLeaf())
    		return true;
    	*/
    	
    	if(a==null)
    		return false;
    	else if(b==null)
    		return true;
    	
    	if(option==1) {
    		if(a.getRadius()>b.getRadius())
    			return true;
    		else {
				return false;
			}
    	}else if(option==2) {
    		if(a.getTotalCoveredPoints()>b.getTotalCoveredPoints())
    			return true;
    		else {
				return false;
			}
    	}
    	return false;
    }
    
    /*
     * using a brute force to compute the distance
     */
    static double computeNodeDistance(indexNode a, indexNode b, int dimension, double[][] dataMatrixa, double[][] dataMatrixb, boolean direction) {
    	// get all the covered points of two nodes
    	ArrayList<Integer> aCoveredPoints = new ArrayList<Integer>();
    	a.getAllCoveredPoints(aCoveredPoints);
    	ArrayList<Integer> bCoveredPoints = new ArrayList<Integer>();
    	b.getAllCoveredPoints(bCoveredPoints);
    	double [][] aDataMatrix = new double[aCoveredPoints.size()][dimension];
    	double [][] bDataMatrix = new double[bCoveredPoints.size()][dimension];
    	int i=0;
    	for(int idx:aCoveredPoints) {
    		aDataMatrix[i++] = dataMatrixa[idx-1];
    	}
    	i=0;
    	for(int idx:bCoveredPoints) {
    		bDataMatrix[i++] = dataMatrixb[idx-1];
    	}
    	return earlyBreaking(aDataMatrix, bDataMatrix, dimension, direction);// take the early breaking idea
    }
    
    
    /*
     * ??? current cost is too high! to be implemented
     * use index to accelerate the processing, we use ball-tree or any index, we can store two vectors max 
     * and min in each dimension, the cost to update is high
     * directed version currently
     */
    public static double HausdorffWithIndexBFS(double [][]point1xys, double [][]point2xys, int dimension, indexNode X, indexNode Y, int splitOption) {
    	double lb = computeInitialLowerBound(point1xys, point2xys, dimension, X, Y, 1);
    	ArrayList<indexNode> checkList = new ArrayList<indexNode>();
    	PriorityQueue<queueForNodes> aHeaps = new PriorityQueue<queueForNodes>();
    	queueForNodes qnode = new queueForNodes(X, Y);
    	X.increaseCounter();
    	aHeaps.add(qnode);
    	int iteration_times = 0;
    	while(!aHeaps.isEmpty()) {
    		queueForNodes nodePair = aHeaps.poll();
    		indexNode anode = nodePair.anode;
    		indexNode bnode = nodePair.bnode;
    		if(checkList.contains(anode))//check whether node is in the checked list, continue if yes.
    			continue;
    		double ub = nodePair.getbound();
    		if(iteration_times++ % 100 == 0)
    			System.out.println(ub+","+lb+","+iteration_times);
    		anode.decreaseCounter();
    		if(ub <= lb) {
    			checkList.add(anode);// add the node into the checked list for future pruning directly
    			anode.getAllCoveredNodes(checkList);// maybe we also need to add all the children node, as they has been prune
    			continue;
    		}
    		
    		boolean stopA = stopSplitCondition(anode, splitOption);
    		boolean stopB = stopSplitCondition(bnode, splitOption);
    		if(stopA && stopB) {
    			if(anode.getMinDis() < nodePair.getLowerBound())//no need to compute the exact
    				continue;
    			double distance = computeNodeDistance(anode, bnode, dimension, point1xys, point2xys, true);// compute the distance directly
    			anode.updateBound(distance);
    			if(anode.getcounter() == 0 && anode.getMinDis() > lb) {// update the lb in a proper case, if the queue does not have any other node
    				lb = anode.getMinDis();
    			}
    		}else {
    			boolean splitPriority = splitPriority(anode, bnode, splitOption);
    			if(splitPriority) {
    				for(indexNode a: anode.getNodelist()) {
    					qnode = new queueForNodes(a, bnode);
    					aHeaps.add(qnode);
    					a.increaseCounter();
    				}
    			}else {
    				for(indexNode b: bnode.getNodelist()) {
    					qnode = new queueForNodes(anode, b);
    					aHeaps.add(qnode);
    					anode.increaseCounter();
    				}
    			}
    		}
    	}
    	return lb;
    }
    
    /*
     * vldb11 method, but with a general index, ball-tree, can work for any types of index tree.
     * An incremental Hausdorff distance calculation algorithm
     * 
     */
    public static double IncrementalDistance(double [][]point1xys, double [][]point2xys, int dimension, indexNode X, indexNode Y, int splitOption) {
    	PriorityQueue<queueMain> aHeaps = new PriorityQueue<queueMain>();
    	PriorityQueue<queueSecond> secondHeaps = new PriorityQueue<queueSecond>();
    	queueSecond secondq = new queueSecond(Y, 0, 0);
    	secondHeaps.add(secondq);
    	double distance = 0;
    	queueMain mainq = new queueMain(X, secondHeaps, Double.MAX_VALUE, 0);
    	aHeaps.add(mainq);
    	while(!aHeaps.isEmpty()) {
    		mainq = aHeaps.poll();
    		secondHeaps = mainq.getQueue();
    		secondq = secondHeaps.peek();
    		indexNode anode = mainq.getIndexNode();
    		double ub = mainq.getbound();
    		indexNode bnode = secondq.getNode();
    		if(anode==null && bnode==null) {//if two nodes are null, i.e., they are points // we can change it to radius threshold
    			distance = ub;
    			break;
    		}else {
    			boolean splitPriority = splitPriority(anode, bnode, splitOption);//
    			if(splitPriority) {
    				traverseX(anode, point1xys, point2xys, secondHeaps, dimension, aHeaps);
    			}else {
    				secondHeaps.poll();
    				ub = traverseY(bnode, point1xys, point2xys, secondHeaps, dimension, aHeaps, ub, anode, mainq.getpointID());
    			}
    		}
    	}
    	return distance;
    }
    
    static void traverseX(indexNode anode, double[][] point1xys, double[][] point2xys, PriorityQueue<queueSecond> secondHeaps, 
    		int dimension, PriorityQueue<queueMain> aHeaps) {
    	ArrayList<double[]> pivtoList = new ArrayList<double[]>();
		ArrayList<Double> radiusArrayList = new ArrayList<Double>();
		if(anode.isLeaf()) {
			for (int a : anode.getpointIdList()) {
				pivtoList.add(point1xys[a-1]);
				radiusArrayList.add(0.0);
			}
		}else {
			for (indexNode a : anode.getNodelist()) {
				pivtoList.add(a.getPivot());
				radiusArrayList.add(a.getRadius());
			}
		}
		int length = pivtoList.size();
		for (int i=0; i<length; i++) {
			PriorityQueue<queueSecond> newsecondHeaps = new PriorityQueue<queueSecond>();
			double minub = Double.MAX_VALUE;
			for (queueSecond aQueueSecond : secondHeaps) {
				indexNode newnodeb = aQueueSecond.getNode();
				double pivot_distance, bradius = 0;
				if (newnodeb != null) {
					pivot_distance = newnodeb.getPivotdistance(pivtoList.get(i));
					bradius = newnodeb.getRadius();
				} else {
					pivot_distance = Util.EuclideanDis(point2xys[aQueueSecond.getPointId()-1], pivtoList.get(i), dimension);
				}
				double newlb = pivot_distance - bradius - radiusArrayList.get(i);
				queueSecond secondqb = new queueSecond(newnodeb, newlb, aQueueSecond.getPointId());
				newsecondHeaps.add(secondqb);
				double newub = pivot_distance + bradius + radiusArrayList.get(i);
				if (newub < minub)
					minub = newub;
			}
			if(anode.isLeaf()) {
				ArrayList<Integer> arrayList = new ArrayList<Integer>(anode.getpointIdList());
				queueMain mainqa = new queueMain(null, newsecondHeaps, minub, arrayList.get(i));
				aHeaps.add(mainqa);
			}else {
				ArrayList<indexNode> arrayList = new ArrayList<indexNode>(anode.getNodelist());
				queueMain mainqa = new queueMain(arrayList.get(i), newsecondHeaps, minub, 0);
				aHeaps.add(mainqa);
			}
		}
    }
    
    static double traverseY(indexNode bnode, double[][] point1xys, double[][] point2xys, PriorityQueue<queueSecond> secondHeaps, 
    		int dimension, PriorityQueue<queueMain> aHeaps, double ub, indexNode anode, int apoint) {
    	ArrayList<double[]> pivtoList = new ArrayList<double[]>();
		ArrayList<Double> radiusArrayList = new ArrayList<Double>();
		if(bnode.isLeaf()) {
			for (int a : bnode.getpointIdList()) {
				pivtoList.add(point2xys[a-1]);
				radiusArrayList.add(0.0);
			}
		}else {
			for (indexNode a : bnode.getNodelist()) {
				pivtoList.add(a.getPivot());
				radiusArrayList.add(a.getRadius());
			}
		}
		int length = pivtoList.size();
		for (int i=0; i<length; i++) {
			double pivot_distance, bradius = 0;
			if (anode != null) {
				pivot_distance = anode.getPivotdistance(pivtoList.get(i));
				bradius = anode.getRadius();
			} else {
				pivot_distance = Util.EuclideanDis(point1xys[apoint-1], pivtoList.get(i), dimension);
			}
			if(bnode.isLeaf()) {
				ArrayList<Integer> arrayList = new ArrayList<Integer>(bnode.getpointIdList());
				queueSecond mainqa = new queueSecond(null, pivot_distance-bradius-radiusArrayList.get(i), arrayList.get(i));
				secondHeaps.add(mainqa);
			}else {
				ArrayList<indexNode> arrayList = new ArrayList<indexNode>(bnode.getNodelist());
				queueSecond mainqa = new queueSecond(arrayList.get(i), pivot_distance-bradius-radiusArrayList.get(i), 0);
				secondHeaps.add(mainqa);
			}
			double newub = pivot_distance + bradius + radiusArrayList.get(i);
			if (newub < ub)
				ub = newub;
		}
		queueMain mainqa = new queueMain(anode, secondHeaps, ub, apoint);
		aHeaps.add(mainqa);
		return ub;
    }
    
    /*
     * depth first, scan another tree again and again
     */
    public static double HausdorffWithIndexDFS(double [][]point1xys, double [][]point2xys, int dimension, indexNode X, 
    		indexNode Y, double cmin, double ub, indexAlgorithm<Object> indexDSS) {
    	if(ub<=cmin) {
    		return cmin;
    	}else if (X.isLeaf()) {
    		for(int pointid:X.getpointIdList()) {
    			double []point = point1xys[pointid-1];
    			double []minDistnearestID = {Double.MAX_VALUE,0};
    			indexDSS.NearestNeighborSearchBall(point, Y, dimension, point2xys, minDistnearestID);
    			if(cmin<minDistnearestID[0])
    				cmin=minDistnearestID[0];
    		}
    	}else {
    		for(indexNode child: X.getNodelist()) {
    			double radius = child.getRadius();
    			double[] pivot = child.getPivot();
    			double []minDistnearestID = {Double.MAX_VALUE,0};
    			indexDSS.NearestNeighborSearchBall(pivot, Y, dimension, point2xys, minDistnearestID);
    			double ub1 = minDistnearestID[0] + radius;
    			//change to basic bound, which is much loose
    			cmin = HausdorffWithIndexDFS(point1xys, point2xys, dimension, child, Y, cmin, ub1, indexDSS);
    		}
    	}
    	return cmin;
    }
    
    /*
     * pami15, early breaking
     */
    public static double earlyBreaking(double [][]point1xys, double [][]point2xys, int dimension, boolean directed) {
    	int x = point1xys.length;
    	int y = point2xys.length;
    	double[][] dist_matrix;
        dist_matrix = new double[x][y];
        double cmax = 0;
    	for(int i=0; i<x; i++) {
    		double cmin= Double.MAX_VALUE;
    		for(int j=0; j<y; j++) {
    			/*
    			 * we add an early breaking on Euclidean distance, good for high-dimension dataset 
    			 */
    			dist_matrix[i][j] = EuclideanDis(point1xys[i], point2xys[j], dimension);
    			if(dist_matrix[i][j]<cmin)
    				cmin = dist_matrix[i][j];
    			if(cmin<cmax)
    				break;
    		}
    		if(cmin>cmax && Double.MAX_VALUE>cmin)
    			cmax = cmin;
    	}
    	
    	if(directed==false)
	    	for(int j=0; j<y; j++) {
	    		double cmin= Double.MAX_VALUE;
	    		for(int i=0; i<x; i++) {
	    			if(dist_matrix[i][j] == 0)
	    				dist_matrix[i][j] = EuclideanDis(point1xys[i], point2xys[j], dimension);
	    			if(dist_matrix[i][j]<cmin)
	    				cmin = dist_matrix[i][j];
	    			if(cmin<cmax)
	    				break;
	    		}
	    		if(cmin>cmax && Double.MAX_VALUE>cmin)
	    			cmax = cmin;
	    	}
    	return cmax;
    }
    
    /*
     * A local start search algorithm to compute exact Hausdorff Distance for arbitrary point sets
     */
    void localSearch() {
    	
    }
    
    /*
     * diffusion search
     */
    void diffusionSearch() {
    	
    }
    
    /*
     * octree method
     * 
     * https://github.com/siddeshpillai/Octree/blob/master/Octree.java
     */
    void octree3D() {
    	
    }
    
    /*
     * index can be various as the pruning rules are general
     * 
     * https://github.com/conversant/rtree, any dimension R-tree
     * 
     * or quad-tree
     * 
     */
}
