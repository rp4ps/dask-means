package edu.nyu.dss.similarity;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import au.edu.rmit.trajectory.clustering.kmeans.indexAlgorithm;
import au.edu.rmit.trajectory.clustering.kmeans.indexNode;

/*
 * this is for top-k dataset search based on similarity search
 * 
 * 3d object dataset: http://yulanguo.me/dataset.html
 * 2d dataset: porto, tdrive, 
 */
public class Search {
	
	void kNN() {
		/*
		 * we evaluate multiple queries
		 */
	}
	
	/*
	 * select random queries to test all the functions
	 */
	void generateQueries() {
		//randomly select a set of datasets, or selecting subset of each dataset to compose new
	}
	
	/*
	 * MBR bound from index
	 */
	void vldb17() {
		
	}
	
	/*
	 * expansion based method with grid index
	 */
	void sigir18() {
		
	}
	
	/*
	 * index-nodes
	 */
	void sigspatial11() {
		
	}
	
	/*
	 * incremental search using bounding box, we can 
	 */
	void vldb11() {
		
	}
	
	/*
	 * lower bounds one by one, integrating a bound into the algorithm, as we know the gap between lower and upper, a+2b, then we 
	 * we can prune by infering the lower bound, need to debug this function.
	 */
	static void hausdorffEarylyAbandonging(Map<Integer, double[][]> dataMap, int querySet, Map<Integer, indexNode> indexMap, 
			int dimension, int limit, boolean topkEarlyBreaking, Map<Integer, indexNode> nodelist, Map<Integer, indexNode> nodelist1) {
		double min_dis = Double.MAX_VALUE;
		int datasetid = 0;
		for(int a:dataMap.keySet()) {
			if(a>limit)
				break;
			Pair<Double, PriorityQueue<queueMain>> aPair;
			if(topkEarlyBreaking)
				aPair = AdvancedHausdorff.IncrementalDistance(dataMap.get(querySet), dataMap.get(a), dimension, indexMap.get(querySet), indexMap.get(a), 1, 2, 0.05, false, min_dis, topkEarlyBreaking, nodelist, nodelist1);
			else
				aPair = AdvancedHausdorff.IncrementalDistance(dataMap.get(querySet), dataMap.get(a), dimension, indexMap.get(querySet), indexMap.get(a), 0, 0, 0.01, false, min_dis, topkEarlyBreaking, nodelist, nodelist1);
			double distance = aPair.getLeft();
			if(min_dis>distance) {
				min_dis = distance;
				datasetid = a;
			}
		}
		System.out.println(min_dis+","+datasetid);
	}
	
	
	
	/*
	 * lower bounds one by one, compute the lower bound first, rank them, then scan one by one
	 */
	static int HausdorffEarylyAbandongingRanking(Map<Integer, double[][]> dataMap, int querySet, Map<Integer, indexNode> indexMap, 
			int dimension, int limit, Map<Integer, Pair<Double, PriorityQueue<queueMain>>> resultPair, 
			Map<Integer, indexNode> nodelist, Map<Integer, indexNode> nodelist1) {
		Map<Integer, Double> setBound = new HashMap<Integer, Double>();
		double radiusThreshold = 0.05; // the error can be bigger for argo
		for(int a:dataMap.keySet()) {
			if(a>limit)
				break;
			Pair<Double, PriorityQueue<queueMain>> aPair = AdvancedHausdorff.IncrementalDistance(dataMap.get(querySet), dataMap.get(a), dimension, indexMap.get(querySet), indexMap.get(a), 1, 2, radiusThreshold, false, 0, false, nodelist,nodelist1);
			setBound.put(a, aPair.getLeft()-radiusThreshold);
		}
		//LinkedHashMap preserve the ordering of elements in which they are inserted
		LinkedHashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
		setBound.entrySet().stream().sorted(Map.Entry.comparingByValue()).forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue()));
		double min_dis = Double.MAX_VALUE;
		int counter = 0;
		int datasetid = 0;
		Pair<Double, PriorityQueue<queueMain>> resPair = null;
		for(int a:sortedMap.keySet()) {
	//		System.out.println(sortedMap.get(a));
			if(min_dis<=setBound.get(a))//the bound
				break;
			Pair<Double, PriorityQueue<queueMain>> aPair = AdvancedHausdorff.IncrementalDistance(dataMap.get(querySet), dataMap.get(a), dimension, indexMap.get(querySet), indexMap.get(a), 0, 1, 0, false, min_dis, false, nodelist, nodelist1);
			double distance = aPair.getLeft();
	//		System.out.println(a+","+distance);
			if(min_dis>distance) {
				min_dis = distance;
				datasetid = a;
				resPair = aPair;
			}
			counter++;
		}
		System.out.println(min_dis+","+counter+","+datasetid);
		resultPair.put(datasetid, resPair);
		return datasetid;
	}
	
	/*
	 * scanning with index, prune if no intersection at all.
	 */
	ArrayList<Integer> rangeQuery(indexNode datalakeRoot, ArrayList<Integer> result, double querymax[], double querymin[], int dim, Map<Integer, indexNode> datalakeIndex) {
		//give a range, find all the dataset that intersect, we just need the data lake tree
		if(datalakeRoot.isrootLeaf()) {
			if (datalakeRoot.intersected(querymax, querymin, dim)) {
				result.add(datalakeRoot.getDatasetID());// add the results into, and ranked by intersected area.
			}
		}else {
			for(indexNode childNode: datalakeRoot.getNodelist(datalakeIndex)) {
				if(childNode.intersected(querymax, querymin, dim)) {
					result = rangeQuery(childNode, result, querymax, querymin, dim, datalakeIndex);
				}
			}
		}
		return result;
	}
	
	/*
	 * scanning with index, compute the intersection, and rank by overlapped range, or number of points.
	 */
	static HashMap<Integer, Double> rangeQueryRankingArea(indexNode datalakeRoot, HashMap<Integer, Double> result, double querymax[], double querymin[], 
			double min_dis, int k, HashMap<Integer, PriorityQueue<queueMain>> queues, int dim, Map<Integer, indexNode> datalakeIndex) {
		//give a range, find all the dataset that intersect, we just need the data lake tree and ranked by intersected areas
		if(datalakeRoot.isrootLeaf()) {
			double distance = -datalakeRoot.intersectedArea(querymax, querymin, dim);
			if (distance < min_dis) {
				result = holdingTopK(result, datalakeRoot.getDatasetID(), distance, k, queues, null);
				if (result.size() == k)
					min_dis = result.entrySet().iterator().next().getValue();
				else
					min_dis = Double.MAX_VALUE;
			}
		}else {
			for(indexNode childNode: datalakeRoot.getNodelist(datalakeIndex)) {
				double bound = childNode.intersectedArea(querymax, querymin, dim);
				if(-bound<min_dis && bound>0) {
					result = rangeQueryRankingArea(childNode, result, querymax, querymin, min_dis, k, queues, dim, datalakeIndex);
				}
			}
		}
		return result;
	}
	
	/*
	 * scanning with index, compute the intersection, and rank by overlapped range, or number of points.
	 */
	static HashMap<Integer, Double> rangeQueryRankingNumberPoints(indexNode datalakeRoot, HashMap<Integer, Double> result, double querymax[], double querymin[], 
			double min_dis, int k, HashMap<Integer, PriorityQueue<queueMain>> queues, int dim, Map<Integer, double[][]> dataMap, Map<Integer, indexNode> datalakeIndex, 
			Map<Integer, Map<Integer, indexNode>> datasetindex, Map<Integer, Integer> argoDataMapping, String indexString, indexAlgorithm<Object> indexDSS)  throws FileNotFoundException, IOException {
		if(datalakeRoot.isrootLeaf()) {
			int datasetid = datalakeRoot.getDatasetID();
			Map<Integer, indexNode> dataindex = indexDSS.restoreSingleIndex(indexString, datasetid, dim);
			indexNode dataseNode = dataindex.get(1);
			double[][] dataset;
			if (dataMap == null)
				dataset = Framework.readSingleFile(argoDataMapping.get(datasetid));
			else
				dataset = dataMap.get(datasetid);
			double distance = -dataseNode.coveredPOints(querymax, querymin, min_dis, dim, dataset, dataindex);
			if (distance < min_dis) {
				result = holdingTopK(result, datasetid, distance, k, queues, null);
				if (result.size() == k)
					min_dis = result.entrySet().iterator().next().getValue();
				else
					min_dis = Double.MAX_VALUE;
			}
		}else {
			for(indexNode childNode: datalakeRoot.getNodelist(datalakeIndex)) {
				int threshold = childNode.coveredPOints(querymax, querymin, dim);
				if(-threshold<min_dis && threshold>0) {
					result = rangeQueryRankingNumberPoints(childNode, result, querymax, querymin, min_dis, k, queues, dim, dataMap, datalakeIndex, datasetindex, argoDataMapping, indexString, indexDSS);
				}
			}
		}
		return result;
	}
	
	
	/*
	 * posting list, similar to term and document
	 */
	void createPostinglists(HashMap<Integer, ArrayList<Integer>> dataset) {
		HashMap<Integer, ArrayList<Integer>> postlisitngHashMap = new HashMap<Integer, ArrayList<Integer>>();
		for(int datasetid: dataset.keySet()) {
			ArrayList<Integer> signatureArrayList = dataset.get(datasetid);
			for(int code: signatureArrayList) {
				ArrayList<Integer> datasetList;
				if(postlisitngHashMap.containsKey(code))
					datasetList = postlisitngHashMap.get(code);
				else
					datasetList = new ArrayList<Integer>();
				if(!datasetList.contains(datasetid))
					datasetList.add(datasetid);
				postlisitngHashMap.put(code, datasetList);
			}
		}
	}
	/*
	 * the similarity is based on the intersection of z-curve code of each datasets.
	 * we design an inverted index to accelerate the pair-wise and top-k, for node, 
	 * we can get a upper bound, and it is more efficient then the inverted index only
	 */
	static HashMap<Integer, Double> gridOverlap(indexNode datalakeRoot, HashMap<Integer, Double> result, int[] query, 
			double min_dis, int k, HashMap<Integer, PriorityQueue<queueMain>> queues, Map<Integer, indexNode> datalakeIndex) {
		// give a query signature, find all the dataset that has a ratio, 
		if (datalakeRoot.isrootLeaf()) {
			int datasetid = datalakeRoot.getDatasetID();
			double distance = datalakeRoot.GridOverlap(query);
			if (distance < min_dis) {
				result = holdingTopK(result, datasetid, distance, k, queues, null);
				if (result.size() == k)
					min_dis = result.entrySet().iterator().next().getValue();
				else
					min_dis = Double.MAX_VALUE;
			}
		} else {
			for (indexNode childNode : datalakeRoot.getNodelist(datalakeIndex)) {
				if (childNode.GridOverlap(query) < min_dis) {
					result = gridOverlap(childNode, result, query, min_dis, k, queues, datalakeIndex);
				}
			}
		}
		return result;
	}
	
	
	/*
	 * access the posting list to filter other datasets, get a candidate list, then call the baseline
	 */
	HashMap<Integer, Double> gridOverlapPosstingList(ArrayList<Integer> query, HashMap<Integer, ArrayList<Integer>> postlisitngHashMap, 
			Map<Integer, ArrayList<Integer>> bitMap, int k) {
		Set<Integer> candidateDatasetArrayList = new HashSet<Integer>();
		HashMap<Integer, Double> result = new HashMap<Integer, Double>();
		HashMap<Integer, PriorityQueue<queueMain>> queues = new HashMap<Integer, PriorityQueue<queueMain>>();
		for(int code: query) {
			if(postlisitngHashMap.containsKey(code))
				candidateDatasetArrayList.addAll(postlisitngHashMap.get(code));
		}
		for(int datasetid: candidateDatasetArrayList) {
			ArrayList<Integer> arrayList = bitMap.get(datasetid);
			int counter = 0;
			for(int queryz: query) {
				if(arrayList.contains(queryz)) {
					counter++;
				}
			}
			double reverseoverlapRatio = 1 - counter/(double)query.size();
			result = Search.holdingTopK(result, datasetid, reverseoverlapRatio, k, queues, null);
		}
		return result;
	}
	
	
	
	/*
	 * we use hashmap to hold top-k
	 */
	static HashMap<Integer, Double> holdingTopK(HashMap<Integer, Double> result, int datasetid, double score, int k, 
			HashMap<Integer, PriorityQueue<queueMain>> queues, PriorityQueue<queueMain> usedQueue) {
		if(result.size()<k) {
			result.put(datasetid, score);
			if(queues!=null)
				queues.put(datasetid, usedQueue);
			if(result.size()==k) {
				HashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
				result.entrySet().stream().sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
						.forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue()));
				result = sortedMap;
			}
		}else {
			Entry<Integer, Double> a = result.entrySet().iterator().next();
			if (score < a.getValue()) {
				result.remove(a.getKey());
				if(queues!=null)
					queues.remove(a.getKey());
				result.put(datasetid, score);
				if(queues!=null)
					queues.put(datasetid, usedQueue);
				HashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
				result.entrySet().stream().sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
						.forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue()));
				result = sortedMap;
			}
		}
		return result;
	}
	
	// we prune those index and return candidates for further purning using sequential methods.
	static HashMap<Integer, Double> pruneByIndex(Map<Integer, double[][]> dataMap, indexNode datalakeRoot, indexNode query, int querySet, int dimension, Map<Integer, indexNode> indexMap,
			Map<Integer, Map<Integer, indexNode>> nodelist1, Map<Integer, indexNode> queryindexmap, Map<Integer, indexNode> datalakeIndex, 
			Map<Integer, Integer> argoDataMap, int k, String indexString) throws FileNotFoundException, IOException {
		PriorityQueue<queueForNode> aForNodes = new PriorityQueue<queueForNode>();
		queueForNode qNodes;
		if(datalakeIndex!=null)
			qNodes = new queueForNode(query, datalakeIndex.get(1));// root node
		else {
			qNodes = new queueForNode(query, datalakeRoot);
		}
		aForNodes.add(qNodes);
		double min_dis = Double.MAX_VALUE;
		int counter = 0;
		indexAlgorithm<Object> indexDSS = new indexAlgorithm<>();
		HashMap<Integer, Double> result = new HashMap<Integer, Double>();
		HashMap<Integer, PriorityQueue<queueMain>> queues = new HashMap<Integer, PriorityQueue<queueMain>>();
		while(!aForNodes.isEmpty()) {
			queueForNode aForNodes2 = aForNodes.poll();
			indexNode aIndexNode = aForNodes2.getNode();
			double lowerbound = aForNodes2.getbound();
			if(lowerbound>min_dis) {
				break;
			}
			if(aIndexNode.isrootLeaf()) {//we get to the dataset level
				int datasetid = aIndexNode.getDatasetID();
				// we can add the bound computation here, before we move on.
				double[][] querydata, dataset;
				if(dataMap!=null) {
					querydata = dataMap.get(querySet);
					dataset = dataMap.get(datasetid);
				}else {
					querydata = Framework.readSingleFile(argoDataMap.get(querySet));//read from files
					dataset = Framework.readSingleFile(argoDataMap.get(datasetid));
				}
				Pair<Double, PriorityQueue<queueMain>> aPair;
				if(nodelist1!=null)
					aPair = AdvancedHausdorff.IncrementalDistance(querydata, dataset, dimension, queryindexmap.get(1), nodelist1.get(datasetid).get(1), 0, 1, 0, false, min_dis, false, queryindexmap, nodelist1.get(datasetid));
				else {
					if(indexMap!=null)
						aPair = AdvancedHausdorff.IncrementalDistance(querydata, dataset, dimension, indexMap.get(querySet), indexMap.get(datasetid), 0, 1, 0, false, min_dis, false, null, null);
					else {
						Map<Integer, indexNode> dataindexMap = indexDSS.restoreSingleIndex(indexString, datasetid, dimension);
						aPair = AdvancedHausdorff.IncrementalDistance(querydata, dataset, dimension, queryindexmap.get(1), dataindexMap.get(1), 0, 1, 0, false, min_dis, false, queryindexmap, dataindexMap);
					}
				}
				double distance = aPair.getLeft();

				result = holdingTopK(result, datasetid, distance, k, queues, aPair.getRight());
				if (result.size() == k)
					min_dis = result.entrySet().iterator().next().getValue();
				else
					min_dis = Double.MAX_VALUE;
				counter++;
			}else {
				for(indexNode childIndexNode: aIndexNode.getNodelist(datalakeIndex)) {
					//	System.out.println(childIndexNode.getRadius());
					qNodes = new queueForNode(query, childIndexNode);
					aForNodes.add(qNodes);
				}
			}
		}
		for(int datastid: result.keySet()) {
			System.out.println(datastid+", "+result.get(datastid));
		}
		System.out.println("\n"+min_dis+", "+counter);
		return result;
	}
	
	// Implement a system which integrate everything, and the outlier version, interface, and think about the outlier of dataset
	
	void detectOutlierRefineIndex() {
		/*
		 * check all the leaf nodes, and see whether there are some big nodes, as index groups similar points together
		 * 
		 * few points with big radius will be judged as outlier
		 * 
		 * Not our research focus, we can just give the idea that putting outier alone will improve the performance.
		 * 
		 * using our own algorithm to detect the own outliers, as it is the still same framework.
		 * 
		 * top-n points distance to their nearest neigbhor, follow the vldb10's definition, and no parameter using our 
		 */
	}
}
