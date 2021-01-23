package skyline0623.balltree;

import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;
import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import java.util.AbstractMap.SimpleEntry;

public class Hypersphere extends Point {
	private double radius;
	private LinkedList<Integer> instances;
	private Hypersphere[] children;
	private double[] sumOfPoints;
	static int COUNT = 0, ALL_COUNT = 0;
	
	Hypersphere(Point center, double r, LinkedList<Integer> ins){
		super(center.pos);
		this.radius = r;
		this.instances = ins;
		sumOfPoints = new double[Process.DIMENSION];
	}
	
	Hypersphere(){
		super(new double[Process.DIMENSION]);
		instances = new LinkedList<Integer>();
		sumOfPoints = new double[Process.DIMENSION];
	}
	
	void addInstance(int index){
		instances.add(index);
		double[] pos = Process.INSTANCES.get(index).getPosition();
		for(int i = 0; i < Process.DIMENSION; i++){
			sumOfPoints[i] += pos[i];
		}
	}
	void addInstance(int index, Map<Integer, double[]> datamapEuc){
		instances.add(index);
		double[] pos = datamapEuc.get(index);
		for(int i = 0; i < Process.DIMENSION; i++){
			sumOfPoints[i] += pos[i];
		}
	}
	
	void endAdding(){
		int size = instances.size();
		for(int i = 0; i < Process.DIMENSION; i++){
			this.pos[i] = this.sumOfPoints[i] / size;
		}
		this.radius = Point.euclideanDistance(this, Process.INSTANCES.get(this.getFarestPoint(this)));
	}
	
	int size(){
		return instances.size();
	}
	
	double maxDistance(Point p){
		return radius + Point.euclideanDistance(p, this);
	}
	
	double minDistance(Point p){
		return Point.euclideanDistance(p, this) - radius;
	}
	
	//��������ڵ�����cluster�У��ͷ���-1 ���򷵻�cluster center��index
	int isInSingleCluster(){
		ALL_COUNT++;
		PriorityQueue<Entry<Integer, Double>> maxpq = new PriorityQueue<Entry<Integer, Double>>(Process.CENTERS.size(), new Comparator<Entry<Integer, Double>>(){
			public int compare(Entry<Integer, Double> e1, Entry<Integer, Double> e2){
				double d1 = e1.getValue(), d2 = e2.getValue();
				if(d1 > d2){
					return 1;
				}
				if(d1 < d2){
					return -1;
				}
				return 0;
			}
		});
		PriorityQueue<Entry<Integer, Double>> minpq = new PriorityQueue<Entry<Integer, Double>>(Process.CENTERS.size(), new Comparator<Entry<Integer, Double>>(){
			public int compare(Entry<Integer, Double> e1, Entry<Integer, Double> e2){
				double d1 = e1.getValue(), d2 = e2.getValue();
				if(d1 > d2){
					return 1;
				}
				if(d1 < d2){
					return -1;
				}
				return 0;
			}
		});
		int index = 0;
		for(ClusteringCenter cen : Process.CENTERS){
			maxpq.add(new SimpleEntry<Integer, Double>(index, this.maxDistance(cen)));
			minpq.add(new SimpleEntry<Integer, Double>(index, this.minDistance(cen)));
			index++;
		}
		Entry<Integer, Double> the = maxpq.poll(), comp;
		index = the.getKey();
		double theDist = the.getValue();
		while((comp = minpq.poll()) != null){
			int ind = comp.getKey();
			double dis = comp.getValue();
			if(theDist < dis){
				if(ind != index){
					COUNT++;
					return index;
				}
				else
					continue;
			}
			else{
				if(ind == index)
					continue;
				return -1;
			}
		}
		return -1;
	}
	
	private int getFarestPoint(Point p){
		double maxDist = 0.0;
		int maxIndex = -1;
		for(int i : this.instances){
			Point pp = Process.INSTANCES.get(i);
			double dist = Point.euclideanDistance(p, pp);
			if(dist >= maxDist){
				maxDist = dist;
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	//split and store it to this node's children field, & return the children.
	Hypersphere[] split(){
		int firstCenter = this.getFarestPoint(this);
		Point fir = Process.INSTANCES.get(firstCenter);
		int secondCenter = this.getFarestPoint(fir);
		Point sec = Process.INSTANCES.get(secondCenter);
		this.children = new Hypersphere[2];
		this.children[0] = new Hypersphere();
		this.children[1] = new Hypersphere();
		this.children[0].addInstance(firstCenter);
		this.children[1].addInstance(secondCenter);
		for(int i : this.instances){
			if(i == firstCenter || i == secondCenter)
				continue;
			Point p = Process.INSTANCES.get(i);
			double dist1 = Point.euclideanDistance(p, fir),
					dist2 = Point.euclideanDistance(p, sec);
			if(dist1 < dist2){
				this.children[0].addInstance(i);
			}
			else{
				this.children[1].addInstance(i);
			}
		}
		this.children[0].endAdding();
		this.children[1].endAdding();
		return this.children;
	}
	
	Hypersphere[] getChildren(){
		return this.children;
	}
	
	static void locateAndAssign(Hypersphere hp){
		int clusterIndex = hp.isInSingleCluster();
		if(clusterIndex != -1){
			ClusteringCenter cc = Process.CENTERS.get(clusterIndex);
			for(int pi : hp.instances){
				cc.addPointToCluster(pi);
			}
			return;
		}
		if(hp.children == null){
			for(int pi : hp.instances){
				Point p = Process.INSTANCES.get(pi);  
				double minDist = Double.MAX_VALUE;
				int minCenIndex = 0, index = 0;
				for(ClusteringCenter cc : Process.CENTERS){
					double dist = Point.euclideanDistance(p, cc);
					if(dist < minDist){
						minDist = dist;
						minCenIndex = index;
					}
					index++;
				}
				ClusteringCenter cen = Process.CENTERS.get(minCenIndex);
				cen.addPointToCluster(pi);
			}
		}
		else{
			for(Hypersphere chp : hp.children){
				Hypersphere.locateAndAssign(chp);
			}
		}
	}

	public int traverseConvert(indexNode rootKmeans, int dimension) {
		rootKmeans.setRadius(radius);
		rootKmeans.setPivot(pos);//
		rootKmeans.setSum(sumOfPoints);
		if(children == null){//for the objects
			Set<Integer> aIntegers = new HashSet<Integer>();
			for(int id : instances) {
				aIntegers.add(id+1);// the pointid				
			}
			rootKmeans.addPoint(aIntegers);		
			rootKmeans.setTotalCoveredPoints(aIntegers.size());
			return aIntegers.size();
		}else {
			int count = 0;
			for(Hypersphere chp : children){
				indexNode childnodekmeans = new indexNode(dimension);
				rootKmeans.addNodes(childnodekmeans);
				count += chp.traverseConvert(childnodekmeans, dimension);
			}
			rootKmeans.setTotalCoveredPoints(count);
			return count;
		}		
	}
}
