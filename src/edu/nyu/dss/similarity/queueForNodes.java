package edu.nyu.dss.similarity;


import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;

/*
 * this class is for the heap method
 */
public class queueForNodes implements Comparable<queueForNodes>{
	double bound;
	indexNode anode;
	indexNode bnode;

	
	public double getbound() {
		return bound;
	}
	
	public double getLowerBound() {
		return bound - 2*anode.getRadius() - 2*bnode.getRadius();
	}
	
	public queueForNodes(indexNode a, indexNode b) {
		this.anode = a;
		this.bnode = b;
		double[] pivota = anode.getPivot();
		double[] pivotb = bnode.getPivot();
		bound = Util.EuclideanDis(pivota, pivotb, pivota.length);
		bound += anode.getRadius() + bnode.getRadius();
	}
	
	@Override
    public int compareTo(queueForNodes other) {
		double gap = this.getbound() - other.getbound();
		
		if(gap>0)
			return 1;
		else if(gap==0)
			return 0;
		else {
			return -1;
		}
    }
	
}
