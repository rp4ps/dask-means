package edu.nyu.dss.similarity;


import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;

/*
 * this class is for the heap method
 */
public class queueForNode implements Comparable<queueForNode>{
	double bound;
	indexNode anode;
	indexNode bnode;

	
	public double getbound() {
		return bound;
	}
	
	public indexNode getNode() {
		return bnode;
	}

	
	public queueForNode(indexNode a, indexNode b) {
		this.anode = a;
		this.bnode = b;
		double[] pivota = anode.getPivot();
		double[] pivotb = bnode.getPivot();
		bound = Util.EuclideanDis(pivota, pivotb, pivota.length);
		bound = Math.max(bound - bnode.getRadius(), 0);
	}
	
	@Override
    public int compareTo(queueForNode other) {
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
