package edu.nyu.dss.similarity;


import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kmeans.indexNode;

/*
 * this class is for the heap method
 */
public class queueSecond implements Comparable<queueSecond>{
	double bound;
	indexNode anode;
	int pointid;

	
	public double getbound() {
		return bound;
	}
	
	public indexNode getNode() {
		return anode;
	}
	
	public int getPointId() {
		return pointid;
	}
	
	public queueSecond(indexNode a, double b, int id) {
		this.anode = a;
		this.bound = b;
		this.pointid = id;
	}
	
	@Override
    public int compareTo(queueSecond other) {
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
