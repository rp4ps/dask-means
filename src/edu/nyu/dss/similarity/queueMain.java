package edu.nyu.dss.similarity;


import java.util.PriorityQueue;

import au.edu.rmit.trajectory.clustering.kmeans.indexNode;

/*
 * this class is for the heap method
 */
public class queueMain implements Comparable<queueMain>{
	double ub;
	indexNode anode;
	int pointid;
	PriorityQueue<queueSecond> aSecond;

	
	public double getbound() {
		return ub;
	}
	
	public indexNode getIndexNode() {
		return anode;
	}
	
	public int getpointID() {
		return pointid;
	}
	
	public PriorityQueue<queueSecond> getQueue(){
		return aSecond;
	}
	
	public queueMain(indexNode a, PriorityQueue<queueSecond> bSecond, double ub, int pointid) {
		this.anode = a;
		this.aSecond = bSecond;
		this.ub = ub;
		this.pointid = pointid;
	}
	
	@Override
    public int compareTo(queueMain other) {
		double gap = this.getbound() - other.getbound();
		
		if(gap>0)
			return -1;
		else if(gap==0)
			return 0;
		else {
			return 1;
		}
    }
	
}
