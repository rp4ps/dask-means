package edu.nyu.unik.expriments;

import java.io.IOException;

import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kmeans.kmeansAlgorithm;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;
import java_cup.internal_error;

public class kmeansEfficiencySynthetic {

	public kmeansEfficiencySynthetic() {
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) throws IOException, KeySizeException, KeyDuplicateException {
		// TODO Auto-generated method stub
		int[] kvalue = new int[]{10, 50, 100, 200, 400, 600, 800, 1000};//10, 100, 1000
		int[] scales = new int[] {100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
		int[] capacities = new int[] {10, 20, 30, 40, 50, 60};//capacity
		int[] dimensions = new int[] {10, 20, 30, 40, 50};
		kvalue = new int[]{10, 100, 400, 700, 1000};//10, 100, 1000
		scales = new int[] {10000};
		capacities = new int[] {30};
		dimensions = new int[] {2, 50};// test the parameters
		int testTime = 1;//test one time
		// compose the args[0] and generate the logs, we get it and 
		
		double means[] = new double[]{-2, -1, 0, 1, 2};
		double vars[] = new double[]{0.01, 0.1, 0.5, 1, 5};
		
		String fileString = args[0];
		String dataname = args[4];
		
	//	for(double mean: means)
	//	for(double var: vars)
		double mean = 0;
		double var = 1;
		/*
		for(int dimension: dimensions)
			for(int kt: kvalue) {
				args[0] =fileString+ "/"+kt+"k"+dimension+"d"+(int)mean+"m"+(int)var+"c.txt";
				args[4] = dataname+mean+"m"+var+"c";
				kmeansAlgorithm<?> runkmeans = new kmeansAlgorithm<>(args);
				for(int capacity: capacities)
					for(int scale: scales) {
						runkmeans.setDimension(dimension);
						runkmeans.setCapacity(capacity);
						runkmeans.setScale(scale);
						int[] tempk = new int[] {kt};
						runkmeans.experiments(tempk, testTime);
					}
			}
		*/
		for(int dimension: dimensions)
			for(double va: vars) {
				int kt = 10;
				if(va>=1)
					args[0] =fileString+ "/"+kt+"k"+dimension+"d"+(int)mean+"m"+(int)va+"c.txt";
				else
					args[0] =fileString+ "/"+kt+"k"+dimension+"d"+(int)mean+"m"+va+"c.txt";
				args[4] = dataname+mean+"m"+va+"c";
				kmeansAlgorithm<?> runkmeans = new kmeansAlgorithm<>(args);
				for(int capacity: capacities)
					for(int scale: scales) {
						runkmeans.setDimension(dimension);
						runkmeans.setCapacity(capacity);
						runkmeans.setScale(scale);
						int[] tempk = new int[] {kt};
						runkmeans.experiments(tempk, testTime);
					}
			}
		
		/*
		for(int dimension: dimensions)
			for(double mea: means) {
				int kt = 10;
				if(mea>=1 || mea<=-1 || mea == 0)
					args[0] =fileString+ "/"+kt+"k"+dimension+"d"+(int)mea+"m"+(int)var+"c.txt";
				else
					args[0] =fileString+ "/"+kt+"k"+dimension+"d"+mea+"m"+(int)var+"c.txt";
				args[4] = dataname+mea+"m"+var+"c";
				kmeansAlgorithm<?> runkmeans = new kmeansAlgorithm<>(args);
				for(int capacity: capacities)
					for(int scale: scales) {
						runkmeans.setDimension(dimension);
						runkmeans.setCapacity(capacity);
						runkmeans.setScale(scale);
						int[] tempk = new int[] {kt};
						runkmeans.experiments(tempk, testTime);
					}
			}*/
	//	runkmeans.staticKmeans(false, true, false);//index sign, bound, scan whole tree again
	//	runkmeans.staticKmeans(true, false, true);//test the functionality		
		//run script to draw the figures using gnuplot
		
	}
}