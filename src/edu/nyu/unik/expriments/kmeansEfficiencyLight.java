package edu.nyu.unik.expriments;

import java.io.IOException;

import au.edu.rmit.trajectory.clustering.kmeans.kmeansAlgorithm;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;
import java_cup.internal_error;

public class kmeansEfficiencyLight {

	public kmeansEfficiencyLight() {
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) throws IOException, KeySizeException, KeyDuplicateException {
		// TODO Auto-generated method stub
		int[] kvalue = new int[]{10, 50, 100, 200, 400, 600, 800, 1000};//10, 100, 1000
		int[] scales = new int[] {100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000}; //1000, 5000, 10000, 50000, 100000, 500000, 1000000
		int[] capacities = new int[] {10, 20, 30, 40, 50, 60};//capacity
		int[] dimensions = new int[] {2, 10, 20, 30, 40, 50};
		
	//	dimensions = new int[] {7};
		kvalue = new int[]{10,100,1000};//10, 100, 1000
	//	kvalue = new int[]{10,20,30,40,50};//10, 100, 1000
		kvalue = new int[]{100, 1000};//10, 100, 1000
		kvalue = new int[]{5000, 10000, 20000, 30000};// test k, change all the short to int if bigger than 32767.
	//	kvalue = new int[]{10, 100, 400, 700, 1000};//10, 100, 1000
		dimensions = new int[] {2, 10, 100, 200, 300, 400, 500}; //for 
		dimensions = new int[] {2};
		scales = new int[] {10000};
	//	scales = new int[] {1000, 5000, 10000, 50000, 100000, 500000, 1000000};
		capacities = new int[] {100, 150, 200};
	//	dimensions = new int[] {2};// test the parameters
		int testTime = 1;//test one time
		kmeansAlgorithm<?> runkmeans = new kmeansAlgorithm<>(args);
		kvalue = new int[]{10000};//100, 1000, 
	//	runkmeans.lightWeightKmeans(false, true, false);//index sign, bound, scan whole tree again
		runkmeans.experiments(kvalue, testTime);// run the configuration that specified, test various k
	//	runkmeans.staticKmeans(true, false, true);//test the functionality		
		//run script to draw the figures using gnuplot
		
		kvalue = new int[]{10000};
		
		// test the fast algorithm itself for various parameters
		/*
		int scaleData = runkmeans.getTraNum();
		for (int capacity : capacities)//test various parameters
		{
			runkmeans.setCapacity(capacity);
			runkmeans.experiments(kvalue, testTime);
		}
		
		runkmeans.setCapacity(30);
		scales = new int[] {(int) (0.01*scaleData), (int) (0.05*scaleData), (int) (0.25*scaleData)};
		for (int scale : scales) {//test various dataset scale
			if (scale > scaleData)
				continue;
			runkmeans.setScale(scale);
			runkmeans.experiments(kvalue, testTime);
		}*/
		
	}
}