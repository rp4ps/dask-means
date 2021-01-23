# pick-means
Fast and Memory-Efficient K-Means on Large-Scale Point Clouds

## Overview
This repo holds the source code and scripts for reproducing the key experiments of fast large-scale point cloud clustering. 

## Build
mvn clean package

## Datasets
Download the six datasets we used in our paper 
https://drive.google.com/drive/folders/1CQCtkGGiKqqWacMN0PyyAHg2rDvkyqg2?usp=sharing

## Run

### Bird

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiencyLight ./dataset/pointcloud/birdfountain_station1_xyz_intensity_rgb.txt 10 40000000 a birdfountain 1 3

### Paris

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiencyLight ./dataset/pointcloud/Paris.csv 10 10000000 a paris  0 2

### L001

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiencyLight ./dataset/pointcloud/Toronto_3D/L1.csv 10 21000000 a toronto1  0 2

### Cathedral

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiencyLight ./dataset/pointcloud/stgallencathedral_station1.txt 10 31000000 a cathedral1  0 2

### Lille

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiencyLight ./dataset/pointcloud/Lille0.csv 10 10000000 a lille0  0 2

### L003

java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiencyLight ./dataset/pointcloud/Toronto_3D/L3.csv 10 39000000 a toronto3  0 2
