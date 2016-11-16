# Clustering


The file program: DCLUSTERING_161115_Manhattan.cpp 

COMPILATION:
compiles with a c++ compilor as g++:
g++ DCLUSTERING_161115_Manhattan.cpp -o

WHAT DOES?
It os thought to analize the cluster structure in a 2d rectangular grid. The input file contain the dicrete coordinates of each node to consider. The program uses the following connectivity criterion:
Two nodes i and j are connected if and only if their relation is that they are first neighbours. Although there are other ways to implement this, I used the Manhattan metric (d=abs(x2-x1)+abs(y2-y1)) and I considered connected only the nodes separated by the integer distance 1.

e(i,j) <====> d=1

For the moment the program outputs two files per configuration. One gives all the clusters and the nodes coordinates associated to each one. the second gives the sizes of the different clusters in descending order.

EXAMPLE
The file FILE_0001.s can be used as an example input file to run and check the program and understand what is doing. I also include a plot of the example configuration to understand the results by eye, FILE_0001.pdf.

