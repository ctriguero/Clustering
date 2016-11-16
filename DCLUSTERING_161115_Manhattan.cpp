///
/// Clustering analisys on a grid using Manhattan metrics
///
/// Version: 161115
///
/// Carles Triguero 2016
///


// To read files i/o
#include <iostream>
#include <fstream>
#include <sstream>
// To store data
#include <vector>
// Sort data
#include <algorithm>
#include <functional>
#include <algorithm>





using namespace std;


int main()
{

unsigned int Files_2_Read=1 ;
unsigned int File_Counter=0 ;
unsigned int i ;
unsigned int j ;
unsigned int k ;
unsigned int l ;
// Open the different files with fifferent names
// Build a name for the different frames

std::cout << endl ;
std::cout << endl ;
std::cout << "\033[1;31mGathering information about clusters associated to a set of nodes\033[0m" << endl ;
std::cout << endl ;
std::cout << endl ;
for (i=0; i<Files_2_Read; i++)
	{
	
		std::cout << endl ;
	
		stringstream ss;
		ifstream input_file;
		++File_Counter ;
		        
		if ( File_Counter <= 9 )
			{
				ss << "FILE_000" << File_Counter << ".s";
				cout << "Openning file: "<< ss.str().c_str() << endl ;
				
				input_file.open(ss.str().c_str(), ios::in ) ;
				// up to here OK
				
				
				//Number of nodes (lines of the file)
				std::string line ;    // variable string to read the line
				unsigned int nodes=0 ;
				while (std::getline(input_file, line))
				{
					nodes++;
				}			
				input_file.close() ;
				std::cout << nodes << " nodes detected in the file: "<< ss.str().c_str() << endl ;


				//Store positions in a vector
				size_t node_number = nodes ; //nodes;    
				std::vector<int> Px(node_number), Py(node_number) ;
				
				nodes=0 ;	
				input_file.open(ss.str().c_str(), ios::in ) ;
				while( std::getline( input_file, line ) )
					{
                    	double x, y ;
						std::stringstream aa(line) ;
						aa >> x >> y ;
						Px[nodes]=x ;
						Py[nodes]=y ;
						nodes++ ;						
					}
				input_file.close() ;
				
				// Clustering now using Manhattan metric d=|x2-x1|+|y2-y1|+|z3-z1|
				// Define a vector that contains the indexes of the cluster
				std::vector<int> C(node_number) ;
				// Initialize: All disconected, each node is a cluster.
				for (j=0; j<nodes; j++) C[j]=j ;
				// INFOOUT for (j=0; j<nodes; j++) cout << "Co(" << j << ")=" << C[j] << endl ;
				
				std::vector<int> auxiliar(node_number,0) ;
				for (j=0; j<nodes-1; j++)
				{
					for (k=j+1; k<nodes; k++)
					{
						if ( abs(Px[k]-Px[j]) + abs(Py[k]-Py[j]) == 1 )
						{
							// INFOOUT cout << "nodes j=" << j << " and k=" << k << " are connected" << endl ;
							// Before connect the nodes identify who also is in the group of k
							// then we will change all the group of the same class of k to j
							
							unsigned int Old_Cluster_Size=0 ;
							for (l=0; l<nodes; l++)
							{
								if ( C[k] == C[l] )
								{
									Old_Cluster_Size++ ;
									auxiliar[Old_Cluster_Size]=l ;
								}
							}
							// We connect k and their cluster mates with j
							for (l=1; l<Old_Cluster_Size+1; l++)
							{
								C[ auxiliar[l] ]=C[j] ;
							}
						}
					}
				}
				
				// Clustering done
				cout << endl ;
				// OUTINFO for (j=0; j<nodes; j++) cout << "C(" << j << ")=" << C[j] << endl ;
				
				
				
			
			
			
			
			
			
				
				
				
				// Cluster Storage
				
				// BEGIN Detect non-empty cluster names without repetition
				// Initialize auxiliar vector
				for (j=0; j<nodes; j++) auxiliar[j]=0 ;
				
				unsigned int Number_clusters=0 ;
				// Loop over possible cluster names j
				for (j=0; j<nodes; j++)
				{
					//cout << "j" << j << endl ;
					int skip=0 ;
					// Loop over all atom indexes k
					for (k=0; k<nodes; k++)
					{
					//cout << "k" << k << endl ;
						// Selects the cluster j, and counts as 1 if not empty 
						// in order to not count multiple we define skip
						// Thus v gives the non empty cluster name (from 1 to Number_clusters)
						if ( ( C[k] == j ) && ( skip == 0 ) )
						{
							Number_clusters++ ;
							auxiliar[Number_clusters]=j ;
							// Once we detect a non-empty cluster skip=1 
							// to go to another cluster name j (exit i loop).
							skip=1 ;
						}
					}
				}
				//OUTINFO for (j=1; j<Number_clusters+1; j++) cout << "aux(" << j << ")=" << auxiliar[j] << endl ;
				cout << endl ;
				cout << "Detected " << Number_clusters << " clusters" << endl ;
				//END Detect non-empty cluster names without repetition
				
				
				// Open a file for this configuration to store all the information about clustering
				ostringstream tt;
				ofstream output_file;
				tt << "Clusters_" << File_Counter << ".dat";
				cout << "Openning file: "<< tt.str().c_str() << endl ;
				output_file.open(tt.str().c_str(), ios::out | ios::trunc ) ;
				
				// Open a file for this configuration to store cluster size distribution
				ostringstream uu;
				ofstream output_sizes;
				uu << "Cluster_sizes_" << File_Counter << ".dat";
				cout << "Openning file: "<< uu.str().c_str() << endl ;
				output_sizes.open(uu.str().c_str(), ios::out | ios::trunc) ;
				output_sizes << "# Sizes of clusters in descending order" << endl ; // Header
				
				// Vector to store sizes that will allow to sort them
				std::vector<int> size(Number_clusters+1,0) ;
				
				unsigned int Cluster_Index=0 ;
				for (j=1; j<Number_clusters+1; j++)
				{
					Cluster_Index++ ;
					output_file << "Cluster number " << Cluster_Index << endl ;
					output_file << "node" << "\t" << "X" << "\t\t" << "Y" << endl ;
					unsigned int cluster_atoms=0 ;
					for (k=0; k<nodes; k++)
					{
						if ( C[k] == auxiliar[j] )
						{
							output_file << k << "\t\t" << Px[k] << "\t\t" << Py[k] << endl ;
							cluster_atoms++ ;
						}
					}
					output_file << "Cluster " << Cluster_Index << " size= " << cluster_atoms << endl ;
					output_file << endl ;
					// Sizes
					size[j]=cluster_atoms ;
				}
				output_file.close() ;
				
				// One can choose Ascending or discending ordering for the sizes of clusters
				// Ascending Sort
////////				std::sort(size.begin(), size.end()); // Ascending order works
////////				for (j=1; j<Number_clusters+1; j++) output_sizes << size[j] << endl ;
////////				cout << "Largest cluster contains " << size[Number_clusters] << " nodes." << endl ;
////////				cout << "Smallest cluster contains " << size[1] << " nodes." << endl ;
				
				// Descending Sort:
				std::vector<int> descending(Number_clusters,0) ;
				for (j=1; j<Number_clusters+1; j++) descending[j-1]=size[j] ;
				std::sort(descending.rbegin(), descending.rend());
				for (j=0; j<Number_clusters; j++) output_sizes << descending[j] << endl ;
				cout << "Largest cluster contains " << descending[0] << " nodes." << endl ;
				cout << "Smallest cluster contains " << descending[Number_clusters-1] << " nodes." << endl ;
				
				output_sizes.close() ;
			}

	}
	std::cout << endl ;
}

//////#include <string>
//////#include <vector>
//////int main(int argc, char *argv[])
//////{
//////  // check if there is more than one argument and use the second one
//////  //  (the first argument is the executable)
//////  if (argc > 1)
//////  {
//////    std::string arg1(argv[1]);
//////    // do stuff with arg1
//////  }

//////  // Or, copy all arguments into a container of strings
//////  std::vector<std::string> allArgs(argv, argv + argc);
//////}

