// classes example
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <unistd.h>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <ctime>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

using namespace std;


///
/// This program calculates the distribution of: 
///
/// (1) Number of first neighbours
/// (2) First neighbours distances
///
/// Based on the connectivity criterion imposed by topological clustering and NOT
/// on the distance based clustering
///



//**********************************************************************************
// ATOM class and methods
//**********************************************************************************
class Atom {
	private:
		int index ;			/// index = index position in the list as given by the file
		int type ;			/// type = type of atom
		int clusterindex ;		/// clusterindex = index telling to which cluster belongs the atom
		double xpos, ypos, zpos ;	/// 3 cartesian coordinate positions of the atoms
	public:
		void SetIndex(int) ;
		int GetIndex() ;
		void SetType(int) ;
		int GetType() ;
		void SetClusterIndex(int) ;
		int GetClusterIndex() ;		
		void SetX(double) ;
		double GetX() ;
		void SetY(double) ;
		double GetY() ;
		void SetZ(double) ;
		double GetZ() ;
};

void Atom::SetIndex (int ii) {
	index = ii ;
}

int Atom::GetIndex () {
	return index ;
}

void Atom::SetType (int tt) {
	type = tt ;
}

int Atom::GetType () {
	return type ;
}

void Atom::SetClusterIndex (int ii)
{
	clusterindex = ii ;
}

int Atom::GetClusterIndex () 
{
	return clusterindex ;
}

void Atom::SetX (double xx) {
	xpos = xx ;
}

double Atom::GetX () {
	return xpos ;
}

void Atom::SetY (double yy) {
	ypos = yy ;
}

double Atom::GetY () {
	return ypos ;
}

void Atom::SetZ (double zz) {
	zpos = zz ;
}

double Atom::GetZ () {
	return zpos ;
}



//**********************************************************************************
// AtomSet class and methods
//**********************************************************************************
class AtomSet
{
	private:
		unsigned int SeaType ;
	public:
		std::vector<Atom> AtomList ;
		std::vector<Atom> NodeList ;
		std::vector<Atom> SeaList ;
		unsigned int NodeType ;
		std::string  FileName ;
		AtomSet (int, std::string) ; // Constructor
		Atom GetAtom (size_t) ;
		Atom GetNode (size_t) ;
		Atom GetSea (size_t) ;
		void ChangeNodeClusterIndex (size_t, int) ;
};

// CONSTRUCTOR Building atom objects and placing them into the AtomSet object (Atom list vector)
AtomSet::AtomSet (int a, std::string b)
{
	NodeType = a ;
	FileName = b ;
	// (0) Pass as arguments Node type and The file name (NodeType, file.xyz) to the class AtomSet
	// (1) Autodetect the number of atoms from the file
	
	// BEGIN pipe in file to analyze
	cout << YELLOW << "    -> Processing file: " << BOLDYELLOW << FileName << RESET << endl ;
	
	std::ifstream InputFile ;
	InputFile.open( FileName.c_str(), ios::in );


	std::string line ;			// variable string to read the line
	unsigned int element ;
	double x, y, z ;

	// Initialize counters
	unsigned int AtomIndex = 0 ;
	unsigned int NodeIndex = 0 ;
	unsigned int SeaIndex = 0 ;
	
	while (getline(InputFile, line))
	{
		//std::getline(file, line) ;
		std::stringstream aa(line) ;
		aa >> element >> x >> y >> z ;
		// creating an Atom object that will contain the type and position
		Atom Atom_object = Atom() ;
		// Assingning the atom [index, type, cluster index, and spatial coordinates] to the Atom object
		Atom_object.SetIndex(AtomIndex) ;
		AtomIndex++ ;
		Atom_object.SetType(element) ;
		Atom_object.SetX(x) ;
		Atom_object.SetY(y) ;
		Atom_object.SetZ(z) ;
		// Assigning this atom object to the i position of the atom list
		AtomList.push_back (Atom_object) ;
		// Assigning this atom object to the node and sea list
		if ( element == NodeType )
		{
			Atom_object.SetClusterIndex (NodeIndex) ;
			NodeList.push_back (Atom_object) ;
			//std::cout << "c index node= " << NodeIndex << endl ;
			NodeIndex++ ;
		}
		else
		{
			Atom_object.SetClusterIndex (SeaIndex) ;
			SeaList.push_back (Atom_object) ;
			//std::cout << "c index sea= " << SeaIndex << endl ;
			SeaIndex++ ;
		}
	}
	InputFile.close() ;
	
	if ( SeaIndex == 0 )
	{
		cout << " " << endl ;
		std::cout << BOLDRED <<"    -> ERROR: No atoms found in the file " << FileName << ". ABORTING" << RESET << std::endl ;
		cout << " " << endl ;
		abort () ;
	}

	size_t dim_AtomSetCard = AtomList.size () ;
	size_t dim_NodeSetCard = NodeList.size () ;
	size_t dim_SeaSetCard = SeaList.size () ;
	
	
	if ( NodeType == 1 ) { SeaType = 2 ; }
	if ( NodeType == 2 ) { SeaType = 1 ; }

	std::cout << YELLOW << "    -> Summary of detected atoms in file " << BOLDYELLOW << FileName << RESET << endl ;
	std::cout << "\t" << BLUE << "    Considering atoms with label " << NodeType << " as node atoms" << RESET << endl ;
	std::cout << "\t" << BLUE << "    Considering atoms with label " << SeaType << " as sea atoms" << RESET << endl ;
	std::cout << "\t" << BLUE << "    * Total number of atoms: " << BOLDBLUE << dim_AtomSetCard << RESET << endl ;
	std::cout << "\t" << BLUE << "    * Node atoms: " << BOLDBLUE << dim_NodeSetCard << RESET << endl ;
	std::cout << "\t" << BLUE << "    * Sea atoms: " << BOLDBLUE << dim_SeaSetCard << RESET << endl ;
}

Atom AtomSet::GetAtom(size_t dim_AtomSetCard)
{
	return AtomList[dim_AtomSetCard] ;
}

Atom AtomSet::GetNode(size_t dim_NodeSetCard)
{
	return NodeList[dim_NodeSetCard] ;
}

Atom AtomSet::GetSea(size_t dim_SeaSetCard)
{
	return SeaList[dim_SeaSetCard] ;
}

// This operator acts over an AtomSet object ( AtomSet_object.ChangeNodeClusterIndex (0,4): Sets to Node cluster index to 4 for node 0 )
void AtomSet::ChangeNodeClusterIndex(size_t dim_NodeSetCard, int value)
{
	NodeList[dim_NodeSetCard].SetClusterIndex (value) ;
}





	

//**********************************************************************************
// main function
//**********************************************************************************
int main( int argc, const char* argv[] )
{
	std::cout << endl ;
	std::cout << BOLDYELLOW << "    _________________________________________________" << std::endl ;
	std::cout << BOLDYELLOW << "            _=_                                      " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "          q(-_-)p                                    " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "          '_) (_`         TOPOLOGICAL CLUSTERING     " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "          /__/  \\         Carles Triguero 2017      " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "        _(<_   / )_                                  " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "       (__\\_\\_|_/__)                               " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "    _________________________________________________" << RESET << std::endl ;
	std::cout << endl ;

	// Check the number of arguments (1 means only the name of the program => No argument)
	if ( argc == 1 )
	{
		cout << "\t Wrong usage. You should execute:" << endl ;
		cout << BOLDRED << "\t " << argv[0] << " [File to analyze.sd]" << RESET << endl ;
		cout << " " << endl ;
		return (0) ;
	}
	
	// Check if the file exists
	ifstream file(argv[1]) ;
	if (!file)
	{
    	std::cout << BOLDRED <<"    -> ERROR: File test.xyz does not exist. ABORTING" << RESET << std::endl ;
    	std::cout << std::endl ;
		abort () ;
	}
	file.close() ;
	
	std::cout << YELLOW << "    -> File " << BOLDYELLOW << argv[1] << RESET << YELLOW << " found" << RESET << std::endl ;

	// Parameters:
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------	
	unsigned int NodeType = 1 ;	// Tells the program what type of atoms will be clustered all other will be consider as sea atoms.
	double R1 = 2.5 ;		// Radius of type 1 atoms.
	double R2 = 1.7 ;		// Radius of type 2 atoms.
	double dT = 10.0 ;		// Distance threshold.
	double hc = 6.0*R2 ;		// Cylinder height (Bin cylinder) //in sea atoms diameter units.
	double rc = 1.5*R2 ;		// Cylinder radius //in sea atoms diameter units.
	double rhoT = 0.01 ;		// Connectivity transition density threshold. 0 => all disconnected 100 => all connected  sea=0.0181
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------		
	//hc = hc*2*R2
	
	const double pi = std::atan(1.0)*4 ;
	double rho ;
	
	double ix, iy, iz, jx, jy, jz, ijx, ijy, ijz, ikx, iky, ikz, modij, Vol, ikpar, modik, ikver, den ;
	unsigned int Atom_count = 0 ;
	unsigned int memo = 0 ;
	unsigned int count = 0 ;
	unsigned int ifriends = 0 ;
	unsigned int n ;
	unsigned int nbins ;
	
	// Count execution time
	int start_s=clock();
	
	// Files to store data
	ofstream Densities ;
	Densities.open ("Densities.dat", ios::out | ios::trunc); // All computed densities inside cylinders
	ofstream Inneratoms ;
	Inneratoms.open ("Inneratoms.dat", ios::out | ios::trunc); // Number of sea atoms inside cylinders
	
	std::string File(argv[1]); // set up a stringstream variable named convert, initialized with the input from argv[1]
	
	// Define the object set of atoms called AtomSet_object with Natoms It is automatically initialized with the data in the file that we defined in the constructor
	AtomSet AtomSet_object(NodeType, File) ;
	
	size_t NodeAtoms = AtomSet_object.NodeList.size() ; // Node atoms
	size_t SeaAtoms = AtomSet_object.SeaList.size() ; // Sea atoms

	std::cout << YELLOW << "    -> Density threshold set to " << BOLDYELLOW << rhoT << RESET << std::endl ;
	std::cout << "\t" << BOLDBLUE   << "    * It is very useful to compute rho_sea before setting rhoT" << RESET << std::endl ;
	std::cout << "\t" << BLUE   << "    * rhoT should be bounded: " << BOLDBLUE << "[rho_sea > rho_T >= 0]" << RESET << std::endl ;
	std::cout << "\t" << BLUE   << "    * rhoT=0 => Strict connectivity (few connexions or zero)" << RESET << std::endl ;
	std::cout << "\t" << BLUE   << "    * rhoT >> rho_sea => All Connected" << RESET << std::endl ;



	// Begin Main loop over node i atoms
	for ( size_t i = 0; i < NodeAtoms-1; ++i )
	{
		ix = AtomSet_object.GetNode(i).GetX() ;
		iy = AtomSet_object.GetNode(i).GetY() ;
		iz = AtomSet_object.GetNode(i).GetZ() ;

		// Begin Main loop over node j atoms
		for ( size_t j = i+1; j < NodeAtoms; ++j )
		{
			jx = AtomSet_object.GetNode(j).GetX() ;
			jy = AtomSet_object.GetNode(j).GetY() ;
			jz = AtomSet_object.GetNode(j).GetZ() ;

			ijx = jx - ix ;
			ijy = jy - iy ;
			ijz = jz - iz ;

			modij = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;


			// Normalized components 
			ijx = ijx/modij ;
			ijy = ijy/modij ;
			ijz = ijz/modij ;
			
			Atom_count = 0 ;
				
			// (1) If the distance between atoms to clusterize is too small to fit a sea atom then they are connected OPTIMIZATION
//			if ( modij < 2*R1 + 2*R2 ) { continue ; }

			// (2) If modij < hc then we just define a single partition: 1 bin cylinder (otherwise nbins=0 and dv-> Infinity) OPTIMIZATION
//			std::cout << "modij= " << modij << std::endl ;
//			std::cout << "hc = " << hc << std::endl ;
			if ( modij < hc )
			{
				// Definition of the number of bins
				nbins = 1 ;
//				std::cout << "nbins= " << nbins << std::endl ;
				Vol = pi*rc*rc*modij ;
			}
			else
			{
				// Definition of the number of bins
				nbins = floor(modij/hc) ;
//				std::cout << "nbins= " << nbins << std::endl ;
				Vol = pi*rc*rc*modij/float(nbins) ;
			}

			size_t dim_nbins = nbins ;
			std::vector<int> bin(dim_nbins+1, 0) ;

			// Loop over sea atoms
			for ( size_t k = 0; k < SeaAtoms; ++k )
			{
//				std::cout << "Inside loop k " << k << std::endl ;
		
				ikx = AtomSet_object.GetSea(k).GetX() - ix ;
				iky = AtomSet_object.GetSea(k).GetY() - iy ;
				ikz = AtomSet_object.GetSea(k).GetZ() - iz ;
				modik = sqrt (ikx*ikx + iky*iky + ikz*ikz) ;						
				ikpar = ijx*ikx + ijy*iky + ijz*ikz ;
				ikver = sqrt (modik*modik - ikpar*ikpar) ;					

				// Check if k atoms are inside the cylinder
				if ( ( ikpar > 0.0 ) && ( ikpar < modij ) && ( ikver < rc ) )
				{
					// Binning along the cyclinder
					n =  int( nbins*ikpar/modij + 1.0 ) ;
					bin[n]++ ;
				}
			}// End k loop


			// Connectivity decision based on the local density in the bins: mini cylinders (hc,rc)
			// If any of the bins of the cylinder has a larger density threshold then the atoms are not connected

			unsigned int connexion = 1 ;

			for ( size_t k = 0; k < nbins; k++ )
			{
				rho = float(bin[k])/Vol ;
				Inneratoms << bin[k] << std::endl ;	// Store the atoms inside each section to analize the system via the distributions
				Densities << rho << std::endl ;		// Store the computed densities to analize the system via the distributions
				if ( rho >= rhoT ) { connexion = 0 ; }
			}

			if ( connexion == 1 )
			{
//				std::cout << "ATOMS " << i << " AND " << j << " ARE CONNECTED" << std::endl ;


				// BEGIN Connecting atoms: Assingning the same cluster index.
	
				// STEP -1- Detecting all atoms already belonging to the cluster C(i) or to the same as i
				std::vector<int> iCluster ;
	
				unsigned int Ci = AtomSet_object.GetNode(i).GetClusterIndex() ;
	
				for ( int m=0; m<NodeAtoms; m++)
				{
					unsigned int Cm = AtomSet_object.GetNode(m).GetClusterIndex() ;
					if ( Ci == Cm ) { iCluster.push_back(m) ; } // iCluster is a list with all node atoms with same cluster index as i node
				}
	
				// STEP -2- Connecting cluster associated to i with cluster associated to j
				unsigned int Cj = AtomSet_object.GetNode(j).GetClusterIndex() ;
				for ( int m=0; m<iCluster.size(); m++) { AtomSet_object.ChangeNodeClusterIndex (iCluster[m],Cj) ; }  //OldCluster i Size = iCluster.size() ;

				// END Connecting atoms: Assingning the same cluster index.
			}


		} // End loop j	
	} // End loop i

//	for ( int i=0; i<NodeAtoms; i++) { std::cout << "Index "<< i << " cluster index " << AtomSet_object.GetNode(i).GetClusterIndex () << std::endl ; }
	
	Densities.close() ;
	Inneratoms.close() ;

	std::cout << YELLOW << "    -> Stored all computed densities in " << BOLDYELLOW << " Densities.dat" << RESET << std::endl ;
	std::cout << YELLOW << "    -> Stored all computed inner atoms in " << BOLDYELLOW << " Inneratoms.dat" << RESET << std::endl ;




////	// BEGIN Clustering recognition
////	for ( int j=0; j<X.size(); j++) auxiliar[j]=0 ;

////	unsigned int ClusterNumber = 0 ;
////	for ( int j=0; j<X.size(); j++)
////	{
////	int skip=0 ;
////	for ( int k=0; k<X.size(); k++)
////	{
////	if ( ( C[k] == j ) && ( skip == 0 ) )
////	{
////	auxiliar[Number_clusters]=j ;
////	Number_clusters++ ;
////	skip=1 ;
////	}
////	}
////	}
////	cout << "\t-> Detected " << BOLDRED << Number_clusters << " clusters" << RESET << endl ;
////	// END Clustering recognition






	
	// End counting time 	
 	int stop_s=clock();
 	
	// Execution time
 	cout << endl ;
	cout << "Execution time: \033[1m" << ( stop_s - start_s )/double(CLOCKS_PER_SEC) << "\033[0m seconds" << endl ;
	cout << endl ;
	
	return (0) ;
}




// Add construction of Histograms?

// count:
////////http://www.cplusplus.com/reference/algorithm/count/
//////// count algorithm example
//////#include <iostream>     // std::cout
//////#include <algorithm>    // std::count
//////#include <vector>       // std::vector

//////int main () {
//////  // counting elements in array:
//////  int myints[] = {10,20,30,30,20,10,10,20};   // 8 elements
//////  int mycount = std::count (myints, myints+8, 10);
//////  std::cout << "10 appears " << mycount << " times.\n";

//////  // counting elements in container:
//////  std::vector<int> myvector (myints, myints+8);
//////  mycount = std::count (myvector.begin(), myvector.end(), 20);
//////  std::cout << "20 appears " << mycount  << " times.\n";

//////  return 0;
//////}

