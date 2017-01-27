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
/// (1) Number of first neigbours
/// (2) First neighbours distances
///
/// Based on the connectivity criterion imposed by topological clustering and NOT
/// on the distance based clustering
///



//**********************************************************************************
// ATOM class and methods
//**********************************************************************************
class Atom {
	int index ; /// index = index position in the list as given by the file
	int type ;  /// type = type of atom
	int clusterindex ;
	double xpos, ypos, zpos ;
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

void Atom::SetClusterIndex (int ii) {
	clusterindex = ii ;
}

int Atom::GetClusterIndex () {
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
//		std::vector<Atom> AtomList ;
//		std::vector<Atom> NodeList ;
//		std::vector<Atom> SeaList ;
//		unsigned int NodeType ;
	public:
		std::vector<Atom> AtomList ;
		std::vector<Atom> NodeList ;
		std::vector<Atom> SeaList ;
		unsigned int NodeType ;
		AtomSet (int) ; // Constructor
		Atom GetAtom (size_t dim_AtomSetCard) ;
		Atom GetNode (size_t dim_AtomSetCard) ;
		Atom GetSea (size_t dim_AtomSetCard) ;
};

// CONSTRUCTOR Building atom objects and placing them into the AtomSet object (Atom list vector)
AtomSet::AtomSet (int a)
{
	NodeType = a ;
	
	// (0) Pass as arguments Node type and The file name (NodeType, file.xyz) to the class AtomSet
	// (1) Autodetect the number of atoms from the file


	
	// Input file
	std::ifstream InputFile ;
	InputFile.open("test.xyz", ios::in) ;


	std::string line ;			// variable string to read the line
	unsigned int element ;
	double x, y, z ;

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
			Atom_object.SetClusterIndex(NodeIndex) ;
			NodeList.push_back (Atom_object) ;
			NodeIndex++ ;
		}
		else
		{
			Atom_object.SetClusterIndex(SeaIndex) ;
			SeaList.push_back (Atom_object) ;
			SeaIndex++ ;
		}
	}
	InputFile.close() ;

	size_t dim_AtomSetCard = AtomList.size () ;
	size_t dim_NodeSetCard = NodeList.size () ;
	size_t dim_SeaSetCard = SeaList.size () ;

	std::cout << BOLDYELLOW << "-> Detected atoms in " << "test.xyz" << RESET << endl ;
	std::cout << "\t" << BLUE << "Total atoms: " << BOLDBLUE << dim_AtomSetCard << RESET << endl ;
	std::cout << "\t" << BLUE << "Node atoms: " << BOLDBLUE << dim_NodeSetCard << RESET << endl ;
	std::cout << "\t" << BLUE << "Sea atoms: " << BOLDBLUE << dim_SeaSetCard << RESET << endl ;

}

// defining AtomSet class methods
//   method get the Atom in a location of the AtomList
Atom AtomSet::GetAtom(size_t dim_AtomSetCard)
{
	return AtomList[dim_AtomSetCard];
}

Atom AtomSet::GetNode(size_t dim_NodeSetCard)
{
	return AtomList[dim_NodeSetCard];
}

Atom AtomSet::GetSea(size_t dim_SeaSetCard)
{
	return AtomList[dim_SeaSetCard];
}


	

//**********************************************************************************
// main function
//**********************************************************************************
int main( int argc, const char* argv[] )
{
	cout << endl ;
	cout << BOLDYELLOW << "-----------------------------" << RESET << endl ;
	cout << BOLDYELLOW << " TOPOLOGICAL CLUSTERING 2017   " << RESET << endl ;
	cout << BOLDYELLOW << "-----------------------------" << RESET << endl ;
	cout << endl ;

	if ( argc == 1 )
	{

		cout << " " << endl ;
		cout << "\t Wrong usage. You should execute:" << endl ;
		cout << BOLDRED << "\t " << argv[0] << " [File to analyze.sd]" << RESET << endl ;
		cout << " " << endl ;
		return (0) ;
	}

	// BEGIN pipe in file to analyze
	std::ifstream InputFile;
	cout << BOLDYELLOW << "-> Processing file " << argv[1] << RESET << endl ;

	// Parameters:
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------	
	unsigned int NodeType = 1 ;	// Tells the program what type of atoms will be clustered all other will be consider as sea atoms.
	double R1 = 2.5 ;		// Radius of type 1 atoms.
	double R2 = 1.7 ;		// Radius of type 2 atoms.
	double dT = 10.0 ;		// Distance threshold.
	double hc = 3.0*R2 ;		// Cylinder height (Bin cylinder) //in sea atoms diameter units.
	double rc = 1.5*R2 ;		// Cylinder radius //in sea atoms diameter units.
	double rhoT = 0.0001 ;		// Connectivity transition density threshold.
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
	
	ofstream Frequencies ;
	Frequencies.open ("Frequencies.dat", ios::out | ios::trunc); // app=append
	ofstream Distances ;
	Distances.open ("Distances.dat", ios::out | ios::trunc); // app=append
	
	// Define the object set of atoms called AtomSet_object with Natoms It is automatically initialized with the data in the file that we defined in the constructor
	AtomSet AtomSet_object(NodeType) ;
	size_t NodeAtoms = AtomSet_object.NodeList.size() ; // Total number of atoms in the file (Node atoms + Sea atoms).

	for ( size_t i = 0; i < NodeAtoms-1; ++i )
	{
		ifriends = 0 ;
		ix = AtomSet_object.GetNode(i).GetX() ;
		iy = AtomSet_object.GetNode(i).GetY() ;
		iz = AtomSet_object.GetNode(i).GetZ() ;

		for ( size_t j = i+1; j < NodeAtoms; ++j )
		{
			std::cout << "i= " << i << " j= " << j << endl ;
			jx = AtomSet_object.GetNode(j).GetX() ;
			jy = AtomSet_object.GetNode(j).GetY() ;
			jz = AtomSet_object.GetNode(j).GetZ() ;

			ijx = jx - ix ;
			ijy = jy - iy ;
			ijz = jz - iz ;

			modij = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;

			// Distance restriction OPTIMIZATION
			if ( modij > dT ) { continue ; }

			// Normalized components 
			ijx = ijx/modij ;
			ijy = ijy/modij ;
			ijz = ijz/modij ;
			
			Atom_count = 0 ;
				
			// (1) If the distance between atoms to clusterize is too small to fit a sea atom then they are connected OPTIMIZATION
			if ( modij < 2*R1 + 2*R2 ) { continue ; }
			// (2) If modij < bc then we just define a single partition: 1 bin cylinder (otherwise nbins=0 and dv-> Infinity) OPTIMIZATION
			if ( modij < hc )
			{
				// Definition of the number of bins
				nbins = 1 ;
				Vol = pi*rc*rc*modij ;
			}
			else
			{
				// Definition of the number of bins
				nbins = floor(modij/hc) ;
				std::cout << "nbins= " << nbins << endl ;
				Vol = pi*rc*rc*modij/float(nbins) ;
			}
////////////					size_t dim_nbins = nbins ;
////////////					std::vector<int> bin ;
////////////					bin.resize(dim_nbins) ;

////////////					for ( size_t k = 0; k < Natoms; ++k )
////////////					{
////////////						unsigned int typek = AtomSet_object.GetAtom(k).GetType() ;
////////////						if ( typek != NodeType )
////////////						{						
////////////							ikx = AtomSet_object.GetAtom(k).GetX() - ix ;
////////////							iky = AtomSet_object.GetAtom(k).GetY() - iy ;	
////////////							ikz = AtomSet_object.GetAtom(k).GetZ() - iz ;	
////////////							modik = sqrt (ikx*ikx + iky*iky + ikz*ikz) ;						
////////////							ikpar = ijx*ikx + ijy*iky + ijz*ikz ;
////////////							ikver = sqrt (modik*modik - ikpar*ikpar) ;						

////////////							// Check if k atoms are inside the cylinder
////////////							if ( ( ikpar > 0.0 ) && ( ikpar < modij ) && ( ikver < rc ))
////////////							{
////////////								// Binning along the cyclinder
////////////								n =  int( nbins*ikpar/modij + 1.0 ) ;
////////////								++bin[n] ;
////////////							}
////////////						}// If k is different from i or j			
////////////					}// Sum over all atoms

////////////					// Connectivity decision based on the local density in the bins: mini cylinders (hc,rc)
////////////					// If any of the bins of the cylinder has a larger density then the atoms are not connected
////////////					for ( size_t k = 0; k < nbins; ++k )
////////////					{
////////////						rho = float(bin[k])/Vol ;
////////////						if ( rho >= rhoT ) { goto ijLoopEnd ; }
////////////					}


////////////					// Detecting all atoms already belonging to the cluster c(i) or to the same as i
////////////					std::vector<int> auxiliar(dim_nbins+1,0) ;
////////////	
////////////					unsigned int OldClusterSize=0 ;
////////////					for ( int m=0; m<Natoms; m++)
////////////					{
////////////						unsigned int Cj = AtomSet_object.GetAtom.GetClusterIndex[j] ;
////////////						unsigned int Ck = AtomSet_object.GetAtom.GetClusterIndex[k] ;
////////////						unsigned int Cm = AtomSet_object.GetAtom.GetClusterIndex[m] ;

////////////						if ( Ck == Cm )
////////////						{
////////////							OldClusterSize++ ;
////////////							auxiliar[OldClusterSize]=m ;
////////////						}
////////////					}
////////////					for ( int m=1; m<OldClusterSize+1; m++)
////////////					{
////////////						AtomSet_object.GetAtom.SetClusterIndex[auxiliar[m]] = Cj ;
////////////					} // END Detecting all atoms already belonging to the cluster c(i) or to the same as i
				//}// i and j are nodes

			//ijLoopEnd:

			} // End j
			Frequencies << ifriends << endl ; // Store in file Frequencies.dat	
			std::cout << "Atom " << i+1 << " has " << ifriends << " first neighours" << endl ;		
		//}	
	} // End i
	
	Frequencies.close() ;
	Distances.close() ;

	
	
	//cout << "\033[1m Results: \033[0m" << endl ;
	cout << endl;
	//cout << "Mean number of first nighbours = \033[1m" << "\033[0m" << endl ;
	cout << endl;
	//cout << "Standard deviation of the number of first nighbours = \033[1m" << "\033[0m" << endl ;		
	
	// End counting time 	
 	int stop_s=clock();
 	
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

