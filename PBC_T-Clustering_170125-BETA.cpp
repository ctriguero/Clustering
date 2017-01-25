// 
// Topological Clustering with periodic boundary conditions
// Carles Triguero 2017
//
// Versiom: 170119
//
// Compilation: g++ -std=c++11 PBC_T-Clustering_170119.cpp
//

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
#include <algorithm>    // std::min
#include <string>

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
/// This program determines the connectivity of the nodes attending to the topological connection criterion
///
/// (1) The number of sea atoms inside the cylinder defined by the nodes
///     Usual procedure when calculating without periodic boundary conditions
///
/// (2) Determines among the 8 possibilities the minimal distance image of the second node of each pair (i,j)
///	min( d(i,IL(j)), d(i,IR(j)), d(i,IT(j)), d(i,IB(j)), d(i,ILB(j)), d(i,ILT(j)), d(i,IRB(j)), d(i,IRT(j)) )
/// (3) The number of sea atoms inside the cylinder defined by the node i and the minimal distance image of node j
///
///  *  The implementation of partition cylinder and comparison with threshold density for both possibilities needs to be implemented
///  *  Select among the two possibilities of connection whether the node is connected or not
///




//**********************************************************************************
// ATOM class and methods
//**********************************************************************************

class Atom
{
	int index ; /// index = index position in the list as given by the file
	int type ;  /// type = type of atom
	double xpos, ypos, zpos ;
	public:
	void SetType(int) ;
	int GetType() ;
	void SetIndex(int) ;
	int GetIndex() ;
	void SetX(double) ;
	double GetX() ;
	void SetY(double) ;
	double GetY() ;
	void SetZ(double) ;
	double GetZ() ;
} ;




void Atom::SetType (int tt) 
{
	type = tt ;
}

int Atom::GetType () 
{
	return type ;
}

void Atom::SetIndex (int ii) 
{
	index = ii ;
}

int Atom::GetIndex () 
{
	return index ;
}

void Atom::SetX (double xx) 
{
	xpos = xx ;
}

double Atom::GetX () 
{
	return xpos ;
}

void Atom::SetY (double yy) 
{
	ypos = yy ;
}

double Atom::GetY () 
{
	return ypos ;
}

void Atom::SetZ (double zz) 
{
	zpos = zz ;
}

double Atom::GetZ () 
{
	return zpos ;
}




//**********************************************************************************
// AtomSet class and methods
//**********************************************************************************

class AtomSet
{
	private:
	size_t AtomSetCard ;
	std::vector<Atom> AtomList ;
	double Lx = 10.0, Ly = 10.0, Lz = 10.0 ; // Size of the box.
	public:
	AtomSet(size_t dim_AtomSetCard) ; // Constructor
	Atom GetAtom(size_t dim_AtomSetCard) ;
} ;




// CONSTRUCTOR Building atom objects and placing them into the AtomSet object (Atom list vector)
AtomSet::AtomSet(size_t dim_AtomSetCard ) : AtomSetCard(dim_AtomSetCard), AtomList()
{
	AtomList.resize(AtomSetCard);

	std::ifstream file( "test_3.xyz" ) ; // input file
	std::string line ;    // variable string to read the line
	unsigned int element ;
	double x, y, z ;

	for ( size_t k = 0; k < AtomSetCard; ++k ) 
	{
		std::getline(file, line) ;
		std::stringstream aa(line) ;
		aa >> element >> x >> y >> z ;
		// creating an Atom object that will contain the type and position
		Atom Atom_object = Atom() ;
		// Assingning the atom type to the Atom object
		Atom_object.SetType(element) ;
		Atom_object.SetIndex(k) ;
		Atom_object.SetX(x) ;
		Atom_object.SetY(y) ;
		Atom_object.SetZ(z) ;
		// assigning this atom object to the i position of the atom list
		AtomList[k] = Atom_object ;
	}
	file.close() ;
}

//   method get the Atom in a location of the AtomList
Atom AtomSet::GetAtom(size_t dim_AtomSetCard) 
{
	return AtomList[dim_AtomSetCard];
}




//**********************************************************************************	std::vector<int> Xp ;
// AtomSet class and methods								Xp.push_back (xpos) ;
//**********************************************************************************

class SeaAtomSet
{
	private:
	size_t SeaAtomSetCard ;
	std::vector<Atom> SeaAtomList ;
	double Lx = 10.0, Ly = 10.0, Lz = 10.0 ; // Size of the box.
	public:
	SeaAtomSet(size_t dim_SeaAtomSetCard) ; // Constructor
	Atom GetSeaAtom(size_t dim_SeaAtomSetCard);
} ;

// Count only sea atoms: discard nodes!
// Superatom set
// CONSTRUCTOR Building atom objects and placing them into the AtomSet object (Atom list vector)
SeaAtomSet::SeaAtomSet(size_t dim_SeaAtomSetCard) : SeaAtomSetCard(dim_SeaAtomSetCard), SeaAtomList()
{
	SeaAtomList.resize(9*SeaAtomSetCard-18);

	std::ifstream file( "test_3.xyz" ) ; // input file
	std::string line ;    // variable string to read the line
	unsigned int element ;
	double x, y, z ;

	unsigned int kk = 0 ;

	for ( size_t k = 0; k < SeaAtomSetCard; ++k ) 
	{
		std::getline(file, line) ;
		std::stringstream aa(line) ;
		aa >> element >> x >> y >> z ;

		if ( element == 2 )
		{
			
			// creating an Atom object that will contain the type and position
			Atom Atom_object = Atom() ;
			Atom_object.SetType(element) ;
			Atom_object.SetIndex(kk) ;
			Atom_object.SetX(x) ;
			Atom_object.SetY(y) ;
			Atom_object.SetZ(z) ;
			SeaAtomList[kk] = Atom_object ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: IL()
			kk++ ;
			Atom Atom_objectL = Atom() ;
			Atom_objectL.SetType(element) ;
			Atom_objectL.SetIndex(kk) ;
			Atom_objectL.SetX(x-Lx) ;
			Atom_objectL.SetY(y) ;
			Atom_objectL.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectL ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: IR()
			kk++ ;
			Atom Atom_objectR = Atom() ;
			Atom_objectR.SetType(element) ;
			Atom_objectR.SetIndex(kk) ;
			Atom_objectR.SetX(x+Lx) ;
			Atom_objectR.SetY(y) ;
			Atom_objectR.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectR ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: IB()
			kk++ ;
			Atom Atom_objectB = Atom() ;
			Atom_objectB.SetType(element) ;
			Atom_objectB.SetIndex(kk) ;
			Atom_objectB.SetX(x) ;
			Atom_objectB.SetY(y-Ly) ;
			Atom_objectB.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectB ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: IT()
			kk++ ;
			Atom Atom_objectT = Atom() ;
			Atom_objectT.SetType(element) ;
			Atom_objectT.SetIndex(kk) ;
			Atom_objectT.SetX(x) ;
			Atom_objectT.SetY(y+Ly) ;
			Atom_objectT.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectT ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: ILB()
			kk++ ;
			Atom Atom_objectLB = Atom() ;
			Atom_objectLB.SetType(element) ;
			Atom_objectLB.SetIndex(kk) ;
			Atom_objectLB.SetX(x-Lx) ;
			Atom_objectLB.SetY(y-Ly) ;
			Atom_objectLB.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectLB ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: ILT()
			kk++ ;
			Atom Atom_objectLT = Atom() ;
			Atom_objectLT.SetType(element) ;
			Atom_objectLT.SetIndex(kk) ;
			Atom_objectLT.SetX(x-Lx) ;
			Atom_objectLT.SetY(y+Ly) ;
			Atom_objectLT.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectLT ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: IRB()
			kk++ ;
			Atom Atom_objectRB = Atom() ;
			Atom_objectRB.SetType(element) ;
			Atom_objectRB.SetIndex(kk) ;
			Atom_objectRB.SetX(x+Lx) ;
			Atom_objectRB.SetY(y-Ly) ;
			Atom_objectRB.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectRB ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Adding Image atoms: IRT()
			kk++ ;
			Atom Atom_objectRT = Atom() ;
			Atom_objectRT.SetType(element) ;
			Atom_objectRT.SetIndex(kk) ;
			Atom_objectRT.SetX(x+Lx) ;
			Atom_objectRT.SetY(y+Ly) ;
			Atom_objectRT.SetZ(z) ;
			SeaAtomList[kk] = Atom_objectRT ;

//			std::cout << "atom index= " << kk << std::endl ;

			// Add 1 to index for the central cell atom
			kk++ ;
		} // end if
	}
	file.close() ;
}

//   method get the Atom in a location of the AtomList
Atom SeaAtomSet::GetSeaAtom(size_t dim_SeaAtomSetCard) 
{
	return SeaAtomList[dim_SeaAtomSetCard];
}

	

//**********************************************************************************
// Main function
//**********************************************************************************

int main()
{


	// Parameters:
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------	
	unsigned int NodeType = 1 ; // Tells the program what type of atoms will be clustered all other will be consider as sea atoms.
	size_t Natoms = 9 ; //8  7 Total number of atoms in the file (Node atoms + Sea atoms).
	size_t SeaNatoms = 63 ; //54   45  9*Natoms-2*9 ; 
	double R1 = 0.5 ; // Radius of type 1 atoms.
	double R2 = 0.5 ; // Radius of type 2 atoms.
	double Lx = 10.0, Ly = 10.0, Lz = 10.0 ; // Size of the box.
	double hc = R1 ; // Size of the partition of the cylinder height
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------		
	
	
	
	const double pi = 3.1415926535897 ;
	double RT ;
	
	double ijx, ijy, ijz, ikx, iky, ikz, modij, Vol, ikpar, modik, ikver, den ;
	double dv ;
	unsigned int Atom_count = 0 ;
	unsigned int memo = 0 ;
	unsigned int count = 0 ;
	unsigned int ifriends = 0 ;
	unsigned int nbins ;
	
	// Count execution time
	int start_s=clock();
	
	ofstream Frequencies ;
	Frequencies.open ("Frequencies.dat", ios::out | ios::trunc); // app=append
	ofstream Distances ;
	Distances.open ("Distances.dat", ios::out | ios::trunc); // app=append
	
	// Define the object set of atoms called AtomSet_object with Natoms
	// It is automatically initialized with the data in the file that we 
	// defined in the constructor

	// Construct Atom set
	AtomSet AtomSet_object(Natoms) ;
	// Construct Sea atom set
	SeaAtomSet SeaAtomSet_object(Natoms) ;

	//std::cout << SeaAtomSet_object.GetSeaAtom(40).GetZ() << std::endl ;

	//return (0) ;
	std::cout << std::endl ;
	std::cout << BOLDYELLOW << "  ------------------------------- " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "  TOPOLOGICAL CLUSTERING WITH PBC " << RESET << std::endl ;
	std::cout << BOLDYELLOW << "  ------------------------------- " << RESET << std::endl ;
	std::cout << std::endl ;

	for ( size_t i = 0; i < Natoms-1; ++i ) {
		if ( AtomSet_object.GetAtom(i).GetType() == NodeType )
		{
			ifriends = 0 ;
			for ( size_t j=i+1; j < Natoms; ++j ) {
				if ( AtomSet_object.GetAtom(j).GetType() == NodeType ) 
				{ 
				
				if ( i != j )
				{
					std::cout << BOLDYELLOW << "* PAIR OF NODES: " << i <<"-"<< j << RESET << std::endl ;

					double ix = AtomSet_object.GetAtom(i).GetX() ;
					double iy = AtomSet_object.GetAtom(i).GetY() ;
					double iz = AtomSet_object.GetAtom(i).GetZ() ;
					double jx = AtomSet_object.GetAtom(j).GetX() ;
					double jy = AtomSet_object.GetAtom(j).GetY() ;
					double jz = AtomSet_object.GetAtom(j).GetZ() ;

					// ---------------------------------------------
					// Edge inside the cell: e(i,j) inside the cell.
					// ---------------------------------------------

					ijx = jx - ix ;
					ijy = jy - iy ;
					ijz = jz - iz ;

					modij = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					Vol = pi*RT*RT*modij ;
					ijx = ijx/modij ;
					ijy = ijy/modij ;
					ijz = ijz/modij ;



//////					// BEGIN Cylinder partittion
//////					//(1) If the distance between atoms to clusterize is to small to fit a sea atom then they are connected


//////					//(2) If |d_ij| < h_c then we just define a single partition: 1 bin cylinder (otherwise nbins=0 and dv-> Infinity).
//////					if ( modij < hc )
//////					{
//////						nbins=1
//////						dv=pi*RT*RT*modij
//////					}
//////					else
//////					{
//////						nbins=int((mod_n-modij % hc)/hc)
//////						dv=pi*RT*RT*modij/Float(nbins)
//////					}
//////					// END Cylinder partittion



				
					Atom_count = 0 ;
				
					// Run over the rest of the atoms
					for ( size_t k = 0; k < Natoms; ++k ) 
					{
						if ( (ikpar > 0) && (ikpar < modij) && (ikver < RT)) // Determines whether the atom k is inside the cyclinder
						{
							if ( (k != i) && (k != j) )
							{						
								ikx = AtomSet_object.GetAtom(k).GetX() - ix ;
								iky = AtomSet_object.GetAtom(k).GetY() - iy ;	
								ikz = AtomSet_object.GetAtom(k).GetZ() - iz ;	
								modik = sqrt (ikx*ikx + iky*iky + ikz*ikz) ;						
								ikpar = ijx*ikx + ijy*iky + ijz*ikz ;		// Parallel component
								ikver = sqrt (modik*modik - ikpar*ikpar) ;	// Vertical component

								//BEGIN Binning along the cylinder axis
								//Once we know this atom is inside the cylinder we calculate to which local density contributes each section.
								//The variable to bin is clearly the projection of the atom position over the axis of the cylinder.
								//Determine the number of the bin data point x(i) falls within 
								n = Int( nbins*ikpar/modij + 1.0 ) ; // Determines in which bin is the atom k
								//Increment the count within that bin
								++bin[n] ;
								//END Binning along the cylinder axis

								// k atoms inside each bin
								unsigned int typei = AtomSet_object.GetAtom(i).GetType() ;
								unsigned int typek = AtomSet_object.GetAtom(k).GetType() ;
								
								if ( (typek ==  typei) ) 
								{ 
									if ( typei == 1)
									{
										RT = R1 + R1 ;
									}
									else
									{
										RT = R2 + R2 ;
									}
								}
								else 
								{
									RT = R1 + R2 ;
								}				
							} // If k is different from i or j			
						} // If atom inside cylinder i-j then count one atom to Atom_count
					} // Sum over all atoms


					
					std::cout << BOLDRED << "\t(Normal Cell) Number of atoms inside the cylinder " << i <<"-"<< j << " is " << Atom_count << RESET << endl ;
























					
					// This If can be removed to speed up just kept for strict consistency
					if ( i != j ) { 
						// Conectivity criterion
						den = Atom_count/Vol ;
						//cout << i << " " << j << " " << den << endl ;
						// computing i-friends
						if ( den < 0.000000001 ) {
						++ifriends ;				
						//cout << " density = " << i+1 << " " << j+1 << " " << ifriends << " " << den << endl ;	
						Distances << modij << endl ; // Store in file Distances.dat
						}// If density is lower than 0.000...
					}


					// ---------------------------------------------
					// Edge outside the cell: minimal distance image of j.
					// ---------------------------------------------

					// min(d(i,IL(j)), d(i,IR(j)), d(i,IT(j)), d(i,IB(j)), d(i,IL(j)))
					// Computation of IL(j):
					ijx = jx - Lx - ix ;
					ijy = jy - iy ;
					ijz = jz - iz ;
					double modiILj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of IR(j):
					ijx = jx + Lx - ix ;
					ijy = jy - iy ;
					ijz = jz - iz ;
					double modiIRj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of IB(j):
					ijx = jx - ix ;
					ijy = jy - Ly - iy ;
					ijz = jz - iz ;
					double modiIBj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of IT(j):
					ijx = jx - ix ;
					ijy = jy + Ly - iy ;
					ijz = jz - iz ;
					double modiITj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of ILB(j):
					ijx = jx - Lx - ix ;
					ijy = jy - Ly - iy ;
					ijz = jz - iz ;
					double modiILBj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of ILT(j):
					ijx = jx - Lx - ix ;
					ijy = jy + Ly - iy ;
					ijz = jz - iz ;
					double modiILTj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of IRB(j):
					ijx = jx + Lx - ix ;
					ijy = jy - Ly - iy ;
					ijz = jz - iz ;
					double modiIRBj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Computation of IRT(j):
					ijx = jx + Lx - ix ;
					ijy = jy + Ly - iy ;
					ijz = jz - iz ;
					double modiIRTj = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					// Check minimal distance and detect what image space to use.
					double ImageDistances[] = {modiILj, modiIRj ,modiIBj ,modiITj, modiILBj, modiILTj, modiIRBj, modiIRTj};
					string Image[] = {"L", "R", "B", "T", "LB", "LT", "RB", "RT"} ;
					double ImVecX[] = {-Lx, +Lx, 0.0, 0.0, -Lx, -Lx, +Lx, +Lx} ; // Traslation vector component X
					double ImVecY[] = {0.0, 0.0, -Ly, +Ly, -Ly, +Ly, -Ly, +Ly} ; // Traslation vector component Y
					double ImVecZ[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} ; // Traslation vector component Z

					double minimal_distance_image= *std::min_element(ImageDistances,ImageDistances+8) ;
					std::cout << "\tMinimal i-j distance image = " << minimal_distance_image << std::endl ;
					std::cout << "\tDistance i-j inside cell = " << modij << std::endl ;
					// Determine which image to use:
					for ( size_t k = 0; k < Natoms+1; ++k ) 
					{
						if ( ImageDistances[k] == minimal_distance_image ) 
						{
							std::cout << std::endl ;
							std::cout << "\tThe correct image for this pair is: I" << Image[k] << "(j)= (" << jx+ImVecX[k] << ", " << jy+ImVecY[k] << ", " << jz+ImVecZ[k] << " )" << std::endl ;
					
							ijx = jx + ImVecX[k] - ix ;
							ijy = jy + ImVecY[k] - iy ;
							ijz = jz + ImVecZ[k] - iz ;

							break ;
						}
					}
					// To detect the atoms inside the cylinde we need the supercell
					// Build here vec(i,min(image_j)) etc

					modij = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					Vol = pi*RT*RT*modij ;
					ijx = ijx/modij ;
					ijy = ijy/modij ;
					ijz = ijz/modij ;

					Atom_count = 0 ;
					
					// Run over the rest of the atoms
					for ( size_t k = 0; k < SeaNatoms; ++k ) 
					{		
						ikx = SeaAtomSet_object.GetSeaAtom(k).GetX() - ix ;
						iky = SeaAtomSet_object.GetSeaAtom(k).GetY() - iy ;	
						ikz = SeaAtomSet_object.GetSeaAtom(k).GetZ() - iz ;	
						modik = sqrt (ikx*ikx + iky*iky + ikz*ikz) ;						
						ikpar = ijx*ikx + ijy*iky + ijz*ikz ;
						ikver = sqrt (modik*modik - ikpar*ikpar) ;						
						// k atoms inside the cylinder
						
						RT = R1 + R2 ;
												
						if ( (ikpar > 0) && (ikpar < modij) && (ikver < RT)) ++Atom_count ;
					}// Sum over all atoms
					
					std::cout << BOLDRED << "\t(Image Cell) Number of atoms inside the image cylinder " << i <<"-"<< j << " is " << Atom_count << RESET << endl ;
					std::cout << std::endl ;
					std::cout << std::endl ;
				}// i and j are nodes (type=1)
			}// i not equal to j
			}// End j
			Frequencies << ifriends << endl ; // Store in file Frequencies.dat		
		}	
	}// End i
	
	Frequencies.close() ;
	Distances.close() ;

	
	
	//cout << "\033[1m Results: \033[0m" << endl ;
	cout << endl;
	//cout << "Mean number of first nighbours = \033[1m" << "\033[0m" << endl ;
	cout << endl;
	//cout << "Standard deviation of the number of first nighbours = \033[1m" << "\033[0m" << endl ;		
	
	// End counting time 	
 	int stop_s=clock();
 	
 	cout << endl;
	cout << "\tExecution time: \033[1m" << ( stop_s - start_s )/double(CLOCKS_PER_SEC) << "\033[0m seconds" << endl;
	cout << endl;
	
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

