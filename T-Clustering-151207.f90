!===========================================================================
!===========================================================================
!
! 2015-12-07 Topological Clustering
! 
! Version: 151207
!
! By Carles Triguero
!
!===========================================================================
!===========================================================================
!
! Brute-force version: That means the program is still very slow and can be optimized.
!
! We need two input files:
!
! Gas.xyz
! Liquid.xyz
!
! And also set the number of atoms in each of this files as:
!
! N_gas:        Number of atoms to be clusterize ("gas atoms" in our case).
! N_liquid:     Background specie disconecting gas domains (This is not clusterized).
!
! (nx,ny,nz)   : Is the vector linking two gas atoms.
! mod_n        : Is the modulus of the vector n.
!
! dn           : Is the subcylinder height. 
! dv           : Is the subcylinder volume. 
!
! Set the size of the simulation box in order to take properly account of the boundary 
! conditions assumed periodic in this program:
!
! Lx, Ly, Lz.
!
! Two parameters control the process:
!
! (1) r_c:      The radius of the cylindrical domain.
! (2) rho_T:    The cutoff density of liquid atoms inside the cylinder.
!  http://www.asciiworld.com/-Geometry-.html


Module Parameters
 
 Implicit None
 
!Numerical precission
 Integer, Parameter :: Prec=selected_real_kind(8,8)

!-------------------------------------------
! BEGIN Main Parameters of the Clustering process
!------------------------------------------- 
 Real(kind=Prec), Parameter :: d_c= 2.47                        !Cylinder diameter in sea atoms diameter units d_g (related with \sigma_{ll}) (domain radial threshold).
 Real(kind=Prec), Parameter :: h_c= 3.0                         !Cylinder height in sea atoms diameter units d_g (to perform the cylinder partition). 
 Real(kind=Prec), Parameter :: rho_T= 0.001                      !Local number density threshold. If set very high we get the real local density distribution.
 Real(kind=Prec), Parameter :: D_T=160.0                        !D-Clustering threshold normaly min[Lx/2,Ly/2,Lz/2].
 Real(kind=Prec), Parameter :: d_l=5.0                          !Liquid atom diameter.
 Real(kind=Prec), Parameter :: d_g=3.4                          !Gas atom diameter.
 Real(kind=prec), Parameter :: Lx=500.0,Ly=500.0,Lz=500.0       !Size of the whole space. Useful when considering periodic boundary conditions.
 Real(kind=prec), Parameter :: Lx2=Lx/2.0_Prec
 Real(kind=prec), Parameter :: Ly2=Ly/2.0_Prec
 Real(kind=prec), Parameter :: Lz2=Lz/2.0_Prec
 Real(kind=prec), Parameter :: Dmax=Sqrt(Lx2*Lx2+Ly2*Ly2+Lz2*Lz2)
!Parameters NOT scaled in d_l: Real ordinary scale.
 Real(kind=Prec), Parameter :: Pi=4.0_Prec*Atan(1.0_Prec)
 Real(kind=Prec), Parameter :: BL=h_c*d_l
 Real(kind=Prec), Parameter :: r_c=d_c*d_l/2.0_Prec
 Real(kind=Prec), Parameter :: r_c2=r_c*r_c
 Real(kind=Prec), Parameter :: PiR_T2=Pi*r_c2
!-------------------------------------------
! END Main Parameters of the Clustering process
!------------------------------------------- 
  
 Integer :: N_Gas
 Integer :: N_Liq
 Integer :: AllocateStatus
 Integer :: DeAllocateStatus
 Integer :: i,j,k,l,Block,Cluster_Index,Number_clusters,skip
 Integer :: cluster_atoms,Old_Cluster_Size
 Integer :: n,nbins
 
 Real(kind=prec) :: dx,dy,dz
 Real(kind=prec) :: xmean,ymean,zmean
 Real(kind=prec) :: xi,xj
 Real(kind=prec) :: yi,yj
 Real(kind=prec) :: zi,zj
 Real(kind=prec) :: nx,ny,nz
 Real(kind=prec) :: nxn,nyn,nzn
 Real(kind=prec) :: mod_nd2,mod_dmnd2 
 Real(kind=prec) :: mod_n,mod_d
 Real(kind=prec) :: d_H,d_P
 Real(kind=prec) :: rho_local,dv
 Character(Len=15) :: number1,number2
 
End Module Parameters






Module Tensorial

 Use Parameters
 Implicit None
 
 Type Vectors
    Real(kind=Prec), Allocatable :: x(:),y(:),z(:)
    Real(kind=Prec), Allocatable :: xk(:),yk(:),zk(:)
    Integer, Allocatable :: c(:),v(:),bin(:)
    Integer :: N_atoms(1:2)
 End Type Vectors

End Module Tensorial


















Program Clustering

 Use Parameters
 Use Tensorial
 Implicit None  
 Type(Vectors) :: V
 
 Call Sets_Cardinality(V)
 
 N_Gas=V%N_atoms(1)
 N_Liq=V%N_atoms(2)
 
 Call Read_Positions(V)
 Call IO

  
  
 
 
      
 !-------------------------------------------------     
 !BEGIN Positions loop over Gas atoms i j to connect
 !-------------------------------------------------

 
 Do i=1,N_Gas
    Do j=i+1,N_Gas
        
       xi=V%x(i)
       yi=V%y(i)
       zi=V%z(i)
     
       xj=V%x(j)
       yj=V%y(j)
       zj=V%z(j)
     
       nx=xj-xi
       ny=yj-yi
       nz=zj-zi
       
! !       !BEGIN PBC (periodic boundary conditions)
! !        If(Abs(nx) > Lx2) nx=Lx-Abs(dx)        
! !        If(Abs(ny) > Ly2) ny=Ly-Abs(dy)
! !        If(Abs(nz) > Lz2) nz=Lz-Abs(dz)  
! !       !END PBC (periodic boundary conditions)
     
       mod_n = Sqrt(nx*nx + ny*ny + nz*nz)
       mod_nd2 = mod_n/2.0_Prec
      
      !-----------------------------------------
      !Begin Distance restriction
      !-----------------------------------------
       If(mod_n > D_T) Goto 50 ! 50 is the jumper for the i & j loop
      !-----------------------------------------
      !End Distance restriction
      !-----------------------------------------
      
      !Normalized components      
       nxn=nx/mod_n
       nyn=ny/mod_n
       nzn=nz/mod_n
   
      
      
  
!----------------------------------------      
!BEGIN Cylinder Partition
!----------------------------------------      
!(1) If the distance between atoms to clusterize is to small to fit a sea atom then they are connected
!
 If(mod_n < d_g + d_l) Then
    Allocate(V%bin(1:2), STAT = AllocateStatus) 
    If (AllocateStatus /= 0) STOP "Not enough memory for BIN vector my friend!"
              !!Write(*,*) 'dij < d_g => connected',i,j
    Goto 60 !Directly considers i-j atoms as connected as they are very very close no sea atoms fit in between.
 Endif

!
!(2) If |d_ij| < h_c then we just define a single partition: 1 bin cylinder (otherwise nbins=0 and dv-> Infinity).
!
 If(mod_n < BL) Then
    nbins=1
    dv=PiR_T2*mod_n
    Else
    nbins=Int((mod_n-Mod(mod_n,BL))/BL)
    dv=PiR_T2*mod_n/Float(nbins)
 Endif
        
 Write(15,*)dv,nbins,mod_n
        
        
 Allocate(V%bin(1:nbins+1), STAT = AllocateStatus) 
 If (AllocateStatus /= 0) STOP "Not enough memory for BIN vector my friend!"
!Initialize
 V%bin=0
!----------------------------------------      
!END Cylinder Partition
!----------------------------------------           
       
       
       
       
       
       
       
       
       
      

       !-----------------------------------------
       !----------------------------------------- 
       !BEGIN Loop over Liquid atoms 
       !-----------------------------------------
       !-----------------------------------------      
        Open(20,File='Liq.xyz')
        Do k=1,N_Liq  
           dx = V%xk(k) - xi 
           dy = V%yk(k) - yi 
           dz = V%zk(k) - zi 
           mod_d = Sqrt( dx*dx + dy*dy + dz*dz )
           
          !---------------------------------------
          !BEGIN Fast check if atom k could be inside cylinder
          !---------------------------------------
           !First exclussion discards all atoms out of a surrounding sphere to the cylinder ij  SIMPLIFY DO NOT NEED SQRT
           !mod_dmnd2 = Sqrt( (dx-0.5*nx)*(dx-0.5*nx) + (dy-0.5*ny)*(dy-0.5*ny) + (dz-0.5*nz)*(dz-0.5*nz) )
           !If( mod_dmnd2 > Sqrt( mod_nd2*mod_nd2 + r_c2 ) ) Goto 30
           mod_dmnd2 = (dx-0.5_Prec*nx)*(dx-0.5_Prec*nx) + &
                       (dy-0.5_Prec*ny)*(dy-0.5_Prec*ny) + &
                       (dz-0.5_Prec*nz)*(dz-0.5_Prec*nz)
                       
           If( mod_dmnd2  >  (mod_nd2*mod_nd2 + r_c2)  ) Goto 30
          !---------------------------------------
          !END Fast check if atom k  could be inside cylinder
          !---------------------------------------                      
           !Horizontal component
           d_H = nxn*dx + nyn*dy + nzn*dz
           !Perpendicular component
           d_P = Sqrt( mod_d*mod_d - d_H*d_H )

           
          !---------------------------------------
          !BEGIN Check if atom k inside cylinder
          !---------------------------------------
          !BEGIN Liquid atoms that are inside the cylinder radius r_c count
           If(d_P > r_c) Goto 30
          !END Liquid atoms that are inside the cylinder radius r_c count          
          !BEGIN Liquid atoms that are inside the cylinder length mod_d count                                ____
           If(d_H < 0.0_Prec)   Goto 30   !----> Ensures the atom it is over the base of the cylinder  -->  o ()__()   ____
           If(d_H > mod_n) Goto 30   !----> Ensures is below the top of the cylinder ------------------------> ()__() o
          !END Liquid atoms that are inside the cylinder length mod_d count
          !!!!!!!!!!IN ONE LINE: If((d_P > r_c).Or.(d_H < 0).Or.(d_H > mod_d ESTA MAL NO ES mod_n)) Goto 30
          !---------------------------------------
          !END Check if atom k inside cylinder
          !---------------------------------------
       
       
          !-----------------------------------------
          !BEGIN Binning along the cylinder axis
          !-----------------------------------------
           !Once we know this atom is inside the cylinder we calculate to which local density contributes each section.
           !The variable to bin is clearly the projection of the atom position over the axis of the cylinder.
           !Determine the number of the bin data point x(i) falls within 
           !n = Int( nbins*( d_H - xmin )/( xmax - xmin ) + 1.0 )
           n = Int( nbins*d_H/mod_n + 1.0 )
           !Increment the count within that bin 
           V%bin(n) = V%bin(n) + 1           
          !-----------------------------------------
          !END Binning along the cylinder axis
          !-----------------------------------------   
          
          
           30 Continue        
        Enddo      
        Close(20)
        
       !-----------------------------------------
       !-----------------------------------------
       !END Loop over Liquid atoms 
       !-----------------------------------------
       !-----------------------------------------
      
      
       !-----------------------------------------
       !BEGIN 
       !-----------------------------------------
        !Now we change the following decision on whether the atoms are connected or not based 
        !on the global density of the cyclinder by the local density.      
        !If any of the sections of the cylinder has a larger density then the atoms are not connected
        Do k=1,nbins
           rho_local=Float(V%bin(k))/dv  
           Write(25,*)rho_local ! It is breakkkk not all the range is shown unless very high threshold densities are set!!!
           Write(35,*)V%bin(k)    ! It is breakkkk not all the range is shown unless very hight threshold densities are set!!!
           If(rho_local >= rho_T) Then
              !Write(*,*)i,j,'NOT connected'
              !If one of the sections has larger density than the threshold.
              !This atom two atom j,j will never be connected.
              !Let us check another pair Goto 50
              Goto 70 !jumper
           Endif
        Enddo
        
        60 Continue
       !-----------------------------------------
       !END 
       !-----------------------------------------  
       
       
       !-----------------------------------------  
       !BEGIN Connected
       !-----------------------------------------  
       
        !Write(*,*)i,j,'YES connected'
       
       !BEGIN Detecting all atoms already belonging to the cluster c(i) or to the same as i
        Allocate(V%v(1:N_Gas), STAT = AllocateStatus)
        If (AllocateStatus /= 0) STOP "Not enough memory for V my friend!"
    
        Old_Cluster_Size=0
        Do l=1,N_Gas
           If(V%c(l) == V%c(i)) Then
              Old_Cluster_Size=Old_Cluster_Size+1
              V%v(Old_Cluster_Size)=l
           Endif
        Enddo
       !END Detecting all atoms already belonging to the cluster c(i) or to the same as i

       !BEGIN Renaming the cluster index c(i) to c(j) in order to add atom j in the cluster: Cluster_Size=Old_Cluster_Size+1
       ! What we basically do is putting all the the atoms in the same cluster i renamed to the index of cluster for atom j
        Do l=1,Old_Cluster_Size
           V%c(V%v(l))=V%c(j)
        Enddo
    
        Deallocate(V%v)     
        !-----------------------------------------  
        !END Connected
        !-----------------------------------------  
        
        70 Continue
        
        Deallocate(V%bin)
        
        50 Continue !Jumper for the i & j loop
         
        
     Enddo
  Enddo
        
  Close(15)
  Close(25)
 !-------------------------------------------------     
 !END Positions loop over Gas atoms i j to connect
 !-------------------------------------------------
  Deallocate(V%xk)
  Deallocate(V%yk)
  Deallocate(V%zk)
 !-----------------      
 !END Clustering   
 !-----------------
  
 !---------------------------    
 !BEGIN Cluster storage
 !---------------------------  
   
 !BEGIN Detect non-empty cluster names without repetition
  Allocate(V%v(1:N_gas), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory for V vector my friend!"

  Number_clusters=0
 !Loop over posible cluster names j
  Do j=1,N_Gas
     skip=0
    !Loop over all atom indexes
     Do i=1,N_Gas
       !Selects the cluster j, and counts as 1 if not empty 
       !in order to not count multiple we define skip
       !Thus v gives the non empty cluster name (from 1 to Number_clusters)
        If((V%c(i) == j).And.(skip == 0)) Then
           Number_clusters=Number_clusters+1
           V%v(Number_clusters)=j
          !Once we detect a non-empty cluster skip=1 to go to another cluster name j (exit i loop).
           skip=1
        Endif           
     Enddo
  Enddo
 !END Detect non-empty cluster names without repetition
 
 !BEGIN Different filenames
  Write (number1, *) Block
  number1 = Adjustl(number1)
  
  Open(20, file='INFORMATION_Clusters_'//Trim(number1)//'.dat')
 
  Write(20,*)"Detected ",Number_clusters," clusters for rho_threshold=",rho_T," clusters for r_threshold=",r_c
  Write(20,*)"Distance_clustering=",D_T," Atomic_diameter=",d_g
  Write(20,*)"RADIUS_CYLINDER      ","LENGTH_CYLINDER      ",Number_clusters
  Write(20,*)"   "       
  Write(*,'(F7.3,5x,I4)')rho_T,Number_clusters

  Open(30, file='Cluster_center_positions_'//Trim(number1)//'.xyz')     
  Open(40, file='Number.dat')   
 !END Different filenames  
     
  Cluster_Index=0  
  Do j=1,Number_clusters
   !BEGIN Different filenames for each cluster
    Write (number2, *) j
    number2 = Adjustl(number2)
    Open(60, file='Cluster_'//Trim(number2)//'.dat')
   !END Different filenames for each cluster
   
    Cluster_Index=Cluster_Index+1                
    Write(20,*)"Cluster number ",Cluster_Index
    Write(20,*)"Atom_Index    X_position    Y_position    Z_position"
        
    cluster_atoms=0
    xmean=0.0
    ymean=0.0
    zmean=0.0
   !Run over all the atoms and classify to which cluster it belongs
    Do i=1,N_Gas
     !We save the positions of all atoms i with cluster name=v(j)
      If(V%c(i) == V%v(j)) Then           
       !In a single file ordered
        Write(20,*)i,V%x(i),V%y(i),V%z(i)
       !In different files for each cluster
        Write(60,*)V%x(i),V%y(i),V%z(i)
       !This counts the number of atoms in each cluster, cluster size.    
        cluster_atoms=cluster_atoms+1
       !This is still a not well defined mean position when boundary conditions act. It only works for bulk atoms.
        xmean=xmean+V%x(i)
        ymean=ymean+V%y(i)
        zmean=zmean+V%z(i)
      Endif           
    Enddo
    Write(*,*) 'Cluster index and number of atoms'
    Write(40,*) cluster_atoms
    Write(30,*)xmean/Float(cluster_atoms),ymean/Float(cluster_atoms),zmean/Float(cluster_atoms)
  Enddo
  
  Deallocate(V%v, STAT = DeAllocateStatus)
  Deallocate(V%x, STAT = DeAllocateStatus)
  Deallocate(V%y, STAT = DeAllocateStatus)
  Deallocate(V%z, STAT = DeAllocateStatus)
     
 !---------------------------    
 !END Cluster storage
 !---------------------------  
 
End Program  Clustering


























Subroutine Sets_Cardinality(V)

 Use Parameters
 Use Tensorial
 Implicit None  
 Type(Vectors) :: V

!Size of the system through bash commands.
 Call execute_command_line ('wc -l < "Gas.xyz" > N_Gas.num')
 Call execute_command_line ('wc -l < "Liq.xyz" > N_Liq.num')
  
!Read number of atoms to clusterize
 Open(10,file="N_Gas.num")
 Read(10,*)V%N_atoms(1)
 Close(10)

!Read number of sea atoms
 Open(10,file="N_Liq.num")
 Read(10,*)V%N_atoms(2)
 Close(10)

!Delete files with the number of atoms
 Call execute_command_line ('rm N_Gas.num N_Liq.num')
   
End Subroutine Sets_Cardinality











Subroutine Read_Positions(V)

 Use Parameters
 Use Tensorial
 Implicit None  
 Type(Vectors) :: V
 
  
!Read positions of atoms to clusterize
 Open(10,File="Gas.xyz")
!BEGIN Allocate position vectors
 Allocate(V%x(1:N_Gas), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for x vector my friend!"
 Allocate(V%y(1:N_Gas), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for y vector my friend!"
 Allocate(V%z(1:N_Gas), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for z vector my friend!"
!END Allocate position vectors
 Block=1 !Do Block=1,Nblock
!We need to keep this vector to then spit out the clusters.
!We do not need that for liquid atoms they are just used to measure the density.
 Do i=(Block-1)*N_Gas+1,(Block-1)*N_Gas+N_Gas
    Read(10,*) V%x(i-(Block-1)*N_gas),V%y(i-(Block-1)*N_gas),V%z(i-(Block-1)*N_gas)   
 Enddo
 Close(10)
!END Read positions of atoms to clusterize



!Read positions of sea atoms 
!Before was done without vectors reading each time!
!BEGIN Read Liq positions
 Open(10,File='Liq.xyz')
!BEGIN Allocate position vectors
 Allocate(V%xk(1:N_Liq), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for xk vector my friend!"
 Allocate(V%yk(1:N_Liq), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for yk vector my friend!"
 Allocate(V%zk(1:N_Liq), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for zk vector my friend!"
!END Allocate position vectors
 Do i=1,N_Liq
    Read(10,*) V%xk(i),V%yk(i),V%zk(i)
 Enddo
 Close(10)
!END Read positions of sea atoms 


!Allocate and initialize cluster matrix each atom with a cluster index
 Allocate(V%c(1:N_Gas), STAT = AllocateStatus)
 If (AllocateStatus /= 0) STOP "Not enough memory for C vector my friend!"
 Do i=1,N_Gas
    V%c(i)=i
 Enddo
 
  Write(*,*)''
  Write(*,*)''
  Write(*,*)' _______      _____ _           _            _             ' 
  Write(*,*)'|__   __|    / ____| |         | |          (_)            '  
  Write(*,*)'   | |______| |    | |_   _ ___| |_ ___ _ __ _ _ __   __ _ '
  Write(*,*)"   | |______| |    | | | | / __| __/ _ \ '__| | '_ \ / _` |"
  Write(*,*)'   | |      | |____| | |_| \__ \ ||  __/ |  | | | | | (_| |'
  Write(*,*)'   |_|       \_____|_|\__,_|___/\__\___|_|  |_|_| |_|\__, |'
  Write(*,*)'                                                      __/ |'
  Write(*,*)'                                                     |___/ '
  Write(*,*)''
  Write(*,*)''
  Write(*,*)'Topological Clustering        '
  Write(*,*)'                               '
  Write(*,*)'Version: 151207               '
  Write(*,*)''
  Write(*,*)'Used parameters:              '
  Write(*,*)''
  Write(*,*)''
  Write(*,'(A,I8)')' Gas Atoms detected: ',N_Gas
  Write(*,'(A,I8)')' Liq Atoms detected: ',N_Liq   
  Write(*,*)''
  Write(*,'(A,3x,F7.4,1x,A)')' d_c =',d_c,'in d_l units'
  Write(*,'(A,1x,F5.2,1x,A)')' h_c =',h_c,'in d_l units'
  Write(*,'(A,1x,F9.4)')' rho_T =',rho_T
  Write(*,'(A,3x,F7.4)')' d_l =',d_l
  Write(*,'(A,3x,F7.4)')' d_g =',d_g
  Write(*,'(A,3x,F8.4)')' D_T =',D_T,' Distance clustering threshold'
  Write(*,*)'--------------------------------------'
 

 
End Subroutine Read_Positions


Subroutine IO
 Open(15,File='Local_volumes.dat')
 Open(25,File='Number_volume_densities.dat')
 Open(35,File='Number_densities.dat')
 Open(27,File='hc_distribution.dat')
End Subroutine IO
 
 

! 
! Subroutine Cylinder_Partition(V,nbins)
! !Partitioning of the cylinder
! 
! !Here we are on a fixed cylinder. We need to allocate the number of bins taking into account
! !The axis length of the cylinder.
! !Allocate vector bin and initiallize
! !Number of bins:  3d_g means that we allow 3 atoms aligned in the axis
! 
!  Use Parameters
!  Use Tensorial
!  Implicit None
!  Type(Vectors) :: V
! 
! 
!         
! End Subroutine Cylinder_Partition
