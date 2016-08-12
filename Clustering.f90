!===========================================================================
!===========================================================================
!
! 2015-02-25 Clustering
! 
! By Carles Triguero
!
!===========================================================================
!===========================================================================


Program Clustering

  Implicit None  
  Real, Dimension(:), Allocatable :: x,y,z
  Real, Dimension(:,:), Allocatable :: d
  Integer, Dimension(:), Allocatable :: c,v
  Integer :: N,AllocateStatus,DeAllocateStatus
  Integer :: i,j,k,l,Block,NBlock,Cluster_Index,Number_clusters,skip
  Integer :: cluster_atoms,d_minloc(2),Old_Cluster_Size
  Real :: dmin,dmax,dx,dy,dz,Lx,Ly,Lz,Lx2,Ly2,Lz2,dthreshold,xmean,ymean,zmean
  Character*15 number1,number2
  
  
  !-------------------
  dthreshold=20.00 !cutoff distance
  !-------------------
  
  Write(*,*) 'D-Clustering running...'
  
  N=2699

  Lx=150.0
  Ly=150.0
  Lz=34.6
  
  Lx2=Lx/2.0
  Ly2=Ly/2.0
  Lz2=Lz/2.0
  
!  NBlock=101
  
  Open(10,File="Gas.xyz")
  
  dmax=Sqrt((Lx/2.0)**2.0 + (Ly/2.0)**2.0 + (Lz/2.0)**2.0)
  
  !Do Block=1,Nblock
  Block=1
  
  
  
 !BEGIN Allocate position vectors
  Allocate(x(1:N), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory my friend!"
  Allocate(y(1:N), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory my friend!"
  Allocate(z(1:N), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory my friend!"
 !END Allocate position vectors
   
   
   
 !BEGIN Read positions
  Do i=(Block-1)*N+1,(Block-1)*N+N
    Read(10,*) x(i-(Block-1)*N),y(i-(Block-1)*N),z(i-(Block-1)*N)   
  Enddo
 !END Read positions
    
    
    
 !BEGIN Distance matrix 
  Allocate(d(1:N,1:N), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory my friend!"      
 !Set all matrix to maximum possible distances
  d=dmax
 !Compute the distances without repeating (only N(N-1)/2 elements due to symmetry)
  Do i=1,N
    Do j=i+1,N
     !BEGIN With Boundary conditions: 
      dx=x(j)-x(i)
      dy=y(j)-y(i)
      dz=z(j)-z(i)
      !If(Abs(dx) .Gt. Lx2) dx=Lx-Abs(dx)        
      !If(Abs(dy) .Gt. Ly2) dy=Ly-Abs(dy)
      !If(Abs(dz) .Gt. Lz2) dz=Lz-Abs(dz)  
     !END With Boundary conditions:
      d(i,j)=Sqrt(dx**2.0 + dy**2.0 + dz**2.0)
    Enddo
  Enddo
 !Simetrize matrix
  Do j=1,N
    Do i=j+1,N
      d(i,j)=d(j,i)
    Enddo
  Enddo
 !END Distance matrix 
    
      
    
 !BEGIN Clustering   
 !Allocate and initialize cluster matrix
  Allocate(c(1:N), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory my friend!"
  Do i=1,N
    c(i)=i
  Enddo
      
  Do k=1,N*(N-1)/2
    dmin=Minval(d)
    
   !This determines the cluster threshold otherwise you end with a single cluster
    If (dmin .Gt. dthreshold) Goto 50
        
    d_minloc=Minloc(d)
    i=d_minloc(1)
    j=d_minloc(2)

   !Detecting all atoms already belonging to the cluster c(i) or to the same as i
    Allocate(v(1:N), STAT = AllocateStatus)
    If (AllocateStatus /= 0) STOP "Not enough memory my friend!"
    
    Old_Cluster_Size=0
    Do l=1,N
      If(c(l) .Eq. c(i)) Then
        Old_Cluster_Size=Old_Cluster_Size+1
        v(Old_Cluster_Size)=l
      Endif
    Enddo

   !Renaming the cluster and adding new element: Cluster_Size=Old_Cluster_Size+1
    Do l=1,Old_Cluster_Size
      c(v(l))=c(j)
    Enddo
    
    Deallocate(v, STAT = DeAllocateStatus)
     
   !In order to do not count again we artificially say that this distance is huge
    d(i,j) = dmax
    d(j,i) = dmax
     
  Enddo
  
  Deallocate(d, STAT = DeAllocateStatus)
 !END Clustering

 
 
  50 Continue

     
     
 !BEGIN Cluster storage
   
   
 !BEGIN Detect non-empty cluster names without repetition
  Allocate(v(1:N), STAT = AllocateStatus)
  If (AllocateStatus /= 0) STOP "Not enough memory my friend!"

  Number_clusters=0
 !Loop over posible cluster names j
  Do j=1,N
    skip=0
   !Loop over all atom indexes
    Do i=1,N
     !Selects the cluster j, and counts as 1 if not empty 
     !in order to not count multiple we define skip
     !Thus v gives the non empty cluster name (from 1 to Number_clusters)
      If((c(i).Eq.j).And.(skip.Eq.0)) Then
        Number_clusters=Number_clusters+1
        v(Number_clusters)=j
       !Once we detect a non-empty cluster skip=1 to go to another cluster name j (exit i loop).
        skip=1
      Endif           
    Enddo
  Enddo
 !END Detect non-empty cluster names without repetition
 !BEGIN Different filenames
  Write (number1, *) Block
  number1 = Adjustl(number1)
  
  Open(20, file='Clusters_'//Trim(number1)//'.dat')
  Write(20,*)"Detected ",Number_clusters," clusters for dthreshold=",dthreshold
  Write(20,*)"   "       

  Open(30, file='Cluster_center_positions_'//Trim(number1)//'.xyz')     
 !END Different filenames  
     
  Cluster_Index=0  
  Do j=1,Number_clusters
   !BEGIN Different filenames for each cluster
    Write (number2, *) j
    number2 = Adjustl(number2)
    Write(*,*) 'Cluster=',number2
    Open(60, file='Cluster_'//Trim(number2)//'.dat')
   !END Different filenames for each cluster
   
    Cluster_Index=Cluster_Index+1                
    Write(20,*)"Cluster number ",Cluster_Index
    Write(20,*)"Atom_Index    X_position    Y_position    Z_position"
        
    cluster_atoms=0.0
    xmean=0.0
    ymean=0.0
    zmean=0.0
   !Run over all the atoms and classify to which cluster it belongs
    Do i=1,N
     !We save the positions of all atoms i with cluster name=v(j)
      If(c(i).Eq.v(j)) Then           
       !In a single file ordered
        Write(20,*)i,x(i),y(i),z(i)
       !In different files for each cluster
        Write(60,*)x(i),y(i),z(i)
       !This counts the number of atoms in each cluster, cluster size.    
        cluster_atoms=cluster_atoms+1
       !This is still a not well defined mean position when boundary conditions act. It only works for bulk atoms.
        xmean=xmean+x(i)
        ymean=ymean+y(i)
        zmean=zmean+z(i)
      Endif           
    Enddo
    
    Write(30,*)xmean/Float(cluster_atoms),ymean/Float(cluster_atoms),zmean/Float(cluster_atoms)
  Enddo
  
  Deallocate(v, STAT = DeAllocateStatus)
  Deallocate(x, STAT = DeAllocateStatus)
  Deallocate(y, STAT = DeAllocateStatus)
  Deallocate(z, STAT = DeAllocateStatus)
     
     ! .....extra info....
     !Open(40,file='Number_of_clusters.dat')
     !Write(40,*)(Float(Block)-1.0)*0.1,Number_clusters
     !Write(*,*) 'Number of points=',N
     !Write(*,*) 'Detected Clusters=',Number_clusters
     
 !END Cluster storage
 
 !Enddo !!BLOCKS
    

End Program  Clustering

  
     

