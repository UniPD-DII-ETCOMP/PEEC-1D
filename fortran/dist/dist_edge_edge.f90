subroutine dist_edge_edge(NN,EP,nEdges,nNodes,n_thread,dist)

implicit none
integer :: EP(2,nEdges), nEdges, nNodes, n_thread
real*8 :: NN(3,nNodes), dist(nEdges,nEdges)
integer :: n,m
real*8, dimension(3) :: centre_n,centre_m

call omp_set_num_threads(n_thread)
dist(:,:)=0.0d0
!$OMP PARALLEL SHARED(nEdges,EP,NN,dist)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(n,m,centre_n,centre_m)
do n=1,nEdges
    !! edge n
    centre_n=0.5*(NN(:,EP(1,n))+NN(:,EP(2,n))) !centre of n-th edge (dual edge)
    do m=n,nEdges
        !! edge m
        centre_m=0.5*(NN(:,EP(1,m))+NN(:,EP(2,m))); !centre of m-th edge (dual edge)
        dist(n,m)=norm2(centre_n-centre_m);
        dist(m,n)=dist(n,m)   
    end do  
end do
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
return
end subroutine dist_edge_edge