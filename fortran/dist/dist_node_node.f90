subroutine dist_node_node(NN,nNodes,n_thread,dist)
implicit none

integer :: nNodes,kk,jj,n_thread
real*8 :: NN(3,nNodes), dist(nNodes,nNodes)

call omp_set_num_threads(n_thread)
dist(:,:)=0.0d0
!$OMP PARALLEL SHARED(nNodes,NN,dist)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(kk,jj)
do kk=1,nNodes
do jj=kk,nNodes
dist(kk,jj)=norm2(NN(1:3,kk)-NN(1:3,jj))
dist(jj,kk)=dist(kk,jj)
end do
end do
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
return
end subroutine dist_node_node