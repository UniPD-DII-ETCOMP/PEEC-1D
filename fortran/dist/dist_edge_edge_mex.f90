!---- never change ------

#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    implicit none

! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
!---- end never change ------
	  
	  
! mwSize stuff for mexing
      mwSize mo,no,siz
	  
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
	  mwPointer NN_pr,EP_pr,nEdges_pr,nNodes_pr,n_thread_pr,dist_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: NN, EP_r
      integer,allocatable,dimension(:,:) :: EP
	  real*8,allocatable,dimension(:,:) :: dist 
	  real*8 :: nEdges_r,nNodes_r,n_thread_r
      integer :: nEdges,nNodes,n_thread
	  character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='log_dist_ee.txt',status='unknown')

!    Check to see INPUTS are numeric.
	  do ii = 1,5
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:fun_compute_Matrix_R:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  if(debu) write(66,*) 'sono arrivato al check numerico'
	  
!     Check that input #3 is INTEGER and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3 must be scalar.')
      endif	  
      siz = m*n
      nEdges_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(nEdges_pr, nEdges_r, siz) ! da double precision a reale
      nEdges=int(nEdges_r) ! da reale a intero
	  if(debu) write(66,*) 'input 3' 
	  
!     Check that input #4 is INTEGER and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be scalar.')
      endif	  
      siz = m*n
      nNodes_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(nNodes_pr, nNodes_r, siz) ! da double precision a reale
      nNodes=int(nNodes_r) ! da reale a intero
	  if(debu) write(66,*) 'input 4' 
	  
!     Check that input #5 is INTEGER and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be scalar.')
      endif	  
      siz = m*n
      n_thread_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(n_thread_pr, n_thread_r, siz) ! da double precision a reale
      n_thread=int(n_thread_r) ! da reale a intero
	  if(debu) write(66,*) 'input 5' 
      
	  
!	  Check that input #1 is REAL matrix and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. 3 .or. n .ne. nNodes) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be 3 x nNodes')
      endif	  
      siz = m*n
      NN_pr = mxGetPr(prhs(1))
	  allocate(NN(3,nNodes))
      call mxCopyPtrToReal8(NN_pr, NN, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 1'
	  
!     Check that input #2 is INTEGER matrix and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 2 .or. n .ne. nEdges) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be 2 x nEdges')
      endif	  
      siz = m*n
      EP_pr = mxGetPr(prhs(2))
	  allocate(EP_r(2,nEdges))
      allocate(EP(2,nEdges))
      call mxCopyPtrToReal8(EP_pr, EP_r, siz) ! da double precision a reale  
	  EP = int(EP_r)
	  if(debu) write(66,*) 'input 2'
      
	  	  
! call the computational subroutine.
      allocate(dist(nEdges,nEdges))
      call dist_edge_edge(NN,EP,nEdges,nNodes,n_thread,dist)

       if(debu) write(66,*) 'ho chiamato la subroutine e creato le matrici'
      
      deallocate(NN,EP_r,EP)
	  if(debu) write(66,*) 'ho deallocato le cose inutili'

! Create a matrix for the return argument 1
      mo=nEdges
      no=nEdges
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)

! Load the output 1 into a MATLAB array.
      dist_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(dist, dist_pr, siz)


      if(debu) write(66,*) 'ho convertito la matrice per matlab'

      deallocate(dist)
      
	  
      if(debu) write(66,*) 'ho deallocato la matrice e ora chiudo il logL.txt'
      if(debu) close(66)
      return
      end

