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
	  mwPointer nSticks_pr,nNodes_pr,NN_pr,G_pr,radius_pr,npg_pr,Wint_pr,n_thread_pr,L_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: NN, G_r, L
      integer,allocatable,dimension(:,:) :: G
	  real*8,allocatable,dimension(:) :: radius 
	  real*8 :: nSticks_r,nNodes_r,npg_r,Wint_r,n_thread_r
      integer :: nSticks,nNodes,npg,Wint,n_thread
	  character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='logL.txt',status='unknown')

!    Check to see INPUTS are numeric.
	  do ii = 1,8
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:fun_compute_Matrix_R:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  if(debu) write(66,*) 'sono arrivato al check numerico'
	  
!     Check that input #1 is INTEGER and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar.')
      endif	  
      siz = m*n
      nSticks_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(nSticks_pr, nSticks_r, siz) ! da double precision a reale
      nSticks=int(nSticks_r) ! da reale a intero
	  if(debu) write(66,*) 'input 1'

!     Check that input #8 is INTEGER and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 8 must be scalar.')
      endif	  
      siz = m*n
      n_thread_pr = mxGetPr(prhs(8))
      call mxCopyPtrToReal8(n_thread_pr, n_thread_r, siz) ! da double precision a reale
      n_thread=int(n_thread_r) ! da reale a intero
	  if(debu) write(66,*) 'input 8' 	  
	  
!     Check that input #2 is INTEGER and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be scalar.')
      endif	  
      siz = m*n
      nNodes_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(nNodes_pr, nNodes_r, siz) ! da double precision a reale
      nNodes=int(nNodes_r) ! da reale a intero
	  if(debu) write(66,*) 'input 2' 
	  
!     Check that input #6 is INTEGER and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 6 must be scalar.')
      endif	  
      siz = m*n
      npg_pr = mxGetPr(prhs(6))
      call mxCopyPtrToReal8(npg_pr, npg_r, siz) ! da double precision a reale
      npg=int(npg_r) ! da reale a intero
	  if(debu) write(66,*) 'input 6' 
	  
!     Check that input #7 is INTEGER and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7 must be scalar.')
      endif	  
      siz = m*n
      Wint_pr = mxGetPr(prhs(7))
      call mxCopyPtrToReal8(Wint_pr, Wint_r, siz) ! da double precision a reale
      Wint=int(Wint_r) ! da reale a intero
	  if(debu) write(66,*) 'input 7'       
	  
!	  Check that input #3 is REAL matrix and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. 3 .or. n .ne. nNodes) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3 must be 3 x nNodes')
      endif	  
      siz = m*n
      NN_pr = mxGetPr(prhs(3))
	  allocate(NN(3,nNodes))
      call mxCopyPtrToReal8(NN_pr, NN, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 3'
	  
!     Check that input #4 is INTEGER matrix and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 2 .or. n .ne. nSticks) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be 2 x nSticks')
      endif	  
      siz = m*n
      G_pr = mxGetPr(prhs(4))
	  allocate(G_r(2,nSticks))
      allocate(G(2,nSticks))
      call mxCopyPtrToReal8(G_pr, G_r, siz) ! da double precision a reale  
	  G = int(G_r)
	  if(debu) write(66,*) 'input 4'
	  
!     Check that input #5 is INTEGER vector and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. nSticks) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be 1 x nSticks')
      endif	  
      siz = m*n
      radius_pr = mxGetPr(prhs(5))
	  allocate(radius(nSticks))
      call mxCopyPtrToReal8(radius_pr, radius, siz) ! da double precision a reale  
	  if(debu) write(66,*) 'input 5'	  
      
	  	  
! call the computational subroutine.
      allocate(L(nSticks,nSticks))
      call L_stick_assemble_par(nSticks,nNodes,NN,G,radius,npg,Wint,n_thread,L)

       if(debu) write(66,*) 'ho chiamato la subroutine e creato le matrici'
      
      deallocate(NN,G_r,G,radius)
	  if(debu) write(66,*) 'ho deallocato le cose inutili'

! Create a matrix for the return argument 1
      mo=nSticks
      no=nSticks
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)

! Load the output 1 into a MATLAB array.
      L_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(L, L_pr, siz)


      if(debu) write(66,*) 'ho convertito la matrice per matlab'

      deallocate(L)
      
	  
      if(debu) write(66,*) 'ho deallocato la matrice e ora chiudo il logL.txt'
      if(debu) close(66)
      return
      end

