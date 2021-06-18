subroutine L_stick_assemble_par(nSticks,nNodes,NN,G,radius,npg,Wint,n_thread,L)
implicit none
integer :: nSticks,nNodes,G(2,nSticks),npg,hh,Wint,kk,ii,jj,n_thread
real*8 :: NN(3,nNodes),L(nSticks,nSticks),radius(nSticks)
real*8, dimension(3,2) :: PPh,PPk
real*8 :: ll_h,ll_k,aa,bb,cc,eps,log_eps,integ,ri,rf
!real*8,dimension(8,15) :: gauss_P,gauss_W
real*8,dimension(npg) :: gauss_P,gauss_W
real*8 :: ut_h(3),ut_k(3),PPgh(3,npg)

L(:,:)=0.0d0
!load gauss points and weights
!100 format(15F20.15)   
!open(77,file='gauss_points_fortran.txt')
!open(88,file='gauss_weights_fortran.txt')
!do jj=1,8
!read(77,100) gauss_P(jj,:)
!read(88,100) gauss_W(jj,:)
!end do
!close(77)
!close(88)

call gauleg(npg,gauss_P,gauss_W) ! calcolo i punti e pesi di Gauss nella linea -[1 1]

call omp_set_num_threads(n_thread)
!$OMP PARALLEL SHARED(nSticks,G,NN,L,radius,npg,Wint,gauss_P,gauss_W)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(hh,PPh,PPgh,ll_h,aa,bb,cc,ut_h,kk,PPk,ll_k,ut_k,integ,ri,rf,eps,log_eps)
do hh=1,nSticks
    !self-inductance
    PPh=NN(1:3,G(1:2,hh))
    call Gauss_line_nvar(PPh,gauss_P,npg,PPgh,ll_h)
    aa=log(ll_h/radius(hh)+sqrt((ll_h/radius(hh))**2+1))
    bb=sqrt(1+(radius(hh)/ll_h)**2)
    cc=radius(hh)/ll_h
    L(hh,hh)=2*1.0d-7*ll_h*(aa-bb+cc+float(Wint)*0.25)
    ut_h=(PPh(1:3,2)-PPh(1:3,1))/ll_h
    !mutual-inductance
    do kk=hh+1,nSticks
        PPk=NN(1:3,G(1:2,kk))
        ll_k=norm2(PPk(1:3,2)-PPk(1:3,1))
        ut_k=(PPk(1:3,2)-PPk(1:3,1))/ll_k
        integ = 0.0d0
        do ii=1,npg
		ri=norm2(PPgh(1:3,ii)-PPk(1:3,1))
		rf=norm2(PPgh(1:3,ii)-PPk(1:3,2))
		eps=ll_k/(ri+rf)
		log_eps=log((1+eps)/(1-eps))
		integ=integ+gauss_W(ii)*log_eps*dot_product(ut_h,ut_k)
        end do
        L(hh,kk)=1.0d-7*0.5*ll_h*integ
	    L(kk,hh)=L(hh,kk)
    end do
end do
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
end subroutine L_stick_assemble_par


!--------------------------------------------------------------------------------
! Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
! integration of polynomial functions.
!      For normalized lower and upper limits of integration -1.0 & 1.0, and
! given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
! containing the abscissas and weights of the Gauss-Legendre n-point quadrature
! formula.
!--------------------------------------------------------------------------------
subroutine gauleg(ngp,xabsc,weig)
implicit none
integer i,j,m
real(kind=8) p1,p2,p3,pp,z,z1
integer, intent(IN) :: ngp !# of Gauss Points
real(kind=8), intent(OUT) :: xabsc(ngp),weig(ngp)
real(kind=8) :: eps2, pi
eps2=3.0d-15
pi=3.141592653589793d0
m=(ngp+1)/2
!Roots are symmetric in the interval so only need to find half of them
do i=1,m
z=cos(pi*(i-0.25d0)/(ngp+0.5d0)) !starting approximation
!Newton's method
100     p1 = 1.0d0
p2 = 0.0d0
!Loop up the recurrence relation to get the Legendre
!polynomial evaluated at z
do j=1,ngp
p3=p2
p2=p1
p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
end do
!p1 is now the desired Legendre polynomial. We next compute pp,
!its derivative, by a standard relation involving also p2, the
!polynomial of one lower order.
pp=ngp*(z*p1-p2)/(z*z-1.0d0)
z1=z
z=z1-p1/pp             !Newton's Method
if (dabs(z-z1) .gt. EPS2) GOTO  100
xabsc(i)=-z                    		! Roots will be bewteen -1.0 & 1.0
xabsc(ngp+1-i)=+z             		! and symmetric about the origin
weig(i)=2.0d0/((1.0d0-z*z)*pp*pp)	! Compute the weight and its
weig(ngp+1-i)=weig(i)               ! symmetric counterpart
end do
end subroutine gauleg

!************** Gauss_line_nvar ******************
subroutine Gauss_line_nvar(NN,gauss_P,np,PPg,ll)
implicit none    
real*8,dimension(np) :: gauss_P
real*8 :: PPg(3,np),NN(3,2),ll
integer :: np
ll = norm2(NN(1:3,2)-NN(1:3,1))
PPg(1,1:np) = spread(NN(1,1),1,np)+0.5*(NN(1,2)-NN(1,1))*(1.0+gauss_P)
PPg(2,1:np) = spread(NN(2,1),1,np)+0.5*(NN(2,2)-NN(2,1))*(1.0+gauss_P)
PPg(3,1:np) = spread(NN(3,1),1,np)+0.5*(NN(3,2)-NN(3,1))*(1.0+gauss_P)
return
end subroutine Gauss_line_nvar
!************** Gauss_line_nvar ****************** DIMITRI
! subroutine Gauss_line_nvar(NN,gauss_P,gauss_W,np,PPg,wg,ll)
! implicit none    
! real*8,dimension(8,15) :: gauss_P,gauss_W
! real*8 :: pg(np),wg(np),PPg(3,np),NN(3,2),ll  
! integer :: np
! ll = norm2(NN(1:3,2)-NN(1:3,1))
! if (mod(np,2)==0) then !even number of points
    ! pg = [gauss_P(1:np/2,np-1),-gauss_P(1:np/2,np-1)]
    ! wg = [gauss_W(1:np/2,np-1),gauss_W(1:np/2,np-1)]
! else if (mod(np,2)==1) then !odd
    ! pg = [gauss_P(1,np-1),gauss_P(2:(np+1)/2,np-1),-gauss_P(2:(np+1)/2,np-1)]
    ! wg = [gauss_W(1,np-1),gauss_W(2:(np+1)/2,np-1),gauss_W(2:(np+1)/2,np-1)]
! end if  
! PPg(1,1:np) = spread(NN(1,1),1,np)+0.5*(NN(1,2)-NN(1,1))*(1+pg)
! PPg(2,1:np) = spread(NN(2,1),1,np)+0.5*(NN(2,2)-NN(2,1))*(1+pg)
! PPg(3,1:np) = spread(NN(3,1),1,np)+0.5*(NN(3,2)-NN(3,1))*(1+pg)
! return
! end subroutine Gauss_line_nvar