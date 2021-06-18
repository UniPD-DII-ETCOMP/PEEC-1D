subroutine P_stick_assemble_par2(nSticks,nNodes,NN,G,radius,npg,ne_cap_max,Cap_Elem,n_thread,P)
implicit none
integer :: nSticks,nNodes,ne_cap_max,ne_cap_hh,ne_cap_kk,npg,n_thread
integer :: Cap_Elem(ne_cap_max+1,nNodes),G(2,nSticks),hh,kk,ii,jj,ff
real*8 :: pi,epsilon0,P_self_jj_jj,P_self_jj_ii,P_mutual_jj_ii,alpha_jj,ll_jj,ll_ii,ll_tot_hh,ll_tot_kk,integ_self,integ_mutual
real*8 :: log_eps,P(nNodes,nNodes),radius(nSticks), glob_P
!real*8,dimension(8,15) :: gauss_P,gauss_W
real*8,dimension(npg) :: gauss_P,gauss_W
real*8 :: idE_hh(ne_cap_max),idE_kk(ne_cap_max),NN(3,nNodes),NN_edge(3,2),NN_jj(3,2),NN_ii(3,2),PPg_jj(3,npg)

!load gauss points and weights
! 100 format(15F20.15)   
! open(77,file='gauss_points_fortran.txt')
! open(88,file='gauss_weights_fortran.txt')
! do jj=1,8
! read(77,100) gauss_P(jj,:)
! read(88,100) gauss_W(jj,:)
! end do
! close(77)
! close(88)

call gauleg(npg,gauss_P,gauss_W) ! calcolo i punti e pesi di Gauss nella linea -[1 1]

call omp_set_num_threads(n_thread)
pi=4.d0*datan(1.d0)
epsilon0=8.85418781762d-12
P(1:nNodes,1:nNodes)=0.0d0
!!!!!!!!!! self-global-capacitance !!!!!!!!!!
!$OMP PARALLEL SHARED(nNodes,G,NN,Cap_Elem,ne_cap_max,radius,gauss_P,gauss_W,npg,pi,epsilon0,P)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(hh,jj,ii,ne_cap_hh,ff,idE_hh,P_self_jj_jj,P_self_jj_ii,ll_tot_hh,NN_edge,NN_jj,alpha_jj,PPg_jj,ll_jj,NN_ii,ll_ii,integ_self,log_eps)
do hh=1,nNodes
    ne_cap_hh=Cap_Elem(1,hh)
    idE_hh=Cap_Elem(2:ne_cap_max+1,hh)
    P_self_jj_jj=0.0d0
    P_self_jj_ii=0.0d0
    ll_tot_hh=0.0d0
    do jj=1,ne_cap_hh
        !self-local-capacitance
        NN_edge=NN(1:3,G(1:2,idE_hh(jj)))
        NN_jj(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)) !costruisco il jj-esimo elemento capacitivo (nodo hh e lato induttivo jj)
        NN_jj(1:3,2) = NN(1:3,hh)
        call fun_1_R_2stick_self(NN_jj,radius(idE_hh(jj)),alpha_jj)
        P_self_jj_jj=P_self_jj_jj+alpha_jj !da dividere per ll_tot^2
!        call Gauss_line_nvar(NN_jj,gauss_P,gauss_W,npg,PPg_jj,whg_jj,ll_jj)
        call Gauss_line_nvar(NN_jj,gauss_P,npg,PPg_jj,ll_jj)
        ll_tot_hh=ll_tot_hh+ll_jj
		!mutual-local-capacitance
        do ii=jj+1,ne_cap_hh
            NN_edge=NN(1:3,G(1:2,idE_hh(ii)))
            NN_ii(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)) !costruisco il ii-esimo elemento capacitivo (nodo hh e lato induttivo ii)
            NN_ii(1:3,2) = NN(1:3,hh)
            ll_ii=norm2(NN_ii(1:3,2)-NN_ii(1:3,1))
            !gauss
            integ_self=0.0d0
            do ff=1,npg
                call fun_1_R_stick(NN_ii,PPg_jj(1:3,ff),ll_ii,log_eps)
                integ_self=integ_self+gauss_W(ff)*log_eps
            end do
            P_self_jj_ii=P_self_jj_ii+0.5*ll_jj*integ_self
        end do
    end do
	P(hh,hh)=(P_self_jj_jj+2*P_self_jj_ii)/(4*pi*epsilon0*ll_tot_hh**2)
end do
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
!!!!!!!!!! mutual-global-capacitance !!!!!!!!!!
!$OMP PARALLEL SHARED(nNodes,G,NN,Cap_Elem,ne_cap_max,gauss_P,gauss_W,npg,pi,epsilon0,P)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(hh,kk,jj,ii,ne_cap_hh,ff,idE_hh,P_mutual_jj_ii,ll_tot_hh,NN_edge,NN_jj,PPg_jj,ll_jj,NN_ii,ll_ii,integ_mutual,log_eps,ne_cap_kk,idE_kk,ll_tot_kk)
do hh=1,nNodes
    ne_cap_hh=Cap_Elem(1,hh) ! numero di elementi del nodo hh
	idE_hh=Cap_Elem(2:ne_cap_max+1,hh)	
	glob_P=0.0
    do kk=hh+1,nNodes
	   ne_cap_kk=Cap_Elem(1,kk)
	   idE_kk=Cap_Elem(2:ne_cap_max+1,kk)
	   ll_tot_hh=0.0	   
	   P_mutual_jj_ii = 0.0
	   do jj = 1,ne_cap_hh
	      NN_edge=NN(1:3,G(1:2,idE_hh(jj))) ! prende gli estremi dell'elemento
          NN_jj(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)) !costruisco il jj-esimo elemento capacitivo (nodo hh e lato induttivo jj)
          NN_jj(1:3,2) = NN(1:3,hh) 
          !call Gauss_line_nvar(NN_jj,gauss_P,gauss_W,npg,PPg_jj,whg_jj,ll_jj) ! metto i punti di gauss
		  call Gauss_line_nvar(NN_jj,gauss_P,npg,PPg_jj,ll_jj) ! metto i punti di gauss
		  ll_tot_hh=ll_tot_hh+ll_jj ! aggiorno la lunghezza
		  ll_tot_kk=0.0		  
	      do ii = 1,ne_cap_kk 
             NN_edge=NN(1:3,G(1:2,idE_kk(ii)))
             NN_ii(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)) !costruisco il ii-esimo elemento capacitivo (nodo kk e lato induttivo ii)
             NN_ii(1:3,2) = NN(1:3,kk)
             ll_ii=norm2(NN_ii(1:3,2)-NN_ii(1:3,1))
             ll_tot_kk=ll_tot_kk+ll_ii 	   			 
			 integ_mutual=0.0;
			 do ff=1,npg
                    call fun_1_R_stick(NN_ii,PPg_jj(1:3,ff),ll_ii,log_eps)
                    integ_mutual=integ_mutual+gauss_W(ff)*log_eps			 
			 end do
	         P_mutual_jj_ii=P_mutual_jj_ii+0.5*ll_jj*integ_mutual
		  end do 
		  !glob_P=glob_P+P_mutual_jj_ii
	   end do  
	   P(hh,kk)=P_mutual_jj_ii/(ll_tot_kk*ll_tot_hh*4*pi*epsilon0)
	   P(kk,hh)=P(hh,kk)
	end do
end do
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
return
end subroutine P_stick_assemble_par2

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
    
!************** fun_1_R_2stick_self ****************
subroutine fun_1_R_2stick_self(NN,radius,alpha)
     
implicit none
real*8 :: NN(3,2),ll,radius,aa,bb,cc,alpha

ll=norm2(NN(1:3,2)-NN(1:3,1))
aa=log(ll/radius+sqrt((ll/radius)**2+1))
bb=sqrt(1+(radius/ll)**2)
cc=radius/ll
alpha=2.0*ll*(aa-bb+cc)
return
end subroutine fun_1_R_2stick_self   
    
!************** fun_1_R_stick ***************
subroutine fun_1_R_stick(PP,PP0,ll,log_eps)
implicit none
real*8 :: PP(3,2),PP0(3)
real*8 :: ll,eps,log_eps,ri,rf

ri=norm2(PP0-PP(1:3,1))
rf=norm2(PP0-PP(1:3,2))
!ll=norm2(PP(1:3,2)-PP(1:3,1))
eps=ll/(ri+rf)
log_eps=log((1+eps)/(1-eps))
return
end subroutine fun_1_R_stick  