! program test_tl_ad
!   implicit none
!   integer, parameter :: n = 100
!   real*8, parameter :: dt = 0.1D0
!   real*8, dimension(3) :: xinit, dyinit, dxinit, x, dx, x1, x2, x3, x_tl, x_ad
!   REAL*8, parameter :: p = 1.3D0, r = 0.6D0, b = 2.0D0 ! constants
!   real*8 :: tl_disp, ad_disp
!   integer :: i,m
  
!   do i=1,3
!      xinit(i) = sin(dble(i))
!   end do
!   x1 = xinit
!   call model(x1,dt,n)
!   do i=1,3
!      dxinit(i) = sin(dble(i + 123))
!   end do
!   dx = dxinit
!   do m=0,11
!      x2 = xinit + dx
!      x3 = xinit
!      x_tl = dx
!      call model(x2,dt,n)
!      call model_tl(x3,x_tl,dt,n)
!      write(*,*) m, (x2 - x1)/x_tl
!      dx = dx*0.1D0
!   end do

!   do i=1,3
!      dyinit(i) = log(cos(sin(dble(i)))**2)
!   end do

!   ! <M \delta x, M \delta x> = <\delta x, M* M \delta x>
!   x = xinit
!   x_tl = dxinit
!   call model_tl(x,x_tl,dt,n)
!   x = xinit
!   x_ad = x_tl
!   call model_ad(x,x_ad,dt,n)
!   tl_disp = sum( x_tl * x_tl )
!   ad_disp = sum( dxinit * x_ad )
!   write(*,*) 'test adj 1',tl_disp, ad_disp, 1.0 - tl_disp / ad_disp

!   ! <M \delta x, \delta y> = <\delta x, M* \delta y>
!   x = xinit
!   x_tl = dxinit
!   call model_tl(x,x_tl,dt,n)
!   x = xinit
!   x_ad = dyinit
!   call model_ad(x,x_ad,dt,n)
!   write(*,*) 'test adj 2',sum(x_tl*dyinit), sum(dxinit*x_ad), log(sum(x_tl*dyinit)/sum(dxinit*x_ad))
  
! contains

!!!!!!!!!!!!!!!!! ORIGINAL VERSION !!!!!!!!!!!!!!!!!!!!!
SUBROUTINE model_lorenz63(x,dt,nstep,params)
  implicit none
  REAL*8,INTENT(INOUT) :: x(3)
  REAL*8,INTENT(IN) :: dt ! constant
  INTEGER, INTENT(IN):: nstep
  real*8, intent(in) :: params(3)
  real*8 :: dxdt(3)
  integer :: j
  DO j = 1,nstep
     CALL calc_dxdt_lorenz63(x,dxdt,params)
     CALL step_lorenz63 (x,dxdt,dt)
  ENDDO
  return
END SUBROUTINE model_lorenz63
! ------------------------------
SUBROUTINE calc_dxdt_lorenz63(x,dxdt,params)
  implicit none
  REAL*8,INTENT(IN) :: x(3)
  REAL*8,INTENT(OUT):: dxdt(3)
  real*8, intent(in) :: params(3)
  !! params: p = 1.3D0, r = 0.6D0, b = 2.0D0 ! constants
  dxdt(1) = -params(1)*x(1)+params(1)*x(2)
  dxdt(2) = x(1)*(params(2)-x(3))-x(2)
  dxdt(3) = x(1)*x(2)-params(3)*x(3)
  return
END SUBROUTINE calc_dxdt_lorenz63
! ---------------------------
SUBROUTINE step_lorenz63(x,dxdt,dt)
  implicit none
  REAL*8,INTENT(INOUT):: x(3)
  REAL*8,INTENT(IN) :: dxdt(3),dt
  integer :: j
  DO j= 1,3
     x(j) = x(j)+dt*dxdt(j)
  ENDDO
  return
END SUBROUTINE step_lorenz63

!!!!!!!!!!!!!!!!! TANGENT LINEAR !!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE model_tl_lorenz63(x,x_tl,dt,nstep,params)
  implicit none
  REAL*8,INTENT(INOUT) :: x(3)
  REAL*8,INTENT(INOUT) :: x_tl(3)
  REAL*8,INTENT(IN) :: dt ! constant
  INTEGER, INTENT(IN):: nstep
  real*8, intent(in) :: params(3)
  real*8 :: dxdt_tl(3), dxdt(3)
  integer :: j

  DO j = 1,nstep
     CALL calc_dxdt_lorenz63(x,dxdt,params)
     CALL calc_dxdt_tl_lorenz63(x,x_tl,dxdt_tl,params)
     CALL step_lorenz63(x,dxdt,dt)
     CALL step_tl_lorenz63(x_tl,dxdt_tl,dt)
  ENDDO
  return
END SUBROUTINE model_tl_lorenz63
! ------------------------------
SUBROUTINE calc_dxdt_tl_lorenz63(x,x_tl,dxdt_tl,params)
  implicit none
  REAL*8,INTENT(IN) :: x(3)
  REAL*8,INTENT(IN) :: x_tl(3)
  REAL*8,INTENT(OUT):: dxdt_tl(3)
  real*8, intent(in) :: params(3)
  dxdt_tl(1) = -params(1)*x_tl(1)+params(1)*x_tl(2)
  dxdt_tl(2) = x_tl(1)*(params(2)-x(3)) - x(1)*x_tl(3) - x_tl(2)
  dxdt_tl(3) = x_tl(1)*x(2) + x(1)*x_tl(2) - params(3)*x_tl(3)
  return
END SUBROUTINE calc_dxdt_tl_lorenz63
! ---------------------------
SUBROUTINE step_tl_lorenz63(x_tl,dxdt_tl,dt)
  implicit none
  REAL*8,INTENT(INOUT):: x_tl(3)
  REAL*8,INTENT(IN) :: dxdt_tl(3),dt
  integer :: j
  DO j= 1,3
     x_tl(j) = x_tl(j)+dt*dxdt_tl(j)
  ENDDO
  return
END SUBROUTINE step_tl_lorenz63

!!!!!!!!!!!!!!!!! ADJOINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE model_ad_lorenz63(x,x_ad,dt,nstep,params)
  implicit none
  REAL*8,INTENT(INOUT) :: x(3)
  REAL*8,INTENT(INOUT) :: x_ad(3)
  REAL*8,INTENT(IN) :: dt ! constant
  INTEGER, INTENT(IN):: nstep
  real*8, intent(in) :: params(3)
  !!
  real*8 :: dxdt_ad(3)
  real*8 :: dxdt(3)
  real*8 :: xstore(3,nstep)
  integer :: j

  do j=1,nstep
     call calc_dxdt_lorenz63(x,dxdt,params)
     xstore(:,j) = x
     call step_lorenz63(x,dxdt,dt)
  end do
  dxdt_ad = 0.0
  DO j = nstep,1,-1
     x = xstore(:,j)
     CALL step_ad_lorenz63(x_ad,dxdt_ad,dt)
     CALL calc_dxdt_ad_lorenz63(x,x_ad,dxdt_ad,params)
  ENDDO
  return
END SUBROUTINE model_ad_lorenz63
! ------------------------------
SUBROUTINE calc_dxdt_ad_lorenz63(x,x_ad,dxdt_ad,params)
  implicit none
  REAL*8,INTENT(IN) :: x(3)
  REAL*8,INTENT(INOUT) :: x_ad(3)
  REAL*8,INTENT(INOUT):: dxdt_ad(3)
  real*8, intent(in) :: params(3)
!  x_tl(1)      | 1                |  x_tl(1)   
!  x_tl(2)    = |        1         |  x_tl(2)   
!  x_tl(3)    = |              1   |  x_tl(3)   
!  dxdt_tl(3)   | x(2)  x(1)  -b 0 |  dxdt_tl(3)
  x_ad(1)    = x_ad(1) + x(2)*dxdt_ad(3)
  x_ad(2)    = x_ad(2) + x(1)*dxdt_ad(3)
  x_ad(3)    = x_ad(3) -    params(3)*dxdt_ad(3)
  dxdt_ad(3) = 0.0

!  x_tl(1)      | 1                    |  x_tl(1)   
!  x_tl(3)      |                 1    |  x_tl(3)   
!  x_tl(2)    = |            1         |  x_tl(2)   
!  dxdt_tl(2)   | (r-x(3))  -1  -x(1) 0|  dxdt_tl(2)!
!  dxdt_tl(2) = x_tl(1)*(r-x(3)) - x(1)*x_tl(3) - x_tl(2)
  x_ad(1)    = x_ad(1) + (params(2)-x(3))*dxdt_ad(2)
  x_ad(3)    = x_ad(3) -   x(1) * dxdt_ad(2)
  x_ad(2)    = x_ad(2) -          dxdt_ad(2)
  dxdt_ad(2) = 0.0

!  x_tl(2)     = |0   1  0| x_tl(2)   
!  x_tl(1)       |1   0  0| x_tl(1)   
!  dxdt_tl(1)    |-p  p  0| dxdt_tl(1) !! p = params(1)
  x_ad(2) = x_ad(2) + params(1)*dxdt_ad(1)
  x_ad(1) = x_ad(1) - params(1)*dxdt_ad(1)
  dxdt_ad(1) = 0.0

  return
END SUBROUTINE calc_dxdt_ad_lorenz63
! ---------------------------
SUBROUTINE step_ad_lorenz63(x_ad,dxdt_ad,dt)
  implicit none
  REAL*8,INTENT(IN):: x_ad(3),dt
  REAL*8,INTENT(INOUT) :: dxdt_ad(3)
  integer :: j
  DO j= 3,1,-1
! dxdt_tl(j) = | 1      | dxdt_tl(j)
! x_tl(j)      | dt   1 | x_tl(j)   
     dxdt_ad(j) = dxdt_ad(j) + dt*x_ad(j)
  ENDDO
  return
END SUBROUTINE step_ad_lorenz63

! end program test_tl_ad

