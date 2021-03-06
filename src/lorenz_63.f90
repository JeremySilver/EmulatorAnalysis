SUBROUTINE model_lorenz63(x,xout,dt,nstep,params)
  implicit none
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(INOUT) :: x(n)
  REAL*8,INTENT(IN) :: dt ! constant
  INTEGER, INTENT(IN):: nstep
  real*8,intent(out) :: xout(n,nstep+1)
  real*8, intent(in) :: params(n)
  real*8 :: dxdt(n)
  integer :: j
  xout(:,1) = x
  DO j = 1,nstep
     CALL calc_dxdt_lorenz63(x,dxdt,params)
     CALL step_lorenz63 (x,dxdt,dt)
     xout(:,j+1) = x
  ENDDO
  return
END SUBROUTINE model_lorenz63
! ------------------------------
SUBROUTINE calc_dxdt_lorenz63(x,dxdt,params)
  implicit none
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(IN) :: x(n)
  REAL*8,INTENT(OUT):: dxdt(n)
  real*8, intent(in) :: params(n)
  !! params: p = 1.3D0, r = 0.6D0, b = 2.0D0 ! constants
  dxdt(1) = -params(1)*x(1)+params(1)*x(2)
  dxdt(2) = x(1)*(params(2)-x(3))-x(2)
  dxdt(3) = x(1)*x(2)-params(3)*x(3)
  return
END SUBROUTINE calc_dxdt_lorenz63
! ---------------------------
SUBROUTINE step_lorenz63(x,dxdt,dt)
  implicit none
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(INOUT):: x(n)
  REAL*8,INTENT(IN) :: dxdt(n),dt
  integer :: j
  DO j= 1,n
     x(j) = x(j)+dt*dxdt(j)
  ENDDO
  return
END SUBROUTINE step_lorenz63

!!!!!!!!!!!!!!!!! TANGENT LINEAR !!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE model_tl_lorenz63(x,x_tl,dt,nstep,params)
  implicit none
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(INOUT) :: x(n)
  REAL*8,INTENT(INOUT) :: x_tl(n)
  REAL*8,INTENT(IN) :: dt ! constant
  INTEGER, INTENT(IN):: nstep
  real*8, intent(in) :: params(n)
  real*8 :: dxdt_tl(n), dxdt(n)
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
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(IN) :: x(n)
  REAL*8,INTENT(IN) :: x_tl(n)
  REAL*8,INTENT(OUT):: dxdt_tl(n)
  real*8, intent(in) :: params(n)
  dxdt_tl(1) = -params(1)*x_tl(1)+params(1)*x_tl(2)
  dxdt_tl(2) = x_tl(1)*(params(2)-x(3)) - x(1)*x_tl(3) - x_tl(2)
  dxdt_tl(3) = x_tl(1)*x(2) + x(1)*x_tl(2) - params(3)*x_tl(3)
  return
END SUBROUTINE calc_dxdt_tl_lorenz63
! ---------------------------
SUBROUTINE step_tl_lorenz63(x_tl,dxdt_tl,dt)
  implicit none
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(INOUT):: x_tl(n)
  REAL*8,INTENT(IN) :: dxdt_tl(n),dt
  integer :: j
  DO j= 1,n
     x_tl(j) = x_tl(j)+dt*dxdt_tl(j)
  ENDDO
  return
END SUBROUTINE step_tl_lorenz63

!!!!!!!!!!!!!!!!! ADJOINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE model_ad_lorenz63(x,x_ad,dt,nstep,params)
  implicit none
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(INOUT) :: x(n)
  REAL*8,INTENT(INOUT) :: x_ad(n)
  REAL*8,INTENT(IN) :: dt ! constant
  INTEGER, INTENT(IN):: nstep
  real*8, intent(in) :: params(n)
  !!
  real*8 :: dxdt_ad(n)
  real*8 :: dxdt(n)
  real*8 :: xstore(n,nstep)
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
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(IN) :: x(n)
  REAL*8,INTENT(INOUT) :: x_ad(n)
  REAL*8,INTENT(INOUT):: dxdt_ad(n)
  real*8, intent(in) :: params(n)
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
  integer, parameter :: n=3 ! number of variables
  REAL*8,INTENT(IN):: x_ad(n),dt
  REAL*8,INTENT(INOUT) :: dxdt_ad(n)
  integer :: j
  DO j= n,1,-1
! dxdt_tl(j) = | 1      | dxdt_tl(j)
! x_tl(j)      | dt   1 | x_tl(j)   
     dxdt_ad(j) = dxdt_ad(j) + dt*x_ad(j)
  ENDDO
  return
END SUBROUTINE step_ad_lorenz63


