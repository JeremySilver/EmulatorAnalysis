subroutine model_lorenz96(x,xout,n,dt,nstep,f)
  implicit none
  real*8,intent(inout) :: x(n)
  integer, intent(in) :: n ! number of variables
  real*8,intent(in) :: dt ! timestep
  real*8,intent(in) :: f ! forcing
  integer, intent(in):: nstep
  real*8,intent(out) :: xout(n,nstep+1)
  integer :: j

  xout(:,1) = x(:)
  do j = 1,nstep
     ! write(*,*) x(1)
     call rk4_lorenz96(x,n,dt,f)
     xout(:,j+1) = x(:)
  enddo

  return
end subroutine model_lorenz96
! ------------------------------


subroutine calc_dxdt_lorenz96(x,n,dxdt,f)
  implicit none
  ! arguments
  real*8,intent(in) :: x(n)
  integer, intent(in) :: n ! number of variables
  real*8,intent(out):: dxdt(n)
  real*8,intent(in) :: f ! forcing
  ! local variables
  integer :: i

  ! initial cases
  i = 1
  dxdt(i) = - x(n-1)*x(n) + x(n)*x(i+1) - x(i) + f
  i = 2
  dxdt(i) = - x(n)*x(i-1) + x(i-1)*x(i+1) - x(i) + f
  ! most cases
  do i=3,n-1
     dxdt(i) = - x(i-2)*x(i-1) + x(i-1)*x(i+1) - x(i) + f
  end do
  ! end case
  i = n
  dxdt(i) = - x(i-2)*x(i-1) + x(i-1)*x(1) - x(i) + f
  return
end subroutine calc_dxdt_lorenz96
! ---------------------------

subroutine rk4_lorenz96(x,n,dt,f)
  implicit none
  real*8,intent(inout) :: x(n)
  integer, intent(in) :: n ! number of variables
  real*8,intent(in) :: dt ! timestep
  real*8,intent(in) :: f ! forcing

  integer :: i
  real*8 :: x2(n), x3(n), x4(n)
  real*8 :: k1(n), k2(n), k3(n), k4(n)
  call calc_dxdt_lorenz96(x,n,k1,f)
  do i=1,n
     x2(i) = x(i) + 0.5d0 * dt * k1(i)
  end do
  call calc_dxdt_lorenz96(x2,n,k2,f)
  do i=1,n
     x3(i) = x(i) + 0.5d0 * dt * k2(i)
  end do
  call calc_dxdt_lorenz96(x3,n,k3,f)
  do i=1,n
     x4(i) = x(i) +         dt * k3(i)
  end do
  call calc_dxdt_lorenz96(x4,n,k4,f)
  
  do i=1,n
     x(i) = x(i) + dt *( k1(i) + 2.0d0 *( k2(i) + k3(i) ) + k4(i) )/6.0d0
  end do
  ! k1 = f (x                  , sigma , rho , beta )
  ! k2 = f (x +0.5d0 * dt * k1 , sigma , rho , beta )
  ! k3 = f (x +0.5d0 * dt * k2 , sigma , rho , beta )
  ! k4 = f (x +        dt * k3 , sigma , rho , beta )
  ! x = x + dt *( k1 + 2.0d0 *( k2 + k3 ) + k4 )/6.0d0
  return
end subroutine rk4_lorenz96

!!!!!!!!!!!!!!!!!! tangent linear !!!!!!!!!!!!!!!!!!!!!!!

subroutine model_tl_lorenz96(x, x_tl, n, dt, nstep, f)
  implicit none
  real*8, intent(inout) :: x(n)
  real*8, intent(inout) :: x_tl(n)
  integer, intent(in) :: n ! number of variables
! timestep
  real*8, intent(in) :: dt
! forcing
  real*8, intent(in) :: f
  integer, intent(in) :: nstep
  integer :: j
  do j=1,nstep
    call rk4_tl_lorenz96(x, x_tl, n, dt, f)
  end do
  return
end subroutine model_tl_lorenz96


subroutine calc_dxdt_tl_lorenz96(x, x_tl, n, dxdt, dxdt_tl, f)
  implicit none
! number of variables
  integer, intent(in) :: n
! arguments
  real*8, intent(in) :: x(n)
  real*8, intent(in) :: x_tl(n)
  real*8, intent(out) :: dxdt(n)
  real*8, intent(out) :: dxdt_tl(n)
! forcing
  real*8, intent(in) :: f
! local variables
  integer :: i
! initial cases
  i = 1
  dxdt_tl = 0.0d0
  dxdt_tl(i) = - x_tl(n-1)*x(n) - x(n-1)*x_tl(n) + x_tl(n)*x(i+1) + x(n)*x_tl(i+1) - x_tl(i)
  dxdt(i) = - x(n-1)*x(n) + x(n)*x(i+1) - x(i) + f
  i = 2
  dxdt_tl(i) = - x_tl(n)*x(i-1) - x(n)*x_tl(i-1) + x_tl(i-1)*x(i+1) + x(i-1)*x_tl(i+1) - x_tl(i)
  dxdt(i) = - x(n)*x(i-1) + x(i-1)*x(i+1) - x(i) + f
! most cases
  do i=3,n-1
    dxdt_tl(i) = - x_tl(i-2)*x(i-1) - x(i-2)*x_tl(i-1) + x_tl(i-1)*x(i+1) + x(i-1)*x_tl(i+1) - x_tl(i)
    dxdt(i) = - x(i-2)*x(i-1) + x(i-1)*x(i+1) - x(i) + f
  end do
! end cases
  i = n
  dxdt_tl(i) = - x_tl(i-2)*x(i-1) - x(i-2)*x_tl(i-1) + x_tl(i-1)*x(1) + x(i-1)*x_tl(1) - x_tl(i)
  dxdt(i) = - x(i-2)*x(i-1) + x(i-1)*x(1) - x(i) + f
  return
end subroutine calc_dxdt_tl_lorenz96

subroutine rk4_tl_lorenz96(x, x_tl, n, dt, f)
  implicit none
! number of variables
  integer, intent(in) :: n
  real*8, intent(inout) :: x(n)
  real*8, intent(inout) :: x_tl(n)
! timestep
  real*8, intent(in) :: dt
! forcing
  real*8, intent(in) :: f
  integer :: i
  real*8 :: x2(n), x3(n), x4(n)
  real*8 :: x2_tl(n), x3_tl(n), x4_tl(n)
  real*8 :: k1(n), k2(n), k3(n), k4(n)
  real*8 :: k1_tl(n), k2_tl(n), k3_tl(n), k4_tl(n)
  call calc_dxdt_tl_lorenz96(x, x_tl, n, k1, k1_tl, f)
  x2_tl = 0.0d0
  do i=1,n
    x2_tl(i) = x_tl(i) + 0.5d0*dt*k1_tl(i)
    x2(i) = x(i) + 0.5d0*dt*k1(i)
  end do
  call calc_dxdt_tl_lorenz96(x2, x2_tl, n, k2, k2_tl, f)
  x3_tl = 0.0d0
  do i=1,n
    x3_tl(i) = x_tl(i) + 0.5d0*dt*k2_tl(i)
    x3(i) = x(i) + 0.5d0*dt*k2(i)
  end do
  call calc_dxdt_tl_lorenz96(x3, x3_tl, n, k3, k3_tl, f)
  x4_tl = 0.0d0
  do i=1,n
    x4_tl(i) = x_tl(i) + dt*k3_tl(i)
    x4(i) = x(i) + dt*k3(i)
  end do
  call calc_dxdt_tl_lorenz96(x4, x4_tl, n, k4, k4_tl, f)
  do i=1,n
    x_tl(i) = x_tl(i) + dt*(k1_tl(i)+2.0d0*(k2_tl(i)+k3_tl(i))+k4_tl(i))/6.0d0
    x(i) = x(i) + dt*(k1(i)+2.0d0*(k2(i)+k3(i))+k4(i))/6.0d0
  end do
! k1 = f (x                  , sigma , rho , beta )
! k2 = f (x +0.5d0 * dt * k1 , sigma , rho , beta )
! k3 = f (x +0.5d0 * dt * k2 , sigma , rho , beta )
! k4 = f (x +        dt * k3 , sigma , rho , beta )
! x = x + dt *( k1 + 2.0d0 *( k2 + k3 ) + k4 )/6.0d0
  return
end subroutine rk4_tl_lorenz96

!!!!!!!!!!!!!!!!!! adjoint !!!!!!!!!!!!!!!!!!!!!!!

subroutine model_ad_lorenz96(x, x_ad, n, dt, nstep, f)
  implicit none
  real*8, intent(inout) :: x(n)
  real*8, intent(inout) :: x_ad(n)
  integer, intent(in) :: n ! number of variables
  real*8, intent(in) :: dt ! timestep
  real*8, intent(in) :: f ! forcing
  integer, intent(in) :: nstep
  integer :: j
  real*8 :: xtraj(n,nstep)
  
  do j=1,nstep
    xtraj(:,j) = x
    call rk4_lorenz96(x, n, dt, f)
  end do
  do j=nstep,1,-1
     x = xtraj(:,j)
     call rk4_ad_lorenz96(x, x_ad, n, dt, f)
  end do
end subroutine model_ad_lorenz96

subroutine rk4_ad_lorenz96(x, x_ad, n, dt, f)
  implicit none
! number of variables
  integer, intent(in) :: n
  real*8, intent(inout) :: x(n)
  real*8, intent(inout) :: x_ad(n)
! timestep
  real*8, intent(in) :: dt
! forcing
  real*8, intent(in) :: f
  integer :: i
  real*8 :: x2(n), x3(n), x4(n)
  real*8 :: x2_ad(n), x3_ad(n), x4_ad(n)
  real*8 :: k1(n), k2(n), k3(n), k4(n)
  real*8 :: k1_ad(n), k2_ad(n), k3_ad(n), k4_ad(n)
  real*8 :: temp_ad
  call calc_dxdt_lorenz96(x, n, k1, f)
  do i=1,n
    x2(i) = x(i) + 0.5d0*dt*k1(i)
  end do
  call calc_dxdt_lorenz96(x2, n, k2, f)
  do i=1,n
    x3(i) = x(i) + 0.5d0*dt*k2(i)
  end do
  call calc_dxdt_lorenz96(x3, n, k3, f)
  do i=1,n
    x4(i) = x(i) + dt*k3(i)
  end do
  k1_ad = 0.0d0
  k2_ad = 0.0d0
  k3_ad = 0.0d0
  k4_ad = 0.0d0
  do i=n,1,-1
    temp_ad = dt*x_ad(i)/6.0d0
    k1_ad(i) = k1_ad(i) + temp_ad
    k2_ad(i) = k2_ad(i) + 2.0d0*temp_ad
    k3_ad(i) = k3_ad(i) + 2.0d0*temp_ad
    k4_ad(i) = k4_ad(i) + temp_ad
  end do
  x4_ad = 0.0d0
  call calc_dxdt_ad_lorenz96(x4, x4_ad, n, k4, k4_ad, f)
  do i=n,1,-1
    x_ad(i) = x_ad(i) + x4_ad(i)
    k3_ad(i) = k3_ad(i) + dt*x4_ad(i)
    x4_ad(i) = 0.0d0
  end do
  x3_ad = 0.0d0
  call calc_dxdt_ad_lorenz96(x3, x3_ad, n, k3, k3_ad, f)
  do i=n,1,-1
    x_ad(i) = x_ad(i) + x3_ad(i)
    k2_ad(i) = k2_ad(i) + dt*0.5d0*x3_ad(i)
    x3_ad(i) = 0.0d0
  end do
  x2_ad = 0.0d0
  call calc_dxdt_ad_lorenz96(x2, x2_ad, n, k2, k2_ad, f)
  do i=n,1,-1
    x_ad(i) = x_ad(i) + x2_ad(i)
    k1_ad(i) = k1_ad(i) + dt*0.5d0*x2_ad(i)
    x2_ad(i) = 0.0d0
  end do
  call calc_dxdt_ad_lorenz96(x, x_ad, n, k1, k1_ad, f)
end subroutine rk4_ad_lorenz96


subroutine calc_dxdt_ad_lorenz96(x, x_ad, n, dxdt, dxdt_ad, f)
  implicit none
! number of variables
  integer, intent(in) :: n
! arguments
  real*8, intent(in) :: x(n)
  real*8 :: x_ad(n)
  real*8 :: dxdt(n)
  real*8 :: dxdt_ad(n)
! forcing
  real*8, intent(in) :: f
! local variables
  integer :: i
! end case
  i = n
  x_ad(i-1) = x_ad(i-1) + (x(1)-x(i-2))*dxdt_ad(i)
  x_ad(1) = x_ad(1) + x(i-1)*dxdt_ad(i)
  x_ad(i-2) = x_ad(i-2) - x(i-1)*dxdt_ad(i)
  x_ad(i) = x_ad(i) - dxdt_ad(i)
  dxdt_ad(i) = 0.0d0
  do i=n-1,3,-1
    x_ad(i-1) = x_ad(i-1) + (x(i+1)-x(i-2))*dxdt_ad(i)
    x_ad(i+1) = x_ad(i+1) + x(i-1)*dxdt_ad(i)
    x_ad(i-2) = x_ad(i-2) - x(i-1)*dxdt_ad(i)
    x_ad(i) = x_ad(i) - dxdt_ad(i)
    dxdt_ad(i) = 0.0d0
  end do
! initial cases
  i = 2
  x_ad(i-1) = x_ad(i-1) + (x(i+1)-x(n))*dxdt_ad(i)
  x_ad(i+1) = x_ad(i+1) + x(i-1)*dxdt_ad(i)
  x_ad(n) = x_ad(n) - x(i-1)*dxdt_ad(i)
  x_ad(i) = x_ad(i) - dxdt_ad(i)
  dxdt_ad(i) = 0.0d0
  i = 1
  x_ad(n) = x_ad(n) + (x(i+1)-x(n-1))*dxdt_ad(i)
  x_ad(i+1) = x_ad(i+1) + x(n)*dxdt_ad(i)
  x_ad(n-1) = x_ad(n-1) - x(n)*dxdt_ad(i)
  x_ad(i) = x_ad(i) - dxdt_ad(i)
end subroutine calc_dxdt_ad_lorenz96

