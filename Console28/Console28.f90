!  Console28.f90 
!
!  FUNCTIONS:
!  Console28 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Console28
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Console28
    parameter nc=50, ny=100, pi=3.141, vp=0.5, vm=-0.5, dy=0.01, sig=10
    real*4 cx(-nc:nc), cy(-nc:nc), fi(-nc:nc,-nc:nc,0:ny), psi(-nc:nc,-nc:nc,0:ny), n(0:ny), T(0:ny), u(0:ny), y(0:ny), nwp, nwm, Twp, Twm, s
    integer*4 i, j, k
    
    do i=-nc, nc
    cx(i)=0.1*i
    write(*,*) cx(i)
        do j=-nc, nc
        cy(j)=0.1*j
            do k=0, ny
            n(k)=1
            T(k)=1
            u(k)=0
            fi(i,j,k)=(n(k)/(pi*T(k)))*exp((-((cx(i)-u(k))**2)-(cy(j)**2))/T(k))
            psi(i,j,k)=(T(k)/2)*(n(k)/(pi*T(k)))*exp((-((cx(i)-u(k))**2)-(cy(j)**2))/T(k))
            enddo
        enddo
    enddo

    nwp=1
    nwm=1
    Twp=1.45
    Twm=0.55
    write(*,*) cy(1), cy(-1)
   do l=1, 200
    do i=-nc, nc
        do j=0, nc
        fi(i,j,0)=(nwp/(pi*Twp))*exp((-((cx(i)-vp)**2)-(cy(j)**2))/Twp)
        psi(i,j,0)=(Twp/2)*(nwp/(pi*Twp))*exp((-((cx(i)-vp)**2)-(cy(j)**2))/Twp)
        fi(i,-j,ny)=(nwm/(pi*Twm))*exp((-((cx(i)-vm)**2)-(cy(-j)**2))/Twm)
        psi(i,-j,ny)=(Twm/2)*(nwm/(pi*Twm))*exp((-((cx(i)-vm)**2)-(cy(-j)**2))/Twm)
        enddo
    enddo
    
    do i=-nc, nc
        do j=0, nc
            do k=1, ny
            fi(i,j,k)=(cy(j)*fi(i,j,k-1)+dy*sig*((n(k)**2)/(pi*sqrt(T(k))))*exp((-((cx(i)-u(k))**2)-(cy(j)**2))/T(k)))/(cy(j)+dy*sig*n(k)*sqrt(T(k)))
            psi(i,j,k)=(cy(j)*psi(i,j,k-1)+dy*sig*(((n(k)**2)*sqrt(T(k)))/(2*pi))*exp((-((cx(i)-u(k))**2)-(cy(j)**2))/T(k)))/(cy(j)+dy*sig*n(k)*sqrt(T(k)))
            enddo
        enddo
   enddo
  
   do i=-nc, nc
        do j=0, nc
            do k=ny-1, 0, -1
            fi(i,-j,k)=(cy(-j)*fi(i,-j,k+1)+(-1)*dy*sig*((n(k)**2)/(pi*sqrt(T(k))))*exp((-((cx(i)-u(k))**2)-(cy(-j)**2))/T(k)))/(cy(-j)+(-1)*dy*sig*n(k)*sqrt(T(k)))
            psi(i,-j,k)=(cy(-j)*psi(i,-j,k+1)+(-1)*dy*sig*(((n(k)**2)*sqrt(T(k)))/(2*pi))*exp((-((cx(i)-u(k))**2)-(cy(-j)**2))/T(k)))/(cy(-j)+(-1)*dy*sig*n(k)*sqrt(T(k)))
            enddo
        enddo
   enddo
   
   do k=0, ny
   n(k)=0
   u(k)=0
   T(k)=0
   y(k)=0
   enddo

   do i=-nc, nc
        do j=-nc, nc
            do k=0, ny
            n(k)=n(k)+0.01*fi(i,j,k)
            !n(k)=n(k)+0.01*(fi(i,j,k)+fi(i+1,j+1,k))/2            
            enddo
        enddo
   enddo
   do i=-nc, nc
        do j=-nc, nc
            do k=0, ny
            u(k)=u(k)+0.01*((1/n(k)))*cx(i)*fi(i,j,k)
            y(k)=y(k)+0.01*((1/n(k)))*cy(j)*fi(i,j,k)
            enddo
        enddo
   enddo
   do i=-nc, nc
        do j=-nc, nc
            do k=0, ny
            T(k)=T(k)+0.01*(2/(3*n(k)))*(psi(i,j,k)+fi(i,j,k)*(((cx(i)-u(k))**2)+(cy(j)**2)))
            enddo
        enddo
   enddo
   nwp=0
   nwm=0
   do i=-nc, nc
        do j=1, nc
        nwp=nwp+2*sqrt(pi)*0.01*cy(j)*fi(i,-j,0)
        enddo
   enddo
   do i=-nc, nc
        do j=1, nc
        nwm=nwm+2*sqrt(pi)*0.01*cy(j)*fi(i,j,ny)
        enddo
   enddo
   s=0
   do k=0, ny
   s=s+n(k)*0.01
   enddo

   do k=0, ny
   n(k)=n(k)/s
   enddo

   write(*,*) s
   
enddo   
       
open(14,file='c.txt',status='unknown')
    
!    do i=-nc, nc
!    do j=-nc,nc
                do k=0, ny
          write(14,*)  '(',dy*k,';',u(k),')'!'значение ',i,j,k,'=',  u(k), T(k)
          enddo
!          enddo
!          enddo
  
    end program Console28

