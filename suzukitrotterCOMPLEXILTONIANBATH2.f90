PROGRAM SuzukiTrotter
  implicit none
  integer, parameter :: M=2,N=4,mcSteps=1500
  complex, parameter :: dtau=(0,1E-2)
  integer :: i,k
  real(kind=8) :: Spin, dumm
  complex :: Ham_sys(2,2,2,2),Ham_int(2,2,2,2)
  complex :: psi_sys(32),psi_cen(16) 
  OPEN(20,FILE='timeevo.dat')
  OPEN(30,FILE='spinevo.dat')
  OPEN(40,FILE='timeevoCEN.dat')

  CALL Initialise(Ham_sys,Ham_int,psi_sys)
    
  do i=1,mcSteps
    
     CALL CalcSpin(psi_sys,Spin)
     CALL UpdateSystem(psi_sys,Ham_sys)
     CALL UpdateInteraction(psi_sys,Ham_int)
     CALL Mapiltonian(psi_sys,psi_cen)
    
     WRITE(20,*) i,REAL(CONJG(Psi_sys)*Psi_sys)
     WRITE(30,*) i,Spin
     WRITE(40,*) i,REAL(CONJG(Psi_cen)*Psi_cen)
  end do
  
CONTAINS
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE Initialise(Ham_sys,Ham_int,psi_sys)
    complex, intent(out) :: Ham_sys(2,2,2,2)
    complex, intent(out) :: Ham_int(2,2,2,2)
    complex :: psi_sys(32)
    integer :: J
    real(kind=8) :: An
    
    psi_sys=(0.d0,0)
    psi_sys(11)=(1.d0,0)!initialization of our state
    psi_sys(27)=(1.d0,0)!initialization of our state
    psi_sys=1.d0/SQRT(REAL(DOT_PRODUCT(psi_sys,psi_sys)))*psi_sys !normalization

    !Ham_sysiltonian Init
    J=20.d0
    Ham_sys=(0.d0,0)
    Ham_sys(2,2,2,2)=exp(-dtau*J/4)
    Ham_sys(1,1,1,1)=Ham_sys(2,2,2,2)
    Ham_sys(1,2,1,2)=cos(abs(dtau)*J/2)*exp(dtau*J/4)
    Ham_sys(2,1,2,1)=Ham_sys(1,2,1,2)
    Ham_sys(2,1,1,2)=(0,1)*sin(abs(dtau)*J/2)*exp(dtau*J/4)
    Ham_sys(1,2,2,1)=Ham_sys(2,1,1,2)

    !Ham_Intiltonian
    An=-0.25d0!-RAND()/2d0
    WRITE(*,*) An
    Ham_int=0.d0
    Ham_int(2,2,2,2)=exp(-dtau*An/4)
    Ham_int(1,1,1,1)=Ham_int(2,2,2,2)
    Ham_int(1,2,1,2)=cos(abs(dtau)*An/2)*exp(dtau*An/4)
    Ham_int(2,1,2,1)=Ham_int(1,2,1,2)
    Ham_int(2,1,1,2)=(0,1)*sin(abs(dtau)*An/2)*exp(dtau*An/4)
    Ham_int(1,2,2,1)=Ham_int(2,1,1,2)  

   
  END SUBROUTINE Initialise
  !0000000000000000000000000000000000000000000000  
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE UpdateSystem(psi,Ham_sys)
    complex, intent(inout) :: psi(32)
    complex, intent(in) :: Ham_sys(2,2,2,2)
    integer :: i,j,b1,b2,b3,b4,c1,c2,d2,d3,e3,e4
    
    do i=1,32

       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)

       !H12
       do c1=0,1
          do c2=0,1
             j=(8*b4)+(4*b3)+(2*c2)+(1*c1)+1 
             psi(j)=psi(j)+Ham_sys(c2+1,c1+1,b2+1,b1+1)*psi(i)
             psi=1d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       !H23
       do d2=0,1
          do d3=0,1
             j=(8*b4)+(4*d3)+(2*d2)+(1*b1)+1
             psi(j)=psi(j)+Ham_sys(d3+1,d2+1,b3+1,b2+1)*psi(i)
             psi=1d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       !H34
       do e3=0,1
          do e4=0,1
             j=(8*e4)+(4*e3)+(2*b2)+(1*b1)+1
             psi(j)=psi(j)+Ham_sys(e4+1,e3+1,b4+1,b3+1)*psi(i)
             psi=1d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
    enddo
    
  END SUBROUTINE UpdateSystem
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE UpdateInteraction(psi,Ham_int)
    complex, intent(inout) :: psi(32)
    complex, intent(in) :: Ham_int(2,2,2,2)
    integer :: i,j,b1,b2,b3,b4,b5,c1,c2,c3,c4,c5
    
    do i=1,32
       
       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)
       b5=MOD((i-1)/16,2)

       !H15
       do c1=0,1
          do c5=0,1
             j=(16*c5)+(8*b4)+(4*b3)+(2*b2)+(1*c1)+1 
             psi(j)=psi(j)+Ham_sys(c5+1,c1+1,b5+1,b1+1)*psi(i)
             psi=1d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       !H23
       do c2=0,1
          do c5=0,1
             j=(16*c5)+(8*b4)+(4*b3)+(2*c2)+(1*b1)+1
             psi(j)=psi(j)+Ham_sys(c5+1,c2+1,b5+1,b2+1)*psi(i)
             psi=1.d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       !H34
       do c3=0,1
          do c5=0,1
             j=(16*c5)+(8*b4)+(4*c3)+(2*b2)+(1*b1)+1 
             psi(j)=psi(j)+Ham_sys(c5+1,c3+1,b5+1,b3+1)*psi(i)
             psi=1.d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       do c4=0,1
          do c5=0,1
             j=(16*c5)+(8*c4)+(4*b3)+(2*b2)+(1*b1)+1 
             psi(j)=psi(j)+Ham_sys(c5+1,c4+1,b5+1,b4+1)*psi(i)
             psi=1d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
    enddo
    
  END SUBROUTINE UPDATEINTERACTION
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE CalcSpin(Psi_sys,Spin)
    complex, intent(inout) :: Psi_sys(32)
    integer :: i
    real(kind=8), intent(out) :: Spin
    
    Spin=0.d0
    do i=1,32
       if(MOD(i-1,2).eq.0) then
          Spin=Spin-0.5d0*REAL(CONJG(Psi_sys(i))*Psi_sys(i))
          else
          Spin=Spin+0.5d0*REAL(CONJG(Psi_sys(i))*Psi_sys(i))   
       endif
    enddo


  END SUBROUTINE CALCSPIN
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000  
  SUBROUTINE Mapiltonian(psi_sys,psi_cen)
    complex, intent(inout) :: psi_sys(32)
    complex, intent(out) :: psi_cen(16)
    integer :: i

    do i=1,16
       psi_cen(i)=REAL(CONJG(Psi_sys(i))*Psi_sys(i))+REAL(CONJG(Psi_sys(i+16))*Psi_sys(i+16))
    enddo
    
  END SUBROUTINE Mapiltonian
END PROGRAM SuzukiTrotter
