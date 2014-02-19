PROGRAM SuzukiTrotter
  implicit none
  integer, parameter :: M=2,N=4,mcSteps=10000
  double complex, parameter :: dtau=(0,1D-2)
  double precision, parameter :: J=1.d0
  double precision :: An
  integer :: i,k
  INTEGER, PARAMETER :: VecLen = 2**N
  real(kind=8) :: Spin
  double complex :: Ham_sys(2,2,2,2),Ham_int(2,2,2,2)
  double complex :: psi_sys(VecLen),psi_cen(VecLen), PsiNew(VecLen)
  OPEN(20,FILE='timeevo.dat')
  OPEN(30,FILE='spinevo.dat')
!  OPEN(40,FILE='timeevoCEN.dat')

  CALL Initialise(Ham_sys,Ham_int,psi_sys)
    
  do i=1,mcSteps

     CALL CalcSpin(psiNew,Spin)
     CALL UpdateSystem(psi_sys,psiNew,Ham_sys)
     psi_sys = psiNew
     WRITE(8,'(8F8.3)') DBLE(PsiNew(1:8))
!     CALL UpdateInteraction(psi_sys)

     print *, i
     CALL Mapiltonian(psi_sys,psi_cen)
    
!     WRITE(20,*) i,REAL(CONJG(Psi_sys)*Psi_sys)
     WRITE(30,*) Spin
!     WRITE(40,*) i,REAL(CONJG(Psi_cen)*Psi_cen)
  end do
  
CONTAINS
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE Initialise(Ham_sys,Ham_int,psi_sys)
    double complex, intent(out) :: Ham_sys(2,2,2,2)
    double complex, intent(out) :: Ham_int(2,2,2,2)
    double complex :: psi_sys(VecLen)
       
    psi_sys=(0.d0,0.D0)
   ! psi_sys(1)=(1.d0,0)!initialization of our state
    psi_sys(7)=(1.d0,0.D0)
    psi_sys=1.d0/SQRT(REAL(DOT_PRODUCT(psi_sys,psi_sys)))*psi_sys !normalization

    !Ham_sysiltonian Init
   
    Ham_sys=(0.d0,0.D0)
    Ham_sys(2,2,2,2)=exp(-dtau*J*0.25D0)
    Ham_sys(1,1,1,1)=Ham_sys(2,2,2,2)
    Ham_sys(1,2,1,2)=cos(abs(dtau)*J*0.5D0)*exp(dtau*J*0.25D0)
    Ham_sys(2,1,2,1)=Ham_sys(1,2,1,2)
    Ham_sys(2,1,1,2)=(0.D0,1.D0)*sin(abs(dtau)*J*0.5D0)*exp(dtau*J*0.25D0)
    Ham_sys(1,2,2,1)=Ham_sys(2,1,1,2)

    !Ham_Intiltonian
   
    Ham_int=0.d0
!!$    Ham_int(2,2,2,2)=exp(-dtau*An/4)
!!$    Ham_int(1,1,1,1)=Ham_int(2,2,2,2)
!!$    Ham_int(1,2,1,2)=cos(abs(dtau)*An/2)*exp(dtau*An/4)
!!$    Ham_int(2,1,2,1)=Ham_int(1,2,1,2)
!!$    Ham_int(2,1,1,2)=(0,1)*sin(abs(dtau)*An/2)*exp(dtau*An/4)
!!$    Ham_int(1,2,2,1)=Ham_int(2,1,1,2)  

   
  END SUBROUTINE Initialise
  !0000000000000000000000000000000000000000000000  
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE UpdateSystem(psiOld,PsiNew,Ham_sys)

    double complex, intent(inout) :: psiOld(VecLen), PsiNew(VecLen)
    double complex, intent(in) :: Ham_sys(2,2,2,2)
    integer :: i,j,b1,b2,b3,b4,b5,c1,c2,d2,d3,e3,e4,f1,f4
    
    PsiNew = 0.D0
    do i=1,VecLen

       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)
       b5=0!MOD((i-1)/16,2)




       !H12
       do c1=0,1
          do c2=0,1
             j=(16*b5)+(8*b4)+(4*b3)+(2*c2)+(1*c1)+1 
             psiNew(j)=psiNew(j)+Ham_sys(c2+1,c1+1,b2+1,b1+1)*psiOld(i)

          enddo
       enddo
     end do
       print *, '1',  SUM(PsiNew*CONJG(PsiNew))
       PsiOld = PsiNew
       PsiNew = 0.D0


       !H34
     DO I=1, VecLen
       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)
       b5=0!MOD((i-1)/16,2)
        do e3=0,1
          do e4=0,1
             j=(16*b5)+(8*e4)+(4*e3)+(2*b2)+(1*b1)+1
             psiNew(j)=psiNew(j)+Ham_sys(e4+1,e3+1,b4+1,b3+1)*psiOld(i)

          enddo
       enddo
     end do

      print *, '3',  SUM(PsiNew*CONJG(PsiNew))
       PsiOld = PsiNew
       PsiNew = 0.D0


       !H23
     do i=1, VecLen
       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)
       b5=0!MOD((i-1)/16,2)
        do d2=0,1
          do d3=0,1
             j=(16*b5)+(8*b4)+(4*d3)+(2*d2)+(1*b1)+1
             psiNew(j)=psiNew(j)+Ham_sys(d3+1,d2+1,b3+1,b2+1)*psiOld(i)

          enddo
       enddo
     END DO
      print *, '2',  SUM(PsiNew*CONJG(PsiNew))
       PsiOld = PsiNew
       PsiNew = 0.D0



       !H14
      DO I=1, VecLen
       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)
       b5=0!MOD((i-1)/16,2)
        do f1=0,1
          do f4=0,1
             j=(16*b5)+(8*f4)+(4*b3)+(2*b2)+(1*f1)+1
             psiNew(j)=psiNew(j)+Ham_sys(f4+1,f1+1,b4+1,b1+1)*psiOld(i)

          enddo
       enddo
     end do




      
      print *, '4',  SUM(PsiNew*CONJG(PsiNew))
      PsiOld = PsiNew

   
    
  END SUBROUTINE UpdateSystem
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE UpdateInteraction(psi)
    double complex, intent(inout) :: psi(VecLen)
    !complex, intent(in) :: Ham_int(2,2,2,2)
    double complex :: Ham
    integer :: i,j,b1,b2,b3,b4,b5,c1,c2,c3,c4,c5
    double precision :: An
    
    
    do i=1,VecLen
       
       b1=MOD((i-1),2)
       b2=MOD((i-1)/2,2)
       b3=MOD((i-1)/4,2)
       b4=MOD((i-1)/8,2)
       b5=MOD((i-1)/16,2)
       
       An=-RAND()*0.5d0
       

       !H15



       do c1=0,1
          do c5=0,1

             j=(16*c5)+(8*b4)+(4*b3)+(2*b2)+(1*c1)+1 
             psi(j)=psi(j)+Ham(c5+1,c1+1,b5+1,b1+1,An)*psi(i)
             psi=1d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo






       !H23
       do c2=0,1
          do c5=0,1
             j=(16*c5)+(8*b4)+(4*b3)+(2*c2)+(1*b1)+1
             psi(j)=psi(j)+Ham(c5+1,c2+1,b5+1,b2+1,An)*psi(i)
             psi=1.d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       !H34
       do c3=0,1
          do c5=0,1
             j=(16*c5)+(8*b4)+(4*c3)+(2*b2)+(1*b1)+1 
             psi(j)=psi(j)+Ham(c5+1,c3+1,b5+1,b3+1,An)*psi(i)
             psi=1.d0/SQRT(REAL(DOT_PRODUCT(psi,psi)))*psi !renormalization
          enddo
       enddo
       do c4=0,1
          do c5=0,1
             j=(16*c5)+(8*c4)+(4*b3)+(2*b2)+(1*b1)+1 
             psi(j)=psi(j)+Ham(c5+1,c4+1,b5+1,b4+1,An)*psi(i)
             WRITE(*,*) Ham(c5+1,c4+1,b5+1,b4+1,An)
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
    double complex, intent(inout) :: Psi_sys(VecLen)
    integer :: i
    real(kind=8), intent(out) :: Spin
    
    Spin=0.d0
    do i=1,VecLen
       if(MOD((i-1),2).eq.0) then
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
    double complex, intent(inout) :: psi_sys(VecLen)
    double complex, intent(out) :: psi_cen(16)
    integer :: i

    do i=1,16
       psi_cen(i)=SQRT(REAL(CONJG(Psi_sys(i))*Psi_sys(i))+REAL(CONJG(Psi_sys(i+16))*Psi_sys(i+16)))
    enddo
    
  END SUBROUTINE Mapiltonian


!!$  COMPLEX FUNCTION Ham(i1,i2,i3,i4,An)
!!$    double complex :: Ham
!!$    double precision, intent(in) :: An
!!$    integer, intent(in) :: i1,i2,i3,i4
!!$    
!!$    Ham=0.D0
!!$    if(i1.eq.i2.and.i1.eq.i3.and.i1.eq.i4)then
!!$       Ham=exp(-dtau*An/4)
!!$    end if
!!$    
!!$    if(i1.eq.i3.and.i2.eq.i4.and.i1.ne.i2)then
!!$       Ham=cos(abs(dtau)*An/2)*exp(dtau*An/4)
!!$    end if
!!$
!!$    if(i1.eq.i4.and.i3.eq.i2.and.i1.ne.i2)then
!!$       Ham=(0,1)*sin(abs(dtau)*An/2)*exp(dtau*An/4)
!!$    end if
!!$       
!!$
!!$  END FUNCTION Ham
    
SUBROUTINE Interactonian(P1,P2,PsiOld,PsiNew)
  integer, intent(in) :: P1, P2
  double complex, intent(inout):: PsiOld(VecLen),PsiNew(VecLen)
  integer :: i,j,b1,b2,b3,b4,b5,c1,c2,d2,d3,e3,e4,f1,f4
  integer :: BitValue,BitArray(N)

  PsiNew = 0.D0
  DO i=1,Veclen
     
     BitValue=i-1
     
     do k=1,N
        BitArray(k)=MOD(BitValue,2)
        BitValue=BitValue/2
     enddo
     
     
     do Q1=0,1
        do Q2=0,1
           
           BitValue2=BitValue+2**(P1-1)*(Q1-BitArray(P1))+2**(P2-1)*(Q2-BitArray(P2))    
           psiNew(BitValue2)=psiNew(BitValue2)+Ham(Q1+1,Q2+1,BitArray(P1)+1,BitArray(P2)+1)*psiOld(i)
           
        enddo
     enddo
  enddo
  PsiOld = PsiNew


END SUBROUTINE INTERACTONIAN       

END PROGRAM SuzukiTrotter





