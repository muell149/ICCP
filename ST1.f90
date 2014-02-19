PROGRAM SuzukiTrotter
  implicit none
  integer, parameter :: M=6,N=2,mcSteps=15000
  double complex, parameter :: dtau=(0,1D-1)
  double precision, parameter :: Jo=16.d0
  double precision :: An
  integer :: i,k,l
  INTEGER, PARAMETER :: VecLen = 2**(N+M)
  real(kind=8) :: Spin
  double complex :: psiARR(100,VecLen)
  real(kind=8) :: SPINarr(100)
  double complex :: psiOld(VecLen),psi_cen(VecLen), PsiNew(VecLen)
  OPEN(20,FILE='timeevo.dat')
  OPEN(30,FILE='spinevo.dat')
  OPEN(40,FILE='AVGstate.dat')
  OPEN(50,FILE='AVGspin.dat')

  CALL Initialise(psiOld)
  L=0
  do i=1,mcSteps
     !Print *, i
     
     CALL CalcSpin(psiOld,Spin)
     CALL UpdateSystem(psiOld,psiNew)
     !CALL UpdateBath(psiOld,psiNew)
     L=L+1
     
     !PSIarr(L,:)=REAL(CONJG(PsiNew)*PsiNew)
     !SPINarr(L)=Spin

     !WRITE(8,'(8F8.3)') DBLE(PsiNew(81:88))
     !WRITE(20,*) i,REAL(CONJG(PsiNew(61:93))*PsiNew(61:93))
     WRITE(30,*) i,Spin

     if(L.eq.10)then
        WRITE(40,*) i,REAL(CONJG(PsiNew)*PsiNew)!SUM(psiARR)/SIZE(psiARR)
        WRITE(50,*) Spin!SUM(SPINarr)/SIZE(SPINarr)
        L=0
     end if

     

  end do
  
CONTAINS
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE Initialise(psiOld)
    double complex :: psiold(VecLen)
       
    psiold=(0.d0,0.D0)
!!$    do i=1,VecLen
!!$       psiold(i)=RAND()*(1.d0,0.d0)
!!$    end do
    psiold(86)=(1.d0,0.D0)
 !   psiold(19)=(1.d0,0.D0)
    psiold=1.d0/SQRT(REAL(DOT_PRODUCT(psiold,psiold)))*psiold !normalization
   
  END SUBROUTINE Initialise
  !0000000000000000000000000000000000000000000000  
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE UpdateSystem(psiOld,PsiNew)
    double complex, intent(inout) :: psiOld(VecLen), PsiNew(VecLen)
    integer :: i,j,b1,b2,b3,b4,b5,c1,c2,d2,d3,e3,e4,f1,f4
    
     CALL Interactonian(1,2,PsiOld,PsiNew,Jo)
    ! print *, '1',  SUM(PsiNew*CONJG(PsiNew))
     PsiOld=PsiNew

!!$     CALL Interactonian(3,4,PsiOld,PsiNew,Jo)
!!$    ! print *, '2',  SUM(PsiNew*CONJG(PsiNew))
!!$     PsiOld=PsiNew
!!$
!!$     CALL Interactonian(2,3,PsiOld,PsiNew,Jo)
!!$    ! print *, '3',  SUM(PsiNew*CONJG(PsiNew))
!!$     PsiOld=PsiNew
!!$
!!$     CALL Interactonian(1,4,PsiOld,PsiNew,Jo)
!!$    ! print *, '4',  SUM(PsiNew*CONJG(PsiNew))
!!$     PsiOld=PsiNew
!!$   
    
  END SUBROUTINE UpdateSystem
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  SUBROUTINE UpdateBath(psiOld,PsiNew)
    double complex, intent(inout) :: psiOld(VecLen), PsiNew(VecLen)
    integer :: i,j
        
    do j=N+1,N+M
       do i=1,N
          An=-RAND()*0.5D0
          CALL Interactonian(i,j,PsiOld,PsiNew,An)
          ! print *, 'bath',  SUM(PsiNew*CONJG(PsiNew))
          PsiOld=PsiNew
       end do
    enddo
   
    
  END SUBROUTINE UpdateBath
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
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000

  DOUBLE COMPLEX FUNCTION Him(i1,i2,i3,i4,An)
    !double complex :: Him
    double precision, intent(in) :: An
    integer, intent(in) :: i1,i2,i3,i4
    
    Him=(0.d0,0.d0)

    if(i1.eq.i2.and.i1.eq.i3.and.i1.eq.i4)then
       Him=exp(-dtau*An/4)
    end if
    
    if(i1.eq.i3.and.i2.eq.i4.and.i1.ne.i2)then
       Him=cos(abs(dtau)*An/2)*exp(dtau*An/4)
    end if

    if(i1.eq.i4.and.i3.eq.i2.and.i1.ne.i2)then
       Him=(0,1)*sin(abs(dtau)*An/2)*exp(dtau*An/4)
    end if
       
    RETURN

  END FUNCTION Him
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000
  !0000000000000000000000000000000000000000000000   
SUBROUTINE Interactonian(P1,P2,PsiOld,PsiNew,An)
  integer, intent(in) :: P1, P2
  double complex, intent(inout):: PsiNew(VecLen)
  integer :: i,k,q1,q2
  integer :: BitValue,BitValue2,BitArray(N+M)
  real(kind=8), intent(in) :: An
  double complex, intent(inout) :: PsiOld(VecLen)

  PsiNew = 0.D0
  DO i=1,Veclen
     
     BitValue=i-1
     
     do k=1,(N+M)
        BitArray(k)=MOD(BitValue,2)
        BitValue=BitValue/2
     enddo
        !  
     do Q1=0,1
        do Q2=0,1
           BitValue2=I+2**(P1-1)*(Q1-BitArray(P1))+2**(P2-1)*(Q2-BitArray(P2))    
           psiNew(BitValue2)=psiNew(BitValue2)+Him(Q1+1,Q2+1,BitArray(P1)+1,BitArray(P2)+1,An)*psiOld(i)
        enddo
     enddo

  enddo

END SUBROUTINE INTERACTONIAN       

END PROGRAM SuzukiTrotter
