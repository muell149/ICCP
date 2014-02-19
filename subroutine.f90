PROGRAM FAKE
implicit none

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
  
END PROGRAM FAKE
