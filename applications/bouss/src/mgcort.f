************************************************************************
      SUBROUTINE MGCORT (NU)
************************************************************************
C
C	Zweck der Uebeung:
C	fuehrt die MG-Korrektur fuer die Temperatur durch.
C	z.Z. wird OMEGAT noch nicht adaptiv bestimmt,
C	sondern aus dem Datenfile eingelesen
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      INCLUDE 'bouss.inc'
      SAVE
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
c
c     do 3 I=1,nu
c3      write (*,*) 'mgcor',DWORK(L(LTOLD)-1+i) , DWORK(KT-1+i)
      AA1=-OMEGAT
      AA2= OMEGAT
      DO 1 I=1,NU
1     DWORK(KT-1+I)=AA1*DWORK(L(LTOLD)-1+I)+AA2*DWORK(KT-1+I)
C
C=======================================================================
C *** Update the solution   T:=TOLD+T
C=======================================================================
C
        AA1= 1.D0
        AA2= 1.D0
c        CALL  LLC1 (DWORK(KTOLD),DWORK(KT),NU,AA1,AA2)
      DO 2 I=1,NU
2     DWORK(KT-1+I)=AA1*DWORK(L(LTOLD)-1+I)+AA2*DWORK(KT-1+I)
C
      END
