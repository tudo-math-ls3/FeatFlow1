      PROGRAM PDIFF
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      OPEN (UNIT=50,FILE='p1_0.1_1d-3_test2')
      OPEN (UNIT=51,FILE='p2_0.1_1d-3_test2')
      WRITE(6,*) 'NNP= ',NNP
      READ (5,*) NNP
C
      OPEN (UNIT=52,FILE='pdiff_0.1_1d-3_test2')
C
C
      DO 10 IP=1,NNP
      READ (50,*) VX,VP1
      READ (51,*) VX,VP2
      write(6,*) VX,VP1,VP2
      WRITE(52,*) VX,VP1-VP2
10    CONTINUE
C
      REWIND(52)
      CLOSE(52)
C
C
      END
