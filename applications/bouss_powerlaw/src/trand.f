***********************************************************
      SUBROUTINE TRAND (KNPR,DDBD,NVT,NMBD,BDFLAG,IMID,IMBD,
     *                 IVAR,PX,PY)
***********************************************************
C
C	DEFINES THE NEUMANN PART OF THE BOUNDARY
C	FOR THE TEMPERATURE. 
C 	BDFLAG=.FALSE.  MEANS ON THIS VERTEX ARE
C	NEUMANN BOUNDARY CONDITIONS IMPOSED
C
C       AUTOMATICALLY IF VELO HAS SET NEUMANN CONDITIONS!!!
C
***********************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      DIMENSION DDBD(*),KNPR(*)
      EXTERNAL PARX,PARY,UE
      SAVE
C
C
      BDFLAG=.FALSE.
C
C
C
      IF (IVAR.EQ.1) THEN
         INPR=KNPR(KNPR(NVT+IMID))
         DPAR=DDBD(IMBD)
C----------------------------------------------------
         IF (INPR.GE.1) THEN
           PX=PARX(DPAR,INPR)
           PY=PARY(DPAR,INPR)
           IF (PX.EQ.0D0) THEN
             BDFLAG=.TRUE.
           ENDIF
           IF (INPR.EQ.2) THEN
             BDFLAG=.TRUE.
           ENDIF
         ENDIF
C----------------------------------------------------
      ENDIF
C
C
      IF (IVAR.NE.1) THEN

C----------------------------------------------------
        DIST=SQRT((PX-0.2D0)**2+(PY-0.2D0)**2)
        IF (PX.EQ.0D0) THEN
          BDFLAG=.TRUE.
        ENDIF
        IF (DIST.LE.0.15D0) THEN
          BDFLAG=.TRUE.
        ENDIF
C----------------------------------------------------
C
      ENDIF
C
C
99999 END


