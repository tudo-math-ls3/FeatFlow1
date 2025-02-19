***********************************************************************
      SUBROUTINE TCALC (MFILE,TMSTEP)
***********************************************************************
C
C	Calculates the temperature via
C	a Boussinesq approximation
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
C-----------------------------------------------------------------------
      EXTERNAL COEFFN
      EXTERNAL PARX,PARY,TMAX,UE
      EXTERNAL E030,E031,EM30,EM31
c
      CALL ZTIME (TTANFA)
      RESI=RE
      RE=1.D0/ALFA
      ISETLV=2
      ILEV=NLEV
      CALL SETLEV(ISETLV)
C
C
      FAKTOR=ALFA*RESI
      THSTEP=THSTEP*FAKTOR
      CALL XMDF1T(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,
     *            KFT,KT,NA,NU,THSTEP,TSTEPH)
C
      THSTEP=THSTEP/FAKTOR
C
      IF (THSTEP.NE.0D0) THEN
C
       THSTEP=-THSTEP
       IF (IUPW.EQ.1) THEN
        CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *              1D0,0D0,DWORK(KT),DWORK(KT),
     *              DWORK(KFT),DWORK(KFT),VWORK(KA1),KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),2,1)
       ELSE
        IF (IELT.EQ.0) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E031,COEFFN,2,1,-1D0)
        IF (IELT.EQ.1) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E030,COEFFN,2,1,-1D0)
        IF (IELT.EQ.2) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM31,COEFFN,2,1,-1D0)
        IF (IELT.EQ.3) THEN
         CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM30,COEFFN,2,1,-1D0)
        ENDIF
       ENDIF
       THSTEP=-THSTEP
      ELSE
       IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
        TOSTEP=THSTEP
        THSTEP=1D0
        IF (IELT.EQ.0) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E031,COEFFN,2,1,0D0)
        IF (IELT.EQ.1) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E030,COEFFN,2,1,0D0)
        IF (IELT.EQ.2) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM31,COEFFN,2,1,0D0)
        IF (IELT.EQ.3) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KT),DWORK(KT),
     *               DWORK(KFT),DWORK(KFT),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM30,COEFFN,2,1,0D0)
       ENDIF
       THSTEP=TOSTEP
      ENDIF
C
       CALL BDSETT (DWORK(KFT),DWORK(KFT),
     *              KWORK(L(KLMBD(ILEV))),DWORK(L(KLDBD(ILEV))),
     *              KWORK(L(LVERT)),KWORK(L(LMID)),
     *              KWORK(L(LNPR)),DWORK(L(LCORVG)),
     *              KNMBD(ILEV),PARX,PARY,UE,0D0)
C
c
C	The vector  KFT now contains the source term of the heat eq,
C       Stokes, Mass and concectif-term of the former time step
C	modified by THSTEP (depending on the time-integration 
C       scheme used). 
C	The boundary values are set!
C 	Now THSTEP is redefined from
C	[\theta-1] k   ==>  \theta k
C
C-----------------------------------------------------------------------
C
      THSTEP=TMSTEP
      CALL ZTIME(TTIM1)
      TTILIN=TTILIN+(TTIM1-TTANFA)
c
      CALL TDEF (MFILE,MSHOW,BSTOP,BNLEND)
C
C-----------------------------------------------------------------------
C
      RE=RESI
      CALL ZTIME (TTENDE)
      TIMTAL=TIMTAL+(TTENDE-TTANFA)

     
99999 END
