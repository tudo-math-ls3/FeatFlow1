************************************************************************
*    M011   (edited from M010)                                         *
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* M011                                                                 *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          multigrid iteration                                         *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LL21 , LLC1                            *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting address of vectors containing the           *
* DB       R*8    solution and the right hand side, DD is used as      *
* DD       R*8    auxiliary vector only                                *
* KOFFX           The actual starting address of DX on level ILEV      *
* KOFFB           is DX(1+KOFFX(ILEV)) (analogously for DB and DD)     *
* KOFFD           Total space required for all vectors is              *
*                 KNEQ(NLMIN)+...+KNEQ(NLMAX)                          *
*                 DX(1+KOFFX(NLMAX)) contains initial solution         *
*                 DB(1+KOFFB(NLMAX)) contains right hand side          *
* KNEQ     I*4    Number of equations for all levels                   *
* NLMAX    I*4    Iteration uses levels NLMIN to NLMAX,                *
* NLMIN    I*4    NLMAX  is the finest level                           *
* NIT      I*4    Maximum number of iterations                         *
*                 Iteration completed after reaching the finest level  *
* EPS1      R*8    Desired precision                                   *
*                 Stop if !!DEF!! < EPS1 !!DEFOLD!!                    *
* EPS2      R*8    Desired precision                                   *
*                 Stop if !!DEF!! < EPS2                               *
* KPRSM    I*4    Number of pre -smoothing steps for all levels        *
* KPOSM    I*4    Number of post-smoothing steps for all levels        *
* ICYCLE   I*4    <0: special cycle types (not yet implemented)        *
*                 =0: F-Cycle                                          *
*                 =1: V-Cycle                                          *
*                 =2: W-Cycle                                          *
*                 >2: Cycle of higher order                            *
* DAX      SUBR   CALL DAX(DX,DAX,NEQ,A1,A2)                           *
*                 Returns  DAX := A1*A*DX+A2*DAX                       *
* DPROL    SUBR   DPROL(DX1,DX2)                                       *
*                 Returns  DX2 := Prolongation(DX1)  to higher level   *
* DSTEP    SUBR   DPROL(DX,DD,DB,DSTEPP)                               *
*                 Returns  DSTEPP := optimal step size for correction  *
* DREST    SUBR   DREST(DD1,DD2)                                       *
*                 Returns  DD2 := Restriction(DD1)  to lower level     *
* DPRSM    SUBR   DPRSM(DX,DB,DD,NEQ,NPRSM)                            *
*                 Returns  DX  after  NPRSM:=KPRSM(ILEV)               *
*                 pre-smoothing steps, DD  is used as auxiliary vector *
* DPOSM    SUBR   Same as above, used for post-smoothing               *
* DEX      SUBR   DEX(DX,DB,DD,NEQ)                                    *
*                 Returns exact solution                               *
* DBC      SUBR   DBC(DX,NEQ)                                          *
*                 Copies boundary data onto components of  DX          *
* KIT0     I*4    auxiliary vectors of length  NLMAX                   *
* KIT      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector on DX(KOFFX(NLMAX))                  *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
* RHOLMG   R*8    Multigrid convergence rate                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE  M011 (DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,
     *                  EPS1,EPS2,DAX,DPROL,DREST,DPRSM,DPOSM,DEX,DEXA,
     *                  DBC,DSTEP,KIT0,KIT,IREL,IDEFMG,RHOLMG,BMGEND)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION DX(*),DB(*),DD(*),KOFFX(*),KOFFB(*),KOFFD(*)
      DIMENSION KNEQ(*),KIT0(*),KIT(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      SAVE 
C
      SUB='M011  '
      IF (ICHECK.GE.997) CALL OTRC('M011  ','05/13/91')
      IER=0
C
C
      BMGEND=.FALSE.
      BREL=IREL.EQ.1
      BMSG2=M.GE.2.OR.MT.GE.2
      MT0=MT
      MT=0
C
      BTIME=IMTIME.GT.0
      IF (BTIME) THEN
       IF (IMTIME.EQ.1) THEN
        TTMG=0D0
        TTS=0D0
        TTE=0D0
        TTD=0D0
        TTP=0D0
        TTR=0D0
       ENDIF
       CALL ZTIME(TTMG0)
      ENDIF
C
C
      NIT0=MAX(ITE,0)
      IF (ICYCLE.LT.0) THEN
       CALL WERR(-181,'M011  ')
       GOTO 99999
      ENDIF
      IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
       CALL WERR(-182,'M011  ')
       GOTO 99999
      ENDIF
      ILEV=NLMAX
C
C *** special case - zero solution
      IF (BTIME) CALL ZTIME(TTD0)
      CALL LLI1(DB(1+KOFFB(NLMAX)),KNEQ(NLMAX),DMAX,INDMAX)
      IF (DMAX.LE.1D-12) THEN
       CALL LCP1(DB(1+KOFFB(NLMAX)),DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
       RHOLMG=0D0
       GOTO 99999
      ENDIF
C
      IF (BTIME) THEN
       CALL ZTIME(TTD1)
       TTD=TTD+TTD1-TTD0
      ENDIF
C
C *** special case - only one level
      IF (NLMIN.EQ.NLMAX) THEN
       IF (BTIME) CALL ZTIME(TTE0)
       CALL DEXA(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),RHOLMG,EPS2,ITECG)
       IF (BTIME) THEN
        CALL ZTIME(TTE1)
        TTE=TTE+TTE1-TTE0
       ENDIF
       ITE=ITECG
       GOTO 99999
      ENDIF
C
C *** level counts
      KIT0(NLMAX)=1
      DO 2 ILEV=NLMIN+1,NLMAX-1
      IF (ICYCLE.EQ.0) THEN
       KIT0(ILEV)=2
      ELSE
       KIT0(ILEV)=ICYCLE
      ENDIF
2     CONTINUE
C
CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0) 
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *            DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
C
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
       DEFOLD=DEF
C *** FD  is considered as initial defect
       FD=DEF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       IF (DEF.LE.EPS2.AND..NOT.BREL) THEN
        ITE=0
        GOTO 1000
       ENDIF
      ENDIF
C
C
C *** Start multigrid iteration
      DO 100 ITE=1,NIT
C
C *** initialize level counts for all levels
      DO 101 ILEV=NLMIN,NLMAX
101   KIT(ILEV)=KIT0(ILEV)
C
      ILEV=NLMAX
110   IF (ILEV.NE.NLMIN) THEN
C
CCC      IF (ILEV.NE.NLMAX) THEN
CCC ***  defect on finest level already available
C
        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPRSM(ILEV).GT.0) 
     *   CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPRSM(ILEV))
C
        IF (BTIME) THEN
         CALL ZTIME(TTS1)
         TTS=TTS+TTS1-TTS0
        ENDIF
C
        IF (BTIME) CALL ZTIME(TTD0)
        CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),
     *           -1D0,1D0)
C
        IF (BTIME) THEN
         CALL ZTIME(TTD1)
         TTD=TTD+TTD1-TTD0
        ENDIF
C
CCC       ENDIF
C
       ILEV=ILEV-1
C
C *** restriction of defect
       IF (BTIME) CALL ZTIME(TTR0)
       CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
C
       IF (BTIME) THEN
        CALL ZTIME(TTR1)
        TTR=TTR+TTR1-TTR0
       ENDIF
C
C ***  choose zero as initial vector on lower level
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
       CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV))
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       GOTO 110
C
      ENDIF
C
C *** exact solution on lowest level
      IF (BTIME) CALL ZTIME(TTE0)
      CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),DD(1+KOFFD(NLMIN)),
     *         KNEQ(NLMIN),RHOLMG)
C
      IF (BTIME) THEN
       CALL ZTIME(TTE1)
       TTE=TTE+TTE1-TTE0
      ENDIF
C
C
130   IF (ILEV.NE.NLMAX) THEN
       ILEV=ILEV+1
C
C *** DPROL  returns  DD:=PROL(DX)
      IF (BTIME) CALL ZTIME(TTP0)
      CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))
      CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV))
C
      CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),DB(1+KOFFB(ILEV)),
     *           DSTEPP)
      CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),DSTEPP,
     *          1D0)
C
      IF (BTIME) THEN
       CALL ZTIME(TTP1)
       TTP=TTP+TTP1-TTP0
      ENDIF
C
C
C *** Post-smoothing
      IF (BTIME) CALL ZTIME(TTS0)
      IF (KPOSM(ILEV).GT.0) 
     * CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),
     *            KNEQ(ILEV),KPOSM(ILEV))
C
      IF (BTIME) THEN
       CALL ZTIME(TTS1)
       TTS=TTS+TTS1-TTS0
      ENDIF
C
       KIT(ILEV)=KIT(ILEV)-1
       IF (KIT(ILEV).EQ.0) THEN
        IF (ICYCLE.EQ.0) THEN
         KIT(ILEV)=1
        ELSE
         KIT(ILEV)=KIT0(ILEV)
        ENDIF
        GOTO 130
       ELSE
        GOTO 110
       ENDIF
C
      ENDIF
C
C
CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0)  
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF
C
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
C
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
       IF (MT.GE.9) WRITE (6,10001) ITE,DEF
CCC       WRITE (6,*) ITE,DEF,(DEF/FD)**(1D0/DBLE(ITE))
       IF (BMSG2) THEN
        MT=MT0
        WRITE (CPARAM,'(I15,D25.16)') ITE,DEF
        CALL OMSG(73,'M011  ')
        WRITE (MPROT,10001) ITE,DEF
        WRITE (*,10001) ITE,DEF
        MT=0
       ENDIF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
C=======================================================================
C ***  Unexpected STOP !!!
C=======================================================================
       IF ((DEF/DEFOLD.GT.1D3).OR.(DEF/FD.GT.5D3)) THEN
        BMGEND=.TRUE.
        GOTO 1000
       ENDIF
C
      DEFOLD=DEF
       IF (BREL) THEN
        IF (DEF.LE.FD*EPS1.AND.DEF.LE.EPS2.AND.ITE.GE.NIT0) GOTO 1000
       ELSE
        IF (DEF.LE.EPS2) GOTO 1000
       ENDIF
      ENDIF
C
100   CONTINUE
C
      ITE=NIT
      IF (IDEFMG.EQ.1) THEN
       MT=MT0
       RHOLMG=DEF/FD
       RHOLMG=RHOLMG**(1D0/DBLE(NIT))
       IF ((RHOLMG).GE.1D0) BMGEND=.TRUE.
       IF (MT.GE.9) WRITE (6,10000) NIT,DEF,FD,DEF/FD,RHOLMG
       IF (BMSG2) THEN
        WRITE (CPARAM,'(I15,4D25.16)') NIT,DEF,FD,DEF/FD,RHOLMG
        CALL OMSG(72,'M011  ')
        WRITE (MPROT,10000) NIT,DEF,FD,DEF/FD,RHOLMG
        WRITE (MPROT,*) 'IER=1 in M011'
       ENDIF
      ELSE
       RHOLMG=DBLE(NIT)
      ENDIF
C
      GOTO 99999
C
1000  IER=0
      MT=MT0
      RHOLMG=0.D0
      IF (FD.GE.1D-70) RHOLMG=(DEF/FD)**(1D0/DBLE(ITE))
      IF (MT.GE.9) WRITE (6,'(I15,4D25.16)') ITE,DEF,FD,DEF/FD,RHOLMG
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,4D25.16)') ITE,DEF,FD,DEF/FD,RHOLMG
       CALL OMSG(72,'M011  ')
       WRITE (CPARAM,'(D25.16)')  RHOLMG
       CALL OMSG(76,'M011  ')
       WRITE (MPROT,10002) ITE,DEF,FD,DEF/FD,RHOLMG
      ENDIF
C
99999 MT=MT0 
      IF (BTIME) THEN
       CALL ZTIME(TTMG1)
       TTMG=TTMG+TTMG1-TTMG0
      ENDIF
C
10000 FORMAT (//' CONVERGENCE OF MG-METHOD FAILED'/
     *        ' NORM OF DEFECT AFTER',I6,' ITERATIONS',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' KONVERGENZRATE          ',D12.3)
10001 FORMAT (' ITERATION',I6,5X,'  RESIDUUM =',D12.3)
10002 FORMAT (' NUMBER OF ITERATIONS    ',I10/
     *        ' NORM OF DEFECT          ',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' KONVERGENZRATE          ',D12.3)
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   PROLU2  (DU1,DV1,DW1,DP1,DU2,DV2,DW2,DP2,
     *                     KAREA1,KADJ1,NAT1,NEL1,
     *                     KAREA2,KADJ2,NAT2,NEL2,IINT)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1,DW1,DP1) 
*                to the fine grid vector (DU2,DV2,DW2,DP2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1,DW1,DP1           - coarse grid vector
*      KAREA1,KADJ1,NAT1,NEL1,  - data of the coarse grid
*      KAREA2,KADJ2,NAT2,NEL2,  - data of the fine grid
*
*    Output:
*      DU2,DV2,DW2,DP2           - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NNAE=6)
      DIMENSION DU1(1),DV1(1),DW1(1),DP1(1),DU2(1),DV2(1),DW2(1),DP2(1),
     *          KAREA1(NNAE,*),KAREA2(NNAE,*),
     *          KADJ1(NNAE,*),KADJ2(NNAE,*)
      SAVE
C-----------------------------------------------------------------------
C
      IF (IINT.GT.0) THEN
       A1= 11D0/24D0
       A2=  7D0/48D0
       A3=- 5D0/48D0
       A4=- 1D0/24D0
       A5= 11D0/24D0
       A6=  1D0/12D0
       A7=- 1D0/24D0
      ELSE
       A1= 0.500D0
       A2= 0.125D0
       A3=-0.125D0
       A4= 0.000D0
       A5= 0.500D0
       A6= 0.000D0
       A7= 0.000D0
      ENDIF
C
C
C
C *** Zero initialization of (DU2,DV2,DW2)
      CALL  LCL1 (DU2,NAT2)
      CALL  LCL1 (DV2,NAT2)
      CALL  LCL1 (DW2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DUH1=DU1(KAREA1(1,IEL1))
      DUH2=DU1(KAREA1(2,IEL1))
      DUH3=DU1(KAREA1(3,IEL1))
      DUH4=DU1(KAREA1(4,IEL1))
      DUH5=DU1(KAREA1(5,IEL1))
      DUH6=DU1(KAREA1(6,IEL1))

C
      DVH1=DV1(KAREA1(1,IEL1))
      DVH2=DV1(KAREA1(2,IEL1))
      DVH3=DV1(KAREA1(3,IEL1))
      DVH4=DV1(KAREA1(4,IEL1))
      DVH5=DV1(KAREA1(5,IEL1))
      DVH6=DV1(KAREA1(6,IEL1))
C
      DWH1=DW1(KAREA1(1,IEL1))
      DWH2=DW1(KAREA1(2,IEL1))
      DWH3=DW1(KAREA1(3,IEL1))
      DWH4=DW1(KAREA1(4,IEL1))
      DWH5=DW1(KAREA1(5,IEL1))
      DWH6=DW1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
      DP2(IELH5)=DPH
      DP2(IELH6)=DPH
      DP2(IELH7)=DPH
      DP2(IELH8)=DPH
C
      B11=A1*DUH1
      B12=A2*DUH1
      B13=A3*DUH1
      B14=A4*DUH1
      B15=A5*DUH1
      B16=A6*DUH1
      B17=A7*DUH1
C
      C11=A1*DUH2
      C12=A2*DUH2
      C13=A3*DUH2
      C14=A4*DUH2
      C15=A5*DUH2
      C16=A6*DUH2
      C17=A7*DUH2
C
      D11=A1*DUH3
      D12=A2*DUH3
      D13=A3*DUH3
      D14=A4*DUH3
      D15=A5*DUH3
      D16=A6*DUH3
      D17=A7*DUH3
C
      E11=A1*DUH4
      E12=A2*DUH4
      E13=A3*DUH4
      E14=A4*DUH4
      E15=A5*DUH4
      E16=A6*DUH4
      E17=A7*DUH4
C
      F11=A1*DUH5
      F12=A2*DUH5
      F13=A3*DUH5
      F14=A4*DUH5
      F15=A5*DUH5
      F16=A6*DUH5
      F17=A7*DUH5
C
      G11=A1*DUH6
      G12=A2*DUH6
      G13=A3*DUH6
      G14=A4*DUH6
      G15=A5*DUH6
      G16=A6*DUH6
      G17=A7*DUH6
      R1=B11+G14
      R2=C11+E14
      R3=D11+F14
      R4=C14+E11
      R5=F11+D14
      R6=B14+G11
      R7=B15+G17
      R8=B16+G16
      R9=B17+G15
C
C
C
C
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DU2(KAREA2(1,IELH1))=DU2(KAREA2(1,IELH1))+
     *                      (R1+C12+D13+
     *                      E13+F12)
       DU2(KAREA2(1,IELH2))=DU2(KAREA2(1,IELH2))+
     *                      (R1+C12+D12+
     *                      E13+F13)
       DU2(KAREA2(1,IELH3))=DU2(KAREA2(1,IELH3))+
     *                      (R1+C13+D12+
     *                      E12+F13)
       DU2(KAREA2(1,IELH4))=DU2(KAREA2(1,IELH4))+
     *                      (R1+C13+D13+
     *                      E12+F12)
      ENDIF
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH1))=DU2(KAREA2(2,IELH1))+
     *                      (B12+R2+D13+
     *                      F12+G13)
       DU2(KAREA2(5,IELH2))=DU2(KAREA2(5,IELH2))+
     *                      (B12+R2+D12+
     *                      F13+G13)
       DU2(KAREA2(5,IELH6))=DU2(KAREA2(5,IELH6))+
     *                      (B13+R2+D12+
     *                      F13+G12)
       DU2(KAREA2(2,IELH5))=DU2(KAREA2(2,IELH5))+
     *                      (B13+R2+D13+
     *                      F12+G12)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH2))=DU2(KAREA2(2,IELH2))+
     *                      (B12+C12+R3+
     *                      E13+G13)
       DU2(KAREA2(5,IELH3))=DU2(KAREA2(5,IELH3))+
     *                      (B12+C13+R3+
     *                      E12+G13)
       DU2(KAREA2(5,IELH7))=DU2(KAREA2(5,IELH7))+
     *                      (B13+C13+R3+
     *                      E12+G12)
       DU2(KAREA2(2,IELH6))=DU2(KAREA2(2,IELH6))+
     *                      (B13+C12+R3+
     *                      E13+G12)
      ENDIF
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH3))=DU2(KAREA2(2,IELH3))+
     *                      (B12+R4+D12+
     *                      F13+G13)
       DU2(KAREA2(5,IELH4))=DU2(KAREA2(5,IELH4))+
     *                      (B12+R4+D13+
     *                      F12+G13)
       DU2(KAREA2(5,IELH8))=DU2(KAREA2(5,IELH8))+
     *                      (B13+R4+D13+
     *                      F12+G12)
       DU2(KAREA2(2,IELH7))=DU2(KAREA2(2,IELH7))+
     *                      (B13+R4+D12+
     *                      F13+G12)
      ENDIF
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH4))=DU2(KAREA2(2,IELH4))+
     *                      (B12+C13+R5+
     *                      E12+G13)
       DU2(KAREA2(5,IELH1))=DU2(KAREA2(5,IELH1))+
     *                      (B12+C12+R5+
     *                      E13+G13)
       DU2(KAREA2(5,IELH5))=DU2(KAREA2(5,IELH5))+
     *                      (B13+C12+R5+
     *                      E13+G12)
       DU2(KAREA2(2,IELH8))=DU2(KAREA2(2,IELH8))+
     *                      (B13+C13+R5+
     *                      E12+G12)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DU2(KAREA2(1,IELH5))=DU2(KAREA2(1,IELH5))+
     *                      (R6+C12+D13+
     *                      E13+F12)
       DU2(KAREA2(1,IELH6))=DU2(KAREA2(1,IELH6))+
     *                      (R6+C12+D12+
     *                      E13+F13)
       DU2(KAREA2(1,IELH7))=DU2(KAREA2(1,IELH7))+
     *                      (R6+C13+D12+
     *                      E12+F13)
       DU2(KAREA2(1,IELH8))=DU2(KAREA2(1,IELH8))+
     *                      (R6+C13+D13+
     *                      E12+F12)
      ENDIF
C
C
      DU2(KAREA2(3,IELH1))=R7+C15+D16+
     *                     E17+F16
      DU2(KAREA2(3,IELH2))=R7+C16+D15+
     *                     E16+F17
      DU2(KAREA2(3,IELH3))=R7+C17+D16+
     *                     E15+F16
      DU2(KAREA2(3,IELH4))=R7+C16+D17+
     *                     E16+F15
C
      DU2(KAREA2(6,IELH1))=R8+C15+D17+
     *                     E17+F15
      DU2(KAREA2(6,IELH2))=R8+C15+D15+
     *                     E17+F17
      DU2(KAREA2(6,IELH3))=R8+C17+D15+
     *                     E15+F17
      DU2(KAREA2(6,IELH4))=R8+C17+D17+
     *                     E15+F15
C
      DU2(KAREA2(3,IELH5))=R9+C15+D16+
     *                     E17+F16
      DU2(KAREA2(3,IELH6))=R9+C16+D15+
     *                     E16+F17
      DU2(KAREA2(3,IELH7))=R9+C17+D16+
     *                     E15+F16
      DU2(KAREA2(3,IELH8))=R9+C16+D17+
     *                     E16+F15
C
      B21=A1*DVH1
      B22=A2*DVH1
      B23=A3*DVH1
      B24=A4*DVH1
      B25=A5*DVH1
      B26=A6*DVH1
      B27=A7*DVH1
C
      C21=A1*DVH2
      C22=A2*DVH2
      C23=A3*DVH2
      C24=A4*DVH2
      C25=A5*DVH2
      C26=A6*DVH2
      C27=A7*DVH2
C
      D21=A1*DVH3
      D22=A2*DVH3
      D23=A3*DVH3
      D24=A4*DVH3
      D25=A5*DVH3
      D26=A6*DVH3
      D27=A7*DVH3
C
      E21=A1*DVH4
      E22=A2*DVH4
      E23=A3*DVH4
      E24=A4*DVH4
      E25=A5*DVH4
      E26=A6*DVH4
      E27=A7*DVH4
C
      F21=A1*DVH5
      F22=A2*DVH5
      F23=A3*DVH5
      F24=A4*DVH5
      F25=A5*DVH5
      F26=A6*DVH5
      F27=A7*DVH5
C
      G21=A1*DVH6
      G22=A2*DVH6
      G23=A3*DVH6
      G24=A4*DVH6
      G25=A5*DVH6
      G26=A6*DVH6
      G27=A7*DVH6
C
C
      S1=B21+G24
      S2=C21+E24
      S3=D21+F24
      S4=C24+E21
      S5=D24+F21
      S6=B24+G21
      S7=B25+G27
      S8=B26+G26
      S9=B27+G25
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DV2(KAREA2(1,IELH1))=DV2(KAREA2(1,IELH1))+
     *                      (S1+C22+D23+
     *                      E23+F22)
       DV2(KAREA2(1,IELH2))=DV2(KAREA2(1,IELH2))+
     *                      (S1+C22+D22+
     *                      E23+F23)
       DV2(KAREA2(1,IELH3))=DV2(KAREA2(1,IELH3))+
     *                      (S1+C23+D22+
     *                      E22+F23)
       DV2(KAREA2(1,IELH4))=DV2(KAREA2(1,IELH4))+
     *                      (S1+C23+D23+
     *                      E22+F22)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH1))=DV2(KAREA2(2,IELH1))+
     *                      (B22+S2+D23+
     *                      F22+G23)
       DV2(KAREA2(5,IELH2))=DV2(KAREA2(5,IELH2))+
     *                      (B22+S2+D22+
     *                      F23+G23)
       DV2(KAREA2(5,IELH6))=DV2(KAREA2(5,IELH6))+
     *                      (B23+S2+D22+
     *                      F23+G22)
       DV2(KAREA2(2,IELH5))=DV2(KAREA2(2,IELH5))+
     *                      (B23+S2+D23+
     *                      F22+G22)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH2))=DV2(KAREA2(2,IELH2))+
     *                      (B22+C22+S3+
     *                      E23+G23)
       DV2(KAREA2(5,IELH3))=DV2(KAREA2(5,IELH3))+
     *                      (B22+C23+S3+
     *                      E22+G23)
       DV2(KAREA2(5,IELH7))=DV2(KAREA2(5,IELH7))+
     *                      (B23+C23+S3+
     *                      E22+G22)
       DV2(KAREA2(2,IELH6))=DV2(KAREA2(2,IELH6))+
     *                      (B23+C22+S3+
     *                      E23+G22)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH3))=DV2(KAREA2(2,IELH3))+
     *                      (B22+S4+D22+
     *                      F23+G23)
       DV2(KAREA2(5,IELH4))=DV2(KAREA2(5,IELH4))+
     *                      (B22+S4+D23+
     *                      F22+G23)
       DV2(KAREA2(5,IELH8))=DV2(KAREA2(5,IELH8))+
     *                      (B23+S4+D23+
     *                      F22+G22)
       DV2(KAREA2(2,IELH7))=DV2(KAREA2(2,IELH7))+
     *                      (B23+S4+D22+
     *                      F23+G22)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH4))=DV2(KAREA2(2,IELH4))+
     *                      (B22+C23+S5+
     *                      E22+G23)
       DV2(KAREA2(5,IELH1))=DV2(KAREA2(5,IELH1))+
     *                      (B22+C22+S5+
     *                      E23+G23)
       DV2(KAREA2(5,IELH5))=DV2(KAREA2(5,IELH5))+
     *                      (B23+C22+S5+
     *                      E23+G22)
       DV2(KAREA2(2,IELH8))=DV2(KAREA2(2,IELH8))+
     *                      (B23+C23+S5+
     *                      E22+G22)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DV2(KAREA2(1,IELH5))=DV2(KAREA2(1,IELH5))+
     *                      (S6+C22+D23+
     *                      E23+F22)
       DV2(KAREA2(1,IELH6))=DV2(KAREA2(1,IELH6))+
     *                      (S6+C22+D22+
     *                      E23+F23)
       DV2(KAREA2(1,IELH7))=DV2(KAREA2(1,IELH7))+
     *                      (S6+C23+D22+
     *                      E22+F23)
       DV2(KAREA2(1,IELH8))=DV2(KAREA2(1,IELH8))+
     *                      (S6+C23+D23+
     *                      E22+F22)
      ENDIF
C
C
      DV2(KAREA2(3,IELH1))=S7+C25+D26+
     *                     E27+F26
      DV2(KAREA2(3,IELH2))=S7+C26+D25+
     *                     E26+F27
      DV2(KAREA2(3,IELH3))=S7+C27+D26+
     *                     E25+F26
      DV2(KAREA2(3,IELH4))=S7+C26+D27+
     *                     E26+F25
C
      DV2(KAREA2(6,IELH1))=S8+C25+D27+
     *                     E27+F25
      DV2(KAREA2(6,IELH2))=S8+C25+D25+
     *                     E27+F27
      DV2(KAREA2(6,IELH3))=S8+C27+D25+
     *                     E25+F27
      DV2(KAREA2(6,IELH4))=S8+C27+D27+
     *                     E25+F25
C
      DV2(KAREA2(3,IELH5))=S9+C25+D26+
     *                     E27+F26
      DV2(KAREA2(3,IELH6))=S9+C26+D25+
     *                     E26+F27
      DV2(KAREA2(3,IELH7))=S9+C27+D26+
     *                     E25+F26
      DV2(KAREA2(3,IELH8))=S9+C26+D27+
     *                     E26+F25
C
      B31=A1*DWH1
      B32=A2*DWH1
      B33=A3*DWH1
      B34=A4*DWH1
      B35=A5*DWH1
      B36=A6*DWH1
      B37=A7*DWH1
C
      C31=A1*DWH2
      C32=A2*DWH2
      C33=A3*DWH2
      C34=A4*DWH2
      C35=A5*DWH2
      C36=A6*DWH2
      C37=A7*DWH2
C
      D31=A1*DWH3
      D32=A2*DWH3
      D33=A3*DWH3
      D34=A4*DWH3
      D35=A5*DWH3
      D36=A6*DWH3
      D37=A7*DWH3
C
      E31=A1*DWH4
      E32=A2*DWH4
      E33=A3*DWH4
      E34=A4*DWH4
      E35=A5*DWH4
      E36=A6*DWH4
      E37=A7*DWH4
C
      F31=A1*DWH5
      F32=A2*DWH5
      F33=A3*DWH5
      F34=A4*DWH5
      F35=A5*DWH5
      F36=A6*DWH5
      F37=A7*DWH5
C
      G31=A1*DWH6
      G32=A2*DWH6
      G33=A3*DWH6
      G34=A4*DWH6
      G35=A5*DWH6
      G36=A6*DWH6
      G37=A7*DWH6
C
C
      T1=B31+G34
      T2=C31+E34
      T3=D31+F34
      T4=C34+E31
      T5=D34+F31
      T6=B34+G31
      T7=B35+G37
      T8=B36+G36
      T9=B37+G35
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DW2(KAREA2(1,IELH1))=DW2(KAREA2(1,IELH1))+
     *                      (T1+C32+D33+
     *                      E33+F32)
       DW2(KAREA2(1,IELH2))=DW2(KAREA2(1,IELH2))+
     *                      (T1+C32+D32+
     *                      E33+F33)
       DW2(KAREA2(1,IELH3))=DW2(KAREA2(1,IELH3))+
     *                      (T1+C33+D32+
     *                      E32+F33)
       DW2(KAREA2(1,IELH4))=DW2(KAREA2(1,IELH4))+
     *                      (T1+C33+D33+
     *                      E32+F32)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH1))=DW2(KAREA2(2,IELH1))+
     *                      (B32+T2+D33+
     *                      F32+G33)
       DW2(KAREA2(5,IELH2))=DW2(KAREA2(5,IELH2))+
     *                      (B32+T2+D32+
     *                      F33+G33)
       DW2(KAREA2(5,IELH6))=DW2(KAREA2(5,IELH6))+
     *                      (B33+T2+D32+
     *                      F33+G32)
       DW2(KAREA2(2,IELH5))=DW2(KAREA2(2,IELH5))+
     *                      (B33+T2+D33+
     *                      F32+G32)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH2))=DW2(KAREA2(2,IELH2))+
     *                      (B32+C32+T3+
     *                      E33+G33)
       DW2(KAREA2(5,IELH3))=DW2(KAREA2(5,IELH3))+
     *                      (B32+C33+T3+
     *                      E32+G33)
       DW2(KAREA2(5,IELH7))=DW2(KAREA2(5,IELH7))+
     *                      (B33+C33+T3+
     *                      E32+G32)
       DW2(KAREA2(2,IELH6))=DW2(KAREA2(2,IELH6))+
     *                      (B33+C32+T3+
     *                      E33+G32)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH3))=DW2(KAREA2(2,IELH3))+
     *                      (B32+T4+D32+
     *                      F33+G33)
       DW2(KAREA2(5,IELH4))=DW2(KAREA2(5,IELH4))+
     *                      (B32+T4+D33+
     *                      F32+G33)
       DW2(KAREA2(5,IELH8))=DW2(KAREA2(5,IELH8))+
     *                      (B33+T4+D33+
     *                      F32+G32)
       DW2(KAREA2(2,IELH7))=DW2(KAREA2(2,IELH7))+
     *                      (B33+T4+D32+
     *                      F33+G32)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH4))=DW2(KAREA2(2,IELH4))+
     *                      (B32+C33+T5+
     *                      E32+G33)
       DW2(KAREA2(5,IELH1))=DW2(KAREA2(5,IELH1))+
     *                      (B32+C32+T5+
     *                      E33+G33)
       DW2(KAREA2(5,IELH5))=DW2(KAREA2(5,IELH5))+
     *                      (B33+C32+T5+
     *                      E33+G32)
       DW2(KAREA2(2,IELH8))=DW2(KAREA2(2,IELH8))+
     *                      (B33+C33+T5+
     *                      E32+G32)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DW2(KAREA2(1,IELH5))=DW2(KAREA2(1,IELH5))+
     *                      (T6+C32+D33+
     *                      E33+F32)
       DW2(KAREA2(1,IELH6))=DW2(KAREA2(1,IELH6))+
     *                      (T6+C32+D32+
     *                      E33+F33)
       DW2(KAREA2(1,IELH7))=DW2(KAREA2(1,IELH7))+
     *                      (T6+C33+D32+
     *                      E32+F33)
       DW2(KAREA2(1,IELH8))=DW2(KAREA2(1,IELH8))+
     *                      (T6+C33+D33+
     *                      E32+F32)
      ENDIF
C
C
      DW2(KAREA2(3,IELH1))=T7+C35+D36+
     *                     E37+F36
      DW2(KAREA2(3,IELH2))=T7+C36+D35+
     *                     E36+F37
      DW2(KAREA2(3,IELH3))=T7+C37+D36+
     *                     E35+F36
      DW2(KAREA2(3,IELH4))=T7+C36+D37+
     *                     E36+F35
C
      DW2(KAREA2(6,IELH1))=T8+C35+D37+
     *                     E37+F35
      DW2(KAREA2(6,IELH2))=T8+C35+D35+
     *                     E37+F37
      DW2(KAREA2(6,IELH3))=T8+C37+D35+
     *                     E35+F37
      DW2(KAREA2(6,IELH4))=T8+C37+D37+
     *                     E35+F35
C
      DW2(KAREA2(3,IELH5))=T9+C35+D36+
     *                     E37+F36
      DW2(KAREA2(3,IELH6))=T9+C36+D35+
     *                     E36+F37
      DW2(KAREA2(3,IELH7))=T9+C37+D36+
     *                     E35+F36
      DW2(KAREA2(3,IELH8))=T9+C36+D37+
     *                     E36+F35
C
C
C
10    CONTINUE
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRM  (DU1,DV1,DW1,DP1,DU2,DV2,DW2,DP2,
     *                     KAREA1,KADJ1,NAT1,NEL1,
     *                     KAREA2,KADJ2,NAT2,NEL2,IINT)
************************************************************************
*    Purpose:    Restricts the fine grid defect vector (DU2,DV2,DW2,DP2)
*                to the coarse grid defect vector (DU1,DV1,DW1,DP1)
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2,DW2,DP2           - fine grid defect vector
*      KAREA1,KADJ1,NAT1,NEL1    - data of the coarse grid
*      KAREA2,KADJ2,NAT2,NEL2    - data of the FIN2 grid
*
*    Output:
*      DU1,DV1,DW1,DP1           - coarse grid defect vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NNAE=6)
      DIMENSION DU1(1),DV1(1),DW1(1),DP1(1),DU2(1),DV2(1),DW2(1),DP2(1),
     *          KAREA2(NNAE,*),KAREA1(NNAE,*),
     *          KADJ1(NNAE,*),KADJ2(NNAE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      IF (IINT.GT.0) THEN
       A1= 11D0/24D0
       A2=  7D0/48D0
       A3=- 5D0/48D0
       A4=- 1D0/24D0
       A5= 11D0/24D0
       A6=  1D0/12D0
       A7=- 1D0/24D0
      ELSE
       A1= 0.500D0
       A2= 0.125D0
       A3=-0.125D0
       A4= 0.000D0
       A5= 0.500D0
       A6= 0.000D0
       A7= 0.000D0
      ENDIF
C
C
C
      CALL LCL1(DU1,NAT1)
      CALL LCL1(DV1,NAT1)
      CALL LCL1(DW1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
C *** Restriction of pressure
C
      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)+
     *           DP2(IELH5)+DP2(IELH6)+DP2(IELH7)+DP2(IELH8)
C
C
      DUH1= DU2(KAREA2(1,IELH1))
      DUH2= DU2(KAREA2(1,IELH2))
      DUH3= DU2(KAREA2(1,IELH3))
      DUH4= DU2(KAREA2(1,IELH4))
      DUH5= DU2(KAREA2(2,IELH1))
      DUH6= DU2(KAREA2(5,IELH2))
      DUH7= DU2(KAREA2(5,IELH6))
      DUH8= DU2(KAREA2(2,IELH5))
      DUH9= DU2(KAREA2(2,IELH2))
      DUH10=DU2(KAREA2(5,IELH3))
      DUH11=DU2(KAREA2(5,IELH7))
      DUH12=DU2(KAREA2(2,IELH6))
      DUH13=DU2(KAREA2(2,IELH3))
      DUH14=DU2(KAREA2(5,IELH4))
      DUH15=DU2(KAREA2(5,IELH8))
      DUH16=DU2(KAREA2(2,IELH7))
      DUH17=DU2(KAREA2(2,IELH4))
      DUH18=DU2(KAREA2(5,IELH1))
      DUH19=DU2(KAREA2(5,IELH5))
      DUH20=DU2(KAREA2(2,IELH8))
      DUH21=DU2(KAREA2(1,IELH5))
      DUH22=DU2(KAREA2(1,IELH6))
      DUH23=DU2(KAREA2(1,IELH7))
      DUH24=DU2(KAREA2(1,IELH8))
      DUH25=DU2(KAREA2(3,IELH1))
      DUH26=DU2(KAREA2(3,IELH2))
      DUH27=DU2(KAREA2(3,IELH3))
      DUH28=DU2(KAREA2(3,IELH4))
      DUH29=DU2(KAREA2(6,IELH1))
      DUH30=DU2(KAREA2(6,IELH2))
      DUH31=DU2(KAREA2(6,IELH3))
      DUH32=DU2(KAREA2(6,IELH4))
      DUH33=DU2(KAREA2(3,IELH5))
      DUH34=DU2(KAREA2(3,IELH6))
      DUH35=DU2(KAREA2(3,IELH7))
      DUH36=DU2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DU1(IA1)= DU1(IA1)
     *          +A1*(DUH1+DUH2+DUH3+DUH4)
     *          +A2*(DUH5+DUH6+DUH9+DUH10+DUH13+DUH14+DUH17+DUH18)
     *          +A3*(DUH7+DUH8+DUH11+DUH12+DUH15+DUH16+DUH19+DUH20)
     *          +A4*(DUH21+DUH22+DUH23+DUH24)
     *          +A5*(DUH25+DUH26+DUH27+DUH28)
     *          +A6*(DUH29+DUH30+DUH31+DUH32)
     *          +A7*(DUH33+DUH34+DUH35+DUH36)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DU1(IA2)= DU1(IA2)
     *          +A1*(DUH5+DUH6+DUH7+DUH8)
     *          +A2*(DUH1+DUH2+DUH9+DUH12+DUH21+DUH22+DUH18+DUH19)
     *          +A3*(DUH3+DUH4+DUH10+DUH11+DUH23+DUH24+DUH17+DUH20)
     *          +A4*(DUH13+DUH14+DUH15+DUH16)
     *          +A5*(DUH25+DUH29+DUH30+DUH33)
     *          +A6*(DUH26+DUH28+DUH34+DUH36)
     *          +A7*(DUH27+DUH31+DUH32+DUH35)     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DU1(IA3)= DU1(IA3)
     *          +A1*(DUH9+DUH10+DUH11+DUH12)
     *          +A2*(DUH2+DUH3+DUH6+DUH7+DUH22+DUH23+DUH13+DUH16)
     *          +A3*(DUH1+DUH4+DUH5+DUH8+DUH21+DUH24+DUH14+DUH15)
     *          +A4*(DUH17+DUH18+DUH19+DUH20)
     *          +A5*(DUH26+DUH30+DUH31+DUH34)
     *          +A6*(DUH25+DUH27+DUH33+DUH35)
     *          +A7*(DUH28+DUH29+DUH32+DUH36)     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DU1(IA4)= DU1(IA4)
     *          +A1*(DUH13+DUH14+DUH15+DUH16)
     *          +A2*(DUH3+DUH4+DUH10+DUH11+DUH23+DUH24+DUH17+DUH20)
     *          +A3*(DUH1+DUH2+DUH9+DUH12+DUH21+DUH22+DUH18+DUH19)
     *          +A4*(DUH5+DUH6+DUH7+DUH8)
     *          +A5*(DUH27+DUH31+DUH32+DUH35)  
     *          +A6*(DUH26+DUH28+DUH34+DUH36)
     *          +A7*(DUH25+DUH29+DUH30+DUH33)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DU1(IA5)= DU1(IA5)
     *          +A1*(DUH17+DUH18+DUH19+DUH20)
     *          +A2*(DUH1+DUH4+DUH5+DUH8+DUH21+DUH24+DUH14+DUH15)
     *          +A3*(DUH2+DUH3+DUH6+DUH7+DUH22+DUH23+DUH13+DUH16)
     *          +A4*(DUH9+DUH10+DUH11+DUH12)
     *          +A5*(DUH28+DUH29+DUH32+DUH36)     
     *          +A6*(DUH25+DUH27+DUH33+DUH35)
     *          +A7*(DUH26+DUH30+DUH31+DUH34)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DU1(IA6)= DU1(IA6)
     *          +A1*(DUH21+DUH22+DUH23+DUH24)
     *          +A2*(DUH7+DUH8+DUH11+DUH12+DUH15+DUH16+DUH19+DUH20)
     *          +A3*(DUH5+DUH6+DUH9+DUH10+DUH13+DUH14+DUH17+DUH18)
     *          +A4*(DUH1+DUH2+DUH3+DUH4)
     *          +A5*(DUH33+DUH34+DUH35+DUH36)     
     *          +A6*(DUH29+DUH30+DUH31+DUH32)
     *          +A7*(DUH25+DUH26+DUH27+DUH28)
      ENDIF
C
C
C
      DVH1= DV2(KAREA2(1,IELH1))
      DVH2= DV2(KAREA2(1,IELH2))
      DVH3= DV2(KAREA2(1,IELH3))
      DVH4= DV2(KAREA2(1,IELH4))
      DVH5= DV2(KAREA2(2,IELH1))
      DVH6= DV2(KAREA2(5,IELH2))
      DVH7= DV2(KAREA2(5,IELH6))
      DVH8= DV2(KAREA2(2,IELH5))
      DVH9= DV2(KAREA2(2,IELH2))
      DVH10=DV2(KAREA2(5,IELH3))
      DVH11=DV2(KAREA2(5,IELH7))
      DVH12=DV2(KAREA2(2,IELH6))
      DVH13=DV2(KAREA2(2,IELH3))
      DVH14=DV2(KAREA2(5,IELH4))
      DVH15=DV2(KAREA2(5,IELH8))
      DVH16=DV2(KAREA2(2,IELH7))
      DVH17=DV2(KAREA2(2,IELH4))
      DVH18=DV2(KAREA2(5,IELH1))
      DVH19=DV2(KAREA2(5,IELH5))
      DVH20=DV2(KAREA2(2,IELH8))
      DVH21=DV2(KAREA2(1,IELH5))
      DVH22=DV2(KAREA2(1,IELH6))
      DVH23=DV2(KAREA2(1,IELH7))
      DVH24=DV2(KAREA2(1,IELH8))
      DVH25=DV2(KAREA2(3,IELH1))
      DVH26=DV2(KAREA2(3,IELH2))
      DVH27=DV2(KAREA2(3,IELH3))
      DVH28=DV2(KAREA2(3,IELH4))
      DVH29=DV2(KAREA2(6,IELH1))
      DVH30=DV2(KAREA2(6,IELH2))
      DVH31=DV2(KAREA2(6,IELH3))
      DVH32=DV2(KAREA2(6,IELH4))
      DVH33=DV2(KAREA2(3,IELH5))
      DVH34=DV2(KAREA2(3,IELH6))
      DVH35=DV2(KAREA2(3,IELH7))
      DVH36=DV2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DV1(IA1)= DV1(IA1)
     *          +A1*(DVH1+DVH2+DVH3+DVH4)
     *          +A2*(DVH5+DVH6+DVH9+DVH10+DVH13+DVH14+DVH17+DVH18)
     *          +A3*(DVH7+DVH8+DVH11+DVH12+DVH15+DVH16+DVH19+DVH20)
     *          +A4*(DVH21+DVH22+DVH23+DVH24)
     *          +A5*(DVH25+DVH26+DVH27+DVH28)
     *          +A6*(DVH29+DVH30+DVH31+DVH32)
     *          +A7*(DVH33+DVH34+DVH35+DVH36)     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DV1(IA2)= DV1(IA2)
     *          +A1*(DVH5+DVH6+DVH7+DVH8)
     *          +A2*(DVH1+DVH2+DVH9+DVH12+DVH21+DVH22+DVH18+DVH19)
     *          +A3*(DVH3+DVH4+DVH10+DVH11+DVH23+DVH24+DVH17+DVH20)
     *          +A4*(DVH13+DVH14+DVH15+DVH16)
     *          +A5*(DVH25+DVH29+DVH30+DVH33)
     *          +A6*(DVH26+DVH28+DVH34+DVH36)
     *          +A7*(DVH27+DVH31+DVH32+DVH35)     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DV1(IA3)= DV1(IA3)
     *          +A1*(DVH9+DVH10+DVH11+DVH12)
     *          +A2*(DVH2+DVH3+DVH6+DVH7+DVH22+DVH23+DVH13+DVH16)
     *          +A3*(DVH1+DVH4+DVH5+DVH8+DVH21+DVH24+DVH14+DVH15)
     *          +A4*(DVH17+DVH18+DVH19+DVH20)
     *          +A5*(DVH26+DVH30+DVH31+DVH34)
     *          +A6*(DVH25+DVH27+DVH33+DVH35)
     *          +A7*(DVH28+DVH29+DVH32+DVH36)     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DV1(IA4)= DV1(IA4)
     *          +A1*(DVH13+DVH14+DVH15+DVH16)
     *          +A2*(DVH3+DVH4+DVH10+DVH11+DVH23+DVH24+DVH17+DVH20)
     *          +A3*(DVH1+DVH2+DVH9+DVH12+DVH21+DVH22+DVH18+DVH19)
     *          +A4*(DVH5+DVH6+DVH7+DVH8)
     *          +A5*(DVH27+DVH31+DVH32+DVH35)  
     *          +A6*(DVH26+DVH28+DVH34+DVH36)
     *          +A7*(DVH25+DVH29+DVH30+DVH33)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DV1(IA5)= DV1(IA5)
     *          +A1*(DVH17+DVH18+DVH19+DVH20)
     *          +A2*(DVH1+DVH4+DVH5+DVH8+DVH21+DVH24+DVH14+DVH15)
     *          +A3*(DVH2+DVH3+DVH6+DVH7+DVH22+DVH23+DVH13+DVH16)
     *          +A4*(DVH9+DVH10+DVH11+DVH12)
     *          +A5*(DVH28+DVH29+DVH32+DVH36)     
     *          +A6*(DVH25+DVH27+DVH33+DVH35)
     *          +A7*(DVH26+DVH30+DVH31+DVH34)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DV1(IA6)= DV1(IA6)
     *          +A1*(DVH21+DVH22+DVH23+DVH24)
     *          +A2*(DVH7+DVH8+DVH11+DVH12+DVH15+DVH16+DVH19+DVH20)
     *          +A3*(DVH5+DVH6+DVH9+DVH10+DVH13+DVH14+DVH17+DVH18)
     *          +A4*(DVH1+DVH2+DVH3+DVH4)
     *          +A5*(DVH33+DVH34+DVH35+DVH36)     
     *          +A6*(DVH29+DVH30+DVH31+DVH32)
     *          +A7*(DVH25+DVH26+DVH27+DVH28)
      ENDIF
C
C
C
      DWH1= DW2(KAREA2(1,IELH1))
      DWH2= DW2(KAREA2(1,IELH2))
      DWH3= DW2(KAREA2(1,IELH3))
      DWH4= DW2(KAREA2(1,IELH4))
      DWH5= DW2(KAREA2(2,IELH1))
      DWH6= DW2(KAREA2(5,IELH2))
      DWH7= DW2(KAREA2(5,IELH6))
      DWH8= DW2(KAREA2(2,IELH5))
      DWH9= DW2(KAREA2(2,IELH2))
      DWH10=DW2(KAREA2(5,IELH3))
      DWH11=DW2(KAREA2(5,IELH7))
      DWH12=DW2(KAREA2(2,IELH6))
      DWH13=DW2(KAREA2(2,IELH3))
      DWH14=DW2(KAREA2(5,IELH4))
      DWH15=DW2(KAREA2(5,IELH8))
      DWH16=DW2(KAREA2(2,IELH7))
      DWH17=DW2(KAREA2(2,IELH4))
      DWH18=DW2(KAREA2(5,IELH1))
      DWH19=DW2(KAREA2(5,IELH5))
      DWH20=DW2(KAREA2(2,IELH8))
      DWH21=DW2(KAREA2(1,IELH5))
      DWH22=DW2(KAREA2(1,IELH6))
      DWH23=DW2(KAREA2(1,IELH7))
      DWH24=DW2(KAREA2(1,IELH8))
      DWH25=DW2(KAREA2(3,IELH1))
      DWH26=DW2(KAREA2(3,IELH2))
      DWH27=DW2(KAREA2(3,IELH3))
      DWH28=DW2(KAREA2(3,IELH4))
      DWH29=DW2(KAREA2(6,IELH1))
      DWH30=DW2(KAREA2(6,IELH2))
      DWH31=DW2(KAREA2(6,IELH3))
      DWH32=DW2(KAREA2(6,IELH4))
      DWH33=DW2(KAREA2(3,IELH5))
      DWH34=DW2(KAREA2(3,IELH6))
      DWH35=DW2(KAREA2(3,IELH7))
      DWH36=DW2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DW1(IA1)= DW1(IA1)
     *          +A1*(DWH1+DWH2+DWH3+DWH4)
     *          +A2*(DWH5+DWH6+DWH9+DWH10+DWH13+DWH14+DWH17+DWH18)
     *          +A3*(DWH7+DWH8+DWH11+DWH12+DWH15+DWH16+DWH19+DWH20)
     *          +A4*(DWH21+DWH22+DWH23+DWH24)
     *          +A5*(DWH25+DWH26+DWH27+DWH28)
     *          +A6*(DWH29+DWH30+DWH31+DWH32)
     *          +A7*(DWH33+DWH34+DWH35+DWH36)     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DW1(IA2)= DW1(IA2)
     *          +A1*(DWH5+DWH6+DWH7+DWH8)
     *          +A2*(DWH1+DWH2+DWH9+DWH12+DWH21+DWH22+DWH18+DWH19)
     *          +A3*(DWH3+DWH4+DWH10+DWH11+DWH23+DWH24+DWH17+DWH20)
     *          +A4*(DWH13+DWH14+DWH15+DWH16)
     *          +A5*(DWH25+DWH29+DWH30+DWH33)
     *          +A6*(DWH26+DWH28+DWH34+DWH36)
     *          +A7*(DWH27+DWH31+DWH32+DWH35)     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DW1(IA3)= DW1(IA3)
     *          +A1*(DWH9+DWH10+DWH11+DWH12)
     *          +A2*(DWH2+DWH3+DWH6+DWH7+DWH22+DWH23+DWH13+DWH16)
     *          +A3*(DWH1+DWH4+DWH5+DWH8+DWH21+DWH24+DWH14+DWH15)
     *          +A4*(DWH17+DWH18+DWH19+DWH20)
     *          +A5*(DWH26+DWH30+DWH31+DWH34)
     *          +A6*(DWH25+DWH27+DWH33+DWH35)
     *          +A7*(DWH28+DWH29+DWH32+DWH36)     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DW1(IA4)= DW1(IA4)
     *          +A1*(DWH13+DWH14+DWH15+DWH16)
     *          +A2*(DWH3+DWH4+DWH10+DWH11+DWH23+DWH24+DWH17+DWH20)
     *          +A3*(DWH1+DWH2+DWH9+DWH12+DWH21+DWH22+DWH18+DWH19)
     *          +A4*(DWH5+DWH6+DWH7+DWH8)
     *          +A5*(DWH27+DWH31+DWH32+DWH35)  
     *          +A6*(DWH26+DWH28+DWH34+DWH36)
     *          +A7*(DWH25+DWH29+DWH30+DWH33)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DW1(IA5)= DW1(IA5)
     *          +A1*(DWH17+DWH18+DWH19+DWH20)
     *          +A2*(DWH1+DWH4+DWH5+DWH8+DWH21+DWH24+DWH14+DWH15)
     *          +A3*(DWH2+DWH3+DWH6+DWH7+DWH22+DWH23+DWH13+DWH16)
     *          +A4*(DWH9+DWH10+DWH11+DWH12)
     *          +A5*(DWH28+DWH29+DWH32+DWH36)     
     *          +A6*(DWH25+DWH27+DWH33+DWH35)
     *          +A7*(DWH26+DWH30+DWH31+DWH34)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DW1(IA6)= DW1(IA6)
     *          +A1*(DWH21+DWH22+DWH23+DWH24)
     *          +A2*(DWH7+DWH8+DWH11+DWH12+DWH15+DWH16+DWH19+DWH20)
     *          +A3*(DWH5+DWH6+DWH9+DWH10+DWH13+DWH14+DWH17+DWH18)
     *          +A4*(DWH1+DWH2+DWH3+DWH4)
     *          +A5*(DWH33+DWH34+DWH35+DWH36)     
     *          +A6*(DWH29+DWH30+DWH31+DWH32)
     *          +A7*(DWH25+DWH26+DWH27+DWH28)
      ENDIF
C
10    CONTINUE
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRU  (DU1,DV1,DW1,DU2,DV2,DW2,
     *                     KVERT1,KAREA1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KAREA2,KADJ2,NEQ2,NEL2,NVT2,
     *                     AVOL)
************************************************************************
*    Purpose:  Restricts the  fine grid solution vector (DU2,DV2,DW2) to
*              the coarse grid solution vector (DU1,DV1,DW1)
*
*    Remark:   i.) RESTRU works for all cases of boundary conditions
*              ii.) weighted interpolation included
*                
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2,DW2               - fine grid solution vector
*      KVERT1,KAREA1,..,NVT1  - data of the coarse grid
*      KVERT2,KAREA2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1,DV1               - coarse grid solution vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL AVOL
C
      PARAMETER (NNVE=8,NNAE=6)
      PARAMETER (A1=0.125D0,A2=0.25D0,A3=0.08333333D0)
      PARAMETER (R1=0.25D0,R2=0.5D0,R3=0.16666666D0)
C
      DIMENSION DU1(*),DV1(*),DW1(*),  DU2(*),DV2(*),DW2(*),
     *          KVERT1(NNVE,*),KAREA1(NNAE,*),KADJ1(NNAE,*),
     *          KVERT2(NNVE,*),KAREA2(NNAE,*),KADJ2(NNAE,*),
     *          AVOL(*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)
C
C
      I1=KAREA2(1,IELH1)
      I2=KAREA2(2,IELH1)
      I3=KAREA2(3,IELH1)
      I4=KAREA2(4,IELH1)
      I5=KAREA2(5,IELH1)
      I6=KAREA2(6,IELH1)
      I7=KAREA2(1,IELH2)
      I8=KAREA2(2,IELH2)
      I9=KAREA2(3,IELH2)
      I10=KAREA2(5,IELH2)
      I11=KAREA2(6,IELH2)
      I12=KAREA2(1,IELH3)
      I13=KAREA2(2,IELH3)
      I14=KAREA2(3,IELH3)
      I15=KAREA2(5,IELH3)
      I16=KAREA2(6,IELH3)
      I17=KAREA2(1,IELH4)
      I18=KAREA2(2,IELH4)
      I19=KAREA2(5,IELH4)
      I20=KAREA2(6,IELH4)
      I21=KAREA2(1,IELH5)
      I22=KAREA2(2,IELH5)
      I23=KAREA2(3,IELH5)
      I24=KAREA2(4,IELH5)
      I25=KAREA2(5,IELH5)
      I26=KAREA2(1,IELH6)
      I27=KAREA2(2,IELH6)
      I28=KAREA2(3,IELH6)
      I29=KAREA2(5,IELH6)
      I30=KAREA2(1,IELH7)
      I31=KAREA2(2,IELH7)
      I32=KAREA2(3,IELH7)
      I33=KAREA2(5,IELH7)
      I34=KAREA2(1,IELH8)
      I35=KAREA2(2,IELH8)
      I36=KAREA2(5,IELH8)
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
      DUH13=DU2(I13)
      DUH14=DU2(I14)
      DUH15=DU2(I15)
      DUH16=DU2(I16)
      DUH17=DU2(I17)
      DUH18=DU2(I18)
      DUH19=DU2(I19)
      DUH20=DU2(I20)
      DUH21=DU2(I21)
      DUH22=DU2(I22)
      DUH23=DU2(I23)
      DUH24=DU2(I24)
      DUH25=DU2(I25)
      DUH26=DU2(I26)
      DUH27=DU2(I27)
      DUH28=DU2(I28)
      DUH29=DU2(I29)
      DUH30=DU2(I30)
      DUH31=DU2(I31)
      DUH32=DU2(I32)
      DUH33=DU2(I33)
      DUH34=DU2(I34)
      DUH35=DU2(I35)
      DUH36=DU2(I36)
C
C
      DVH1= DV2(I1)
      DVH2= DV2(I2)
      DVH3= DV2(I3)
      DVH4= DV2(I4)
      DVH5= DV2(I5)
      DVH6= DV2(I6)
      DVH7= DV2(I7)
      DVH8= DV2(I8)
      DVH9= DV2(I9)
      DVH10=DV2(I10)
      DVH11=DV2(I11)
      DVH12=DV2(I12)
      DVH13=DV2(I13)
      DVH14=DV2(I14)
      DVH15=DV2(I15)
      DVH16=DV2(I16)
      DVH17=DV2(I17)
      DVH18=DV2(I18)
      DVH19=DV2(I19)
      DVH20=DV2(I20)
      DVH21=DV2(I21)
      DVH22=DV2(I22)
      DVH23=DV2(I23)
      DVH24=DV2(I24)
      DVH25=DV2(I25)
      DVH26=DV2(I26)
      DVH27=DV2(I27)
      DVH28=DV2(I28)
      DVH29=DV2(I29)
      DVH30=DV2(I30)
      DVH31=DV2(I31)
      DVH32=DV2(I32)
      DVH33=DV2(I33)
      DVH34=DV2(I34)
      DVH35=DV2(I35)
      DVH36=DV2(I36)
C
      DWH1= DW2(I1)
      DWH2= DW2(I2)
      DWH3= DW2(I3)
      DWH4= DW2(I4)
      DWH5= DW2(I5)
      DWH6= DW2(I6)
      DWH7= DW2(I7)
      DWH8= DW2(I8)
      DWH9= DW2(I9)
      DWH10=DW2(I10)
      DWH11=DW2(I11)
      DWH12=DW2(I12)
      DWH13=DW2(I13)
      DWH14=DW2(I14)
      DWH15=DW2(I15)
      DWH16=DW2(I16)
      DWH17=DW2(I17)
      DWH18=DW2(I18)
      DWH19=DW2(I19)
      DWH20=DW2(I20)
      DWH21=DW2(I21)
      DWH22=DW2(I22)
      DWH23=DW2(I23)
      DWH24=DW2(I24)
      DWH25=DW2(I25)
      DWH26=DW2(I26)
      DWH27=DW2(I27)
      DWH28=DW2(I28)
      DWH29=DW2(I29)
      DWH30=DW2(I30)
      DWH31=DW2(I31)
      DWH32=DW2(I32)
      DWH33=DW2(I33)
      DWH34=DW2(I34)
      DWH35=DW2(I35)
      DWH36=DW2(I36)
      AVOL11=AVOL(IEL1)
C
C *** The area IA1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DU1(IA1)=(A1*(DUH1+DUH7+DUH12+DUH17)+
     *    A2*(DUH3+DUH14+DUH4+DUH9)-
     *    A3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+
     *        DUH5+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11
        DV1(IA1)=(A1*(DVH1+DVH7+DVH12+DVH17)+
     *    A2*(DVH3+DVH14+DVH4+DVH9)-
     *    A3*(DVH2+DVH10+DVH8+DVH15+DVH13+DVH19+DVH18+
     *        DVH5+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11
        DW1(IA1)=(A1*(DWH1+DWH7+DWH12+DWH17)+
     *    A2*(DWH3+DWH14+DWH4+DWH9)-
     *    A3*(DWH2+DWH10+DWH8+DWH15+DWH13+DWH19+DWH18+
     *        DWH5+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11
      ELSE
        DU1(IA1)=DU1(IA1)+
     *    (A1*(DUH1+DUH7+DUH12+DUH17)+A2*(DUH3+DUH14+DUH4+DUH9)-
     *    A3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+
     *        DUH5+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11
        DV1(IA1)=DV1(IA1)+
     *    (A1*(DVH1+DVH7+DVH12+DVH17)+A2*(DVH3+DVH14+DVH4+DVH9)-
     *    A3*(DVH2+DVH10+DVH8+DVH15+DVH13+DVH19+DVH18+
     *        DVH5+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11
        DW1(IA1)=DW1(IA1)+
     *    (A1*(DWH1+DWH7+DWH12+DWH17)+A2*(DWH3+DWH14+DWH4+DWH9)-
     *    A3*(DWH2+DWH10+DWH8+DWH15+DWH13+DWH19+DWH18+
     *        DWH5+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11
      ENDIF
      ELSE
C     case of a boundary edge
        DU1(IA1)=R1*(DUH1+DUH7+DUH12+DUH17)+
     *     R2*(DUH3+DUH14+DUH4+DUH9)-
     *    R3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+
     *        DUH5+DUH6+DUH11+
     *        DUH16+DUH20)
        DV1(IA1)=R1*(DVH1+DVH7+DVH12+DVH17)+
     *     R2*(DVH3+DVH14+DVH4+DVH9)-
     *    R3*(DVH2+DVH10+DVH8+DVH15+DVH13+DVH19+DVH18+
     *        DVH5+DVH6+DVH11+
     *        DVH16+DVH20)
        DW1(IA1)=R1*(DWH1+DWH7+DWH12+DWH17)+
     *    R2*(DWH3+DWH14+DWH4+DWH9)-
     *    R3*(DWH2+DWH10+DWH8+DWH15+DWH13+DWH19+DWH18+
     *        DWH5+DWH6+DWH11+
     *        DWH16+DWH20)
      ENDIF
C
C *** The area IA2
C
      IF (KADJ1(2,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DU1(IA2)=(A1*(DUH2+DUH10+DUH29+DUH22)+
     *     A2*(DUH3+DUH23+DUH6+DUH11)-
     *     A3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+
     *         DUH5+DUH4+DUH9+
     *         DUH24+DUH28))*AVOL11  
        DV1(IA2)=(A1*(DVH2+DVH10+DVH29+DVH22)+
     *     A2*(DVH3+DVH23+DVH6+DVH11)-
     *     A3*(DVH1+DVH7+DVH8+DVH27+DVH26+DVH21+DVH25+
     *         DVH5+DVH4+DVH9+
     *         DVH24+DVH28))*AVOL11  
        DW1(IA2)=(A1*(DWH2+DWH10+DWH29+DWH22)+
     *     A2*(DWH3+DWH23+DWH6+DWH11)-
     *     A3*(DWH1+DWH7+DWH8+DWH27+DWH26+DWH21+DWH25+
     *         DWH5+DWH4+DWH9+
     *         DWH24+DWH28))*AVOL11  
      ELSE
        DU1(IA2)=DU1(IA2)+
     *     (A1*(DUH2+DUH10+DUH29+DUH22)+A2*(DUH3+DUH23+DUH6+DUH11)-
     *     A3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+
     *         DUH5+DUH4+DUH9+
     *         DUH24+DUH28))*AVOL11  
        DV1(IA2)=DV1(IA2)+
     *     (A1*(DVH2+DVH10+DVH29+DVH22)+A2*(DVH3+DVH23+DVH6+DVH11)-
     *     A3*(DVH1+DVH7+DVH8+DVH27+DVH26+DVH21+DVH25+
     *         DVH5+DVH4+DVH9+
     *         DVH24+DVH28))*AVOL11  
        DW1(IA2)=DW1(IA2)+
     *     (A1*(DWH2+DWH10+DWH29+DWH22)+A2*(DWH3+DWH23+DWH6+DWH11)-
     *     A3*(DWH1+DWH7+DWH8+DWH27+DWH26+DWH21+DWH25+
     *         DWH5+DWH4+DWH9+
     *         DWH24+DWH28))*AVOL11
       ENDIF
       ELSE  
C     case of a boundary area
        DU1(IA2)=R1*(DUH2+DUH10+DUH29+DUH22)+
     *     R2*(DUH3+DUH23+DUH6+DUH11)-
     *     R3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+
     *         DUH5+DUH4+DUH9+
     *         DUH24+DUH28)  
        DV1(IA2)=R1*(DVH2+DVH10+DVH29+DVH22)+
     *     R2*(DVH3+DVH23+DVH6+DVH11)-
     *     R3*(DVH1+DVH7+DVH8+DVH27+DVH26+DVH21+DVH25+
     *         DVH5+DVH4+DVH9+
     *         DVH24+DVH28)  
        DW1(IA2)=R1*(DWH2+DWH10+DWH29+DWH22)+
     *     R2*(DWH3+DWH23+DWH6+DWH11)-
     *     R3*(DWH1+DWH7+DWH8+DWH27+DWH26+DWH21+DWH25+
     *         DWH5+DWH4+DWH9+
     *         DWH24+DWH28)  
      ENDIF
C
C *** The edge IA3
C
      IF (KADJ1(3,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DU1(IA3)=(A1*(DUH8+DUH15+DUH33+DUH27)+
     *     A2*(DUH11+DUH16+DUH9+DUH28)-
     *     A3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+
     *         DUH31+DUH30+DUH26+
     *          DUH29+DUH10))*AVOL11  
        DV1(IA3)=(A1*(DVH8+DVH15+DVH33+DVH27)+
     *     A2*(DVH11+DVH16+DVH9+DVH28)-
     *     A3*(DVH3+DVH14+DVH32+DVH23+DVH7+DVH12+DVH13+
     *         DVH31+DVH30+DVH26+
     *          DVH29+DVH10))*AVOL11  
        DW1(IA3)=(A1*(DWH8+DWH15+DWH33+DWH27)+
     *     A2*(DWH11+DWH16+DWH9+DWH28)-
     *     A3*(DWH3+DWH14+DWH32+DWH23+DWH7+DWH12+DWH13+
     *         DWH31+DWH30+DWH26+
     *          DWH29+DWH10))*AVOL11  
      ELSE
        DU1(IA3)= DU1(IA3)+
     *     (A1*(DUH8+DUH15+DUH33+DUH27)+A2*(DUH11+DUH16+DUH9+DUH28)-
     *     A3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+
     *         DUH31+DUH30+DUH26+
     *          DUH29+DUH10))*AVOL11  
        DV1(IA3)=DV1(IA3)+
     *     (A1*(DVH8+DVH15+DVH33+DVH27)+A2*(DVH11+DVH16+DVH9+DVH28)-
     *     A3*(DVH3+DVH14+DVH32+DVH23+DVH7+DVH12+DVH13+
     *         DVH31+DVH30+DVH26+
     *          DVH29+DVH10))*AVOL11  
        DW1(IA3)=DW1(IA3)+
     *     (A1*(DWH8+DWH15+DWH33+DWH27)+A2*(DWH11+DWH16+DWH9+DWH28)-
     *     A3*(DWH3+DWH14+DWH32+DWH23+DWH7+DWH12+DWH13+
     *         DWH31+DWH30+DWH26+
     *          DWH29+DWH10))*AVOL11
        ENDIF
        ELSE  
C     case of a boundary area 
        DU1(IA3)=R1*(DUH8+DUH15+DUH33+DUH27)+
     *     R2*(DUH11+DUH16+DUH9+DUH28)-
     *     R3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+
     *         DUH31+DUH30+DUH26+
     *          DUH29+DUH10)  
        DV1(IA3)=R1*(DVH8+DVH15+DVH33+DVH27)+
     *     R2*(DVH11+DVH16+DVH9+DVH28)-
     *     R3*(DVH3+DVH14+DVH32+DVH23+DVH7+DVH12+DVH13+
     *         DVH31+DVH30+DVH26+
     *          DVH29+DVH10)  
        DW1(IA3)=R1*(DWH8+DWH15+DWH33+DWH27)+
     *     R2*(DWH11+DWH16+DWH9+DWH28)-
     *     R3*(DWH3+DWH14+DWH32+DWH23+DWH7+DWH12+DWH13+
     *         DWH31+DWH30+DWH26+
     *          DWH29+DWH10)  
      ENDIF
C
C *** The area IA4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DU1(IA4)=(A1*(DUH13+DUH31+DUH36+DUH19)+
     *   A2*(DUH32+DUH14+DUH16+DUH20)-
     *   A3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+
     *       DUH18+DUH4+DUH9+
     *       DUH28+DUH24))*AVOL11  
        DV1(IA4)=(A1*(DVH13+DVH31+DVH36+DVH19)+
     *   A2*(DVH32+DVH14+DVH16+DVH20)-
     *   A3*(DVH17+DVH12+DVH15+DVH33+DVH30+DVH34+DVH35+
     *       DVH18+DVH4+DVH9+
     *       DVH28+DVH24))*AVOL11  
        DW1(IA4)=(A1*(DWH13+DWH31+DWH36+DWH19)+
     *   A2*(DWH32+DWH14+DWH16+DWH20)-
     *   A3*(DWH17+DWH12+DWH15+DWH33+DWH30+DWH34+DWH35+
     *       DWH18+DWH4+DWH9+
     *       DWH28+DWH24))*AVOL11  
      ELSE
        DU1(IA4)=DU1(IA4)+
     *   (A1*(DUH13+DUH31+DUH36+DUH19)+A2*(DUH32+DUH14+DUH16+DUH20)-
     *   A3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+
     *       DUH18+DUH4+DUH9+
     *       DUH28+DUH24))*AVOL11  
        DV1(IA4)= DV1(IA4)+
     *   (A1*(DVH13+DVH31+DVH36+DVH19)+A2*(DVH32+DVH14+DVH16+DVH20)-
     *   A3*(DVH17+DVH12+DVH15+DVH33+DVH30+DVH34+DVH35+
     *       DVH18+DVH4+DVH9+
     *       DVH28+DVH24))*AVOL11  
        DW1(IA4)=DW1(IA4)+
     *   (A1*(DWH13+DWH31+DWH36+DWH19)+A2*(DWH32+DWH14+DWH16+DWH20)-
     *   A3*(DWH17+DWH12+DWH15+DWH33+DWH30+DWH34+DWH35+
     *       DWH18+DWH4+DWH9+
     *       DWH28+DWH24))*AVOL11
        ENDIF
        ELSE  
C     case of a boundary area
        DU1(IA4)=R1*(DUH13+DUH31+DUH36+DUH19)+
     *   R2*(DUH32+DUH14+DUH16+DUH20)-
     *   R3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+
     *       DUH18+DUH4+DUH9+
     *       DUH28+DUH24)  
        DV1(IA4)=R1*(DVH13+DVH31+DVH36+DVH19)+
     *   R2*(DVH32+DVH14+DVH16+DVH20)-
     *   R3*(DVH17+DVH12+DVH15+DVH33+DVH30+DVH34+DVH35+
     *       DVH18+DVH4+DVH9+
     *       DVH28+DVH24)  
        DW1(IA4)=R1*(DWH13+DWH31+DWH36+DWH19)+
     *   R2*(DWH32+DWH14+DWH16+DWH20)-
     *   R3*(DWH17+DWH12+DWH15+DWH33+DWH30+DWH34+DWH35+
     *       DWH18+DWH4+DWH9+
     *       DWH28+DWH24)  
      ENDIF
C
C *** The area IA5
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(5,IEL1).GT.IEL1) THEN
        DU1(IA5)=(A1*(DUH5+DUH18+DUH35+DUH25)+
     *    A2*(DUH6+DUH20+DUH4+DUH24)-
     *    A3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+
     *        DUH2+DUH3+DUH14+
     *        DUH32+DUH23))*AVOL11  
        DV1(IA5)=(A1*(DVH5+DVH18+DVH35+DVH25)+
     *    A2*(DVH6+DVH20+DVH4+DVH24)-
     *    A3*(DVH1+DVH17+DVH19+DVH36+DVH34+DVH21+DVH22+
     *        DVH2+DVH3+DVH14+
     *        DVH32+DVH23))*AVOL11  
        DW1(IA5)=(A1*(DWH5+DWH18+DWH35+DWH25)+
     *    A2*(DWH6+DWH20+DWH4+DWH24)-
     *    A3*(DWH1+DWH17+DWH19+DWH36+DWH34+DWH21+DWH22+
     *        DWH2+DWH3+DWH14+
     *        DWH32+DWH23))*AVOL11  
      ELSE
        DU1(IA5)=DU1(IA5)+
     *    (A1*(DUH5+DUH18+DUH35+DUH25)+A2*(DUH6+DUH20+DUH4+DUH24)-
     *    A3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+
     *        DUH2+DUH3+DUH14+
     *        DUH32+DUH23))*AVOL11  
        DV1(IA5)=DV1(IA5)+
     *    (A1*(DVH5+DVH18+DVH35+DVH25)+A2*(DVH6+DVH20+DVH4+DVH24)-
     *    A3*(DVH1+DVH17+DVH19+DVH36+DVH34+DVH21+DVH22+
     *        DVH2+DVH3+DVH14+
     *        DVH32+DVH23))*AVOL11  
        DW1(IA5)= DW1(IA5)+
     *    (A1*(DWH5+DWH18+DWH35+DWH25)+A2*(DWH6+DWH20+DWH4+DWH24)-
     *    A3*(DWH1+DWH17+DWH19+DWH36+DWH34+DWH21+DWH22+
     *        DWH2+DWH3+DWH14+
     *        DWH32+DWH23))*AVOL11
       ENDIF
       ELSE  
C     case of a boundary area
        DU1(IA5)=R1*(DUH5+DUH18+DUH35+DUH25)+
     *    R2*(DUH6+DUH20+DUH4+DUH24)-
     *    R3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+
     *        DUH2+DUH3+DUH14+
     *        DUH32+DUH23)  
        DV1(IA5)=R1*(DVH5+DVH18+DVH35+DVH25)+
     *    R2*(DVH6+DVH20+DVH4+DVH24)-
     *    R3*(DVH1+DVH17+DVH19+DVH36+DVH34+DVH21+DVH22+
     *        DVH2+DVH3+DVH14+
     *        DVH32+DVH23)  
        DW1(IA5)=R1*(DWH5+DWH18+DWH35+DWH25)+
     *    R2*(DWH6+DWH20+DWH4+DWH24)-
     *    R3*(DWH1+DWH17+DWH19+DWH36+DWH34+DWH21+DWH22+
     *        DWH2+DWH3+DWH14+
     *        DWH32+DWH23)  
      ENDIF
C
C *** The area IA6
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(6,IEL1).GT.IEL1) THEN
        DU1(IA6)=(A1*(DUH26+DUH30+DUH34+DUH21)+
     *    A2*(DUH23+DUH32+DUH24+DUH28)-
     *    A3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+
     *        DUH25+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11  
        DV1(IA6)=(A1*(DVH26+DVH30+DVH34+DVH21)+
     *    A2*(DVH23+DVH32+DVH24+DVH28)-
     *    A3*(DVH22+DVH29+DVH27+DVH33+DVH31+DVH36+DVH21+
     *        DVH25+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11  
        DW1(IA6)=(A1*(DWH26+DWH30+DWH34+DWH21)+
     *    A2*(DWH23+DWH32+DWH24+DWH28)-
     *    A3*(DWH22+DWH29+DWH27+DWH33+DWH31+DWH36+DWH21+
     *        DWH25+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11  
      ELSE
        DU1(IA6)= DU1(IA6)+
     *    (A1*(DUH26+DUH30+DUH34+DUH21)+A2*(DUH23+DUH32+DUH24+DUH28)-
     *    A3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+
     *        DUH25+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11  
        DV1(IA6)=DV1(IA6)+
     *    (A1*(DVH26+DVH30+DVH34+DVH21)+A2*(DVH23+DVH32+DVH24+DVH28)-
     *    A3*(DVH22+DVH29+DVH27+DVH33+DVH31+DVH36+DVH21+
     *        DVH25+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11  
        DW1(IA6)=DW1(IA6)+
     *    (A1*(DWH26+DWH30+DWH34+DWH21)+A2*(DWH23+DWH32+DWH24+DWH28)-
     *    A3*(DWH22+DWH29+DWH27+DWH33+DWH31+DWH36+DWH21+
     *        DWH25+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11 
        ENDIF
        ELSE 
C     case of a boundary area
        DU1(IA6)=R1*(DUH26+DUH30+DUH34+DUH21)+
     *    R2*(DUH23+DUH32+DUH24+DUH28)-
     *    R3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+
     *        DUH25+DUH6+DUH11+
     *        DUH16+DUH20)  
        DV1(IA6)=R1*(DVH26+DVH30+DVH34+DVH21)+
     *    R2*(DVH23+DVH32+DVH24+DVH28)-
     *    R3*(DVH22+DVH29+DVH27+DVH33+DVH31+DVH36+DVH21+
     *        DVH25+DVH6+DVH11+
     *        DVH16+DVH20)  
        DW1(IA6)=R1*(DWH26+DWH30+DWH34+DWH21)+
     *    R2*(DWH23+DWH32+DWH24+DWH28)-
     *    R3*(DWH22+DWH29+DWH27+DWH33+DWH31+DWH36+DWH21+
     *        DWH25+DWH6+DWH11+
     *        DWH16+DWH20)  
      ENDIF
C
10    CONTINUE
C
      DO 555 IEL1=1,NEL1
      DO 666 I=1,6
      IADJ=KADJ1(I,IEL1)
      IA=KAREA1(I,IEL1)
      IF (IADJ.EQ.0) THEN
      GOTO 666
      ENDIF
      IF (IADJ.GT.IEL1) THEN
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DU1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DU1(IA)
      DV1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DV1(IA)
      DW1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DW1(IA)
      ELSE
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DU1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DU1(IA)
      DV1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DV1(IA)
      DW1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DW1(IA)
      ENDIF
666   CONTINUE
555   CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE   PROLU  (DU1,DV1,DW1,DP1,DU2,DV2,DW2,DP2,
     *                     KAREA1,KADJ1,KVERT1,NAT1,NEL1,AVOL1,
     *                     KAREA2,KADJ2,KVERT2,NAT2,NEL2,AVOL2,
     *                     VPL1,VPL2,VAUX,NVT1,NVT2,IINT)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1,DW1,DP1) 
*                to the fine grid vector (DU2,DV2,DW2,DP2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1,DW1,DP1           - coarse grid vector
*      KAREA1,KADJ1,KVERT1,NAT1,NEL1,NVT1,AVOL1 - data of the coarse grid
*      KAREA2,KADJ2,KVERT2,NAT2,NEL2,NVT2,AVOL2 - data of the fine grid
*
*    Output:
*      DU2,DV2,DW2,DP2           - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL AVOL1,AVOL2
      REAL VPL1,VPL2,VAUX
C
      PARAMETER (NNAE=6,NNVE=8)
      DIMENSION DU1(1),DV1(1),DW1(1),DP1(1),DU2(1),DV2(1),DW2(1),DP2(1),
     *          KAREA1(NNAE,*),KAREA2(NNAE,*),KVERT1(NNVE,*),
     *          KVERT2(NNVE,*),VPL1(*),VPL2(*),VAUX(*),
     *          KADJ1(NNAE,*),KADJ2(NNAE,*),AVOL1(*),AVOL2(*)
      SAVE
C-----------------------------------------------------------------------
C
      IF (IINT.GT.0) THEN
       A1= 11D0/24D0
       A2=  7D0/48D0
       A3=- 5D0/48D0
       A4=- 1D0/24D0
       A5= 11D0/24D0
       A6=  1D0/12D0
       A7=- 1D0/24D0
      ELSE
       A1= 0.500D0
       A2= 0.125D0
       A3=-0.125D0
       A4= 0.000D0
       A5= 0.500D0
       A6= 0.000D0
       A7= 0.000D0
      ENDIF
C
C
      IF ((IINT.EQ.2).OR.(IINT.EQ.-2)) THEN
C     
C     LINEAR Prolongation of the pressure
C
      DO 111 IEL1=1,NEL1
      VPIEL=REAL(DP1(IEL1))
      VAVOL=AVOL1(IEL1)
C
      IV1=KVERT1(1,IEL1)
      IV2=KVERT1(2,IEL1)
      IV3=KVERT1(3,IEL1)
      IV4=KVERT1(4,IEL1)
      IV5=KVERT1(5,IEL1)
      IV6=KVERT1(6,IEL1)
      IV7=KVERT1(7,IEL1)
      IV8=KVERT1(8,IEL1)
C
      VPL1(IV1)=VPL1(IV1)+0.125E0*VAVOL*VPIEL
      VPL1(IV2)=VPL1(IV2)+0.125E0*VAVOL*VPIEL
      VPL1(IV3)=VPL1(IV3)+0.125E0*VAVOL*VPIEL
      VPL1(IV4)=VPL1(IV4)+0.125E0*VAVOL*VPIEL
      VPL1(IV5)=VPL1(IV5)+0.125E0*VAVOL*VPIEL
      VPL1(IV6)=VPL1(IV6)+0.125E0*VAVOL*VPIEL
      VPL1(IV7)=VPL1(IV7)+0.125E0*VAVOL*VPIEL
      VPL1(IV8)=VPL1(IV8)+0.125E0*VAVOL*VPIEL
C
      VAUX(IV1)=VAUX(IV1)+0.125E0*VAVOL
      VAUX(IV2)=VAUX(IV2)+0.125E0*VAVOL
      VAUX(IV3)=VAUX(IV3)+0.125E0*VAVOL
      VAUX(IV4)=VAUX(IV4)+0.125E0*VAVOL
      VAUX(IV5)=VAUX(IV5)+0.125E0*VAVOL
      VAUX(IV6)=VAUX(IV6)+0.125E0*VAVOL
      VAUX(IV7)=VAUX(IV7)+0.125E0*VAVOL
      VAUX(IV8)=VAUX(IV8)+0.125E0*VAVOL
111   CONTINUE
C
      DO 20 IVT1=1,NVT1
20    VPL1(IVT1)=VPL1(IVT1)/VAUX(IVT1)
C
      CALL LCP2(VPL1,VPL2,NVT1)
C
      DO 112 IEL1=1,NEL1
C
      VXH1=0.5D0*VPL1(KVERT1(1,IEL1))
      VXH2=0.5D0*VPL1(KVERT1(2,IEL1))
      VXH3=0.5D0*VPL1(KVERT1(3,IEL1))
      VXH4=0.5D0*VPL1(KVERT1(4,IEL1))
      VXH5=0.5D0*VPL1(KVERT1(5,IEL1))
      VXH6=0.5D0*VPL1(KVERT1(6,IEL1))
      VXH7=0.5D0*VPL1(KVERT1(7,IEL1))
      VXH8=0.5D0*VPL1(KVERT1(8,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      IF ((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)) 
     *   VPL2(KVERT2(3,IELH1))=0.5D0*(VXH1+VXH2+VXH3+VXH4)
C
      IF ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0)) 
     *   VPL2(KVERT2(6,IELH1))=0.5D0*(VXH1+VXH2+VXH5+VXH6)
C
      IF ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0)) 
     *   VPL2(KVERT2(6,IELH2))=0.5D0*(VXH2+VXH3+VXH6+VXH7)
C
      IF ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0)) 
     *   VPL2(KVERT2(6,IELH3))=0.5D0*(VXH3+VXH4+VXH7+VXH8)
C
      IF ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0)) 
     *   VPL2(KVERT2(6,IELH4))=0.5D0*(VXH1+VXH4+VXH5+VXH8)
C
      IF ((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)) 
     *   VPL2(KVERT2(3,IELH5))=0.5D0*(VXH5+VXH6+VXH7+VXH8)
C
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH1))=VXH1+VXH2
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH2))=VXH2+VXH3
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH3))=VXH3+VXH4
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH4))=VXH1+VXH4
C
      IF (((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0)).AND.
     *    ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0))) 
     *   VPL2(KVERT2(5,IELH1))=VXH1+VXH5
C
      IF (((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0)).AND.
     *    ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0))) 
     *   VPL2(KVERT2(5,IELH2))=VXH2+VXH6
C
      IF (((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0)).AND.
     *    ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0))) 
     *   VPL2(KVERT2(5,IELH3))=VXH3+VXH7
C
      IF (((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0)).AND.
     *    ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0))) 
     *   VPL2(KVERT2(5,IELH4))=VXH4+VXH8
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH5))=VXH5+VXH6
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH6))=VXH6+VXH7
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH7))=VXH7+VXH8
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0))) 
     *   VPL2(KVERT2(2,IELH8))=VXH5+VXH8
C
      VPL2(KVERT2(7,IELH1))=0.25D0*(VXH1+VXH2+VXH3+VXH4+
     *                             VXH5+VXH6+VXH7+VXH8)
C
112   CONTINUE
C
      DO 113 IEL2=1,NEL2
      V1=VPL2(KVERT2(1,IEL2))
      V2=VPL2(KVERT2(2,IEL2))
      V3=VPL2(KVERT2(3,IEL2))
      V4=VPL2(KVERT2(4,IEL2))
      V5=VPL2(KVERT2(5,IEL2))
      V6=VPL2(KVERT2(6,IEL2))
      V7=VPL2(KVERT2(7,IEL2))
      V8=VPL2(KVERT2(8,IEL2))
      DP2(IEL2)=0.125D0*(V1+V2+V3+V4+V5+V6+V7+V8)
C
113   CONTINUE
C
      ENDIF
C
C *** Zero initialization of (DU2,DV2,DW2)
      CALL  LCL1 (DU2,NAT2)
      CALL  LCL1 (DV2,NAT2)
      CALL  LCL1 (DW2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DUH1=DU1(KAREA1(1,IEL1))
      DUH2=DU1(KAREA1(2,IEL1))
      DUH3=DU1(KAREA1(3,IEL1))
      DUH4=DU1(KAREA1(4,IEL1))
      DUH5=DU1(KAREA1(5,IEL1))
      DUH6=DU1(KAREA1(6,IEL1))

C
      DVH1=DV1(KAREA1(1,IEL1))
      DVH2=DV1(KAREA1(2,IEL1))
      DVH3=DV1(KAREA1(3,IEL1))
      DVH4=DV1(KAREA1(4,IEL1))
      DVH5=DV1(KAREA1(5,IEL1))
      DVH6=DV1(KAREA1(6,IEL1))
C
      DWH1=DW1(KAREA1(1,IEL1))
      DWH2=DW1(KAREA1(2,IEL1))
      DWH3=DW1(KAREA1(3,IEL1))
      DWH4=DW1(KAREA1(4,IEL1))
      DWH5=DW1(KAREA1(5,IEL1))
      DWH6=DW1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      IF ((IINT.EQ.1).OR.(IINT.EQ.-1)) THEN
C
C ***  constant Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
      DP2(IELH5)=DPH
      DP2(IELH6)=DPH
      DP2(IELH7)=DPH
      DP2(IELH8)=DPH
C
      ENDIF
C
      B11=A1*DUH1
      B12=A2*DUH1
      B13=A3*DUH1
      B14=A4*DUH1
      B15=A5*DUH1
      B16=A6*DUH1
      B17=A7*DUH1
C
      C11=A1*DUH2
      C12=A2*DUH2
      C13=A3*DUH2
      C14=A4*DUH2
      C15=A5*DUH2
      C16=A6*DUH2
      C17=A7*DUH2
C
      D11=A1*DUH3
      D12=A2*DUH3
      D13=A3*DUH3
      D14=A4*DUH3
      D15=A5*DUH3
      D16=A6*DUH3
      D17=A7*DUH3
C
      E11=A1*DUH4
      E12=A2*DUH4
      E13=A3*DUH4
      E14=A4*DUH4
      E15=A5*DUH4
      E16=A6*DUH4
      E17=A7*DUH4
C
      F11=A1*DUH5
      F12=A2*DUH5
      F13=A3*DUH5
      F14=A4*DUH5
      F15=A5*DUH5
      F16=A6*DUH5
      F17=A7*DUH5
C
      G11=A1*DUH6
      G12=A2*DUH6
      G13=A3*DUH6
      G14=A4*DUH6
      G15=A5*DUH6
      G16=A6*DUH6
      G17=A7*DUH6
      R1=B11+G14
      R2=C11+E14
      R3=D11+F14
      R4=C14+E11
      R5=F11+D14
      R6=B14+G11
      R7=B15+G17
      R8=B16+G16
      R9=B17+G15
C
C
C
      AVOL11=AVOL1(IEL1)
C
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DU2(KAREA2(1,IELH1))=DU2(KAREA2(1,IELH1))+
     *                      (R1+C12+D13+
     *                      E13+F12)*AVOL11
       DU2(KAREA2(1,IELH2))=DU2(KAREA2(1,IELH2))+
     *                      (R1+C12+D12+
     *                      E13+F13)*AVOL11
       DU2(KAREA2(1,IELH3))=DU2(KAREA2(1,IELH3))+
     *                      (R1+C13+D12+
     *                      E12+F13)*AVOL11
       DU2(KAREA2(1,IELH4))=DU2(KAREA2(1,IELH4))+
     *                      (R1+C13+D13+
     *                      E12+F12)*AVOL11
      ENDIF
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH1))=DU2(KAREA2(2,IELH1))+
     *                      (B12+R2+D13+
     *                      F12+G13)*AVOL11
       DU2(KAREA2(5,IELH2))=DU2(KAREA2(5,IELH2))+
     *                      (B12+R2+D12+
     *                      F13+G13)*AVOL11
       DU2(KAREA2(5,IELH6))=DU2(KAREA2(5,IELH6))+
     *                      (B13+R2+D12+
     *                      F13+G12)*AVOL11
       DU2(KAREA2(2,IELH5))=DU2(KAREA2(2,IELH5))+
     *                      (B13+R2+D13+
     *                      F12+G12)*AVOL11
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH2))=DU2(KAREA2(2,IELH2))+
     *                      (B12+C12+R3+
     *                      E13+G13)*AVOL11
       DU2(KAREA2(5,IELH3))=DU2(KAREA2(5,IELH3))+
     *                      (B12+C13+R3+
     *                      E12+G13)*AVOL11
       DU2(KAREA2(5,IELH7))=DU2(KAREA2(5,IELH7))+
     *                      (B13+C13+R3+
     *                      E12+G12)*AVOL11
       DU2(KAREA2(2,IELH6))=DU2(KAREA2(2,IELH6))+
     *                      (B13+C12+R3+
     *                      E13+G12)*AVOL11
      ENDIF
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH3))=DU2(KAREA2(2,IELH3))+
     *                      (B12+R4+D12+
     *                      F13+G13)*AVOL11
       DU2(KAREA2(5,IELH4))=DU2(KAREA2(5,IELH4))+
     *                      (B12+R4+D13+
     *                      F12+G13)*AVOL11
       DU2(KAREA2(5,IELH8))=DU2(KAREA2(5,IELH8))+
     *                      (B13+R4+D13+
     *                      F12+G12)*AVOL11
       DU2(KAREA2(2,IELH7))=DU2(KAREA2(2,IELH7))+
     *                      (B13+R4+D12+
     *                      F13+G12)*AVOL11
      ENDIF
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH4))=DU2(KAREA2(2,IELH4))+
     *                      (B12+C13+R5+
     *                      E12+G13)*AVOL11
       DU2(KAREA2(5,IELH1))=DU2(KAREA2(5,IELH1))+
     *                      (B12+C12+R5+
     *                      E13+G13)*AVOL11
       DU2(KAREA2(5,IELH5))=DU2(KAREA2(5,IELH5))+
     *                      (B13+C12+R5+
     *                      E13+G12)*AVOL11
       DU2(KAREA2(2,IELH8))=DU2(KAREA2(2,IELH8))+
     *                      (B13+C13+R5+
     *                      E12+G12)*AVOL11
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DU2(KAREA2(1,IELH5))=DU2(KAREA2(1,IELH5))+
     *                      (R6+C12+D13+
     *                      E13+F12)*AVOL11
       DU2(KAREA2(1,IELH6))=DU2(KAREA2(1,IELH6))+
     *                      (R6+C12+D12+
     *                      E13+F13)*AVOL11
       DU2(KAREA2(1,IELH7))=DU2(KAREA2(1,IELH7))+
     *                      (R6+C13+D12+
     *                      E12+F13)*AVOL11
       DU2(KAREA2(1,IELH8))=DU2(KAREA2(1,IELH8))+
     *                      (R6+C13+D13+
     *                      E12+F12)*AVOL11
      ENDIF
C
C
      DU2(KAREA2(3,IELH1))=R7+C15+D16+
     *                     E17+F16
      DU2(KAREA2(3,IELH2))=R7+C16+D15+
     *                     E16+F17
      DU2(KAREA2(3,IELH3))=R7+C17+D16+
     *                     E15+F16
      DU2(KAREA2(3,IELH4))=R7+C16+D17+
     *                     E16+F15
C
      DU2(KAREA2(6,IELH1))=R8+C15+D17+
     *                     E17+F15
      DU2(KAREA2(6,IELH2))=R8+C15+D15+
     *                     E17+F17
      DU2(KAREA2(6,IELH3))=R8+C17+D15+
     *                     E15+F17
      DU2(KAREA2(6,IELH4))=R8+C17+D17+
     *                     E15+F15
C
      DU2(KAREA2(3,IELH5))=R9+C15+D16+
     *                     E17+F16
      DU2(KAREA2(3,IELH6))=R9+C16+D15+
     *                     E16+F17
      DU2(KAREA2(3,IELH7))=R9+C17+D16+
     *                     E15+F16
      DU2(KAREA2(3,IELH8))=R9+C16+D17+
     *                     E16+F15
C
      B21=A1*DVH1
      B22=A2*DVH1
      B23=A3*DVH1
      B24=A4*DVH1
      B25=A5*DVH1
      B26=A6*DVH1
      B27=A7*DVH1
C
      C21=A1*DVH2
      C22=A2*DVH2
      C23=A3*DVH2
      C24=A4*DVH2
      C25=A5*DVH2
      C26=A6*DVH2
      C27=A7*DVH2
C
      D21=A1*DVH3
      D22=A2*DVH3
      D23=A3*DVH3
      D24=A4*DVH3
      D25=A5*DVH3
      D26=A6*DVH3
      D27=A7*DVH3
C
      E21=A1*DVH4
      E22=A2*DVH4
      E23=A3*DVH4
      E24=A4*DVH4
      E25=A5*DVH4
      E26=A6*DVH4
      E27=A7*DVH4
C
      F21=A1*DVH5
      F22=A2*DVH5
      F23=A3*DVH5
      F24=A4*DVH5
      F25=A5*DVH5
      F26=A6*DVH5
      F27=A7*DVH5
C
      G21=A1*DVH6
      G22=A2*DVH6
      G23=A3*DVH6
      G24=A4*DVH6
      G25=A5*DVH6
      G26=A6*DVH6
      G27=A7*DVH6
C
C
      S1=B21+G24
      S2=C21+E24
      S3=D21+F24
      S4=C24+E21
      S5=D24+F21
      S6=B24+G21
      S7=B25+G27
      S8=B26+G26
      S9=B27+G25
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DV2(KAREA2(1,IELH1))=DV2(KAREA2(1,IELH1))+
     *                      (S1+C22+D23+
     *                      E23+F22)*AVOL11
       DV2(KAREA2(1,IELH2))=DV2(KAREA2(1,IELH2))+
     *                      (S1+C22+D22+
     *                      E23+F23)*AVOL11
       DV2(KAREA2(1,IELH3))=DV2(KAREA2(1,IELH3))+
     *                      (S1+C23+D22+
     *                      E22+F23)*AVOL11
       DV2(KAREA2(1,IELH4))=DV2(KAREA2(1,IELH4))+
     *                      (S1+C23+D23+
     *                      E22+F22)*AVOL11
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH1))=DV2(KAREA2(2,IELH1))+
     *                      (B22+S2+D23+
     *                      F22+G23)*AVOL11
       DV2(KAREA2(5,IELH2))=DV2(KAREA2(5,IELH2))+
     *                      (B22+S2+D22+
     *                      F23+G23)*AVOL11
       DV2(KAREA2(5,IELH6))=DV2(KAREA2(5,IELH6))+
     *                      (B23+S2+D22+
     *                      F23+G22)*AVOL11
       DV2(KAREA2(2,IELH5))=DV2(KAREA2(2,IELH5))+
     *                      (B23+S2+D23+
     *                      F22+G22)*AVOL11
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH2))=DV2(KAREA2(2,IELH2))+
     *                      (B22+C22+S3+
     *                      E23+G23)*AVOL11
       DV2(KAREA2(5,IELH3))=DV2(KAREA2(5,IELH3))+
     *                      (B22+C23+S3+
     *                      E22+G23)*AVOL11
       DV2(KAREA2(5,IELH7))=DV2(KAREA2(5,IELH7))+
     *                      (B23+C23+S3+
     *                      E22+G22)*AVOL11
       DV2(KAREA2(2,IELH6))=DV2(KAREA2(2,IELH6))+
     *                      (B23+C22+S3+
     *                      E23+G22)*AVOL11
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH3))=DV2(KAREA2(2,IELH3))+
     *                      (B22+S4+D22+
     *                      F23+G23)*AVOL11
       DV2(KAREA2(5,IELH4))=DV2(KAREA2(5,IELH4))+
     *                      (B22+S4+D23+
     *                      F22+G23)*AVOL11
       DV2(KAREA2(5,IELH8))=DV2(KAREA2(5,IELH8))+
     *                      (B23+S4+D23+
     *                      F22+G22)*AVOL11
       DV2(KAREA2(2,IELH7))=DV2(KAREA2(2,IELH7))+
     *                      (B23+S4+D22+
     *                      F23+G22)*AVOL11
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH4))=DV2(KAREA2(2,IELH4))+
     *                      (B22+C23+S5+
     *                      E22+G23)*AVOL11
       DV2(KAREA2(5,IELH1))=DV2(KAREA2(5,IELH1))+
     *                      (B22+C22+S5+
     *                      E23+G23)*AVOL11
       DV2(KAREA2(5,IELH5))=DV2(KAREA2(5,IELH5))+
     *                      (B23+C22+S5+
     *                      E23+G22)*AVOL11
       DV2(KAREA2(2,IELH8))=DV2(KAREA2(2,IELH8))+
     *                      (B23+C23+S5+
     *                      E22+G22)*AVOL11
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DV2(KAREA2(1,IELH5))=DV2(KAREA2(1,IELH5))+
     *                      (S6+C22+D23+
     *                      E23+F22)*AVOL11
       DV2(KAREA2(1,IELH6))=DV2(KAREA2(1,IELH6))+
     *                      (S6+C22+D22+
     *                      E23+F23)*AVOL11
       DV2(KAREA2(1,IELH7))=DV2(KAREA2(1,IELH7))+
     *                      (S6+C23+D22+
     *                      E22+F23)*AVOL11
       DV2(KAREA2(1,IELH8))=DV2(KAREA2(1,IELH8))+
     *                      (S6+C23+D23+
     *                      E22+F22)*AVOL11
      ENDIF
C
C
      DV2(KAREA2(3,IELH1))=S7+C25+D26+
     *                     E27+F26
      DV2(KAREA2(3,IELH2))=S7+C26+D25+
     *                     E26+F27
      DV2(KAREA2(3,IELH3))=S7+C27+D26+
     *                     E25+F26
      DV2(KAREA2(3,IELH4))=S7+C26+D27+
     *                     E26+F25
C
      DV2(KAREA2(6,IELH1))=S8+C25+D27+
     *                     E27+F25
      DV2(KAREA2(6,IELH2))=S8+C25+D25+
     *                     E27+F27
      DV2(KAREA2(6,IELH3))=S8+C27+D25+
     *                     E25+F27
      DV2(KAREA2(6,IELH4))=S8+C27+D27+
     *                     E25+F25
C
      DV2(KAREA2(3,IELH5))=S9+C25+D26+
     *                     E27+F26
      DV2(KAREA2(3,IELH6))=S9+C26+D25+
     *                     E26+F27
      DV2(KAREA2(3,IELH7))=S9+C27+D26+
     *                     E25+F26
      DV2(KAREA2(3,IELH8))=S9+C26+D27+
     *                     E26+F25
C
      B31=A1*DWH1
      B32=A2*DWH1
      B33=A3*DWH1
      B34=A4*DWH1
      B35=A5*DWH1
      B36=A6*DWH1
      B37=A7*DWH1
C
      C31=A1*DWH2
      C32=A2*DWH2
      C33=A3*DWH2
      C34=A4*DWH2
      C35=A5*DWH2
      C36=A6*DWH2
      C37=A7*DWH2
C
      D31=A1*DWH3
      D32=A2*DWH3
      D33=A3*DWH3
      D34=A4*DWH3
      D35=A5*DWH3
      D36=A6*DWH3
      D37=A7*DWH3
C
      E31=A1*DWH4
      E32=A2*DWH4
      E33=A3*DWH4
      E34=A4*DWH4
      E35=A5*DWH4
      E36=A6*DWH4
      E37=A7*DWH4
C
      F31=A1*DWH5
      F32=A2*DWH5
      F33=A3*DWH5
      F34=A4*DWH5
      F35=A5*DWH5
      F36=A6*DWH5
      F37=A7*DWH5
C
      G31=A1*DWH6
      G32=A2*DWH6
      G33=A3*DWH6
      G34=A4*DWH6
      G35=A5*DWH6
      G36=A6*DWH6
      G37=A7*DWH6
C
C
      T1=B31+G34
      T2=C31+E34
      T3=D31+F34
      T4=C34+E31
      T5=D34+F31
      T6=B34+G31
      T7=B35+G37
      T8=B36+G36
      T9=B37+G35
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DW2(KAREA2(1,IELH1))=DW2(KAREA2(1,IELH1))+
     *                      (T1+C32+D33+
     *                      E33+F32)*AVOL11
       DW2(KAREA2(1,IELH2))=DW2(KAREA2(1,IELH2))+
     *                      (T1+C32+D32+
     *                      E33+F33)*AVOL11
       DW2(KAREA2(1,IELH3))=DW2(KAREA2(1,IELH3))+
     *                      (T1+C33+D32+
     *                      E32+F33)*AVOL11
       DW2(KAREA2(1,IELH4))=DW2(KAREA2(1,IELH4))+
     *                      (T1+C33+D33+
     *                      E32+F32)*AVOL11
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH1))=DW2(KAREA2(2,IELH1))+
     *                      (B32+T2+D33+
     *                      F32+G33)*AVOL11
       DW2(KAREA2(5,IELH2))=DW2(KAREA2(5,IELH2))+
     *                      (B32+T2+D32+
     *                      F33+G33)*AVOL11
       DW2(KAREA2(5,IELH6))=DW2(KAREA2(5,IELH6))+
     *                      (B33+T2+D32+
     *                      F33+G32)*AVOL11
       DW2(KAREA2(2,IELH5))=DW2(KAREA2(2,IELH5))+
     *                      (B33+T2+D33+
     *                      F32+G32)*AVOL11
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH2))=DW2(KAREA2(2,IELH2))+
     *                      (B32+C32+T3+
     *                      E33+G33)*AVOL11
       DW2(KAREA2(5,IELH3))=DW2(KAREA2(5,IELH3))+
     *                      (B32+C33+T3+
     *                      E32+G33)*AVOL11
       DW2(KAREA2(5,IELH7))=DW2(KAREA2(5,IELH7))+
     *                      (B33+C33+T3+
     *                      E32+G32)*AVOL11
       DW2(KAREA2(2,IELH6))=DW2(KAREA2(2,IELH6))+
     *                      (B33+C32+T3+
     *                      E33+G32)*AVOL11
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH3))=DW2(KAREA2(2,IELH3))+
     *                      (B32+T4+D32+
     *                      F33+G33)*AVOL11
       DW2(KAREA2(5,IELH4))=DW2(KAREA2(5,IELH4))+
     *                      (B32+T4+D33+
     *                      F32+G33)*AVOL11
       DW2(KAREA2(5,IELH8))=DW2(KAREA2(5,IELH8))+
     *                      (B33+T4+D33+
     *                      F32+G32)*AVOL11
       DW2(KAREA2(2,IELH7))=DW2(KAREA2(2,IELH7))+
     *                      (B33+T4+D32+
     *                      F33+G32)*AVOL11
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH4))=DW2(KAREA2(2,IELH4))+
     *                      (B32+C33+T5+
     *                      E32+G33)*AVOL11
       DW2(KAREA2(5,IELH1))=DW2(KAREA2(5,IELH1))+
     *                      (B32+C32+T5+
     *                      E33+G33)*AVOL11
       DW2(KAREA2(5,IELH5))=DW2(KAREA2(5,IELH5))+
     *                      (B33+C32+T5+
     *                      E33+G32)*AVOL11
       DW2(KAREA2(2,IELH8))=DW2(KAREA2(2,IELH8))+
     *                      (B33+C33+T5+
     *                      E32+G32)*AVOL11
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DW2(KAREA2(1,IELH5))=DW2(KAREA2(1,IELH5))+
     *                      (T6+C32+D33+
     *                      E33+F32)*AVOL11
       DW2(KAREA2(1,IELH6))=DW2(KAREA2(1,IELH6))+
     *                      (T6+C32+D32+
     *                      E33+F33)*AVOL11
       DW2(KAREA2(1,IELH7))=DW2(KAREA2(1,IELH7))+
     *                      (T6+C33+D32+
     *                      E32+F33)*AVOL11
       DW2(KAREA2(1,IELH8))=DW2(KAREA2(1,IELH8))+
     *                      (T6+C33+D33+
     *                      E32+F32)*AVOL11
      ENDIF
C
C
      DW2(KAREA2(3,IELH1))=T7+C35+D36+
     *                     E37+F36
      DW2(KAREA2(3,IELH2))=T7+C36+D35+
     *                     E36+F37
      DW2(KAREA2(3,IELH3))=T7+C37+D36+
     *                     E35+F36
      DW2(KAREA2(3,IELH4))=T7+C36+D37+
     *                     E36+F35
C
      DW2(KAREA2(6,IELH1))=T8+C35+D37+
     *                     E37+F35
      DW2(KAREA2(6,IELH2))=T8+C35+D35+
     *                     E37+F37
      DW2(KAREA2(6,IELH3))=T8+C37+D35+
     *                     E35+F37
      DW2(KAREA2(6,IELH4))=T8+C37+D37+
     *                     E35+F35
C
      DW2(KAREA2(3,IELH5))=T9+C35+D36+
     *                     E37+F36
      DW2(KAREA2(3,IELH6))=T9+C36+D35+
     *                     E36+F37
      DW2(KAREA2(3,IELH7))=T9+C37+D36+
     *                     E35+F36
      DW2(KAREA2(3,IELH8))=T9+C36+D37+
     *                     E36+F35
C
C
C
10    CONTINUE
C
      DO 555 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
      AVOL22=AVOL1(IEL1)
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       AVOL23=AVOL1(KADJ1(1,IEL1))
       DU2(KAREA2(1,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH1))
       DU2(KAREA2(1,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH2))
       DU2(KAREA2(1,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH3))
       DU2(KAREA2(1,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH4))
       DV2(KAREA2(1,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH1))
       DV2(KAREA2(1,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH2))
       DV2(KAREA2(1,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH3))
       DV2(KAREA2(1,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH4))
       DW2(KAREA2(1,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH1))
       DW2(KAREA2(1,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH2))
       DW2(KAREA2(1,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH3))
       DW2(KAREA2(1,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH4))
      ENDIF
      IF (KADJ1(2,IEL1).NE.0) THEN 
       AVOL23=AVOL1(KADJ1(2,IEL1))
       DU2(KAREA2(2,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH1))
       DU2(KAREA2(5,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH2))
       DU2(KAREA2(5,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH6))
       DU2(KAREA2(2,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH5))
       DV2(KAREA2(2,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH1))
       DV2(KAREA2(5,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH2))
       DV2(KAREA2(5,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH6))
       DV2(KAREA2(2,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH5))
       DW2(KAREA2(2,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH1))
       DW2(KAREA2(5,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH2))
       DW2(KAREA2(5,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH6))
       DW2(KAREA2(2,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH5))
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       AVOL23=AVOL1(KADJ1(3,IEL1))
       DU2(KAREA2(2,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH2))
       DU2(KAREA2(5,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH3))
       DU2(KAREA2(5,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH7))
       DU2(KAREA2(2,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH6))
       DV2(KAREA2(2,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH2))
       DV2(KAREA2(5,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH3))
       DV2(KAREA2(5,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH7))
       DV2(KAREA2(2,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH6))
       DW2(KAREA2(2,IELH2))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH2))
       DW2(KAREA2(5,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH3))
       DW2(KAREA2(5,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH7))
       DW2(KAREA2(2,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH6))
      ENDIF
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       AVOL23=AVOL1(KADJ1(4,IEL1))
       DU2(KAREA2(2,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH3))
       DU2(KAREA2(5,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH4))
       DU2(KAREA2(5,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH8))
       DU2(KAREA2(2,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH7))
       DV2(KAREA2(2,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH3))
       DV2(KAREA2(5,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH4))
       DV2(KAREA2(5,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH8))
       DV2(KAREA2(2,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH7))
       DW2(KAREA2(2,IELH3))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH3))
       DW2(KAREA2(5,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH4))
       DW2(KAREA2(5,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH8))
       DW2(KAREA2(2,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH7))
      ENDIF
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       AVOL23=AVOL1(KADJ1(5,IEL1))
       DU2(KAREA2(2,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH4))
       DU2(KAREA2(5,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH1))
       DU2(KAREA2(5,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(5,IELH5))
       DU2(KAREA2(2,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(2,IELH8))
       DV2(KAREA2(2,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH4))
       DV2(KAREA2(5,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH1))
       DV2(KAREA2(5,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(5,IELH5))
       DV2(KAREA2(2,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(2,IELH8))
       DW2(KAREA2(2,IELH4))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH4))
       DW2(KAREA2(5,IELH1))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH1))
       DW2(KAREA2(5,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(5,IELH5))
       DW2(KAREA2(2,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(2,IELH8))
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       AVOL23=AVOL1(KADJ1(6,IEL1))
       DU2(KAREA2(1,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH5))
       DU2(KAREA2(1,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH6))
       DU2(KAREA2(1,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH7))
       DU2(KAREA2(1,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DU2(KAREA2(1,IELH8))
       DV2(KAREA2(1,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH5))
       DV2(KAREA2(1,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH6))
       DV2(KAREA2(1,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH7))
       DV2(KAREA2(1,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DV2(KAREA2(1,IELH8))
       DW2(KAREA2(1,IELH5))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH5))
       DW2(KAREA2(1,IELH6))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH6))
       DW2(KAREA2(1,IELH7))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH7))
       DW2(KAREA2(1,IELH8))=SQRT(2D0)/(SQRT(AVOL22+AVOL23))*
     *                      DW2(KAREA2(1,IELH8))
      ENDIF
C
555   CONTINUE 
      END
c
************************************************************************
      SUBROUTINE  RESTRD  (DU1,DV1,DW1,DP1,DU2,DV2,DW2,DP2,
     *                     KAREA1,KADJ1,KVERT1,NAT1,NEL1,
     *                     KAREA2,KADJ2,KVERT2,NAT2,NEL2,IINT,
     *                     NVT1,NVT2,VPL1,VPL2,VAUX,AVOL,AVOL3)
************************************************************************
*    Purpose:    Restricts the fine grid defect vector (DU2,DV2,DW2,DP2)
*                to the coarse grid defect vector (DU1,DV1,DW1,DP1)
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2,DW2,DP2                - fine grid defect vector
*      KAREA1,KADJ1,KVERT1,NAT1,NEL1,NVT1,AVOL - data of the coarse grid
*      KAREA2,KADJ2,KVERT2,NAT2,NEL2,NVT2,AVOL3  - data of the fine grid
*
*    Output:
*      DU1,DV1,DW1,DP1           - coarse grid defect vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL AVOL,AVOL3
      REAL VPL1,VPL2,VAUX
C
      PARAMETER (NNAE=6,NNVE=8)
      DIMENSION DU1(1),DV1(1),DW1(1),DP1(1),DU2(1),DV2(1),DW2(1),DP2(1),
     *          KAREA2(NNAE,*),KAREA1(NNAE,*),KVERT1(NNVE,*),
     *          KVERT2(NNVE,*),VPL1(*),VPL2(*),VAUX(*),
     *          KADJ1(NNAE,*),KADJ2(NNAE,*),AVOL(*),AVOL3(*)
C
      SAVE
C-----------------------------------------------------------------------
C
      IF (IINT.GT.0) THEN
       A1= 11D0/24D0
       A2=  7D0/48D0
       A3=- 5D0/48D0
       A4=- 1D0/24D0
C       A5= 11D0/24D0
       A6=  1D0/12D0
C       A7=- 1D0/24D0
      ELSE
       A1= 0.500D0
       A2= 0.125D0
       A3=-0.125D0
       A4= 0.000D0
C       A5= 0.500D0
       A6= 0.000D0
C       A7= 0.000D0
      ENDIF
C
C
      IF ((IINT.EQ.2).OR.(IINT.EQ.-2)) THEN
C
C *** linear Restriction of pressure
C
      DO 111 IEL2=1,NEL2
      VPIEL=REAL(DP2(IEL2))
      VAVOL=AVOL3(IEL2)
C
      IV1=KVERT2(1,IEL2)
      IV2=KVERT2(2,IEL2)
      IV3=KVERT2(3,IEL2)
      IV4=KVERT2(4,IEL2)
      IV5=KVERT2(5,IEL2)
      IV6=KVERT2(6,IEL2)
      IV7=KVERT2(7,IEL2)
      IV8=KVERT2(8,IEL2)
C
      VPL2(IV1)=VPL2(IV1)+0.125E0*VAVOL*VPIEL
      VPL2(IV2)=VPL2(IV2)+0.125E0*VAVOL*VPIEL
      VPL2(IV3)=VPL2(IV3)+0.125E0*VAVOL*VPIEL
      VPL2(IV4)=VPL2(IV4)+0.125E0*VAVOL*VPIEL
      VPL2(IV5)=VPL2(IV5)+0.125E0*VAVOL*VPIEL
      VPL2(IV6)=VPL2(IV6)+0.125E0*VAVOL*VPIEL
      VPL2(IV7)=VPL2(IV7)+0.125E0*VAVOL*VPIEL
      VPL2(IV8)=VPL2(IV8)+0.125E0*VAVOL*VPIEL
C
      VAUX(IV1)=VAUX(IV1)+0.125E0*VAVOL
      VAUX(IV2)=VAUX(IV2)+0.125E0*VAVOL
      VAUX(IV3)=VAUX(IV3)+0.125E0*VAVOL
      VAUX(IV4)=VAUX(IV4)+0.125E0*VAVOL
      VAUX(IV5)=VAUX(IV5)+0.125E0*VAVOL
      VAUX(IV6)=VAUX(IV6)+0.125E0*VAVOL
      VAUX(IV7)=VAUX(IV7)+0.125E0*VAVOL
      VAUX(IV8)=VAUX(IV8)+0.125E0*VAVOL
111   CONTINUE
C
      DO 13 IVT2=1,NVT2
13    VPL2(IVT2)=VPL2(IVT2)/VAUX(IVT2)
C
      CALL LCP2(VPL2,VPL1,NVT1)
C
      DO 112 IEL2=1,NEL2
C
      I1=KVERT2(1,IEL2)
      I3=KVERT2(3,IEL2)
      I6=KVERT2(6,IEL2)
      I7=KVERT2(7,IEL2)
      I8=KVERT2(8,IEL2)
C
      IF (KADJ2(1,IEL2).EQ.0) THEN
       AA3=0.25D0
      ELSE
       AA3=0.125D0
      ENDIF
C
      IF (KADJ2(2,IEL2).EQ.0) THEN
       AA6=0.25D0
      ELSE
       AA6=0.125D0
      ENDIF
C
      IF (KADJ2(5,IEL2).EQ.0) THEN
       AA8=0.25D0
      ELSE
       AA8=0.125D0
      ENDIF
C
      VPL1(I1)=VPL1(I1)+AA3*VPL2(I3)+AA6*VPL2(I6)+0.125D0*VPL2(I7)+
     *         AA8*VPL2(I8)
C
112   CONTINUE
C
      IVT=NVT1
C
      DO 113 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      I1 =KVERT2(2,IELH1)
      I2 =KVERT2(2,IELH2)
      I3 =KVERT2(2,IELH3)
      I4 =KVERT2(2,IELH4)
      I5 =KVERT2(5,IELH1)
      I6 =KVERT2(5,IELH2)
      I7 =KVERT2(5,IELH3)
      I8 =KVERT2(5,IELH4)
      I9 =KVERT2(2,IELH5)
      I10=KVERT2(2,IELH6)
      I11=KVERT2(2,IELH7)
      I12=KVERT2(2,IELH8)
C
      J1=KVERT1(1,IEL1)
      J2=KVERT1(2,IEL1)
      J3=KVERT1(3,IEL1)
      J4=KVERT1(4,IEL1)
      J5=KVERT1(5,IEL1)
      J6=KVERT1(6,IEL1)
      J7=KVERT1(7,IEL1)
      J8=KVERT1(8,IEL1)
C
      IF (I1.GT.IVT) THEN
       VPL1(J1)=VPL1(J1)+0.5D0*VPL2(I1)
       VPL1(J2)=VPL1(J2)+0.5D0*VPL2(I1)
       IVT=IVT+1
      ENDIF
C
      IF (I2.GT.IVT) THEN
       VPL1(J2)=VPL1(J2)+0.5D0*VPL2(I2)
       VPL1(J3)=VPL1(J3)+0.5D0*VPL2(I2)
       IVT=IVT+1
      ENDIF
C
      IF (I3.GT.IVT) THEN
       VPL1(J3)=VPL1(J3)+0.5D0*VPL2(I3)
       VPL1(J4)=VPL1(J4)+0.5D0*VPL2(I3)
       IVT=IVT+1
      ENDIF
C
      IF (I4.GT.IVT) THEN
       VPL1(J1)=VPL1(J1)+0.5D0*VPL2(I4)
       VPL1(J4)=VPL1(J4)+0.5D0*VPL2(I4)
       IVT=IVT+1
      ENDIF
C
      IF (I5.GT.IVT) THEN
       VPL1(J1)=VPL1(J1)+0.5D0*VPL2(I5)
       VPL1(J5)=VPL1(J5)+0.5D0*VPL2(I5)
       IVT=IVT+1
      ENDIF
C
      IF (I6.GT.IVT) THEN
       VPL1(J2)=VPL1(J2)+0.5D0*VPL2(I6)
       VPL1(J6)=VPL1(J6)+0.5D0*VPL2(I6)
       IVT=IVT+1
      ENDIF
C
      IF (I7.GT.IVT) THEN
       VPL1(J3)=VPL1(J3)+0.5D0*VPL2(I7)
       VPL1(J7)=VPL1(J7)+0.5D0*VPL2(I7)
       IVT=IVT+1
      ENDIF
C
      IF (I8.GT.IVT) THEN
       VPL1(J4)=VPL1(J4)+0.5D0*VPL2(I8)
       VPL1(J8)=VPL1(J8)+0.5D0*VPL2(I8)
       IVT=IVT+1
      ENDIF
C
      IF (I9.GT.IVT) THEN
       VPL1(J5)=VPL1(J5)+0.5D0*VPL2(I9)
       VPL1(J6)=VPL1(J6)+0.5D0*VPL2(I9)
       IVT=IVT+1
      ENDIF
C
      IF (I10.GT.IVT) THEN
       VPL1(J6)=VPL1(J6)+0.5D0*VPL2(I10)
       VPL1(J7)=VPL1(J7)+0.5D0*VPL2(I10)
       IVT=IVT+1
      ENDIF
C
      IF (I11.GT.IVT) THEN
       VPL1(J7)=VPL1(J7)+0.5D0*VPL2(I11)
       VPL1(J8)=VPL1(J8)+0.5D0*VPL2(I11)
       IVT=IVT+1
      ENDIF
C
      IF (I12.GT.IVT) THEN
       VPL1(J5)=VPL1(J5)+0.5D0*VPL2(I12)
       VPL1(J8)=VPL1(J8)+0.5D0*VPL2(I12)
       IVT=IVT+1
      ENDIF
C
113    CONTINUE
C
C
      ENDIF
C
      CALL LCL1(DU1,NAT1)
      CALL LCL1(DV1,NAT1)
      CALL LCL1(DW1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      V1=VPL1(KVERT1(1,IEL1))
      V2=VPL1(KVERT1(2,IEL1))
      V3=VPL1(KVERT1(3,IEL1))
      V4=VPL1(KVERT1(4,IEL1))
      V5=VPL1(KVERT1(5,IEL1))
      V6=VPL1(KVERT1(6,IEL1))
      V7=VPL1(KVERT1(7,IEL1))
      V8=VPL1(KVERT1(8,IEL1))
C
C
C *** Restriction of pressure
C
      IF ((IINT.EQ.2).OR.(IINT.EQ.-2)) THEN
C *** linear restriction of pressure
       DP1(IEL1)=0.125D0*(V1+V2+V3+V4+V5+V6+V7+V8)
      ENDIF
C
      IF ((IINT.EQ.1).OR.(IINT.EQ.-1)) THEN
C *** constant restriction
       DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)+
     *           DP2(IELH5)+DP2(IELH6)+DP2(IELH7)+DP2(IELH8)
      ENDIF
C
C
C
      DUH1= DU2(KAREA2(1,IELH1))
      DUH2= DU2(KAREA2(1,IELH2))
      DUH3= DU2(KAREA2(1,IELH3))
      DUH4= DU2(KAREA2(1,IELH4))
      DUH5= DU2(KAREA2(2,IELH1))
      DUH6= DU2(KAREA2(5,IELH2))
      DUH7= DU2(KAREA2(5,IELH6))
      DUH8= DU2(KAREA2(2,IELH5))
      DUH9= DU2(KAREA2(2,IELH2))
      DUH10=DU2(KAREA2(5,IELH3))
      DUH11=DU2(KAREA2(5,IELH7))
      DUH12=DU2(KAREA2(2,IELH6))
      DUH13=DU2(KAREA2(2,IELH3))
      DUH14=DU2(KAREA2(5,IELH4))
      DUH15=DU2(KAREA2(5,IELH8))
      DUH16=DU2(KAREA2(2,IELH7))
      DUH17=DU2(KAREA2(2,IELH4))
      DUH18=DU2(KAREA2(5,IELH1))
      DUH19=DU2(KAREA2(5,IELH5))
      DUH20=DU2(KAREA2(2,IELH8))
      DUH21=DU2(KAREA2(1,IELH5))
      DUH22=DU2(KAREA2(1,IELH6))
      DUH23=DU2(KAREA2(1,IELH7))
      DUH24=DU2(KAREA2(1,IELH8))
      DUH25=DU2(KAREA2(3,IELH1))
      DUH26=DU2(KAREA2(3,IELH2))
      DUH27=DU2(KAREA2(3,IELH3))
      DUH28=DU2(KAREA2(3,IELH4))
      DUH29=DU2(KAREA2(6,IELH1))
      DUH30=DU2(KAREA2(6,IELH2))
      DUH31=DU2(KAREA2(6,IELH3))
      DUH32=DU2(KAREA2(6,IELH4))
      DUH33=DU2(KAREA2(3,IELH5))
      DUH34=DU2(KAREA2(3,IELH6))
      DUH35=DU2(KAREA2(3,IELH7))
      DUH36=DU2(KAREA2(3,IELH8))
C
C
      B1=DUH1+DUH2+DUH3+DUH4
      B2=DUH5+DUH6+DUH9+DUH10+DUH13+DUH14+DUH17+DUH18
      B3=DUH7+DUH8+DUH11+DUH12+DUH15+DUH16+DUH19+DUH20
      B4=DUH21+DUH22+DUH23+DUH24
      B5=DUH25+DUH26+DUH27+DUH28
      B6=DUH29+DUH30+DUH31+DUH32
      B7=DUH33+DUH34+DUH35+DUH36
      B8=DUH5+DUH6+DUH7+DUH8
      B9=DUH1+DUH2+DUH9+DUH12+DUH21+DUH22+DUH18+DUH19
      B10=DUH3+DUH4+DUH10+DUH11+DUH23+DUH24+DUH17+DUH20
      B11=DUH13+DUH14+DUH15+DUH16
      B12=DUH25+DUH29+DUH30+DUH33
      B13=DUH26+DUH28+DUH34+DUH36
      B14=DUH27+DUH31+DUH32+DUH35     
      B15=DUH9+DUH10+DUH11+DUH12
      B16=DUH2+DUH3+DUH6+DUH7+DUH22+DUH23+DUH13+DUH16
      B17=DUH1+DUH4+DUH5+DUH8+DUH21+DUH24+DUH14+DUH15
      B18=DUH17+DUH18+DUH19+DUH20
      B19=DUH26+DUH30+DUH31+DUH34
      B20=DUH25+DUH27+DUH33+DUH35
      B21=DUH28+DUH29+DUH32+DUH36
C
      AVOL11=AVOL(IEL1)
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DU1(IA1)= DU1(IA1)+(A1*(B1+B5)+A2*B2+A3*B3+
     *           A4*(B4+B7)+A6*B6)*AVOL11
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DU1(IA2)= DU1(IA2)+(A1*(B8+B12)+A2*B9+A3*B10+
     *           A4*(B11+B14)+A6*B13)*AVOL11
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DU1(IA3)= DU1(IA3)+(A1*(B15+B19)+A2*B16+A3*B17+
     *           A4*(B18+B21)+A6*B20)*AVOL11
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DU1(IA4)= DU1(IA4)+(A1*(B11+B14)+A2*B10+A3*B9+
     *           A4*(B8+B12)+A6*B13)*AVOL11
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DU1(IA5)= DU1(IA5)+(A1*(B18+B21)+A2*B17+A3*B16+
     *           A4*(B15+B19)+A6*B20)*AVOL11
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DU1(IA6)= DU1(IA6)+(A1*(B4+B7)+A2*B3+A3*B2+
     *           A4*(B1+B5)+A6*B6)*AVOL11
      ENDIF
C
C
C
      DVH1= DV2(KAREA2(1,IELH1))
      DVH2= DV2(KAREA2(1,IELH2))
      DVH3= DV2(KAREA2(1,IELH3))
      DVH4= DV2(KAREA2(1,IELH4))
      DVH5= DV2(KAREA2(2,IELH1))
      DVH6= DV2(KAREA2(5,IELH2))
      DVH7= DV2(KAREA2(5,IELH6))
      DVH8= DV2(KAREA2(2,IELH5))
      DVH9= DV2(KAREA2(2,IELH2))
      DVH10=DV2(KAREA2(5,IELH3))
      DVH11=DV2(KAREA2(5,IELH7))
      DVH12=DV2(KAREA2(2,IELH6))
      DVH13=DV2(KAREA2(2,IELH3))
      DVH14=DV2(KAREA2(5,IELH4))
      DVH15=DV2(KAREA2(5,IELH8))
      DVH16=DV2(KAREA2(2,IELH7))
      DVH17=DV2(KAREA2(2,IELH4))
      DVH18=DV2(KAREA2(5,IELH1))
      DVH19=DV2(KAREA2(5,IELH5))
      DVH20=DV2(KAREA2(2,IELH8))
      DVH21=DV2(KAREA2(1,IELH5))
      DVH22=DV2(KAREA2(1,IELH6))
      DVH23=DV2(KAREA2(1,IELH7))
      DVH24=DV2(KAREA2(1,IELH8))
      DVH25=DV2(KAREA2(3,IELH1))
      DVH26=DV2(KAREA2(3,IELH2))
      DVH27=DV2(KAREA2(3,IELH3))
      DVH28=DV2(KAREA2(3,IELH4))
      DVH29=DV2(KAREA2(6,IELH1))
      DVH30=DV2(KAREA2(6,IELH2))
      DVH31=DV2(KAREA2(6,IELH3))
      DVH32=DV2(KAREA2(6,IELH4))
      DVH33=DV2(KAREA2(3,IELH5))
      DVH34=DV2(KAREA2(3,IELH6))
      DVH35=DV2(KAREA2(3,IELH7))
      DVH36=DV2(KAREA2(3,IELH8))
C
      C1=DVH1+DVH2+DVH3+DVH4
      C2=DVH5+DVH6+DVH9+DVH10+DVH13+DVH14+DVH17+DVH18
      C3=DVH7+DVH8+DVH11+DVH12+DVH15+DVH16+DVH19+DVH20
      C4=DVH21+DVH22+DVH23+DVH24
      C5=DVH25+DVH26+DVH27+DVH28
      C6=DVH29+DVH30+DVH31+DVH32
      C7=DVH33+DVH34+DVH35+DVH36
      C8=DVH5+DVH6+DVH7+DVH8
      C9=DVH1+DVH2+DVH9+DVH12+DVH21+DVH22+DVH18+DVH19
      C10=DVH3+DVH4+DVH10+DVH11+DVH23+DVH24+DVH17+DVH20
      C11=DVH13+DVH14+DVH15+DVH16
      C12=DVH25+DVH29+DVH30+DVH33
      C13=DVH26+DVH28+DVH34+DVH36
      C14=DVH27+DVH31+DVH32+DVH35     
      C15=DVH9+DVH10+DVH11+DVH12
      C16=DVH2+DVH3+DVH6+DVH7+DVH22+DVH23+DVH13+DVH16
      C17=DVH1+DVH4+DVH5+DVH8+DVH21+DVH24+DVH14+DVH15
      C18=DVH17+DVH18+DVH19+DVH20
      C19=DVH26+DVH30+DVH31+DVH34
      C20=DVH25+DVH27+DVH33+DVH35
      C21=DVH28+DVH29+DVH32+DVH36     
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DV1(IA1)= DV1(IA1)+(A1*(C1+C5)+A2*C2+A3*C3+
     *           A4*(C4+C7)+A6*C6)*AVOL11
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DV1(IA2)= DV1(IA2)+(A1*(C8+C12)+A2*C9+A3*C10+
     *           A4*(C11+C14)+A6*C13)*AVOL11
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DV1(IA3)= DV1(IA3)+(A1*(C15+C19)+A2*C16+A3*C17+
     *           A4*(C18+C21)+A6*C20)*AVOL11
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DV1(IA4)= DV1(IA4)+(A1*(C11+C14)+A2*C10+A3*C9+
     *           A4*(C8+C12)+A6*C13)*AVOL11
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DV1(IA5)= DV1(IA5)+(A1*(C18+C21)+A2*C17+A3*C16+
     *           A4*(C15+C19)+A6*C20)*AVOL11
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DV1(IA6)= DV1(IA6)+(A1*(C4+C7)+A2*C3+A3*C2+
     *           A4*(C1+C5)+A6*C6)*AVOL11
      ENDIF
C
C
C
      DWH1= DW2(KAREA2(1,IELH1))
      DWH2= DW2(KAREA2(1,IELH2))
      DWH3= DW2(KAREA2(1,IELH3))
      DWH4= DW2(KAREA2(1,IELH4))
      DWH5= DW2(KAREA2(2,IELH1))
      DWH6= DW2(KAREA2(5,IELH2))
      DWH7= DW2(KAREA2(5,IELH6))
      DWH8= DW2(KAREA2(2,IELH5))
      DWH9= DW2(KAREA2(2,IELH2))
      DWH10=DW2(KAREA2(5,IELH3))
      DWH11=DW2(KAREA2(5,IELH7))
      DWH12=DW2(KAREA2(2,IELH6))
      DWH13=DW2(KAREA2(2,IELH3))
      DWH14=DW2(KAREA2(5,IELH4))
      DWH15=DW2(KAREA2(5,IELH8))
      DWH16=DW2(KAREA2(2,IELH7))
      DWH17=DW2(KAREA2(2,IELH4))
      DWH18=DW2(KAREA2(5,IELH1))
      DWH19=DW2(KAREA2(5,IELH5))
      DWH20=DW2(KAREA2(2,IELH8))
      DWH21=DW2(KAREA2(1,IELH5))
      DWH22=DW2(KAREA2(1,IELH6))
      DWH23=DW2(KAREA2(1,IELH7))
      DWH24=DW2(KAREA2(1,IELH8))
      DWH25=DW2(KAREA2(3,IELH1))
      DWH26=DW2(KAREA2(3,IELH2))
      DWH27=DW2(KAREA2(3,IELH3))
      DWH28=DW2(KAREA2(3,IELH4))
      DWH29=DW2(KAREA2(6,IELH1))
      DWH30=DW2(KAREA2(6,IELH2))
      DWH31=DW2(KAREA2(6,IELH3))
      DWH32=DW2(KAREA2(6,IELH4))
      DWH33=DW2(KAREA2(3,IELH5))
      DWH34=DW2(KAREA2(3,IELH6))
      DWH35=DW2(KAREA2(3,IELH7))
      DWH36=DW2(KAREA2(3,IELH8))
C
      D1=DWH1+DWH2+DWH3+DWH4
      D2=DWH5+DWH6+DWH9+DWH10+DWH13+DWH14+DWH17+DWH18
      D3=DWH7+DWH8+DWH11+DWH12+DWH15+DWH16+DWH19+DWH20
      D4=DWH21+DWH22+DWH23+DWH24
      D5=DWH25+DWH26+DWH27+DWH28
      D6=DWH29+DWH30+DWH31+DWH32
      D7=DWH33+DWH34+DWH35+DWH36
      D8=DWH5+DWH6+DWH7+DWH8
      D9=DWH1+DWH2+DWH9+DWH12+DWH21+DWH22+DWH18+DWH19
      D10=DWH3+DWH4+DWH10+DWH11+DWH23+DWH24+DWH17+DWH20
      D11=DWH13+DWH14+DWH15+DWH16
      D12=DWH25+DWH29+DWH30+DWH33
      D13=DWH26+DWH28+DWH34+DWH36
      D14=DWH27+DWH31+DWH32+DWH35     
      D15=DWH9+DWH10+DWH11+DWH12
      D16=DWH2+DWH3+DWH6+DWH7+DWH22+DWH23+DWH13+DWH16
      D17=DWH1+DWH4+DWH5+DWH8+DWH21+DWH24+DWH14+DWH15
      D18=DWH17+DWH18+DWH19+DWH20
      D19=DWH26+DWH30+DWH31+DWH34
      D20=DWH25+DWH27+DWH33+DWH35
      D21=DWH28+DWH29+DWH32+DWH36     
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DW1(IA1)= DW1(IA1)+(A1*(D1+D5)+A2*D2+A3*D3+
     *           A4*(D4+D7)+A6*D6)*AVOL11
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DW1(IA2)= DW1(IA2)+(A1*(D8+D12)+A2*D9+A3*D10+
     *           A4*(D11+D14)+A6*D13)*AVOL11
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DW1(IA3)= DW1(IA3)+(A1*(D15+D19)+A2*D16+A3*D17+
     *           A4*(D18+D21)+A6*D20)*AVOL11
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DW1(IA4)= DW1(IA4)+(A1*(D11+D14)+A2*D10+A3*D9+
     *           A4*(D8+D12)+A6*D13)*AVOL11
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DW1(IA5)= DW1(IA5)+(A1*(D18+D21)+A2*D17+A3*D16+
     *           A4*(D15+D19)+A6*D20)*AVOL11
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DW1(IA6)= DW1(IA6)+(A1*(D4+D7)+A2*D3+A3*D2+
     *           A4*(D1+D5)+A6*D6)*AVOL11
      ENDIF
C
10    CONTINUE
C
      DO 555 IEL1=1,NEL1
      DO 666 I=1,6
      IADJ=KADJ1(I,IEL1)
      IA=KAREA1(I,IEL1)
      IF (IADJ.EQ.0) THEN
      GOTO 666
      ELSE
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DU1(IA)=(SQRT(2D0)/SQRT(AVOL1+AVOL2))*DU1(IA)
      DV1(IA)=(SQRT(2D0)/SQRT(AVOL1+AVOL2))*DV1(IA)
      DW1(IA)=(SQRT(2D0)/SQRT(AVOL1+AVOL2))*DW1(IA)
      ENDIF
666   CONTINUE
555   CONTINUE
      END
c
************************************************************************
      SUBROUTINE MP031(DX1,DX2,KAREA1,KAREA2,KADJ1,KADJ2,NAT1,NAT2,
     *                 NEL1,NEL2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 11D0/24D0,A2= 7D0/48D0,A3=-5D0/48D0,A4=-1D0/24D0)
      PARAMETER (A5= 11D0/24D0,A6= 1D0/12D0,A7=-1D0/24D0)
      DIMENSION DX1(1),DX2(1),KAREA1(NNAE,1),KAREA2(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DX2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DY1=DX1(KAREA1(1,IEL1))
      DY2=DX1(KAREA1(2,IEL1))
      DY3=DX1(KAREA1(3,IEL1))
      DY4=DX1(KAREA1(4,IEL1))
      DY5=DX1(KAREA1(5,IEL1))
      DY6=DX1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+
     *                      A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A4*DY6
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+
     *                      A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A4*DY6
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+
     *                      A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A4*DY6
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+
     *                      A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A4*DY6
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+
     *                      A2*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+
     *                      A2*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+
     *                      A3*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A2*DY6
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+
     *                      A3*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+
     *                      A2*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A3*DY6
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+
     *                      A2*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A3*DY6
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+
     *                      A3*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A2*DY6
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+
     *                      A3*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+
     *                      A2*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+
     *                      A2*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+
     *                      A3*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A2*DY6
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+
     *                      A3*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+
     *                      A2*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+
     *                      A2*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+
     *                      A3*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A2*DY6
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+
     *                      A3*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+
     *                      A4*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+
     *                      A4*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+
     *                      A4*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+
     *                      A4*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6
      ENDIF
C
C
      DX2(KAREA2(3,IELH1))=A5*DY1+A5*DY2+A6*DY3+A7*DY4+A6*DY5+A7*DY6
      DX2(KAREA2(3,IELH2))=A5*DY1+A6*DY2+A5*DY3+A6*DY4+A7*DY5+A7*DY6
      DX2(KAREA2(3,IELH3))=A5*DY1+A7*DY2+A6*DY3+A5*DY4+A6*DY5+A7*DY6
      DX2(KAREA2(3,IELH4))=A5*DY1+A6*DY2+A7*DY3+A6*DY4+A5*DY5+A7*DY6
C
      DX2(KAREA2(6,IELH1))=A6*DY1+A5*DY2+A7*DY3+A7*DY4+A5*DY5+A6*DY6
      DX2(KAREA2(6,IELH2))=A6*DY1+A5*DY2+A5*DY3+A7*DY4+A7*DY5+A6*DY6
      DX2(KAREA2(6,IELH3))=A6*DY1+A7*DY2+A5*DY3+A5*DY4+A7*DY5+A6*DY6
      DX2(KAREA2(6,IELH4))=A6*DY1+A7*DY2+A7*DY3+A5*DY4+A5*DY5+A6*DY6
C
      DX2(KAREA2(3,IELH5))=A7*DY1+A5*DY2+A6*DY3+A7*DY4+A6*DY5+A5*DY6
      DX2(KAREA2(3,IELH6))=A7*DY1+A6*DY2+A5*DY3+A6*DY4+A7*DY5+A5*DY6
      DX2(KAREA2(3,IELH7))=A7*DY1+A7*DY2+A6*DY3+A5*DY4+A6*DY5+A5*DY6
      DX2(KAREA2(3,IELH8))=A7*DY1+A6*DY2+A7*DY3+A6*DY4+A5*DY5+A5*DY6
C
      IF (KADJ1(1,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+2D0*
     *                   (A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A4*DY6)
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+2D0*
     *                   (A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A4*DY6)
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+2D0*
     *                   (A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A4*DY6)
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+2D0*
     *                   (A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A4*DY6)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+2D0*
     *                   (A2*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+2D0*
     *                   (A2*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+2D0*
     *                   (A3*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A2*DY6)
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+2D0*
     *                   (A3*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+2D0*
     *                   (A2*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A3*DY6)
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+2D0*
     *                   (A2*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A3*DY6)
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+2D0*
     *                   (A3*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A2*DY6)
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+2D0*
     *                   (A3*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+2D0*
     *                   (A2*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+2D0*
     *                   (A2*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+2D0*
     *                   (A3*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A2*DY6)
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+2D0*
     *                   (A3*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+2D0*
     *                   (A2*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+2D0*
     *                   (A2*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+2D0*
     *                   (A3*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A2*DY6)
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+2D0*
     *                   (A3*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+2D0*
     *                   (A4*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6)
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+2D0*
     *                   (A4*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+2D0*
     *                   (A4*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+2D0*
     *                   (A4*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6)
      ENDIF
C
c
10    CONTINUE
C
C
      END
C
************************************************************************
      SUBROUTINE MP030(DX1,DX2,KAREA1,KAREA2,KADJ1,KADJ2,NAT1,NAT2,
     *                 NEL1,NEL2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 0.5D0,A2= 0.125D0,A3=-0.125D0,A4= 0D0)
      PARAMETER (A5=0.5D0,A6=0D0,A7=0D0)
      DIMENSION DX1(1),DX2(1),KAREA1(NNAE,1),KAREA2(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DX2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DY1=DX1(KAREA1(1,IEL1))
      DY2=DX1(KAREA1(2,IEL1))
      DY3=DX1(KAREA1(3,IEL1))
      DY4=DX1(KAREA1(4,IEL1))
      DY5=DX1(KAREA1(5,IEL1))
      DY6=DX1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+
     *                      A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+
     *                      A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+
     *                      A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+
     *                      A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+
     *                      A2*DY1+A1*DY2+A3*DY3+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+
     *                      A2*DY1+A1*DY2+A2*DY3+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+
     *                      A3*DY1+A1*DY2+A2*DY3+A3*DY5+A2*DY6
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+
     *                      A3*DY1+A1*DY2+A3*DY3+A2*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+
     *                      A2*DY1+A2*DY2+A1*DY3+A3*DY4+A3*DY6
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+
     *                      A2*DY1+A3*DY2+A1*DY3+A2*DY4+A3*DY6
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+
     *                      A3*DY1+A3*DY2+A1*DY3+A2*DY4+A2*DY6
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+
     *                      A3*DY1+A2*DY2+A1*DY3+A3*DY4+A2*DY6
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+
     *                      A2*DY1+A2*DY3+A1*DY4+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+
     *                      A2*DY1+A3*DY3+A1*DY4+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+
     *                      A3*DY1+A3*DY3+A1*DY4+A2*DY5+A2*DY6
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+
     *                      A3*DY1+A2*DY3+A1*DY4+A3*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+
     *                      A2*DY1+A3*DY2+A2*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+
     *                      A2*DY1+A2*DY2+A3*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+
     *                      A3*DY1+A2*DY2+A3*DY4+A1*DY5+A2*DY6
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+
     *                      A3*DY1+A3*DY2+A2*DY4+A1*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+
     *                      A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+
     *                      A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+
     *                      A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+
     *                      A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6
      ENDIF
C
C
      DX2(KAREA2(3,IELH1))=A5*(DY1+DY2)
      DX2(KAREA2(3,IELH2))=A5*(DY1+DY3)
      DX2(KAREA2(3,IELH3))=A5*(DY1+DY4)
      DX2(KAREA2(3,IELH4))=A5*(DY1+DY5)
C
      DX2(KAREA2(6,IELH1))=A5*(DY2+DY5)
      DX2(KAREA2(6,IELH2))=A5*(DY2+DY3)
      DX2(KAREA2(6,IELH3))=A5*(DY3+DY4)
      DX2(KAREA2(6,IELH4))=A5*(DY4+DY5)
C
      DX2(KAREA2(3,IELH5))=A5*(DY2+DY6)
      DX2(KAREA2(3,IELH6))=A5*(DY3+DY6)
      DX2(KAREA2(3,IELH7))=A5*(DY4+DY6)
      DX2(KAREA2(3,IELH8))=A5*(DY5+DY6)
C
      IF (KADJ1(1,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+2D0*
     *                   (A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5)
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+2D0*
     *                   (A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5)
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+2D0*
     *                   (A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5)
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+2D0*
     *                   (A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+2D0*
     *                   (A2*DY1+A1*DY2+A3*DY3+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+2D0*
     *                   (A2*DY1+A1*DY2+A2*DY3+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+2D0*
     *                   (A3*DY1+A1*DY2+A2*DY3+A3*DY5+A2*DY6)
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+2D0*
     *                   (A3*DY1+A1*DY2+A3*DY3+A2*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+2D0*
     *                   (A2*DY1+A2*DY2+A1*DY3+A3*DY4+A3*DY6)
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+2D0*
     *                   (A2*DY1+A3*DY2+A1*DY3+A2*DY4+A3*DY6)
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+2D0*
     *                   (A3*DY1+A3*DY2+A1*DY3+A2*DY4+A2*DY6)
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+2D0*
     *                   (A3*DY1+A2*DY2+A1*DY3+A3*DY4+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+2D0*
     *                   (A2*DY1+A2*DY3+A1*DY4+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+2D0*
     *                   (A2*DY1+A3*DY3+A1*DY4+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+2D0*
     *                   (A3*DY1+A3*DY3+A1*DY4+A2*DY5+A2*DY6)
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+2D0*
     *                   (A3*DY1+A2*DY3+A1*DY4+A3*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+2D0*
     *                   (A2*DY1+A3*DY2+A2*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+2D0*
     *                   (A2*DY1+A2*DY2+A3*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+2D0*
     *                   (A3*DY1+A2*DY2+A3*DY4+A1*DY5+A2*DY6)
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+2D0*
     *                   (A3*DY1+A3*DY2+A2*DY4+A1*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+2D0*
     *                   (A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6)
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+2D0*
     *                   (A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+2D0*
     *                   (A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+2D0*
     *                   (A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6)
      ENDIF
C
C
10    CONTINUE
C
C
      END
C
C
************************************************************************
      SUBROUTINE MR031(DF2,DF1,KAREA2,KAREA1,KADJ2,KADJ1,NAT2,NAT1,
     *                 NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 11D0/24D0,A2= 7D0/48D0,A3=-5D0/48D0,A4=-1D0/24D0)
      PARAMETER (A5= 11D0/24D0,A6= 1D0/12D0,A7=-1D0/24D0)
      DIMENSION DF1(1),DF2(1),KAREA2(NNAE,1),KAREA1(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DF1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      DY1= DF2(KAREA2(1,IELH1))
      DY2= DF2(KAREA2(1,IELH2))
      DY3= DF2(KAREA2(1,IELH3))
      DY4= DF2(KAREA2(1,IELH4))
      DY5= DF2(KAREA2(2,IELH1))
      DY6= DF2(KAREA2(5,IELH2))
      DY7= DF2(KAREA2(5,IELH6))
      DY8= DF2(KAREA2(2,IELH5))
      DY9= DF2(KAREA2(2,IELH2))
      DY10=DF2(KAREA2(5,IELH3))
      DY11=DF2(KAREA2(5,IELH7))
      DY12=DF2(KAREA2(2,IELH6))
      DY13=DF2(KAREA2(2,IELH3))
      DY14=DF2(KAREA2(5,IELH4))
      DY15=DF2(KAREA2(5,IELH8))
      DY16=DF2(KAREA2(2,IELH7))
      DY17=DF2(KAREA2(2,IELH4))
      DY18=DF2(KAREA2(5,IELH1))
      DY19=DF2(KAREA2(5,IELH5))
      DY20=DF2(KAREA2(2,IELH8))
      DY21=DF2(KAREA2(1,IELH5))
      DY22=DF2(KAREA2(1,IELH6))
      DY23=DF2(KAREA2(1,IELH7))
      DY24=DF2(KAREA2(1,IELH8))
      DY25=DF2(KAREA2(3,IELH1))
      DY26=DF2(KAREA2(3,IELH2))
      DY27=DF2(KAREA2(3,IELH3))
      DY28=DF2(KAREA2(3,IELH4))
      DY29=DF2(KAREA2(6,IELH1))
      DY30=DF2(KAREA2(6,IELH2))
      DY31=DF2(KAREA2(6,IELH3))
      DY32=DF2(KAREA2(6,IELH4))
      DY33=DF2(KAREA2(3,IELH5))
      DY34=DF2(KAREA2(3,IELH6))
      DY35=DF2(KAREA2(3,IELH7))
      DY36=DF2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DF1(IA1)= DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A4*(DY21+DY22+DY23+DY24)
     *          +A5*(DY25+DY26+DY27+DY28)
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY33+DY34+DY35+DY36)     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DF1(IA2)= DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A4*(DY13+DY14+DY15+DY16)
     *          +A5*(DY25+DY29+DY30+DY33)
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY27+DY31+DY32+DY35)     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DF1(IA3)= DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A4*(DY17+DY18+DY19+DY20)
     *          +A5*(DY26+DY30+DY31+DY34)
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY28+DY29+DY32+DY36)     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DF1(IA4)= DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A4*(DY5+DY6+DY7+DY8)
     *          +A5*(DY27+DY31+DY32+DY35)  
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY25+DY29+DY30+DY33)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DF1(IA5)= DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A4*(DY9+DY10+DY11+DY12)
     *          +A5*(DY28+DY29+DY32+DY36)     
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY26+DY30+DY31+DY34)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DF1(IA6)= DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A4*(DY1+DY2+DY3+DY4)
     *          +A5*(DY33+DY34+DY35+DY36)     
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY25+DY26+DY27+DY28)
      ENDIF
C
      IF (KADJ1(1,IEL1).EQ.0) THEN
       DF1(IA1)= 2D0*(DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A4*(DY21+DY22+DY23+DY24)
     *          +A5*(DY25+DY26+DY27+DY28)
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY33+DY34+DY35+DY36))     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN
       DF1(IA2)= 2D0*(DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A4*(DY13+DY14+DY15+DY16)
     *          +A5*(DY25+DY29+DY30+DY33)
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY27+DY31+DY32+DY35))     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN
       DF1(IA3)= 2D0*(DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A4*(DY17+DY18+DY19+DY20)
     *          +A5*(DY26+DY30+DY31+DY34)
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY28+DY29+DY32+DY36))     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN
       DF1(IA4)= 2D0*(DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A4*(DY5+DY6+DY7+DY8)
     *          +A5*(DY27+DY31+DY32+DY35)  
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY25+DY29+DY30+DY33))
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN
       DF1(IA5)= 2D0*(DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A4*(DY9+DY10+DY11+DY12)
     *          +A5*(DY28+DY29+DY32+DY36)     
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY26+DY30+DY31+DY34))
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN
       DF1(IA6)= 2D0*(DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A4*(DY1+DY2+DY3+DY4)
     *          +A5*(DY33+DY34+DY35+DY36)     
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY25+DY26+DY27+DY28))
      ENDIF
C
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE MR030(DF2,DF1,KAREA2,KAREA1,KADJ2,KADJ1,NAT2,NAT1,
     *                 NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 0.5D0,A2= 0.125D0,A3=-0.125D0,A4= 0D0)
      PARAMETER (A5=0.5D0,A6=0D0,A7=0D0)
      DIMENSION DF1(1),DF2(1),KAREA2(NNAE,1),KAREA1(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DF1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      DY1= DF2(KAREA2(1,IELH1))
      DY2= DF2(KAREA2(1,IELH2))
      DY3= DF2(KAREA2(1,IELH3))
      DY4= DF2(KAREA2(1,IELH4))
      DY5= DF2(KAREA2(2,IELH1))
      DY6= DF2(KAREA2(5,IELH2))
      DY7= DF2(KAREA2(5,IELH6))
      DY8= DF2(KAREA2(2,IELH5))
      DY9= DF2(KAREA2(2,IELH2))
      DY10=DF2(KAREA2(5,IELH3))
      DY11=DF2(KAREA2(5,IELH7))
      DY12=DF2(KAREA2(2,IELH6))
      DY13=DF2(KAREA2(2,IELH3))
      DY14=DF2(KAREA2(5,IELH4))
      DY15=DF2(KAREA2(5,IELH8))
      DY16=DF2(KAREA2(2,IELH7))
      DY17=DF2(KAREA2(2,IELH4))
      DY18=DF2(KAREA2(5,IELH1))
      DY19=DF2(KAREA2(5,IELH5))
      DY20=DF2(KAREA2(2,IELH8))
      DY21=DF2(KAREA2(1,IELH5))
      DY22=DF2(KAREA2(1,IELH6))
      DY23=DF2(KAREA2(1,IELH7))
      DY24=DF2(KAREA2(1,IELH8))
      DY25=DF2(KAREA2(3,IELH1))
      DY26=DF2(KAREA2(3,IELH2))
      DY27=DF2(KAREA2(3,IELH3))
      DY28=DF2(KAREA2(3,IELH4))
      DY29=DF2(KAREA2(6,IELH1))
      DY30=DF2(KAREA2(6,IELH2))
      DY31=DF2(KAREA2(6,IELH3))
      DY32=DF2(KAREA2(6,IELH4))
      DY33=DF2(KAREA2(3,IELH5))
      DY34=DF2(KAREA2(3,IELH6))
      DY35=DF2(KAREA2(3,IELH7))
      DY36=DF2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DF1(IA1)= DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A5*(DY25+DY26+DY27+DY28)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DF1(IA2)= DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A5*(DY25+DY29+DY30+DY33)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DF1(IA3)= DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A5*(DY26+DY30+DY31+DY34)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DF1(IA4)= DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A5*(DY27+DY31+DY32+DY35)  
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DF1(IA5)= DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A5*(DY28+DY29+DY32+DY36)     
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DF1(IA6)= DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A5*(DY33+DY34+DY35+DY36)     
      ENDIF
C
      IF (KADJ1(1,IEL1).EQ.0) THEN
       DF1(IA1)= 2D0*(DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A5*(DY25+DY26+DY27+DY28))
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN
       DF1(IA2)= 2D0*(DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A5*(DY25+DY29+DY30+DY33))
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN
       DF1(IA3)= 2D0*(DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A5*(DY26+DY30+DY31+DY34))
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN
       DF1(IA4)= 2D0*(DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A5*(DY27+DY31+DY32+DY35))  
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN
       DF1(IA5)= 2D0*(DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A5*(DY28+DY29+DY32+DY36))     
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN
       DF1(IA6)= 2D0*(DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A5*(DY33+DY34+DY35+DY36))    
      ENDIF
C
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   YAX (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
*
*            - DX,DAX  have the structure  D=(D1,D2,D3,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DAX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C=======================================================================
C     Getting all parameters for MATML
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      I3=1+NU+NU
      IP=I3+NU
C
      KABD =L(KLABD(ILEV))
      NABD= KNABD(ILEV)
C
C=======================================================================
C
      CALL MATML(DAX(1),DAX(I2),DAX(I3),DAX(IP),
     *            DX(1),DX(I2),DX(I3),DX(IP),A1,A2,
     *            VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            VWORK(KB1),VWORK(KB2),VWORK(KB3),KWORK(KCOLB),
     *            KWORK(KLDB),NU,NP,KWORK(KABD),NABD,INEUM)
C
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   YDBC  (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX=(D1,D2,D3) to 0
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C=======================================================================
C     Getting all parameters for DX
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      I3=1+NU+NU
C
      KABD =L(KLABD(ILEV))
      NABD= KNABD(ILEV)
C
      CALL  BDRY0 (DX(1),DX(I2),DX(I3),KWORK(KABD),NABD)
C
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   YEX (DX,DB,DD,NEQ,RHO)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controlled by variables on
*              COMMON blocks 
*  
*            - DX,DB,DD have the structure  D=(D1,D2,D3,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      I3=1+NU+NU
      IP=I3+NU
C
      KVERT=L(LVERT)
      KAREA =L(LAREA)
      KNPR =L(LNPR )
      KABD =L(KLABD(ILEV))
      NABD= KNABD(ILEV)
C
      KVOL=L(KLVOL(ILEV))
C
C=======================================================================
C
      IF (ISL.EQ.1) THEN
C
       DO 11 ITE=1,NSL
C
       CALL  VANCAE (DX(1),DX(I2),DX(I3),DX(IP),DD(1),DD(I2),DD(I3),
     *               DD(IP),DB(1),DB(I2),DB(I3),DB(IP),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               VWORK(KB1),VWORK(KB2),VWORK(KB3),KWORK(KCOLB),
     *               KWORK(KLDB),NU,NP,KWORK(KABD),KWORK(KVERT),
     *               KWORK(KAREA),KWORK(KNPR),RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KVOL),NP,INEUM)
C
       RES=MAX(DMAXU,DMAXP)
       IF (ITE.EQ.1) RESINI=RES
       RHO=(RES/RESINI)**(1D0/DBLE(ITE))
ccc       write(6,*) ite,res,DMAXU,DMAXP,DMPSL*RESINI,RHO
       IF ((RES.LT.EPSSL).AND.(RES.LT.DMPSL*RESINI)) GOTO 99999
C
11     CONTINUE
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE   YEXA (DX,DB,DD,NEQ,RHO,EPS,ITE)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
*            - DX,DB,DD have the structure  D=(D1,D2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C
C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      I3=1+NU+NU
      IP=I3+NU
C
      KVERT=L(LVERT)
      KAREA =L(LAREA)
      KNPR =L(LNPR )
      KABD =L(KLABD(ILEV))
      NABD= KNABD(ILEV)
C
      KVOL=L(KLVOL(ILEV))
C
C=======================================================================
C
      IF (ISL.EQ.1) THEN
C
       DO 11 ITE=1,NSL
C
       CALL  VANCAE (DX(1),DX(I2),DX(I3),DX(IP),DD(1),DD(I2),DD(I3),
     *               DD(IP),DB(1),DB(I2),DB(I3),DB(IP),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               VWORK(KB1),VWORK(KB2),VWORK(KB3),KWORK(KCOLB),
     *               KWORK(KLDB),NU,NP,KWORK(KABD),KWORK(KVERT),
     *               KWORK(KAREA),KWORK(KNPR),RLXSL,DMAXU,DMAXP)
C
      IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KVOL),NP,INEUM)
C
       RES=MAX(DMAXU,DMAXP)
       IF (ITE.EQ.1) RESINI=RES
       RHO=(RES/RESINI)**(1D0/DBLE(ITE))
ccc       write(*,*) ite,res,DMAXU,DMAXP,DMPSL*RESINI,RHO
       IF ((RES.LT.EPSSL).AND.(RES.LT.DMPSL*RESINI)) GOTO 99999
C
11     CONTINUE
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE   YPROL (DUC,DUF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=ALPHA*p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
*  
*            - DUF and DUC have the structure  DU=(DU1,DU2,DU3,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DUF(*),DUC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C=======================================================================
C     Getting all parameters for PROLU
C=======================================================================
C
C *** addresses for the fine level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2F=1+NU
      I3F=1+NU+NU
      IPF=I3F+NU
C
      KAREAF=L(LAREA)
      KADJF =L(LADJ )
      KVERTF=L(LVERT)
C
C *** addresses for the coarse level ILEV-1
      I1=ILEV-1
      NUC=KNU(I1)
      NPC=KNP(I1)
      NATC=KNAT(I1)
      NVTC=KNVT(I1)
      I2C=1+NUC
      I3C=1+NUC+NUC
      IPC=I3C+NUC
C
      KAREAC=L(KLAREA(I1))
      KADJC =L(KLADJ (I1))
      KVERTC=L(KLVERT(I1))
C
C=======================================================================
C
      CALL  PROLU2(DUC(1),DUC(I2C),DUC(I3C),DUC(IPC), 
     *             DUF(1),DUF(I2F),DUF(I3F),DUF(IPF),
     *             KWORK(KAREAC),KWORK(KADJC),NATC,NPC,
     *             KWORK(KAREAF),KWORK(KADJF),NU,NP,IINT)
C
      IF (ABS(IINT).EQ.2) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
C
       CALL C2N2DM(DUC(IPC),DWORK(KPLC),KWORK(KAREAC),KWORK(KADJC),
     *             NPC,NATC,NVTC,0)
C
       IF (IINT.GT.0) THEN
        CALL MP031(DWORK(KPLC),DWORK(KPLF),KWORK(KAREAC),KWORK(KAREAF),
     *             KWORK(KADJC),KWORK(KADJF),NATC,NAT,NPC,NP)
       ELSE
        CALL MP030(DWORK(KPLC),DWORK(KPLF),KWORK(KAREAC),KWORK(KAREAF),
     *             KWORK(KADJC),KWORK(KADJF),NATC,NAT,NPC,NP)
       ENDIF
C
       CALL C2N2DM(DUF(IPF),DWORK(KPLF),KWORK(KAREAF),KWORK(KADJF),
     *             NP,NAT,NVT,1)
      ENDIF
C
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   YREST (DDF,DDC)  
************************************************************************
*
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*  
*            - DDF and DDC have the structure  DD=(DD1,DD2,DDP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DDF(*),DDC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C=======================================================================
C     Getting all parameters for RESTRD
C=======================================================================
C
C *** addresses for the coarse level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2C=1+NU
      I3C=1+NU+NU
      IPC=I3C+NU
C
      KAREAC=L(LAREA)
      KADJC =L(LADJ )
      KVERTC=L(LVERT)
C
C *** addresses for the fine level ILEV+1
      I1=ILEV+1
      NUF=KNU(I1)
      NPF=KNP(I1)
      NATF=KNAT(I1)
      NVTF=KNVT(I1)
      I2F=1+NUF
      I3F=1+NUF+NUF
      IPF=I3F+NUF
C
      KAREAF=L(KLAREA(I1))
      KADJF =L(KLADJ (I1))
      KVERTF=L(KLVERT(I1))
C
C=======================================================================
C
      CALL RESTRM (DDC(1),DDC(I2C),DDC(I3C),DDC(IPC), 
     *             DDF(1),DDF(I2F),DDF(I3F),DDF(IPF),
     *             KWORK(KAREAC),KWORK(KADJC),NAT,NP,
     *             KWORK(KAREAF),KWORK(KADJF),NATF,NPF,IINT)
C
      IF (ABS(IINT).EQ.2) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
C
       CALL C2N2DM(DDF(IPF),DWORK(KPLF),KWORK(KAREAF),KWORK(KADJF),
     *             NPF,NATF,NVTF,0)
C
       IF (IINT.GT.0) THEN
       CALL MR031(DWORK(KPLF),DWORK(KPLC),KWORK(KAREAF),KWORK(KAREAC),
     *            KWORK(KADJF),KWORK(KADJC),NATF,NAT,NPF,NP)
       ELSE
       CALL MR030(DWORK(KPLF),DWORK(KPLC),KWORK(KAREAF),KWORK(KAREAC),
     *            KWORK(KADJF),KWORK(KADJC),NATF,NAT,NPF,NP)
       ENDIF
C
       CALL C2N2DM(DDC(IPC),DWORK(KPLC),KWORK(KAREAC),KWORK(KADJC),
     *             NP,NAT,NVT,1)
      ENDIF
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   YSM (DX,DB,DD,NEQ,NSMO)  
************************************************************************
*
*   Purpose: - performs NSMO smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
*            - DX,DB,DD have the structure  D=(D1,D2,D3,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      I3=1+NU+NU
      IP=I3+NU
C
      KVERT=L(LVERT)
      KAREA=L(LAREA)
      KNPR =L(LNPR )
      KABD =L(KLABD(ILEV))
      NABD =KNABD(ILEV)
C
      KVOL =L(KLVOL(ILEV))
C
C=======================================================================
C
      IF (ISM.EQ.1) THEN
       DO 11  ITE=1,NSMO
       CALL VANCAS (DX(1),DX(I2),DX(I3),DX(IP),DD(1),DD(I2),DD(I3),
     *              DD(IP),DB(1),DB(I2),DB(I3),DB(IP),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              VWORK(KB1),VWORK(KB2),VWORK(KB3),KWORK(KCOLB),
     *              KWORK(KLDB),NU,NP,KWORK(KABD),KWORK(KVERT),
     *              KWORK(KAREA),KWORK(KNPR),DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KVOL),NP,INEUM)
C
11     CONTINUE
      ENDIF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   YSTEP (DX,DD,DB,ALPHA)  
************************************************************************
*
*   Purpose: - performs step size control for prolongation with
*                  DX    - old fine solution vector on level ILEV
*                  DD    - fine correction vector on level ILEV
*                  DB    - fine right hand side vector on level ILEV
*                  ALPHA - relaxation parameter according to some
*                          optimization criterion but within the limits
*                          ALFMIN and ALFMAX (COMMON /RPARM/)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DD(*),DB(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C=======================================================================
C     Getting all parameters for PROLU
C=======================================================================
C
C
C *** Setting of given ALPHA
      IF (AMINMG.EQ.AMAXMG) THEN
       ALPHA=AMINMG
       RETURN
      ENDIF
C
C
C *** Calculation of optimal ALPHA
      IF (AMINMG.NE.AMAXMG) THEN
C
       ISETLV=2
       CALL SETLEV (ISETLV)
C
       CALL LCP1 (DB,DWORK(L(LD1)),NUP)
       CALL YAX(DX,DWORK(L(LD1)),NUP,-1D0,1D0)
       DO 100 IEQ=NUP,NUP-NEL+1,-1
100    DWORK(L(LD1)+IEQ-1)=-DWORK(L(LD1)+IEQ-1)
       CALL LSP1(DD,DWORK(L(LD1)),NUP,DBY)
C
C
       CALL YAX(DD,DWORK(L(LD1)),NUP,1D0,0D0)
       DO 110 IEQ=NUP,NUP-NEL+1,-1
110    DWORK(L(LD1)+IEQ-1)=-DWORK(L(LD1)+IEQ-1)
       CALL LSP1(DD,DWORK(L(LD1)),NUP,DBX)
C
       ALPHA=DBY/DBX
       IF (ALPHA.LT.AMINMG) ALPHA=AMINMG
       IF (ALPHA.GT.AMAXMG) ALPHA=AMAXMG
CCC       WRITE(*,*)'ALPHA=====',ALPHA,DBY,DBX,ILEV
C
       RETURN
C
      ENDIF
C
C
      END
