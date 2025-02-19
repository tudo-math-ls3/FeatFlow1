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
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
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
      SUBROUTINE   PROLU  (DU1,DV1,DP1,DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1,DP1) to
*                the fine grid vector (DU2,DV2,DP2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1,DP1           - coarse grid vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU2,DV2,DP2           - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.46875D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DP1(*),  DU2(*),DV2(*),DP2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C *** Zero initialization of (DU2,DV2)
      CALL  LCL1 (DU2,NEQ2)
      CALL  LCL1 (DV2,NEQ2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      DVH1=DV1(IM1)
      DVH2=DV1(IM2)
      DVH3=DV1(IM3)
      DVH4=DV1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
       DV2(IA)=DV2(IA)+   A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4
       DV2(IB)=DV2(IB)+   A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
      DV2(IC)=A5*DVH1+A6*(DVH2+DVH4)+A7*DVH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
       DV2(IA)=DV2(IA)+   A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1
       DV2(IB)=DV2(IB)+   A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
      DV2(IC)=A5*DVH2+A6*(DVH3+DVH1)+A7*DVH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
       DV2(IA)=DV2(IA)+   A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2
       DV2(IB)=DV2(IB)+   A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
      DV2(IC)=A5*DVH3+A6*(DVH4+DVH2)+A7*DVH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
       DV2(IA)=DV2(IA)+   A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3
       DV2(IB)=DV2(IB)+   A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
      DV2(IC)=A5*DVH4+A6*(DVH1+DVH3)+A7*DVH2
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   PROLUM (DU1,DV1,DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1) to
*                the fine grid vector (DU2,DV2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1               - coarse grid vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU2,DV2               - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.46875D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DU2(*),DV2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C *** Zero initialization of (DU2,DV2)
      CALL  LCL1 (DU2,NEQ2)
      CALL  LCL1 (DV2,NEQ2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      DVH1=DV1(IM1)
      DVH2=DV1(IM2)
      DVH3=DV1(IM3)
      DVH4=DV1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
       DV2(IA)=DV2(IA)+   A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4
       DV2(IB)=DV2(IB)+   A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
      DV2(IC)=A5*DVH1+A6*(DVH2+DVH4)+A7*DVH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
       DV2(IA)=DV2(IA)+   A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1
       DV2(IB)=DV2(IB)+   A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
      DV2(IC)=A5*DVH2+A6*(DVH3+DVH1)+A7*DVH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
       DV2(IA)=DV2(IA)+   A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2
       DV2(IB)=DV2(IB)+   A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
      DV2(IC)=A5*DVH3+A6*(DVH4+DVH2)+A7*DVH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
       DV2(IA)=DV2(IA)+   A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3
       DV2(IB)=DV2(IB)+   A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
      DV2(IC)=A5*DVH4+A6*(DVH1+DVH3)+A7*DVH2
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRD  (DU1,DV1,DP1,DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Restricts the  fine grid defect vector (DU2,DV2,DP2)
*                to the coarse grid defect vector (DU1,DV1,DP1)
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2,DP2           - fine grid defect vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1,DP1           - coarse grid defect vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.9375D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DP1(*),  DU2(*),DV2(*),DP2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Restriction of pressure
C
      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)
C
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
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
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)= A1*(DVH1+DVH2)+A2*(DVH4+DVH7)
     *                 +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                   +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)=DV1(IM1)+A2*(DVH4+DVH7)
     *                        +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                          +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
       DV1(IM1)=     A1*(DVH1+DVH2)+2D0*A2*(DVH4+DVH7)
     *          +2D0*A3*(DVH5+DVH6)+2D0*A4*(DVH3+DVH8)
     *          +    A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)= A1*(DVH3+DVH4)+A2*(DVH6+DVH1)
     *                 +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                   +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)=DV1(IM2)+A2*(DVH6+DVH1)
     *                        +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                          +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
       DV1(IM2)=     A1*(DVH3+DVH4)+2D0*A2*(DVH6+DVH1)
     *          +2D0*A3*(DVH7+DVH8)+2D0*A4*(DVH5+DVH2)
     *          +    A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)= A1*(DVH5+DVH6)+A2*(DVH8+DVH3)
     *                 +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                   +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)=DV1(IM3)+A2*(DVH8+DVH3)
     *                        +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                          +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
       DV1(IM3)=     A1*(DVH5+DVH6)+2D0*A2*(DVH8+DVH3)
     *          +2D0*A3*(DVH1+DVH2)+2D0*A4*(DVH7+DVH4)
     *          +    A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)= A1*(DVH7+DVH8)+A2*(DVH2+DVH5)
     *                 +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                   +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)=DV1(IM4)+A2*(DVH2+DVH5)
     *                        +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                          +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
       DV1(IM4)=     A1*(DVH7+DVH8)+2D0*A2*(DVH2+DVH5)
     *          +2D0*A3*(DVH3+DVH4)+2D0*A4*(DVH1+DVH6)
     *          +    A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
      ENDIF
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRM  (DU1,DV1,DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Restricts the  fine grid defect vector (DU2,DV2)
*                to the coarse grid defect vector (DU1,DV1)
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2               - fine grid defect vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1               - coarse grid defect vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.9375D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DU2(*),DV2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
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
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)= A1*(DVH1+DVH2)+A2*(DVH4+DVH7)
     *                 +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                   +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)=DV1(IM1)+A2*(DVH4+DVH7)
     *                        +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                          +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
       DV1(IM1)=     A1*(DVH1+DVH2)+2D0*A2*(DVH4+DVH7)
     *          +2D0*A3*(DVH5+DVH6)+2D0*A4*(DVH3+DVH8)
     *          +    A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)= A1*(DVH3+DVH4)+A2*(DVH6+DVH1)
     *                 +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                   +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)=DV1(IM2)+A2*(DVH6+DVH1)
     *                        +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                          +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
       DV1(IM2)=     A1*(DVH3+DVH4)+2D0*A2*(DVH6+DVH1)
     *          +2D0*A3*(DVH7+DVH8)+2D0*A4*(DVH5+DVH2)
     *          +    A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)= A1*(DVH5+DVH6)+A2*(DVH8+DVH3)
     *                 +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                   +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)=DV1(IM3)+A2*(DVH8+DVH3)
     *                        +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                          +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
       DV1(IM3)=     A1*(DVH5+DVH6)+2D0*A2*(DVH8+DVH3)
     *          +2D0*A3*(DVH1+DVH2)+2D0*A4*(DVH7+DVH4)
     *          +    A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)= A1*(DVH7+DVH8)+A2*(DVH2+DVH5)
     *                 +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                   +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)=DV1(IM4)+A2*(DVH2+DVH5)
     *                        +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                          +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
       DV1(IM4)=     A1*(DVH7+DVH8)+2D0*A2*(DVH2+DVH5)
     *          +2D0*A3*(DVH3+DVH4)+2D0*A4*(DVH1+DVH6)
     *          +    A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
      ENDIF
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRU  (DU1,DV1,DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2 )
************************************************************************
*    Purpose:  Restricts the  fine grid solution vector (DU2,DV2) to
*              the coarse grid solution vector (DU1,DV1)
*
*    Remark:   RESTRU works for all cases of boundary conditions
*                
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2               - fine grid solution vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1               - coarse grid solution vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.1875D0, A2=0.375D0, A3=-0.0625D0)
      PARAMETER (R1=0.375D0, R2=0.75D0, R3=-0.125D0)
C
      DIMENSION DU1(*),DV1(*),  DU2(*),DV2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
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
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DU1(IM1)=A1*(DUH1+DUH2)  +A2*DUH9   +A3*(DUH8+DUH3+DUH10+DUH12)
        DV1(IM1)=A1*(DVH1+DVH2)  +A2*DVH9   +A3*(DVH8+DVH3+DVH10+DVH12)
       ELSE
        DU1(IM1)=DU1(IM1) +
     *           A1*(DUH1+DUH2)  +A2*DUH9   +A3*(DUH8+DUH3+DUH10+DUH12)
        DV1(IM1)=DV1(IM1) +
     *           A1*(DVH1+DVH2)  +A2*DVH9   +A3*(DVH8+DVH3+DVH10+DVH12)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM1)=R1*(DUH1+DUH2)  +R2*DUH9   +R3*(DUH8+DUH3+DUH10+DUH12)
        DV1(IM1)=R1*(DVH1+DVH2)  +R2*DVH9   +R3*(DVH8+DVH3+DVH10+DVH12)
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DU1(IM2)=A1*(DUH3+DUH4)  +A2*DUH10  +A3*(DUH2+DUH5+DUH9 +DUH11)
        DV1(IM2)=A1*(DVH3+DVH4)  +A2*DVH10  +A3*(DVH2+DVH5+DVH9 +DVH11)
       ELSE
        DU1(IM2)=DU1(IM2) +
     *           A1*(DUH3+DUH4)  +A2*DUH10  +A3*(DUH2+DUH5+DUH9 +DUH11)
        DV1(IM2)=DV1(IM2) +
     *           A1*(DVH3+DVH4)  +A2*DVH10  +A3*(DVH2+DVH5+DVH9 +DVH11)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM2)=R1*(DUH3+DUH4)  +R2*DUH10  +R3*(DUH2+DUH5+DUH9 +DUH11)
        DV1(IM2)=R1*(DVH3+DVH4)  +R2*DVH10  +R3*(DVH2+DVH5+DVH9 +DVH11)
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DU1(IM3)=A1*(DUH5+DUH6)  +A2*DUH11  +A3*(DUH4+DUH7+DUH10+DUH12)
        DV1(IM3)=A1*(DVH5+DVH6)  +A2*DVH11  +A3*(DVH4+DVH7+DVH10+DVH12)
       ELSE
        DU1(IM3)=DU1(IM3) +
     *           A1*(DUH5+DUH6)  +A2*DUH11  +A3*(DUH4+DUH7+DUH10+DUH12)
        DV1(IM3)=DV1(IM3) +
     *           A1*(DVH5+DVH6)  +A2*DVH11  +A3*(DVH4+DVH7+DVH10+DVH12)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM3)=R1*(DUH5+DUH6)  +R2*DUH11  +R3*(DUH4+DUH7+DUH10+DUH12)
        DV1(IM3)=R1*(DVH5+DVH6)  +R2*DVH11  +R3*(DVH4+DVH7+DVH10+DVH12)
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DU1(IM4)=A1*(DUH7+DUH8)  +A2*DUH12  +A3*(DUH6+DUH1+DUH9 +DUH11)
        DV1(IM4)=A1*(DVH7+DVH8)  +A2*DVH12  +A3*(DVH6+DVH1+DVH9 +DVH11)
       ELSE
        DU1(IM4)=DU1(IM4) +
     *           A1*(DUH7+DUH8)  +A2*DUH12  +A3*(DUH6+DUH1+DUH9 +DUH11)
        DV1(IM4)=DV1(IM4) +
     *           A1*(DVH7+DVH8)  +A2*DVH12  +A3*(DVH6+DVH1+DVH9 +DVH11)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM4)=R1*(DUH7+DUH8)  +R2*DUH12  +R3*(DUH6+DUH1+DUH9 +DUH11)
        DV1(IM4)=R1*(DVH7+DVH8)  +R2*DVH12  +R3*(DVH6+DVH1+DVH9 +DVH11)
      ENDIF
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MR010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NNVE=4)
      DIMENSION DP1(*),DP2(*)
      DIMENSION KADJ1(NNVE,*),KADJ2(NNVE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Restriction of pressure
C
      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MP030(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *                 KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.5D0,A2=-0.0625D0,A3=0D0,A4=0.0625D0)
      PARAMETER (A5=0.625D0,A6=0.125D0,A7=0.125D0,A8=0.125D0)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KVERT2(NNVE,*),
     *          KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *          KADJ2(NNVE,*)
C
C
C
C *** Zero initialization of (DU2)
      CALL  LCL1 (DU2,NMT2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
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
      SUBROUTINE MP031(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *                 KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.46875D0,A2=-0.09375D0,A3=-0.03125D0,A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KVERT2(NNVE,*),
     *          KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *          KADJ2(NNVE,*)
C
C
C
C *** Zero initialization of (DU2)
      CALL  LCL1 (DU2,NMT2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
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
      SUBROUTINE MR030(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *                 KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=1D0,A2=-0.0625D0,A3=0D0,A4=0.0625D0)
      PARAMETER (A5=0.625D0,A6=0.125D0,A7=0.125D0,A8=0.125D0)
      DIMENSION DU1(1),DU2(1),KVERT2(NNVE,1),KVERT1(NNVE,1),
     *          KMID1(NNVE,1),KMID2(NNVE,1),
     *          KADJ1(NNVE,1),KADJ2(NNVE,1)
C
C
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
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
C
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF
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
      SUBROUTINE MP010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION DP1(*),DP2(*)
      DIMENSION KADJ1(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MR031(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *                 KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.9375D0,A2=-0.09375D0,A3=-0.03125D0,A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
      DIMENSION DU1(1),DU2(1),KVERT2(NNVE,1),KVERT1(NNVE,1),
     *          KMID1(NNVE,1),KMID2(NNVE,1),
     *          KADJ1(NNVE,1),KADJ2(NNVE,1)
C
C
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
C
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
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
C
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF
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
      SUBROUTINE  YAX (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
*
*            - DX,DAX  have the structure  D=(D1,D2,DP)
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
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
      IP=I2+NU
C
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
C=======================================================================
C
      CALL MATML(DAX(1),DAX(I2),DAX(IP),DX(1),DX(I2),DX(IP),A1,A2,
     *            VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *            NU,NP,KWORK(KMBD),NMBD,INEUM)
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE  YDBC (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX=(D1,D2) to 0
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
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
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
C
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      CALL  BDRY0 (DX(1),DX(I2),KWORK(KMBD),NMBD)
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   YEX (DX,DB,DD,NEQ,RHO)  
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
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      IF (ISL.EQ.1) THEN
C
       DO 11 ITE=1,NSL
C
       CALL  VANCAE (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *               DB(1),DB(I2),DB(IP),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *               NU,NP,KWORK(KMBD),KWORK(KVERT),
     *               KWORK(KMID),KWORK(KNPR),NMBD,RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
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
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      IF (ISL.EQ.1) THEN
C
       DO 11 ITE=1,NSL
C
       CALL  VANCAE (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *               DB(1),DB(I2),DB(IP),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *               NU,NP,KWORK(KMBD),KWORK(KVERT),
     *               KWORK(KMID),KWORK(KNPR),NMBD,RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
C
       RES=MAX(DMAXU,DMAXP)
       IF (ITE.EQ.1) RESINI=RES
       RHO=(RES/RESINI)**(1D0/DBLE(ITE))
ccc       write(6,*) ite,res,DMAXU,DMAXP,DMPSL*RESINI,RHO
       IF ((RES.LT.EPSMG).AND.(RES.LT.DMPMG*RESINI)) GOTO 99999
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
      SUBROUTINE  YPROL (DUC,DUF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
*  
*            - DUF and DUC have the structure  DU=(DU1,DU2,DP)
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
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
      IPF=I2F+NU
C
      KVERTF=L(LVERT)
      KMIDF =L(LMID )
      KADJF =L(LADJ )
C
C *** addresses for the coarse level ILEV-1
      I1=ILEV-1
      NUC=KNU(I1)
      NPC=KNP(I1)
      NVTC=KNVT(I1)
      NMTC=KNMT(I1)
      I2C=1+NUC
      IPC=I2C+NUC
C
      KVERTC=L(KLVERT(I1))
      KMIDC =L(KLMID (I1))
      KADJC =L(KLADJ (I1))
C
C=======================================================================
C
      CALL PROLU (DUC(1),DUC(I2C),DUC(IPC),DUF(1),DUF(I2F),DUF(IPF),
     *            KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),NUC,NPC,NVTC,
     *            KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),NU,NP,NVT)
C
      IF (ABS(IINT).EQ.2) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
C
       CALL C2N2DM(DUC(IPC),DWORK(KPLC),KWORK(KMIDC),
     *             KWORK(KADJC),NPC,NMTC,NVTC,0)
C
       IF (IINT.GT.0) THEN
        CALL MP031(DWORK(KPLC),DWORK(KPLF),KWORK(KVERTC),KWORK(KVERTF),
     *             KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
     *             NVTC,NVT,NPC,NP,NMT)
       ELSE
        CALL MP030(DWORK(KPLC),DWORK(KPLF),KWORK(KVERTC),KWORK(KVERTF),
     *             KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
     *             NVTC,NVT,NPC,NP,NMT)
       ENDIF
C
       CALL C2N2DM(DUF(IPF),DWORK(KPLF),KWORK(KMIDF),
     *             KWORK(KADJF),NP,NMT,NVT,1)
      ENDIF
C
C
C
      END
C
C
C
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
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
      IPC=I2C+NU
C
      KVERTC=L(LVERT)
      KMIDC =L(LMID )
      KADJC =L(LADJ )
C
C *** addresses for the fine level ILEV+1
      I1=ILEV+1
      NUF=KNU(I1)
      NPF=KNP(I1)
      NVTF=KNVT(I1)
      NMTF=KNMT(I1)
      I2F=1+NUF
      IPF=I2F+NUF
C
      KVERTF=L(KLVERT(I1))
      KMIDF =L(KLMID (I1))
      KADJF =L(KLADJ (I1))
C
C=======================================================================
C
      CALL RESTRD(DDC(1),DDC(I2C),DDC(IPC),DDF(1),DDF(I2F),DDF(IPF),
     *            KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),NU,NP,NVT,
     *            KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),NUF,NPF,NVTF)
C
      IF (ABS(IINT).EQ.2) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
C
       CALL C2N2DM(DDF(IPF),DWORK(KPLF),KWORK(KMIDF),KWORK(KADJF),
     *             NPF,NMTF,NVTF,0)
C
       IF (IINT.GT.0) THEN
        CALL MR031(DWORK(KPLF),DWORK(KPLC),KWORK(KVERTF),KWORK(KVERTC),
     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
     *             NVTF,NVT,NPF,NP)
       ELSE
        CALL MR030(DWORK(KPLF),DWORK(KPLC),KWORK(KVERTF),KWORK(KVERTC),
     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
     *             NVTF,NVT,NPF,NP)
       ENDIF
C
       CALL C2N2DM(DDC(IPC),DWORK(KPLC),KWORK(KMIDC),KWORK(KADJC),
     *             NP,NMT,NVT,1)
      ENDIF
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE  YSM (DX,DB,DD,NEQ,NSMO)  
************************************************************************
*
*   Purpose: - performs NSMO smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
*            - DX,DB,DD have the structure  D=(D1,D2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
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
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      IF (ISM.EQ.1) THEN
       DO 11  ITE=1,NSMO
       CALL VANCAS (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *              DB(1),DB(I2),DB(IP),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *              NU,NP,KWORK(KMBD),KWORK(KVERT),KWORK(KMID),
     *              KWORK(KNPR),NMBD)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
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
*                          AMINMG and AMAXMG (COMMON /RPARM/)
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
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
ccc       WRITE(*,*)'ALPHA=====',ALPHA,DBY,DBX,ILEV
C
       RETURN
C
      ENDIF
C
C
      END
