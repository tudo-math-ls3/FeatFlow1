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
* EPS1      R*8   Desired precision                                    *
*                 Stop if !!DEFn!! < EPS1 !!DEF0!!                     *
* EPS2      R*8   Desired precision                                    *
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
     *                  EPS1,EPS2,DEF,DAX,DPROL,DREST,DPRSM,DPOSM,
     *                  DEX,DEXA,DBC,DSTEP,KIT0,KIT,IREL,IDEFMG,RHOLMG,
     *                  BMGEND)
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
c
c     Wir fangen mit null an!
c
ccccc      CALL LCL1(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
c
c
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
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX (DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *           -1D0,1D0)
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
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
C ***  restriction of defect
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
      CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),
     *        DD(1+KOFFD(NLMIN)),KNEQ(NLMIN),RHOLMG)
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
      CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),
     *           DB(1+KOFFB(ILEV)),KNEQ(ILEV),DSTEPP)
      CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),
     *          DSTEPP,1D0)
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
     * CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *            DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPOSM(ILEV))
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
c       IF (MT.GE.0) WRITE (6,10001) ITE,DEF
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
      SUBROUTINE   PROLU  (DU1,DV1,DP1,  DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,  NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,  NEQ2,NEL2,NVT2 )
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
      SUBROUTINE  RESTRD  (DU1,DV1,DP1,  DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,  NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,  NEQ2,NEL2 ,NVT2)
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
      SUBROUTINE  RESTRU  (DU1,DV1,  DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,  NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,  NEQ2,NEL2,NVT2,
     *                     ICOMP )
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
*      ICOMP                 - number of vector components
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
      IF (ICOMP.NE.1) THEN
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
      ENDIF
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DU1(IM1)=A1*(DUH1+DUH2)  +A2*DUH9   +A3*(DUH8+DUH3+DUH10+DUH12)
        IF (ICOMP.NE.1) 
     * DV1(IM1)=A1*(DVH1+DVH2)  +A2*DVH9   +A3*(DVH8+DVH3+DVH10+DVH12)
       ELSE
        DU1(IM1)=DU1(IM1) +
     *           A1*(DUH1+DUH2)  +A2*DUH9   +A3*(DUH8+DUH3+DUH10+DUH12)
        IF (ICOMP.NE.1) 
     *  DV1(IM1)=DV1(IM1) +
     *           A1*(DVH1+DVH2)  +A2*DVH9   +A3*(DVH8+DVH3+DVH10+DVH12)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM1)=R1*(DUH1+DUH2)  +R2*DUH9   +R3*(DUH8+DUH3+DUH10+DUH12)
        IF (ICOMP.NE.1) 
     *   DV1(IM1)=R1*(DVH1+DVH2) +R2*DVH9   +R3*(DVH8+DVH3+DVH10+DVH12)
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DU1(IM2)=A1*(DUH3+DUH4)  +A2*DUH10  +A3*(DUH2+DUH5+DUH9 +DUH11)
        IF (ICOMP.NE.1) 
     *   DV1(IM2)=A1*(DVH3+DVH4)  +A2*DVH10  +A3*(DVH2+DVH5+DVH9 +DVH11)
       ELSE
        DU1(IM2)=DU1(IM2) +
     *           A1*(DUH3+DUH4)  +A2*DUH10  +A3*(DUH2+DUH5+DUH9 +DUH11)
        IF (ICOMP.NE.1) 
     *   DV1(IM2)=DV1(IM2) +
     *           A1*(DVH3+DVH4)  +A2*DVH10  +A3*(DVH2+DVH5+DVH9 +DVH11)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM2)=R1*(DUH3+DUH4)  +R2*DUH10  +R3*(DUH2+DUH5+DUH9 +DUH11)
        IF (ICOMP.NE.1) 
     *   DV1(IM2)=R1*(DVH3+DVH4)  +R2*DVH10  +R3*(DVH2+DVH5+DVH9 +DVH11)
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DU1(IM3)=A1*(DUH5+DUH6)  +A2*DUH11  +A3*(DUH4+DUH7+DUH10+DUH12)
        IF (ICOMP.NE.1) 
     *   DV1(IM3)=A1*(DVH5+DVH6)  +A2*DVH11  +A3*(DVH4+DVH7+DVH10+DVH12)
       ELSE
        DU1(IM3)=DU1(IM3) +
     *           A1*(DUH5+DUH6)  +A2*DUH11  +A3*(DUH4+DUH7+DUH10+DUH12)
        IF (ICOMP.NE.1) 
     *   DV1(IM3)=DV1(IM3) +
     *           A1*(DVH5+DVH6)  +A2*DVH11  +A3*(DVH4+DVH7+DVH10+DVH12)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM3)=R1*(DUH5+DUH6)  +R2*DUH11  +R3*(DUH4+DUH7+DUH10+DUH12)
        IF (ICOMP.NE.1) 
     *   DV1(IM3)=R1*(DVH5+DVH6)  +R2*DVH11  +R3*(DVH4+DVH7+DVH10+DVH12)
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DU1(IM4)=A1*(DUH7+DUH8)  +A2*DUH12  +A3*(DUH6+DUH1+DUH9 +DUH11)
        IF (ICOMP.NE.1) 
     *   DV1(IM4)=A1*(DVH7+DVH8)  +A2*DVH12  +A3*(DVH6+DVH1+DVH9 +DVH11)
       ELSE
        DU1(IM4)=DU1(IM4) +
     *           A1*(DUH7+DUH8)  +A2*DUH12  +A3*(DUH6+DUH1+DUH9 +DUH11)
        IF (ICOMP.NE.1) 
     *   DV1(IM4)=DV1(IM4) +
     *           A1*(DVH7+DVH8)  +A2*DVH12  +A3*(DVH6+DVH1+DVH9 +DVH11)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM4)=R1*(DUH7+DUH8)  +R2*DUH12  +R3*(DUH6+DUH1+DUH9 +DUH11)
        IF (ICOMP.NE.1) 
     *   DV1(IM4)=R1*(DVH7+DVH8)  +R2*DVH12  +R3*(DVH6+DVH1+DVH9 +DVH11)
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
      SUBROUTINE MR010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
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
C
C
C
************************************************************************
      SUBROUTINE YAXU (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
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
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DAX,NEQ)
      ENDIF
C
C
      CALL LAX37(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *           KWORK(L(KLLDA(ILEV))),NEQ,DX,DAX,A1,A2)
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DAX,NEQ)
      ENDIF
C
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE YDBCU (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX to zero
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
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      CALL  BDRY0(DX,DX,KWORK(KMBD),NMBD)
C
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE YEXU (DX,DB,DD,NEQ,RHO)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX37,YIA137,YID137,YIF137
C
      SAVE
C
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.1) THEN
       CALL IC037(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSLU,ITE,DMPUSL,
     *            RLXSLU)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.2) THEN
C
      OMEGA=RLXSLU
C
      IREQ=5*NEQ
      BNOCON=OMEGA.GT.0D0
      IREQ=MAX(IREQ,5)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      L5=L4+NEQ
C
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE IF (OMEGA.LT.2D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YID137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,I000,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.3) THEN
       CALL IF027(VWORK(L(KLA(ILEV))),VWORK(L(KLAILU(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDA(ILEV))),
     *            DX,DB,DD,NEQ,NSLU,ITE,DMPUSL,RLXSLU,RHO)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.4) THEN
C
      OMEGA=RLXSLU
C
      IREQ=5*NEQ
      BNOCON=.TRUE.
      IREQ=MAX(IREQ,5)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      L5=L4+NEQ
C
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      KXYPAR(4)=L(KLAILU(ILEV))
      DXYPAR(1)=RLXSLU
C
      CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIF137,BNOCON,
     *           DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *           DWORK(L5),1,RHO)
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YEXAU (DX,DB,DD,NEQ,RHO,EPS,ITECG)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX37,YIA137,YID137,YIF137
C
      SAVE
C
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.1) THEN
       CALL IC037(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSLU,ITE,DMPUSL,
     *            RLXSLU)
       RHO  =1D0
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.2) THEN
C
      OMEGA=RLXSLU
C
      IREQ=1*NEQ
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
C
      L1=L(LWORK)
      L2=KAUX1
      L3=KAUX2
      L4=L(LD1)
      L5=L(LD2)
C
      BNOCON=OMEGA.GT.0D0
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE IF (OMEGA.LT.2D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YID137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,I000,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ENDIF
      ITECG=ITE
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
      IF (ISLU.EQ.3) THEN
       CALL IF027(VWORK(L(KLA(ILEV))),VWORK(L(KLAILU(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDA(ILEV))),
     *            DX,DB,DD,NEQ,NSLU,ITE,DMPUSL,RLXSLU,RHO)
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.4) THEN
C
      IREQ=1*NEQ
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
C
      L1=L(LWORK)
      L2=KAUX1
      L3=KAUX2
      L4=L(LD1)
      L5=L(LD2)
C
      BNOCON=.TRUE.
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      KXYPAR(4)=L(KLAILU(ILEV))
      DXYPAR(1)=RLXSLU
C
      CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIF137,BNOCON,
     *           DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *           DWORK(L5),1,RHO)
      ITECG=ITE
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YPROLU (DUC,DUF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DUF(*),DUC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
      LLV1=L(KLVERT(ILEV-1))
      LLV2=L(KLVERT(ILEV))
      LLM1=L(KLMID(ILEV-1))
      LLM2=L(KLMID(ILEV))
      LLA1=L(KLADJ(ILEV-1))
      LLA2=L(KLADJ(ILEV))
      NVT1=KNVT(ILEV-1)
      NVT2=KNVT(ILEV)
      NEL1=KNEL(ILEV-1)
      NEL2=KNEL(ILEV)
      NMT1=KNMT(ILEV-1)
      NMT2=KNMT(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTU.EQ.1) THEN
       CALL MP031(DUC,DUF,KWORK(LLV1),KWORK(LLV2),
     *            KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),
     *            KWORK(LLA2),NVT1,NVT2,NEL1,NEL2,NMT2)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.2) THEN
       CALL MP030(DUC,DUF,KWORK(LLV1),KWORK(LLV2),
     *            KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),
     *            KWORK(LLA2),NVT1,NVT2,NEL1,NEL2,NMT2)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YRESTU(DDF,DDC)  
************************************************************************
*
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DDF(*),DDC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON block
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      LLV2=L(KLVERT(ILEV+1))
      LLV1=L(KLVERT(ILEV))
      LLA2=L(KLADJ(ILEV+1))
      LLA1=L(KLADJ(ILEV))
      LLM2=L(KLMID(ILEV+1))
      LLM1=L(KLMID(ILEV))
      NVT2=KNVT(ILEV+1)
      NVT1=KNVT(ILEV)
      NEL2=KNEL(ILEV+1)
      NEL1=KNEL(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTU.EQ.1) THEN
       CALL MR031(DDF,DDC,KWORK(LLV2),KWORK(LLV1),
     *            KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),
     *            KWORK(LLA1),NVT2,NVT1,NEL2,NEL1)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.2) THEN
       CALL MR030(DDF,DDC,KWORK(LLV2),KWORK(LLV1),
     *            KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),
     *            KWORK(LLA1),NVT2,NVT1,NEL2,NEL1)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YSMU (DX,DB,DD,NEQ,NSM)  
************************************************************************
*
*   Purpose: - performs NSM smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISMU.EQ.1) THEN
C
      CALL IA237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *           KWORK(L(KLLDA(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMU)
      GOTO 999
C
      ENDIF
C-----------------------------------------------------------------------
      IF (ISMU.EQ.2) THEN
C
      IF (RLXSMU.EQ.1D0) THEN
       CALL IB237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSM)
      ELSE
       CALL IC237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSM,RLXSMU)
      ENDIF
      GOTO 999
C
      ENDIF
C-----------------------------------------------------------------------
      IF (ISMU.EQ.3) THEN
       CALL ID237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSM,RLXSMU)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
      IF (ISMU.EQ.4) THEN
       CALL IF227(VWORK(L(KLA(ILEV))),VWORK(L(KLAILU(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDA(ILEV))),
     *            DX,DB,DD,NEQ,NSM,RLXSMU)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C
C
999   IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C 
99999 END
c
c
c
************************************************************************
      SUBROUTINE YSTEPU(DX,DD,DB,NEQ,ALPHA)  
************************************************************************
*
*   Purpose: - performs step size control for prolongation with
*                  DX    - old fine solution vector on level ILEV
*                  DD    - fine correction vector on level ILEV
*                  DB    - fine right hand side vector on level ILEV
*                  ALPHA - relaxation parameter according to some
*                          optimization criterion but within the limits
*                          AMINU and AMAXU (COMMON /RPARM/)
*  
*            - adaptive strategie not yet implemented
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DD(*),DB(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
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
C *** Setting of given ALPHA
      IF (AMINU.EQ.AMAXU) THEN
       ALPHA=AMINU
       RETURN
      ENDIF
C
C
C *** Calculation of optimal ALPHA
      IF (AMINU.NE.AMAXU) THEN
C
       IF (ISORTU.GT.0) THEN
        KTRA1=L(KLTRA1(ILEV))
        KTRA2=L(KLTRA2(ILEV))
        CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
        CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
        CALL VECSRT(DD,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(L(LD1)),DD,NEQ)
       ENDIF
C
       CALL VYAX7(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),NEQ,DD,DX,DD,DYAX1,DYAX2)
       CALL LSP1(DB,DD,NEQ,DBX)
       ALPHA=(DBX-DYAX1)/DYAX2
       IF (ALPHA.LT.AMINU) ALPHA=AMINU
       IF (ALPHA.GT.AMAXU) ALPHA=AMAXU
C
       IF (ISORTU.GT.0) THEN
        KTRA1=L(KLTRA1(ILEV))
        KTRA2=L(KLTRA2(ILEV))
        CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
        CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
        CALL VECSRT(DD,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(L(LD1)),DD,NEQ)
       ENDIF
C
      ENDIF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE YAXP (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
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
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C-----------------------------------------------------------------------
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DAX,NEQ)
      ENDIF
C
C
      CALL LAX17(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *           KWORK(L(KLLDC(ILEV))),NEQ,DX,DAX,A1,A2)
      CALL TOL20A(DAX,VWORK(L(KLAREA(ILEV))),NEQ,INEUM)
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DAX,NEQ)
      ENDIF
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YDBCP (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX to zero
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
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      END
c
c
c
************************************************************************
      SUBROUTINE YEXP (DX,DB,DD,NEQ,RHO)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX17,YIA117,YID117,YIF137
C
      SAVE
C
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.1) THEN
       CALL IC017(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSLP,ITE,DMPPSL,
     *            RLXSLP)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.2) THEN
C
       OMEGA=RLXSLP
C
       IREQ=4*NEQ
       BNOCON=OMEGA.LT.0D0
       IF (BNOCON) IREQ=3*NEQ
       IREQ=MAX(IREQ,4)
       CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
       IF (IER.NE.0) GOTO 99999
       L1=L(LWORK)
       L2=L1+NEQ
       L3=L2+NEQ
       L4=L3+NEQ
       IF (BNOCON) L4=L1
C
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       DXYPAR(1)=OMEGA
C
       IF (OMEGA.EQ.0D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,YIA117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
       ELSE IF (OMEGA.LT.2D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,YID117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
       ELSE
        CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,I000,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
       ENDIF
C
       IER1=IER
       CALL ZDISP(0,LWORK,'WORKCG')
       IER=IER1
       GOTO 999
C
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.3) THEN
       CALL IF037(DWORK(L(KLC(ILEV))),VWORK(L(KLCILU(ILEV))),
     *            KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *            DX,DB,DD,NEQ,NSLP,ITE,DMPPSL,RLXSLP,RHO)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.4) THEN
C
       IREQ=4*NEQ
       BNOCON=.FALSE.
       IREQ=MAX(IREQ,4)
       CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
       IF (IER.NE.0) GOTO 99999
       L1=L(LWORK)
       L2=L1+NEQ
       L3=L2+NEQ
       L4=L3+NEQ
C
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       KXYPAR(4)=L(KLCILU(ILEV))
       DXYPAR(1)=RLXSLP
C
       CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,YIF137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
C
       IER1=IER
       CALL ZDISP(0,LWORK,'WORKCG')
       IER=IER1
       GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
      CALL TOL20A(DX,VWORK(L(KLAREA(ILEV))),NEQ,INEUM)
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YEXAP (DX,DB,DD,NEQ,RHO,EPS,ITECG)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX17,YIA117,YID117,YIF137
C
      SAVE
C
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.1) THEN
       CALL IC017(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSLP,ITE,EPS,RLXSLP)
       RHO  =1D0
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.2) THEN
C
       OMEGA=RLXSLP
C
       L1=KAUX1
       L2=KAUX2
       L3=L(LD1)
       L4=L(LD2)
C
       BNOCON=OMEGA.LT.0D0
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       DXYPAR(1)=OMEGA
C
       IF (OMEGA.EQ.0D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,YIA117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ELSE IF (OMEGA.LT.2D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,YID117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ELSE
        CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,I000,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ENDIF
       ITECG=ITE
C
       GOTO 999
C
      ENDIF 
C
C-----------------------------------------------------------------------
C
      IF (ISLP.EQ.3) THEN
       CALL IF037(DWORK(L(KLC(ILEV))),VWORK(L(KLCILU(ILEV))),
     *            KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *            DX,DB,DD,NEQ,NSLP,ITE,EPS,RLXSLP,RHO)
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
C
      IF (ISLP.EQ.4) THEN
C
       L1=KAUX1
       L2=KAUX2
       L3=L(LD1)
       L4=L(LD2)
C
       BNOCON=.FALSE.
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       KXYPAR(4)=L(KLCILU(ILEV))
       DXYPAR(1)=RLXSLP
C
       CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,YIF137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ITECG=ITE
C
       GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
      CALL TOL20A(DX,VWORK(L(KLAREA(ILEV))),NEQ,INEUM)
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YPROLP (DPC,DPF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DPF(*),DPC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
      LLV1=L(KLVERT(ILEV-1))
      LLV2=L(KLVERT(ILEV))
      LLM1=L(KLMID(ILEV-1))
      LLM2=L(KLMID(ILEV))
      LLA1=L(KLADJ(ILEV-1))
      LLA2=L(KLADJ(ILEV))
      NVT1=KNVT(ILEV-1)
      NVT2=KNVT(ILEV)
      NEL1=KNEL(ILEV-1)
      NEL2=KNEL(ILEV)
      NMT1=KNMT(ILEV-1)
      NMT2=KNMT(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.1) THEN
C
       CALL MP010(DPC,DPF,KWORK(LLA1),KWORK(LLA2),NEL1,NEL2)
       GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.2) THEN
C
      KPLC=L(LD1)
      KPLF=L(LD2)
      CALL ZNEW(NVT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DPC,DWORK(KPLC),VWORK(L(KLAREA(ILEV-1))),
     *           VWORK(L(LAUXCL)),KWORK(LLV1),NEL1,NVT1,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MP011(DWORK(KPLC),DWORK(KPLF),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLA1),KWORK(LLA2),NVT1,NEL1)
C
      CALL ZNEW(NVT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DPF,DWORK(KPLF),VWORK(L(KLAREA(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLV2),NEL2,NVT2,1)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.3) THEN
C
      CALL ZNEW(NMT1,1,LPLC  ,'DPLC  ')
      CALL ZNEW(NMT2,1,LPLF  ,'DPLF  ')
      CALL ZNEW(NMT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DPC,DWORK(L(LPLC)),VWORK(L(KLAREA(ILEV-1))),
     *           VWORK(L(LAUXCL)),KWORK(LLM1),NEL1,NMT1,NVT1,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MP031(DWORK(L(LPLC)),DWORK(L(LPLF)),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           NVT1,NVT2,NEL1,NEL2,NMT2)
C
      CALL ZNEW(NMT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DPF,DWORK(L(LPLF)),VWORK(L(KLAREA(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLM2),NEL2,NMT2,NVT2,1)
C     
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      CALL ZDISP (0,LPLF  ,'DPLF  ')
      CALL ZDISP (0,LPLC  ,'DPLC  ')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.4) THEN
C
      KPLC=L(LD1)
      KPLF=L(LD2)
C
      CALL C2N2DM(DPC,DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *            NEL1,NMT1,NVT1,0)
C
      CALL MP031(DWORK(KPLC),DWORK(KPLF),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           NVT1,NVT2,NEL1,NEL2,NMT2)
C
      CALL C2N2DM(DPF,DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *            NEL2,NMT2,NVT2,1)
C
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C
C
99999 CALL TOL20A(DPF,VWORK(L(KLAREA(ILEV))),NEL2,INEUM)
C
      END
c
c
c
************************************************************************
      SUBROUTINE YRESTP (DDF,DDC)  
************************************************************************
*
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DDF(*),DDC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON block
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      LLV2=L(KLVERT(ILEV+1))
      LLV1=L(KLVERT(ILEV))
      LLA2=L(KLADJ(ILEV+1))
      LLA1=L(KLADJ(ILEV))
      LLM2=L(KLMID(ILEV+1))
      LLM1=L(KLMID(ILEV))
      NVT2=KNVT(ILEV+1)
      NVT1=KNVT(ILEV)
      NMT2=KNMT(ILEV+1)
      NMT1=KNMT(ILEV)
      NEL2=KNEL(ILEV+1)
      NEL1=KNEL(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.1) THEN
C
       CALL MR010(DDC,DDF,KWORK(LLA1),KWORK(LLA2),NEL1,NEL2)
       GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.2) THEN
C
      KPLC=L(LD1)
      KPLF=L(LD2)
      CALL ZNEW(NVT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DDF,DWORK(KPLF),VWORK(L(KLAREA(ILEV+1))),
     *           VWORK(L(LAUXCL)),KWORK(LLV2),NEL2,NVT2,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MR011(DWORK(KPLF),DWORK(KPLC),KWORK(LLV2),
     *           KWORK(LLA2),NVT1,NEL2)
C
      CALL ZNEW(NVT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DDC,DWORK(KPLC),VWORK(L(KLAREA(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLV1),NEL1,NVT1,1)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.3) THEN
C
      CALL ZNEW(NMT1,1,LPLC  ,'DPLC  ')
      CALL ZNEW(NMT2,1,LPLF  ,'DPLF  ')
      CALL ZNEW(NMT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DDF,DWORK(L(LPLF)),VWORK(L(KLAREA(ILEV+1))),
     *           VWORK(L(LAUXCL)),KWORK(LLM2),NEL2,NMT2,NVT2,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MR031(DWORK(L(LPLF)),DWORK(L(LPLC)),KWORK(LLV2),KWORK(LLV1),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           NVT2,NVT1,NEL2,NEL1)
C     
      CALL ZNEW(NMT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DDC,DWORK(L(LPLC)),VWORK(L(KLAREA(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLM1),NEL1,NMT1,NVT1,1)
C     
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      CALL ZDISP (0,LPLF  ,'DPLF  ')
      CALL ZDISP (0,LPLC  ,'DPLC  ')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.4) THEN
C
      KPLC=L(LD1)
      KPLF=L(LD2)
C
      CALL C2N2DM(DDF,DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *            NEL2,NMT2,NVT2,0)
C
      CALL MR031(DWORK(KPLF),DWORK(KPLC),KWORK(LLV2),KWORK(LLV1),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           NVT2,NVT1,NEL2,NEL1)
C
      CALL C2N2DM(DDC,DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *            NEL1,NMT1,NVT1,1)
C
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C
C
99999 CALL TOL20A(DDC,VWORK(L(KLAREA(ILEV))),NEL1,INEUM)
C
      END
c
c
c
************************************************************************
      SUBROUTINE YSMP (DX,DB,DD,NEQ,NSM)  
************************************************************************
*
*   Purpose: - performs NSM smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.1) THEN
       CALL IA217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMP)
       GOTO 999
      ENDIF
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.2) THEN
       IF (RLXSMP.EQ.1D0) THEN
        CALL IB217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSM)
       ELSE
        CALL IC217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSM,RLXSMP)
       ENDIF
       GOTO 999
      ENDIF
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.3) THEN
       CALL ID217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSM,RLXSMP)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.4) THEN
       CALL IF237(DWORK(L(KLC(ILEV))),VWORK(L(KLCILU(ILEV))),
     *            KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *            DX,DB,DD,NEQ,NSM,RLXSMP)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C
C
999   IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
      CALL TOL20A(DX,VWORK(L(KLAREA(ILEV))),NEQ,INEUM)
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YSTEPP (DX,DD,DB,NEQ,ALPHA)  
************************************************************************
*
*   Purpose: - performs step size control for prolongation with
*                  DX    - old fine solution vector on level ILEV
*                  DD    - fine correction vector on level ILEV
*                  DB    - fine right hand side vector on level ILEV
*                  ALPHA - relaxation parameter according to some
*                          optimization criterion but within the limits
*                          AMINP and AMAXP (COMMON /RPARM/)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DD(*),DB(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
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
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C
C *** Setting of given ALPHA
      IF (AMINP.EQ.AMAXP) THEN
       ALPHA=AMINP
       RETURN
      ENDIF
C
C
C *** Calculation of optimal ALPHA
      IF (AMINP.NE.AMAXP) THEN
C
       IF (ISORTP.GT.0) THEN
        KTRC1=L(KLTRC1(ILEV))
        KTRC2=L(KLTRC2(ILEV))
        CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
        CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
        CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
        CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
        CALL VECSRT(DD ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
        CALL LCP1  (DWORK(L(LDP)),DD ,NEQ)
       ENDIF
C
       CALL DYAX7(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),NEQ,DD,DX,DD,DYAX1,DYAX2)
       CALL LSP1(DB,DD,NEQ,DBX)
       ALPHA=(DBX-DYAX1)/DYAX2
       IF (ALPHA.LT.AMINP) ALPHA=AMINP
       IF (ALPHA.GT.AMAXP) ALPHA=AMAXP
C
       IF (ISORTP.GT.0) THEN
        KTRC1=L(KLTRC1(ILEV))
        KTRC2=L(KLTRC2(ILEV))
        CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
        CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
        CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
        CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
        CALL VECSRT(DD ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
        CALL LCP1  (DWORK(L(LDP)),DD ,NEQ)
       ENDIF
C
       RETURN
C
      ENDIF
C
C
      END
