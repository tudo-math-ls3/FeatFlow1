************************************************************************
      SUBROUTINE  XOPTCN (KU1,KU2,KU1OLD,KU2OLD,KF1,KF2,KD1,KD2,
     *                    KAUX1,KAUX2,KA1,KCOLA,KLDA,KST1,KM1,KMASS1,
     *                    NA,NU,DELU1,DELU2,OMEGA,KMBD,NMBD,INEUM)
************************************************************************
*    Purpose: - Computes for given vectors U and UOLD the optimal 
*               weighted correction OMEGA*V with V:=U-UOLD
*
*    Input:
*      - vectors U,UOLD,F
*      - OMEGA for the new calculation of the nonlinear block A
*      - matrices and pointers A,...,ST,NA
*      - number of equations NU
*    Output:
*      - updated vector U
*      - optimal relaxation parameter OMEGA
*      - maximum relative changes  DELU 
*      - note that D is changed to  D=K*V with correction V=U-UOLD 
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION KMBD(*)
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
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
      common /locvis/ klny(NNLEV),kny
      SAVE 
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Coefficient of stiffness matrix
      EXTERNAL COEFNN
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30
C
      EXTERNAL MATML
C
C-----------------------------------------------------------------------
C
      IF (OMGMIN.EQ.OMGMAX) THEN
       IF (OMGMIN.LT.0D0) THEN
        OMEGA=ABS(OMGMIN)
        DELU1=0D0
        DELU2=0D0
        GOTO 99999
       ELSE
        OMEGA=OMGMIN
        GOTO 999
       ENDIF
      ENDIF
C
C=======================================================================
C *** Calculate on AUX the linearization point: UOLD+OMEGA*U
C=======================================================================
C
      CALL ZTIME(TTT0)
      AA1=1.0D0
      AA2=OMEGA
      CALL  LCP1 (DWORK(KU1),DWORK(KAUX1),NU)
      CALL  LCP1 (DWORK(KU2),DWORK(KAUX2),NU)
      CALL  LLC1 (DWORK(KU1OLD),DWORK(KAUX1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU2OLD),DWORK(KAUX2),NU,AA1,AA2)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Calculate the new nonlinear block A at the point AUX
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL XMADF3(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      IF (IUPW.EQ.1) THEN
            if (IPRECA.eq.4) then
             IF (IELT.EQ.0) 
     *         CALL ADDSTP(DWORK(KAUX1),DWORK(KAUX2),
     *              DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *         DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *         KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *         DWORK(L(LCORVG)),E031,COEFNN,0,3,1D0,DWORK(kny))
             IF (IELT.EQ.1) 
     *         CALL ADDSTP(DWORK(KAUX1),DWORK(KAUX2),
     *            DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *         DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *         KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *         DWORK(L(LCORVG)),E030,COEFNN,0,3,1D0,DWORK(kny))
             IF (IELT.EQ.2) 
     *         CALL ADDSTN(DWORK(KAUX1),DWORK(KAUX2),
     *            DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *         DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *         KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *         DWORK(L(LCORVG)),EM31,COEFNN,0,3,1D0,DWORK(kny))
             IF (IELT.EQ.3) 
     *         CALL ADDSTN(DWORK(KAUX1),DWORK(KAUX2),
     *            DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *         DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *         KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *         DWORK(L(LCORVG)),EM30,COEFNN,0,3,1D0,DWORK(kny))
            endif
            if(ISTOK.NE.1)
     *  CALL GUPWD (DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),0,3,DWORK(kny))
      ELSE
       IF (IELT.EQ.0) 
     *  CALL SUPWDG(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E031,COEFNN,0,3,1D0,DWORK(kny))
       IF (IELT.EQ.1) 
     *  CALL SUPWDG(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E030,COEFNN,0,3,1D0,DWORK(kny))
       IF (IELT.EQ.2) 
     *  CALL SUPWNP(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM31,COEFNN,0,3,1D0,DWORK(kny))
       IF (IELT.EQ.3) 
     *  CALL SUPWNP(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),VWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM30,COEFNN,0,3,1D0,DWORK(kny))
      ENDIF
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL BDRYA  (VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),KMBD,NMBD)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
C=======================================================================
C *** Calculate the defect  D=F-K*UOLD
C=======================================================================
      CALL ZTIME(TTT0)
      CALL LCP1 (DWORK(KF1),DWORK(KD1),NU)
      CALL LCP1 (DWORK(KF2),DWORK(KD2),NU)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL MATML(DWORK(KD1),DWORK(KD2),DWORK(KU1OLD),DWORK(KU2OLD),
     *     -1D0,1D0,VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C *** Calculate the value  AUX=K*U
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL MATML(DWORK(KAUX1),DWORK(KAUX2),DWORK(KU1),DWORK(KU2),
     *            1D0,0D0,VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C *** Calculate   SKV1:= (K*U,D)   = (AUX,D)
C *** Calculate   SKV2:= (K*U,K*U) = (AUX,AUX)
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL  LSP1 (DWORK(KAUX1),DWORK(KD1)  ,2*NU,SKV1)
      CALL  LSP1 (DWORK(KAUX1),DWORK(KAUX1),2*NU,SKV2)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Calculate the optimal relaxation parameter OMEGA
C=======================================================================
C
      IF (SKV2.LT. 1.0D-40) THEN
          WRITE(MTERM,*) 'ERROR in OPTCOR: SKV2 is nearly zero'
          STOP
      ENDIF
C
      OMEGA=SKV1/SKV2
      IF (OMEGA.LT.ABS(OMGMIN)) OMEGA=ABS(OMGMIN)
      IF (OMEGA.GT.ABS(OMGMAX)) OMEGA=ABS(OMGMAX)
C
C=======================================================================
C *** Calculate the optimal correction  U:=UOLD+OMEGA*U
C=======================================================================
C
999   CALL ZTIME(TTT0)
      CALL  LLI1 (DWORK(KU1),NU,DELU1,INDU1)
      CALL  LLI1 (DWORK(KU2),NU,DELU2,INDU2)
C
      AA1= 1.D0
      AA2= OMEGA
      CALL  LLC1 (DWORK(KU1OLD),DWORK(KU1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU2OLD),DWORK(KU2),NU,AA1,AA2)
C
C=======================================================================
C *** relative maximum changes   DELU
C=======================================================================
C
      CALL  LLI1 (DWORK(KU1),NU,DELT1,INDT1)
      CALL  LLI1 (DWORK(KU2),NU,DELT2,INDT2)
      DELT=MAX(DELT1,DELT2)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELU1=DELU1/DELT
      DELU2=DELU2/DELT
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C
99999 END
