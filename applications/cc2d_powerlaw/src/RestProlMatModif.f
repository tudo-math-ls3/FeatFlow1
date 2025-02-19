************************************************************************
      SUBROUTINE RESTU2 (DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,KADJ2,KADJ1,
     *                   DCORVG,VAREA1,ASPRAT,
     *                   NVT2,NVT1,NEL2,NEL1,NMT1,NMT2)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KMID1(NNVE,*),
     *          KADJ1(NNVE,*),KVERT2(NNVE,*),KMID2(NNVE,*),
     *          KADJ2(NNVE,*),DCORVG(2,*),VAREA1(*), ASPRAT(*)
C
C
C
      CALL  LCL1 (DU1,NMT1)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELA1=KADJ1(1,IEL1)
      IELA2=KADJ1(2,IEL1)
      IELA3=KADJ1(3,IEL1)
      IELA4=KADJ1(4,IEL1)
C
C
      CALL ARCALC(ARIEL ,IEL1 ,KVERT1,DCORVG)
c$$$      ARIEL=ASPRAT(iel1)
      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL
      DAREA=DBLE(VAREA1(IEL1))
C
C
      IF (IELA1.NE.0) THEN
       CALL ARCALC(ARIEL1,IELA1,KVERT1,DCORVG)
c$$$      ARIEL1=ASPRAT(iela1)
       IF (ARIEL1.LT.1D0) ARIEL1=1D0/ARIEL1
       DAREA1=DBLE(VAREA1(IELA1))
      ELSE
       ARIEL1=0D0
       DAREA1=0D0
      ENDIF
C
      IF (IELA2.NE.0) THEN
       CALL ARCALC(ARIEL2,IELA2,KVERT1,DCORVG)
c$$$      ARIEL2=ASPRAT(iela2)
       IF (ARIEL2.LT.1D0) ARIEL2=1D0/ARIEL2
       DAREA2=DBLE(VAREA1(IELA2))
      ELSE
       ARIEL2=0D0
       DAREA2=0D0
      ENDIF
C
      IF (IELA3.NE.0) THEN
       CALL ARCALC(ARIEL3,IELA3,KVERT1,DCORVG)
c$$$      ARIEL3=ASPRAT(iela3)
       IF (ARIEL3.LT.1D0) ARIEL3=1D0/ARIEL3
       DAREA3=DBLE(VAREA1(IELA3))
      ELSE
       ARIEL3=0D0
       DAREA3=0D0
      ENDIF
C
      IF (IELA4.NE.0) THEN
       CALL ARCALC(ARIEL4,IELA4,KVERT1,DCORVG)
c$$$      ARIEL4=ASPRAT(iela4)
       IF (ARIEL4.LT.1D0) ARIEL4=1D0/ARIEL4
       DAREA4=DBLE(VAREA1(IELA4))
      ELSE
       ARIEL4=0D0
       DAREA4=0D0
      ENDIF
C
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
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA1.NE.0D0) THEN
       DQAREA=DAREA/DAREA1
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINTU).EQ.1).AND.(ARIEL.GT.PSBARR)) IGRID=1
      IF ((ABS(IINTU).EQ.2).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL1.GT.PSBARR))) IGRID=1
      IF ((ABS(IINTU).EQ.3).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL1.GT.PSBARR).OR.
     *     (DQAREA.GT.PSBVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       R1= 0.75D0
       R2=-0.25D0
      ELSE
       R1=0.5D0
       R2=0.0D0
      ENDIF
C
      WEIGHT=0.5D0*DAREA/(DAREA+DAREA1)
      DU1(IM1)=DU1(IM1)+WEIGHT*( R1*(DUH1+DUH2)+2D0*R1*DUH9
     *                          +R2*(DUH8+DUH3+DUH10+DUH12))
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA2.NE.0D0) THEN
       DQAREA=DAREA/DAREA2
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINTU).EQ.1).AND.(ARIEL.GT.PSBARR)) IGRID=1
      IF ((ABS(IINTU).EQ.2).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL2.GT.PSBARR))) IGRID=1
      IF ((ABS(IINTU).EQ.3).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL2.GT.PSBARR).OR.
     *     (DQAREA.GT.PSBVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       R1= 0.75D0
       R2=-0.25D0
      ELSE
       R1=0.5D0
       R2=0.0D0
      ENDIF
C
      WEIGHT=0.5D0*DAREA/(DAREA+DAREA2)
      DU1(IM2)=DU1(IM2)+WEIGHT*( R1*(DUH3+DUH4)+2D0*R1*DUH10  
     *                          +R2*(DUH2+DUH5+DUH9+DUH11))
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA3.NE.0D0) THEN
       DQAREA=DAREA/DAREA3
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINTU).EQ.1).AND.(ARIEL.GT.PSBARR)) IGRID=1
      IF ((ABS(IINTU).EQ.2).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL3.GT.PSBARR))) IGRID=1
      IF ((ABS(IINTU).EQ.3).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL3.GT.PSBARR).OR.
     *     (DQAREA.GT.PSBVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       R1= 0.75D0
       R2=-0.25D0
      ELSE
       R1=0.5D0
       R2=0.0D0
      ENDIF
C
      WEIGHT=0.5D0*DAREA/(DAREA+DAREA3)
      DU1(IM3)=DU1(IM3)+WEIGHT*( R1*(DUH5+DUH6)+2D0*R1*DUH11  
     *                          +R2*(DUH4+DUH7+DUH10+DUH12))
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA4.NE.0D0) THEN
       DQAREA=DAREA/DAREA4
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINTU).EQ.1).AND.(ARIEL.GT.PSBARR)) IGRID=1
      IF ((ABS(IINTU).EQ.2).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL4.GT.PSBARR))) IGRID=1
      IF ((ABS(IINTU).EQ.3).AND.
     *    ((ARIEL.GT.PSBARR).OR.(ARIEL4.GT.PSBARR).OR.
     *     (DQAREA.GT.PSBVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       R1= 0.75D0
       R2=-0.25D0
      ELSE
       R1=0.5D0
       R2=0.0D0
      ENDIF
C
      WEIGHT=0.5D0*DAREA/(DAREA+DAREA4)
      DU1(IM4)=DU1(IM4)+WEIGHT*( R1*(DUH7+DUH8)+2D0*R1*DUH12  
     *                          +R2*(DUH6+DUH1+DUH9+DUH11))
C
C
10    CONTINUE
C
C
      END
c
************************************************************************
      SUBROUTINE YREST2(DDF,DDC)  
************************************************************************
C
C-----------------------------------------------------------------------
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*            - DDF and DDC have the structure  DD=(DD1,DD2,DDP)
C  Entstanden aus "YRESTA" von cp2d
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      DIMENSION DDF(*),DDC(*)
C
C
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
      NMT2=KNMT(ILEV+1)
      NMT1=KNMT(ILEV)
C
      LLAR2=L(KLAREA(ILEV+1))
      LLAR1=L(KLAREA(ILEV))
C
      I2C=  1+NMT1
      IPC=I2C+NMT1
      I2F=  1+NMT2
      IPF=I2F+NMT2
C
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.0) THEN
       CALL MRUC (DDF(1),DDC(1),KWORK(LLV2),KWORK(LLV1),
     *            KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *            VWORK(LLAR1),NVT2,NVT1,NEL2,NEL1,NMT1,IINTU)
       CALL MRUC (DDF(I2F),DDC(I2C),KWORK(LLV2),KWORK(LLV1),
     *            KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *            VWORK(LLAR1),NVT2,NVT1,NEL2,NEL1,NMT1,IINTU)
      ELSE
       CALL MRU (DDF(1),DDC(1),KWORK(LLV2),KWORK(LLV1),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           DWORK(L(LCORVG)),VWORK(LLAR1),
     *           NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,
     *           IINTU,PSBARR,PSBVOL)
       CALL MRU (DDF(I2F),DDC(I2C),KWORK(LLV2),KWORK(LLV1),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           DWORK(L(LCORVG)),VWORK(LLAR1),
     *           NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,
     *           IINTU,PSBARR,PSBVOL)
      ENDIF
C

C-----------------------------------------------------------------------
C
      IF (IINTP.EQ.0) THEN
       CALL MR010(DDC(IPC),DDF(IPF),KWORK(LLA1),KWORK(LLA2),NEL1,NEL2)
      ELSEIF (ABS(IINTP).EQ.1) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
       CALL C2N2DM(DDF(IPF),DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *             NEL2,NMT2,NVT2,0)
       CALL MRUC (DWORK(KPLF),DWORK(KPLC),KWORK(LLV2),KWORK(LLV1),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *            VWORK(LLAR1),NVT2,NVT1,NEL2,NEL1,NMT1,2)
       CALL C2N2DM(DDC(IPC),DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *             NEL1,NMT1,NVT1,1)
      ELSE
       KPLC=L(LD1)
       KPLF=L(LD2)
       CALL C2N2DM(DDF(IPF),DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *             NEL2,NMT2,NVT2,0)
       CALL MRU (DWORK(KPLF),DWORK(KPLC),KWORK(LLV2),KWORK(LLV1),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           DWORK(L(LCORVG)),VWORK(LLAR1),
     *           NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,
     *           IINTP,PSBARR,PSBVOL)
       CALL C2N2DM(DDC(IPC),DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *             NEL1,NMT1,NVT1,1)
      ENDIF
C
c      write (*,*) (i,ddc(ipc-1+i),'\n',i=i,10)
C-----------------------------------------------------------------------
C
99999 END
c
************************************************************************
      SUBROUTINE YPROL2 (DUC,DUF)  
************************************************************************
C
C-----------------------------------------------------------------------
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
*            - DDF and DDC have the structure  DD=(DD1,DD2,DDP)
C  entstanden aus "YPROLA" aus cp2d
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
CC
      DIMENSION DUF(*),DUC(*)
C
C
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
      LLAR1=L(KLAREA(ILEV-1))
      LLAR2=L(KLAREA(ILEV))
C
      I2C=  1+NMT1
      IPC=I2C+NMT1
      I2F=  1+NMT2
      IPF=I2F+NMT2
C
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.0) THEN
       CALL MPUC (DUC(1),DUF(1),KWORK(LLV1),KWORK(LLV2),
     *            KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *            VWORK(LLAR1),NVT1,NVT2,NEL1,NEL2,NMT2,IINTU)
       CALL MPUC (DUC(I2C),DUF(I2F),KWORK(LLV1),KWORK(LLV2),
     *            KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *            VWORK(LLAR1),NVT1,NVT2,NEL1,NEL2,NMT2,IINTU)
      ELSE
       CALL MPU (DUC(1),DUF(1),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           DWORK(L(LCORVG)),VWORK(LLAR1),
     *           NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *           IINTU,PSBARP,PSBVOL)
       CALL MPU (DUC(I2C),DUF(I2F),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           DWORK(L(LCORVG)),VWORK(LLAR1),
     *           NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *           IINTU,PSBARP,PSBVOL)
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF (IINTP.EQ.0) THEN
       CALL MP010(DUC(IPC),DUF(IPF),KWORK(LLA1),KWORK(LLA2),NEL1,NEL2)
      ELSEIF (abs(IINTP).eq.1) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
       CALL C2N2DM(DUC(IPC),DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *             NEL1,NMT1,NVT1,0)
       CALL MPUC (DWORK(KPLC),DWORK(KPLF),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *            VWORK(LLAR1),NVT1,NVT2,NEL1,NEL2,NMT2,2)
       CALL C2N2DM(DUF(IPF),DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *             NEL2,NMT2,NVT2,1)
      ELSE
       KPLC=L(LD1)
       KPLF=L(LD2)
       CALL C2N2DM(DUC(IPC),DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *             NEL1,NMT1,NVT1,0)
       CALL MPU (DWORK(KPLC),DWORK(KPLF),KWORK(LLV1),KWORK(LLV2),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           DWORK(L(LCORVG)),VWORK(LLAR1),
     *           NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *           IINTP,PSBARP,PSBVOL)
       CALL C2N2DM(DUF(IPF),DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *             NEL2,NMT2,NVT2,1)
      ENDIF
C
C-----------------------------------------------------------------------
C
99999 END
C
************************************************************************
      SUBROUTINE MPU(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,KADJ1,KADJ2,
     *               DCORVG,VAREA1,NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *               IINT1,PARP,PVOL)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KVERT2(NNVE,*),
     *          KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),KADJ2(NNVE,*),
     *          DCORVG(2,*),VAREA1(*)
C
C
C
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
      IELA1=KADJ1(1,IEL1)
      IELA2=KADJ1(2,IEL1)
      IELA3=KADJ1(3,IEL1)
      IELA4=KADJ1(4,IEL1)
C
C
      CALL ARCALC(ARIEL ,IEL1 ,KVERT1,DCORVG)
      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL
      DAREA=DBLE(VAREA1(IEL1))
C
C
      IF (IELA1.NE.0) THEN
       CALL ARCALC(ARIEL1,IELA1,KVERT1,DCORVG)
       IF (ARIEL1.LT.1D0) ARIEL1=1D0/ARIEL1
       DAREA1=DBLE(VAREA1(IELA1))
      ELSE
       ARIEL1=0D0
       DAREA1=0D0
      ENDIF
C
      IF (IELA2.NE.0) THEN
       CALL ARCALC(ARIEL2,IELA2,KVERT1,DCORVG)
       IF (ARIEL2.LT.1D0) ARIEL2=1D0/ARIEL2
       DAREA2=DBLE(VAREA1(IELA2))
      ELSE
       ARIEL2=0D0
       DAREA2=0D0
      ENDIF
C
      IF (IELA3.NE.0) THEN
       CALL ARCALC(ARIEL3,IELA3,KVERT1,DCORVG)
       IF (ARIEL3.LT.1D0) ARIEL3=1D0/ARIEL3
       DAREA3=DBLE(VAREA1(IELA3))
      ELSE
       ARIEL3=0D0
       DAREA3=0D0
      ENDIF
C
      IF (IELA4.NE.0) THEN
       CALL ARCALC(ARIEL4,IELA4,KVERT1,DCORVG)
       IF (ARIEL4.LT.1D0) ARIEL4=1D0/ARIEL4
       DAREA4=DBLE(VAREA1(IELA4))
      ELSE
       ARIEL4=0D0
       DAREA4=0D0
      ENDIF
C
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
C
      IA=KMID2(1,IELH1)-NVT2
      IB=KMID2(4,IELH2)-NVT2
      IC=KMID2(2,IELH1)-NVT2
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA1.NE.0D0) THEN
       DQAREA=DAREA/DAREA1
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARP)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL1.GT.PARP))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL1.GT.PARP).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA1)
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C
C
      IA=KMID2(1,IELH2)-NVT2
      IB=KMID2(4,IELH3)-NVT2
      IC=KMID2(2,IELH2)-NVT2
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA2.NE.0D0) THEN
       DQAREA=DAREA/DAREA2
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARP)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL2.GT.PARP))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL2.GT.PARP).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA2)
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C
C
      IA=KMID2(1,IELH3)-NVT2
      IB=KMID2(4,IELH4)-NVT2
      IC=KMID2(2,IELH3)-NVT2
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA3.NE.0D0) THEN
       DQAREA=DAREA/DAREA3
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARP)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL3.GT.PARP))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL3.GT.PARP).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA3)
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C
C
      IA=KMID2(1,IELH4)-NVT2
      IB=KMID2(4,IELH1)-NVT2
      IC=KMID2(2,IELH4)-NVT2
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA4.NE.0D0) THEN
       DQAREA=DAREA/DAREA4
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARP)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL4.GT.PARP))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARP).OR.(ARIEL4.GT.PARP).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA4)
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
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
      SUBROUTINE MRU(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,KADJ2,KADJ1,
     *               DCORVG,VAREA1,NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,
     *               IINT1,PARR,PVOL)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'

C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      DIMENSION DU1(1),DU2(1),KVERT2(NNVE,1),KVERT1(NNVE,1),
     *          KMID1(NNVE,1),KMID2(NNVE,1),KADJ1(NNVE,1),KADJ2(NNVE,1),
     *          DCORVG(2,*),VAREA1(*)
C
C
      CALL  LCL1 (DU1,NMT1)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELA1=KADJ1(1,IEL1)
      IELA2=KADJ1(2,IEL1)
      IELA3=KADJ1(3,IEL1)
      IELA4=KADJ1(4,IEL1)
C
C
      CALL ARCALC(ARIEL ,IEL1 ,KVERT1,DCORVG)
      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL
      DAREA=DBLE(VAREA1(IEL1))
C
C
      IF (IELA1.NE.0) THEN
       CALL ARCALC(ARIEL1,IELA1,KVERT1,DCORVG)
       IF (ARIEL1.LT.1D0) ARIEL1=1D0/ARIEL1
       DAREA1=DBLE(VAREA1(IELA1))
      ELSE
       ARIEL1=0D0
       DAREA1=0D0
      ENDIF
C
      IF (IELA2.NE.0) THEN
       CALL ARCALC(ARIEL2,IELA2,KVERT1,DCORVG)
       IF (ARIEL2.LT.1D0) ARIEL2=1D0/ARIEL2
       DAREA2=DBLE(VAREA1(IELA2))
      ELSE
       ARIEL2=0D0
       DAREA2=0D0
      ENDIF
C
      IF (IELA3.NE.0) THEN
       CALL ARCALC(ARIEL3,IELA3,KVERT1,DCORVG)
       IF (ARIEL3.LT.1D0) ARIEL3=1D0/ARIEL3
       DAREA3=DBLE(VAREA1(IELA3))
      ELSE
       ARIEL3=0D0
       DAREA3=0D0
      ENDIF
C
      IF (IELA4.NE.0) THEN
       CALL ARCALC(ARIEL4,IELA4,KVERT1,DCORVG)
       IF (ARIEL4.LT.1D0) ARIEL4=1D0/ARIEL4
       DAREA4=DBLE(VAREA1(IELA4))
      ELSE
       ARIEL4=0D0
       DAREA4=0D0
      ENDIF
C
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
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA1.NE.0D0) THEN
       DQAREA=DAREA/DAREA1
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARR)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL1.GT.PARR))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL1.GT.PARR).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA1)
      DU1(IM1)=DU1(IM1)+WEIGHT*( A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                          +A3*(DUH5+DUH6)+A4*(DUH3+DUH8))
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA2.NE.0D0) THEN
       DQAREA=DAREA/DAREA2
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARR)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL2.GT.PARR))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL2.GT.PARR).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA2)
      DU1(IM2)=DU1(IM2)+WEIGHT*( A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                          +A3*(DUH7+DUH8)+A4*(DUH5+DUH2))
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA3.NE.0D0) THEN
       DQAREA=DAREA/DAREA3
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARR)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL3.GT.PARR))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL3.GT.PARR).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA3)
      DU1(IM3)=DU1(IM3)+WEIGHT*( A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                          +A3*(DUH1+DUH2)+A4*(DUH7+DUH4))
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
C
C
      IGRID=0
      DQAREA=0D0
      IF (DAREA4.NE.0D0) THEN
       DQAREA=DAREA/DAREA4
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
      ENDIF
      IF ((ABS(IINT1).EQ.1).AND.(ARIEL.GT.PARR)) IGRID=1
      IF ((ABS(IINT1).EQ.2).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL4.GT.PARR))) IGRID=1
      IF ((ABS(IINT1).EQ.3).AND.
     *    ((ARIEL.GT.PARR).OR.(ARIEL4.GT.PARR).OR.
     *     (DQAREA.GT.PVOL))) IGRID=1
C
      IF (IGRID.EQ.0) THEN
       IF (IINT1.LT.0) THEN
        A1=1D0
        A2=-0.125D0
        A3=0D0
        A4=0.125D0
        A5=0.625D0
        A6=0.125D0
        A7=0.125D0
        A8=0.125D0
       ELSE
        A1=0.9375D0
        A2=-0.1875D0
        A3=-0.0625D0
        A4=0.3125D0
        A5=0.5625D0
        A6=0.1875D0
        A7=0.0625D0
        A8=0.1875D0
       ENDIF
      ELSE
       A1=1D0
       A2=0D0
       A3=0D0
       A4=0D0
       A5=1D0
       A6=0D0
       A7=0D0
       A8=0D0
      ENDIF
C
      WEIGHT=DAREA/(DAREA+DAREA4)
      DU1(IM4)=DU1(IM4)+WEIGHT*( A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                          +A3*(DUH3+DUH4)+A4*(DUH1+DUH6))
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
C
C
10    CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE MRUC (DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,KADJ2,KADJ1,
     *                 VAREA1,NVT2,NVT1,NEL2,NEL1,NMT1,IINT1)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:  KONSTANTE Restriktion , Mittelung der Nachbarknoten
C           Wird aufgerufen, falls IINTU=0
C           Warum die Faelle IINTU<> 0 abgefragt werden
c           weiss ich nicht
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'

C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      PARAMETER (A1=1D0,A2=0D0,A3=0D0,A4=0D0)
      PARAMETER (A5=1D0,A6=0D0,A7=0D0,A8=0D0)
      DIMENSION DU1(1),DU2(1),KVERT2(NNVE,1),KVERT1(NNVE,1),
     *          KMID1(NNVE,1),KMID2(NNVE,1),KADJ1(NNVE,1),KADJ2(NNVE,1),
     *          VAREA1(*)
C
C
C
c      write (*,*) 'in mruc, iintu= ', iintu, iintp, iint1

      CALL  LCL1 (DU1,NMT1)
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
      DAREA=DBLE(VAREA1(IEL1))
C
C
      IADJ1=KADJ1(1,IEL1)
      IF (IADJ1.NE.0) THEN
       DAREA1=DBLE(VAREA1(IADJ1))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA1)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA1/(DAREA+DAREA1)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU1(IM1)=DU1(IM1)
     *         +WEIGHT*( A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                  +A3*(DUH5+DUH6)+A4*(DUH3+DUH8))
     *         +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
C
C
      IADJ2=KADJ1(2,IEL1)
      IF (IADJ2.NE.0) THEN 
       DAREA2=DBLE(VAREA1(IADJ2))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA2)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA2/(DAREA+DAREA2)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU1(IM2)=DU1(IM2)
     *         +WEIGHT*( A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                  +A3*(DUH7+DUH8)+A4*(DUH5+DUH2))
     *         +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
C
C
      IADJ3=KADJ1(3,IEL1)
      IF (IADJ3.NE.0) THEN 
       DAREA3=DBLE(VAREA1(IADJ3))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA3)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA3/(DAREA+DAREA3)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU1(IM3)=DU1(IM3)
     *         +WEIGHT*( A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                  +A3*(DUH1+DUH2)+A4*(DUH7+DUH4))
     *         +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
C
C
      IADJ4=KADJ1(4,IEL1)
      IF (IADJ4.NE.0) THEN 
       DAREA4=DBLE(VAREA1(IADJ4))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA4)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA4/(DAREA+DAREA4)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU1(IM4)=DU1(IM4)
     *         +WEIGHT*( A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                  +A3*(DUH3+DUH4)+A4*(DUH1+DUH6))
     *         +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
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
      SUBROUTINE MPUC (DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,KADJ1,KADJ2,
     *                 VAREA1,NVT1,NVT2,NEL1,NEL2,NMT2,IINT1)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      PARAMETER (A1=1D0,A2=0D0,A3=0D0,A4=0D0)
      PARAMETER (A5=1D0,A6=0D0,A7=0D0,A8=0D0)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KVERT2(NNVE,*),
     *          KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),KADJ2(NNVE,*),
     *          VAREA1(*)
C
C
C
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
c
c$$$      DUH1=1d0!DU1(IM1)
c$$$      DUH2=1d0!DU1(IM2)
c$$$      DUH3=1d0!DU1(IM3)
c$$$      DUH4=1d0!DU1(IM4)
c
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      DAREA=DBLE(VAREA1(IEL1))
C
C
      IA=KMID2(1,IELH1)-NVT2
      IB=KMID2(4,IELH2)-NVT2
      IC=KMID2(2,IELH1)-NVT2
      IADJ1=KADJ1(1,IEL1)
      IF (IADJ1.NE.0) THEN 
       DAREA1=DBLE(VAREA1(IADJ1))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA1)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA1/(DAREA+DAREA1)
c       if (weight.ne.0.5d0) write (*,*) 'weight', weight
      ELSE
       WEIGHT=1D0
      ENDIF
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C
      IA=KMID2(1,IELH2)-NVT2
      IB=KMID2(4,IELH3)-NVT2
      IC=KMID2(2,IELH2)-NVT2
      IADJ2=KADJ1(2,IEL1)
      IF (IADJ2.NE.0) THEN 
       DAREA2=DBLE(VAREA1(IADJ2))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA2)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA2/(DAREA+DAREA2)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C
      IA=KMID2(1,IELH3)-NVT2
      IB=KMID2(4,IELH4)-NVT2
      IC=KMID2(2,IELH3)-NVT2
      IADJ3=KADJ1(3,IEL1)
      IF (IADJ3.NE.0) THEN 
       DAREA3=DBLE(VAREA1(IADJ3))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA3)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA3/(DAREA+DAREA3)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C
      IA=KMID2(1,IELH4)-NVT2
      IB=KMID2(4,IELH1)-NVT2
      IC=KMID2(2,IELH4)-NVT2
      IADJ4=KADJ1(4,IEL1)
      IF (IADJ4.NE.0) THEN 
       DAREA4=DBLE(VAREA1(IADJ4))
       IF (ABS(IINT1).LE.1) WEIGHT=0.5D0
       IF (ABS(IINT1).EQ.2) WEIGHT=DAREA /(DAREA+DAREA4)
       IF (ABS(IINT1).EQ.3) WEIGHT=DAREA4/(DAREA+DAREA4)
      ELSE
       WEIGHT=1D0
      ENDIF
      DU2(IA)=DU2(IA)+WEIGHT*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
      DU2(IB)=DU2(IB)+WEIGHT*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
C
   
C
10    CONTINUE
C
C
c      write (*,*) (du2(i),'\n',i=1,nmt2)

      END
C
C
C
************************************************************************
C     Restriktion der Matrix
************************************************************************
      SUBROUTINE MAREST(KVERT1,KVERT2,KMID1,KMID2,KADJ1,KADJ2,VA1,VA2,
     *                  KLDA1,KLDA2,KCOLA1,KCOLA2,DCORVG,VAREA1,
     *                  NEL1,NEL2,NVT1,NVT2)
************************************************************************
C
C-----------------------------------------------------------------------
*   Purpose:  
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      REAL VA1,VA2
C
c      PARAMETER (NNVE=4)
      DIMENSION KVERT1(NNVE,*),KVERT2(NNVE,*),KMID1(NNVE,*),
     *          KMID2(NNVE,*),KADJ1(NNVE,*),KADJ2(NNVE,*)
      DIMENSION VA1(*),VA2(*),KLDA1(*),KLDA2(*),KCOLA1(*),KCOLA2(*)
      DIMENSION DCORVG(2,*),VAREA1(*)
      DIMENSION KIEL(NNVE)
C
C
CCC      CALL LCL2(VA1,NA)
C
      DO 10 IEL=1,NEL1
C
      CALL ARCALC(ARIEL,IEL,KVERT1,DCORVG)
      IF (ARIEL.LT.1D0) ARIEL=1D0/ARIEL
      DAREA=DBLE(VAREA1(IEL))
C
      KIEL(1)=IEL
      KIEL(2)=KADJ2(2,KIEL(1))
      KIEL(3)=KADJ2(2,KIEL(2))
      KIEL(4)=KADJ2(2,KIEL(3))
C
      DO 20 IVE=1,4
      IVE1=IVE
      IVE2=MOD(IVE1,4)+1
      IVE4=MOD(IVE1+2,4)+1
      IADJ1=KADJ1(IVE1,IEL)
      IADJ2=KADJ1(IVE2,IEL)
      IADJ4=KADJ1(IVE4,IEL)
      IEL1=KIEL(IVE1)
      IEL2=KIEL(IVE2)
C
      IMID1=0
      IMID2=0
      IMID3=0
      IMID4=0
      IMID5=0
C
      DVAL1=0D0
      DVAL2=0D0
      DVAL3=0D0
      DVAL4=0D0
      DVAL5=0D0
C
      IF ((IADJ1.LT.IEL).AND.(IADJ1.NE.0)) GOTO 20
C
      IF (IADJ1.NE.0) THEN
       CALL ARCALC(ARIELA,IADJ1,KVERT1,DCORVG)
       IF (ARIELA.LT.1D0) ARIELA=1D0/ARIELA
       DAREAA=DBLE(VAREA1(IADJ1))
       DQAREA=DAREA/DAREAA
       IF (DQAREA.LT.1D0) DQAREA=1D0/DQAREA
ccc       WRITE(6,*) IEL,IVE,IADJ1,ARIEL,ARIELA,DQAREA
       IF ((ARIEL.LT.PSBARM).AND.(ARIELA.LT.PSBARM).AND.
     *     (DQAREA.LT.PSBVOL)) GOTO 20
      ELSE
ccc       WRITE(6,*) IEL,IVE,IADJ1,ARIEL
       IF (ARIEL.LT.PSBARM) GOTO 20
      ENDIF
C
      IMID1 =KMID1(IVE1,IEL)-NVT1
      IMID2 =KMID1(IVE2,IEL)-NVT1
      IMID3 =KMID1(IVE4,IEL)-NVT1
ccc      WRITE(6,*) 'IMID ',IMID1,IMID2,IMID3,IEL1,IEL2
C
      IM1 =KMID2(1,IEL1)-NVT2
      IM2 =KMID2(4,IEL2)-NVT2
      IM3 =KMID2(2,IEL1)-NVT2
      IM4 =KMID2(4,IEL1)-NVT2
      IM5 =KMID2(1,IEL2)-NVT2
      IM6 =KMID2(3,IEL1)-NVT2
      IM7 =KMID2(2,IEL2)-NVT2
ccc      WRITE(6,*) 'IM ',IM1,IM2,IM3,IM4,IM5,IM6,IM7

C
      IF (IADJ1.NE.0) THEN
       DO 30 JVE=1,4
30     IF (KADJ1(JVE,IADJ1).EQ.IEL) GOTO 32
C
32     JVE1=JVE
       JVE2=MOD(JVE1,4)+1
       JVE4=MOD(JVE1+2,4)+1
       JADJ2=KADJ1(JVE2,IADJ1)
       JADJ4=KADJ1(JVE4,IADJ1)
       JEL1=KADJ2(4,IEL2)
       JEL2=KADJ2(1,IEL1)
C
       IMID4 =KMID1(JVE2,IADJ1)-NVT1
       IMID5 =KMID1(JVE4,IADJ1)-NVT1
C
       IM8 =KMID2(2,JEL1)-NVT2
       IM9 =KMID2(4,JEL1)-NVT2
       IM10=KMID2(1,JEL2)-NVT2
       IM11=KMID2(3,JEL1)-NVT2
       IM12=KMID2(2,JEL2)-NVT2
C
      ELSE
C
       IM8 =0
C
      ENDIF
C
C
ccc      write(6,*) ILD,IM1,IM3,KLDA2(IM1),ICOL,IEL
      DO 101 ILD=KLDA2(IM3)+1,KLDA2(IM3+1)-1
      ICOL=KCOLA2(ILD)
ccc      write(6,*) ILD,IM1,IM3,KLDA2(IM1),ICOL,IEL
101   IF (ICOL.EQ.IM1) GOTO 102
      WRITE(6,*) 'WRONG 101'
C
102   PV1=DBLE(VA2(KLDA2(IM1))+VA2(ILD))
C
      IF (IM8.NE.0) THEN
       DO 104 ILD=KLDA2(IM8)+1,KLDA2(IM8+1)-1
       ICOL=KCOLA2(ILD)
104    IF (ICOL.EQ.IM1) GOTO 105
       WRITE(6,*) 'WRONG 104'
C
105    PV1=PV1+DBLE(VA2(ILD))
      ENDIF
C
C
      DO 111 ILD=KLDA2(IM3)+1,KLDA2(IM3+1)-1
      ICOL=KCOLA2(ILD)
111   IF (ICOL.EQ.IM2) GOTO 112
      WRITE(6,*) 'WRONG 111'
C
112   PV2=DBLE(VA2(KLDA2(IM2))+VA2(ILD))
C
      IF (IM8.NE.0) THEN
       DO 114 ILD=KLDA2(IM8)+1,KLDA2(IM8+1)-1
       ICOL=KCOLA2(ILD)
114    IF (ICOL.EQ.IM2) GOTO 115
       WRITE(6,*) 'WRONG 114'
C
115    PV2=PV2+DBLE(VA2(ILD))
      ENDIF
C
C
      DO 121 ILD=KLDA2(IM1)+1,KLDA2(IM1+1)-1
      ICOL=KCOLA2(ILD)
121   IF (ICOL.EQ.IM3) GOTO 122
      WRITE(6,*) 'WRONG 121'
C
122   PV3=DBLE(VA2(KLDA2(IM3))+VA2(ILD))
C
      DO 124 ILD=KLDA2(IM2)+1,KLDA2(IM2+1)-1
      ICOL=KCOLA2(ILD)
124   IF (ICOL.EQ.IM3) GOTO 125
      WRITE(6,*) 'WRONG 124'
C
125   PV3=PV3+DBLE(VA2(ILD))
C
C
      DO 131 ILD=KLDA2(IM1)+1,KLDA2(IM1+1)-1
      ICOL=KCOLA2(ILD)
131   IF (ICOL.EQ.IM4) GOTO 132
      WRITE(6,*) 'WRONG 131'
C
132   PV4=DBLE(VA2(ILD))
C
      DO 134 ILD=KLDA2(IM3)+1,KLDA2(IM3+1)-1
      ICOL=KCOLA2(ILD)
134   IF (ICOL.EQ.IM4) GOTO 135
      WRITE(6,*) 'WRONG 134'
C
135   PV4=PV4+DBLE(VA2(ILD))
C
C
      DO 141 ILD=KLDA2(IM2)+1,KLDA2(IM2+1)-1
      ICOL=KCOLA2(ILD)
141   IF (ICOL.EQ.IM5) GOTO 142
      WRITE(6,*) 'WRONG 141'
C
142   PV5=DBLE(VA2(ILD))
C
      DO 144 ILD=KLDA2(IM3)+1,KLDA2(IM3+1)-1
      ICOL=KCOLA2(ILD)
144   IF (ICOL.EQ.IM5) GOTO 145
      WRITE(6,*) 'WRONG 144'
C
145   PV5=PV5+DBLE(VA2(ILD))
C
C
      DO 151 ILD=KLDA2(IM1)+1,KLDA2(IM1+1)-1
      ICOL=KCOLA2(ILD)
151   IF (ICOL.EQ.IM6) GOTO 152
      WRITE(6,*) 'WRONG 151'
C
152   PV6=DBLE(VA2(ILD))
C
      DO 154 ILD=KLDA2(IM3)+1,KLDA2(IM3+1)-1
      ICOL=KCOLA2(ILD)
154   IF (ICOL.EQ.IM6) GOTO 155
      WRITE(6,*) 'WRONG 154'
C
155   PV6=PV6+DBLE(VA2(ILD))
C
C
      DO 161 ILD=KLDA2(IM2)+1,KLDA2(IM2+1)-1
      ICOL=KCOLA2(ILD)
161   IF (ICOL.EQ.IM7) GOTO 162
      WRITE(6,*) 'WRONG 161'
C
162   PV7=DBLE(VA2(ILD))
C
      DO 164 ILD=KLDA2(IM3)+1,KLDA2(IM3+1)-1
      ICOL=KCOLA2(ILD)
164   IF (ICOL.EQ.IM7) GOTO 165
      WRITE(6,*) 'WRONG 164'
C
165   PV7=PV7+DBLE(VA2(ILD))
C
C
      IF (IM8.NE.0) THEN
C
       DO 171 ILD=KLDA2(IM1)+1,KLDA2(IM1+1)-1
       ICOL=KCOLA2(ILD)
171    IF (ICOL.EQ.IM8) GOTO 172
       WRITE(6,*) 'WRONG 171'
C
172    PV8=DBLE(VA2(KLDA2(IM8))+VA2(ILD))
C
       DO 174 ILD=KLDA2(IM2)+1,KLDA2(IM2+1)-1
       ICOL=KCOLA2(ILD)
174    IF (ICOL.EQ.IM8) GOTO 175
       WRITE(6,*) 'WRONG 174'
C
175    PV8=PV8+DBLE(VA2(ILD))
C
C
       DO 181 ILD=KLDA2(IM8)+1,KLDA2(IM8+1)-1
       ICOL=KCOLA2(ILD)
181    IF (ICOL.EQ.IM9) GOTO 182
       WRITE(6,*) 'WRONG 181'
C
182    PV9=DBLE(VA2(ILD))
C
       DO 184 ILD=KLDA2(IM2)+1,KLDA2(IM2+1)-1
       ICOL=KCOLA2(ILD)
184    IF (ICOL.EQ.IM9) GOTO 185
       WRITE(6,*) 'WRONG 184'
C
185    PV9=PV9+DBLE(VA2(ILD))
C
C
       DO 191 ILD=KLDA2(IM8)+1,KLDA2(IM8+1)-1
       ICOL=KCOLA2(ILD)
191    IF (ICOL.EQ.IM10) GOTO 192
       WRITE(6,*) 'WRONG 191'
C
192    PV10=DBLE(VA2(ILD))
C
       DO 194 ILD=KLDA2(IM1)+1,KLDA2(IM1+1)-1
       ICOL=KCOLA2(ILD)
194    IF (ICOL.EQ.IM10) GOTO 195
       WRITE(6,*) 'WRONG 194'
C
195    PV10=PV10+DBLE(VA2(ILD))
C
C
       DO 201 ILD=KLDA2(IM8)+1,KLDA2(IM8+1)-1
       ICOL=KCOLA2(ILD)
201    IF (ICOL.EQ.IM11) GOTO 202
       WRITE(6,*) 'WRONG 201'
C
202    PV11=DBLE(VA2(ILD))
C
       DO 204 ILD=KLDA2(IM2)+1,KLDA2(IM2+1)-1
       ICOL=KCOLA2(ILD)
204    IF (ICOL.EQ.IM11) GOTO 205
       WRITE(6,*) 'WRONG 204'
C
205    PV11=PV11+DBLE(VA2(ILD))
C
C
       DO 211 ILD=KLDA2(IM8)+1,KLDA2(IM8+1)-1
       ICOL=KCOLA2(ILD)
211    IF (ICOL.EQ.IM12) GOTO 212
       WRITE(6,*) 'WRONG 211'
C
212    PV12=DBLE(VA2(ILD))
C
       DO 214 ILD=KLDA2(IM1)+1,KLDA2(IM1+1)-1
       ICOL=KCOLA2(ILD)
214    IF (ICOL.EQ.IM12) GOTO 215
       WRITE(6,*) 'WRONG 214'
C
215    PV12=PV12+DBLE(VA2(ILD))
C
      ENDIF
C
C
C
      IF (IADJ1.EQ.0) THEN
       DVAL1=1D0/1D0*(PV1+PV2+PV3)
      ELSE
       DVAL1=1D0/1D0*(PV1+PV2+PV3+PV8)
      ENDIF
C
      IF (IADJ2.EQ.0) THEN
       DVAL2=1D0/1D0*(PV5+PV7)
      ELSE
       DVAL2=1D0/1D0*(PV5+PV7)
      ENDIF
C
      IF (IADJ4.EQ.0) THEN
       DVAL3=1D0/1D0*(PV4+PV6)
      ELSE
       DVAL3=1D0/1D0*(PV4+PV6)
      ENDIF
C
      IF (IADJ1.NE.0) THEN
C
       IF (JADJ2.EQ.0) THEN
        DVAL4=1D0/1D0*(PV10+PV12)
       ELSE
        DVAL4=1D0/1D0*(PV10+PV12)
       ENDIF
C
       IF (JADJ4.EQ.0) THEN
        DVAL5=1D0/1D0*(PV9+PV11)
       ELSE
        DVAL5=1D0/1D0*(PV9+PV11)
       ENDIF
C
      ENDIF
C
C
ccc      WRITE(6,*) IEL,IVE,IMID1,IMID2,IMID3,IMID4,IMID5
ccc      WRITE(6,*) IEL,IVE,DVAL1,DVAL2,DVAL3,DVAL4,DVAL5
C
C
ccc      GOTO 20
C
      DO 1010 ILD=KLDA1(IMID1),KLDA1(IMID1+1)-1
1010  VA1(ILD)=0E0
C
      VA1(KLDA1(IMID1))=REAL(DVAL1)
C
      DO 1020 ILD=KLDA1(IMID1)+1,KLDA1(IMID1+1)-1
      ICOL=KCOLA1(ILD)
1020  IF (ICOL.EQ.IMID2) GOTO 1022
      WRITE(6,*) 'WRONG 1020'
C
1022  VA1(ILD)=DVAL2
C
C
      DO 1030 ILD=KLDA1(IMID1)+1,KLDA1(IMID1+1)-1
      ICOL=KCOLA1(ILD)
1030  IF (ICOL.EQ.IMID3) GOTO 1032
      WRITE(6,*) 'WRONG 1030'
C
1032  VA1(ILD)=DVAL3
C
C
      IF (IADJ1.NE.0) THEN
C
       DO 1040 ILD=KLDA1(IMID1)+1,KLDA1(IMID1+1)-1
       ICOL=KCOLA1(ILD)
1040   IF (ICOL.EQ.IMID4) GOTO 1042
       WRITE(6,*) 'WRONG 1040'
C
1042   VA1(ILD)=DVAL4
C
C
       DO 1050 ILD=KLDA1(IMID1)+1,KLDA1(IMID1+1)-1
       ICOL=KCOLA1(ILD)
1050   IF (ICOL.EQ.IMID5) GOTO 1052
       WRITE(6,*) 'WRONG 1050'
C
1052   VA1(ILD)=DVAL5
C
      ENDIF
C
C
20    CONTINUE
C
10    CONTINUE
C
C
C
C
      END
