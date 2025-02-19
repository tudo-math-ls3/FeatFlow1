      SUBROUTINE   PROLU (DU1,DV1,DP1,DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,NMT1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NMT2,NEL2,NVT2,
     *                     DCORVG, AREA)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1,DP1) to
*                the fine grid vector (DU2,DV2,DP2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1,DP1           - coarse grid vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*      DCORVG                - Coarse grid vertex coordinates
*      AREA                  - Array with area sizes of elements
*
*    Output:
*      DU2,DV2,DP2           - interpolated fine grid vector
*-----------------------------------------------------------------------
      
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      INCLUDE 'cinidat.inc'
      INCLUDE 'cinidat2.inc'

C parameters

      DOUBLE PRECISION DU1(*),DV1(*),DP1(*),  DU2(*),DV2(*),DP2(*)
      DOUBLE PRECISION DCORVG(2,*)
      DOUBLE PRECISION AREA(*)      
      INTEGER KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*)
      INTEGER KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)

      INTEGER NMT1, NEL1, NVT1
      INTEGER NMT2, NEL2, NVT2
      
C local variables

      INTEGER IADPR1, IAVRT2

      IF ((IINT.EQ.1).OR.(IINT.EQ.2)) THEN
        CALL MP031 (DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *              KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
        CALL MP031 (DV1,DV2,KVERT1,KVERT2,KMID1,KMID2,
     *              KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
     
      ELSE IF ((IINT.EQ.-1).OR.(IINT.EQ.-2)) THEN
        CALL MP030 (DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *              KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
        CALL MP030 (DV1,DV2,KVERT1,KVERT2,KMID1,KMID2,
     *              KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
     
      ELSE IF ((IINT.EQ.3).OR.(IINT.EQ.4)) THEN
C Evaluate bitfield for prolongation/restriction
        IADPR1=0
        IF(IAND(IAPRM,1).EQ.1) IADPR1=1
        IF(IAND(IAPRM,9).EQ.9) IADPR1=2
C Evaluate bitfield for weighting
        IAVRT2=0
        IF(IAND(IAVPR,1).EQ.1) IAVRT2=1
        IF(IAND(IAVPR,5).EQ.5) IAVRT2=2
        CALL MP031X(DU1,DU2,KVERT1,KVERT2,
     *              KMID1,KMID2,KADJ1,KADJ2,
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
        CALL MP031X(DV1,DV2,KVERT1,KVERT2,
     *              KMID1,KMID2,KADJ1,KADJ2,
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
     
      ELSE IF ((IINT.EQ.-3).OR.(IINT.EQ.-4)) THEN
C Evaluate bitfield for prolongation/restriction
        IADPR1=0
        IF(IAND(IAPRM,1).EQ.1) IADPR1=1
        IF(IAND(IAPRM,9).EQ.9) IADPR1=2
C Evaluate bitfield for weighting
        IAVRT2=0
        IF(IAND(IAVPR,1).EQ.1) IAVRT2=1
        IF(IAND(IAVPR,5).EQ.5) IAVRT2=2
        CALL MP030X(DU1,DU2,KVERT1,KVERT2,
     *              KMID1,KMID2,KADJ1,KADJ2,
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
        CALL MP030X(DV1,DV2,KVERT1,KVERT2,
     *              KMID1,KMID2,KADJ1,KADJ2,
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
      END IF
      
      CALL MP010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
     
      END 

      SUBROUTINE  RESTRD (DU1,DV1,DP1,DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,NMT1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NMT2,NEL2,NVT2,
     *                     DCORVG, AREA)
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
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cinidat2.inc'
      
C parameters

      DOUBLE PRECISION DU1(*),DV1(*),DP1(*),  DU2(*),DV2(*),DP2(*)
      DOUBLE PRECISION DCORVG(2,*)
      DOUBLE PRECISION AREA(*)      
      INTEGER KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*)
      INTEGER KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)

      INTEGER NMT1, NEL1, NVT1
      INTEGER NMT2, NEL2, NVT2

C local variables

      INTEGER IADPR1, IAVRT2

      IF ((IINT.EQ.1).OR.(IINT.EQ.2)) THEN
        CALL MR031(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *           KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
        CALL MR031(DV2,DV1,KVERT2,KVERT1,KMID2,KMID1,
     *           KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
     
      ELSE IF ((IINT.EQ.-1).OR.(IINT.EQ.-2)) THEN
        CALL MR030(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *           KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
        CALL MR030(DV2,DV1,KVERT2,KVERT1,KMID2,KMID1,
     *           KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
     
      ELSE IF ((IINT.EQ.3).OR.(IINT.EQ.4)) THEN
C Evaluate Bitfield for prol/rest
        IADPR1=0
        IF(IAND(IAPRM,2).EQ.2) IADPR1=1
        IF(IAND(IAPRM,18).EQ.18) IADPR1=2
C Evaluate Bitfield for weighting
        IAVRT2=0
        IF(IAND(IAVPR,2).EQ.2) IAVRT2=1
        IF(IAND(IAVPR,10).EQ.10) IAVRT2=2
        CALL MR031X(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *              KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1,
     *              NMT2,NMT1,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
        CALL MR031X(DV2,DV1,KVERT2,KVERT1,KMID2,KMID1,
     *              KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1,
     *              NMT2,NMT1,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
     
      ELSE IF ((IINT.EQ.-3).OR.(IINT.EQ.-4)) THEN
C Evaluate Bitfield for prol/rest
        IADPR1=0
        IF(IAND(IAPRM,2).EQ.2) IADPR1=1
        IF(IAND(IAPRM,18).EQ.18) IADPR1=2
C Evaluate Bitfield for weighting
        IAVRT2=0
        IF(IAND(IAVPR,2).EQ.2) IAVRT2=1
        IF(IAND(IAVPR,10).EQ.10) IAVRT2=2
        CALL MR030X(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *              KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1,
     *              NMT2,NMT1,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
        CALL MR030X(DV2,DV1,KVERT2,KVERT1,KMID2,KMID1,
     *              KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1,
     *              NMT2,NMT1,DCORVG,
     *              AREA,IAVRT2,DPREP,IADPR1)
      END IF

      CALL MR010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)

      END

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
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C constants
      
      DOUBLE PRECISION A1, A2, A3, R1, R2, R3

      PARAMETER (A1=0.1875D0, A2=0.375D0, A3=-0.0625D0)
      PARAMETER (R1=0.375D0, R2=0.75D0, R3=-0.125D0)

C parameters

      DOUBLE PRECISION DU1(*),DV1(*),  DU2(*),DV2(*)
      INTEGER KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*)
      INTEGER KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)

      INTEGER NEQ1, NEL1, NVT1
      INTEGER NEQ2, NEL2, NVT2

C local variables

      INTEGER IEL1, IM1, IM2, IM3, IM4, IELH1, IELH2, IELH3, IELH4
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
      DOUBLE PRECISION DUH1, DUH2, DUH3, DUH4, DUH5, DUH6, DUH7, DUH8
      DOUBLE PRECISION DUH9, DUH10, DUH11, DUH12
      DOUBLE PRECISION DVH1, DVH2, DVH3, DVH4, DVH5, DVH6, DVH7, DVH8
      DOUBLE PRECISION DVH9, DVH10, DVH11, DVH12

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

************************************************************************
* Routines for prolongation/restriction of E03x/EM3x can be found in
* MGROUT30.F/MGROUT31.F !
* Routines for prolongation/restriction of E010 (pressure) can be
* found in MGROUT10.F !
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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
C
C parameters

      DOUBLE PRECISION DX(*),DAX(*),A1,A2
      INTEGER NEQ

C externals

      EXTERNAL MTMUL
      
C local variables

      INTEGER ISETLV,I2,IP,KMBD,NMBD
      
C=======================================================================
C     Getting all parameters for MTMUL
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
      CALL MTMUL(DAX(1),DAX(I2),DAX(IP),DX(1),DX(I2),DX(IP),A1,A2,
     *            DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'

C parameters

      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ)

C externals

      EXTERNAL BDRY0

C local variables

      INTEGER KMBD, NMBD, ISETLV, I2

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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'

C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ), RHO
      
C local variables

      INTEGER ISETLV, I2, IP
      INTEGER KVERT,KMID,KNPR,KMBD,NMBD,KAREA,ITE
      DOUBLE PRECISION DMAXU,DMAXP,RES,RESINI

C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================

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
       CALL  VANCAE (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *               DB(1),DB(I2),DB(IP),
     *               DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *               NU,NP,KWORK(KMBD),KWORK(KVERT),
     *               KWORK(KMID),KWORK(KNPR),NMBD,
     *               NVT,NEL,RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),DWORK(KAREA),NP,INEUM)
C
       RES=MAX(DMAXU,DMAXP)
       IF (ITE.EQ.1) RESINI=RES
       RHO=(RES/RESINI)**(1D0/DBLE(ITE))
ccc       
C       write(6,*) ite,res
C       ,DMAXU,DMAXP,DMPSL*RESINI,RHO
       IF ((RES.LT.EPSSL).AND.(RES.LT.DMPSL*RESINI)) GOTO 99999
       
       IF (RES.GT.1E6*RESINI) THEN
         WRITE (MTERM,*) 'Warning: coarse grid solver not convergent!'
         GOTO 99999
       END IF
C
11     CONTINUE
C
      ENDIF
C
C
C
99999  END
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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'

C parameters
      
      INTEGER NEQ,ITE
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ), RHO,EPS
      
C local variables

      INTEGER ISETLV, I2, IP
      INTEGER KVERT,KMID,KNPR,KMBD,NMBD,KAREA
      DOUBLE PRECISION DMAXU,DMAXP,RES,RESINI

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
     *               DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *               NU,NP,KWORK(KMBD),KWORK(KVERT),
     *               KWORK(KMID),KWORK(KNPR),NMBD,
     *               NVT,NEL,RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),DWORK(KAREA),NP,INEUM)
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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'
      
C parameters

      DOUBLE PRECISION DUF(*),DUC(*)

C local variables
      
      INTEGER ISETLV, I2F, IPF, KVERTF, KMIDF, KADJF, I1, KPLF
      INTEGER NUC, NPC, NVTC, NMTC, I2C, IPC, KVERTC, KMIDC, KADJC, KPLC
      INTEGER KCRVGC,KAREAC
      
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
      KCRVGC=L(KLCVG (I1))
      KAREAC=L(KLAREA(I1))
C
C=======================================================================
C
      CALL PROLU (DUC(1),DUC(I2C),DUC(IPC),DUF(1),DUF(I2F),DUF(IPF),
     *            KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),NUC,NPC,NVTC,
     *            KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),NU,NP,NVT,
     *            DWORK(KCRVGC),DWORK(KAREAC))
C
      IF ((ABS(IINT).EQ.2).OR.(ABS(IINT).EQ.4)) THEN

C Use defect vector space as auxiliary vector. LD1 has size
C 2*NVT(NLMAX)+NEL(NLMAX) and is therefore large enough to carry the
C linear pressure on the coarse- and fine grid of level <= NLMAX.

       KPLC=L(LD1)
       KPLF=L(LD1)+NMTC

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

      END

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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'
      
C parameters

      DOUBLE PRECISION DDF(*),DDC(*)

C local variables
      
      INTEGER ISETLV, I2F, KVERTF, KMIDF, KADJF, I1, KPLF
      INTEGER NUF, NPF, NVTF, NMTF, IPF
      INTEGER I2C, IPC, KVERTC, KMIDC, KADJC, KPLC
      INTEGER KCRVGC,KAREAC

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
      KCRVGC=L(KLCVG (I1))
      KAREAC=L(KLAREA(I1))

C=======================================================================

      CALL RESTRD(DDC(1),DDC(I2C),DDC(IPC),DDF(1),DDF(I2F),DDF(IPF),
     *            KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),NU,NP,NVT,
     *            KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),NUF,NPF,NVTF,
     *            DWORK(KCRVGC),DWORK(KAREAC))

      IF ((ABS(IINT).EQ.2).OR.(ABS(IINT).EQ.4)) THEN

C Use defect vector space as auxiliary vector. LD1 has size
C 2*NVT(NLMAX)+NEL(NLMAX) and is therefore large enough to carry the
C linear pressure on the coarse- and fine grid of level <= NLMAX.

       KPLF=L(LD1)
       KPLC=L(LD1)+NMTF

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
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'
      
C parameters

      INTEGER NEQ, NSMO
      DOUBLE PRECISION DX(NEQ), DB(NEQ), DD(NEQ)

C local variables
      
      INTEGER ISETLV, I2, IP, KVERT, KMID, KNPR, KMBD, NMBD, KAREA
      INTEGER ITE

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
     *              DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *              NU,NP,KWORK(KMBD),KWORK(KVERT),KWORK(KMID),
     *              KWORK(KNPR),NMBD,NVT,NEL,RLXSM)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),DWORK(KAREA),NP,INEUM)
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
      SUBROUTINE   YSTEP (DX,DD,DB,NEQ,ALPHA)  
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
      IMPLICIT NONE
      
C standard COMMON blocks
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'

C parameters

      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ)
      DOUBLE PRECISION ALPHA

C local variables

      INTEGER ISETLV, IEQ
      DOUBLE PRECISION DBX, DBY

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
