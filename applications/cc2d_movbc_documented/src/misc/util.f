************************************************************************
*
*    Purpose:  - Interpolates the solution vector (DU1,DU2) to
*                the vector (DL1,DL2) of dimension NVT with
*                values in the vertices
*              - if BNODBD=FALSE,
*                the values of vertices at the boundary are computed
*                via the exact function UE; otherwise these
*                values are computed by interpolation like inner points
*                (can be used for debugging to prevent writing
*                out of hardcoded boundary-conditions)
*
************************************************************************

      SUBROUTINE INTUVD (DU1,DU2,DL1,DL2,DAUX,NVT,NEL,NVBD,KNPR,
     *                   KMID,KVERT,KVBD,KMBD,DCORVG,UE,BNODBD)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters      

      DOUBLE PRECISION DU1(*),DU2(*),DL1(*),DL2(*),DAUX(*),DCORVG(2,*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KMBD(*),KNPR(*)
      INTEGER NVT,NEL,NVBD
      LOGICAL BNODBD
      
      DOUBLE PRECISION UE
      EXTERNAL UE

C local variables

      INTEGER IVT,IEL
      INTEGER IM1,IM2,IM3,IM4,IV1,IV2,IV3,IV4,IV
      INTEGER I,I1,I0,IM0
      INTEGER INPR1,INPR2,INPR3,INPR4
      DOUBLE PRECISION DUH1,DUH2,DUH3,DUH4,DVH1,DVH2,DVH3,DVH4
      DOUBLE PRECISION X,Y

C-----------------------------------------------------------------------
      DO 1 IVT=1,NVT
      DL1 (IVT)=0D0
      DL2 (IVT)=0D0
1     DAUX(IVT)=0D0
C
C
      DO 10 IEL=1,NEL
C
      IM1=KMID(1,IEL)-NVT
      IM2=KMID(2,IEL)-NVT
      IM3=KMID(3,IEL)-NVT
      IM4=KMID(4,IEL)-NVT
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      DVH1=DU2(IM1)
      DVH2=DU2(IM2)
      DVH3=DU2(IM3)
      DVH4=DU2(IM4)
C
      DAUX(IV1)=DAUX(IV1)+1D0
      DAUX(IV2)=DAUX(IV2)+1D0
      DAUX(IV3)=DAUX(IV3)+1D0
      DAUX(IV4)=DAUX(IV4)+1D0
C
      DL1(IV1)=DL1(IV1) + 0.75D0*(DUH1+DUH4) - 0.25D0*(DUH2+DUH3)
      DL1(IV2)=DL1(IV2) + 0.75D0*(DUH2+DUH1) - 0.25D0*(DUH3+DUH4)
      DL1(IV3)=DL1(IV3) + 0.75D0*(DUH3+DUH2) - 0.25D0*(DUH4+DUH1)
      DL1(IV4)=DL1(IV4) + 0.75D0*(DUH4+DUH3) - 0.25D0*(DUH1+DUH2)
C
      DL2(IV1)=DL2(IV1) + 0.75D0*(DVH1+DVH4) - 0.25D0*(DVH2+DVH3)
      DL2(IV2)=DL2(IV2) + 0.75D0*(DVH2+DVH1) - 0.25D0*(DVH3+DVH4)
      DL2(IV3)=DL2(IV3) + 0.75D0*(DVH3+DVH2) - 0.25D0*(DVH4+DVH1)
      DL2(IV4)=DL2(IV4) + 0.75D0*(DVH4+DVH3) - 0.25D0*(DVH1+DVH2)
C
10    CONTINUE
C
      DO 20  IV=1,NVT
      DL1(IV)=DL1(IV)/DAUX(IV)
      DL2(IV)=DL2(IV)/DAUX(IV)
 20   CONTINUE
 
C Skip boundary value recomputation if this should not be done:

      IF (BNODBD) GOTO 99999
C
C=======================================================================
C *** boundary values
C=======================================================================
C
      DO 30  I=1,NVBD
      I1=I
      IF (I.EQ.1) THEN 
       I0=NVBD
      ELSE
       I0=I-1
      ENDIF
C
      IV =KVBD(I)
      IM1=KMBD(I1)
      IM0=KMBD(I0)
      IF ((IM1.LT.0).OR.(IM0.LT.0)) GOTO 30
C
      X=DCORVG(1,IV)
      Y=DCORVG(2,IV)
      DL1(IV)=UE(X,Y,1)
      DL2(IV)=UE(X,Y,2)
c      DL1(IV)=0D0
c      DL2(IV)=0D0
30    CONTINUE
C
99999 END

************************************************************************
* Interpolates the solution pressure DP to the vector DPL of dimension 
* NVT with values in the vertices.

* In:
*  DP     - Pressure to be interpolated. Array [1..NVT] of double
*  DAUX   - Auxiliary vector. Array [1..NVT] of double
*  AREA   - Array with area of elements. Array [1..NVT] of double
*
* Out:
*  DPL    - Interpolated pressure. Array [1..NVT]
************************************************************************

      SUBROUTINE   INTPV (DP,DPL,DAUX,AREA,KVERT,KNPR)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
C parameters
      
      DOUBLE PRECISION DP(*),DPL(*),DAUX(*)
      INTEGER KVERT(NNVE,*),KNPR(*)
      DOUBLE PRECISION AREA(*)

C local variables

      INTEGER IVT,IEL
      INTEGER IV1,IV2,IV3,IV4
      INTEGER INPR1,INPR2,INPR3,INPR4

      DOUBLE PRECISION DPIEL, DAREA

C-----------------------------------------------------------------------
c
      DO 1 IVT=1,NVT
      DPL (IVT)=0D0
1     DAUX(IVT)=0D0
C
C
      DO 10 IEL=1,NEL
C
      DPIEL=DP(IEL)
      DAREA=AREA(IEL)
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
C
      DPL(IV1)=DPL(IV1)+0.25D0*DAREA*DPIEL
      DPL(IV2)=DPL(IV2)+0.25D0*DAREA*DPIEL
      DPL(IV3)=DPL(IV3)+0.25D0*DAREA*DPIEL
      DPL(IV4)=DPL(IV4)+0.25D0*DAREA*DPIEL
C
      DAUX(IV1)=DAUX(IV1)+0.25D0*DAREA
      DAUX(IV2)=DAUX(IV2)+0.25D0*DAREA
      DAUX(IV3)=DAUX(IV3)+0.25D0*DAREA
      DAUX(IV4)=DAUX(IV4)+0.25D0*DAREA
C
10    CONTINUE
C
C
      DO 20 IVT=1,NVT
20    DPL(IVT)=DPL(IVT)/DAUX(IVT)
C
C
C      VPH=0E0
C      DO 30 IVT=1,NVT
C30    VPH=VPH+VPL(IVT)
C      VMWP=VPH/REAL(NVT)
C
C
C      DO 40 IVT=1,NVT
C40    VPL(IVT)=VPL(IVT)-VMWP
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   SETARE  (AREA,NEL,KVERT,DCORVG)
************************************************************************
*
*   Purpose: - writes on  AREA(IEL)  the area of the element IEL,
*              IEL=1,...,NEL
*            - writes on  AREA(NEL+1) the sum of all  AREA(IEL)
*            - KVERT,DCORVG are the usual FEAT arrays
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C parameters

      INTEGER NEL
      DOUBLE PRECISION  AREA(*)
      INTEGER KVERT(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      
C local variables

      DOUBLE PRECISION SUM,AAA
      INTEGER IEL,I1,I2,I3,I4
      DOUBLE PRECISION X1,X2,X3,X4,Y1,Y2,Y3,Y4
      
C=======================================================================
      SUM=0.D0
      DO  11  IEL=1,NEL
C
      I1=KVERT(1,IEL)
      I2=KVERT(2,IEL)
      I3=KVERT(3,IEL)
      I4=KVERT(4,IEL)
C
      X1=DCORVG(1,I1)
      X2=DCORVG(1,I2)
      X3=DCORVG(1,I3)
      X4=DCORVG(1,I4)
C
      Y1=DCORVG(2,I1)
      Y2=DCORVG(2,I2)
      Y3=DCORVG(2,I3)
      Y4=DCORVG(2,I4)
C
      AAA=0.5D0*(  DABS((X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2))
     *            +DABS((X1-X4)*(Y3-Y4)-(Y1-Y4)*(X3-X4)) )
      AREA(IEL)=AAA
      SUM=SUM+AAA
  11  CONTINUE
C
      AREA(NEL+1)=SUM
C
      END

************************************************************************
* Transform vector to L2_0
*
* This routine transforms the vector P into the space L2_0.
* For this purpose, the vector AREA with the areas of all elements
* is used. AREA(NEL+1) is assumed to be the sum of all AREA(IEL).
*
* In:
*   NEL    - Number of elements in the geometry
*   P      - array [1..NEL] of double
*            Pressure vector in the space P0, one pressure value
*            per element
*   AREA   - array [1..NEL+1] of double
*            Area of all elements; AREA[NEL+1] is the sum of the
*            areas of all elements
*   INEUM  - =0: No Neumann boundary edges in the geometry
*            =1: There are Neumann boundary edges in the geometry.
*
* Out:
*   P      - the modified vector
*
* If INEUM=1, the vector is assumed to be already in L2_0, so nothing
* is done.
************************************************************************

      SUBROUTINE TOL20A (P,AREA,NEL,INEUM)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C parameters

      INTEGER INEUM, NEL
      DOUBLE PRECISION AREA(*)
      DOUBLE PRECISION P(*)
      
C local variables

      DOUBLE PRECISION PINT,C
      INTEGER IEL

      IF (INEUM.EQ.1) RETURN

C     Build the integral
C       int_Omega p dx
C     This is approximated by
C       PINT = SUM_Elements P(Element)*Volume(Element)

      PINT=0D0
      DO IEL=1,NEL
        PINT=PINT+P(IEL)*AREA(IEL)
      END DO

C     Divide PINT by the volume of the domain; this gives the integral
C     mean value of the pressure:

      C=PINT/AREA(NEL+1)

C     Subtract the integral mean value C of the pressure from all
C     pressure components. Afterwards, P has integral mean value = 0.

      DO IEL=1,NEL
        P(IEL)=P(IEL)-C
      END DO

      END

************************************************************************
* Set all data for current level ILEV (from /MGPAR/)
*
* This routine changes the values in the /TRIAx/ COMMON blocks
* according to the level ILEV in the multigrid COMMON block.
* The parameter ISETLV decides on which information is initialised.
*
* In:
*   ISETLV - >=1  - update of /TRIAA/,TRIAD/   
*            >=2  - update of /LEVDIM/,/ADRFLD/
*
* In (from COMMON block /MGPAR/):
*   ILEV   - current level number (from /MGPAR/)
*
* Out:
*   The COMMON blocks /TRIAx/ are set according to ISETLV and ILEV.
************************************************************************

      SUBROUTINE SETLEV (ISETLV)

      IMPLICIT NONE

C include the necessary COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'

C parameters

      INTEGER ISETLV

C-----------------------------------------------------------------------
C *** elementary check
      IF ( (ILEV.GT.NLMAX) .OR. (ISETLV.LT.1) .OR. (ISETLV.GT.2) .OR.
     *     ((ISETLV.EQ.2).AND.(ILEV.LT.NLMIN)) )  THEN
        WRITE(MTERM,*) 'ERROR in SETLEV: ILEV or ISETLV is wrong'
        STOP
      ENDIF
C
C *** update of /TRIAD/,/TRIAA/
C
      NEL =KNEL  (ILEV)
      NVT =KNVT  (ILEV)
      NMT =KNMT  (ILEV)
      NVEL=KNVEL (ILEV)
      NVBD=KNVBD (ILEV)
C
      LCORVG=KLCVG  (ILEV)
      LCORMG=KLCMG  (ILEV)
      LVERT =KLVERT (ILEV)
      LMID  =KLMID  (ILEV)
      LADJ  =KLADJ  (ILEV)
      LVEL  =KLVEL  (ILEV)
      LMEL  =KLMEL  (ILEV)
      LNPR  =KLNPR  (ILEV)
      LMM   =KLMM   (ILEV)
      LVBD  =KLVBD  (ILEV)
      LEBD  =KLEBD  (ILEV)
      LBCT  =KLBCT  (ILEV)
      LVBDP =KLVBDP (ILEV)
      LMBDP =KLMBDP (ILEV)
C
C *** update of /LEVDIM/,/ADRFLD/ if  ISETLV=2
C
      IF (ISETLV.EQ.2)  THEN
C
         NA =KNA  (ILEV)
         NB =KNB  (ILEV)
         NU =KNU  (ILEV)
         NP =KNP  (ILEV)
         NUP=KNUP (ILEV)
C
         KA1  =L(KLA    (ILEV))
         KST1 =L(KLST   (ILEV))
         KM1  =L(KLM    (ILEV))
         KCOLA=L(KLCOLA (ILEV))
         KLDA =L(KLLDA  (ILEV))
         KB1  =L(KLB1   (ILEV))
         KB2  =L(KLB2   (ILEV))
         KCOLB=L(KLCOLB (ILEV))
         KLDB =L(KLLDB  (ILEV))
         KU1  =L(KLUP   (ILEV))
         KU2  =KU1+NU
         KP   =KU2+NU
         KF1  =L(KLF12P (ILEV))
         KF2  =KF1+NU
         KFP  =KF2+NU
         KAUX1=L(KLAUX  (ILEV))
         KAUX2=KAUX1+NU
         KAUXP=KAUX2+NU
C
      ENDIF
C
      END

************************************************************************
* Calculate streamfunction
*
* Based on a (U,V) solution vector, this routine routine calculates
* the (scalar) streamfunction in all vertices of the triangulation.
* U,V are expected as point values in the corners of the elements
* in a given triangulation.
*
* In:
*   DU1,
*   DU2    : array [1..NMT] of double
*            X- and Y-velocity field, midpoint/edge-based
*   DCORVG,
*   KVERT,
*   KMID,
*   KADJ,
*   NVT,
*   NVE,
*   NEL    : Usual geometry information; must correspond to TRIA!
*   DVIND  : array [1..NVT] of double
*            Auxiliary array.
*
* Out:
*   DX     : array [1..NVT] of double precision
*            The stream function in all (corner-) vertices.
************************************************************************

      SUBROUTINE U2ISO (DCORVG,KVERT,KMID,KADJ,NVT,NEL,NVE,
     *                  DVIND,DX,DU1,DU2)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C     parameters

      DOUBLE PRECISION DCORVG(2,*),DVIND(*),DX(*),DU1(*),DU2(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*)
      INTEGER NVT,NEL,NVE

C local variables

      INTEGER IVT, IEH, IEL, IH, IVE, JVE, IHV, IND, INDH, IVTH, IMID
      INTEGER IME, IELH
      DOUBLE PRECISION PX1,PY1,PX2,PY2,DN1,DN2,DXH

C     Clear the solution and auxiliary array

      CALL LCL1(DX,NVT)
      CALL LCL1(DVIND,NVT)

C     Start with element 1 and its first vertex.
C     Set the streamfunction there to 0.0. The streamfunction
C     in the other vertices are then relatively calculated 
C     to this basis.
C     DVIND is a marker which is set to 1 for every node that 
C     is finished.

      DX(KVERT(1,1))    = 0D0
      DVIND(KVERT(1,1)) = 1D0

C     Loop over the elements:

      DO IEH=1,NEL

C       We set the current element IEL to IEH:

        IEL=IEH
        IH=0

C       On the current element, loop over the vertices of
C       that element. Add the DVIND-values of all vertices of the
C       current element together. So, IH gets the number of marked
C       vertices and IND will be the number of the "last" marked
C       vertex of the four on the element.

        DO IVE = 1,NVE
          JVE=KVERT(IVE,IEL)
          IHV=INT(DVIND(JVE))
          IH=IH+IHV
          IF (IHV.GE.1) IND=IVE
        END DO

C       If all four vertices are marked, there's nothing to calculate
C       on the element. If no vertex is marked, we can't calculate
C       anything on the current element. In both cases skip the
C       computation and search for a better element:

        IF ((IH.GE.NVE).OR.(IH.EQ.0)) GOTO 20

C       Ok, here we found an element where some of the vertices are
C       marked and some not. Here we can calculate a part of the
C       streamfunction.

13      CONTINUE  
      
C       Loop over the vertices on the element:
      
        DO IVE=1,NVE-1

C         IND is the index of a marked vertex. Calculate the "next"
C         vertex following IND and its vertex number into IVH.

          INDH=MOD(IND,NVE)+1
          IVTH=KVERT(INDH,IEL)

C         If that vertex is not marked, the streamfunction is not
C         calculated there. Otherwise we are just looking at two 
C         marked neighboured vertices, so there's nothing to gain here. 

          IF (DVIND(IVTH).LT.1D0) THEN
          
C           Vertex IVT (corresponding to IND) is marked, vertex IVTH 
C           is not. 
C
C               x---x IVTH
C               |   |
C               x---O IVT
C
C           Mark vertex IVTH to indicate that the streamfunction
C           is not being calculated there:
          
            DVIND(IVTH)=1D0

C           and calculate the streamfunction in INTH.

            IVT =KVERT(IND,IEL)
            IMID=KMID (IND,IEL)-NVT
            
C           IMID is now the midpoint number following IVT and thus
C           the number of the DOF in the FE function.
C           Calculate the normal vector of the current edge into
C           N=(DN1,DN2) - not normalized.
C
C               x-------x IVTH
C               |       |
C               |  IMID x--> (DN1,DN2)
C               |       |
C               x-------O IVT

            PX1=DCORVG(1,IVT)
            PY1=DCORVG(2,IVT)
            PX2=DCORVG(1,IVTH)
            PY2=DCORVG(2,IVTH)
            DN1 = PY2-PY1
            DN2 =-PX2+PX1
            
C           Calculate the streamfunction in IVTH from the value
C           in IVT by:
C
C           sfc(IVTH) = sfc(IVT) + U(IMID)*N
C
C           which is the "amount of flow over the edge (IVT,IVTH)".

            DX(IVTH)=DX(IVT)+(DU1(IMID)*DN1+DU2(IMID)*DN2)
          
          END IF ! (DVIND(IVTH) < 1D0)

C         Go on to the next vertex on the element to look if that one
C         has a not marked neighbour.

          IND=INDH
            
        END DO ! IVE
        
C       Now on the current element IEL, on all (corner) vertices the
C       streamfunction is calculated. We go on looking to the adjacent
C       elements of IEL if there's an element where the streamfunction
C       is not calculated in all vertices...

20      CONTINUE

C       Look onto the adjacent elements of the current element if there's
C       a suitable neighbour element where we can continue the calculation.
C
C       Loop over the edges of the current element

        DO IME=1,NVE

C         Get the neighbour element adjacent to the current element

          IELH=KADJ(IME,IEL)
          
          IF (IELH.NE.0) THEN
          
C           Now we have the number of the neighbour element in IELH.
C           Loop about the vertices of the element and sum up the
C           markers into IH.
          
            IH=0
            DO IVE=1,NVE
              JVE=KVERT(IVE,IELH)
              IHV=INT(DVIND(JVE))
              IH=IH+IHV
              IF (IHV.GE.1) INDH=IVE
            END DO
            
C           If there is at least one but not all markers set, the
C           element can be used for further calculation.
C           Switch the current element IEL to that one and
C           continue the calculation there.

            IF ((IH.LT.NVE).AND.(IH.GT.0)) THEN
              IEL=IELH
              IND=INDH
              GOTO 13
            END IF
            
          END IF ! IELH <> 0

        END DO ! IME
      
      END DO ! IEH

C     At last, normalize the streamfunction such that vertex 1
C     has value 0.0. Remember that we assigned a value of 0.0
C     to the first vertex of element 1, which is usually not
C     vertex 1 of the triangulation!

      DXH=DX(1)
      DO IVT=1,NVT
        DX(IVT)=DX(IVT)-DXH
      END DO

      END

************************************************************************
* Constant <-> nonconforming interpolation
* 
* This converts between constant and nonconforming piecewise linear
* FE pressure functions.
* IPAR=0 : convert piecewise constant DPC to piecewise linear DPL
* IPAR=1 : convert piecewise linear DPL to piecewise constant DPC
*
* DPC    : array [1..NEL] of double
* DPL    : array [1..NMT] of double
************************************************************************

      SUBROUTINE C2N2DM (DPC,DPL,KMID,KADJ,NEL,NMT,NVT,IPAR)
************************************************************************
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DPC(*),DPL(*)
      INTEGER KMID(NNVE,*),KADJ(NNVE,*)
      INTEGER NEL,NMT,NVT,IPAR

C local variables

      DOUBLE PRECISION DPH
      INTEGER IEL, IVE, IADJ, IMID

C-----------------------------------------------------------------------
C
C *** constant to nonconforming = 0
      IF (IPAR.EQ.0) THEN
C
       DO 10 IEL=1,NEL
       DPH  =DPC(IEL)
C
       DO 20 IVE=1,4
       IADJ=KADJ(IVE,IEL)
       IMID=KMID(IVE,IEL)-NVT
C
       IF (IADJ.EQ.0)   DPL(IMID)=DPH
       IF (IADJ.GT.IEL) DPL(IMID)=0.5D0*(DPH+DPC(IADJ))
C
20     CONTINUE
10     CONTINUE
C
      ELSE
C
       DO 110 IEL=1,NEL
       DPC(IEL)=0.25D0*( DPL(KMID(1,IEL)-NVT)+DPL(KMID(2,IEL)-NVT)
     *                  +DPL(KMID(3,IEL)-NVT)+DPL(KMID(4,IEL)-NVT))
110    CONTINUE
C
      ENDIF      

      END

************************************************************************
* Grid disturbance
*
* This routine stochastically disturbes the current grid.
* Every gridpoint is moved DIEPS percent in a - more or less random -
* direction. The grid is assumed to be uniform, the minimum cell
* size is computed with the help of NVT!
*
* In:
*  DCORVG, KNPR, NVT - as usual
*  DIEPS             - Rate of disturbing. 0.2D0 = 20%
*
* The current grid will be modified directly.
************************************************************************
      SUBROUTINE GRDIST(DCORVG,KNPR,NVT,DIEPS)
      
      IMPLICIT NONE
      
C parameters
      
      DOUBLE PRECISION DCORVG(2,*),DIEPS
      INTEGER KNPR(*),NVT
      
C local variables

      DOUBLE PRECISION H, HDIST
      INTEGER IVT
      
      H=1D0/(SQRT(DBLE(NVT))-1)
      HDIST=DIEPS*H
C
      DO 1 IVT=1,NVT
      IF (KNPR(IVT).EQ.0) THEN
C       I17=17
C       J1=MOD(IVT,I17)
C       C1=DBLE(J1)
C       C2=(-1D0)**C1   ! Be careful with this statement: -Real**Real=not working on some machines
C       DCORVG(1,IVT)=DCORVG(1,IVT)+C2*HDIST
       DCORVG(1,IVT)=DCORVG(1,IVT)+DBLE((-1)**MOD(IVT,17))*HDIST
C       C1=DBLE(IVT)
C       C2=(-1D0)**C1   ! Be careful with this statement: -Real**Real=not working on some machines
C       DCORVG(2,IVT)=DCORVG(2,IVT)+C2*HDIST
       DCORVG(2,IVT)=DCORVG(2,IVT)+DBLE((-1)**IVT)*HDIST
      ENDIF
1     CONTINUE
C
      END
C
C
C
************************************************************************
C Change coordinates
C
C The follwing function corrects the grid on domains where non-linear/
C curved boundary segments are used. 
C When a boundary segment of the domain is a line, new boundary nodes
C are automatically positioned on the boundary. But when the boundary
C is a curve, a new boundary vertex is not automatically on the
C boundary, but it first has tobe moved to there. This is already
C performed in the refinement routine XSB0X. Unfortunately this
C procedure can lead to very anisotropic elements near the boundary,
C depending on how sharp the curved boundary is.
C
C CHCOOR now tries to reduce these effects of anisotropy. The element
C midpoint of the coarser element (which is at the same time the vertex where
C all the four finer elements meet) is taken as the average of the four
C edge midpoints that arise from natural refinement.
C
C o----------o----------o
C \          \          |
C |\         \          | 
C | \      -> \         |
C |  o---------o--------o   Averaging of the element midpoint by interpolation
C | /      -> /         |
C |/         /          |
C /          /          |
C o----------o----------o

      SUBROUTINE CHCOOR(DCORVG,KVERT,KADJ,NEL)
      
      IMPLICIT NONE
      
C parameters
      
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KVERT(4,*),KADJ(4,*)
      INTEGER NEL

C local variables

      INTEGER IADJ3, IADJ4
      INTEGER IVT1, IVT2, IVT3, IVT4
      INTEGER IVTM, IEL
      DOUBLE PRECISION PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,PXM,PYM,PX,PY
C
      DO 10 IEL=1,NEL/4
C
      IADJ3=KADJ(2,IEL)
      IADJ4=KADJ(3,IEL)
      IVT1=KVERT(2,IEL)
      IVT2=KVERT(4,IEL)
      IVT3=KVERT(2,IADJ3)
      IVT4=KVERT(4,IADJ4)
      IVTM=KVERT(3,IEL)
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PX3=DCORVG(1,IVT3)
      PX4=DCORVG(1,IVT4)
C
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PY3=DCORVG(2,IVT3)
      PY4=DCORVG(2,IVT4)
C
      PXM=DCORVG(1,IVTM)
      PYM=DCORVG(2,IVTM)
C
      PX=0.25D0*(PX1+PX2+PX3+PX4)
      PY=0.25D0*(PY1+PY2+PY3+PY4)
C
      DCORVG(1,IVTM)=PX
      DCORVG(2,IVTM)=PY
C
10    CONTINUE
C
      END

************************************************************************
* Calculate stopping criterion for adaptive time stepping
*
* Based on the configuration, CRITAD will compute a bound EPSAD
* that can be used in the adaptive time stepping as stopping criterion.
*
* In:
*   TIMEIN  - Length of startup phase
*   TIMENS  - Current simulation time
*   TIMEST  - Initial simulation time
*   EPSADI  - error limit during startup phase
*   EPSADL  - standard error limit after startup phase
*   IADIN   - Identifier for the error control in the startup phase
*             =0: use EPSADI during startup phase
*             =1: use linear blending from EPSADI to EPSADL during
*                 startup phase
*             =2: use logarithmic blending from EPSADI to EPSADL during
*                 startup phase
*
* Out:
*   EPSAD   - Stopping criterion for adaptive time step control.
*             During the startup phase (time T in the range 
*             TIMEST..TIMEST+TIMEIN), the stopping criterion
*             will be calculated according to IADIN. After the
*             startup phase, EPSADL will be used.
************************************************************************

      SUBROUTINE CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C parameters

      DOUBLE PRECISION TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD
      INTEGER IADIN

C local variables

      DOUBLE PRECISION TDIFF
             
C     Standard result: EPSADL
             
      EPSAD=EPSADL

C     Do we have a startup phase at all? Would be better, otherwise
C     Navier-Stokes might stagnate with decreasing time steps at the
C     beginning...

      IF (TIMEIN.GT.0) THEN
      
C       Calculate the "current position" of the simulation time in
C       the interval TIMEST..TIMEST+TIMEIN
      
        TDIFF=TIMENS-TIMEST

C       Standard error control in the startup phase

        IF (IADIN.EQ.0) THEN
        
C         Use EPSADI during startup phase and EPSADL afterwards:
        
          IF (TDIFF.LE.TIMEIN) THEN
            EPSAD=EPSADI
          ELSE
            EPSAD=EPSADL
          ENDIF
          
        ENDIF

C       Linear blending as control in the startup phase.
C       Blend linearly between EPSADI and EPSADL from the initial
C       simulation time TIMNEST to the end of the startup phase
C       TIMEST+TIMEIN.
C       After the startup phase, use EPSADL.

        IF (IADIN.EQ.1) THEN
          IF (TDIFF.LE.TIMEIN) THEN
            EPSAD=EPSADI+TDIFF/TIMEIN*(EPSADL-EPSADI)
          ELSE
            EPSAD=EPSADL
          ENDIF
        ENDIF

C       Logarithmic blending as control in the startup phase.
C       Blend logarithmically between EPSADI and EPSADL from the initial
C       simulation time TIMNEST to the end of the startup phase
C       TIMEST+TIMEIN.
C       After the startup phase, use EPSADL.

        IF (IADIN.EQ.2) THEN
          IF (TDIFF.LE.TIMEIN) THEN
            EPSAD=EPSADI**(1D0-TDIFF/TIMEIN)*EPSADL**(TDIFF/TIMEIN)
          ELSE
            EPSAD=EPSADL
          ENDIF
        ENDIF

      ENDIF

      END
