************************************************************************
* Update/mark DIRICHLET- and NEUMANN components in the geometry vectors
*
* This routine modifies the geometry vectors KNPR,KMBD to mark which
* points are dirichlet-nodes and which points are Neumann nodes.
*
* Result:
*  KNPR  - modified nodal-property vector
*  KMBD  - 
*  NMBD  - Number of moving boundary nodes
*  INEUM - >0, if there are Neumann boundaries in the geometry
*
* The given geometry vectors are updated for Dirichlet- and Neumann-
* components as defined in the input parameters
* (-> indat2d.f, fctbdry.f).
************************************************************************

************************************************************************
* Update boundary information
*
* This routine updates/marks DIRICHLET and NEUMANN boundary components
* in a triangulation structure. There is a loop through all vertices
* and edges that tests which are DIRICHLET elements, NEUMANN elements
* and fictitious boundary elements.
*
* In:
*   KVBD,
*   KEBD,
*   KVERT,
*   KMID,
*   KNPR
*   DMBDP,
*   DCORVG,
*   KADJ,
*   NVBD,
*   NVT,
*   NEL,   - Usual geometry information
*
* Out:
*   INEUM  - >0, if there are Neumann boundaries in the geometry
* 
*   KNPR   - Is build according to the geometry information:
*            The KNPR values of all Neumann-type edges on the
*            real boundary are overwritten by 0. The KNPR values
*            of all fictitious boundary edges in the geometry
*            are set to a value <> 0.
*   NMBD   - Number of edges on the real boundary
*   KMBD   - array [1..NMT+NVBD] of integer
*            "Shortcut nodal property" array. This array stores all
*            edges on the boundary. The number is >= for all
*            Dirichlet edges and < 0 for all edges that should be
*            treated as Neumann edges.
*            The array has a maximum size on NMT, as it stores not
*            only edges on the real boundary, but also probably
*            fictitious boundary edges inside of the geometry - so
*            in the worst case all edges!
*   DDBD   - array [1..NVBD] of double
*            The parameter value of each midpoint on the boundary.
************************************************************************

      SUBROUTINE BDRNEU (KMBD,KVBD,KEBD,KVERT,KMID,KNPR,DDBD,DMBDP,
     *                   DCORVG,KADJ,NVBD,NVT,NEL,
     *                   NMBD,INEUM)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cbasictria.inc'

      INCLUDE 'cnsparfrac.inc'

C constants

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
C parameters

      INTEGER KMBD(*),KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KNPR(*),KADJ(4,*)
      DOUBLE PRECISION DDBD(*),DMBDP(*),DCORVG(2,*)
      INTEGER NVBD, NVT, NEL
      INTEGER NMBD, INEUM
      
C externals

      INTEGER ISFBDY,ISFBDM,NFBDYC
      EXTERNAL ISFBDY,ISFBDM,NFBDYC

C local variables

      INTEGER IVBD, IMID, INPR, INPART, INPRN, NPART
      INTEGER IV1, IV2, IEL, IVE, IVT1, IVT2, IMBCP, IVE1, IM1
      DOUBLE PRECISION DPAR, DPARN1, DPARN2, XX, YY, XX1, XX2, YY1, YY2

      NMBD =0
      INEUM=0

      DO IVBD=1,NVBD
      
C       Get all information about the edge on the boundary:
      
        CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
      
        IF (IMID.NE.0) THEN

C         Get the parameter value of the edge midpoint

          DPAR=DMBDP(IVBD)
          
C         Add the edge to KMBD - with a positive value to indicate
C         a Dirichlet edge at first.
          
          NMBD=NMBD+1
          KMBD(NMBD)=IMID
          DDBD(NMBD)=DPAR

C         At first ask the user-defined routines how many
C         Neumann boundary components exist:

          INPART=0
          CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
          NPART=INPART

C         NPART is now the number of Neumann components.
C
C         Loop through the NPARTS boundary parts and check, if the
C         vertex/edge belongs to a Neumann segment:
          
          DO INPART=1,NPART

C           Get the parameter values DPARN1,DPARN2 of the boundary segment
C           from the user-defined routines

            CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)

C           Check if DPARV is inside of that segment. 

            IF ((DPAR.GT.DPARN1).AND.(DPAR.LT.DPARN2)
     *                          .AND.(INPR.EQ.INPRN)) THEN

C             Yes, we have at least one Neumann boundary segment.
C             Set KNPR appropriately and add the edge to KMBD -
C             with a "negative" value to indicate a Neumann edge.

              INEUM=1
              KNPR(NVT+IMID)=0
              KMBD(NMBD)=-KMBD(NMBD)
              
            ENDIF
          END DO
        END IF

      END DO

C
C *MK*: Handling of moving boundaries
C
C If there are no fict. boundary components, we skip this...

      IF (NFBDYC().GT.0) THEN

C Loop through all elements and all edges to find the nodes 
C belonging to elements where a moving boundary goes through...
     
        DO IEL=1,NEL
        
         DO IVE=1,4
         
          IF (KADJ(IVE,IEL).LT.IEL) THEN
      
C take the X/Y-coordinate of the current vertex
      
            IVT1=KVERT(IVE,IEL)
            XX=DCORVG(1,IVT1)
            YY=DCORVG(2,IVT1)
        
C Test, if the point is in a fictitious boundary component. If yes,
C replace KNPR by a number generated from the appropriate boundary
C component.
C By convention the boundary component numbers of all fictitious
C boundary components are NBCT+1..NBCT+NFBDYC. Theoretically we
C could store this value directly, but then we woould loose
C the information that this point is in reality an inner node
C (as all real boundary nodes have KNPR-values > 0). So to
C avoint loosing this information we store a negative value to
C the KNPR entry of the points in a fictitious boundary component.
C 
C More precisely, if this point is in a fictitious boundary component, 
C set KNPR of the appropriate midpoint from 0 to 
C -(NEL+number of fict. boundary component) to indicate this.
C The standard FEAT documentation uses INPR >= -NEL, therefore
C this value describes a vertex in the inner of a fictitious boundary
C component uniquely, and we don't loose the information that the node
C is in reality an "inner" node - we can find these nodes and set the
C KNPR-value back to 0 at any time.
C
C Remark: We first handle the (corner) vertices, then the midpoints.
C Only (corner) vertices are exported later to GMV files.
C Using the KNPR-array he GMV output routine can decide which points
C belong to fictitious boundary components to "stamp out" this region.
C This has nothing to do with any calculation! As we are using
C nonconforming finite elements that are midpoint oriented, only
C the KNPR-information of the midpoints are concerned during
C the modification the matrix!

            IMBCP = ISFBDY (XX,YY,0)

            IF (IMBCP.GT.0) THEN
C              Handlung of vertices disabled, we word with
C              midpoint-oriented elements!
C              KNPR(IVT1) = IMBCP
              KNPR(IVT1) = -(NEL+IMBCP)
            ENDIF

C Take a look to the midpoint of the current element, following the
C current vertex. 

            IVE1=IVE+1
            IF(IVE1.EQ.5) IVE1=1
            IVT2=KVERT(IVE1,IEL)
            XX1=DCORVG(1,IVT2)
            YY1=DCORVG(2,IVT2)
            XX2=0.5D0*(XX+XX1)
            YY2=0.5D0*(YY+YY1)
      
C Test also for this point if it's in a fictitious boundary

            IM1=KMID(IVE,IEL)    

            IMBCP = ISFBDY (XX2,YY2,0)

C If this point is in a fictitious boundary component, set
C KNPR of the appropriate midpoint from 0 to -(NEL+1) at least, too.
C Don't add edges that are part of the real boundary. These have
C been added above and may cause an array overflow in KMBD if too
C many edges are of fictitious boundary type!
C
C Remark: The INPR >= -NEL convention in the FEAT documentation does
C only hold for (corver) vertices as only these vertices can become
C "irregular" (i.e. hanging nodes). But to avoid any confusion
C we don't exploit this fact but treat midpoints the same way.

            IM1=KMID(IVE,IEL)
            
            IF ((IMBCP.GT.0).AND.(KNPR(IM1).LE.0)) THEN

              KNPR(IM1) = -(NEL+IMBCP)
              
C Increment the number of boundary points and store the number of 
C the "new" boundary node

              NMBD=NMBD+1
              KMBD(NMBD)=IM1-NVT
            END IF
            
          END IF
         END DO
        END DO
      END IF

      END

************************************************************************
* Implement PRESSURE DROP values into the velocity RHS-vectors DF1/DF2.
*
* In:
*   KVBD,
*   KEBD,
*   KVERT,
*   KMID,
*   KNPR,
*   DCORVG,
*   DMBDP,
*   NVBD   - Usual geometry information
*   DF1,
*   DF2    - array [1..*] of double
*            X- any Y-RHS vectors which are to be modified. These
*            are assumed to be assembled by E030/EM30/E031/E031 finite
*            element.
*   TSTEPB - Length of current time step in time stepping scheme.
*
* Out:
*   DF1,
*   DF2    - Modified RHS vectors
************************************************************************

      SUBROUTINE PDSET (KVBD,KEBD,KVERT,KMID,KNPR,DCORVG,DMBDP,DF1,DF2,
     *                  TSTEPB,NVBD)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'

      INCLUDE 'cnsparfrac.inc'

      INCLUDE 'cinidat.inc'

C parameters 

      INTEGER KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
      DOUBLE PRECISION DCORVG(2,*),DMBDP(*),DF1(*),DF2(*)
      INTEGER NVBD

C local variables

      INTEGER IVBD, IMID, IV1, IV2, INPR, INPART, INPRN, NPART
      DOUBLE PRECISION DPAR, DPARN1, DPARN2, PX1, PY1, PX2, PY2
      DOUBLE PRECISION DN1, DN2, PMEAN, FDATIN, TSTEPB
      EXTERNAL FDATIN
      INTEGER TRIA

      DO IVBD=1,NVBD

        CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
        
        IF (IMID.NE.0) THEN

          DPAR=DMBDP(IVBD)

          INPART=0
          CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
          NPART=INPART
          
          DO INPART=1,NPART
        
C           Get the parameter values DPARN1,DPARN2 of the boundary segment
C           from the user-defined routines

            CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
            IF ((DPAR.GT.DPARN1).AND.(DPAR.LT.DPARN2)
     *                          .AND.(INPR.EQ.INPRN)) THEN
              PX1  =DCORVG(1,IV1)
              PY1  =DCORVG(2,IV1)
              PX2  =DCORVG(1,IV2)
              PY2  =DCORVG(2,IV2)
              DN1  =-PY2+PY1
              DN2  = PX2-PX1

C             Call the user defined callback routine to get the
C             mean pressure value:

              PMEAN= FDATIN(7,INPR,DPAR,DPAR,TIMENS,RE,TSTEP,
     *               0,TRIA,0,0)

C             Include that into the RHS vectors - cf p. 257 (235) in 
C             Turek's book.
C             The pressure drop condition is a time independent
C             right hand side condition! Therefore it's not dependent
C             on the THETA-value of any THETA-Scheme, but only
C             dependent on the length of the current time step.

              DF1(IMID)=DF1(IMID)+PMEAN*DN1*TSTEPB
              DF2(IMID)=DF2(IMID)+PMEAN*DN2*TSTEPB
            ENDIF
          END DO
        END IF

      END DO

      END

************************************************************************
* Implement DIRICHLET values into the vector (DX1,DX2).
*
* This routine can be used for updating the velocity solution
* vectors DU1/DU2 and velocity RHS vectors DF1/DF2 with the fixed
* values on the Dirichlet boundary.
*
* It's assumed that DF1/DF2 is assembled with E030/E031/EM30/EM31.
*
* In:
*   DX1,
*   DX2    - array [1..*] of double
*            X- any Y-velocity vectors which are to be modified.
*            Both vectors are assumed to correspond to values in the
*            midpoints of edges!
*   PARX,
*   PARY   - Subroutines of the used parametrization to retrieve
*            X/Y-coordinate of a parameter value.
*   UE     - Subroutine that must return the Dirichlet function value
*            in a point X/Y.
*   KNPR,
*   NVT,
*   DCORVG,
*   KADJ,
*   NEL,
*   KMID,
*   KVERT  - Usual geometry information
*
* Return:
*   DX1, 
*   DX2    - Modified vectors. The Dirichlet values are implemented.
************************************************************************

      SUBROUTINE    BDRSET  (DX1,DX2,KNPR,NVT,
     *                       PARX,PARY,UE,DCORVG,KADJ,NEL,KMID,KVERT)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'

C parameters      
      
      DOUBLE PRECISION DX1(*),DX2(*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KNPR(*),KADJ(4,*),KMID(NNVE,*),KVERT(NNVE,*)
      INTEGER NVT, NEL
      
C externals

      DOUBLE PRECISION PARX, PARY, UE, FBINDT
      EXTERNAL PARX,PARY,UE,FBINDT 

C local variables

      INTEGER INPR, IEL, IVE, IVT1, IVT2, IVE1, IM1
      DOUBLE PRECISION U1, U2, XX, YY, XX1, YY1, XX2, YY2

C Process all vertices -> loop through all elements and on
C each element through all vertices.

      DO IEL=1,NEL
        DO IVE=1,4
          IF (KADJ(IVE,IEL).LT.IEL) THEN
            IVT1=KVERT(IVE,IEL)
            XX=DCORVG(1,IVT1)
            YY=DCORVG(2,IVT1)
            IVE1=IVE+1
            IF(IVE1.EQ.5) IVE1=1
            
            IVT2=KVERT(IVE1,IEL)
            XX1=DCORVG(1,IVT2)
            YY1=DCORVG(2,IVT2)
            XX2=0.5D0*(XX+XX1)
            YY2=0.5D0*(YY+YY1)
            IM1=KMID(IVE,IEL)
            
            INPR=KNPR(IM1)

            IF (INPR.LT.-NEL) THEN
C Fictitious boundary midpoint; has INPR = -(NEL+fict. boundary number)
              U1 = FBINDT (1,XX2,YY2)
              U2 = FBINDT (2,XX2,YY2)
              DX1(IM1-NVT)=U1
              DX2(IM1-NVT)=U2
            END IF
            IF (INPR.GE.1) THEN
C Vertex on a boundary
              U1=UE(XX2,YY2,1)
              U2=UE(XX2,YY2,2)
              DX1(IM1-NVT)=U1
              DX2(IM1-NVT)=U2
            END IF
          END IF
        END DO
      END DO

      END

************************************************************************
* Multiply the NEUMANN-components of the vector DX with A1
*
* This will search the geometry for all nodes that belong to Neumann
* boundary of real boundary components. Fictitious boundary Neumann
* nodes are ignored. The corresponding entries in the vector
* DX (which is assumed to be assembled with E030/E031/EM30/EM31) are
* multiplied by A1.
*
* In:
*   DX     - source vector
*   NMBD   - Number of edges on the real boundary
*   KMBD   - array [1..NMT] of integer
*            "Shortcut nodal property" array. This array stores all
*            edges on the boundary. The number is >= for all
*            Dirichlet edges and < 0 for all edges that should be
*            treated as Neumann edges.
*   A1     - double
*
* Out:
*   DX     - modified vector
************************************************************************

      SUBROUTINE BDRDEF (DX,KMBD,NMBD,A1)

      IMPLICIT NONE
      
C parameter

      INTEGER NMBD
      DOUBLE PRECISION DX(*),A1
      INTEGER KMBD(*)
      
C local variables

      INTEGER IMBD, IMID     

C     loop over all boundary vertices

      DO IMBD=1,NMBD
      
C       Get the shortcut nodal property

        IMID=KMBD(IMBD)

C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges

        IF (IMID.LE.0) THEN
          DX(-IMID)=A1*DX(-IMID)
        END IF
        
      END DO
        
      END

************************************************************************
* Updates the matrix entries for all DIRICHLET boundary nodes.
*
* Replaces all matrix lines corresponding to DIRICHLET nodes by
* unit vectors.
*
* In:
*   DA,
*   KCOL,
*   KLD    - System matrix, to be modified
*   NMBD   - Number of edges on the real boundary
*   KMBD   - array [1..NMT] of integer
*            "Shortcut nodal property" array. This array stores all
*            edges on the boundary. The number is >= for all
*            Dirichlet edges and < 0 for all edges that should be
*            treated as Neumann edges.
*
* Out:
*   DA     - modified matrix
************************************************************************

      SUBROUTINE BDRYA (DA,KCOL,KLD,KMBD,NMBD)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'

C parameters 

      DOUBLE PRECISION DA(*)
      INTEGER KCOL(*),KLD(*),KMBD(*),NMBD

C local variables

      INTEGER IMBD, IMID, ICOL

      DO IMBD=1,NMBD
      
C       Get the shortcut nodal property

        IMID=KMBD(IMBD)

C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges
        
        IF (IMID.GE.0) THEN

C         The number IMID is the number of the edge, similar to KMBD.
C
C         The diagonal element is set to 1, the other elements
C         are set to 0.
C         We don't use LCL1 here, since the number of elements
C         to set to 0 is usually so small, that even an LCL1 would
C         take too much time.

          DA(KLD(IMID)) = 1D0

          DO ICOL=KLD(IMID)+1,KLD(IMID+1)-1
            DA(ICOL) = 0D0
          END DO

        END IF
      END DO
      
      END

************************************************************************
* Set the DIRICHLET-components of the vector (D1,D2) to zero
*
* This is typically used to force entries in the defect vector 
* corresponding to Dirichlet nodes to zero.
*
* In:
*   D1     - array [1..*] of double
*   D2     - array [1..*] of double
*   NMBD   - Number of edges on the real boundary
*   KMBD   - array [1..NMT] of integer
*            "Shortcut nodal property" array. This array stores all
*            edges on the boundary. The number is >= for all
*            Dirichlet edges and < 0 for all edges that should be
*            treated as Neumann edges.
*
* Out:
*   D1,
*   D2     - modified vectors
************************************************************************

      SUBROUTINE    BDRY0  (D1,D2,KMBD,NMBD)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'

C parameters

      DOUBLE PRECISION D1(*),D2(*)
      INTEGER KMBD(*),NMBD

C local variables

      INTEGER IMBD, IMID

      DO IMBD=1,NMBD

C       Get the shortcut nodal property

        IMID=KMBD(IMBD)

C       Values > 0 are Dirichlet, values < 0 Neumann boundary edges

        IF (IMID.GE.0) THEN
          D1(IMID)=0D0
          D2(IMID)=0D0
        END IF
        
      END DO
      
      END

************************************************************************
* Calculate information about a boundary edge.
*
* This routine calculates the number of the boundary edge IMID, 
* the global numbers IV1, IV2 of the vertices adjacent to that edge
* as well as the nodal property of the edge (which is the number
* of the "first" adjacent vertex on the edge).
*
* In:
*   IVBD   - Index of the boundary edge in the KEBD array
*   KVBD,
*   KEBD,
*   KVERT,
*   KMID,
*   KNPR   - Usual geometry information
*
* Out:
*   IMID   - Number of the edge (range: 1..NMT)
*   INPR   - Nodal property of edhe IMID.
************************************************************************

      SUBROUTINE GETMBD (IMID,IVT1,IVT2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,
     *                   INPR)
      
      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      INCLUDE 'cout.inc'
      
C parameters

      INTEGER KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)

C local variables
      
      INTEGER IVT1, IVT2, IVBD, INPR, IEL, IMID, I, IVERT, IVERT2

C     Fetch global number of the IVBDth boundary point into IVT1
C     (see manual, P. 15)

      IVT1=KVBD(IVBD)
      INPR=KNPR(IVT1)
      IEL=KEBD(IVBD)
      IF (IEL.EQ.0) THEN
        IMID=0
        GOTO 99999
      ENDIF

C     Loop through the four vertices on the element
C     to find the "neighbour" of IVT1.

      DO I=1,4
        IVERT=I
        IF (KVERT(I,IEL).EQ.IVT1) GOTO 2
      END DO

C *** Error
      WRITE(MTERM,*) 'ERROR in GETMBD: vertice not found'
      RETURN

2     CONTINUE

      IMID=KMID(IVERT,IEL)-NVT
      IVERT2=IVERT+1
      IF (IVERT2.GT.NNVE) IVERT2=1
      IVT2=KVERT(IVERT2,IEL)

99999 END

************************************************************************
* Calculates integral boundary pressure
*
* This routine calculates the integral over the pressure on the
* boundary between the two user defined boundary segments
* [DPI(1,1)..DPI(2,1)] and [DPI(1,2)..DPI(2,2)] of the boundary
* components KPI(1) and KPI(2), respectively.
*
* In:
*   DP      - array [1..NEL] of double
*             Pressure vector, constant on each element
*   KVERT,
*   KNPR,
*   KVBD,
*   KMM,
*   DCORVG,
*   DVBDP,
*   NVBD    - Usual geometry information
*
*  In (from COMMON blocks):
*   DPI,
*   KPI     - definition of the boundary segments where to integrate
*
* Out:
*   P1,
*   P2      - Integrap boundary pressure on the two user defined
*             boundary components.
************************************************************************

      SUBROUTINE BDPRES(DP,KVERT,KNPR,KVBD,KMM,DCORVG,DVBDP,NVBD,P1,P2)

      IMPLICIT NONE
      
C main COMMON blocks

      INCLUDE 'cbasictria.inc'
      INCLUDE 'cnspts.inc'

C parameters

      DOUBLE PRECISION DP(*),DCORVG(2,*),DVBDP(*)
      INTEGER KVERT(NNVE,*),KNPR(*),KVBD(*),NVBD,KMM(2,*)
      DOUBLE PRECISION P1, P2

C local variables
      
      DOUBLE PRECISION DLEN, DPAR, PX1, PY1, PX2, PY2
      DOUBLE PRECISION DPARH, DL
      INTEGER IVBD, IVT1, IVT2, IVBDH

      P1=0D0
      DLEN=0D0

      IF (KPI(1).EQ.0) GOTO 19

C     Loop over all points on the boundary.

      DO IVBD=1,NVBD
      
C       Only look at those points on boundary component KPI(1)
C       between parameter values DPI(1,1) and DPI(2,1).
      
        DPAR=DVBDP(IVBD)
        IF ((DPAR.GE.DPI(1,1)).AND.(DPAR.LT.DPI(2,1))
     *                        .AND.(KNPR(KVBD(IVBD)).EQ.KPI(1))) THEN
          IVT1=KVBD(IVBD)
          PX1=DCORVG(1,IVT1)
          PY1=DCORVG(2,IVT1)

C         Get the "next" neighbour of vertex IVT1.
C         If IVT1 is the last point on the boundary component,
C         the "next" point is the first one.

          IF (IVT1.EQ.KMM(2,KPI(1))) THEN
            IVT2=KMM(1,KPI(1))
            PX2=DCORVG(1,IVT2)
            PY2=DCORVG(2,IVT2)
          ELSE
            IVBDH=IVBD+1
            DPARH=DVBDP(IVBDH)
            IVT2=KVBD(IVBDH)
            PX2=DCORVG(1,IVT2)
            PY2=DCORVG(2,IVT2)
          ENDIF

C         Calculate the length of the segment

          DL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
          DLEN=DLEN+DL
          
C         Simple midpoint-rule for calculating the pressure integral
          
          P1=P1+0.5D0*DL*(DP(IVT1)+DP(IVT2))
        
        ENDIF

      END DO

C     Divide by the length to get the integral value.

      IF (DLEN.NE.0) THEN
        P1=P1/DLEN
      ELSE
        P1=0D0
      END IF

C     Second user defined boundary component:
      
19    P2=0D0
      DLEN=0D0

      IF (KPI(2).EQ.0) RETURN

C     Loop over all points on the boundary.

      DO IVBD=1,NVBD
      
C       Only look at those points on boundary component KPI(1)
C       between parameter values DPI(1,1) and DPI(2,1).

        DPAR=DVBDP(IVBD)
        IF ((DPAR.GE.DPI(1,2)).AND.(DPAR.LT.DPI(2,2))
     *                        .AND.(KNPR(KVBD(IVBD)).EQ.KPI(2))) THEN
          IVT1=KVBD(IVBD)
          PX1=DCORVG(1,IVT1)
          PY1=DCORVG(2,IVT1)
          
C         Get the "next" neighbour of vertex IVT1.
C         If IVT1 is the last point on the boundary component,
C         the "next" point is the first one.

          IF (IVT1.EQ.KMM(2,KPI(2))) THEN
            IVT2=KMM(1,KPI(2))
            PX2=DCORVG(1,IVT2)
            PY2=DCORVG(2,IVT2)
          ELSE
            IVBDH=IVBD+1
            DPARH=DVBDP(IVBDH)
            IVT2=KVBD(IVBDH)
            PX2=DCORVG(1,IVT2)
            PY2=DCORVG(2,IVT2)
          ENDIF

C         Calculate the length of the segment

          DL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
          DLEN=DLEN+DL

C         Simple midpoint-rule for calculating the pressure integral

          P2=P2+0.5D0*DL*(DP(IVT1)+DP(IVT2))
        ENDIF

      END DO

C     Divide by the length to get the integral value.

      IF (DLEN.NE.0) THEN
        P2=P2/DLEN
      ELSE
        P2 = 0D0
      END IF

      END
