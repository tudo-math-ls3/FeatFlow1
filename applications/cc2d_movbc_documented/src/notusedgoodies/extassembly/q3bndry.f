************************************************************************
* Q3 boundary vertex collection routine
*
* This routine collects all DOF's on the boundary into
* one array KBDVC. The array is created on the heap, a handle is
* returned. This array can be used for the implementation of boundary
* conditions, as it represents a replacement for KVBD for all
* non-corner vertices.
*
* In:
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure of the underlying grid
* Out:
*   LBDOF : Handle to array [1..NBDVC] of integer
*           Array with the global DOF's on the boundary
*   LBDVT : Handle to array [1..NBDVC] of integer
*           Array with node number corresponding to the DOF on the
*           boundary
*   NBDVC : Number of global DOF's on the boundary.
*
* Remark: Node numbering
*  KBDOF(I) contains the number of a global DOF on the boundary.
*    This may be a corner vertex or an inner-edge vertex
*    (cp. X---x---x---X on an edge) as Q3 has 2 DOF's on one edge.
*  KBDVT(I) contains the node-number in the triangulation
*    corresponding to KBDVC(I). This may be the number of a
*    corner vertex or the number of an edge.
*  So the pair (KBDVC(I),KBDND(I)) gives an association between
*  the global DOF and the position of the DOF inside of the
*  triangulation. This is necessary as boundary-value implementation
*  routines usually check the triangulation for the type of
*  boundary conditions (e.g. "is it a Neumann edge") and don't know
*  anything about the global DOF's of the element that is
*  behind that.
************************************************************************

      SUBROUTINE Q3BDVC (TRIA,LBDOF,LBDVT,NBDVC)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      
      INTEGER TRIA(SZTRIA),LBDOF,NBDVC,LBDVT
      
C     local variables

      INTEGER KVBD,NVBD,KEBD,IELI,IEL,IELTYP,KVERT,KMID,IV,IVIDX
      INTEGER IV1,IV2,IBDVC,KBDVC,KBDND,KMBD,IEDG,NEQ,I
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
C     externals

      INTEGER NDFGX
      EXTERNAL NDFGX
      
C     This routine only works for Q3 elements:

      IELTYP = 14
      
C     How many boundary vertices to we have?
C     = Number of corners on the boundary + 2*(number of edges).
C     Where: Number of edges = Number of corner vertices on the boundary.

      NVBD = TRIA(ONVBD)
      NBDVC = NVBD + 2*NVBD
      
C     Get the number of global DOF's

      NEQ = NDFGX(IELTYP,TRIA)
      
C     Allocate memory

      CALL ZNEW(NBDVC,3,LBDOF,'KBDVC ')
      CALL ZNEW(NBDVC,3,LBDVT,'KBDND ')
      
      KBDVC = L(LBDOF)
      KBDND = L(LBDVT)
      
C     Copy KVBD to KBDVC. This is the first part of boundary vertices.
      
      CALL LCP3(KWORK(L(TRIA(OLVBD))),KWORK(KBDVC),NVBD)
      
C     The first NVBD vertex numbers coincide with the DOF's:

      CALL LCP3(KWORK(L(TRIA(OLVBD))),KWORK(KBDND),NVBD)
      
C     IBDVC counts how many nodes in KBDVC are used already.
      
      IBDVC = TRIA(ONVBD)
      IDFL = 16
      
C     Now we have to collect the inner-edge vertices.
C     For this, we loop through the elements on the boundary.

      KVBD  = L(TRIA(OLVBD))
      KMBD  = L(TRIA(OLMBD))
      KEBD  = L(TRIA(OLEBD))
      KVERT = L(TRIA(OLVERT))
      KMID  = L(TRIA(OLMID))
      
      DO IELI = 1,NVBD

C       We have two DOF's per edge which do not coincide with any 
C       points on the edge in the triangulation! Q3 is numbered 
C       locally as:
C
C           4 --10---9-- 3
C           |            |
C          11   16  15   8
C           |            |
C          12   13  14   7
C           |            |
C           1 ---5---6-- 2
C
C       We are at boundary-element index IELI - the real element
C       number of the boundary element can be obtained with KEBD:

        IEL = KWORK(KEBD+IELI-1)
        
C       and the corresponding vertex before the edge on the 
C       boundary is:
        
        IV = KWORK(KVBD+IELI-1)
        
C       Also fetch the edge number:

        IEDG = KWORK(KMBD+IELI-1)
        
C       Get the local and global DOF's on this element - sorted.

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KWORK(KVERT),KWORK(KMID),
     *              KDFG,KDFL)
     
C       Sort for the local DOF

        CALL NGLSD(KDFL,KDFG,IDFL)
     
C       The corner vertex IV corresponds to the global DOF!
C       Search in KDFG for IV to find its local DOF.

        DO IVIDX = 1,IDFL
          IF (KDFG(IVIDX).EQ.IV) GOTO 10
        END DO
        WRITE (*,*) 'Q3BDVC: Vertex not found!'
        STOP
10      CONTINUE

C       The DOF's are sorted, so KFGL(I)=I for I=1,NDFL!
C       We need the vertices on the edge following IV.
C       The local DOF's are simply calculated from
            
        IV1 = 5 + 2*(IVIDX-1)
        IV2 = 5 + 2*(IVIDX-1)+1
          
C       (see the picture above). Add the corresponding global DOF's
C       to the KBDVC-array.

        KWORK(KBDVC+IBDVC)   = KDFG(IV1)
        KWORK(KBDVC+IBDVC+1) = KDFG(IV2)
        
C       Associate these boundary nodes with the edge in the 
C       triangulation.

        KWORK(KBDND+IBDVC)   = IEDG
        KWORK(KBDND+IBDVC+1) = IEDG
        
        IBDVC = IBDVC+2
        
      END DO
        
C     That's it.
        
      END

