***********************************************************************
* This file provides drag-/lift calculation routines like the file
* draglift.f. In contrast to draglift.f the routines here
* build up upon the new calculation method with the general
* volume integration method VOLINT in bdvol.f.
***********************************************************************

***********************************************************************
* Initialise ALPHA-vector
*
* This routine prepares the ALPHA-vector DALPHA such that it is 1
* on all entries belonging to nodes in the (fictitious) boundary
* component IFBC and 0 everywhere else. If IFBC=0, ALPHA is prepared
* such that it is 1 on all nodes belonging to any fictitious boundary
* component.
* This routine is designed to work with midpoints, which are
* nodes if EM30/EM31 is used. In this case there must be NU=NMT.
*
* In:
*  DALPHA  - array [1..NU] of double; ALPHA-vector to be filled
*  KNPR,
*  KMID    - information from the triangulation
*  IFBC    - Number of boundary component to prepare the ALPHA
*            vector for:
*            0 = all fictitious boundary components
*            NBCT+1..NBCT+NFBDYC: prepare for fictitious boundary
*                component IFBC-NBCT (=1..NFBDYC)
*
* Out:
*  DALPHA  - the ALPHA-vector filled with 1/0
***********************************************************************

      SUBROUTINE INALPH (DALPHA,KNPR,KMID,NEL,NVT,NBCT,NU,IFBC)

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C parameters

      INTEGER KNPR(*),KMID(NNVE,*),IFBC,NU,NEL,NVT,NBCT
      DOUBLE PRECISION DALPHA(NU)

C local variables

      INTEGER IEL0,IVE,IM1,IVT
      
      CALL  LCL1(DALPHA, NU) 
      
      DO IEL0=1,NEL
        DO IVE=1,4
          
          IM1=KMID(IVE,IEL0)
          
          IF (IFBC.EQ.0) THEN
C Handle the case of all fictitious boundaries:
            IF(KNPR(IM1).LT.-NEL) THEN
              DALPHA(IM1-NVT) = 1.0D0
            END IF
          ELSE
            IF ((IFBC.GT.NBCT).AND.(KNPR(IM1).EQ.-(NEL+IFBC))) THEN
C A specific fictitious boundary:
              IVT=KNPR(IM1)
              DALPHA(IM1-NVT) = 1.0D0
            ELSE IF (KNPR(IM1).EQ.IFBC) THEN
C Undocumented, not tested: non-fictitious boundary component IFBC>0
              IVT=KNPR(IM1)
              DALPHA(IM1-NVT) = 1.0D0
            ENDIF
          ENDIF
       
        ENDDO
      ENDDO

      END      

