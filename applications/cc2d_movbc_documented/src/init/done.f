***********************************************************************
* The following code is used for deallocation of allocated vectors.
* The first part cleans up all initializations that are performed
* by INIT1.
***********************************************************************

C ---------------------------------------------------------------------
C Clean up the parametrization of the domain (->RDPARM)
C ---------------------------------------------------------------------

      SUBROUTINE DISPAR 
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cparqdata.inc'
      
      CALL ZDISP(0,LNCOMP,'NCOMP ')
      CALL ZDISP(0,LICPTR,'ICPTR ')
      CALL ZDISP(0,LITYP,'ITYP  ')
      CALL ZDISP(0,LNSPLN,'NSPLN ')
      CALL ZDISP(0,LNPAR,'NPAR  ')
      CALL ZDISP(0,LIPPTR,'IPPTR ')
      CALL ZDISP(0,LXPAR,'DXPAR ')
      CALL ZDISP(0,LYPAR,'DYPAR ')

      END

C ---------------------------------------------------------------------
C Clean up the triangulation (-> XMSB2)
C ---------------------------------------------------------------------

      SUBROUTINE DISORS ()
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmg.inc'

      INCLUDE 'cinidat.inc'

      INTEGER I,J
      
      DO I=NLMAX,NLMIN,-1
      
        CALL ZDISP(0,KLCMG(I) ,'KLCMG'); 
        CALL ZDISP(0,KLVERT(I),'KLVERT');
        IF (KLMID(I).NE.0) CALL ZDISP(0,KLMID(I) ,'KLMID');
        IF (KLADJ(I).NE.0) CALL ZDISP(0,KLADJ(I) ,'KLADJ');
        CALL ZDISP(0,KLVEL(I) ,'KLVEL');
        IF (KLMEL(I).NE.0) CALL ZDISP(0,KLMEL(I) ,'KLMEL');
        CALL ZDISP(0,KLNPR(I) ,'KLNPR');
        CALL ZDISP(0,KLMM(I)  ,'KLMM');
        CALL ZDISP(0,KLVBD(I) ,'KLVBD');
        CALL ZDISP(0,KLEBD(I) ,'KLEBD');
        CALL ZDISP(0,KLBCT(I) ,'KLBCT');
        CALL ZDISP(0,KLVBDP(I),'KLVBDP');
        CALL ZDISP(0,KLMBDP(I),'KLMBDP');

      END DO

C The grid coordinates in DCORVG only exist once!
C The handles in KLCVG(...) are (normally) all the same!
C
C So we have to make sure every handle is only released once.
C All other KLCVG-entries are set to 0 like the ZDISP function
C would do.

      DO I=NLMAX,1,-1
      
        IF (KLCVG(I).NE.0) THEN

          DO J=I-1,1,-1
            IF (KLCVG(J).EQ.KLCVG(I)) KLCVG(J)=0
          END DO

          CALL ZDISP(0,KLCVG(I) ,'KLCVG');
          
        END IF
      
      END DO
      
      END

C ---------------------------------------------------------------------
C Clean up backup of KNPR if necessary
C (in case of time dependent boundary components)
C ---------------------------------------------------------------------

      SUBROUTINE DISKNP ()
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmg.inc'

      INCLUDE 'cinidat.inc'

      INTEGER I
      
      INTEGER NFBDYC
      EXTERNAL NFBDYC
      
      IF ((IBDR.GE.2).OR.(NFBDYC().GT.0)) THEN
        DO I=NLMAX,NLMIN,-1
          CALL ZDISP(0,KLNPRO(I),'KNPRO ')
        END DO
      END IF

      END


C ---------------------------------------------------------------------
C Clean up the structures describing the midpoints on the
C boundary:
C   KMBD=no. of midpoint on boundary
C   DDBD=parameters of midpoints on boundary
C ---------------------------------------------------------------------

      SUBROUTINE DISMBD ()
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'cinidat.inc'

      INTEGER I

      DO I=NLMIN,NLMAX

        CALL ZDISP(0,KLMBD(I),'KMBD  ')
        CALL ZDISP(0,KLDBD(I),'DDBD  ')

      END DO
      
      END 

C ---------------------------------------------------------------------
C Clean up the matrix structure of the system matrices
C (--> normally shared between Stokes matrices, mass matrices,...)
C ---------------------------------------------------------------------

      SUBROUTINE DISMTS
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgpar.inc'
      
      INTEGER I

C *** Release of matrix structure

      DO I=NLMAX,NLMIN,-1
        CALL ZDISP (0,KLCOLA(I),'KLCOL ')
        CALL ZDISP (0,KLLDA(I),'KLLD  ')
      END DO

      END 


C ---------------------------------------------------------------------
C Clean up the Stokes matrix content
C
C If ILV=0, all matrices will be released, otherwise only the matrix
C on level ILV
C ---------------------------------------------------------------------

      SUBROUTINE DISSTM (ILV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgpar.inc'
      
      INTEGER I,ILV

C *** Release stokes matrix

      IF (ILV.EQ.0) THEN
        DO I=NLMAX,NLMIN,-1
          CALL ZDISP (0,KLST(I),'KLST  ')
        END DO
      ELSE
        CALL ZDISP (0,KLST(ILV),'KLST  ')
      END IF

      END 
      
      
C ---------------------------------------------------------------------
C Clean up the mass matrices
C
C If ILV=0, all matrices will be released, otherwise only the matrix
C on level ILV
C ---------------------------------------------------------------------

      SUBROUTINE DISMM (ILV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INTEGER I,ILV

C *** Release mass matrices

      IF (ILV.EQ.0) THEN
        DO I=NLMAX,NLMIN,-1
          CALL ZDISP (0,KLM(I),'KLM   ')
        END DO
      ELSE
        CALL ZDISP (0,KLM(ILV),'KLM   ')
      END IF
        
      END 


C ---------------------------------------------------------------------
C Clean up the matrix blocks B1,B2
C
C If ILV=0, all matrices will be released, otherwise only the matrix
C on level ILV
C ---------------------------------------------------------------------

      SUBROUTINE DISMTB (ILV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INTEGER I,ILV

C *** Release matrix blocks B1,B2
 
      IF (ILV.EQ.0) THEN
        DO I=NLMAX,NLMIN,-1
          CALL ZDISP (0,KLB2(I),'KLB2  ')
          CALL ZDISP (0,KLB1(I),'KLB1  ')
          CALL ZDISP (0,KLCOLB(I),'KCOLB ')
          CALL ZDISP (0,KLLDB (I),'KLLDB ')
        END DO
      ELSE 
        CALL ZDISP (0,KLB2(ILV),'KLB2  ')
        CALL ZDISP (0,KLB1(ILV),'KLB1  ')
        CALL ZDISP (0,KLCOLB(ILV),'KCOLB ')
        CALL ZDISP (0,KLLDB (ILV),'KLLDB ')
      END IF
 
      END 
      

C ---------------------------------------------------------------------
C Clean up the solution vectors and vectors of right hand side.
C
C In:
C  BRHS   - release right hand side handles
C  BSOL   - release solution vector handles 
C  BSLD   - release defect correction vector handles
C  BSLOLD - release handles of backup of defect correction vector 
C ---------------------------------------------------------------------

      SUBROUTINE DISSOL (BRHS, BSOL, BSLD, BSLOLD)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cinidat.inc'
      
      LOGICAL BRHS, BSOL, BSLD, BSLOLD
      INTEGER I

C *** Release solution vectors, UOLD on NLMAX

      DO I=NLMAX,NLMIN,-1
        IF (BRHS) CALL ZDISP(0,KLF12P(I),'DF12P ')
        IF (BSOL) CALL ZDISP(0,KLUP(I),'KLUP  ')
      END DO

      I=2
      ILEV=NLMAX
      CALL SETLEV(I)

      IF (BSLOLD) THEN
        CALL ZDISP(0,LU1OLD,'DU1OLD')
      END IF
      
      IF (BSLD) THEN
        CALL ZDISP(0,LD1,'DD12P ')
      END IF
      
      END 
      

C ---------------------------------------------------------------------
C Clean up the system matrix content
C ---------------------------------------------------------------------

      SUBROUTINE DISMTA
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INTEGER I

C *** Release system matrix

      DO I=NLMAX,NLMIN,-1
        CALL ZDISP(0,KLA(I),'KLA   ')
      END DO

      END 


C ---------------------------------------------------------------------
C Clean up auxiliary vectors
C ---------------------------------------------------------------------

      SUBROUTINE DISAUX
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INTEGER I

C *** Release auxil. vectors

      DO I=NLMAX,NLMIN,-1
        CALL ZDISP(0,KLAUX(I),'KLAUX ')
      END DO

      END 

C ---------------------------------------------------------------------
C Clean up vectors with area of elements
C ---------------------------------------------------------------------

      SUBROUTINE DISARE
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INTEGER I

C *** Release vector of element areas and auxil. vectors

      DO I=NLMAX,NLMIN,-1
        CALL ZDISP(0,KLAREA(I),'KLA   ')
      END DO

      END 

C ---------------------------------------------------------------------
C Clean up current level information
C (I.e. all handles,... that are still stored in /TRIAA/, /TRIAD/,...
C
C This is done by a call to SETLEV for setting the current level to 1,
C so before this routine is called, all other information in the 
C MG common blocks must have been cleared!
C ---------------------------------------------------------------------

      SUBROUTINE DISCLV ()
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INTEGER I
      
C Set the current level to NLMIN to delete any informations in the
C "current level" block

      I = 2
      ILEV = NLMIN
      
      CALL SETLEV (I)

      END

***********************************************************************
* DONE1
*
* Deallocates all vectors that have been allocated during INIT1.
* Sets the current multigrid level to 1 and clears information abount
* the current level as far as they have to be cleared.
*
* The following constants can be specified in a bitmask INTDIS to 
* prevent the deallocation of certain structures:
*
*       1 - Don't release system matrix A (only content)
*       2 - Don't release system matrix B (only content)
*       4 - Don't release mass matrices (if allocated)
*       8 - Don't release Stokes matrices S1,S2 (KLST)
*      16 - Don't release matrix structure of system matrix A
*      32 - Don't release boundary midp. structures KLMBD,DDBD
*      64 - Don't release solution vectors
*     128 - Don't release auxiliary vectors 
*     256 - Don't release vector with area of elements
*     512 - Don't release multigrid triangulation
*    1024 - Don't release backup of KNPR-arrays (if allocated)
*    2048 - Don't release domain parametrisation
*
* So:
*      CALL DONE1(0) clears everything
*      CALL DONE1(65535) clears nothing but set's the current MG-level
*                        to 1
***********************************************************************

      SUBROUTINE DONE1 (INTDIS)

      IMPLICIT NONE
      
      INCLUDE 'cnsparfrac.inc'
      INCLUDE 'cinidat.inc'
      INCLUDE 'cfiles.inc'
      
      INTEGER INTDIS
      
C Release vectors allocated in the standard implementation of CC2D:      
      
      IF (IAND(INTDIS,    1).EQ.0) CALL DISMTA
      IF ((IAND(INTDIS,    2).EQ.0).AND. 
     *             (IPRECB.NE.2)) CALL DISMTB (0)
C Release lumped or real mass matrices
      IF ((IAND(INTDIS,    4).EQ.0).AND.
     *              (ISTAT.NE.0)) CALL DISMM (0)
      IF (IAND(INTDIS,    8).EQ.0) CALL DISSTM (0)
      IF (IAND(INTDIS,   16).EQ.0) CALL DISMTS
      IF (IAND(INTDIS,   32).EQ.0) CALL DISMBD
      IF (IAND(INTDIS,   64).EQ.0) THEN
C Memory of backup vectors only has to be released if allocated:
        CALL DISSOL (.TRUE.,.TRUE.,.TRUE.,
     *               (OMGMIN.GT.0D0).OR.(OMGMAX.GT.0D0))
      END IF
      IF (IAND(INTDIS,  128).EQ.0) CALL DISAUX
      IF (IAND(INTDIS,  256).EQ.0) CALL DISARE
      IF (IAND(INTDIS,  512).EQ.0) CALL DISORS
      IF (IAND(INTDIS, 1024).EQ.0) CALL DISKNP
      IF (IAND(INTDIS, 2048).EQ.0) THEN
        IF (IMESH1.GT.0) CALL DISPAR
      END IF

C Switch to level 1 to update information in the Common blocks
C descibing the current level

      CALL DISCLV
      
      END
 
 
      
