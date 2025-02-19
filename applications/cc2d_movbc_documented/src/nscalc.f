************************************************************************
* (Re)initialise the Navier Stokes solver for a new run.
*
* In:
*  BPRT  - false=perform a complete reinitialisation
*          true =perform a partial reinitialisation, don't destroy
*                current solution vectors
*
* Generally the routine performs the following:
*  - regeneration of the mass matrices, Stokes matrices, RHS
************************************************************************

      SUBROUTINE ININSO (BPRT)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cinigeometry.inc'
      
C parameters

      LOGICAL BPRT
      
C local variables

      INTEGER ILV
      INTEGER LF1,LF2,LF12P

C Reinit all Navier-Stokes structures.

C Release the mass-/Stokes matrix entries,...

      IF (ISTAT.NE.0) THEN
        CALL DISMM (0)
      END IF
      
      IF (IPRECA.NE.4) THEN
        CALL DISSTM (0)
      ENDIF
      
      CALL DISMTB (0)

C At this point it's able to include grid deformation...
C This will be included in a later FeatFlow version.
C For now we continue rebuilding everything necessary.

C Restore the KNPR-arrays and rebuild them. 

      DO ILV=NLMIN,NLMAX
        CALL XLCP3(KLNPRO(ILV),KLNPR(ILV),KNVT(ILV)+KNMT(ILV))
        CALL IMPBDG (ILV)
      END DO
        
C Rebuild the mass/Stokes matrices, depending on our
C current grid

      IF (IPRECA.NE.4) THEN
        CALL GENSTM 
      ENDIF

      IF (ISTAT.NE.0) THEN
        CALL GENMM (IMASS.EQ.0,IMASSL.NE.0)
      END IF

      CALL GENMTB(IPRECB.GE.2)
          
C Rebuild right hand side like in INIT1 for finest level

      CALL GENRHS(NLMAX,LF1,LF2)

      LF12P = KLF12P(NLMAX)
      CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LF12P)),KNEQA(NLMAX))
      CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LF12P)+KNEQA(NLMAX)), 
     *                  KNEQA(NLMAX))
      CALL  ZDISP (0, LF1, 'LF1TMP' )
      CALL  ZDISP (0, LF2, 'LF2TMP' )

C Immediately implement boundary Dirichlet/Neumann conditions.

      CALL IMPSLR (NLMAX,(IBDR.EQ.1).AND.(ISTAT.EQ.0),
     *                 .TRUE.,.FALSE.)
     
C Reinitialise the current solution vectors if necessary

      IF (.NOT.BPRT) THEN
        CALL LCL1 (DWORK(L(KLUP(NLMAX))),KNUP(NLMAX))
      END IF
     
C Implement boundary conditions into solution. This is always necessary,
C as for a changed geometry we have at least different Dirichlet solution
C values. Although this would not be necessary in the first iteration
C or if the solution vectors contain precalculated solutions,
C we do this here for "safetyness" and "simplicity", since it's quickly
C be done...

      CALL IMPSLR (NLMAX,.FALSE.,.FALSE.,.TRUE.)

      END
            