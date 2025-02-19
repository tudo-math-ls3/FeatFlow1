************************************************************************
* This file contains additional COMMON blocks for the Navier Stokes
* solver - to be more exact: the main parameters about the time 
* discretization.
************************************************************************

      DOUBLE PRECISION TSTEP,THETA,THSTEP,TIMENS,EPSNS
      INTEGER NITNS,ITNS,ISTAT

      DOUBLE PRECISION THETAP,FALPHA,FBETA
      INTEGER IFRSTP

      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP

      SAVE /NSPAR/,/NSFRAC/
      