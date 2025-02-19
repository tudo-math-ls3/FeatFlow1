************************************************************************
* This file contains the multigrid COMMON block /PAR/ 
* that are very often used in multigrid routines. These describe
* masic multigrid parameters, matrices, right hand sides,...
*
* The file CBASICMG.FOR must be included prior to this file!
************************************************************************

      INTEGER ILEV,NLEV,NLMIN,NLMAX
      INTEGER ICYCLE,KPRSM,KPOSM

      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)

      SAVE /MGPAR/
