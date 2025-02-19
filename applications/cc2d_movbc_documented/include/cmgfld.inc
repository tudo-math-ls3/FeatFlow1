************************************************************************
* This file contains the multigrid COMMON block /FLD/
* that are very often used in multigrid routines. These describe
* masic multigrid parameters, matrices, right hand sides,...
*
* The file CBASICMG.FOR must be included prior to this file!
************************************************************************

      INTEGER KLA,KLST,KLM
      INTEGER KLCOLA,KLLDA
      INTEGER KLB1,KLB2
      INTEGER KLCOLB,KLLDB
      INTEGER KLUP,KLF12P,KLAUX
      INTEGER KLAREA,LU1OLD,LD1
      INTEGER KNEQA,KNEQB

      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),
     *                KNEQA(NNLEV),KNEQB(NNLEV),
     *                LU1OLD,LD1

C KNEQA describes the number of equations in matrix A;
C this is normally the same as KNU, but build in the matrix
C building routine.
C KNEQB describes the number of equations in matrix A;
C this is normally the same as KNP, but build in the matrix
C building routine.

      SAVE /MGFLD/
