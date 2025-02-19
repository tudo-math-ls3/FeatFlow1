************************************************************************
* This file contains COMMON blocks important for describing the basic
* structures in multigrid.
*
* The file CBASICMG.FOR must be included prior to this file!
************************************************************************

      DOUBLE PRECISION TTMG,TTS,TTE,TTD,TTP,TTR
      INTEGER IMTIME
      INTEGER KLMBD,KLDBD,KNMBD
      INTEGER KLNPRO,INEUM
      INTEGER KNA,KNB,KNU,KNP
      INTEGER KNUP

      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)

      SAVE /MGTIME/,/MGBDRY/,/MGDIM/

C KNA    - Number of entries in each Stokes-matrix-part
C KNB    - Number of entries in each B-matrix
C KNU    - Number of unknowns per velocity; this is normally the same
C          as KNEQA, except for that KNEQ is build by the matrix building
C          routine whereas KNU is computed manually
C KNP    - Number of unknowns in pressure;; this is normally the same
C          as KNEQB, except for that KNEQ is build by the matrix building
C          routine whereas KNP is computed manually
C KNUP   - Length of whole solution vector: 2*velocity+pressure
