************************************************************************
* This file contains COMMON blocks holding the handles for the
* triangulations that are used in multigrid.
*
* The file CBASICMG.FOR must be included prior to this file!
************************************************************************

      INTEGER KNEL,KNVT,KNMT
      INTEGER KNVEL,KNVBD
      INTEGER KLCVG,KLCMG,KLVERT
      INTEGER KLMID,KLADJ,KLVEL
      INTEGER KLMEL,KLNPR,KLMM 
      INTEGER KLVBD,KLEBD,KLBCT
      INTEGER KLVBDP,KLMBDP

      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)

      SAVE /MGTRD/,/MGTRA/
