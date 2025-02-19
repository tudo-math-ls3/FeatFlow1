************************************************************************
* This file contains additional COMMON blocks with multigrid data.
* These are mostly auxiliary vectors and variables denoting information
* on the current level.
************************************************************************

      INTEGER NA,NB,NU,NP,NUP
      INTEGER KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB
      INTEGER KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP

      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP

      SAVE /LEVDIM/,/ADRFLD/
      