************************************************************************
* This file contains the variable representing the user provided
* information in the DAT file of the application.
************************************************************************

C Access each variable in the following COMMON block separately or 
C access all variables in this COMMON block at once (e.g. for quick
C saving for future use)

      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   

C Array VIPARM can be used to make a snapshot
C of the whole configuration

      INTEGER IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *        IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *        INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *        NSM,NSL,NSMFAC,IRMESH

      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IRMESH

C Double precision parameters

      DOUBLE PRECISION VRPARM
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          

C Array VRPARM can be used to make a
C snapshot of the whole configuration.

      DOUBLE PRECISION NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                 EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                 RLXSM,RLXSL,AMINMG,AMAXMG

      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG

      SAVE /IPARM/, /RPARM/
      