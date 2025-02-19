      BLOCK DATA ZVALUE1
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*1
      PARAMETER (NNARR=299,NNLEV=9)
      PARAMETER (NBLOCA=1,NBLOCB=2,NBLOCF=2)
      PARAMETER (NBLA1=NBLOCA*NNLEV,NBLB1=NBLOCB*NNLEV,
     *           NBLF1=NBLOCF*NNLEV)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM(120)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
      DATA M/2/,MT/2/,MKEYB/5/,MTERM/6/,IER/0/,ICHECK/1/
      DATA MERR/11/,MPROT/12/,MSYS/13/,MTRC/14/,IRECL8/512/
      DATA SUB/'MAIN  '/
      DATA FMT/'(3D25.16)','(5E16.7)','(6I12)'/
      DATA CPARAM/120*' '/
      DATA NEL/0/,NVT/0/,NET/0/,NAT/0/,NVE/0/,NEE/0/,NAE/0/,NVEL/0/,
     *     NEEL/0/,NVED/0/,NVAR/0/,NEAR/0/,NBCT/0/,NVBD/0/,NEBD/0/,
     *     NABD/0/
      DATA LCORVG/0/,LCORMG/0/,LCORAG/0/,LVERT/0/,LEDGE/0/,LAREA/0/,
     *     LADJ/0/,LVEL/0/,LEEL/0/,LAEL/0/,LVED/0/,LAED/0/,LVAR/0/,
     *     LEAR/0/,LEVE/0/,LAVE/0/,LNPR/0/,LBCT/0/,LVBD/0/,LEBD/0/,
     *     LABD/0/
      DATA KTYPE/NNARR*0/,KLEN/NNARR*0/,KLEN8/NNARR*0/,IFLAG/0/
      DATA NMAEL/0/,NMAVT/0/,NMAEDG/0/,NMAVE/0/,NMAVEL/0/,NMABCT/0/,
     *     NMAVBD/0/
      DATA LMACVG/0/,LMACMG/0/,LMAVT/0/,LMAMID/0/,LMAADJ/0/,LMAVEL/0/,
     *     LMAMEL/0/,LMANPR/0/,LMAMM/0/,LMAVBD/0/,LMAEBD/0/,LMABCT/0/,
     *     LMAVBP/0/,LMAMBP/0/,LMAVE/0/
      DATA KNEL/NNLEV*0/,KNVT/NNLEV*0/,KNET/NNLEV*0/,
     *     KNAT/NNLEV*0/,KNVE/NNLEV*0/,KNEE/NNLEV*0/,    
     *     KNAE/NNLEV*0/,KNVEL/NNLEV*0/,KNEEL/NNLEV*0/, 
     *     KNVED/NNLEV*0/,KNVAR/NNLEV*0/,KNEAR/NNLEV*0/,
     *     KNBCT/NNLEV*0/,KNVBD/NNLEV*0/,KNEBD/NNLEV*0/, 
     *     KNABD/NNLEV*0/          
      DATA KLCVG/NNLEV*0/,KLCMG/NNLEV*0/,KLCAG/NNLEV*0/,
     *     KLVERT/NNLEV*0/,KLEDGE/NNLEV*0/,KLAREA/NNLEV*0/,
     *     KLADJ/NNLEV*0/,KLVEL/NNLEV*0/,KLEEL/NNLEV*0/,
     *     KLAEL/NNLEV*0/,KLVED/NNLEV*0/,KLAED/NNLEV*0/,
     *     KLVAR/NNLEV*0/,KLEAR/NNLEV*0/,KLEVE/NNLEV*0/,
     *     KLAVE/NNLEV*0/,KLNPR/NNLEV*0/,KLBCT/NNLEV*0/,
     *     KLVBD/NNLEV*0/,KLEBD/NNLEV*0/,KLABD/NNLEV*0/
      DATA ILEV/0/,NLEV/0/,NLMIN/0/,NLMAX/0/,
     *     ICYCLE/0/,KPRSM/NNLEV*0/,KPOSM/NNLEV*0/
      DATA KLA/NNLEV*0/,KLST/NBLA1*0/,KLMASS/NBLA1*0/,KLM/NBLA1*0/,
     *     KLCOLA/NNLEV*0/,KLLDA/NNLEV*0/,
     *     KLB1/NNLEV*0/,KLB2/NNLEV*0/,KLB3/NNLEV*0/
     *     KLCOLB/NNLEV*0/,KLLDB/NNLEV*0/,
     *     KLUP/NNLEV*0/,KLF12P/NNLEV*0/,KLAUX/NNLEV*0/,KLVOL/NNLEV*0/,
     *     LU1OLD/0/,LU2OLD/0/,LU3OLD/0/,
     *     LPOLD/0/,LD1/0/,LD2/0/,LD3/0/,LDP/0/
      DATA KNA/NNLEV*0/,KNB/NNLEV*0/,KNU/NNLEV*0/,KNP/NNLEV*0/,
     *     KNUP/NNLEV*0/
      DATA KLC/NNLEV*0/,KLCOLC/NNLEV*0/,KLLDC/NNLEV*0/,KNC/NNLEV*0/
      DATA KLAILU/NNLEV*0/,KLTRA1/NNLEV*0/,KLTRA2/NNLEV*0/
      DATA KLCILU/NNLEV*0/,KLTRC1/NNLEV*0/,KLTRC2/NNLEV*0/
      DATA INEUM/0/,LELBD/0/,KELBD/NNLEV*0/,KLNPRO/NNLEV*0/
      DATA NA/0/,NB/0/,NU/0/,NP/0/,NUP/0/
      DATA KA1/0/,KST1/0/,KMASS1/0/,KM1/0/,KCOLA/0/,KLDA/0/,KB1/0/,
     *     KB2/0/,KB3/0/,KCOLB/0/,KLDB/0/,KU1/0/,KU2/0/,
     *     KU3/0/,KP/0/,
     *     KF1/0/,KF2/0/,KF3/0/,KFP/0/,KAUX1/0/,KAUX2/0/,KAUX3/0/,
     *     KAUXP/0/
C
      END
      
