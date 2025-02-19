************************************************************************
c
c     Unterprogramme zum Themenkreis
c     Anisotrope Gitterverfeinerung
c
************************************************************************

************************************************************************
      SUBROUTINE  NEUTST (KMBD,KVBD,KEBD,KVERT,KMID,KNPR,DDBD,DMBDP,
     *                    NMBD)
************************************************************************
C
C-----------------------------------------------------------------------
*    Purpose:  sets the DIRICHLET- and NEUMANN components
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      DIMENSION KMBD(*),KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*)
      DIMENSION KNPR(*),DDBD(*),DMBDP(*)
C
C
C
      NMBD =0
      INEUM=0
      DO 1 IVBD=1,NVBD
      CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
      IF (IMID.EQ.0) GOTO 1
C
      DPAR=DMBDP(IVBD)
c      NMBD=NMBD+1
c      KMBD(NMBD)=IMID
c      DDBD(NMBD)=DPAR
C
      INPART=0
      CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
      NPART=INPART
C      
      DO 10 INPART=1,NPART
      CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
      IF ((DPAR.GT.DPARN1).AND.(DPAR.LT.DPARN2)
     *                    .AND.(INPR.EQ.INPRN)) THEN
c       INEUM=1
       KNPR(NVT+IMID)=0
c       KMBD(NMBD)=-KMBD(NMBD)
      ENDIF
10    CONTINUE
C
1     CONTINUE
C
      END

************************************************************************
      SUBROUTINE REFANI (DCORVG,KVERT,KWAND,NELS,i1,i2,i3,i4)
*
************************************************************************
c    adjusts the 'Wall-elements' of the regularly refined grid
c    according to a finer resolution on the wall.


       INCLUDE 'common.inc'
       INCLUDE 'dwork.inc'
       INCLUDE 'bouss.inc'
       INCLUDE 'block.inc'
 
      DIMENSION x(4),y(4),INP(4)
      DIMENSION DCORVG(2,*), KVERT(4,*), KWAND (*)
      DOUBLE PRECISION LAENGE, HOEHE,OFFSET,Parame
      EXTERNAL   PARX,PARY,TMAX
      SAVE
c
C
         CALL ZNEW (NELS,3,LCHECK,'LCHECK')
         CALL ZNEW (NVT ,3,LCHCK2,'LCHCK2')
c
c===4=============================================
c$$$c     Benchmark:
c$$$c     Dimensionen
c$$$      Laenge=1.1d0
c$$$      HOEHE =0.17d0
c$$$      Parame=1100d0
c$$$      OFFSET=-0.3d0
c$$$c     
c$$$c     Eckkenparameter
c$$$      TECK1 = 1d0
c$$$      TECK2 = 2d0
c$$$      TECK3 = 3d0
c$$$      TECK4 = 4d0
c$$$c
c$$$c     T1-T4 
c$$$      T1=1185d0
c$$$      T2=1185d0
c$$$      T3=2455D0
c$$$      T4=2455D0
c================================================
c     Einheitsquadrat:
c     Dimensionen
      Laenge=1d0
      Parame=1d0
      OFFSET=-0.0d0
      HOEHE =1.0d0
c     
c     Eckkenparameter
      TECK1 = 1d0
      TECK2 = 2d0
      TECK3 = 3d0
      TECK4 = 4d0
c
c      mittlere Randparameter T1-T4 
      T1=1.5d0
      T2=1.5d0
      T3=3.5D0
      T4=3.5D0
c================================================

         DO 10 IEL=1,NELS
            IELNR=abs(KWAND(IEL))
            IF (KWORK(L(LCHECK)-1+IELNR).NE.0) GOTO 10
            KWORK(L(LCHECK)-1+IELNR)=1
C

            INP(1)=KVERT(1,IELNR)
            INP(2)=KVERT(2,IELNR)
            INP(3)=KVERT(3,IELNR)
            INP(4)=KVERT(4,IELNR)
c

            X(1)=DCORVG(1,INP(1))
            Y(1)=DCORVG(2,INP(1))
            X(2)=DCORVG(1,INP(2))
            Y(2)=DCORVG(2,INP(2))
            X(3)=DCORVG(1,INP(3))
            Y(3)=DCORVG(2,INP(3))
            X(4)=DCORVG(1,INP(4))
            Y(4)=DCORVG(2,INP(4))
c
            YMAX=MAX(Y(1),Y(2),Y(3),Y(4))
            YMIN=MIN(Y(1),Y(2),Y(3),Y(4))
c

            IF (YMAX.EQ.HOEHE) THEN
               DO 11 J=1,4
                 IF (Y(J).NE.HOEHE) then
                    IF (KWORK(L(LCHCK2)-1+INP(j)).EQ.0)
     *              DCORVG(2,INP(J))=Y(J)+(HOEHE-Y(J))*OME1
                 endif

 11            CONTINUE
            ELSE
               BBOTTN=.FALSE.
               DO 22 J=1,4
                  xjpar=(((X(J)-OFFSET)/Laenge))*Parame
                  IF (abs(Y(J)-PARY(xjpar)).le.1d-9)
     *            BBOTTN=.TRUE.
 22               CONTINUE
               IF (BBOTTN) THEN
c         wenn das Element in der Schicht am Boden ist,dann
c         verschiebe die Knoten, die NICHT an der Wand liegen!
c
               
                  DO 24 J=1,4
                     XJPAR=(((X(J)-OFFSET)/Laenge))*Parame
c
                     IF (abs(y(j)-PARY(XJPAR)).le.1d-9)
     *                  KWORK(L(LCHCK2)-1+INP(j))=1
c        dann liegt der Knoten auf dem Rand
c
                     IF (KWORK(L(LCHCK2)-1+INP(j)).EQ.0) THEN
c
                        XPARP=(((X(MOD(4,J)+1)-OFFSET)/Laenge))*Parame
                        IF (abs(Y(MOD(4,J)+1)-PARY(XPARP)).le.1d-9) THEN
                           J2=(MOD(4,J)+1)
                        ELSE
                           J2=(MOD(4,J+3))
                           XPARP=(((X(J2)-OFFSET)/Laenge))*Parame
                        ENDIF
c
c       J2 ist der Index des naechsten Wandpunktes
c
                          DCORVG(1,INP(J))=x(j)-
     *                   (x(J)-X(j2))*OME2
                          DCORVG(2,INP(J))=y(j)-
     *                   (Y(J)-PARY(xparp))*OME2
                  ENDIF

 24               CONTINUE
               ENDIF
            ENDIF
C
            DO 101 I=1,4
               KWORK(L(LCHCK2)-1+INP(I))=1
 101           CONTINUE
C

 10         CONTINUE
c
c     Korrektur der Parameterwerte fuer Randknoten und -mittelpunkte
c
c
            DO 30 I=1,NVBD
               TI=DWORK(L(LVBDP)-1+i)
               IF ((TI.GT.TECK1).AND.(TI.LT.T1)) then
                  T1=TI
                  I1=i
                  endif
               IF ((TI.LT.TECK2).AND.(TI.GT.T2)) then 
                  T2=TI
                  I2=i
                  endif
               IF ((TI.GT.TECK3).AND.(TI.LT.T3)) then
                  T3=TI
                  I3=i
                  endif
               IF ((TI.LT.TECK4).AND.(TI.GT.T4)) then 
                  T4=TI
                  I4=i
                  endif
 30            CONTINUE
c
               DWORK(L(LVBDP)-1+I1  )=T1-(T1-TECK1)*OME2
               DWORK(L(LMBDP)-1+I1-1)=(DWORK(L(LVBDP)-1+I1)+TECK1)/2d0
               DWORK(L(LMBDP)-1+I1  )=(DWORK(L(LVBDP)-1+I1)+
     *                                 DWORK(L(LVBDP)+I1))/2d0
c
               DWORK(L(LVBDP)-1+I2)  =T2+(TECK2-T2)*OME1
               DWORK(L(LMBDP)-1+I2)  =(TECK2+DWORK(L(LVBDP)-1+I2))/2d0
               DWORK(L(LMBDP)-1+I2-1)=(DWORK(L(LVBDP)-1+I2)+
     *                                 DWORK(L(LVBDP)-2+I2))/2d0
c
               DWORK(L(LVBDP)-1+I3)  =T3+(TECK3-T3)*OME1
               DWORK(L(LMBDP)-1+I3-1)=(DWORK(L(LVBDP)-1+I3)+TECK3)/2d0
               DWORK(L(LMBDP)-1+I3)  =(DWORK(L(LVBDP)-1+I3)+
     *                               DWORK(L(LVBDP)  +I3))/2d0

c
               DWORK(L(LVBDP)-1+I4  )=T4-(T4-TECK4)*OME2
               DWORK(L(LMBDP)-1+I4  )=(TECK4+DWORK(L(LVBDP)-1+I4))/2d0
               DWORK(L(LMBDP)-1+I4-1)=(DWORK(L(LVBDP)-1+I4)+
     *                                 DWORK(L(LVBDP)-1+I4-1))/2d0
c

         CALL ZDISP (0,LCHECK,'LCHECK')
         CALL ZDISP (0,LCHCK2,'LCHECK')

c
c
      END

**********************************************************
      SUBROUTINE RATEST (DCORVG,KVERT,KMID,KVBD,KEBD,KRAND)
***********************************************************
C
C     Bestimmt den Abstand der einzelnen Elemente von der
C     naechsten WAND (DIRICHLET RAND=0) und speichert diese
C     auf den Vektor LDIST  (DIMENSION NEL).
C	Als Abstand wir der Abstand der Elementmitte
C	Zu den Mitten der Kanten der Randelemente berechnet
C     Zusaetzlich wird in LWAND die Nummer des naechsten
C	Wandelementes gespeichert. 
C	NUR EINE RANDKOMPONENTE
C
c      PARAMETER (NNARR=299,NNLEV=9,NNAB=21,NNWORK=1)
C
      INCLUDE 'common.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
C
C *** Standard dimensioning for workspace concept
      DIMENSION DCORVG(2,*),KVERT(4,*),KVBD(*),KMID(4,*),KEBD(*)
      DIMENSION KRAND(*)
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** Standard COMMON blocks
c      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
C***********************************************************
      SAVE
C
C	Welche ECKEN liegen auf dem Neumannrand?
C
      DO 10 I=1,NVBD
         VERTNR=KVBD(I)
         ELNR  =KEBD(I)
         DO 11 ILOC=1,4
            VTLOC=ILOC
            IF (KVERT(ILOC,ELNR).EQ.VERTNR) GOTO 12
11       CONTINUE
         WRITE (*,*) 'SCHEISSE, WAS IST LOS? FEHLER IN RATEST'
12       MIDNR=KMID(VTLOC,ELNR)
         IF (KWORK(L(LNPR)-1+MIDNR).EQ.0) then
                KWORK(L(LNPR)-1+VERTNR)=-1        
         endif
10    CONTINUE 
C
C***********************************************************
C
C
      DO 1 IEL=1,NEL
         neck1=KVERT(1,IEL)
         neck2=KVERT(2,IEL)
         neck3=KVERT(3,IEL)
         neck4=KVERT(4,IEL)
c
c	KOORDINATEN DER ECKEN
c
         x1=DCORVG(1,neck1)
         y1=DCORVG(2,neck1)
         x2=DCORVG(1,neck2)
         y2=DCORVG(2,neck2)
         x3=DCORVG(1,neck3)
         y3=DCORVG(2,neck3)
         x4=DCORVG(1,neck4)
         y4=DCORVG(2,neck4)
c
         XX=(x1+x2+x3+x4)/4D0
         YY=(y1+y2+y3+y4)/4d0
         DELTA=1D99
C

         DO 2 IVBD=1,NVBD
            NRVERT=KVBD(IVBD)
c
            IF (KWORK(L(LNPR)-1+NRVERT).LE.0) GOTO 2  
c 		Dann ist die folgende Kante Neumann
            IF (IVBD.EQ.NVBD) THEN
              IV1=1
            ELSE
              IV1=IVBD+1
            ENDIF

            NRVER2=KVBD(IV1)

            XVERT1=(DCORVG(1,NRVERT))
            YVERT1=(DCORVG(2,NRVERT))
            XVERT2=(DCORVG(1,NRVER2))
            YVERT2=(DCORVG(2,NRVER2))
c
           XMID=(XVERT1+XVERT2)/2D0
           YMID=(YVERT1+YVERT2)/2D0
C
            VRAND=FDATIN(1,1,XMID,YMID,0d0,RE)
            IF (VRAND.NE.0d0) GOTO 2
C              Dann ist der Knoten keine Wand,
C              Sondern Einströmrand
CC            ABST1=SQRT((xvert1-xx)**2+(yvert1-yy)**2)
CC            ABST2=SQRT((xvert2-xx)**2+(yvert2-yy)**2)
C
            ABST=SQRT((XMID-xx)**2+(YMID-yy)**2)
C
CC            IF ((ABST1+ABST2)/2d0.lt.DELTA) THEN
            IF (ABST.lt.DELTA) THEN
                DELTA =ABST
                IELRND=KEBD(IVBD)
            ENDIF
2        CONTINUE
C
	 DWORK(L(LDIST)-1+IEL)=DELTA
         KRAND(IEL) =IELRND
1     CONTINUE
C
      DO 4 I=1,NVT
4     IF (KWORK(L(LNPR)-1+I).EQ.-1) KWORK(L(LNPR)-1+I)=1
C
C
      END
