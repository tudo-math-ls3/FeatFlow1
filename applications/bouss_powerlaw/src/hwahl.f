************************************************************************
      SUBROUTINE HWAHL (HLOCAL, UNORM,  XBETA1, 
     *                  XBETA2, IEL,KVERT,KMID,DCORVG)
************************************************************************
C
C	WAEHLT EIN LOKALES H
C
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL LAMBDA
C
      PARAMETER (NNBAS=21,NNLEV=9,NNARR=299,NNWORK=1,NNVE=4)
C
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
C
      SAVE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	NUMMERN DER ECKEN
c
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
c     Skalieren
c
c      skal=max(xbeta1,xbeta2)
      xbeta1=xbeta1
      xbeta2=xbeta2
      alpmax=0d0
c    
c
c     es werden alle 4 Ecken durchgegangen und getestet,
c     ob die Gerade mit Steigung Beta durch diese Ecke
c     eine der Viereckskanten schneidet.
c     MERKE: es muessen pro Ecke nur zwei Kanten getestet werden,
c     naemlich die gegenueberliegenden.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          ERSTE ECKE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call schnitt(x1,y1,alpha,xbeta1,xbeta2,
     *            x3,y3,lambda,x2,y2)
       ALPMAX=MAX(alpha,alpmax)
c
      call schnitt(x1,y1,alpha,xbeta1,xbeta2,
     *            x3,y3,lambda,x4,y4)
       ALPMAX=MAX(alpha,alpmax)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          ZWEITE ECKE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call schnitt(x2,y2,alpha,xbeta1,xbeta2,
     *            x4,y4,lambda,x1,y1)
       ALPMAX=MAX(alpha,alpmax)
c
      call schnitt(x2,y2,alpha,xbeta1,xbeta2,
     *            x4,y4,lambda,x3,y3)
       ALPMAX=MAX(alpha,alpmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          DRITTE ECKE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call schnitt(x3,y3,alpha,xbeta1,xbeta2,
     *            x1,y1,lambda,x2,y2)
       ALPMAX=MAX(alpha,alpmax)
c
      call schnitt(x3,y3,alpha,xbeta1,xbeta2,
     *            x1,y1,lambda,x4,y4)
       ALPMAX=MAX(alpha,alpmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          VIERTE ECKE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call schnitt(x4,y4,alpha,xbeta1,xbeta2,
     *            x2,y2,lambda,x1,y1)
       ALPMAX=MAX(alpha,alpmax)
c
      call schnitt(x4,y4,alpha,xbeta1,xbeta2,
     *            x2,y2,lambda,x3,y3)
       ALPMAX=MAX(alpha,alpmax)
c
      HLOCAL=ALPMAX*4D0*unorm
c
C
      END

*******************************************************************
      SUBROUTINE SCHNITT (xo,yo,alpha,beta1,beta2
     *                   ,xa,ya,lambda,xb,yb)
C
C   Schnitt von zwei Geraden im R^2
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)  
      DOUBLE PRECISION beta1, beta2, alpha
      REAL lambda, skal
c
      skal=beta2*(xb-xa)-beta1*(yb-ya)
C      
      if (skal.eq.0D0) then
C  
C        beta und der Richtungsvektor sind parallel
C
         alpha=0D0
c         write(*,*) 'eins'
         bflag=.false.
      else  
           lambda=real((beta1*(ya-yo)-beta2*(xa-xo))/skal)
           bflag=.true.      
      endif
C
C     Ist der Schnittpunkt innerhalb des Elements?
C
      if (bflag) then
      if ((lambda.ge.-1e-1).and.(lambda.le.1.11e0)) then
         if (beta1.ne.0D0) then
            alpha=((xa-xo)+lambda*(xb-xa))/beta1
         else
            if (beta2.ne.0D0) then
               alpha=((ya-yo)+lambda*(yb-ya))/beta2
            else
               alpha=0D0
           endif
         endif
      else
c         write(*,*) 'drei'
         alpha=0D0
      endif
      endif
      END
