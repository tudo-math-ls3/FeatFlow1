***********************************************************************
* Calculates the drag/lift values.
*
* In:
*  ILV    - The level where tha calculation takes place. 
*           Normally NLMAX
*  ICMETH - Calculation method.
*           0=no force calculation, only calculate approx. volume
*           1=Body-Force
*           2=Volume integration, constant pressure (standard)
*           3=Volume integration, linear pressure
*           4=Volume integration, variable alpha
*  BKNPUP - true =update KNPR-arrays before computation
*           false=don't update KNPR-arrays
*
* Out:
*  DFWX - Drag coefficient
*  DFWY - Lift coefficient
*  DVOL - approximative volume of fictitious boundary components
*
* The level will be changed if necessary.
*
* If BKNPUP=true, the KNPR-arrays will be restored to the backup KNPRO 
* and rebuit according to the current setting of the fictitious boundary
* implementation. This can be necessary for a proper calculation of
* the drag/lift values if the radii of the obstacle(s) in the geometry
* have been relaxated by a modification of DCRRLX/DCRRLY.
* The KNPR-array is not restored to it's state before the call,
* so e.g. if the caller resets DCRRLX/DCRRLY to 0D0, it has to take
* care of resetting the KNPR-arrays too.
***********************************************************************

      SUBROUTINE DLCALC(ILV,ICMETH,BKNPUP,DFWX,DFWY,DVOL)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      INCLUDE 'cnspts.inc'
      
C *** user COMMON blocks
      INCLUDE 'cinidat.inc'
      
      INCLUDE 'cfiles.inc'
      
C parameters

      LOGICAL BKNPUP
      INTEGER ILV,ICMETH,I
      DOUBLE PRECISION DFWX,DFWY

C externals

      EXTERNAL EM30, EM31, E030, E031

C local variables

      INTEGER LAL,LALFA,LMALF,ISETLV,LTMP,KPL,KAUXM,LAREA
      DOUBLE PRECISION DVOL
      
      DFWX = 0D0
      DFWY = 0D0
      DVOL = 0D0

C Switch the level

      ILEV=ILV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      
      IF (BKNPUP) THEN

C Update KNPR-arrays on all levels to current configuration
C of the fictitious boundary, so that the integration routines
C (which rely on KNPR) work correctly.

        DO I=NLMIN,NLMAX
          CALL XLCP3(KLNPRO(I),KLNPR(I),KNVT(I)+KNMT(I))
          CALL IMPBDG (I)
        END DO
      END IF
      
C allocate some temporary space

      CALL ZNEW (2*NVT,1,LTMP,'KTMP  ')

C calculate pressure vector if we have to calculate drag/lift

      KPL  =L(LTMP)
      KAUXM=L(LTMP)+NVT
      LAREA=KLAREA(NLEV)

      IF (ICMETH.NE.0) THEN
        CALL  INTPV (DWORK(KP),DWORK(KPL),DWORK(KAUXM),DWORK(L(LAREA)),
     *               KWORK(L(LVERT)),KWORK(L(LNPR)))
      END IF

C=======================================================================
C *** lift and drag, featflow original
C=======================================================================
      
      IF (ICMETH.EQ.1) THEN
        IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *    CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EM31,DFWX,DFWY)

        IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *    CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EM30,DFWX,DFWY)
      END IF

C=======================================================================
C *** lift and drag by volume integration, constant pressure
C=======================================================================

      IF (ICMETH.EQ.2) THEN
        CALL ZNEW (NU,1,LAL,'AL   ')

        IF (IELT.EQ.0) 
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),E031,
     *    DWORK(L(lal)),DFWX,DFWY,0)
        IF (IELT.EQ.1) 
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),E030,
     *    DWORK(L(lal)),DFWX,DFWY,0)
        IF (IELT.EQ.2) 
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),EM31,
     *    DWORK(L(lal)),DFWX,DFWY,0)
        IF (IELT.EQ.3)
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),EM30,
     *    DWORK(L(lal)),DFWX,DFWY,0)

        CALL ZDISP (0,LAL,'al   ')
        
      END IF

C=======================================================================
C *** lift and drag by volume integration, linear pressure
C=======================================================================

      IF (ICMETH.EQ.3) THEN
        CALL ZNEW (NU,1,LAL,'al   ')

        IF (IELT.EQ.0) 
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KPL),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),E031,
     *    DWORK(L(lal)),DFWX,DFWY,1)
        IF (IELT.EQ.1) 
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KPL),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),E030,
     *    DWORK(L(lal)),DFWX,DFWY,1)
        IF (IELT.EQ.2) 
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KPL),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),EM31,
     *    DWORK(L(lal)),DFWX,DFWY,1)

        IF (IELT.EQ.3)
     *    CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KPL),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),EM30,
     *    DWORK(L(lal)),DFWX,DFWY,1)

        CALL ZDISP (0,LAL,'al   ')
      END IF

      
C=======================================================================
C *** lift and drag by volume, L2-projected alpha-vector
C=======================================================================

      IF (ICMETH.EQ.4) THEN
        CALL ZNEW(NMT,1,LALFA,'DALFA  ')
        CALL ZNEW(NMT,1,LMALF,'DMALF  ') 

        IF (IELT.EQ.0) 
     *    CALL BDFVL2 (DWORK(KU1),DWORK(KU2),DWORK(KP),NY,
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *                 DWORK(L(LCORVG)),DWORK(L(LCORMG)),  
     *                 E031,E031,IELT,IELT,8,DFWX,DFWY,
     *                 DWORK(L(LALFA)),DWORK(L(LMALF)))
        IF (IELT.EQ.1) 
     *    CALL BDFVL2 (DWORK(KU1),DWORK(KU2),DWORK(KP),NY,
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *                 DWORK(L(LCORVG)),DWORK(L(LCORMG)),  
     *                 E030,E030,IELT,IELT,8,DFWX,DFWY,
     *                 DWORK(L(LALFA)),DWORK(L(LMALF)))     
        IF (IELT.EQ.2) 
     *    CALL BDFVL2 (DWORK(KU1),DWORK(KU2),DWORK(KP),NY,
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *                 DWORK(L(LCORVG)),DWORK(L(LCORMG)),  
     *                 EM31,EM31,IELT,IELT,8,DFWX,DFWY,
     *                 DWORK(L(LALFA)),DWORK(L(LMALF)))

        IF (IELT.EQ.3)
     *    CALL BDFVL2 (DWORK(KU1),DWORK(KU2),DWORK(KP),NY,
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *                 DWORK(L(LCORVG)),DWORK(L(LCORMG)),  
     *                 EM30,EM30,IELT,IELT,8,DFWX,DFWY,
     *                 DWORK(L(LALFA)),DWORK(L(LMALF)))

        CALL ZDISP(0,LALFA,'DALFA  ')
        CALL ZDISP(0,LMALF,'DMALF  ') 
      END IF
      
      CALL ZDISP (0,LTMP,'KTMP  ')

      END

************************************************************************
      SUBROUTINE BDFORC(DU1,DU2,DP,KVERT,KMID,KVBD,KEBD,KMM,DCORVG,ELE,
     *                  DFW,DAW)
************************************************************************
*    Purpose:  Calculates lift (DFW) and drag (DAW)
*-----------------------------------------------------------------------
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cnspts.inc'

C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DFW,DAW
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*),KMM(2,*)
      EXTERNAL ELE

C externals

      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER IELTYP,IW1,IW2,IVBD,IVT,IVT1,IVT2,ISTOP,II,IVBD1
      INTEGER IVE,JP,JDFL
      DOUBLE PRECISION DLEN,PX1,PX2,PY1,PY2,PXM,PYM,DLH
      DOUBLE PRECISION DTX,DTY,DNX,DNY,DPCONT,DJ1,DJ2,DJ3,DJ4,XX,YY,DUT

      DFW=0D0
      DAW=0D0
      DLEN=0D0
C
C
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
C
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      NCUBP=1
      ICUBP=1
C
      IW1=1
      IW2=1
C
      DO 10 IVBD=1,NVBD
      IVT=KVBD(IVBD)
      IF (IVT.EQ.KMM(1,2)) GOTO 15
10    CONTINUE
C
      WRITE(6,*) 'ERROR 1 IN BDFORC'
      RETURN
C
15    IVBD1=IVBD
      ISTOP=0
C
      DO 100 IVBD=IVBD1,NVBD
      IF (ISTOP.EQ.1) GOTO 1000
      IEL =KEBD(IVBD)
      IVT1=KVBD(IVBD)
      IF (IVT1.EQ.KMM(2,2)) THEN
       IVT2=KVBD(IVBD1)
       ISTOP=1
      ELSE
       IVT2=KVBD(IVBD+1)
      ENDIF
C
      DO 101 II=1,4
101   IF (KVERT(II,IEL).EQ.IVT1) GOTO 102
      WRITE(6,*) 'WRONG 1'
      IW1=IW1+1
102   DO 103 II=1,4
103   IF (KVERT(II,IEL).EQ.IVT2) GOTO 104
      WRITE(6,*) 'WRONG 2'
      IW2=IW2+1
104   CONTINUE
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PXM=0.5D0*(PX1+PX2)
      PYM=0.5D0*(PY1+PY2)
      DLH=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
C
      DTX= (PX2-PX1)/DLH
      DTY= (PY2-PY1)/DLH
      DNX=-(PY2-PY1)/DLH
      DNY= (PX2-PX1)/DLH
      DPCONT=DP(IEL)
      DLEN=DLEN+DLH
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
      XX=PXM
      YY=PYM
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(XX,YY,-3)
C
      DUT=0
      DO 130 JDFL=1,IDFL
      DUT=DUT+DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTX*DNX
     *       +DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTY*DNX
     *       +DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTX*DNY
     *       +DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTY*DNY
130   CONTINUE
C
      DFW=DFW+DLH*(DPF(1)*DUT*DNY-DPCONT*DNX)
      DAW=DAW-DLH*(DPF(1)*DUT*DNX+DPCONT*DNY)
C
100   CONTINUE
C
1000  DFW=2D0*DFW/DPF(2)      
      DAW=2D0*DAW/DPF(2) 
C
C
C
99999 END
C

c#######################################################################
c#   new lift and drag evaluations .... 
c#######################################################################

c
************************************************************************
      SUBROUTINE BDFVOL(DU1,DU2,DP,KVERT,KNPR,KMID,
     *    DCORVG,ELE,DAL,DFWX,DFWY,IP)
************************************************************************
*       purpose: calculates lift (dfwy) and drag (dfwx) 
*       volume integration
*       ip=0,1 ... for constant, linear pressure
*-----------------------------------------------------------------------
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cmgadr.inc'

      INCLUDE 'cinidat.inc'
      INCLUDE 'cnspts.inc'

C parameters
      DOUBLE PRECISION du1(*),du2(*),dp(*),dcorvg(2,*),dal(*)
      INTEGER knpr(*),kvert(nnve,*),kmid(nnve,*)
      integer ip
      double precision dfwx,dfwy

C externals
      EXTERNAL ELE
      
C local variables

      INTEGER iel0,ibd,ive,im1,ivt
      save
c       
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
      
      
C number of FIXED boundary part where to compute the forces
      !ibd=2  

C set to negative to compute forces on the ficitious boundary only
      ibd=-1

c       prepare vector dal such that it is 1 at the boundary where the
c       force should be computed and 0 otherwise 
      
      call  lcl1( dal, nu) 
      
      do iel0=1,nel
        do ive=1,4
          
          im1=kmid(ive,iel0)
          if((knpr(im1).lt.0).and.(ibd.lt.0)) then
            dal(im1-nvt) = 1.0d0
          else
            if((knpr(im1).gt.0).and.(ibd.lt.0)) then
              ivt=knpr(im1)
              if(knpr(ivt).eq.ibd) then
                dal(im1-nvt) = 1.0d0
              endif
            endif
          endif
c       
        enddo
      ENDDO
      
c       use the vector dal as a test function in the weak formulation
c       (only the difusion part is taken, should take even the
c       convection but in featflow it seems to be ok) ... 
c       this needs to be checked in timedependent case
      
      dfwx=0d0
      dfwy=0d0

c       ip=0 ... constant pressure , ip=1 ... linear pressure
      if(ip.eq.0) then
        call cstgforce(du1,du2,dp,kvert,kmid,
     *      dcorvg,ele,dfwx,dfwy,dal,ielt,8,dpf(1),dpf(2))
      else
        call lstgforce(du1,du2,dp,kvert,kmid,
     *      dcorvg,ele,dfwx,dfwy,dal,ielt,8,dpf(1),dpf(2))
      endif        
      if (dpf(2).ne.0d0) then
        dfwx=2d0*dfwx/dpf(2)
        dfwy=2d0*dfwy/dpf(2)
      endif      
      
99999 end

c=============================================================================
c  the volume integration with constant pressure dp(1..nel)

      SUBROUTINE CSTGFORCE(DU1,DU2,DP,KVERT,KMID,
     *    DCORVG,ELE,DFWX,DFWY,DAL,IELT,ICUB,DPF1,DPF2)
c       
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'

      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cnspts.inc'

C parameters

      DOUBLE PRECISION du1(*),du2(*),dp(*),dcorvg(2,*),dal(*)
      double precision dfwx,dfwy,dpf1,dpf2
      INTEGER kvert(nnve,*),kmid(nnve,*)
      integer ielt, icub

C externals

      external ele
      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER IELTYP,I,IVE,JP,IG
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DAV,DAX,DAY,DPV
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,OM
      DOUBLE PRECISION XX,YY,XI1,XI2,DN1,DN2,AH1,AH2

c       *** preparation - evaluation of parameters
      ier=0

c       *** which derivatives of basis functions are needed?
      do  i = 1,nnder
        bder(i)=.false.
      enddo
      
      bder(1)=.true.
      bder(2)=.true.
      bder(3)=.true.
      
c       *** dummy call of ele sets number of element
      ieltyp=-1
      call ele(0d0,0d0,ieltyp)

      idfl=ndfl(ieltyp)
      call cb2q(icub)
      if (ier.ne.0) goto 99999
c       
c       *** dummy call - ele may save arithmetic operations
      icubp=icub
      call ele(0d0,0d0,-2)
      
      dfwx=0d0
      dfwy=0d0
c       
c       *** loop over all elements
      do iel=1,nel
        call ndfgl(iel,1,ieltyp,kvert,kmid,kdfg,kdfl)
        if (ier.lt.0) goto 99999
c       
c       *** evaluation of coordinates of the vertices
        do ive = 1, nve
          jp=kvert(ive,iel)
          kve(ive)=jp
          dx(ive)=dcorvg(1,jp)
          dy(ive)=dcorvg(2,jp)
        enddo
c       
        dj1=0.5d0*(-dx(1)-dx(2)+dx(3)+dx(4))
        dj2=0.5d0*( dx(1)-dx(2)+dx(3)-dx(4))
        dj3=0.5d0*(-dy(1)+dy(2)-dy(3)+dy(4))
        dj4=0.5d0*(-dy(1)+dy(2)+dy(3)-dy(4))
c       
        call ele(0d0,0d0,-2)
c       *** loop over all cubature points
        do icubp = 1, ncubp
          
c       *** cubature point on the reference element
          xi1=dxi(icubp,1)
          xi2=dxi(icubp,2)
          
c       *** jacobian of the bilinear mapping onto the reference element
          djac(1,1)=0.5d0*(dx(2)-dx(1)+dj2)+0.5d0*dj2*xi2
          djac(1,2)=0.5d0*dj1+0.5d0*dj2*xi1
          djac(2,1)=0.5d0*dj4-0.5d0*dj3*xi2
          djac(2,2)=0.5d0*(dy(3)-dy(1)-dj4)-0.5d0*dj3*xi1
          detj=djac(1,1)*djac(2,2)-djac(1,2)*djac(2,1)
          om=domega(icubp)*detj
          
c       *** cubature point on the real element
          xx=0.5d0*(dx(1)+dx(2)+dj1)+0.5d0*(dx(2)-dx(1)+dj2)*xi1
     *        +0.5d0*dj1*xi2+0.5d0*dj2*xi1*xi2
          yy=0.5d0*(dy(1)+dy(3)+dj3)+0.5d0*dj4*xi1
     *        +0.5d0*(dy(3)-dy(1)-dj4)*xi2-0.5d0*dj3*xi1*xi2
          
c       *** evaluate the basis functions in the cubature point
c       for the velocities
          if((ielt.eq.2).or.(ielt.eq.3)) then
            call ele(xx,yy,-3)
          else
            call ele(xi1,xi2,-3)
          endif
          if (ier.lt.0) goto 99999

c       evaluate the solution values and derivatives in the cubature point     
          du1v=0d0 ! value
          du1x=0d0 ! x dreiv.
          du1y=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du1v=du1v+du1(ig)*dbas(kdfl(i),1)
            du1x=du1x+du1(ig)*dbas(kdfl(i),2)
            du1y=du1y+du1(ig)*dbas(kdfl(i),3)
          enddo
          
          du2v=0d0 ! value
          du2x=0d0 ! x dreiv.
          du2y=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du2v=du2v+du2(ig)*dbas(kdfl(i),1)
            du2x=du2x+du2(ig)*dbas(kdfl(i),2)
            du2y=du2y+du2(ig)*dbas(kdfl(i),3)
          enddo
          
          dav=0d0 ! value
          dax=0d0 ! x dreiv.
          day=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            dav=dav+dal(ig)*dbas(kdfl(i),1)
            dax=dax+dal(ig)*dbas(kdfl(i),2)
            day=day+dal(ig)*dbas(kdfl(i),3)
          enddo

c       form the integrand

          dn1=-dax
          dn2=-day

          dpv=dp(iel)
c
c          ah1=du1v*du1x*dn1+du2v*du1y*dn1
c     *        -dpv*dn1+dpf1*(2.0d0*du1x*dn1+(du2x+du1y)*dn2)
c          ah2=du1v*du2x*dn1+du2v*du2y*dn1
c     *        -dpv*dn2+dpf1*((du2x+du1y)*dn1+2.0d0*du2y*dn2)
c
          ah1=-dpv*dn1+dpf1*(du1x*dn1+du1y*dn2)
          ah2=-dpv*dn2+dpf1*(du2x*dn1+du2y*dn2)

          dfwx=dfwx+ah1*om
          dfwy=dfwy+ah2*om
          
        enddo
      enddo
c       
c       
99999 END

C=============================================================================
c  the volume integration with linear pressure dp(1..nel)

      SUBROUTINE LSTGFORCE(DU1,DU2,DP,KVERT,KMID,
     *    DCORVG,ELE,DFWX,DFWY,DAL,IELT,ICUB,DPF1,DPF2)
c       
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cnspts.inc'

C parameters

      DOUBLE PRECISION du1(*),du2(*),dp(*),dcorvg(2,*),dal(*)
      double precision dfwx,dfwy,dpf1,dpf2
      INTEGER kvert(nnve,*),kmid(nnve,*)
      integer ielt

C externals

      EXTERNAL ELE,E011
      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER ITYP,I,ICUB,IVE,JP,IG
      INTEGER ITYP1, IDFL1, KDFG1(NNBAS), KDFL1(NNBAS)   
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DAV,DAX,DAY,DPV,OM
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XX,YY,DN1,DN2,AH1,AH2,XI1,XI2

c       *** preparation - evaluation of parameters
      ier=0

c       *** which derivatives of basis functions are needed?
      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
c       *** dummy call of ele sets number of element
      ityp=-1
      call ele(0d0,0d0,ityp)
      idfl=ndfl(ityp)

      ityp1=-1
      call e011(0d0,0d0,ityp1)
      idfl1=ndfl(ityp1)

      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
c       
c       *** dummy call - ele may save arithmetic operations
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
      
      dfwx=0d0
      dfwy=0d0
c       
c       *** loop over all elements
      DO IEL=1,NEL

        call ndfgl(iel,1,ityp,kvert,kmid,kdfg,kdfl)
        if (ier.lt.0) goto 99999
c       
        call ndfgl(iel,1,ityp1,kvert,kmid,kdfg1,kdfl1)
        IF (IER.LT.0) GOTO 99999
c       
c       *** evaluation of coordinates of the vertices
        DO IVE = 1, NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        ENDDO
C       
        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C       
        CALL ELE(0D0,0D0,-2)
c       *** loop over all cubature points
        DO ICUBP = 1, NCUBP
          
c       *** cubature point on the reference element
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
c       *** jacobian of the bilinear mapping onto the reference element
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ
          
c       *** cubature point on the real element
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *        +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *        +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
          
c       *** evaluate the basis functions in the cubature point
c       for the velocities
          if((ielt.eq.2).or.(ielt.eq.3)) then
            CALL ELE(XX,YY,-3)
          ELSE
            CALL ELE(XI1,XI2,-3)
          ENDIF
          IF (IER.LT.0) GOTO 99999
          
c       evaluate the solution values and derivatives in the cubature point     
          DU1V=0D0 ! VALUE
          DU1X=0D0 ! X DREIV.
          DU1Y=0D0 ! Y DERIV
          DO I=1,IDFL
            IG=KDFG(I)
            du1v=du1v+du1(ig)*dbas(kdfl(i),1)
            du1x=du1x+du1(ig)*dbas(kdfl(i),2)
            du1y=du1y+du1(ig)*dbas(kdfl(i),3)
          ENDDO
          
          DU2V=0D0 ! VALUE
          DU2X=0D0 ! X DREIV.
          DU2Y=0D0 ! Y DERIV
          DO I=1,IDFL
            IG=KDFG(I)
            du2v=du2v+du2(ig)*dbas(kdfl(i),1)
            du2x=du2x+du2(ig)*dbas(kdfl(i),2)
            du2y=du2y+du2(ig)*dbas(kdfl(i),3)
          ENDDO
          
          DAV=0D0 ! VALUE
          DAX=0D0 ! X DREIV.
          DAY=0D0 ! Y DERIV
          DO I=1,IDFL
            IG=KDFG(I)
            dav=dav+dal(ig)*dbas(kdfl(i),1)
            dax=dax+dal(ig)*dbas(kdfl(i),2)
            day=day+dal(ig)*dbas(kdfl(i),3)
          enddo

c       *** evaluate the basis functions in the cubature point
c       for the pressure
          call e011(xi1,xi2,-3)
          if (ier.lt.0) goto 99999
          
c       evaluate the solution values and derivatives in the cubature point     
          dpv=0d0 ! value
          do i=1,idfl1
            ig=kdfg1(i)
            dpv=dpv+dp(ig)*dbas(kdfl1(i),1)
          ENDDO

c       form the integrand

          DN1=-DAX
          DN2=-DAY
c
c          ah1=du1v*du1x*dn1+du2v*du1y*dn1
c     *        -dpv*dn1+dpf1*(2.0d0*du1x*dn1+(du2x+du1y)*dn2)
c          ah2=du1v*du2x*dn1+du2v*du2y*dn1
c     *        -dpv*dn2+dpf1*((du2x+du1y)*dn1+2.0d0*du2y*dn2)

          ah1=-dpv*dn1+dpf1*(du1x*dn1+du1y*dn2)
          ah2=-dpv*dn2+dpf1*(du2x*dn1+du2y*dn2)

          dfwx=dfwx+ah1*om
          dfwy=dfwy+ah2*om
          
        ENDDO
      ENDDO
c       
c       
99999 END

***********************************************************************
* Drag and Lift calculation
* Non-constant ALPHA
*
* In:
*  DFALF  - Workspace vector of length NMT for ALPHA-vectors
*  DM     - Workspace vector of length NMT for transformation with
*           mass matrix
*  IELT   - Type of element (0-3) for ALPHA discretization
*  ELE    - Element function for ALPHA discretization
*  IELT1  - Type of element (0-3) for evaluation of Drag/Lift
*  ELE1   - Element function for evaluating Drag/Lift
*
* Out:
*  DFXVIF - Drag
*  DFYVIF - Lift
*
* The method used here is a projection of the piecewise constant
* funtion ALPHA into the finite element space determined by ELE.
* The L2-projection used here is 
*
*    (FE-projected alpha,phi_i) = (alpha,phi_i)
*
* Translating this into matrix-vector  means solving the discrete
* system
*                           Ma~ = b
* 
* with b the RHS-vector that results by building the right hand side
* using the piecewise constant function ALPHA. M is the mass matrix
* and a~ the FE-vector of the projected ALPHA.
* The mass matrix is lumped here for simplicity.
***********************************************************************

      SUBROUTINE BDFVL2 (DU1,DU2,DP,DNY,KVERT,KMID,KNPR,DCORVG,
     *                   DCORMG,ELE,ELE1,IELT,IELT1,ICUB,
     *                   DFXVIF,DFYVIF,DALFAF,DM)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'

      INCLUDE 'cnspts.inc'

C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DALFAF(*),DM(*)
      DOUBLE PRECISION DNY, DCORVG(2,*),DCORMG(2,*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
      INTEGER IELT, IELT1, ICUB
      DOUBLE PRECISION DFXVIF, DFYVIF
     
C externals

      EXTERNAL ELE,ELE1
      INTEGER NDFL
      EXTERNAL NDFL
C      INTEGER ISFBDY
C      EXTERNAL ISFBDY
      
C local variables

      INTEGER IELTYP,I,IVE,JP,IG,JDOFEL
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,OM
      DOUBLE PRECISION XX,YY,XI1,XI2,DN1,DN2,AH1,AH2
      DOUBLE PRECISION DALV,DALX,DALY,DAL1X,DAL1Y
      DOUBLE PRECISION DBI1,DBI2,DBI3,DPP

C initialise the result variables

      DFXVIF=0D0
      DFYVIF=0D0

C initialise the ALPHA vector and lumped mass matrix vector

      DO I=1,NMT
        DALFAF(I)=0D0
        DM(I)=0D0                                             
      END DO

C Build the right hand side vector of the system.

      IF(IELT1.EQ.0) THEN
        CALL BALRHS(DALFAF,DM,KVERT,KMID,DCORVG,
     *                ELE1,3,IELT1)
      ELSE
        CALL BALRHS(DALFAF,DM,KVERT,KMID,DCORVG,
     *               ELE1,10,IELT1)
       END IF

C Solve the system by multiplying by M^-1.

      DO I=1,NMT                    
        DALFAF(I)=DALFAF(I)/DM(I)                       

CCC To test if the rest of the routine works correctly, comment out the
CCC previous line and comment in the next one. The calculated numbers
CCC should be the same like in CGSTFRCE.
CCC     IF (ISFBDY (DCORMG(1,I),DCORMG(2,I),0).GT.0) DALFAF(I)=1.0D0
      END DO
         
C Ok, the ALPHA-vector is prepared. Now we can continue calculating
C the drag-/lift-forces like in the routines above...

C----------------------------------------------------------------

c       *** Preparation - evaluation of parameters

      IER=0
      
c       *** Which derivatives of basis functions are needed?

      DO I = 1,NNDER
       BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.            
      
c       *** Dummy call of ele sets number of element

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
       
c       *** Dummy call - ele may save arithmetic operations

      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)

C============================================================================
c       *** Loop over all elements

      JDOFEL=1
       
C-----------------------------------------------------------------------
      
      DO IEL=1,NEL        !! LOOP FOR EVERY ELEMENT
       
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

c       *** evaluation of coordinates of the vertices

        DO IVE = 1, NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999        

c       *** Loop over all cubature points        

        DO ICUBP = 1, NCUBP    !!  LOOP FOR EVERY CUBTURE POINT
          
c       *** Cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
c       *** Jacobian of the bilinear mapping onto the reference element

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ
          
c       *** Cubature point on the real element
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *      +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2

C-----------------------------------------------------------------

c       *** Evaluate the basis functions in the cubature point
          IF((IELT.EQ.2).OR.(IELT.EQ.3)) THEN
            CALL ELE(XX,YY,-3)
          ELSE
            CALL ELE(XI1,XI2,-3)           
          END IF
          IF (IER.LT.0) GOTO 99999
          
c       Evaluate the solution values and derivatives in the cubature point     

          DU1V=0D0     ! U1 value
          DU1X=0D0     ! U1 x dreiv.
          DU1Y=0D0     ! U1 y deriv
          DU2V=0D0     ! U2 value
          DU2X=0D0     ! U2 x dreiv.
          DU2Y=0D0     ! U2 y deriv
          DALV=0D0     ! ALFA value
          DALX=0D0     ! ALFA x dreiv.
          DALY=0D0     ! ALFA y deriv
          DAL1X=0D0     ! ALFA x dreiv.
          DAL1Y=0D0     ! ALFA y deriv

          DO I=1,IDFL
            IG=KDFG(I)
            DBI1=DBAS(KDFL(I),1)
            DBI2=DBAS(KDFL(I),2)            
            DBI3=DBAS(KDFL(I),3)            
C---------------FOR U1----------------
            DU1V=DU1V+DU1(IG)*DBI1
            DU1X=DU1X+DU1(IG)*DBI2            
            DU1Y=DU1Y+DU1(IG)*DBI3
C---------------FOR U2----------------
            DU2V=DU2V+DU2(IG)*DBI1
            DU2X=DU2X+DU2(IG)*DBI2            
            DU2Y=DU2Y+DU2(IG)*DBI3                        
C---------------FOR ALFA----------------            
            DALV=DALV+DALFAF(IG)*DBI1
            DALX=DALX+DALFAF(IG)*DBI2            
            DALY=DALY+DALFAF(IG)*DBI3                            
          END DO
C--------------------------------------------------------

	    DPP=DP(IEL)
C	    DPP=0D0
          
c-----------Form the integrand------------------

          DN1=-DALX
          DN2=-DALY

c-------------------Acting force-------------------------          

          AH1=-DPP*DN1+DNY*(DU1X*DN1+DU1Y*DN2)           
          AH2=-DPP*DN2+DNY*(DU2X*DN1+DU2Y*DN2)

c
c        AH1=-DPP*DN1+DNY*(2D0*DU1X*DN1+(DU1Y+DU2X)*DN2)           
c        AH2=-DPP*DN2+DNY*((DU1Y+DU2X)*DN1+2D0*DU2Y*DN2)

c--------------------------------------------------------
c        AH1=-DP(IEL)*DN1+DNY*(DU1X*DN1+DU1Y*DN2)
c     *                  +DNY*(DU2X*DN2+DU2Y*DN1)            
c        AH2=-DP(IEL)*DN2+DNY*(DU2X*DN1+DU2Y*DN2)
c     *                  +DNY*(DU1X*DN2+DU1Y*DN1)         
                
          DFXVIF=DFXVIF+AH1*OM
          DFYVIF=DFYVIF+AH2*OM

        END DO     !! LOOP FOR CUBTURE POINT

c------------------------------------------------------------------

      END DO      !! LOOP FOR EVERY ELEMENT
       
      DFXVIF=2D0*DFXVIF/DPF(2)      
      DFYVIF=2D0*DFYVIF/DPF(2)

99999 END
      
C ---------------------------------------------------------------------
C Build right hand side for ALPHA-projection
C
C The above routine BDFVL2 performes a L2-projection of the piecewise
C constant function ALPHA into our finite element space Q1~.
C This routine builds the discrete right hand side vector for the
C system Ma=RHS that arises from the formula
C
C          (FE-projected alpha,phi_i) = (discrete alpha,phi_i)
C ---------------------------------------------------------------------
      
      SUBROUTINE BALRHS(DF,DM,KVERT,KMID,DCORVG,ELE,ICUB,IELT)

      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'

      INCLUDE 'casmbly.inc'
      
C parameters
      
      DOUBLE PRECISION DF(*),DM(*),DCORVG(2,*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),ICUB,IELT
      EXTERNAL ELE
      
C externals

      INTEGER NDFL,NFBDYC,ISFBDY
      EXTERNAL NDFL,NFBDYC,ISFBDY
      
C local variables

      INTEGER I, J, IELTYP, JDOFE1, IVE, JP, IG, INPR1, ISOLID
      DOUBLE PRECISION DJ1, DJ2, DJ3, DJ4, XI1, XI2, OM, XX, YY
      DOUBLE PRECISION DBI1, DBI2, DBI3, DBJ1, DBJ2, DBJ3

C *** Preparation - evaluation of parameters

      IER=0
          
C *** Which derivatives of basis functions are needed?

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
C     *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C     
C************************************************************************
C     *** Calculation of the matrix - storage technique 7 or 8
C************************************************************************
C     *** Dummy call - ELE may save arithmetic operations
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
     
C     *** Loop over all elements
      JDOFE1=1
      DO IEL=1,NEL
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C     *** Evaluation of coordinates of the vertices

        DO IVE = 1, NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO
     
        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
     
        CALL ELE(0D0,0D0,-2)

C     *** Loop over all cubature points

        DO ICUBP = 1, NCUBP

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C     *** Jacobian of the bilinear mapping onto the reference element

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM = DOMEGA(ICUBP)*DETJ

C     *** ELE needs the information ICUBP because of preceeding
            
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *           +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *           (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2

c       *** Evaluate the basis functions in the cubature point
      
          IF((IELT.EQ.2).OR.(IELT.EQ.3)) THEN
            CALL ELE(XX,YY,-3)
          ELSE
            CALL ELE(XI1,XI2,-3)     
          ENDIF

          IF (IER.LT.0) GOTO 99999                                    

          DO I=1,IDFL

            IG=KDFG(I)
            DBI1=DBAS(KDFL(I),1)
            DBI2=DBAS(KDFL(I),2)
            DBI3=DBAS(KDFL(I),3)
            INPR1=0             

            DO j=1,IDFL
              DBJ1=DBAS(KDFL(J),1)
              DBJ2=DBAS(KDFL(J),2)
              DBJ3=DBAS(KDFL(J),3)
              DM(IG)=DM(IG)+DBI1*DBJ1*OM
            END DO

            IF ((NFBDYC().GT.0).AND.(ISFBDY(XX,YY,0).GT.0)) THEN
              ISOLID=1
            ELSE
              ISOLID=0
            END IF
     
            DF(IG)=DF(IG)+DBLE(ISOLID)*DBI1*OM

          END DO

        END DO

      END DO

99999 END
      
