*CMZ :  4.00/15 27/04/2022  15.20.19  by  Michael Scheer
*CMZ :  3.06/00 14/01/2019  17.27.04  by  Michael Scheer
*CMZ :  3.03/02 18/01/2016  13.02.27  by  Michael Scheer
*CMZ :  3.02/03 30/10/2014  17.17.28  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random_gauss_omp(n,g,rr)

      implicit none

c Based on textbook "Numerical Recipies"
c The subroutine util_random is called, which can be initialized
c and controlled by util_random_init, util_random_set_seed, and
c util_random_get_seed.

      real g(n),r,v1,v2,fac,rr(2)
      integer n,i

c      print*,n

      if (n.eq.1) then
1       call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 1
        fac=sqrt(-2.*log(r)/r)
        g(1)=v1*fac
      else
        do i=1,(n/2*2),2
11      call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 11
        fac=sqrt(-2.*log(r)/r)
        g(i)=v1*fac
        g(i+1)=v2*fac
        enddo
      endif

      if (mod(n,2).ne.0) then
12      call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 12
        fac=sqrt(-2.*log(r)/r)
        g(n)=v1*fac
      endif

      return
      end
*CMZ :  3.03/02 29/02/2016  16.24.34  by  Michael Scheer
*CMZ :  3.02/03 05/09/2014  12.29.51  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random(n,r)
*KEEP,random.
      include 'random.cmn'
*KEND.

      real r(n)

      call random_number(r)

      irancalls=irancalls+n

      return
      end
*CMZ :  4.00/11 28/05/2021  09.17.01  by  Michael Scheer
*CMZ :  3.05/05 12/07/2018  13.12.16  by  Michael Scheer
*CMZ :  3.02/00 24/09/2014  13.51.08  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.24.14  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/19 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/08 15/12/2010  14.05.16  by  Michael Scheer
*CMZ : 00.00/07 07/05/2008  14.28.10  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_INTER(XA,YA,Y2A,N,X,Y,MODE)
*KEEP,GPLHINT.
*KEND.

C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       YA2:  ARRAY SPLINE COEFFICIENTS
C-       X: Y(X) IS CALCULATED
C-       MODE: CONTROL FLAG:
C-             MODE.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             MODE.LT.0: NEW INITIALIZATION

C--   OUTPUT:

C-       Y: Y(X) IS CALCULATED

      IMPLICIT NONE

      INTEGER NOLD,N,KLO,KHI,KLOLD,K,MODE,NORDER

      double precision Y,X,XA1OLD,XANOLD,H,A,B

      double precision XA(*),YA(*),Y2A(*),EPS,XX

c      save klold,nold,xa1old,xanold

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      if (klold.lt.1.or.klold.gt.n-1) klold=1
      if (nold.lt.1.or.nold.gt.n-1) nold=1
      if (xa1old.lt.xa(1).or.xa1old.gt.xa(n-1)) xa1old=xa(1)
      if (xanold.lt.xa(2).or.xanold.gt.xa(n)) xanold=xa(n)

      EPS=ABS(XA(N)-XA(1))/1.0D10
      XX=X

      IF(XA(1).LT.XA(N)) THEN

        IF(XX.LT.XA(1).AND.XX.GT.XA(1)-EPS) THEN
          XX=XA(1)
        ELSE IF(XX.GT.XA(N).AND.XX.LT.XA(N)+EPS) THEN
          XX=XA(N)
        ENDIF

        IF(XX.LT.XA(1).OR.XX.GT.XA(N)) THEN
          WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
          WRITE(6,*)'X:'
          WRITE(6 ,*)'***ERROR IN UTIL_SPLINE_INTER: X OUT OF RANGE ***'
          STOP
        ENDIF

      ELSE

        IF(XX.LT.XA(N).AND.XX.GT.XA(N)-EPS) THEN
          XX=XA(N)
        ELSE IF(XX.GT.XA(1).AND.XX.LT.XA(N)+EPS) THEN
          XX=XA(1)
        ENDIF

        IF(XX.LT.XA(N).OR.XX.GT.XA(1)) THEN
          WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
          WRITE(6,*)'X:',X
          WRITE(6 ,*)'***ERROR IN UTIL_SPLINE_INTER: X OUT OF RANGE ***'
          STOP
        ENDIF

      ENDIF

      norder=1
      if (xa(n).lt.xa(1)) then
        norder=-1
      endif

      if (norder.eq.1) then

        IF (MODE.LT.0.OR.KLOLD.GE.N) THEN
          KLO=1
        ELSE IF(NOLD.EQ.N
     &      .AND. XA(1).EQ.XA1OLD
     &      .AND. XA(N).EQ.XANOLD
     &      .AND. XX.GT.XA(KLOLD)
     &      ) THEN
          KLO=KLOLD
        ELSE
          KLO=1
        ENDIF

        IF (XX.LT.XA(KLO+1)) THEN
          KHI=KLO+1
          GOTO 2
        ENDIF

        KHI=N
1       IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.XX)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 1
        ENDIF

2       H=XA(KHI)-XA(KLO)

        IF (H.le.0.0D0) THEN
          WRITE(6 ,*)'*** ERROR IN UTIL_SPLINE_INTER: BAD INPUT ***'
          STOP
        ENDIF

        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     &    (A*(A+1.D0)*(A-1.D0)*Y2A(KLO)+B*(B+1.D0)*
     &    (B-1.D0)*Y2A(KHI))*(H**2)/6.D0

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)

      else !(norder.eq.1) then

        IF (MODE.LT.0.or.nold.ne.n) THEN
          KLO=1
        ELSE IF(
     &      XA(1).EQ.XA1OLD
     &      .AND. XA(N).EQ.XANOLD
     &      .AND. XX.lt.XA(KLOLD)
     &      ) THEN
          KLO=KLOLD
        ELSE
          KLO=1
        ENDIF

        IF (XX.gt.XA(KLO+1)) THEN
          KHI=KLO+1
          GOTO 21
        ENDIF

        KHI=N
11      IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).LT.XX)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 11
        ENDIF

21      H=XA(KHI)-XA(KLO)

        IF (H.ge.0.0D0) THEN
          WRITE(6 ,*)'*** ERROR IN UTIL_SPLINE_INTER: BAD INPUT ***'
          STOP
        ENDIF

        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     &    (A*(A+1.D0)*(A-1.D0)*Y2A(KLO)+B*(B+1.D0)*
     &    (B-1.D0)*Y2A(KHI))*(H**2)/6.D0

        KLOLD=KLO
        NOLD=N
        XA1OLD=XA(1)
        XANOLD=XA(N)

      endif !(norder.eq.1) then

      RETURN
      END
*CMZ :  4.00/11 28/05/2021  09.18.09  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.56.50  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.67/05 16/05/2012  12.45.37  by  Michael Scheer
*CMZ : 00.00/07 12/10/2009  12.17.45  by  Michael Scheer
*CMZ : 00.00/02 14/04/2003  12.46.09  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.48  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_COEF(X,Y,N,YP1,YPN,Y2,AA,BB,CC,C)
*KEEP,gplhint.
*KEND.

C--- CALCULATES SPLINE COEFFICIENTS

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       X: ARRAY OF X-VALUES
C-       Y: ARRAY OF Y-VALUES
C-       YP1:  SECOND DERIVATIVE AT FIRST X-VALUE
C-       YPN:  SECOND DERIVATIVE AT LAST X-VALUE

C--   OUPUT:

C-       Y2:   SPLINE-COEFFICIENTS

C--   WORKINGSPACE: AA(N),BB(N),CC(N),C(N)


      IMPLICIT NONE

      INTEGER N,J
      double precision  X(N),Y(N),Y2(N),AA(N),BB(N),CC(N),C(N)

      double precision YP1,YPN

      double precision xx(3),yy(3),a(3),yp(3),xopt,yopt
      INTEGER ifail

      IF (N.LT.3) then
        if (abs(yp1).eq.9999.0d0) then
          y2(1)=0.0d0
        else
          y2(1)=yp1
        endif
        if (abs(ypn).eq.9999.0d0) then
          y2(n)=0.0d0
        else
          y2(n)=ypn
        endif
        RETURN
      endif

      if (abs(yp1).eq.9999.0d0) then
        xx=x(1:3)
        yy=y(1:3)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(1)=2.0d0*a(3)
        else
          y2(1)=0.0d0
        endif
      else
        Y2(1)=YP1
      endif

      if (abs(ypn).eq.9999.0d0) then
        xx=x(n-2:n)
        yy=y(n-2:n)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(n)=2.0d0*a(3)
        else
          y2(N)=0.0d0
        endif
      else
        Y2(N)=YPN
      endif

      C(1)=Y2(1)
      C(N)=y2(n)

      BB(1)=1.D0
      CC(1)=0.D0
      CC(N)=1.D0

      DO J=2,N-1
        if(x(j+1).eq.x(j)) then
          write(6,*)
          write(6,*)
     &      '*** Error in util_spline_coef: Intervall of zero length'
          write(6,*)'j, x(j), x(j+1):',j,x(j),x(j+1)
          write(6,*)
          stop
        endif
          AA(J)=(X(J  )-X(J-1))/6.D0
          BB(J)=(X(J+1)-X(J-1))/3.D0
          CC(J)=(X(J+1)-X(J  ))/6.D0
          C(J)=(Y(J+1)-Y(J  ))/(X(J+1)-X(J  ))
     &          -(Y(J  )-Y(J-1))/(X(J  )-X(J-1))
      ENDDO !J

      DO J=2,N-1

          BB(J)=BB(J)-AA(J)*CC(J-1)
           C(J)= C(J)-AA(J)* C(J-1)
C          AA(J)=AA(J)-AA(J)*BB(J-1)

          CC(J)=CC(J)/BB(J)
           C(J)= C(J)/BB(J)
          BB(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
         Y2(J)=C(J)-CC(J)*Y2(J+1)
      ENDDO

      RETURN
      END
*CMZ :          31/10/2022  17.05.45  by  Michael Scheer
*CMZ : 00.00/16 19/03/2014  12.30.26  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.27.13  by  Michael Scheer
*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_PARABEL(Xin,Yin,A,YP,XOPT,yopt,IFAIL)

C--- CALCULATES A(1),A(2),A(3), THE DERIVATIVES YP(X(1)),YP(X(2)),YP(X(3)),
C    AND THE EXTREMUM (XOPT,A(XOPT)) OF PARABOLA Y=A1+A2*X+A3*X**2
C    FROM COORDINATES OF THE THREE POINTS (X(1),Y(1)),(X(2),Y(2)),(X(3),Y(3))
C

      IMPLICIT NONE

      INTEGER IFAIL

      REAL*8 A(3),X(3),Y(3),DXM,DXP,x0,a1,a2,dxm2,dxp2,dxmax,dymax,
     &  xin(3),yin(3),DET,YP(3),XOPT,yopt,a22,fm,fp,f0

      IFAIL=0

c calculate f=a0+a1*(x-x0)+a2*(x-x0)**2
c  = a0 + a1*x - a0*x0 + a2*x**2 - 2*a2*x*x0 + a2*x0**2
c  = a0 + (a2*x0 -a0)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

c change system: (x0,s0)->(0,0), i.e.
c calculate f=a1*dx+a2*dx**2
c  df/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

      x=xin
      y=yin

      call util_sort_func(3,x,y)

c      if (x(1).le.x(2).and.x(1).le.x(3)) then
c      endif

      x0=x(2)
      f0=y(2)

      fm=y(1)-f0
      fp=y(3)-f0

      dxm=x(1)-x0
      dxp=x(3)-x0

c fm=a1*dxm+a2*dxm**2
c fp=a1*dxp+a2*dxp**2

c (dxm dxm2) (a1) = (y(1))
c (dxp dxp2) (a2) = (y(3))

      dxm2=dxm*dxm
      dxp2=dxp*dxp

      det=dxm*dxp2-dxp*dxm2

      if (det.ne.0.0d0) then
        a1=(fm*dxp2-fp*dxm2)/det
        a2=(fp*dxm-fm*dxp)/det
      else
        ifail=1
        return
      endif

      if (a2.ne.0.0d0) then
        dxmax=-a1/(2.0d0*a2)
        dymax=(a1+a2*dxmax)*dxmax
        xopt=x0+dxmax
        yopt=f0+dymax
      endif

c calculate f=f0+a1*dx+a2*dx**2
c = a1*x - a1*x0 + a2*x**2 + a2*x0**2 - 2*a2*x*x0
c  f = f0 + (a2*x0 -a1)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

      a22=2.0d0*a2

      a(1)=f0 + (a2*x0 -a1)*x0
      a(2)=a1 - a22*x0
      a(3)=a2

c calculate yp=a1+2*a2*dx

      yp(1)=a1+a22*dxm
      yp(2)=a1
      yp(3)=a1+a22*dxp

      RETURN
      END
*CMZ :  4.00/11 28/05/2021  09.19.32  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.18.57  by  Michael Scheer
*CMZ :  2.68/05 03/09/2012  09.27.27  by  Michael Scheer
*CMZ : 00.00/11 11/02/2011  15.34.09  by  Michael Scheer
*-- Author :    Michael Scheer   11/02/2011
      SUBROUTINE util_spline_running_integral(X,Y,N,RESULT
     &                                 ,COEF,WORK1,WORK2,WORK3,WORK4)
*KEEP,gplhint.
*KEND.

C---  CALCULATES RUNNING INTERGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N
      double precision X(N),Y(N),RESULT(n)
      double precision COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF(X,Y,N,-9999.0d0,-9999.0d0,COEF,
     &  WORK1,WORK2,WORK3,WORK4)

C--- INTEGRATION

      RESULT(1)=0.0D0
      DO I=1,N-1

        RESULT(i+1)=RESULT(i)
     &          +(X(I+1)-X(I))*0.5D0
     &          *(Y(I)+Y(I+1))
     &          -(X(I+1)-X(I))**3/24.D0
     &          *(COEF(I)+COEF(I+1))

      ENDDO

      RETURN
      END
*CMZ :  3.05/05 12/07/2018  13.00.23  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.51.04  by  Michael Scheer
*CMZ :  3.02/05 27/03/2015  15.15.15  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.19.21  by  Michael Scheer
*CMZ :  2.68/05 03/09/2012  09.26.37  by  Michael Scheer
*-- Author :    Michael Scheer   10/05/2012
      subroutine util_g1_static(y,g1)
*KEEP,gplhint.
*KEND.

c calculates G1(y) with an estimated precision of about 1.5e-3 for y<=30,
c and 1.5e-2 for y>30.

      implicit none

      integer ical,npoi,ipoi,npoilow,npoihigh,ndatp
      parameter(ndatp=29)

      double precision y,g1,g1_30,c_30,g1_5em5,c_5em5,y_30,
     &  y_5em5,ydum,g1dum,r1

      double precision
     &  ywlow(ndatp),ywhigh(ndatp),coeflow(ndatp),coefhigh(ndatp),
     &  g1wlow(ndatp),g1whigh(ndatp),g1walow(ndatp),
     &  g1wahigh(ndatp),r1low(ndatp),r1high(ndatp),
     &  w1(ndatp),w2(ndatp),w3(ndatp),w4(ndatp),
     &  ystat(ndatp),g1stat(ndatp)

      data ical/0/
      data y_30/30.0d0/
      data y_5em5/5.0d-5/
      data g1_30/6.580794488121591d-013/ !WAVE
      data g1_5em5/7.909860755922665E-002/ !WAVE

* Numerisch mit WAVE berechnet 11.5.2012 (ISPECDIP=2)
        data ystat/
     &    5.0D-005, 7.0D-005, 2.0D-004, 5.0D-004, 1.0D-003,
     &    2.0D-003, 5.0D-003, 1.0D-002, 2.0D-002, 5.0D-002,
     &    0.10D0, 0.20D0, 0.50D0, 1.0D0, 2.0D0,
     &    3.0D0, 4.0D0, 5.0D0, 6.0D0, 7.0D0,
     &    8.0D0, 9.0D0, 10.00D0, 20.00D0, 30.00D0,
     &    40.0D0, 50.0D0, 60.0D0, 70.0D0
     &    /

        data g1stat/
     &    7.909860755922665D-002, 8.846122555950733D-002,
     &    0.125342415210665d0, 0.169701277907238d0, 0.2131391d0,
     &    0.2671962d0, 0.3584969d0, 0.4449725d0, 0.5472394d0,
     &    0.7015719d0, 0.8181855d0, 0.9033860d0, 0.8708191d0,
     &    0.65142282d0, 0.30163590d0, 0.128565710002655d0,
     &    5.282739666852105D-002,
     &    2.12481297D-002, 8.426079715722744D-003,
     &    3.307610970763407D-003, 1.288451614441198D-003,
     &    4.988932935072772D-004, 1.92238264D-004,
     &    1.19686345D-008, 6.58079455D-013, 3.42988745D-017,
     &    1.73478519D-021, 8.60693915D-026, 4.21333348D-030
     &    /

      save ical,c_5em5,c_30,
     &  ywlow,r1low,coeflow,npoilow,
     &  ywhigh,r1high,coefhigh,npoihigh

      if (ical.eq.0) then

        c_5em5=g1_5em5/y_5em5**(1.0d0/3.0d0)
        c_30=g1_30/(sqrt(y_30)*exp(-y_30))

        npoilow=0
        npoihigh=0

        do npoi=1,ndatp
          ydum=ystat(npoi)
          if (ydum.ge.y_5em5.and.ydum.le.4.0d0) then !zwei Abfragen wegen 4.0
            npoilow=npoilow+1
          endif
          if (ydum.ge.4.0d0.and.ydum.le.y_30) then
            npoihigh=npoihigh+1
          endif
        enddo

        npoi=ndatp

        npoilow=0
        npoihigh=0
        do ipoi=1,npoi
          ydum=ystat(ipoi)
          g1dum=g1stat(ipoi)
          if (ydum.ge.y_5em5.and.ydum.le.4.0d0) then
            npoilow=npoilow+1
            ywlow(npoilow)=ydum
            g1wlow(npoilow)=g1dum
            g1walow(npoilow)=
     &        391.8d0 * ydum**(1.0d0/3.0d0) * exp(-ydum*0.8307d0)
     &        -192.0d0 * sqrt(ydum) * exp(-ydum*0.7880d0)
          endif
          if (ydum.ge.4.0d0.and.ydum.le.y_30) then
            npoihigh=npoihigh+1
            ywhigh(npoihigh)=ydum
            g1whigh(npoihigh)=g1dum
            g1wahigh(npoihigh)=164.0d0*sqrt(ydum)*EXP(-ydum)
          endif
        enddo

        r1low(1:npoilow)=g1wlow(1:npoilow)/g1walow(1:npoilow)
        r1high(1:npoihigh)=g1whigh(1:npoihigh)/g1wahigh(1:npoihigh)

        call util_spline_coef(ywlow,r1low,npoilow,0.0d0,0.0d0,coeflow,
     &    w1,w2,w3,w4)
        call util_spline_coef(ywhigh,r1high,npoihigh,0.0d0,0.0d0,
     &    coefhigh,w1,w2,w3,w4)

        ical=1
      endif

      if (y.le.5.0d-5) then
        g1=c_5em5*y**(1.0d0/3.0d0)
      else if (y.ge.30.0d0) then
        g1=c_30*sqrt(y)*exp(-y)
      else

        if (y.ge.y_5em5.and.y.lt.4.0d0) then
          call util_spline_inter(ywlow,r1low,coeflow,npoilow,y,r1,-1)
          g1=r1*(
     &      391.8d0 * y**(1.0d0/3.0d0) * exp(-y*0.8307d0)
     &      -192.0d0 * sqrt(y) * exp(-y*0.7880d0))
        else if (y.ge.4.0d0.and.y.le.y_30) then
          call util_spline_inter(
     &      ywhigh,r1high,coefhigh,npoihigh,y,r1,-1)
          g1=r1*(164.0d0*sqrt(y)* EXP(-y))
        endif

      endif

      return
9999  stop '*** File wave-g1.dat not found ***'
      end
*CMZ : 00.00/15 05/01/2012  13.52.39  by  Michael Scheer
*CMZ : 00.00/00 11/01/95  11.41.04  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SORT_FUNC(N,RA,YA)

C--- HEAPSORT ROUTINE; SEE NUMERICAL RECIPES 8.2 S 231
C--- ARRAY YA IS FUNCTION OF RA AND SORTED ACCORDINGLY

      IMPLICIT NONE

      INTEGER N,L,IR,I,J

      REAL*8 RA(N),RRA
      REAL*8 YA(N),YYA

      if (n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          YYA=YA(L)
        ELSE
          RRA=RA(IR)
          YYA=YA(IR)
          RA(IR)=RA(1)
          YA(IR)=YA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            YA(1)=YYA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            YA(I)=YA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        YA(I)=YYA
      GO TO 10
      END
*CMZ :          15/04/2024  15.07.21  by  Michael Scheer
*CMZ :  4.01/02 11/05/2023  12.10.27  by  Michael Scheer
*-- Author :    Michael Scheer   08/05/2023
      subroutine util_get_electron(xbeta,betah,alphah,betav,alphav,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  espread,bunchlen,xelec,yelec,zelec,ypelec,zpelec,deelec,modebunch)

      implicit none

      double precision xbeta,betah,alphah,betav,alphav,emith,emitv,
     &  ebeam,xelec,yelec,zelec,ypelec,zpelec,deelec,
     &  s0h,beta0h,gamma0h,espread,bunchlen,
     &  s1h,beta1h,betap1h,alpha1h,gamma1h,phase1h,
     &  s2h,beta2h,betap2h,alpha2h,gamma2h,phase2h,
     &  s0v,beta0v,gamma0v,
     &  s1v,beta1v,betap1v,alpha1v,gamma1v,pvase1v,
     &  s2v,beta2v,betap2v,alpha2v,gamma2v,pvase2v,
     &  sigz,sigzp,sigy,sigyp,
     &  disph,dispph,dispv,disppv

      integer,parameter :: ng=6
      real rr(2),g(ng)

      integer :: modebunch

      if (emith.ne.0.0d0) then

        s1h=xbeta
        s2h=xelec
        beta1h=betah
        betap1h=-2.0d0*alphah

        call util_beta_function_drift(
     &    s0h,beta0h,gamma0h,
     &    s1h,beta1h,betap1h,alpha1h,gamma1h,phase1h,
     &    s2h,beta2h,betap2h,alpha2h,gamma2h,phase2h)

        sigz=sqrt(emith*beta0h)
        sigzp=sqrt(emith/beta0h)

      else

        sigz=0.0d0
        sigzp=0.0d0

      endif

      if (emitv.ne.0.0d0) then

        s1v=xbeta
        s2v=xelec
        beta1v=betav
        betap1v=-2.0d0*alphav

        call util_beta_function_drift(
     &    s0v,beta0v,gamma0v,
     &    s1v,beta1v,betap1v,alpha1v,gamma1v,pvase1v,
     &    s2v,beta2v,betap2v,alpha2v,gamma2v,pvase2v)

        sigy=sqrt(emitv*beta0v)
        sigyp=sqrt(emitv/beta0v)

      else

        sigy=0.0d0
        sigyp=0.0d0

      endif

      call util_random_gauss(ng,g,rr)

      deelec=g(1)*espread
      xelec=xelec+g(2)*bunchlen

      zelec=g(3)*sigz+deelec*disph
      zpelec=g(4)*sigzp+deelec*dispph

      yelec=g(5)*sigy+deelec*dispv
      ypelec=g(6)*sigyp+deelec*disppv

      zelec=zelec+(xelec-s0h)*zpelec
      yelec=yelec+(xelec-s0v)*ypelec

      return
      end
*CMZ : 00.00/15 09/10/2012  13.45.25  by  Michael Scheer
*-- Author :    Michael Scheer   05/10/2012
      subroutine util_beta_function_drift(
     &  s0,beta0,gamma0,
     &  s1,beta1,betap1,alpha1,gamma1,phase1,
     &  s2,beta2,betap2,alpha2,gamma2,phase2)

      real*8 s1,beta1,betap1,alpha1,gamma1,phase1,phase2,
     &  s2,beta2,betap2,alpha2,gamma2,beta0,gamma0,s0

c Calculates beta(s2) etc. from beta(s1) and betap(s1)
c beta(s)=beta0+(s-s0)**2/beta(0)

      alpha1=-betap1/2.0d0
      gamma1=(1.0d0+alpha1**2)/beta1
      s0=s1+alpha1/gamma1
      beta0=1.0d0/gamma1
      gamma0=1.0d0/beta0
      phase1=atan((s1-s0)/beta0)

      beta2=beta0+(s2-s0)**2/beta0
      betap2=2.0d0*(s2-s0)/beta0
      alpha2=-(s2-s0)/beta0
      gamma2=(1.0d0+alpha2**2)/beta2
      phase2=atan((s2-s0)/beta0)

      return
      end
*CMZ :  4.00/15 27/04/2022  08.09.29  by  Michael Scheer
*CMZ :  3.06/00 14/01/2019  17.27.04  by  Michael Scheer
*CMZ :  3.03/02 18/01/2016  13.02.27  by  Michael Scheer
*CMZ :  3.02/03 30/10/2014  17.17.28  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random_gauss(n,g,rr)

      implicit none

c Based on textbook "Numerical Recipies"
c The subroutine util_random is called, which can be initialized
c and controlled by util_random_init, util_random_set_seed, and
c util_random_get_seed.

      real g(n),r,v1,v2,fac,rr(2)
      integer n,i

      if (n.eq.1) then
1       call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 1
        fac=sqrt(-2.*log(r)/r)
        g(1)=v1*fac
      else
        do i=1,(n/2*2),2
11      call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 11
        fac=sqrt(-2.*log(r)/r)
        g(i)=v1*fac
        g(i+1)=v2*fac
        enddo
      endif

      if (mod(n,2).ne.0) then
12      call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 12
        fac=sqrt(-2.*log(r)/r)
        g(n)=v1*fac
      endif

      return
      end
*CMZ :  4.01/03 17/05/2023  11.24.58  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  11.49.33  by  Michael Scheer
*CMZ : 00.00/16 21/11/2014  14.53.59  by  Michael Scheer
*-- Author :    Michael Scheer   21/11/2014
      subroutine util_spline_integral_2d(nx,ny,x,y,f,result,istat,kalloc)
*KEEP,gplhint.
*KEND.

      implicit none

      double precision x(nx),y(ny),f(nx,ny),result
      integer istat,nx,ny,ix,iy,kstat,kalloc

      double precision, allocatable :: fb(:),fb2(:),coef(:),
     &  w1(:),w2(:),w3(:),w4(:)

      save

      if (kalloc.gt.0) then
        allocate(fb(max(nx,ny)))
        allocate(fb2(max(nx,ny)))
        allocate(coef(max(nx,ny)))
        allocate(w1(max(nx,ny)))
        allocate(w2(max(nx,ny)))
        allocate(w3(max(nx,ny)))
        allocate(w4(max(nx,ny)))
      else if (kalloc.lt.0) then
        deallocate(fb)
        deallocate(fb2)
        deallocate(coef)
        deallocate(w1)
        deallocate(w2)
        deallocate(w3)
        deallocate(w4)
        return
      endif

      kstat=0

      if (ny.gt.nx) then
        do ix=1,nx
          fb(1:ny)=f(ix,1:ny)
          call util_spline_integral_stat(y,fb,ny,fb2(ix)
     &      ,coef,w1,w2,w3,w4,istat)
          kstat=kstat+istat
        enddo
        call util_spline_integral_stat(x,fb2,nx,result
     &    ,coef,w1,w2,w3,w4,istat)
        kstat=kstat+istat
      else !nx.gt.ny?
        do iy=1,ny
          fb(1:nx)=f(1:nx,iy)
          call util_spline_integral_stat(x,fb,nx,fb2(iy)
     &      ,coef,w1,w2,w3,w4,istat)
          kstat=kstat+istat
        enddo
        call util_spline_integral_stat(y,fb2,ny,result
     &    ,coef,w1,w2,w3,w4,istat)
        kstat=kstat+istat
      endif !nx.gt.ny

      return
      end
*CMZ :  4.01/03 16/05/2023  19.38.31  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_INTEGRAL_STAT(X,Y,N,RESULT
     &                                 ,COEF,WORK1,WORK2,WORK3,WORK4,ISTAT)

C---  CALCULATES INTERGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N,ISTAT
      REAL*8 X(N),Y(N),RESULT
      REAL*8 COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF_STATus(X,Y,N,-9999.0d0,-9999.0d0,COEF,
     &  WORK1,WORK2,WORK3,WORK4,ISTAT)

C--- INTEGRATION

      RESULT=0.0D0
      DO I=1,N-1

      RESULT=RESULT
     &          +(X(I+1)-X(I))*0.5D0
     &          *(Y(I)+Y(I+1))
     &          -(X(I+1)-X(I))**3/24.D0
     &          *(COEF(I)+COEF(I+1))

      ENDDO

      RETURN
      END
*CMZ :  4.01/03 17/05/2023  11.21.58  by  Michael Scheer
*CMZ : 00.00/20 18/11/2016  15.04.18  by  Michael Scheer
*CMZ : 00.00/19 19/11/2015  13.56.50  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.67/05 16/05/2012  12.45.37  by  Michael Scheer
*CMZ : 00.00/07 12/10/2009  12.17.45  by  Michael Scheer
*CMZ : 00.00/02 14/04/2003  12.46.09  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.48  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_COEF_STATUS(X,Y,N,YP1,YPN,Y2,AA,BB,CC,C,istat)

C--- CALCULATES SPLINE COEFFICIENTS

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       X: ARRAY OF X-VALUES
C-       Y: ARRAY OF Y-VALUES
C-       YP1:  SECOND DERIVATIVE AT FIRST X-VALUE
C-       YPN:  SECOND DERIVATIVE AT LAST X-VALUE

C--   OUPUT:

C-       Y2:   SPLINE-COEFFICIENTS
c:       istat: STATUS

C--   WORKINGSPACE: AA(N),BB(N),CC(N),C(N)


      IMPLICIT NONE

      INTEGER N,J
      REAL*8  X(N),Y(N),Y2(N),AA(N),BB(N),CC(N),C(N)

      REAL*8 YP1,YPN

      double precision xx(3),yy(3),a(3),yp(3),xopt,yopt
      INTEGER ifail,istat

      istat=0

      IF (N.LT.3) then
        if (abs(yp1).eq.9999.0d0) then
          y2(1)=0.0d0
        else
          y2(1)=yp1
        endif
        if (abs(ypn).eq.9999.0d0) then
          y2(n)=0.0d0
        else
          y2(n)=ypn
        endif
        RETURN
      endif

      if (abs(yp1).eq.9999.0d0) then
        xx=x(1:3)
        yy=y(1:3)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(1)=2.0d0*a(3)
        else
          print*,"*** Warning in util_spline_coef: Calculation of yp1 by util_parabel failed ***"
          print*,"*** Setting second derivative of first point to zero ***"
          y2(1)=0.0d0
        endif
      else
        Y2(1)=YP1
      endif

      if (abs(ypn).eq.9999.0d0) then
        xx=x(n-2:n)
        yy=y(n-2:n)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(n)=2.0d0*a(3)
        else
          print*,"*** Warning in util_spline_coef: Calculation of ypn by util_parabel failed ***"
          print*,"*** Setting second derivative of last point to zero ***"
          y2(N)=0.0d0
        endif
      else
        Y2(N)=YPN
      endif

      C(1)=Y2(1)
      C(N)=y2(n)

      BB(1)=1.D0
      CC(1)=0.D0
      CC(N)=1.D0

      DO J=2,N-1
        if(x(j+1).le.x(j)) then
          write(6,*)
          write(6,*)
     &      '*** Error in util_spline_coef: Intervall of zero length or bad ordering of data'
          write(6,*)'j, x(j), x(j+1):',j,x(j),x(j+1)
          write(6,*)
          istat=-1
          return
        endif

        AA(J)=(X(J  )-X(J-1))/6.D0
        BB(J)=(X(J+1)-X(J-1))/3.D0
        CC(J)=(X(J+1)-X(J  ))/6.D0
        C(J)=(Y(J+1)-Y(J  ))/(X(J+1)-X(J  ))
     &    -(Y(J  )-Y(J-1))/(X(J  )-X(J-1))
      ENDDO !J

      DO J=2,N-1

        BB(J)=BB(J)-AA(J)*CC(J-1)
        C(J)= C(J)-AA(J)* C(J-1)
C          AA(J)=AA(J)-AA(J)*BB(J-1)

        CC(J)=CC(J)/BB(J)
        C(J)= C(J)/BB(J)
        BB(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
        Y2(J)=C(J)-CC(J)*Y2(J+1)
        if (Abs(y2(j)).lt.1.0d-15) y2(j)=0.0d0
      ENDDO

      RETURN
      END
*CMZ :  3.03/04 06/11/2017  16.07.30  by  Michael Scheer
*CMZ :  3.03/02 20/01/2016  10.34.45  by  Michael Scheer
*CMZ :  3.02/03 27/10/2014  10.39.03  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random_set_seed(isize,iseed)

c The actual size of the seed array is at least 64
c If iseed(1:n) are set zero, the behaviour is unclear to me, but it seems
c to reduce the used seed array to iseed(n+1:64), but what else is changed?
c So be careful!

      use iso_fortran_env !, only: int64

      implicit none

      integer isize
      integer iseed(isize),myseed(64)

      if (isize.lt.64) then
        myseed(1:64-isize)=0
        myseed(isize+1:64)=iseed(1:isize)
        call random_seed(put=myseed)
      else
        call random_seed(put=iseed)
      endif

      return
      end
*CMZ :  4.01/04 22/01/2018  16.49.42  by  Michael Scheer
*CMZ : 00.00/20 19/08/2016  15.08.19  by  Michael Scheer
*CMZ : 00.00/15 04/01/2013  12.22.07  by  Michael Scheer
*CMZ : 00.00/05 27/02/2007  16.32.04  by  Michael Scheer
*CMZ : 00.00/02 04/08/2006  14.56.41  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.58  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_zeit_kommentar(lun,comment)

c to determine date and time and write it to logical unit lun

      implicit none

      integer lun

      character(*) comment

      character spacer(50)
      character(10) dtday,dttime,dtzone
      integer idatetime(8),ilast

      data spacer/50*' '/

      call date_and_time(dtday,dttime,dtzone,idatetime)

      ilast=len_trim(comment)

      write(lun,*)
      write(lun,*)
      if (ilast.gt.0) then
        write(lun,*)comment(1:ilast),spacer,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      else
        write(lun,*)spacer,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      endif
      write(lun,*)

      return
      end
*CMZ :          16/04/2024  12.38.45  by  Michael Scheer
*CMZ :  2.05/04 16/12/2023  12.06.33  by  Michael Scheer
*CMZ :  2.04/19 16/09/2023  16.33.50  by  Michael Scheer
*CMZ :  2.03/00 26/07/2022  07.55.50  by  Michael Scheer
*CMZ :  2.02/02 01/07/2022  17.30.28  by  Michael Scheer
*CMZ :  2.02/00 29/03/2021  09.26.44  by  Michael Scheer
*CMZ :  1.00/00 19/08/2016  18.24.11  by  Michael Scheer
*CMZ : 00.00/15 04/01/2013  12.22.07  by  Michael Scheer
*CMZ : 00.00/05 27/02/2007  16.32.04  by  Michael Scheer
*CMZ : 00.00/02 04/08/2006  14.56.41  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.58  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_zeit_kommentar_delta(lun,comment,iset)

c to determine date and time and write it to logical unit lun

      implicit none

      integer lun,ilast

      character(*) comment

      character spacer(30)
      character(10) dtday,dttime,dtzone
      character(32) c32
      character(2048) :: cblank=''
      character(2048) cline

      real :: seconds=0.0,secondso=0.0
      integer idatetime(8),iyear,imonth,iday,ihour,iminute,isec
      integer iyearo,imontho,idayo,ihouro,iminuteo,iseco,kday,khour,kminit,ksec

      integer :: ical=0,linlen,iset,iseconds,isecondso

      data spacer/30*' '/

      save

      call date_and_time(dtday,dttime,dtzone,idatetime)

      iyear=idatetime(1)
      imonth=idatetime(2)
      iday=idatetime(3)
      ihour=idatetime(5)
      iminute=idatetime(6)
      isec=idatetime(7)

      ilast=len_trim(comment)

      write(lun,*)
      if (ilast.gt.0) then
        write(cline,*)comment(1:ilast),spacer,dttime(1:2),':',dttime(3:4),':',
     &    dttime(5:6),' ',dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      else
        write(cline,*)spacer,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      endif

      ilast=len_trim(cline)

      if (ilast.lt.64) then
        write(lun,'(a)') cblank(1:64-ilast) // cline(1:ilast)
      else
        write(lun,'(a)') trim(cline)
      endif

      if (ical.gt.0.and.iset.eq.0) then
        cline='  delta time: '
        seconds=secnds(secondso)
        iseconds=nint(seconds)
        kday=0
        khour=0
        kminit=0
        ksec=0
        if (iseconds.ge.3600*24) then
          kday=iseconds/(3600*24)
          iseconds=iseconds-kday*3600*24
          write(c32,*)kday
          cline=trim(cline) // trim(c32) // ' days '
        endif
        if (iseconds.ge.3600) then
          khour=iseconds/3600
          iseconds=iseconds-khour*3600
          write(c32,*)khour
          cline=trim(cline) // trim(c32) // ' hours '
        endif
        if (iseconds.ge.60) then
          kminit=iseconds/60
          iseconds=iseconds-kminit*60
          write(c32,*)kminit
          cline=trim(cline) // trim(c32) // ' minutes '
        endif
        write(c32,*)iseconds
        cline=trim(cline) // trim(c32) // ' seconds'
        write(lun,*) trim(cline)
      endif

      write(lun,*)''

      if (ical.eq.0.or.iset.ne.0) then
        secondso=secnds(0.0)
        isecondso=nint(secondso)
        iyearo=iyear
        imontho=imonth
        idayo=iday
        ihouro=ihour
        iminuteo=iminute
        iseco=isec
      endif

      ical=1

      return
      end
*CMZ :  3.03/02 17/11/2015  14.46.41  by  Michael Scheer
*CMZ :  3.02/03 05/09/2014  12.38.49  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random_get_seed(isize,iseed)

      use iso_fortran_env !, only: int64

      implicit none

      integer isize
      integer iseed(isize)

      call random_seed(get=iseed)

      return
      end
*CMZ :  3.06/00 22/01/2019  16.22.02  by  Michael Scheer
*CMZ :  3.03/02 17/11/2015  14.46.41  by  Michael Scheer
*CMZ :  3.02/03 05/09/2014  12.37.33  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random_init(isize,iseed)

c isize contains size if seed-array to be provided by the user
c For a simple use of the generated array call random_seed(isize), where
c isize is an integer

C IT SEEMS NOT TO WORK AS EXPECTED, MAYBE THE FIRST HALF OF THE SEEDS MUST BE
C SAVED TO THE SECOND HALF OF THE ARRAY!?

      use iso_fortran_env !, only: int64

      implicit none

      integer isize,n
      integer iseed(isize)

      n=isize
      call random_seed(size=isize)

C      if (isize.ne.64) then
      if (n.ne.64) then
        iseed(1)=-1
        print*,
     &    '*** Error in util_random_init: Dimension of seed-array must be 64:'
c     &    ,isize
        return
      endif

      call random_seed(get=iseed)
      iseed(33:64)=iseed(1:32)

      isize=n

      return
      end
*CMZ :          12/01/2024  16.28.35  by  Michael Scheer
*-- Author :    Michael Scheer   11/11/2023

      SUBROUTINE util_CFT(A,B,NTOT,N,NSPAN,ISN,maxf,maxp,at,ck,bt,sk,np)
      ! based on cft.F of cernlib D702

C
C     MULTIVARIATE COMPLEX FOURIER TRANSFORM, COMPUTED IN PLACE
C     USING MIXED-RADIX FAST FOURIER TRANSFORM ALGORITHM.
C     BY R. C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968
C     ARRAYS A AND B ORIGINALLY HOLD THE REAL AND IMAGINARY
C     COMPONENTS OF THE DATA, AND RETURN THE REAL AND
C     IMAGINARY COMPONENTS OF THE RESULTING FOURIER COEFFICIENTS.
C     MULTIVARIATE DATA IS INDEXED ACCORDING TO THE FORTRAN
C     ARRAY ELEMENT SUCCESSOR FUNCTION, WITHOUT LIMIT
C     ON THE NUMBER OF IMPLIED MULTIPLE SUBSCRIPTS.
C     THE SUBROUTINE IS CALLED ONCE FOR EACH VARIATE.
C     THE CALLS FOR A MULTIVARIATE TRANSFORM MAY BE IN ANY ORDER.
C     NTOT IS THE TOTAL NUMBER OF COMPLEX DATA VALUES.
C     N IS THE DIMENSION OF THE CURRENT VARIABLE.
C     NSPAN/N IS THE SPACING OF CONSUCUTIVE DATA VALUES
C     WHILE INDEXING THE CURRENT VARIABLE.
C     THE SIGN OF ISN DETERMINES THE SIGN OF THE COMPLEX
C     EXPONENTIAL, AND THE MAGNITUDE OF ISN IS NORMALLY ONE.
C
C     FOR A SINGLE-VARIATE TRANSFORM,
C     NTOT = N = NSPAN = (NUMBER OF COMPLEX DATA VALUES), F.G.
C     CALL CFT(A,B,N,N,N,1)
C
C     A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3), B(N1,N2,N3)
C     IS COMPUTED BY
C     CALL CFT(A,B,N1*N2*N3,N1,N1,1)
C     CALL CFT(A,B,N1*N2*N3,N2,N1*N2,1)
C     CALL CFT(A,B,N1*N2*N3,N3,N1*N2*N3,1)
C
C     THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
C     ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
C     GIVE THE CORRECT INDEXING INCREMENT AND THE SECOND PARAMETER
C     USED TO PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF
C     IMAGINARY VALUES, E.G.
C
C        REAL S(2)
C        EQUIVALENCE (A,S)
C        ....
C        ....
C        CALL CFT(A,S(2),NTOT,N,NSPAN,2)
C
C     ARRAYS AT(MAXF), CK(MAXF), BT(MAXF), SK(MAXF), AND NP(MAXP)
C     ARE USED FOR TEMPORARY STORAGE. IF THE AVAILABLE STORAGE
C     IS INSUFFICIENT, THE PROGRAM IS TERMINATED BY A STOP.
C     MAXF MUST BE .GE. THE MAXIMUM PRIME FACTOR OF N.
C     MAXP MUST BE .GT. THE NUMBER OF PRIME FACTORS OF N.
C     IN ADDITION, IF THE SQUARE-FREE PORTION K CF N HAS TWO OR
C     MORE PRIME FACTORS, THEN MAXP MUST BE .GE. K-1.
C     ARRAY STORAGE IN NFAC FOR A MAXIMUM OF 11 FACTORS OF N.
C     IF N HAS MORE THAN ONE SQUARE-FREE FACTOR, THE PRODUCT OF THE
C     SQUARE-FREE FACTORS MUST BE .LE. 210
C
      implicit none

      double precision A(ntot),B(ntot)
      double precision  AT(maxf),CK(maxf),BT(maxf),SK(maxf),
     &  rad,s72,c72,s120,sd,s3,aa,aj,ajm,ajp,ak,akm,akp,bb,bj,bjm,bjp,bk,bkm,
     &  bkp,c1,cd,c2,c3,radf,s1,s2
cmsh      integer NFAC(11),NP(209)
      integer maxp,maxf,n,inc,isn,nspan,i,ii,ntot,j,jc,jf,jj,k,k1,k2,k3,
     &  k4,kk,ks,kspan,kspnn,kt,m,nn,nt
      integer NFAC(11),NP(maxp)
C     ARRAY STORAGE FOR MAXIMUM PRIME FACTOR OF 23
cmsh      double precision  AT(23),CK(23),BT(23),SK(23)
      EQUIVALENCE (I,II)
C     THE FOLLOWING TWO CONSTANTS SHOULD AGREE WITH THE ARRAY DIMENSIONS
cmsh      MAXF=23
cmsh      MAXP=209
      IF(N .LT. 2) RETURN
      INC=ISN
C     THE FOLLOWING CONSTANTS ARE RAD = 2.*PI , S72 = SIN(0.4*PI) ,
C     C72 = COS(0.4*PI) AND S120 = SQRT(0.75)
      RAD = 6.2831853071796
      S72 = 0.95105651629515
      C72 = 0.30901699437495
      S120 = 0.86602540378444
      IF(ISN .GE. 0) GO TO 10
      S72=-S72
      S120=-S120
      RAD=-RAD
      INC=-INC
   10 NT=INC*NTOT
      KS=INC*NSPAN
      KSPAN=KS
      NN=NT-INC
      JC=KS/N
      RADF=RAD*JC*0.5
      I=0
      JF=0
C     DETERMINE THE FACTORS OF N
      M=0
      K=N
      GO TO 20
   15 M=M+1
      NFAC(M)=4
      K=K/16
   20 IF(K-(K/16)*16 .EQ. 0) GO TO 15
      J=3
      JJ=9
      GO TO 30
   25 M=M+1
      NFAC(M)=J
      K=K/JJ
   30 IF(MOD(K,JJ) .EQ. 0) GO TO 25
      J=J+2
      JJ=J**2
      IF(JJ .LE. K) GO TO 30
      IF(K .GT. 4) GO TO 40
      KT=M
      NFAC(M+1)=K
      IF(K .NE. 1) M=M+1
      GO TO 80
   40 IF(K-(K/4)*4 .NE. 0) GO TO 50
      M=M+1
      NFAC(M)=2
      K=K/4
   50 KT=M
      J=2
   60 IF(MOD(K,J) .NE. 0) GO TO 70
      M=M+1
      NFAC(M)=J
      K=K/J
   70 J=((J+1)/2)*2+1
      IF(J .LE. K) GO TO 60
   80 IF(KT .EQ. 0) GO TO 100
      J=KT
   90 M=M+1
      NFAC(M)=NFAC(J)
      J=J-1
      IF(J .NE. 0) GO TO 90
C     COMPUTE FOURIER TRANSFORM
  100 SD=RADF/KSPAN
      CD=2.0*SIN(SD)**2
      SD=SIN(SD+SD)
      KK=1
      I=I+1
      IF(NFAC(I) .NE. 2) GO TO 400
C     TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)
      KSPAN=KSPAN/2
      K1=KSPAN+2
  210 K2=KK+KSPAN
      AK=A(K2)
      BK=B(K2)
      A(K2)=A(KK)-AK
      B(K2)=B(KK)-BK
      A(KK)=A(KK)+AK
      B(KK)=B(KK)+BK
      KK=K2+KSPAN
      IF(KK .LE. NN) GO TO 210
      KK=KK-NN
      IF(KK .LE. JC) GO TO 210
      IF(KK .GT. KSPAN) GO TO 800
  220 C1=1.0-CD
      S1=SD
  230 K2=KK+KSPAN
      AK=A(KK)-A(K2)
      BK=B(KK)-B(K2)
      A(KK)=A(KK)+A(K2)
      B(KK)=B(KK)+B(K2)
      A(K2)=C1*AK-S1*BK
      B(K2)=S1*AK+C1*BK
      KK=K2+KSPAN
      IF(KK .LT. NT) GO TO 230
      K2=KK-NT
      C1=-C1
      KK=K1-K2
      IF(KK .GT. K2) GO TO 230
      AK=C1-(CD*C1+SD*S1)
      S1=(SD*C1-CD*S1)+S1
C     THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C     ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY BE DELETED.
C     C1=0.5/(AK**2+S1**2)+0.5
C     S1=C1*S1
C     C1=C1*AK
C     NEXT STATEMENT SHOULD BE DELETED IF NON-ROUNDED ARITHMETIC IS USED
      C1=AK
      KK=KK+JC
      IF(KK .LT. K2) GO TO 230
      K1=K1+INC+INC
      KK=(K1-KSPAN)/2+JC
      IF(KK .LE. JC+JC) GO TO 220
      GO TO 100
C     TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)
  320 K1=KK+KSPAN
      K2=K1+KSPAN
      AK=A(KK)
      BK=B(KK)
      AJ=A(K1)+A(K2)
      BJ=B(K1)+B(K2)
      A(KK)=AK+AJ
      B(KK)=BK+BJ
      AK=-0.5*AJ+AK
      BK=-0.5*BJ+BK
      AJ=(A(K1)-A(K2))*S120
      BJ=(B(K1)-B(K2))*S120
      A(K1)=AK-BJ
      B(K1)=BK+AJ
      A(K2)=AK+BJ
      B(K2)=BK-AJ
      KK=K2+KSPAN
      IF(KK .LT. NN) GO TO 320
      KK=KK-NN
      IF(KK .LE. KSPAN) GO TO 320
      GO TO 700
C     TRANSFORM FOR FACTOR OF 4
  400 IF(NFAC(I) .NE. 4) GO TO 600
      KSPNN=KSPAN
      KSPAN=KSPAN/4
  410 C1=1.0
      S1=0
  420 K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
      AJP=AKP-AJP
      BKP=B(KK)+B(K2)
      BKM=B(KK)-B(K2)
      BJP=B(K1)+B(K3)
      BJM=B(K1)-B(K3)
      B(KK)=BKP+BJP
      BJP=BKP-BJP
      IF(ISN .LT. 0) GO TO 450
      AKP=AKM-BJM
      AKM=AKM+BJM
      BKP=BKM+AJM
      BKM=BKM-AJM
      IF(S1 .EQ. 0.0) GO TO 460
  430 A(K1)=AKP*C1-BKP*S1
      B(K1)=AKP*S1+BKP*C1
      A(K2)=AJP*C2-BJP*S2
      B(K2)=AJP*S2+BJP*C2
      A(K3)=AKM*C3-BKM*S3
      B(K3)=AKM*S3+BKM*C3
      KK=K3+KSPAN
      IF(KK .LE. NT) GO TO 420
  440 C2=C1-(CD*C1+SD*S1)
      S1=(SD*C1-CD*S1)+S1
C     THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C     ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY BE DELETED.
C     C1=0.5/(C2**2+S1**2)+0.5
C     S1=C1*S1
C     C1=C1*C2
C     NEXT STATEMENT SHOULD BE DELETED IF NON-ROUNDED ARITHMETIC IS USED
      C1=C2
      C2=C1**2-S1**2
      S2=2.0*C1*S1
      C3=C2*C1-S2*S1
      S3=C2*S1+S2*C1
      KK=KK-NT+JC
      IF(KK .LE. KSPAN) GO TO 420
      KK=KK-KSPAN+INC
      IF(KK .LE. JC) GO TO 410
      IF(KSPAN .EQ. JC) GO TO 800
      GO TO 100
  450 AKP=AKM+BJM
      AKM=AKM-BJM
      BKP=BKM-AJM
      BKM=BKM+AJM
      IF(S1 .NE. 0.0) GO TO 430
  460 A(K1)=AKP
      B(K1)=BKP
      A(K2)=AJP
      B(K2)=BJP
      A(K3)=AKM
      B(K3)=BKM
      KK=K3+KSPAN
      IF(KK .LE. NT) GO TO 420
      GO TO 440
C     TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)
  510 C2=C72**2-S72**2
      S2=2.0*C72*S72
  520 K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
      K4=K3+KSPAN
      AKP=A(K1)+A(K4)
      AKM=A(K1)-A(K4)
      BKP=B(K1)+B(K4)
      BKM=B(K1)-B(K4)
      AJP=A(K2)+A(K3)
      AJM=A(K2)-A(K3)
      BJP=B(K2)+B(K3)
      BJM=B(K2)-B(K3)
      AA=A(KK)
      BB=B(KK)
      A(KK)=AA+AKP+AJP
      B(KK)=BB+BKP+BJP
      AK=AKP*C72+AJP*C2+AA
      BK=BKP*C72+BJP*C2+BB
      AJ=AKM*S72+AJM*S2
      BJ=BKM*S72+BJM*S2
      A(K1)=AK-BJ
      A(K4)=AK+BJ
      B(K1)=BK+AJ
      B(K4)=BK-AJ
      AK=AKP*C2+AJP*C72+AA
      BK=BKP*C2+BJP*C72+BB
      AJ=AKM*S2-AJM*S72
      BJ=BKM*S2-BJM*S72
      A(K2)=AK-BJ
      A(K3)=AK+BJ
      B(K2)=BK+AJ
      B(K3)=BK-AJ
      KK=K4+KSPAN
      IF(KK .LT. NN) GO TO 520
      KK=KK-NN
      IF(KK .LE. KSPAN) GO TO 520
      GO TO 700
C     TRANSFORM FOR ODD FACTORS
  600 K=NFAC(I)
      KSPNN=KSPAN
      KSPAN=KSPAN/K
      IF(K .EQ. 3) GO TO 320
      IF(K .EQ. 5) GO TO 510
      IF(K .EQ. JF) GO TO 640
      JF=K
      S1=RAD/K
      C1=COS(S1)
      S1=SIN(S1)
      IF(JF .GT. MAXF) GO TO 998
      CK(JF)=1.0
      SK(JF)=0.0
      J=1
  630 CK(J)=CK(K)*C1+SK(K)*S1
      SK(J)=CK(K)*S1-SK(K)*C1
      K=K-1
      CK(K)=CK(J)
      SK(K)=-SK(J)
      J=J+1
      IF(J .LT. K) GO TO 630
  640 K1=KK
      K2=KK+KSPNN
      AA=A(KK)
      BB=B(KK)
      AK=AA
      BK=BB
      J=1
      K1=K1+KSPAN
  650 K2=K2-KSPAN
      J=J+1
      AT(J)=A(K1)+A(K2)
      AK=AT(J)+AK
      BT(J)=B(K1)+B(K2)
      BK=BT(J)+BK
      J=J+1
      AT(J)=A(K1)-A(K2)
      BT(J)=B(K1)-B(K2)
      K1=K1+KSPAN
      IF(K1 .LT. K2) GO TO 650
      A(KK)=AK
      B(KK)=BK
      K1=KK
      K2=KK+KSPNN
      J=1
  660 K1=K1+KSPAN
      K2=K2-KSPAN
      JJ=J
      AK=AA
      BK=BB
      AJ=0.0
      BJ=0.0
      K=1
  670 K=K+1
      AK=AT(K)*CK(JJ)+AK
      BK=BT(K)*CK(JJ)+BK
      K=K+1
      AJ=AT(K)*SK(JJ)+AJ
      BJ=BT(K)*SK(JJ)+BJ
      JJ=JJ+J
      IF(JJ .GT. JF) JJ=JJ-JF
      IF(K .LT. JF) GO TO 670
      K=JF-J
      A(K1)=AK-BJ
      B(K1)=BK+AJ
      A(K2)=AK+BJ
      B(K2)=BK-AJ
      J=J+1
      IF(J .LT. K) GO TO 660
      KK=KK+KSPNN
      IF(KK .LE. NN) GO TO 640
      KK=KK-NN
      IF(KK .LE. KSPAN) GO TO 640
C     MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4)
  700 IF(I .EQ. M) GO TO 800
      KK=JC+1
  710 C2=1.0-CD
      S1=SD
  720 C1=C2
      S2=S1
      KK=KK+KSPAN
  730 AK=A(KK)
      A(KK)=C2*AK-S2*B(KK)
      B(KK)=S2*AK+C2*B(KK)
      KK=KK+KSPNN
      IF(KK .LE. NT) GO TO 730
      AK=S1*S2
      S2=S1*C2+C1*S2
      C2=C1*C2-AK
      KK=KK-NT+KSPAN
      IF(KK .LE. KSPNN) GO TO 730
      C2=C1-(CD*C1+SD*S1)
      S1=S1+(SD*C1-CD*S1)
C     THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C     ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY
C     BE DELETED.
C     C1=0.5/(C2**2+S1**2)+0.5
C     S1=C1*S1
C     C2=C1*C2
      KK=KK-KSPNN+JC
      IF(KK .LE. KSPAN) GO TO 720
      KK=KK-KSPAN+JC+INC
      IF(KK .LE. JC+JC) GO TO 710
      GO TO 100
C     PERMUTE THE RESULTS TO NORMAL ORDER---DONE IN TWO STAGES
C     PERMUTATION FOR SQUARE FACTORS OF N
  800 NP(1)=KS
      IF(KT .EQ. 0) GO TO 890
      K=KT+KT+1
      IF(M .LT. K) K=K-1
      J=1
      NP(K+1)=JC
  810 NP(J+1)=NP(J)/NFAC(J)
      NP(K)=NP(K+1)*NFAC(J)
      J=J+1
      K=K-1
      IF(J .LT. K) GO TO 810
      K3=NP(K+1)
      KSPAN=NP(2)
      KK=JC+1
      K2=KSPAN+1
      J=1
      IF(N .NE. NTOT) GO TO 850
C     PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE)
  820 AK=A(KK)
      A(KK)=A(K2)
      A(K2)=AK
      BK=B(KK)
      B(KK)=B(K2)
      B(K2)=BK
      KK=KK+INC
      K2=KSPAN+K2
      IF(K2 .LT. KS) GO TO 820
  830 K2=K2-NP(J)
      J=J+1
      K2=NP(J+1)+K2
      IF(K2 .GT. NP(J)) GO TO 830
      J=1
  840 IF(KK .LT. K2) GO TO 820
      KK=KK+INC
      K2=KSPAN+K2
      IF(K2 .LT. KS) GO TO 840
      IF(KK .LT. KS) GO TO 830
      JC=K3
      GO TO 890
C     PERMUTATION FOR MULTIVARIATE TRANSFORM
  850 K=KK+JC
  860 AK=A(KK)
      A(KK)=A(K2)
      A(K2)=AK
      BK=B(KK)
      B(KK)=B(K2)
      B(K2)=BK
      KK=KK+INC
      K2=K2+INC
      IF(KK .LT. K) GO TO 860
      KK=KK+KS-JC
      K2=K2+KS-JC
      IF(KK .LT. NT) GO TO 850
      K2=K2-NT+KSPAN
      KK=KK-NT+JC
      IF(K2 .LT. KS) GO TO 850
  870 K2=K2-NP(J)
      J=J+1
      K2=NP(J+1)+K2
      IF(K2 .GT. NP(J)) GO TO 870
      J=1
  880 IF(KK .LT. K2) GO TO 850
      KK=KK+JC
      K2=KSPAN+K2
      IF(K2 .LT. KS) GO TO 880
      IF(KK .LT. KS) GO TO 870
      JC=K3
  890 IF(2*KT+1 .GE. M) RETURN
      KSPNN=NP(KT+1)
C     PERMUTATION FOR SQUARE-FREE FACTORS OF N
      J=M-KT
      NFAC(J+1)=1
  900 NFAC(J)=NFAC(J)*NFAC(J+1)
      J=J-1
      IF(J .NE. KT) GO TO 900
      KT=KT+1
      NN=NFAC(KT)-1
      IF(NN .GT. MAXP) GO TO 998
      JJ=0
      J=0
      GO TO 906
  902 JJ=JJ-K2
      K2=KK
      K=K+1
      KK=NFAC(K)
  904 JJ=KK+JJ
      IF(JJ .GE. K2) GO TO 902
      NP(J)=JJ
  906 K2=NFAC(KT)
      K=KT+1
      KK=NFAC(K)
      J=J+1
      IF(J .LE. NN) GO TO 904
C     DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1
      J=0
      GO TO 914
  910 K=KK
      KK=NP(K)
      NP(K)=-KK
      IF(KK .NE. J) GO TO 910
      K3=KK
  914 J=J+1
      KK=NP(J)
      IF(KK .LT. 0) GO TO 914
      IF(KK .NE. J) GO TO 910
      NP(J)=-J
      IF(J .NE. NN) GO TO 914
      MAXF=INC*MAXF
C     REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES
      GO TO 950
  924 J=J-1
      IF(NP(J) .LT. 0) GO TO 924
      JJ=JC
  926 KSPAN=JJ
      IF(JJ .GT. MAXF) KSPAN=MAXF
      JJ=JJ-KSPAN
      K=NP(J)
      KK=JC*K+II+JJ
      K1=KK+KSPAN
      K2=0
  928 K2=K2+1
      AT(K2)=A(K1)
      BT(K2)=B(K1)
      K1=K1-INC
      IF(K1 .NE. KK) GO TO 928
  932 K1=KK+KSPAN
      K2=K1-JC*(K+NP(K))
      K=-NP(K)
  936 A(K1)=A(K2)
      B(K1)=B(K2)
      K1=K1-INC
      K2=K2-INC
      IF(K1 .NE. KK) GO TO 936
      KK=K2
      IF(K .NE. J) GO TO 932
      K1=KK+KSPAN
      K2=0
  940 K2=K2+1
      A(K1)=AT(K2)
      B(K1)=BT(K2)
      K1=K1-INC
      IF(K1 .NE. KK) GO TO 940
      IF(JJ .NE. 0) GO TO 926
      IF(J .NE. 1) GO TO 924
  950 J=K3+1
      NT=NT-KSPNN
      II=NT-INC+1
      IF(NT .GE. 0) GO TO 924
      RETURN
C     ERROR FINISH, INSUFFICIENT ARRAY STORAGE
  998 ISN=0
      WRITE(6,999)
  999 FORMAT('*** Error in util_cft: ARRAY BOUNDS EXCEEDED, change number of points to avoid too big a prime factor')
      END
*CMZ :          12/01/2024  18.00.46  by  Michael Scheer
*-- Author :    Michael Scheer   11/11/2023
      subroutine util_cft_2d(nx,ny,arin,aiin,arout,aiout,isn)
      implicit none

      integer, parameter :: maxfp=23,maxpp=209

      integer :: nx,ny,np(maxpp),maxf=maxfp,maxp=maxpp,isn,ksn
      double precision :: arin(nx,ny),aiin(nx,ny),arout(nx,ny),aiout(nx,ny)
      double precision :: at(maxfp),ck(maxfp),bt(maxfp),sk(maxfp)

      if (isn.ge.0) then
        ksn=1
      else
        ksn=-1 !backward transformation
      endif

      arout=arin
      aiout=aiin

      ! Check cernlib.pdf, D702
      call util_cft(arout,aiout,nx*ny,nx,nx,ksn,maxf,maxp,at,ck,bt,sk,np)
      if (ksn.eq.0) then
        isn=ksn
        return
      endif
      call util_cft(arout,aiout,nx*ny,ny,nx*ny,ksn,maxf,maxp,at,ck,bt,sk,np)
      isn=ksn

      return
      end
