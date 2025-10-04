*CMZ :          10/04/2019  09.43.55  by  aaa_license
*-- Author :    Michael Scheer   14/12/2017

! This software is under the GNU General Public License:

!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
!    The this software contains parts of the CERN PROGRAM LIBRARY.
!    See http://www.cern.ch
!******************************************************************************
*CMZ :          07/10/2014  13.43.31  by  Michael Scheer
*CMZ :  1.16/04 17/04/2014  11.31.05  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*
* $Id: abend.F,v 1.1.1.1 1996/02/15 17:50:37 mclareni Exp $
*
* $Log: abend.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*

      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  '*** Aborted by abend.f ***' !msh
      END
*CMZ :          07/10/2014  14.02.53  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
*
* $Id: bsir364.F,v 1.1.1.1 1996/04/01 15:02:07 mclareni Exp $
*
* $Log: bsir364.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*

      FUNCTION DBSIR3(X,NU)

*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*80 ERRTXT
      CHARACTER NAMEI*(*),NAMEK*(*),NAMEIE*(*),NAMEKE*(*)

      PARAMETER (NAMEI = 'BSIR3/DBSIR3', NAMEIE = 'EBSIR3/DEBIR3')
      PARAMETER (NAMEK = 'BSKR3/DBSKR3', NAMEKE = 'EBSKR3/DEBKR3')

      LOGICAL LEX

      DIMENSION BC(0:23,2),CC(0:15,2),PP(-2:2),GG(-2:2)

      PARAMETER (EPS = 1D-15)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (W3 = 1.73205 08075 68877 29D0)
      PARAMETER (G1 = 2.67893 85347 07747 63D0)
      PARAMETER (G2 = 1.35411 79394 26400 42D0)
      PARAMETER (PIH = PI/2, RPIH = 2/PI, RPI2 = 1/(2*PI))
      PARAMETER (C1 = 2*PI/(3*W3))
      PARAMETER (GM = 3*(1/G2-3/G1)/2, GP = (3/G1+1/G2)/2)

      DATA PP(-2) /-0.66666 66666 66666 67D0/
      DATA PP(-1) /-0.33333 33333 33333 33D0/
      DATA PP( 1) / 0.33333 33333 33333 33D0/
      DATA PP( 2) / 0.66666 66666 66666 67D0/

      DATA GG(-2) / 0.37328 21739 07395 23D0/
      DATA GG(-1) / 0.73848 81116 21648 31D0/
      DATA GG( 1) / 1.11984 65217 22185 68D0/
      DATA GG( 2) / 1.10773 21674 32472 47D0/

      DATA BC( 0,1) / 1.00458 61710 93207 35D0/
      DATA BC( 1,1) / 0.00467 34791 99873 60D0/
      DATA BC( 2,1) / 0.00009 08034 04815 04D0/
      DATA BC( 3,1) / 0.00000 37262 16110 59D0/
      DATA BC( 4,1) / 0.00000 02520 73237 90D0/
      DATA BC( 5,1) / 0.00000 00227 82110 77D0/
      DATA BC( 6,1) / 0.00000 00012 91332 28D0/
      DATA BC( 7,1) /-0.00000 00006 11915 16D0/
      DATA BC( 8,1) /-0.00000 00003 75616 85D0/
      DATA BC( 9,1) /-0.00000 00001 16415 46D0/
      DATA BC(10,1) /-0.00000 00000 14443 25D0/
      DATA BC(11,1) / 0.00000 00000 05373 69D0/
      DATA BC(12,1) / 0.00000 00000 03074 27D0/
      DATA BC(13,1) / 0.00000 00000 00297 66D0/
      DATA BC(14,1) /-0.00000 00000 00265 20D0/
      DATA BC(15,1) /-0.00000 00000 00091 37D0/
      DATA BC(16,1) / 0.00000 00000 00015 52D0/
      DATA BC(17,1) / 0.00000 00000 00014 12D0/
      DATA BC(18,1) /-0.00000 00000 00000 23D0/
      DATA BC(19,1) /-0.00000 00000 00001 98D0/
      DATA BC(20,1) /-0.00000 00000 00000 13D0/
      DATA BC(21,1) / 0.00000 00000 00000 29D0/
      DATA BC(22,1) / 0.00000 00000 00000 03D0/
      DATA BC(23,1) /-0.00000 00000 00000 05D0/

      DATA BC( 0,2) / 0.99363 49867 16925 14D0/
      DATA BC( 1,2) /-0.00646 71526 00616 03D0/
      DATA BC( 2,2) /-0.00010 60188 22351 55D0/
      DATA BC( 3,2) /-0.00000 41406 57716 24D0/
      DATA BC( 4,2) /-0.00000 02916 95418 21D0/
      DATA BC( 5,2) /-0.00000 00365 71574 33D0/
      DATA BC( 6,2) /-0.00000 00075 81590 37D0/
      DATA BC( 7,2) /-0.00000 00019 23008 52D0/
      DATA BC( 8,2) /-0.00000 00004 20438 80D0/
      DATA BC( 9,2) /-0.00000 00000 39372 04D0/
      DATA BC(10,2) / 0.00000 00000 19007 44D0/
      DATA BC(11,2) / 0.00000 00000 10137 64D0/
      DATA BC(12,2) / 0.00000 00000 01331 30D0/
      DATA BC(13,2) /-0.00000 00000 00676 92D0/
      DATA BC(14,2) /-0.00000 00000 00311 72D0/
      DATA BC(15,2) / 0.00000 00000 00011 87D0/
      DATA BC(16,2) / 0.00000 00000 00040 21D0/
      DATA BC(17,2) / 0.00000 00000 00004 78D0/
      DATA BC(18,2) /-0.00000 00000 00004 74D0/
      DATA BC(19,2) /-0.00000 00000 00001 16D0/
      DATA BC(20,2) / 0.00000 00000 00000 59D0/
      DATA BC(21,2) / 0.00000 00000 00000 21D0/
      DATA BC(22,2) /-0.00000 00000 00000 08D0/
      DATA BC(23,2) /-0.00000 00000 00000 03D0/

      DATA CC( 0,1) / 0.99353 64122 76093 39D0/
      DATA CC( 1,1) /-0.00631 44392 60798 63D0/
      DATA CC( 2,1) / 0.00014 30095 80961 13D0/
      DATA CC( 3,1) /-0.00000 57870 60592 03D0/
      DATA CC( 4,1) / 0.00000 03265 50333 20D0/
      DATA CC( 5,1) /-0.00000 00231 23231 95D0/
      DATA CC( 6,1) / 0.00000 00019 39555 14D0/
      DATA CC( 7,1) /-0.00000 00001 85897 89D0/
      DATA CC( 8,1) / 0.00000 00000 19868 42D0/
      DATA CC( 9,1) /-0.00000 00000 02326 79D0/
      DATA CC(10,1) / 0.00000 00000 00294 68D0/
      DATA CC(11,1) /-0.00000 00000 00039 95D0/
      DATA CC(12,1) / 0.00000 00000 00005 75D0/
      DATA CC(13,1) /-0.00000 00000 00000 87D0/
      DATA CC(14,1) / 0.00000 00000 00000 14D0/
      DATA CC(15,1) /-0.00000 00000 00000 02D0/

      DATA CC( 0,2) / 1.00914 95380 72789 40D0/
      DATA CC( 1,2) / 0.00897 12068 42483 60D0/
      DATA CC( 2,2) /-0.00017 13895 98261 54D0/
      DATA CC( 3,2) / 0.00000 65547 92549 82D0/
      DATA CC( 4,2) /-0.00000 03595 19190 49D0/
      DATA CC( 5,2) / 0.00000 00250 24412 19D0/
      DATA CC( 6,2) /-0.00000 00020 74924 13D0/
      DATA CC( 7,2) / 0.00000 00001 97223 67D0/
      DATA CC( 8,2) /-0.00000 00000 20946 47D0/
      DATA CC( 9,2) / 0.00000 00000 02440 93D0/
      DATA CC(10,2) /-0.00000 00000 00307 91D0/
      DATA CC(11,2) / 0.00000 00000 00041 61D0/
      DATA CC(12,2) /-0.00000 00000 00005 97D0/
      DATA CC(13,2) / 0.00000 00000 00000 91D0/
      DATA CC(14,2) /-0.00000 00000 00000 14D0/
      DATA CC(15,2) / 0.00000 00000 00000 02D0/

      LEX=.FALSE.
      GO TO 8





      ENTRY DEBIR3(X,NU)

      LEX=.TRUE.

    8 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .OR. NU .LT. 0 .AND. X .LE. 0
     1   .OR. NU .GT. 0 .AND. X .LT. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEI ,'C340.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEIE,'C340.1',ERRTXT)
      ELSEIF(X .EQ. 0) THEN
       S=0
      ELSEIF(X .LT. 8) THEN
       Y=(HF*X)**2
       XN=PP(NU)
       XL=XN+2
       A0=1
       A1=1+2*Y/((XL+1)*(XN+1))
       A2=1+Y*(4+3*Y/((XL+2)*(XN+2)))/((XL+3)*(XN+1))
       B0=1
       B1=1-Y/(XL+1)
       B2=1-Y*(1-Y/(2*(XL+2)))/(XL+3)
       T1=3+XL
       V1=3-XL
       V3=XL-1
       V2=V3+V3
       C=0
       DO 33 N = 3,30
       C0=C
       T1=T1+2
       T2=T1-1
       T3=T2-1
       T4=T3-1
       T5=T4-1
       T6=T5-1
       V1=V1+1
       V2=V2+1
       V3=V3+1
       U1=N*T4
       E=V3/(U1*T3)
       U2=E*Y
       F1=1+Y*V1/(U1*T1)
       F2=(1+Y*V2/(V3*T2*T5))*U2
       F3=-Y*Y*U2/(T4*T5*T5*T6)
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 34
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   33  CONTINUE
   34  S=GG(NU)*(HF*X)**PP(NU)*C
       IF(LEX) S=EXP(-X)*S
      ELSE
       R=1/X
       W=SQRT(RPI2*R)
       H=16*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 23,0,-1
       B0=BC(I,MU)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=W*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(X)*S
       T=0
       IF(NU .LT. 0) THEN
        H=10*R-1
        ALFA=H+H
        B1=0
        B2=0
        DO 3 I = 15,0,-1
        B0=CC(I,MU)+ALFA*B1-B2
        B2=B1
    3   B1=B0
        R=EXP(-X)
        T=W3*W*R*(B0-H*B2)
        IF(LEX) T=R*T
       END IF
       S=S+T
      END IF
      GO TO 99





      ENTRY DBSKR3(X,NU)

      LEX=.FALSE.
      GO TO 9





      ENTRY DEBKR3(X,NU)

      LEX=.TRUE.

    9 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .OR. X .LE. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEK ,'C340.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEKE,'C340.1',ERRTXT)
      ELSEIF(X .LE. 1) THEN
       A0=PP(-1)
       B=HF*X
       D=-LOG(B)
       F=A0*D
       E=EXP(F)
       G=(GM*A0+GP)*E
       BK=C1*(HF*GM*(E+1/E)+GP*D*SINH(F)/F)
       F=BK
       E=A0**2
       P=HF*C1*G
       Q=HF/G
       C=1
       D=B**2
       BK1=P
       DO 11 N = 1,15
       FN=N
       F=(FN*F+P+Q)/(FN**2-E)
       C=C*D/FN
       P=P/(FN-A0)
       Q=Q/(FN+A0)
       G=C*(P-FN*F)
       H=C*F
       BK=BK+H
       BK1=BK1+G
       IF(H*BK1+ABS(G)*BK .LE. EPS*BK*BK1) GO TO 12
   11  CONTINUE
   12  S=BK
       IF(MU .EQ. 2) S=BK1/B
       IF(LEX) S=EXP(X)*S
      ELSEIF(X .LE. 5) THEN
       XN=4*PP(MU)**2
       A=9-XN
       B=25-XN
       C=768*X**2
       C0=48*X
       A0=1
       A1=(16*X+7+XN)/A
       A2=(C+C0*(XN+23)+XN*(XN+62)+129)/(A*B)
       B0=1
       B1=(16*X+9-XN)/A
       B2=(C+C0*B)/(A*B)+1
       C=0
       DO 24 N = 3,30
       C0=C
       FN=N
       FN2=FN+FN
       FN1=FN2-1
       FN3=FN1/(FN2-3)
       FN4=12*FN**2-(1-XN)
       FN5=16*FN1*X
       RAN=1/((FN2+1)**2-XN)
       F1=FN3*(FN4-20*FN)+FN5
       F2=28*FN-FN4-8+FN5
       F3=FN3*((FN2-5)**2-XN)
       A=(F1*A2+F2*A1+F3*A0)*RAN
       B=(F1*B2+F2*B1+F3*B0)*RAN
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 25
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   24  CONTINUE
   25  S=C/SQRT(RPIH*X)
       IF(.NOT.LEX) S=EXP(-X)*S
      ELSE
       R=1/X
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 13 I = 15,0,-1
       B0=CC(I,MU)+ALFA*B1-B2
       B2=B1
   13  B1=B0
       S=SQRT(PIH*R)*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(-X)*S
      END IF




   99 DBSIR3=S

      RETURN
  101 FORMAT('INCORRECT ARGUMENT OR INDEX, X = ',1P,E15.6,' NU = ',I5)
      END
*CMZ :          21/11/2017  14.35.24  by  Michael Scheer
*-- Author :
# 1 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F"
*
* $Id: bsja64.F,v 1.1.1.1 1996/04/01 15:02:08 mclareni Exp $
*
* $Log: bsja64.F,v $
* Revision 1.1.1.1 1996/04/01 15:02:08 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F" 2




      SUBROUTINE BSJA(X,A,NMAX,ND,B)


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* imp64.inc
*

      IMPLICIT REAL (A-H,O-Z)
# 18 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F" 2
      REAL SX,D,T,Q,U,V,TC(11)
      CHARACTER*80 ERRTXT
      CHARACTER NAMEJ*(*),NAMEI*(*)






      PARAMETER (NAMEJ = 'BSJA', NAMEI = 'BSIA')
      EXTERNAL GAMMA

      LOGICAL LJA,LIA,LEV,LER
      DIMENSION B(0:*),BA(0:100),RR(0:100)

      PARAMETER (Z1 = 1, HF = Z1/2, Z10 = 10)

      DATA TC / 5.7941 E-5,-1.76148E-3, 2.08645E-2,-1.29013E-1,
     1 8.5777 E-1, 1.0125 E+0, 7.75 E-1, 2.3026 E+0,
     2 1.3863 E+0, 7.3576 E-1, 1.3591 E+0/

      LJA=.TRUE.
      LIA=.FALSE.
      SGN=-1
      GO TO 9





      ENTRY BSIA(X,A,NMAX,ND,B)

      LJA=.FALSE.
      LIA=.TRUE.
      SGN=1

    9 LER=.FALSE.
      IF(X .LE. 0) THEN
       WRITE(ERRTXT,101) X
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.1',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.1',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. A .AND. A .LT. 1)) THEN
       WRITE(ERRTXT,102) A
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.2',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.2',ERRTXT)
       LER=.TRUE.
      ELSEIF(ABS(NMAX) .GT. 100) THEN
       WRITE(ERRTXT,103) NMAX
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.3',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.3',ERRTXT)
       LER=.TRUE.
      END IF
      IF(LER) RETURN
      EPS=HF*Z10**(-ND)
      NMX=ABS(NMAX)
      IF(NMAX .LE. 0) NMX=1
      DO 5 N = 0,NMX
      RR(N)=0
    5 BA(N)=0
      D=TC(8)*ND+TC(9)
      SX=X
      Q=0
      IF(NMX .GT. 0) THEN
       V=0.5*D/NMX
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 6 I = 2,6
    6 T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       ENDIF
       Q=NMX*T
      ENDIF




      F=(HF*X)**A/GAMMA(1+A)

      T=1
      V=TC(10)*D/SX
      IF(LIA) THEN
       F=EXP(X)*F
       V=V-TC(10)
      ENDIF
      IF(LJA .OR. LIA .AND. X .LT. D) THEN
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 7 I = 2,6
    7 T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       ENDIF
      ENDIF
      NU=1+MAX(Q,TC(11)*SX*T)

      MU=-1
    2 MU=MU+1
      AL=1
      IF(LJA) THEN
       DO 3 N = 1,NU/2
       XN=N
    3 AL=AL*(XN+A)/(XN+1)
       R=0
       S=0
       LEV=.TRUE.
       DO 4 N = 2*(NU/2),1,-1
       XN=N
       XA=XN+A
       R=1/(2*XA/X-R)
       IF(N .LE. NMX) RR(N-1)=R
       IF(LEV) THEN
        AL=AL*(XN+2)/(XA+A)
        S=R*(AL*XA+S)
       ELSE
        S=R*S
       ENDIF
       LEV=.NOT.LEV
    4 CONTINUE
      ELSE
       DO 23 N = 1,NU
       XN=N
   23 AL=AL*(XN+2*A)/(XN+1)
       R=0
       S=0
       DO 24 N = NU,1,-1
       XN=N
       XA=XN+A
       XA2=XA+XA
       R=1/(XA2/X+R)
       IF(N .LE. NMX) RR(N-1)=R
       AL=AL*(XN+1)/(XA+A)
       S=R*(XA2*AL+S)
   24 CONTINUE
      ENDIF
      B(0)=F/(1+S)
      DO 10 N = 0,NMX-1
   10 B(N+1)=RR(N)*B(N)
      DO 11 N = 0,NMX
      IF(ABS(B(N)-BA(N)) .GT. EPS*ABS(B(N))) THEN
       DO 12 M = 0,NMX
   12 BA(M)=B(M)
       NU=NU+5
       IF(MU .LE. 50) GO TO 2
       WRITE(ERRTXT,104) X,A
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.4',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.4',ERRTXT)
       RETURN
      ENDIF
   11 CONTINUE
      IF(NMAX .LT. 0) THEN
       AL=2/X
       B(1)=AL*A*B(0)+SGN*B(1)
       DO 13 N = 1,-NMAX-1
   13 B(N+1)=AL*(A-N)*B(N)+SGN*B(N-1)
      ENDIF
      RETURN
  101 FORMAT('ILLEGAL ARGUMENT X = ',1P,D15.8)
  102 FORMAT('ILLEGAL ORDER A = ',1P,D15.8)
  103 FORMAT('ILLEGAL NMAX = ',I5)
  104 FORMAT('NO CONVERGENCE FOR X = ',1P,D15.8,' A = ',D15.8,
     1 ' TRY SMALLER ND')
      END
*CMZ :          28/08/2014  12.00.37  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
*
* $Id: r8dp.F,v 1.1.1.1 1996/04/01 15:01:57 mclareni Exp $
*
* $Log: r8dp.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:57  mclareni
* Mathlib gen
*
*

cmsh #include "gen/pilot.h"
#define CERNLIB_DOUBLE

#if defined(CERNLIB_DOUBLE)
      FUNCTION C309R8(Z,ACC)
      COMPLEX*16 C309R8,Z
      DOUBLE PRECISION ACC
      DOUBLE PRECISION X,Y,AX,AY,A

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
      X=DREAL(Z)
      Y=DIMAG(Z)
      AX=ABS(X)
      AY=ABS(Y)
      A=5*ACC*(AX+AY)
      IF(AX .LT. A) X=0
      IF(AY .LT. A) Y=0
      C309R8=DCMPLX(X,Y)
      RETURN
      END
#endif
*CMZ :          28/08/2014  12.25.01  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
# 1 "cdigam64.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "cdigam64.F"
*
* $Id: cdigam64.F,v 1.1.1.1 1996/04/01 15:01:56 mclareni Exp $
*
* $Log: cdigam64.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:56 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "cdigam64.F" 2



cmsh: Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX cdigam64.F


      FUNCTION WDIGAM(Z)
# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
# 18 "cdigam64.F" 2
# 1 "/usr/include/gen/defc64.inc" 1 3 4
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* defc64.inc
*







      COMPLEX*16
# 19 "cdigam64.F" 2
     + WDIGAM

# 1 "/usr/include/gen/defc64.inc" 1 3 4
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* defc64.inc
*







      COMPLEX*16
# 22 "cdigam64.F" 2
     + Z,U,V,H,R,P
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT




      PARAMETER (NAME = 'CDIGAM/WDIGAM')

      DIMENSION C(6)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)

# 1 "/usr/include/gen/gcmpfun.inc" 1 3 4
*
* $Id: gcmpfun.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: gcmpfun.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
*
* gcmpfun.inc
*

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 14 "/usr/include/gen/gcmpfun.inc" 2 3 4
     + GREAL,GIMAG,XARG,YARG
# 1 "/usr/include/gen/defc64.inc" 1 3 4
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* defc64.inc
*







      COMPLEX*16
# 16 "/usr/include/gen/gcmpfun.inc" 2 3 4
     + ZARG,GCONJG,GCMPLX
      GREAL( ZARG)=DREAL( ZARG)
      GIMAG( ZARG)=DIMAG( ZARG)
      GCONJG(ZARG)=DCONJG(ZARG)
      GCMPLX(XARG,YARG)=DCMPLX(XARG,YARG)
# 37 "cdigam64.F" 2
CSEQ,GCMPLX.

      DATA C(1) / 8.33333 33333 33333 33D-2/
      DATA C(2) /-8.33333 33333 33333 33D-3/
      DATA C(3) / 3.96825 39682 53968 25D-3/
      DATA C(4) /-4.16666 66666 66666 67D-3/
      DATA C(5) / 7.57575 75757 57575 76D-3/
      DATA C(6) /-2.10927 96092 79609 28D-2/

      U=Z
      X=U
      A=ABS(X)
      IF(GIMAG(U) .EQ. 0 .AND. -A .EQ. INT(X)) THEN
       H=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C307.1',ERRTXT)
      ELSE
       IF(X .LT. 0) U=-U
       V=U
       H=0
       IF(A .LT. 15) THEN
        H=1/V
        DO 1 I = 1,14-INT(A)
        V=V+1
    1 H=H+1/V
        V=V+1
       END IF
       R=1/V**2
       P=R*C(1)
       DO 2 I = 6,1,-1
    2 P=R*(C(I)+P)
       H=LOG(V)-HF/V-P-H
       IF(X .LT. 0) THEN
        V=PI*U
        X=V
        A=SIN(X)
        X=COS(X)
        Y=TANH(GIMAG(V))
        H=H+1/U+PI*GCMPLX(X,-A*Y)/GCMPLX(A,X*Y)
       END IF
      ENDIF

      WDIGAM=H




      RETURN
  101 FORMAT(1X,'ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END
*CMZ :          03/12/2024  10.40.21  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX cfft.F

# 1 "cfft.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "cfft.F"
*
* $Id: cfft.F,v 1.1.1.1 1996/02/15 17:48:48 mclareni Exp $
*
* $Log: cfft.F,v $
* Revision 1.1.1.1 1996/02/15 17:48:48 mclareni
* Kernlib
*
*
# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 10 "cfft.F" 2
      SUBROUTINE CFFT(A,MSIGN)
cmsh      COMPLEX A(1),U,W,T
      COMPLEX A(2),U,W,T
      IF(MSIGN.EQ.0) RETURN
      M=IABS(MSIGN)
      N=2**M
      NV2=N/2
      NM1=N-1
      J=1
      DO 7 I=1,NM1
        IF(I.GE.J) GO TO 5
        T=A(J)
        A(J)=A(I)
        A(I)=T
 5      K=NV2
 6      IF(K.GE.J) GO TO 7
        J=J-K
        K=K/2
        GO TO 6
 7      J=J+K
        DO 8 I=1,N,2
          T=A(I+1)
          A(I+1)=A(I)-T
 8        A(I )=A(I)+T
          IF(M.EQ.1) RETURN
          C=0.
          S=ISIGN(1,MSIGN)
          LE=2
          DO 20 L=2,M
            W=CMPLX(C,S)
            U=W
            C=SQRT(C*.5+.5)
            S=AIMAG(W)/(C+C)
            LE1=LE
            LE=LE1+LE1
            DO 9 I=1,N,LE
              IP=I+LE1
              T=A(IP)
              A(IP)=A(I)-T
 9          A(I) =A(I)+T
            DO 20 J=2,LE1
              DO 10 I=J,N,LE
                  IP=I+LE1
                  T=A(IP)*U
                  A(IP)=A(I)-T
 10             A(I) =A(I)+T
 20           U=U*W
      RETURN
      END
*CMZ :  1.16/04 16/04/2014  14.22.39  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
*
* $Id: cfstft.F,v 1.2 1997/12/15 16:18:42 mclareni Exp $
*
* $Log: cfstft.F,v $
* Revision 1.2  1997/12/15 16:18:42  mclareni
* Changes for the Portland Group f77 compiler inside cpp define CERNLIB_QFPGF77
*
* Revision 1.1  1996/04/16 15:57:01  mclareni
* The name of cfstft was mistyped and also becomes D706, not 705
*
*
*
      SUBROUTINE CFSTFT(MS,A)

      COMPLEX A(0:*),U,W,T

      IF(MS .EQ. 0) GO TO 3
      M=ABS(MS)
      N=2**M
      J=0
      DO 7 I = 0,N-2
      IF(I .LT. J) THEN
       T=A(J)
       A(J)=A(I)
       A(I)=T
      ENDIF
      K=N/2
    6 IF(K .LE. J) THEN
       J=J-K
       K=K/2
       GO TO 6
      ENDIF
    7 J=J+K
      DO 8 I = 0,N-1,2
      T=A(I+1)
      A(I+1)=A(I)-T
      A(I)=A(I)+T
    8 CONTINUE
      C=0
      S=SIGN(1,MS)
      LE=2
      DO 2 L = 2,M
      W=CMPLX(C,S)
      U=W
      C=SQRT(0.5*C+0.5)
      S=AIMAG(W)/(C+C)
      LE1=LE
      LE=LE1+LE1
      DO 9 I = 0,N-1,LE
      T=A(I+LE1)
      A(I+LE1)=A(I)-T
      A(I)=A(I)+T
    9 CONTINUE
      DO 2 J = 2,LE1
      DO 1 I = J-1,N-1,LE
      T=A(I+LE1)*U
      A(I+LE1)=A(I)-T
      A(I)=A(I)+T
    1 CONTINUE
      U=U*W
    2 CONTINUE
    3 RETURN
      END
*CMZ :          28/08/2014  12.41.21  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
# 1 "clogam64.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "clogam64.F"
*
* $Id: clogam64.F,v 1.1.1.1 1996/04/01 15:01:55 mclareni Exp $
*
* $Log: clogam64.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:55 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "clogam64.F" 2



cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX clogam64.F
cmsh equivalent to : cpp -E -DCERNLIB_DOUBLE -DCERNLIB_LINUX
cms for this routine

      FUNCTION WLGAMA(Z)
# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
# 19 "clogam64.F" 2
# 1 "/usr/include/gen/defc64.inc" 1 3 4
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* defc64.inc
*







      COMPLEX*16
# 20 "clogam64.F" 2
     + WLGAMA
     + ,WLOGAM

# 1 "/usr/include/gen/defc64.inc" 1 3 4
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* defc64.inc
*







      COMPLEX*16
# 24 "clogam64.F" 2
C + Z,W,U,V,H,P,R,GCONJG,GCMPLX
     + Z, U,V,H,P,R
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT




      PARAMETER (NAME = 'CLGAMA/WLGAMA')

      DIMENSION C(10)

      PARAMETER (Z1 = 1, HF = Z1/2)

# 1 "/usr/include/gen/gcmpfun.inc" 1 3 4
*
* $Id: gcmpfun.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: gcmpfun.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
*
* gcmpfun.inc
*

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 14 "/usr/include/gen/gcmpfun.inc" 2 3 4
     + GREAL,GIMAG,XARG,YARG
# 1 "/usr/include/gen/defc64.inc" 1 3 4
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* defc64.inc
*







      COMPLEX*16
# 16 "/usr/include/gen/gcmpfun.inc" 2 3 4
     + ZARG,GCONJG,GCMPLX
      GREAL( ZARG)=DREAL( ZARG)
      GIMAG( ZARG)=DIMAG( ZARG)
      GCONJG(ZARG)=DCONJG(ZARG)
      GCMPLX(XARG,YARG)=DCMPLX(XARG,YARG)
# 39 "clogam64.F" 2

      DATA PI /3.14159 26535 89793 24D+0/
      DATA C1 /9.18938 53320 46727 42D-1/
      DATA C2 /1.14472 98858 49400 17D+0/

      DATA C( 1) / 8.33333 33333 33333 33D-2/
      DATA C( 2) /-2.77777 77777 77777 78D-3/
      DATA C( 3) / 7.93650 79365 07936 51D-4/
      DATA C( 4) /-5.95238 09523 80952 38D-4/
      DATA C( 5) / 8.41750 84175 08417 51D-4/
      DATA C( 6) /-1.91752 69175 26917 53D-3/
      DATA C( 7) / 6.41025 64102 56410 26D-3/
      DATA C( 8) /-2.95506 53594 77124 18D-2/
      DATA C( 9) / 1.79644 37236 88305 73D-1/
      DATA C(10) /-1.39243 22169 05901 12D+0/
C GREAL(U)=DREAL(U)
C GIMAG(U)=DIMAG(U)
C GCONJG(U)=DCONJG(U)
C GCMPLX(X,Y)=DCMPLX(X,Y)





      ENTRY WLOGAM(Z)


      X=Z
      Y=GIMAG(Z)
      IF(Y .EQ. 0 .AND. -ABS(X) .EQ. INT(X)) THEN
       H=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C306.1',ERRTXT)
      ELSE
       YA=ABS(Y)
       U=GCMPLX(X,YA)
       IF(X .LT. 0) U=1-U
       H=0
       UR=U
       IF(UR .LT. 7) THEN
        UI=GIMAG(U)
        A=ATAN2(UI,UR)
        H=U
        DO 1 I = 1,6-INT(UR)
        UR=UR+1
        U=GCMPLX(UR,UI)
        H=H*U
    1 A=A+ATAN2(UI,UR)
        H=GCMPLX(HF*LOG(GREAL(H)**2+GIMAG(H)**2),A)
        U=U+1
       ENDIF
       R=1/U**2
       P=R*C(10)
       DO 2 I = 9,2,-1
    2 P=R*(C(I)+P)
       H=C1+(U-HF)*LOG(U)-U+(C(1)+P)/U-H
       IF(X .LT. 0) THEN
        UR=INT(X)-1
        UI=PI*(X-UR)
        X=PI*YA
        T=EXP(-X-X)
        A=SIN(UI)
        T=X+HF*LOG(T*A**2+(HF*(1-T))**2)
        A=ATAN2(COS(UI)*TANH(X),A)-UR*PI
        H=C2-GCMPLX(T,A)-H
       ENDIF
       IF(Y .LT. 0) H=GCONJG(H)
      ENDIF

      WLGAMA=H




      RETURN
  101 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END
*CMZ :          02/05/2017  14.44.57  by  Michael Scheer
*-- Author :

*KEEP,CMSH.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 17 "cvcpy.F" 2

# 1 "d501l1.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501l1.F"
*
* $Id: d501l1.F,v 1.2 2003/09/02 12:41:10 mclareni Exp $
*
* $Log: d501l1.F,v $
* Revision 1.2  2003/09/02 12:41:10  mclareni
* Version corrected by D.A. and C.H. (Aug 2003).
* After column pivoting the components of the covariance matrix were not
* restored in the correct order by using the JPVT vector. This resulted in a
* quasi-random reshuffling of the errors in output.
*
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
* Version corrected by D.A. and C.H. (Aug 2003)
* After coloumn pivoting the components of the covariance
* matrix were not restored in the correct order by using the
* JPVT vector. This resulted in a quasi-random reshuffling of
* the errors in output
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 22 "d501l1.F" 2
      SUBROUTINE D501L1(VERS,SUB,K,N,X,NX,Y,SY,MODE,EPS,MAXIT,
     1                  IPRT,M,A,AL,AU,PHI1,DPHI,IAFR,MFR,
     2                  COV,NC,STD,P,LAMU,DSCAL,W1,W2,W3,TAU,
     3                  COPYF,COPYDF,R2,R1,F,DF,JPVT,NERROR)


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 28 "d501l1.F" 2

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 29 "d501l1.F" 2
     +   JP2,LAMBDA,LAMU,LK,MY
      CHARACTER VERS*6
      LOGICAL LFN,LID,LRP,LPR
      DIMENSION X(*),Y(*),SY(*),A(*),AL(*),AU(*),DPHI(*),F(*),DF(N,*)
      DIMENSION STD(*),P(*),LAMU(*),DSCAL(*),W1(*),W2(*),W3(*),COV(NC,*)
      DIMENSION IAFR(*),TAU(*),JPVT(*),COPYF(*),COPYDF(N+M,*)
      DIMENSION R1(M,*),R2(M,*)
      DIMENSION W64(64)

      PARAMETER (Z0 = 0, Z1 = 1, HALF = Z1/2, R3 = Z1/3, R10 = Z1/10)
      PARAMETER (SIG1 = R10, SIG2 = 11*R10, COEF = R10**3, STEP = Z1)
      PARAMETER (RHO1 = R10**4, RHO2 = Z1/4, RHO3 = 3*Z1/4)

      EXTERNAL SUB

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE IS ONE OF THE MAIN ROUTINE OF THE LEAMAX PACKAGE.
*   IT SOLVES TWO DIFFERENT PROBLEMS DEPENDING ON THE VALUE
*   OF THE PARAMETER VERS.
*   ( VERS = DSUMSQ : GENERAL NONLINEAR LEAST SQUARES PROBLEM
*     VERS = DFUNFT : LEAST SQUARES DATA FITTING PROBLEM      )
*
*   IN ALL CASES BOUNDS ON THE VARIABLES MAY BE SET.
*
************************************************************************

************************************************************************
*   COMPUTE AN APPROXIMATION  EPS0  TO THE RELATIVE MACHINE PRECISION
************************************************************************

      EPS0=Z1
    5 EPS0=EPS0/10
      IF (Z1+EPS0 .NE. Z1) GO TO 5
      EPS0=10*EPS0

************************************************************************
*   CHECK THE VAUES OF INPUT PARAMETERS
************************************************************************

      NERROR=0

cmsh      CALL D501P1(K,N,NC,X,NX,Y,SY,MODE,EPS0,EPS,MAXIT,IPRT,M,A,AL,AU,
      CALL D501P1(K,N,NC,X(1),NX,Y(1),SY,MODE,EPS0,EPS,MAXIT,IPRT,M,A,AL,AU,
     +            NERROR,VERS)

      IF (NERROR .NE. 0) RETURN

************************************************************************
*   SET INITIAL VALUES
************************************************************************

      EPS1=10*EPS0

      LFN=.FALSE.
      LID=.FALSE.
      LRP=.FALSE.
      LPR=IPRT .NE. 0

      ITER=0

      CALL DVSET(M,Z1,DSCAL(1),DSCAL(2))
      CALL DVSET(M,Z0,LAMU(1),LAMU(2))
      CALL DVSET(M,Z0,STD(1),STD(2))

************************************************************************
*   COMPUTE INITIAL VALUE  PHI1  OF OBJECTIVE FUNCTION
************************************************************************

       CALL D501SF(VERS,SUB,0,M,A,N,F,DF,K,NX,X,Y,SY,W2,NERROR)
       IF (NERROR .NE. 0) RETURN
       PHI1=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

************************************************************************
*   COMPUTE  F, DF, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,W2,DPHI,DSCAL,LAMU,F,DF,IAFR,
     +            MFR,SUB,EPS0,EPS1,MODE,VERS,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF MFR = 0 MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

       IF(MFR .EQ. 0) GO TO 230

************************************************************************
*   ITERATION BEGINS
************************************************************************

      DELTA=0
      LAMBDA=0

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=0
      DO 10 I=1,MFR
   10 DPHINO=DPHINO+DPHI(IAFR(I))**2
      DPHINO=SQRT(DPHINO)

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,VERS)

      DA=0
      DO 20 I=1,MFR
   20 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)

************************************************************************
*   ITERATION WITH GAUSS-NEWTON STEP
************************************************************************

   30 LAMBDA=0

************************************************************************
*   COPY  F  AND  DF
*   COMPUTE THE QR FACTORIZATION WITH COLUMN PIVOTING OF  DF
*   AND SOLVE THE LINEAR LEAST SQUARES PROBLEM USING
*   LAPACK ROUTINES  DGEQPF , DORMQR , DTRTRS
************************************************************************

      CALL DVSCL(N,-Z1,F(1),F(2),COPYF(1),COPYF(2))
      CALL DMCPY(N,MFR,DF(1,1),DF(1,2),DF(2,1),COPYDF(1,1),
     1           COPYDF(1,2),COPYDF(2,1))
C**** KSK 25.07.95
      DO 31 NN = 1,MFR
   31 JPVT(NN)=0
C     CALL DVSET(MFR,Z0,JPVT(1),JPVT(2))
C**** KSK 25.07.95

      CALL DGEQPF(N,MFR,COPYDF,N+M,JPVT,TAU,W3,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DORMQR('L','T',N,1,MFR,COPYDF,N+M,TAU,COPYF,N,W64,64,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DTRTRS('U','N','N',MFR,1,COPYDF,N+M,COPYF,N,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      DO 40 I=1,MFR
   40 P(JPVT(I))=COPYF(I)

************************************************************************
*   COMPUTE THE MATRIX  R2
************************************************************************

      DO  50 I=1,MFR
      DO  50 J=1,MFR
      IF (I .GT. J) THEN
       R1(I,J)=0
      ELSE
       R1(I,J)=COPYDF(I,J)
      ENDIF
   50 CONTINUE

      DO  60 I=1,MFR
      DO  60 J=1,MFR
C   60 R2(I,J)=R1(I,JPVT(J))
C D.A. & C.H. Aug 2003
   60 R2(I,JPVT(J))=R1(I,J)

      CALL DMMLT(MFR,MFR,MFR,R2(1,1),R2(2,1),R2(1,2),R2(1,1),R2(1,2),
     +           R2(2,1),COV(1,1),COV(1,2),COV(2,1),W2)

************************************************************************
*   COMPUTE THE L2-NORM OF THE SCALED VECTOR  P
************************************************************************

      DP=0
      DO 70 I=1,MFR
   70 DP=DP+(DSCAL(I)*P(I))**2
      DP=SQRT(DP)

************************************************************************
*   COMPUTE THE STEP SIZE  ALFA
************************************************************************

      ALFA=1
      DO 80 I=1,MFR
      IF(P(I) .NE. 0) THEN
       IF(P(I) .GT. 0) THEN
        ALFA1=AU(IAFR(I))
       ELSE
        ALFA1=AL(IAFR(I))
       ENDIF
       ALFA1=(ALFA1-A(IAFR(I)))/P(I)
       IF(ALFA1 .EQ. 0) THEN
        P(I)=0
       ELSE
        ALFA=MIN(ALFA,ALFA1)
       ENDIF
      ENDIF
   80 CONTINUE

************************************************************************
*   COMPUTE INITIAL DELTA IF NECESSARY
************************************************************************

      IF(.NOT.LID) THEN
       DELTA=STEP*MAX(DA,DP/SIG2)
       LID=.TRUE.
      ENDIF
      IF(DELTA .LE. EPS*DA) GO TO 230

************************************************************************
*   CONTINUATION WITH GAUSS-NEWTON OR SWITCHING TO LEVENBERG-MARQUARDT?
************************************************************************

      IF(DP .GT. SIG2*DELTA) THEN

***********************************************************************
*   DO THE LEVENBERG - MARQUARDT STEP, (HEBDEN'S METHOD).
*   - COMPUTE THE LM - PARAMETER LAMBDA
*   - COMPUTE THE CORRESPONDING STEP P, ITS DP AND ALFA.
***********************************************************************

       DO 90 I=1,MFR
   90  STD(I)=-DPHI(IAFR(I))

       UK=0
       DO 100 I=1,MFR
  100  UK=UK+(DPHI(IAFR(I))/DSCAL(I))**2
       UK=SQRT(UK)/DELTA

************************************************************************
*   COMPUTE INITIAL LAMBDA
************************************************************************

       LAMBDA=COEF*UK

       LK=0
       ITERA=0

  110  ITERA=ITERA+1
       IF(ITERA .GE. 50) GO TO 230

************************************************************************
*   RESET LAMBDA IF NECESSARY
************************************************************************

       IF(LK .GE. LAMBDA .OR. LAMBDA .GE. UK)
     +    LAMBDA=MAX(COEF*UK,SQRT(LK*UK))

************************************************************************
*   COMPUTE NEW P FOR NEW LAMBDA
************************************************************************

************************************************************************
*   COPY F AND DF, AND EXTEND  DF  BY  SQRT(LAMBDA) * DIAG(DSCAL(I))
*   COMPUTE THE QR FACTORIZATION WITH COLUMN PIVOTING OF  EXTENDED  DF
*   AND SOLVE THE LINEAR LEAST SQUARES PROBLEM USING  LAPACK  ROUTINES
*   DGEQPF , DORMQR , DTRTRS
************************************************************************

      CALL DVSET(N+MFR,Z0,COPYF(1),COPYF(2))
      CALL DVSCL(N,-Z1,F(1),F(2),COPYF(1),COPYF(2))
      CALL DMSET(N+MFR,MFR,Z0,COPYDF(1,1),COPYDF(1,2),COPYDF(2,1))
      CALL DMCPY(N,MFR,DF(1,1),DF(1,2),DF(2,1),
     +           COPYDF(1,1),COPYDF(1,2),COPYDF(2,1))
      DO 120 I=1,MFR
  120 COPYDF(N+I,I)=SQRT(LAMBDA)*DSCAL(I)
C**** KSK 25.07.95
      DO 121 NN = 1,MFR
  121 JPVT(NN)=0
C     CALL DVSET(MFR,Z0,JPVT(1),JPVT(2))
C**** KSK 25.07.95

      CALL DGEQPF(N+MFR,MFR,COPYDF,N+M,JPVT,TAU,W3,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DORMQR('L','T',N+MFR,1,MFR,COPYDF,N+M,TAU,COPYF,
     +            N+MFR,W64,64,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DTRTRS('U','N','N',MFR,1,COPYDF,N+M,COPYF,N+MFR,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      DO 130 I=1,MFR
  130 P(JPVT(I))=COPYF(I)


************************************************************************
*   STOP ITERATION?
************************************************************************

       DP=0
       DO 140 I=1,MFR
  140  DP=DP+(DSCAL(I)*P(I))**2
       DP=SQRT(DP)

       IF(SIG1*DELTA .GT. DP .OR. DP .GT. SIG2*DELTA) THEN

************************************************************************
*   CONTINUE ITERATION FOR LAMBDA
************************************************************************

        IF(DP .LE. 0) GO TO 230
        P1=DP-DELTA
        DO 150 I=1,MFR
  150   W1(I)=DSCAL(I)**2*P(I)

************************************************************************
*   COMPUTE THE MATRIX  R1
************************************************************************

      DO 160 I=1,MFR
      DO 160 J=1,MFR
      IF (I .GT. J) THEN
       R1(I,J)=0
      ELSE
       R1(I,J)=COPYDF(I,J)
      ENDIF
  160 CONTINUE

      DO 170 I=1,MFR
      DO 170 J=1,MFR
C  170 R2(I,J)=R1(I,JPVT(J))
C D.A. & C.H. Aug 2003
  170 R2(I,JPVT(J))=R1(I,J)

      CALL DMMLT(MFR,MFR,MFR,R2(1,1),R2(2,1),R2(1,2),R2(1,1),R2(1,2),
     +           R2(2,1),R1(1,1),R1(1,2),R1(2,1),W2)

      CALL DSINV(MFR,R1,M,NERROR)
      IF (NERROR .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      P1P=-DMBIL(MFR,W1(1),W1(2),R1(1,1),R1(1,2),R1(2,1),
     +               W1(1),W1(2))/DP

************************************************************************
*   UPDATE LK, UK, LAMBDA
************************************************************************

        IF(P1 .LT. 0) UK=LAMBDA
        LK=MAX(LK,LAMBDA-P1/P1P)
        IF(LK .GE. UK) UK=2*LK
        LAMBDA=LAMBDA-(DP/DELTA)*(P1/P1P)
        GO TO 110
       ENDIF
      ENDIF

************************************************************************
*   END OF LEVENBERG - MARQUARDT STEP
************************************************************************

      ALFA=1
      DO 180 I=1,MFR
      IF(P(I) .NE. 0) THEN
       IF(P(I) .GT. 0) THEN
        ALFA1=AU(IAFR(I))
       ELSE
        ALFA1=AL(IAFR(I))
       ENDIF
       ALFA1=(ALFA1-A(IAFR(I)))/P(I)
       IF(ALFA1 .EQ. 0) THEN
        P(I)=0
       ELSE
        ALFA=MIN(ALFA,ALFA1)
       ENDIF
      ENDIF
  180 CONTINUE

************************************************************************
*   COMPUTE   A + ALPHA * P
************************************************************************

      CALL DVCPY(M,A(1),A(2),W1(1),W1(2))
      DO 190 I=1,MFR
  190 W1(IAFR(I))=A(IAFR(I))+ALFA*P(I)

************************************************************************
*   COMPUTE VALUE  PHI2  OF THE OBJECTIVE FUNCTION
************************************************************************

      CALL D501SF(VERS,SUB,0,M,W1,N,F,DF,K,NX,X,Y,SY,W2,NERROR)
      IF (NERROR .NE. 0) RETURN

      PHI2=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

      PHMAXI=1
      IF(PHI1 .GT. 0) PHMAXI=1/SQRT(PHI1)
      CALL DVSCL(MFR,PHMAXI,P(1),P(2),W2(1),W2(2))
      JP2=DMBIL(MFR,W2(1),W2(2),COV(1,1),COV(1,2),COV(2,1),W2(1),W2(2))

************************************************************************
*   COMPUTE THE APPROXIMATION MEASURE  RHO  AND THE UPDATING FACTOR  MY
*   FOR DELTA
************************************************************************

       IF(PHI1 .LE. 0) THEN
        RHO=1
        MY=HALF
       ELSE
        S2=LAMBDA*DP**2/PHI1
        S3=1-PHI2/PHI1
        S4=HALF*JP2+S2
        IF(S4 .EQ. 0) THEN
         RHO=1
         MY=R10
        ELSE
         RHO=0
         IF(S3 .GT. 0) RHO=S3/S4
         MY=-HALF*(JP2+S2)
         S2=2*MY+S3
         IF(S2 .EQ. 0) THEN
          MY=R10
         ELSEIF(S3 .EQ. 0) THEN
          MY=HALF
         ELSE
          MY=MIN(MAX(MY/S2,R10),HALF)
         ENDIF
        ENDIF
       ENDIF

************************************************************************
*   END OF COMPUTATTION OF RHO AND MY
************************************************************************

************************************************************************
*   IF RHO .LE. RHO1, REDUCE DELTA BY FACTOR MY AND MAKE NEW LEVENBERG-
*   MARQUARDT STEP, OTHERWISE ACCEPT P
************************************************************************

      IF(RHO .LE. RHO1) THEN
       DELTA=MY*DELTA
       DA=0
       DO 200 I=1,MFR
  200  DA=DA+(DSCAL(I)*A(IAFR(I)))**2
       DA=SQRT(DA)
       GO TO 30
      ENDIF
      CALL DVCPY(M,W1(1),W1(2),A(1),A(2))
      DA=0
      DO 210 I=1,MFR
  210 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)
      MFROLD=MFR

************************************************************************
*   COMPUTE  F, DF, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,W2,DPHI,DSCAL,LAMU,F,DF,IAFR,
     +            MFR,SUB,EPS0,EPS1,MODE,VERS,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF MFR = 0  MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

      IF(MFR .EQ. 0) THEN
       ITER=ITER+1
       GO TO 230
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=0
      DO 220 I=1,MFR
  220 DPHINO=DPHINO+DPHI(IAFR(I))**2
      DPHINO=SQRT(DPHINO)

************************************************************************
*   TERMINATION CRITERION
************************************************************************

      IF (     PHI2      .LE. PHI1
     1   .AND. PHI1-PHI2 .LE. EPS*(1+ABS(PHI2))
     2   .AND. DP        .LE. SQRT(EPS)*(1+DA)
     3   .AND. DPHINO    .LE. EPS**R3*(1+ABS(PHI2))) LFN=.TRUE.

      ITER=ITER+1
      PHI1=PHI2

      IF(.NOT.LFN) THEN
       CALL D501SF(VERS,SUB,0,M,A,N,F,DF,K,NX,X,Y,SY,W2,NERROR)
       IF (NERROR .NE. 0) RETURN
       PHI1=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

       IF(LPR) THEN
          IF((MOD(ITER,IPRT) .EQ. 0  .OR.  ITER .GE. MAXIT))
     1       CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     2                   MODE,VERS)
       ENDIF

       IF(ITER .GE. MAXIT) THEN
        NERROR=2
        GO TO 230
       ENDIF

************************************************************************
*   UPDATE DELTA AND GO BACK TO GAUSS-NEWTON STEP
************************************************************************

       IF(MFROLD .NE. MFR) LID=.FALSE.
       IF(RHO .LE. RHO2) THEN
        DELTA=MY*DELTA
       ELSE IF(RHO .GE. RHO3 .OR. LAMBDA .EQ. 0) THEN
        DELTA=2*DP
       ENDIF
       GO TO 30
      ENDIF

************************************************************************
*   END OF ITERATION
************************************************************************

  230 LFN=.TRUE.

************************************************************************
*   COMPUTE  F, DF, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,W2,DPHI,DSCAL,LAMU,F,DF,IAFR,
     +            MFR,SUB,EPS0,EPS1,MODE,VERS,MERROR)
      IF(MERROR .NE. 0) THEN
       NERROR=MERROR
       RETURN
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=0
      DO 240 I=1,MFR
  240 DPHINO=DPHINO+DPHI(IAFR(I))**2
      DPHINO=SQRT(DPHINO)

************************************************************************
*   COMPUTE THE VALUE PHI1 OF THE OBJECTIVE FUNCTION
************************************************************************

      PHI1=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

************************************************************************
*   COMPUTE THE COVARIANCE MATRIX  COV  AND THE STANDARD DEVIATION  STD
*   FOR THE FREE VARIABLES
************************************************************************

      CALL DVSET(M,Z0,STD(1),STD(2))

      IF(MFR .GT. 0) THEN

       CALL DSINV(MFR,COV,NC,MERROR)
       IF(MERROR .NE. 0) THEN
        NERROR=4
        RETURN
       ENDIF

       S=2*PHI1
       IF(N .NE. MFR) S=S/(N-MFR)
       CALL DMSCL(MFR,MFR,S,COV(1,1),COV(1,2),COV(2,1),
     +                      COV(1,1),COV(1,2),COV(2,1))
       DO 250 I=1,MFR
  250  STD(IAFR(I))=SQRT(COV(I,I))
      ENDIF

************************************************************************
*   PRINT LAST ITERATION RESULTS
************************************************************************

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,VERS)

      RETURN

      END
*CMZ :          02/05/2017  14.44.29  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

# 1 "d501l2.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501l2.F"
*
* $Id: d501l2.F,v 1.1.1.1 1996/04/01 15:02:19 mclareni Exp $
*
* $Log: d501l2.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "d501l2.F" 2
      SUBROUTINE D501L2(K,N,X,NX,MODE,EPS,MAXIT,IPRT,M,A,AL,AU,
     1                  PHI1,DPHI,IAFR,MFR,STD,P,LAMU,DSCAL,B,W1,W2,
     2                  AM,COV,SUB,NERROR)


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 15 "d501l2.F" 2

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 16 "d501l2.F" 2
     +    JP2,LAMBDA,LAMU,LK,MY
      LOGICAL LFN,LID,LRP,LPR
      DIMENSION X(*),A(*),AL(*),AU(*),DPHI(*),STD(*),P(*),LAMU(*)
      DIMENSION DSCAL(*),B(*),W1(*),W2(*),AM(M,*),COV(M,*),IAFR(*)
      PARAMETER (Z0 = 0, Z1 = 1, HALF = Z1/2, R3 = Z1/3, R10 = Z1/10)
      PARAMETER (SIG1 = R10, SIG2 = 11*R10, COEF = R10**3, STEP = Z1)
      PARAMETER (RHO1 = R10**4, RHO2 = Z1/4, RHO3 = 3*Z1/4)
      dimension y(1),sy(1) !cmsh
      EXTERNAL SUB

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE IS ONE OF THE MAIN ROUTINES OF THE  LEAMAX  PACKAGE
*   FOR SOLVING THE PROBLEM OF MAXIMUM LIKELIHOOD ESTIMATION.
*   BOUNDS ON THE VARIABLES MAY BE SET.
*
************************************************************************

************************************************************************
*   COMPUTE AN APPROXIMATION  EPS0  TO THE RELATIVE MACHINE PRECISION
************************************************************************

      EPS0=Z1
   10 EPS0=EPS0/10
      IF (Z1+EPS0 .NE. Z1) GO TO 10
      EPS0=10*EPS0

************************************************************************
*   CHECK THE VALUES OF INPUT PARAMETERS
************************************************************************

      NERROR=0

      CALL D501P1(K,N,M,X,NX,Y,SY,MODE,EPS0,EPS,MAXIT,IPRT,M,A,AL,AU,
     +            NERROR,'DMAXLK')

      IF(NERROR .NE. 0) RETURN

************************************************************************
*   SET INITIAL VALUES
************************************************************************

      EPS1=10*EPS0

      LFN=.FALSE.
      LID=.FALSE.
      LRP=.FALSE.
      LPR=IPRT .NE. 0

      ITER=0

      CALL DVSET(M,Z1,DSCAL(1),DSCAL(2))
      CALL DVSET(M,Z0,LAMU(1),LAMU(2))
      CALL DVSET(M,Z0,STD(1),STD(2))

************************************************************************
*   COMPUTE INITIAL VALUE  PHI1  OF OBJECTIVE FUNCTION
************************************************************************

      PHI1=0
      IX=1
      DO 20 I=1,N
      CALL SUB(K,X(IX),M,A,F0,ZZ,0,NERROR)
      IF(F0 .LE. 0  .OR.  NERROR .NE. 0) THEN
       NERROR=3
       RETURN
      ENDIF
      PHI1=PHI1-LOG(F0)
   20 IX=IX+NX

************************************************************************
*   COMPUTE J(TRANS) X J, B, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N2(K,N,M,A,AL,AU,X,NX,W1,B,DPHI,DSCAL,LAMU,AM,COV,
     +            IAFR,MFR,SUB,EPS0,EPS1,MODE,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF  MFR = 0  MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

      IF(MFR .EQ. 0) GO TO 190

************************************************************************
*   ITERATION BEGINS
************************************************************************

      DELTA=0
      LAMBDA=0

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=SQRT(DVMPY(MFR,B(1),B(2),B(1),B(2)))

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,'DMAXLK')

      DA=0
      DO 40 I=1,MFR
   40 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)

************************************************************************
*   ITERATION WITH GAUSS-NEWTON STEP
************************************************************************

   50 LAMBDA=0

************************************************************************
*   SOLVE NORMAL EQUATIONS
************************************************************************

      CALL DVSCL(MFR,-Z1,B(1),B(2),W1(1),W1(2))
      CALL DSINV(MFR,AM,M,NERROR)

      IF(NERROR .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF
      CALL DMMPY(MFR,MFR,AM(1,1),AM(1,2),AM(2,1),W1(1),W1(2),P(1),P(2))


************************************************************************
*   COMPUTE THE L2-NORM OF THE SCALED VECTOR  P
************************************************************************

      DP=0
      DO 60 I=1,MFR
   60 DP=DP+(DSCAL(I)*P(I))**2
      DP=SQRT(DP)

************************************************************************
*   COMPUTE THE STEP SIZE  ALFA
************************************************************************

      ALFA=1
      DO 70 I=1,MFR
      IF(P(I) .NE. 0) THEN
       IF(P(I) .GT. 0) THEN
        ALFA1=AU(IAFR(I))
       ELSE
        ALFA1=AL(IAFR(I))
       ENDIF
       ALFA1=(ALFA1-A(IAFR(I)))/P(I)
       IF(ALFA1 .EQ. 0) THEN
        P(I)=0
       ELSE
        ALFA=MIN(ALFA,ALFA1)
       ENDIF
      ENDIF
   70 CONTINUE

************************************************************************
*   COMPUTE INITIAL DELTA IF NECESSARY
************************************************************************

      IF(.NOT.LID) THEN
       DELTA=STEP*MAX(DA,DP/SIG2)
       LID=.TRUE.
      ENDIF
      IF(DELTA .LE. EPS*DA) GO TO 190

************************************************************************
*   CONTINUATION WITH GAUSS-NEWTON OR SWITCHING TO LEVENBERG-MARQUARDT?
************************************************************************

      IF(DP .GT. SIG2*DELTA) THEN

************************************************************************
*   DO THE LEVENBERG - MARQUARDT STEP, (HEBDEN'S METHOD).
*   - COMPUTE THE LM - PARAMETER LAMBDA
*   - COMPUTE THE CORRESPONDING STEP P, ITS DP AND ALFA.
************************************************************************

       CALL DVSCL(MFR,-Z1,B(1),B(2),STD(1),STD(2))
       UK=0
       DO 80 I=1,MFR
   80  UK=UK+(DPHI(IAFR(I))/DSCAL(I))**2
       UK=SQRT(UK)/DELTA

************************************************************************
*   COMPUTE INITIAL LAMBDA
************************************************************************

       LAMBDA=COEF*UK

       LK=0
       ITERA=0

   90  ITERA=ITERA+1
       IF(ITERA .GE. 50) GO TO 190

************************************************************************
*   RESET LAMBDA IF NECESSARY
************************************************************************

       IF(LK .GE. LAMBDA .OR. LAMBDA .GE. UK)
     +    LAMBDA=MAX(COEF*UK,SQRT(LK*UK))

************************************************************************
*   COMPUTE NEW P FOR NEW LAMBDA
************************************************************************

       CALL DMCPY(MFR,MFR,COV(1,1),COV(1,2),COV(2,1),
     +                    AM(1,1),AM(1,2),AM(2,1))
       DO 100 I=1,MFR
  100  AM(I,I)=AM(I,I)+LAMBDA*DSCAL(I)**2

************************************************************************
*   SOLVE NORMAL EQUATIONS
************************************************************************

       CALL DSINV(MFR,AM,M,NERROR)
       IF(NERROR .NE. 0) THEN
        NERROR=4
        RETURN
       ENDIF
       CALL DMMPY(MFR,MFR,AM(1,1),AM(1,2),AM(2,1),STD(1),STD(2),
     +                    P(1),P(2))

************************************************************************
*   COMPUTE THE L2-NORM OF THE SCALED VECTOR  P
************************************************************************

       DP=0
       DO 110 I=1,MFR
  110  DP=DP+(DSCAL(I)*P(I))**2
       DP=SQRT(DP)

************************************************************************
*   STOP ITERATION IN THE CASE OF NORM EQUAL TO ZERO
************************************************************************

       IF (DP .LE. 0) GO TO 190

       IF(SIG1*DELTA .GT. DP .OR. DP .GT. SIG2*DELTA) THEN

************************************************************************
*   CONTINUE ITERATION FOR LAMBDA
************************************************************************

        P1=DP-DELTA
        DO 120 I=1,MFR
  120   W1(I)=DSCAL(I)**2*P(I)
        P1P=-DMBIL(MFR,W1(1),W1(2),AM(1,1),AM(1,2),AM(2,1),
     +                 W1(1),W1(2))/DP

************************************************************************
*   UPDATE LK, UK, LAMBDA
************************************************************************

        IF(P1 .LT. 0) UK=LAMBDA
        LK=MAX(LK,LAMBDA-P1/P1P)
        IF(LK .GE. UK) UK=2*LK
        LAMBDA=LAMBDA-(DP/DELTA)*(P1/P1P)
        GO TO 90
       ENDIF
      ENDIF

************************************************************************
*   END OF LEVENBERG - MARQUARDT STEP
************************************************************************

      ALFA=1
      DO 130 I=1,MFR
      IF(P(I) .NE. 0) THEN
       IF(P(I) .GT. 0) THEN
        ALFA1=AU(IAFR(I))
       ELSE
        ALFA1=AL(IAFR(I))
       ENDIF
       ALFA1=(ALFA1-A(IAFR(I)))/P(I)
       IF(ALFA1 .EQ. 0) THEN
        P(I)=0
       ELSE
        ALFA=MIN(ALFA,ALFA1)
       ENDIF
      ENDIF
  130 CONTINUE

************************************************************************
*   COMPUTE  A + ALPHA * P
************************************************************************
      CALL DVCPY(M,A(1),A(2),W1(1),W1(2))
      DO 140 I=1,MFR
  140 W1(IAFR(I))=A(IAFR(I))+ALFA*P(I)

************************************************************************
*   COMPUTE VALUE  PHI2  OF OBJECTIVE FUNCTION
************************************************************************

      PHI2=0
      IX=1
      DO 150 I=1,N
      CALL SUB(K,X(IX),M,W1,F0,ZZ,0,NERROR)
      IF(F0 .LE. 0  .OR.  NERROR .NE. 0) THEN
       NERROR=3
       RETURN
      ENDIF
      PHI2=PHI2-LOG(F0)
  150 IX=IX+NX

      PHMAXI=1
      AAU=MAX(ABS(PHI1),ABS(PHI2))
      IF(AAU .GT. 0) PHMAXI=1/AAU

      CALL DVSCL(MFR,PHMAXI,P(1),P(2),W2(1),W2(2))
      JP2=DMBIL(MFR,W2(1),W2(2),COV(1,1),COV(1,2),COV(2,1),W2(1),W2(2))

************************************************************************
*   COMPUTE THE APPROXIMATION MEASURE  RHO  AND THE UPDATING FACTOR  MY
*   FOR  DELTA
************************************************************************

       IF(PHI1 .LE. PHI2) THEN
        RHO=0
        MY=R10
       ELSE
        S2=2*(PHI1-PHI2)*PHMAXI**2
        S3=LAMBDA*(DP*PHMAXI)**2
        RHO=S2/(JP2+2*S3)
        MY=-JP2-S3
        S2=S2+2*MY
        IF(S2 .EQ. 0) THEN
         MY=R10
        ELSE
         MY=MIN(MAX(MY/S2,R10),HALF)
        ENDIF
       ENDIF

************************************************************************
*   END OF COMPUTATTION OF RHO AND MY
************************************************************************

************************************************************************
*   IF RHO .LE. RHO1, REDUCE DELTA BY FACTOR MY AND MAKE NEW LEVENBERG-
*   MARQUARDT STEP, OTHERWISE ACCEPT P
************************************************************************

      IF(RHO .LE. RHO1) THEN
       DELTA=MY*DELTA
       DA=0
       DO 160 I=1,MFR
  160  DA=DA+(DSCAL(I)*A(IAFR(I)))**2
       DA=SQRT(DA)
       GO TO 50
      ENDIF
      CALL DVCPY(M,W1(1),W1(2),A(1),A(2))
      DA=0
      DO 170 I=1,MFR
  170 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)
      MFROLD=MFR

************************************************************************
*   COMPUTE J(TRANS) X J, B, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N2(K,N,M,A,AL,AU,X,NX,W1,B,DPHI,DSCAL,LAMU,AM,COV,
     +            IAFR,MFR,SUB,EPS0,EPS1,MODE,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF  MFR = 0  MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

      IF(MFR .EQ. 0) THEN
       ITER=ITER+1
       GO TO 190
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=SQRT(DVMPY(MFR,B(1),B(2),B(1),B(2)))

************************************************************************
*   TERMINATION CRITERION
************************************************************************

      IF (     PHI2      .LE. PHI1
     1   .AND. PHI1-PHI2 .LE. EPS*(1+ABS(PHI2))
     2   .AND. DP        .LE. SQRT(EPS)*(1+DA)
     3   .AND. DPHINO    .LE. EPS**R3*(1+ABS(PHI2)))    LFN=.TRUE.

      ITER=ITER+1
      PHI1=PHI2

      IF(.NOT.LFN) THEN
       IF(ITER .GE. MAXIT) THEN
        IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                      MODE,'DMAXLK')
        NERROR=2
        GO TO 190
       ENDIF

       IF(LPR) THEN
          IF(MOD(ITER,IPRT) .EQ. 0)
     1       CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     2                   MODE,'DMAXLK')
       ENDIF

************************************************************************
*   UPDATE DELTA AND GO BACK TO GAUSS-NEWTON STEP
************************************************************************

       IF(MFROLD .NE. MFR) LID=.FALSE.
       IF(RHO .LE. RHO2) THEN
        DELTA=MY*DELTA
       ELSE IF(RHO .GE. RHO3 .OR. LAMBDA .EQ. 0) THEN
        DELTA=2*DP
       ENDIF

       GO TO 50

      ENDIF

************************************************************************
*   END OF ITERATION
************************************************************************

  190 LFN=.TRUE.

************************************************************************
*   COMPUTE J(TRANS) X J, B, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N2(K,N,M,A,AL,AU,X,NX,W1,B,DPHI,DSCAL,LAMU,AM,COV,
     +            IAFR,MFR,SUB,EPS0,EPS1,MODE,MERROR)
      IF(MERROR .NE. 0) THEN
       NERROR=MERROR
       RETURN
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=SQRT(DVMPY(MFR,B(1),B(2),B(1),B(2)))

************************************************************************
*   COMPUTE THE VALUE  PHI1  OF THE OBJECTIVE FUNCTION
************************************************************************

      PHI1=0
      IX=1
      DO 200 I=1,N
       CALL SUB(K,X(IX),M,A,F0,ZZ,0,MERROR)
       IF(F0 .LE. 0  .OR.  MERROR .NE. 0) THEN
        NERROR=3
        RETURN
       ENDIF
       PHI1=PHI1-LOG(F0)
  200 IX=IX+NX


************************************************************************
*   PRINT LAST ITERATION RESULTS
************************************************************************

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,'DMAXLK')

      RETURN

      END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "d501n1.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501n1.F"
*
* $Id: d501n1.F,v 1.1.1.1 1996/04/01 15:02:19 mclareni Exp $
*
* $Log: d501n1.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "d501n1.F" 2
      SUBROUTINE D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,WORK,DPHI,DSCAL,LAMU,
     +                  F,DF,IAFR,MFR,SUB,EPS0,EPS,MODE,VERS,NERROR)

*************************************************************************
*   LEAMAX, VERSION: 15.03.1993
*************************************************************************
*
*   THIS ROUTINE COMPUTES FUNCTION VALUES, DERIVATIVES, THE GRADIENT,
*   AND THE SCALING PARAMETERS. IT ALSO DETERMINES THE ACTIVE SET OF
*   CONSTRAINTS AND THE LAGRANGE MULTIPLIER.
*
************************************************************************


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 24 "d501n1.F" 2

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 25 "d501n1.F" 2
     +    LAMU
      CHARACTER VERS*6
      DIMENSION A(*),AL(*),AU(*),X(*),Y(*),SY(*),WORK(*),DPHI(*)
      DIMENSION DSCAL(*),LAMU(*),F(*),DF(N,*),IAFR(*)
      EXTERNAL SUB

      PARAMETER (Z0 = 0)

************************************************************************
*   COMPUTE INITIAL VALUES
************************************************************************

      HREL=SQRT(EPS0)
      HABS=10*EPS0

      NERROR=0

************************************************************************
*   COMPUTE FUNCTION VALUES AND DERIVATIVES (IF MODE NOTEQUAL ZERO)
************************************************************************

      CALL D501SF(VERS,SUB,MODE,M,A,N,F,DF,K,NX,X,Y,SY,WORK(N+1),NERROR)
      IF(NERROR .NE. 0) RETURN

      IF(MODE .EQ. 0) THEN

************************************************************************
*    APPROXIMATE DERIVATIVES
************************************************************************

       DO 10 J=1,M
       H =ABS(A(J))*HREL+HABS
       IF (A(J)+H .GT. AU(J)) H=-H
       A(J)=A(J)+H
       CALL D501SF
     +      (VERS,SUB,MODE,M,A,N,WORK,DF,K,NX,X,Y,SY,WORK(N+1),NERROR)
       IF(NERROR .NE. 0) RETURN
       A(J)=A(J)-H
       CALL DVSUB(N,WORK(1),WORK(2),F(1),F(2),DF(1,J),DF(2,J))
   10  CALL DVSCL(N,1/H,DF(1,J),DF(2,J),DF(1,J),DF(2,J))
      ENDIF

************************************************************************
*   COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION
************************************************************************

      CALL DMMPY(M,N,DF(1,1),DF(2,1),DF(1,2),F(1),F(2),DPHI(1),DPHI(2))

************************************************************************
*   DETERMINE THE DIAGONAL MATRIX   DSCAL   FOR SCALING THE PROBLEM
************************************************************************

      DO 30 I=1,M
      AI=0
      DO 20 J=1,N
   20 AI=AI+DF(J,I)**2
   30 DSCAL(I)=MAX(DSCAL(I),SQRT(AI))

************************************************************************
*     DETERMINE FREE VARIABLES AND STORE THEIR INDECES IN IAFR
*     DETERMINE LAGRANGE-MULTIPLIER   LAMU
************************************************************************

      GR=0
      DO 40 I=1,MFR
   40 GR=GR+(DSCAL(I)*A(IAFR(I)))**2
      GR=HREL*SQRT(GR)

      CALL DVSET(M,Z0,LAMU(1),LAMU(2))

      MFR=0

      DO 50 I=1,M
      IF(AU(I)-AL(I) .LT. EPS*(ABS(AU(I))+ABS(AL(I)))+2*HABS) THEN
        A(I)=AU(I)
        LAMU(I)=DPHI(I)
      ELSE
       IF(A(I) .GE. AU(I)-(EPS * ABS(AU(I)) + HABS )) THEN
        A(I)=AU(I)
        IF(DPHI(I) .GT. -GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE IF(A(I) .LE. AL(I)+(EPS * ABS(AL(I)) + HABS )) THEN
        A(I)=AL(I)
        IF(DPHI(I) .LT. GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE
        MFR=MFR+1
        IAFR(MFR)=I
       ENDIF
      ENDIF

   50 CONTINUE

************************************************************************
*   DELETE ROWS OF  DSCAL  AND COLUMNS  OF  DF
*   WHICH BELONG TO NON-FREE VARIABLES
************************************************************************

       DO 60 I=1,MFR
       DSCAL(I)=DSCAL(IAFR(I))
       DO 60 L=1,N
   60  DF(L,I)=DF(L,IAFR(I))

      RETURN
      END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "d501n2.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501n2.F"
*
* $Id: d501n2.F,v 1.1.1.1 1996/04/01 15:02:19 mclareni Exp $
*
* $Log: d501n2.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "d501n2.F" 2
      SUBROUTINE D501N2(K,N,M,A,AL,AU,X,NX,WORK,B,DPHI,DSCAL,LAMU,
     1                  AM,COV,IAFR,MFR,SUB,EPS0,EPS,MODE,NERROR)


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 14 "d501n2.F" 2

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 15 "d501n2.F" 2
     +   LAMU
      DIMENSION A(*),AL(*),AU(*),X(*),WORK(*),B(*),DPHI(*),DSCAL(*)
      DIMENSION LAMU(*),AM(M,*),COV(M,*),IAFR(*)
      PARAMETER (Z0 = 0)

*************************************************************************
*   LEAMAX, VERSION: 15.03.1993
*************************************************************************
*
*   THIS ROUTINE COMPUTES THE GRADIENT, THE JACOBIAN, AND IT SETS UP
*   THE MATRIX FOR THE NORMAL EQUATIONS. IT ALSO DETERMINES THE ACTIVE
*   SET OF CONSTRAINTS AND THE LAGRANGE-MULTIPLIER.
*
************************************************************************

************************************************************************
*   SET INITIAL VALUES
************************************************************************

      HREL=SQRT(EPS0)
      HABS=10*EPS0

************************************************************************
*   COMPUTE THE GRADIENT   B  OF THE OBJECTIVE FUNCTION
*   COMPUTE AN APPROXIMATION  AM  OF THE SECOND DERIVATIVE (THE HESSIAN)
*   OF THE OBJECTIVE FUNCTION
************************************************************************

      NERROR=0
      CALL DVSET(M,Z0,B(1),B(2))
      CALL DMSET(M,M,Z0,AM(1,1),AM(1,2),AM(2,1))
      IX=1

      DO 30 I=1,N

      CALL SUB(K,X(IX),M,A,F0,WORK,MODE,NERROR)
      IF(NERROR .NE. 0  .OR.  F0 .LE. 0) THEN
       NERROR=3
       RETURN
      ENDIF

      IF(MODE .EQ. 0) THEN

************************************************************************
*   APPROXIMATE DERIVATIVES
************************************************************************

       DO 10 J=1,M
        H =ABS(A(J))*HREL+HABS
        IF (A(J)+H .GT. AU(J)) H =-H
        A(J)=A(J)+H
        CALL SUB(K,X(IX),M,A,FH,WORK,MODE,NERROR)
        IF(NERROR .NE. 0) THEN
         NERROR=3
         RETURN
        ENDIF
        A(J)=A(J)-H
   10   WORK(J)=(FH-F0)/H
       ENDIF

       CALL DVSCL(M,1/F0,WORK(1),WORK(2),WORK(1),WORK(2))
       CALL DVSUB(M,B(1),B(2),WORK(1),WORK(2),B(1),B(2))

       DO 20 L=1,M
       DO 20 J=L,M
   20  AM(L,J)=AM(L,J)+WORK(L)*WORK(J)

   30  IX=IX+NX

       CALL DMUTL(M,AM(1,1),AM(1,2),AM(2,1))

************************************************************************
*   COPY THE GRADIENT OF THE OBJECTIVE FUNCTION TO  DPHI
************************************************************************

      CALL DVCPY(M,B(1),B(2),DPHI(1),DPHI(2))

************************************************************************
*   DETERMINE THE DIAGONAL MATRIX  DSCAL  FOR SCALING THE PROBLEM
************************************************************************

      DO 40 I=1,M
   40 DSCAL(I)=MAX(DSCAL(I),SQRT(AM(I,I)))

************************************************************************
*   DETERMINE FREE VARIABLES AND STORE THEIR INDICES IN IAFR
*   DETERMINE LAGRANGE MULTIPLIER  LAMU
************************************************************************

      GR=0
      DO 50 I=1,MFR
   50 GR=GR+(DSCAL(I)*A(IAFR(I)))**2
      GR=HREL*SQRT(GR)
      CALL DVSET(M,Z0,LAMU(1),LAMU(2))

      MFR=0

      DO 60 I=1,M
      IF(AU(I)-AL(I) .LT. EPS*(ABS(AU(I))+ABS(AL(I)))+2*HABS) THEN
        A(I)=AU(I)
        LAMU(I)=DPHI(I)
      ELSE
       IF(A(I) .GE. AU(I)-(EPS * ABS(AU(I)) + HABS)) THEN
        A(I)=AU(I)
        IF(DPHI(I) .GT. -GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE IF(A(I) .LE. AL(I)+(EPS * ABS(AL(I)) + HABS)) THEN
        A(I)=AL(I)
        IF(DPHI(I) .LT. GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE
        MFR=MFR+1
        IAFR(MFR)=I
       ENDIF
      ENDIF

   60 CONTINUE

***********************************************************************
*   DELETE ROWS AND COLUMNS OF  AM  AND  B  WHICH BELONG TO NON-FREE
*   VARIABLES
************************************************************************

      IF(MFR .EQ. 0 .OR. MFR .EQ. M) THEN
       MFC=M
      ELSE
       MFC=MFR
       DO 70 I =1,MFR
       B(I)=B(IAFR(I))
       DSCAL(I)=DSCAL(IAFR(I))
       DO 70 L = 1,M
   70  AM(L,I)=AM(L,IAFR(I))
       DO 80 I=1,MFR
       DO 80 L=1,M
   80  AM(I,L)=AM(IAFR(I),L)
      ENDIF

      CALL DMCPY(MFC,MFC,AM(1,1),AM(1,2),AM(2,1),
     +                   COV(1,1),COV(1,2),COV(2,1))
      RETURN

      END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "d501p1.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501p1.F"
*
* $Id: d501p1.F,v 1.1.1.1 1996/04/01 15:02:19 mclareni Exp $
*
* $Log: d501p1.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "d501p1.F" 2
      SUBROUTINE D501P1(K,M,NC,X,NX,Y,SY,MODE,EPS0,EPS,MAXIT,IPRT,
     +                  N,A,AL,AU,NERROR,VERS)

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE CHECKS THE VALUES OF INPUT PARAMETERS OF THE
*   SUBROUTINES  DSUMSQ, DFUNFT, DMAXLK  DEPENDING ON THE VALUE OF
*   THE PARAMETER  VERS.
*   IF  IPRT < 0  ALL VALUES OF THE INPUT PARAMETERS ARE PRINTED.
*
*************************************************************************


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 25 "d501p1.F" 2
      CHARACTER VERS*6,TIT1(3)*35,TIT2(0:1)*35
      DIMENSION X(*),Y(*),SY(*),A(*),AL(*),AU(*)
      PARAMETER (Z1 = 1, R10 = Z1/10)

      DATA TIT1(1) /'MINIMIZATION OF A SUM OF SQUARES'/
      DATA TIT1(2) /'LEAST-SQUARES DATA FITTING'/
      DATA TIT1(3) /'MAXIMUM LIKELIHOOD ESTIMATION'/

      DATA TIT2(0) /'APPROXIMATE DERIVATIVES (MODE = 0)'/
      DATA TIT2(1) /'ANALYTICAL DERIVATIVES (MODE = 1)'/

      IF(IPRT .NE. 0) THEN
       IF(VERS .EQ. 'DSUMSQ') IV=1
       IF(VERS .EQ. 'DFUNFT') IV=2
       IF(VERS .EQ. 'DMAXLK') IV=3

       WRITE(6,1000) VERS,TIT1(IV),TIT2(MODE)
      ENDIF

************************************************************************
*   PRINT INPUT PARAMETERS (IF IPRT .LT. 0)
************************************************************************

      IF(IPRT .LT. 0) THEN
       WRITE(6,1010) VERS,M,N
       IF(VERS .NE. 'DMAXLK') WRITE(6,1020) NC
       IF(VERS .NE. 'DSUMSQ') WRITE(6,1030) K,NX
       WRITE(6,1040) MAXIT,MODE,IPRT,EPS
      ENDIF

************************************************************************
*   CHECK VALUES OF INPUT PARAMETERS, AND PRINT THEM (IF IPRT .LT. 0)
************************************************************************

      NERROR=0

      IF(     MAXIT .LT. 1
     1   .OR. K     .LT. 1
     2   .OR. N     .LT. 1
     3   .OR. M     .LT. N
     4   .OR. NC    .LT. N
     5   .OR. NX    .LT. K ) THEN
       NERROR=1
       RETURN
      ENDIF

      IF(IPRT .LT. 0) THEN
       WRITE(6,1050) (AL(I), I=1,N)
       WRITE(6,1060) (AU(I), I=1,N)
       WRITE(6,1070) (A(I),  I=1,N)
      IF(VERS .NE. 'DSUMSQ')WRITE(6,1080)((X(I),I=L,M*NX,NX),L=1,K)
       IF(VERS .EQ. 'DFUNFT') THEN
        WRITE(6,1090) (Y(I), I=1,M)
        WRITE(6,1100) (SY(I),I=1,M)
       ENDIF
      ENDIF

      DO 10 I=1,N
      IF(AL(I) .GT. AU(I)) THEN
       NERROR=1
       RETURN
      ENDIF
   10 CONTINUE

************************************************************************
*   IF VALUES OF THE PARAMETERS A, SY, MODE OR EPS ARE NOT PRACTICABLE
*   SET RECOMMENDED VALUES FOR THIS PARAMETERS
************************************************************************

      DO 20 I=1,N
      IF(A(I) .GT. AU(I)) A(I)=AU(I)
      IF(A(I) .LT. AL(I)) A(I)=AL(I)
   20 CONTINUE

      IF (VERS .EQ. 'DFUNFT') THEN
       DO 30 I = 1,M
       IF(SY(I) .LE. 0) THEN
        NERROR=1
        RETURN
       ENDIF
   30  CONTINUE
      ENDIF

CC    IF(STEP .LE. 0) STEP=1
      IF(MODE .NE. 1) MODE=0
      IF (EPS .LT. EPS0  .OR.  EPS .GT. R10) EPS=10*EPS0

      RETURN

 1000 FORMAT(7(/),30X,'MATHLIB PACKAGE   D501   VERSION 15.03.93'//
     1       30X,'PACKAGE LEAMAX  ****  ROUTINE ',A6,' ****'///
     2       15X,A35,A35//)
 1010 FORMAT(' INPUT  OF  ',A6,' :'//'  M :',I5,6X,'N :',I5)
 1020 FORMAT('  NC:',I5)
 1030 FORMAT('  K :',I5,6X,'NX:',I5)
 1040 FORMAT('  MAXIT :',I5,8X,'MODE :',I5,8X,'IPRT :',I5/
     +       '  EPS   :',1PD11.1/)
 1050 FORMAT(/'  AL :',/(5(1PD15.5)))
 1060 FORMAT( '  AU :',/(5(1PD15.5)))
 1070 FORMAT( '  A :', /(5(1PD15.5)))
 1080 FORMAT( '  X :', /(5(1PD15.5)))
 1090 FORMAT( '  Y :', /(5(1PD15.5)))
 1100 FORMAT( '  SY :',/(5(1PD15.5)))

      END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "d501p2.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501p2.F"
*
* $Id: d501p2.F,v 1.1.1.1 1996/04/01 15:02:20 mclareni Exp $
*
* $Log: d501p2.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "d501p2.F" 2
      SUBROUTINE D501P2(LRP,N,A,B,C,LAMU,PHI,DPHINO,ITER,LFN,MODE,VERS)

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE CONTROLS THE PRINTING OF THE PACKAGE LEAMAX.
*
************************************************************************


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 21 "d501p2.F" 2

# 1 "/usr/include/gen/def64.inc" 1 3 4
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*







      DOUBLE PRECISION
# 22 "d501p2.F" 2
     +   LAMBDA,LAMU
      LOGICAL LFN,LRP
      CHARACTER VERS*6,TIT(2)*18
      DIMENSION A(*),B(*),C(*),LAMU(*)


      DATA TIT(1),TIT(2) /' ','STANDARD DEVIATION'/

      IF(.NOT.LRP) THEN

       WRITE(6,1030)
       LRP=.TRUE.
      ENDIF

      IF(LFN) THEN
       WRITE(6,1010) 'END:',ITER,PHI,DPHINO
       IF(VERS .EQ. 'DFUNFT') THEN
        WRITE(6,1020) TIT(2)(1:8),TIT(2)(10:18)
       ELSE
        WRITE(6,1020) TIT(1)(1:8),TIT(1)(10:18)
       ENDIF
      ELSE
       WRITE(6,1010) '    ',ITER,PHI,DPHINO
       WRITE(6,1020) TIT(1)(1:8),TIT(1)(10:18)
      ENDIF

      IF(LFN .AND. VERS .EQ. 'DFUNFT') THEN
       WRITE(6,1040) (I,A(I),B(I),LAMU(I),C(I), I=1,N)
      ELSE
       WRITE(6,1050) (I,A(I),B(I),LAMU(I), I=1,N)
      ENDIF

      RETURN

 1010 FORMAT(/6X,A4,' ITERATION',I5,3X,'PHI = ',1PD12.5,6X,
     1       'GNO = ',1PD12.5/)
 1020 FORMAT(12X,'PARAMETER',7X,'PARAMETER',9X,'GRADIENT',
     1       10X,'LAGRANGE',8X,A8/
     2       14X,'NUMBER',10X,'VALUE',28X,'MULTIPLIER',7X,A9/)
 1030 FORMAT(//' ITERATION'//11X,'PHI = VALUE OF OBJECTIVE FUNCTION',
     1        10X,'GNO = NORM OF GRADIENT')
 1040 FORMAT (15X,I3,4X,1PD17.5,1PD17.5,1PD17.5,1PD17.5)
 1050 FORMAT (15X,I3,4X,1PD17.5,1PD17.5,1PD17.5)

      END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "d501sf.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "d501sf.F"
*
* $Id: d501sf.F,v 1.1.1.1 1996/04/01 15:02:20 mclareni Exp $
*
* $Log: d501sf.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "d501sf.F" 2
      SUBROUTINE D501SF (VERS,SUB,MODE,M,A,N,F,DF,K,NX,X,Y,SY,W,NERROR)

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE COMPUTES FUNCTION VALUES AND DERIVATIVES DEPENDING ON
*   THE VALUE OF THE PARAMETER  VERS.
*
*************************************************************************


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 22 "d501sf.F" 2
      DIMENSION A(*),F(*),DF(N,*),X(*),Y(*),SY(*),W(*)
      CHARACTER VERS*6

      NERROR=0

      IF (VERS .EQ. 'DSUMSQ') THEN
       CALL SUB (M,A,N,F,DF,MODE,NERROR)
       IF (NERROR .NE. 0) NERROR=3
       RETURN
      ENDIF

      IF (VERS .EQ. 'DFUNFT') THEN
       IX=1
       DO 20 I=1,N
        CALL SUB (K,X(IX),M,A,SF,W,MODE,NERROR)
        IF (NERROR .NE. 0) THEN
         NERROR=3
         RETURN
        ENDIF
        F(I)=(Y(I)-SF)/SY(I)
        IX=IX+NX
       IF (MODE .EQ. 0) GOTO 20
       DO 10 J=1,M
   10  DF(I,J)=-W(J)/SY(I)
   20  CONTINUE
       RETURN
      ENDIF

      END
*CMZ :          21/11/2017  14.35.24  by  Michael Scheer
*-- Author :
# 1 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F"
*
* $Id: bsja64.F,v 1.1.1.1 1996/04/01 15:02:08 mclareni Exp $
*
* $Log: bsja64.F,v $
* Revision 1.1.1.1 1996/04/01 15:02:08 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F" 2

      SUBROUTINE DBSJA(X,A,NMAX,ND,B)





# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
# 18 "/opt/cern/pro/src/mathlib/gen/c/bsja64.F" 2
      REAL SX,D,T,Q,U,V,TC(11)
      CHARACTER*80 ERRTXT
      CHARACTER NAMEJ*(*),NAMEI*(*)

      PARAMETER (NAMEJ = 'BSJA/DBSJA',
     1 NAMEI = 'BSIA/DBSIA')
      EXTERNAL DGAMMA





      LOGICAL LJA,LIA,LEV,LER
      DIMENSION B(0:*),BA(0:100),RR(0:100)

      PARAMETER (Z1 = 1, HF = Z1/2, Z10 = 10)

      DATA TC / 5.7941 E-5,-1.76148E-3, 2.08645E-2,-1.29013E-1,
     1 8.5777 E-1, 1.0125 E+0, 7.75 E-1, 2.3026 E+0,
     2 1.3863 E+0, 7.3576 E-1, 1.3591 E+0/

      LJA=.TRUE.
      LIA=.FALSE.
      SGN=-1
      GO TO 9


      ENTRY DBSIA(X,A,NMAX,ND,B)




      LJA=.FALSE.
      LIA=.TRUE.
      SGN=1

    9 LER=.FALSE.
      IF(X .LE. 0) THEN
       WRITE(ERRTXT,101) X
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.1',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.1',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. A .AND. A .LT. 1)) THEN
       WRITE(ERRTXT,102) A
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.2',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.2',ERRTXT)
       LER=.TRUE.
      ELSEIF(ABS(NMAX) .GT. 100) THEN
       WRITE(ERRTXT,103) NMAX
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.3',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.3',ERRTXT)
       LER=.TRUE.
      END IF
      IF(LER) RETURN
      EPS=HF*Z10**(-ND)
      NMX=ABS(NMAX)
      IF(NMAX .LE. 0) NMX=1
      DO 5 N = 0,NMX
      RR(N)=0
    5 BA(N)=0
      D=TC(8)*ND+TC(9)
      SX=X
      Q=0
      IF(NMX .GT. 0) THEN
       V=0.5*D/NMX
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 6 I = 2,6
    6 T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       ENDIF
       Q=NMX*T
      ENDIF

      F=(HF*X)**A/DGAMMA(1+A)




      T=1
      V=TC(10)*D/SX
      IF(LIA) THEN
       F=EXP(X)*F
       V=V-TC(10)
      ENDIF
      IF(LJA .OR. LIA .AND. X .LT. D) THEN
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 7 I = 2,6
    7 T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       ENDIF
      ENDIF
      NU=1+MAX(Q,TC(11)*SX*T)

      MU=-1
    2 MU=MU+1
      AL=1
      IF(LJA) THEN
       DO 3 N = 1,NU/2
       XN=N
    3 AL=AL*(XN+A)/(XN+1)
       R=0
       S=0
       LEV=.TRUE.
       DO 4 N = 2*(NU/2),1,-1
       XN=N
       XA=XN+A
       R=1/(2*XA/X-R)
       IF(N .LE. NMX) RR(N-1)=R
       IF(LEV) THEN
        AL=AL*(XN+2)/(XA+A)
        S=R*(AL*XA+S)
       ELSE
        S=R*S
       ENDIF
       LEV=.NOT.LEV
    4 CONTINUE
      ELSE
       DO 23 N = 1,NU
       XN=N
   23 AL=AL*(XN+2*A)/(XN+1)
       R=0
       S=0
       DO 24 N = NU,1,-1
       XN=N
       XA=XN+A
       XA2=XA+XA
       R=1/(XA2/X+R)
       IF(N .LE. NMX) RR(N-1)=R
       AL=AL*(XN+1)/(XA+A)
       S=R*(XA2*AL+S)
   24 CONTINUE
      ENDIF
      B(0)=F/(1+S)
      DO 10 N = 0,NMX-1
   10 B(N+1)=RR(N)*B(N)
      DO 11 N = 0,NMX
      IF(ABS(B(N)-BA(N)) .GT. EPS*ABS(B(N))) THEN
       DO 12 M = 0,NMX
   12 BA(M)=B(M)
       NU=NU+5
       IF(MU .LE. 50) GO TO 2
       WRITE(ERRTXT,104) X,A
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.4',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.4',ERRTXT)
       RETURN
      ENDIF
   11 CONTINUE
      IF(NMAX .LT. 0) THEN
       AL=2/X
       B(1)=AL*A*B(0)+SGN*B(1)
       DO 13 N = 1,-NMAX-1
   13 B(N+1)=AL*(A-N)*B(N)+SGN*B(N-1)
      ENDIF
      RETURN
  101 FORMAT('ILLEGAL ARGUMENT X = ',1P,D15.8)
  102 FORMAT('ILLEGAL ORDER A = ',1P,D15.8)
  103 FORMAT('ILLEGAL NMAX = ',I5)
  104 FORMAT('NO CONVERGENCE FOR X = ',1P,D15.8,' A = ',D15.8,
     1 ' TRY SMALLER ND')
      END
*CMZ :          27/11/2017  14.29.31  by  Michael Scheer
*-- Author :
      ! TAKEN from CERNLIB C340
      FUNCTION  BSIR3(X,NU)
      CHARACTER*80 ERRTXT
      CHARACTER NAMEI*(*),NAMEK*(*),NAMEIE*(*),NAMEKE*(*)
      PARAMETER (NAMEI = 'BSIR3', NAMEIE = 'EBSIR3')
      PARAMETER (NAMEK = 'BSKR3', NAMEKE = 'EBSKR3')
      LOGICAL LEX

      DIMENSION BC(0:23,2),CC(0:15,2),PP(-2:2),GG(-2:2)

      PARAMETER (EPS = 1D-15)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (W3 = 1.73205 08075 68877 29D0)
      PARAMETER (G1 = 2.67893 85347 07747 63D0)
      PARAMETER (G2 = 1.35411 79394 26400 42D0)
      PARAMETER (PIH = PI/2, RPIH = 2/PI, RPI2 = 1/(2*PI))
      PARAMETER (C1 = 2*PI/(3*W3))
      PARAMETER (GM = 3*(1/G2-3/G1)/2, GP = (3/G1+1/G2)/2)

      DATA PP(-2) /-0.66666 66666 66666 67D0/
      DATA PP(-1) /-0.33333 33333 33333 33D0/
      DATA PP( 1) / 0.33333 33333 33333 33D0/
      DATA PP( 2) / 0.66666 66666 66666 67D0/

      DATA GG(-2) / 0.37328 21739 07395 23D0/
      DATA GG(-1) / 0.73848 81116 21648 31D0/
      DATA GG( 1) / 1.11984 65217 22185 68D0/
      DATA GG( 2) / 1.10773 21674 32472 47D0/

      DATA BC( 0,1) / 1.00458 61710 93207 35D0/
      DATA BC( 1,1) / 0.00467 34791 99873 60D0/
      DATA BC( 2,1) / 0.00009 08034 04815 04D0/
      DATA BC( 3,1) / 0.00000 37262 16110 59D0/
      DATA BC( 4,1) / 0.00000 02520 73237 90D0/
      DATA BC( 5,1) / 0.00000 00227 82110 77D0/
      DATA BC( 6,1) / 0.00000 00012 91332 28D0/
      DATA BC( 7,1) /-0.00000 00006 11915 16D0/
      DATA BC( 8,1) /-0.00000 00003 75616 85D0/
      DATA BC( 9,1) /-0.00000 00001 16415 46D0/
      DATA BC(10,1) /-0.00000 00000 14443 25D0/
      DATA BC(11,1) / 0.00000 00000 05373 69D0/
      DATA BC(12,1) / 0.00000 00000 03074 27D0/
      DATA BC(13,1) / 0.00000 00000 00297 66D0/
      DATA BC(14,1) /-0.00000 00000 00265 20D0/
      DATA BC(15,1) /-0.00000 00000 00091 37D0/
      DATA BC(16,1) / 0.00000 00000 00015 52D0/
      DATA BC(17,1) / 0.00000 00000 00014 12D0/
      DATA BC(18,1) /-0.00000 00000 00000 23D0/
      DATA BC(19,1) /-0.00000 00000 00001 98D0/
      DATA BC(20,1) /-0.00000 00000 00000 13D0/
      DATA BC(21,1) / 0.00000 00000 00000 29D0/
      DATA BC(22,1) / 0.00000 00000 00000 03D0/
      DATA BC(23,1) /-0.00000 00000 00000 05D0/

      DATA BC( 0,2) / 0.99363 49867 16925 14D0/
      DATA BC( 1,2) /-0.00646 71526 00616 03D0/
      DATA BC( 2,2) /-0.00010 60188 22351 55D0/
      DATA BC( 3,2) /-0.00000 41406 57716 24D0/
      DATA BC( 4,2) /-0.00000 02916 95418 21D0/
      DATA BC( 5,2) /-0.00000 00365 71574 33D0/
      DATA BC( 6,2) /-0.00000 00075 81590 37D0/
      DATA BC( 7,2) /-0.00000 00019 23008 52D0/
      DATA BC( 8,2) /-0.00000 00004 20438 80D0/
      DATA BC( 9,2) /-0.00000 00000 39372 04D0/
      DATA BC(10,2) / 0.00000 00000 19007 44D0/
      DATA BC(11,2) / 0.00000 00000 10137 64D0/
      DATA BC(12,2) / 0.00000 00000 01331 30D0/
      DATA BC(13,2) /-0.00000 00000 00676 92D0/
      DATA BC(14,2) /-0.00000 00000 00311 72D0/
      DATA BC(15,2) / 0.00000 00000 00011 87D0/
      DATA BC(16,2) / 0.00000 00000 00040 21D0/
      DATA BC(17,2) / 0.00000 00000 00004 78D0/
      DATA BC(18,2) /-0.00000 00000 00004 74D0/
      DATA BC(19,2) /-0.00000 00000 00001 16D0/
      DATA BC(20,2) / 0.00000 00000 00000 59D0/
      DATA BC(21,2) / 0.00000 00000 00000 21D0/
      DATA BC(22,2) /-0.00000 00000 00000 08D0/
      DATA BC(23,2) /-0.00000 00000 00000 03D0/

      DATA CC( 0,1) / 0.99353 64122 76093 39D0/
      DATA CC( 1,1) /-0.00631 44392 60798 63D0/
      DATA CC( 2,1) / 0.00014 30095 80961 13D0/
      DATA CC( 3,1) /-0.00000 57870 60592 03D0/
      DATA CC( 4,1) / 0.00000 03265 50333 20D0/
      DATA CC( 5,1) /-0.00000 00231 23231 95D0/
      DATA CC( 6,1) / 0.00000 00019 39555 14D0/
      DATA CC( 7,1) /-0.00000 00001 85897 89D0/
      DATA CC( 8,1) / 0.00000 00000 19868 42D0/
      DATA CC( 9,1) /-0.00000 00000 02326 79D0/
      DATA CC(10,1) / 0.00000 00000 00294 68D0/
      DATA CC(11,1) /-0.00000 00000 00039 95D0/
      DATA CC(12,1) / 0.00000 00000 00005 75D0/
      DATA CC(13,1) /-0.00000 00000 00000 87D0/
      DATA CC(14,1) / 0.00000 00000 00000 14D0/
      DATA CC(15,1) /-0.00000 00000 00000 02D0/

      DATA CC( 0,2) / 1.00914 95380 72789 40D0/
      DATA CC( 1,2) / 0.00897 12068 42483 60D0/
      DATA CC( 2,2) /-0.00017 13895 98261 54D0/
      DATA CC( 3,2) / 0.00000 65547 92549 82D0/
      DATA CC( 4,2) /-0.00000 03595 19190 49D0/
      DATA CC( 5,2) / 0.00000 00250 24412 19D0/
      DATA CC( 6,2) /-0.00000 00020 74924 13D0/
      DATA CC( 7,2) / 0.00000 00001 97223 67D0/
      DATA CC( 8,2) /-0.00000 00000 20946 47D0/
      DATA CC( 9,2) / 0.00000 00000 02440 93D0/
      DATA CC(10,2) /-0.00000 00000 00307 91D0/
      DATA CC(11,2) / 0.00000 00000 00041 61D0/
      DATA CC(12,2) /-0.00000 00000 00005 97D0/
      DATA CC(13,2) / 0.00000 00000 00000 91D0/
      DATA CC(14,2) /-0.00000 00000 00000 14D0/
      DATA CC(15,2) / 0.00000 00000 00000 02D0/

      LEX=.FALSE.
      GO TO 8

      ENTRY  EBSIR3(X,NU)
      LEX=.TRUE.

    8 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .OR. NU .LT. 0 .AND. X .LE. 0
     1   .OR. NU .GT. 0 .AND. X .LT. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEI ,'C340.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEIE,'C340.1',ERRTXT)
      ELSEIF(X .EQ. 0) THEN
       S=0
      ELSEIF(X .LT. 8) THEN
       Y=(HF*X)**2
       XN=PP(NU)
       XL=XN+2
       A0=1
       A1=1+2*Y/((XL+1)*(XN+1))
       A2=1+Y*(4+3*Y/((XL+2)*(XN+2)))/((XL+3)*(XN+1))
       B0=1
       B1=1-Y/(XL+1)
       B2=1-Y*(1-Y/(2*(XL+2)))/(XL+3)
       T1=3+XL
       V1=3-XL
       V3=XL-1
       V2=V3+V3
       C=0
       DO 33 N = 3,30
       C0=C
       T1=T1+2
       T2=T1-1
       T3=T2-1
       T4=T3-1
       T5=T4-1
       T6=T5-1
       V1=V1+1
       V2=V2+1
       V3=V3+1
       U1=N*T4
       E=V3/(U1*T3)
       U2=E*Y
       F1=1+Y*V1/(U1*T1)
       F2=(1+Y*V2/(V3*T2*T5))*U2
       F3=-Y*Y*U2/(T4*T5*T5*T6)
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 34
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   33  CONTINUE
   34  S=GG(NU)*(HF*X)**PP(NU)*C
       IF(LEX) S=EXP(-X)*S
      ELSE
       R=1/X
       W=SQRT(RPI2*R)
       H=16*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 23,0,-1
       B0=BC(I,MU)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=W*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(X)*S
       T=0
       IF(NU .LT. 0) THEN
        H=10*R-1
        ALFA=H+H
        B1=0
        B2=0
        DO 3 I = 15,0,-1
        B0=CC(I,MU)+ALFA*B1-B2
        B2=B1
    3   B1=B0
        R=EXP(-X)
        T=W3*W*R*(B0-H*B2)
        IF(LEX) T=R*T
       END IF
       S=S+T
      END IF
      GO TO 99

      ENTRY  BSKR3(X,NU)
      LEX=.FALSE.
      GO TO 9

      ENTRY  EBSKR3(X,NU)
      LEX=.TRUE.

    9 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .OR. X .LE. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEK ,'C340.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEKE,'C340.1',ERRTXT)
      ELSEIF(X .LE. 1) THEN
       A0=PP(-1)
       B=HF*X
       D=-LOG(B)
       F=A0*D
       E=EXP(F)
       G=(GM*A0+GP)*E
       BK=C1*(HF*GM*(E+1/E)+GP*D*SINH(F)/F)
       F=BK
       E=A0**2
       P=HF*C1*G
       Q=HF/G
       C=1
       D=B**2
       BK1=P
       DO 11 N = 1,15
       FN=N
       F=(FN*F+P+Q)/(FN**2-E)
       C=C*D/FN
       P=P/(FN-A0)
       Q=Q/(FN+A0)
       G=C*(P-FN*F)
       H=C*F
       BK=BK+H
       BK1=BK1+G
       IF(H*BK1+ABS(G)*BK .LE. EPS*BK*BK1) GO TO 12
   11  CONTINUE
   12  S=BK
       IF(MU .EQ. 2) S=BK1/B
       IF(LEX) S=EXP(X)*S
      ELSEIF(X .LE. 5) THEN
       XN=4*PP(MU)**2
       A=9-XN
       B=25-XN
       C=768*X**2
       C0=48*X
       A0=1
       A1=(16*X+7+XN)/A
       A2=(C+C0*(XN+23)+XN*(XN+62)+129)/(A*B)
       B0=1
       B1=(16*X+9-XN)/A
       B2=(C+C0*B)/(A*B)+1
       C=0
       DO 24 N = 3,30
       C0=C
       FN=N
       FN2=FN+FN
       FN1=FN2-1
       FN3=FN1/(FN2-3)
       FN4=12*FN**2-(1-XN)
       FN5=16*FN1*X
       RAN=1/((FN2+1)**2-XN)
       F1=FN3*(FN4-20*FN)+FN5
       F2=28*FN-FN4-8+FN5
       F3=FN3*((FN2-5)**2-XN)
       A=(F1*A2+F2*A1+F3*A0)*RAN
       B=(F1*B2+F2*B1+F3*B0)*RAN
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 25
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   24  CONTINUE
   25  S=C/SQRT(RPIH*X)
       IF(.NOT.LEX) S=EXP(-X)*S
      ELSE
       R=1/X
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 13 I = 15,0,-1
       B0=CC(I,MU)+ALFA*B1-B2
       B2=B1
   13  B1=B0
       S=SQRT(PIH*R)*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(-X)*S
      END IF
   99  BSIR3=S
      RETURN
  101 FORMAT('INCORRECT ARGUMENT OR INDEX, X = ',1P,E15.6,' NU = ',I5)
      END
*CMZ :          02/05/2017  15.32.17  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
*CMZ :          18/03/2015  10.18.08  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX

# 1 "deqinv.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "deqinv.F"
*
* $Id: deqinv.F,v 1.1.1.1 1996/02/15 17:48:48 mclareni Exp $
*
* $Log: deqinv.F,v $
* Revision 1.1.1.1 1996/02/15 17:48:48 mclareni
* Kernlib
*
*
# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 10 "deqinv.F" 2
cmsh      SUBROUTINE DEQINV(N,A,IDIM,R,IFAIL,K,B)
      SUBROUTINE DEQINV(N,A,IDIM,ir,IFAIL,K,B)
cmsh      REAL R(N),T1,T2,T3
      integer ir(n)
      real T1,T2,T3
      DOUBLE PRECISION A(IDIM,N),B(IDIM,K),DET,TEMP,S,
     $ B1,B2,C11,C12,C13,C21,C22,C23,C31,C32,C33
      CHARACTER*6 NAME
      DATA NAME/'DEQINV'/,KPRNT/1/
C
C ******************************************************************
C
C REPLACES B BY THE SOLUTION X OF A*X=B, AND REPLACES A BY ITS IN-
C VERSE.
C
C N ORDER OF THE SQUARE MATRIX IN ARRAY A.
C
C A (DOUBLE PRECISION) TWO-DIMENSIONAL ARRAY CONTAINING
C AN N BY N MATRIX.
C
C IDIM FIRST DIMENSION PARAMETER OF ARRAYS A AND B.
C
C R (REAL) WORKING VECTOR OF LENGTH NOT LESS THAN N.
C
C IFAIL OUTPUT PARAMETER. IFAIL= 0 ... NORMAL EXIT.
C IFAIL=-1 ... SINGULAR MATRIX.
C
C K NUMBER OF COLUMNS OF THE MATRIX IN ARRAY B.
C
C B (DOUBLE PRECISION) TWO-DIMENSIONAL ARRAY CONTAINING
C AN N BY K MATRIX.
C
C CALLS ... DFACT, DFINV, F010PR, ABEND.
C
C ******************************************************************
C
C TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM).OR.(K.LT.1)) GO TO 10
C
C TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 9
      IFAIL=0
      IF(N.LT.3) GO TO 5
C
C N=3 CASE.
C
C COMPUTE COFACTORS.
      C11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      C12=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      C13=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      C21=A(3,2)*A(1,3)-A(3,3)*A(1,2)
      C22=A(3,3)*A(1,1)-A(3,1)*A(1,3)
      C23=A(3,1)*A(1,2)-A(3,2)*A(1,1)
      C31=A(1,2)*A(2,3)-A(1,3)*A(2,2)
      C32=A(1,3)*A(2,1)-A(1,1)*A(2,3)
      C33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      T1=ABS(SNGL(A(1,1)))
      T2=ABS(SNGL(A(2,1)))
      T3=ABS(SNGL(A(3,1)))
C
C (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      IF(T1.GE.T2) GO TO 1
         IF(T3.GE.T2) GO TO 2
C (PIVOT IS A21)
            TEMP=A(2,1)
            DET=C13*C32-C12*C33
            GO TO 3
    1 IF(T3.GE.T1) GO TO 2
C (PIVOT IS A11)
         TEMP=A(1,1)
         DET=C22*C33-C23*C32
         GO TO 3
C (PIVOT IS A31)
    2 TEMP=A(3,1)
         DET=C23*C12-C22*C13
C
C SET ELEMENTS OF INVERSE IN A.
    3 IF(DET.EQ.0D0) GO TO 11
      S=TEMP/DET
      A(1,1)=S*C11
      A(1,2)=S*C21
      A(1,3)=S*C31
      A(2,1)=S*C12
      A(2,2)=S*C22
      A(2,3)=S*C32
      A(3,1)=S*C13
      A(3,2)=S*C23
      A(3,3)=S*C33
C
C REPLACE B BY AINV*B.
      DO 4 J=1,K
         B1=B(1,J)
         B2=B(2,J)
         B(1,J)=A(1,1)*B1+A(1,2)*B2+A(1,3)*B(3,J)
         B(2,J)=A(2,1)*B1+A(2,2)*B2+A(2,3)*B(3,J)
         B(3,J)=A(3,1)*B1+A(3,2)*B2+A(3,3)*B(3,J)
    4 CONTINUE
      RETURN
C
    5 IF(N.LT.2) GO TO 7
C
C N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DET.EQ.0D0) GO TO 11
      S=1D0/DET
      C11 =S*A(2,2)
      A(1,2)=-S*A(1,2)
      A(2,1)=-S*A(2,1)
      A(2,2)=S*A(1,1)
      A(1,1)=C11
      DO 6 J=1,K
         B1=B(1,J)
         B(1,J)=C11*B1+A(1,2)*B(2,J)
         B(2,J)=A(2,1)*B1+A(2,2)*B(2,J)
    6 CONTINUE
      RETURN
C
C N=1 CASE.
C
    7 IF(A(1,1).EQ.0D0) GO TO 11
      A(1,1)=1D0/A(1,1)
      DO 8 J=1,K
         B(1,J)=A(1,1)*B(1,J)
    8 CONTINUE
      RETURN
C
C N.GT.3 CASES. FACTORIZE MATRIX, INVERT AND SOLVE SYSTEM.
C
cmsh    9 CALL DFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
    9 CALL DFACT(N,A,IDIM,ir,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
cmsh      CALL DFEQN(N,A,IDIM,R,K,B)
cmsh      CALL DFINV(N,A,IDIM,R)
      CALL DFEQN(N,A,IDIM,ir,K,B)
      CALL DFINV(N,A,IDIM,ir)
      RETURN
C
C ERROR EXITS.
C
   10 IFAIL=+1
      CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
   11 IFAIL=-1
      RETURN
C
      END
*CMZ :  1.16/04 17/04/2014  12.50.22  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*
* $Id: deqn.F,v 1.1.1.1 1996/02/15 17:48:49 mclareni Exp $
*
* $Log: deqn.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:49  mclareni
* Kernlib
*
*
cmsh      SUBROUTINE DEQN(N,A,IDIM,R,IFAIL,K,B)
cmsh      REAL R(N),T1,T2,T3
      SUBROUTINE DEQN(N,A,IDIM,IR,IFAIL,K,B)
      integer ir(n)
      REAL T1,T2,T3
      DOUBLE PRECISION A(IDIM,N),B(IDIM,K),DET,S,TEMP,
     $                 B1,Y1,Y2,L11,L21,L22,L31,L32,L33,U12,U13,U23
      CHARACTER*6 NAME
      DATA NAME/'DEQN'/,KPRNT/1/
C
C     ******************************************************************
C
C     REPLACES B BY THE SOLUTION X OF A*X=B, AFTER WHICH A IS UNDEFINED.
C
C     (PARAMETERS AS FOR DEQINV.)
C
C     CALLS ... DFACT, DFEQN, F010PR, ABEND.
C
C     ******************************************************************
C
C  TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM).OR.(K.LT.1)) GO TO 11
C
C  TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 10
      IFAIL=0
      IF(N.LT.3) GO TO 6
C
C  N=3 CASE.
C
C     FACTORIZE MATRIX A=L*U.
C     (FIRST PIVOT SEARCH)
      T1=ABS(SNGL(A(1,1)))
      T2=ABS(SNGL(A(2,1)))
      T3=ABS(SNGL(A(3,1)))
      IF(T1.GE.T2) GO TO 1
         IF(T3.GE.T2) GO TO 2
C        (PIVOT IS A21)
            M1=2
            M2=1
            M3=3
            GO TO 3
    1 IF(T3.GE.T1) GO TO 2
C     (PIVOT IS A11)
         M1=1
         M2=2
         M3=3
         GO TO 3
C     (PIVOT IS A31)
    2    M1=3
         M2=2
         M3=1
    3 TEMP=A(M1,1)
      IF(TEMP.EQ.0D0) GO TO 10
      L11=1D0/TEMP
      U12=L11*A(M1,2)
      U13=L11*A(M1,3)
      L22=A(M2,2)-A(M2,1)*U12
      L32=A(M3,2)-A(M3,1)*U12
C     (SECOND PIVOT SEARCH)
      IF( ABS(SNGL(L22)) .GE. ABS(SNGL(L32)) )  GO TO 4
         I=M2
         M2=M3
         M3=I
         TEMP=L22
         L22=L32
         L32=TEMP
    4 L21=A(M2,1)
      L31=A(M3,1)
      IF(L22.EQ.0D0) GO TO 10
      L22=1D0/L22
      U23=L22*(A(M2,3)-L21*U13)
      TEMP=A(M3,3)-L31*U13-L32*U23
      IF(TEMP.EQ.0D0) GO TO 10
      L33=1D0/TEMP
C
C     SOLVE L*Y=B AND U*X=Y.
      DO 5 J=1,K
         Y1=L11*B(M1,J)
         Y2=L22*(B(M2,J)-L21*Y1)
         B(3,J)=L33*(B(M3,J)-L31*Y1-L32*Y2)
         B(2,J)=Y2-U23*B(3,J)
         B(1,J)=Y1-U12*B(2,J)-U13*B(3,J)
    5 CONTINUE
      RETURN
C
    6 IF(N.LT.2) GO TO 8
C
C  N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DET.EQ.0D0) GO TO 12
      S=1D0/DET
      DO 7 J=1,K
         B1=B(1,J)
         B(1,J)=S*(A(2,2)*B1-A(1,2)*B(2,J))
         B(2,J)=S*(-A(2,1)*B1+A(1,1)*B(2,J))
    7 CONTINUE
      RETURN
C
C  N=1 CASE.
C
    8 IF(A(1,1).EQ.0D0) GO TO 12
      S=1D0/A(1,1)
      DO 9 J=1,K
         B(1,J)=S*B(1,J)
    9 CONTINUE
      RETURN
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND SOLVE SYSTEM.
C
cmsh   10 CALL DFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
   10 CALL DFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
cmsh      CALL DFEQN(N,A,IDIM,R,K,B)
      CALL DFEQN(N,A,IDIM,IR,K,B)
      RETURN
C
C  ERROR EXITS.
C
   11 IFAIL=+1
      CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
   12 IFAIL=-1
      RETURN
C
      END
*CMZ :          25/08/2014  16.19.01  by  Michael Scheer
*CMZ :  1.16/04 17/04/2014  12.54.46  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*# 1 "dfact.F"
*# 1 "<command-line>"
*# 1 "dfact.F"
*
* $Id: dfact.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfact.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*

*# 1 "kernnum/pilot.h" 1
*# 21 "kernnum/pilot.h"

*# 33 "kernnum/pilot.h"

*# 10 "dfact.F" 2
          SUBROUTINE          DFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
          INTEGER             IR(*),    IPAIRF
          DOUBLE PRECISION    A(IDIM,*),DET,      ZERO,     ONE,X,Y,TF
          REAL                G1,       G2
          REAL                PIVOTF,   P,        Q,        SIZEF,  T
          DOUBLE PRECISION    S11, S12, DOTF
          CHARACTER*6         HNAME
          IPAIRF(J,K)  =  J*2**12 + K
          PIVOTF(X)    =  ABS(SNGL(X))
          SIZEF(X)     =  ABS(SNGL(X))
          DOTF(X,Y,S11)  =  X * Y + S11
*# 31 "dfact.F"
          DATA      G1, G2              /  1.E-19,  1.E19  /




          DATA      HNAME               /  ' DFACT'  /
          DATA      ZERO, ONE           /  0.D0, 1.D0  /
          DATA      NORMAL, IMPOSS      /  0, -1  /
          DATA      JRANGE, JOVER, JUNDER  /  0, +1, -1  /

*# 1 "fact.inc" 1
*
* $Id: fact.inc,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: fact.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
*
* fact.inc
*
          IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 110
             CALL TMPRNT(HNAME,N,IDIM,0)
             RETURN
 110      IFAIL  =  NORMAL
          JFAIL  =  JRANGE
          NXCH   =  0
          DET    =  ONE
          DO 144    J  =  1, N
cmsh 120         K  =  J
             K  =  J
             P  =  PIVOTF(A(J,J))
             IF(J .EQ. N)  GOTO 122
             JP1  =  J+1
             DO 121    I  =  JP1, N
                Q  =  PIVOTF(A(I,J))
                IF(Q .LE. P)  GOTO 121
                   K  =  I
                   P  =  Q
 121            CONTINUE
             IF(K .NE. J)  GOTO 123
 122         IF(P .GT. 0.)  GOTO 130
                DET    =  ZERO
                IFAIL  =  IMPOSS
                JFAIL  =  JRANGE
                RETURN
 123         DO 124    L  =  1, N
                TF      =  A(J,L)
                A(J,L)  =  A(K,L)
                A(K,L)  =  TF
 124            CONTINUE
             NXCH      =  NXCH + 1
             IR(NXCH)  =  IPAIRF(J,K)
 130         DET     =  DET * A(J,J)
             A(J,J)  =  ONE / A(J,J)
             T  =  SIZEF(DET)
             IF(T .LT. G1)  THEN
                DET    =  ZERO
                IF(JFAIL .EQ. JRANGE)  JFAIL  =  JUNDER
             ELSEIF(T .GT. G2)  THEN
                DET    =  ONE
                IF(JFAIL .EQ. JRANGE)  JFAIL  =  JOVER
             ENDIF
             IF(J .EQ. N)  GOTO 144
             JM1  =  J-1
             JP1  =  J+1
             DO 143   K  =  JP1, N
                S11  =  -A(J,K)
                S12  =  -A(K,J+1)
                IF(J .EQ. 1)  GOTO 142
                DO 141  I  =  1, JM1
                   S11  =  DOTF(A(I,K),A(J,I),S11)
                   S12  =  DOTF(A(I,J+1),A(K,I),S12)
 141               CONTINUE
 142            A(J,K)    =  -S11 * A(J,J)
                A(K,J+1)  =  -DOTF(A(J,J+1),A(K,J),S12)
 143            CONTINUE
 144         CONTINUE
cmsh 150      IF(MOD(NXCH,2) .NE. 0)  DET  =  -DET
          IF(MOD(NXCH,2) .NE. 0)  DET  =  -DET
          IF(JFAIL .NE. JRANGE)   DET  =  ZERO
          IR(N)  =  NXCH
*# 41 "dfact.F" 2
          RETURN
          END
*CMZ :  1.16/04 17/04/2014  11.20.34  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*# 1 "dfeqn.F"
*# 1 "<command-line>"
*# 1 "dfeqn.F"
*
* $Id: dfeqn.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfeqn.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*

*# 1 "kernnum/pilot.h" 1
*# 21 "kernnum/pilot.h"

*# 33 "kernnum/pilot.h"

*# 10 "dfeqn.F" 2
          SUBROUTINE          DFEQN(N,A,IDIM,IR,K,B)
          INTEGER             IR(*)
          DOUBLE PRECISION    A(IDIM,*),B(IDIM,*),X,Y,TE
          DOUBLE PRECISION    S21, S22, DOTF
          CHARACTER*6         HNAME
          DOTF(X,Y,S21)  =  X*Y + S21
          DATA      HNAME               /  ' DFEQN'  /

*# 1 "feqn.inc" 1
*
* $Id: feqn.inc,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: feqn.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
*
* feqn.inc
*
          IF(IDIM .GE. N  .AND.  N .GT. 0  .AND.  K .GT. 0)  GOTO 210
          CALL TMPRNT(HNAME,N,IDIM,K)
          RETURN
 210      NXCH  =  IR(N)
          IF(NXCH .EQ. 0)  GOTO 220
          DO 212    M  =  1, NXCH
             IJ  =  IR(M)
             I   =  IJ / 4096
             J   =  MOD(IJ,4096)
             DO 211   L  =  1, K
                TE      =  B(I,L)
                B(I,L)  =  B(J,L)
                B(J,L)  =  TE
 211            CONTINUE
 212         CONTINUE
 220      DO 221    L  =  1, K
             B(1,L)  =  A(1,1)*B(1,L)
 221         CONTINUE
          IF(N .EQ. 1)  GOTO 299
          DO 243    L  =  1, K
             DO 232   I  =  2, N
                IM1  =  I-1
                S21  =  - B(I,L)
                DO 231   J  =  1, IM1
                   S21  =  DOTF(A(I,J),B(J,L),S21)
 231               CONTINUE
                B(I,L)  =  - A(I,I)*S21
 232            CONTINUE
             NM1  =  N-1
             DO 242   I  =  1, NM1
                NMI  =  N-I
                S22  =  - B(NMI,L)
                DO 241   J  =  1, I
                   NMJP1  =  N - J+1
                   S22    =  DOTF(A(NMI,NMJP1),B(NMJP1,L),S22)
 241               CONTINUE
                B(NMI,L)  =  - S22
 242            CONTINUE
 243         CONTINUE
 299      CONTINUE
*# 18 "dfeqn.F" 2
          RETURN
          END
*CMZ :          28/08/2014  14.01.08  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX dfinv.F

# 1 "dfinv.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "dfinv.F"
*
* $Id: dfinv.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfinv.F,v $
* Revision 1.1.1.1 1996/02/15 17:49:03 mclareni
* Kernlib
*
*
# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 10 "dfinv.F" 2
          SUBROUTINE DFINV(N,A,IDIM,IR)
          INTEGER IR(*)
          DOUBLE PRECISION A(IDIM,*),ZERO, X, Y, TI
          DOUBLE PRECISION S31, S32, S33, S34, DOTF
          CHARACTER*6 HNAME
          DATA HNAME / ' DFINV' /
          DOTF(X,Y,S31) = X*Y + S31
          DATA ZERO / 0.D0 /
# 1 "finv.inc" 1
*
* $Id: finv.inc,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: finv.inc,v $
* Revision 1.1.1.1 1996/02/15 17:49:03 mclareni
* Kernlib
*
*
*
* finv.inc
*
          IF(IDIM .GE. N .AND. N .GT. 0) GOTO 310
             CALL TMPRNT(HNAME,N,IDIM,0)
             RETURN
 310     IF(N .EQ. 1) RETURN
          A(2,1) = -A(2,2) * DOTF(A(1,1),A(2,1),ZERO)
          A(1,2) = -A(1,2)
          IF(N .EQ. 2) GOTO 330
          DO 314 I = 3, N
             IM2 = I-2
             DO 312 J = 1, IM2
                S31 = ZERO
                S32 = A(J,I)
                DO 311 K = J, IM2
                   S31 = DOTF(A(K,J),A(I,K),S31)
                   S32 = DOTF(A(J,K+1),A(K+1,I),S32)
 311      CONTINUE
                A(I,J) = -A(I,I) * DOTF(A(I-1,J),A(I,I-1),S31)
                A(J,I) = -S32
 312      CONTINUE
             A(I,I-1) = -A(I,I) * DOTF(A(I-1,I-1),A(I,I-1),ZERO)
             A(I-1,I) = -A(I-1,I)
 314      CONTINUE
 330      NM1 = N-1
          DO 335 I = 1, NM1
             NMI = N-I
             DO 332 J = 1, I
                S33 = A(I,J)
                DO 331 K = 1, NMI
                   S33 = DOTF(A(I+K,J),A(I,I+K),S33)
 331      CONTINUE
                A(I,J) = S33
 332      CONTINUE
             DO 334 J = 1, NMI
                S34 = ZERO
                DO 333 K = J, NMI
                   S34 = DOTF(A(I+K,I+J),A(I,I+K),S34)
 333      CONTINUE
                A(I,I+J) = S34
 334      CONTINUE
 335      CONTINUE
          NXCH = IR(N)
          IF(NXCH .EQ. 0) RETURN
            DO 342 M = 1, NXCH
             K = NXCH - M+1
             IJ = IR(K)
             I = IJ / 4096
             J = MOD(IJ,4096)
             DO 341 K = 1, N
                TI = A(K,I)
                A(K,I) = A(K,J)
                A(K,J) = TI
 341         CONTINUE
 342       CONTINUE
# 19 "dfinv.F" 2
          RETURN
          END
*CMZ :          02/05/2017  15.30.22  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
*CMZ :          02/05/2017  15.19.40  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
*CMZ :          02/05/2017  10.18.43  by  Michael Scheer
*-- Author :
      SUBROUTINE DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
*
*  -- LAPACK test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  This routine is deprecated and has been replaced by routine DGEQP3.
*
*  DGEQPF computes a QR factorization with column pivoting of a
*  real M-by-N matrix A: A*P = Q*R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the upper triangle of the array contains the
*          min(M,N)-by-N upper triangular matrix R; the elements
*          below the diagonal, together with the array TAU,
*          represent the orthogonal matrix Q as a product of
*          min(m,n) elementary reflectors.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
*          to the front of A*P (a leading column); if JPVT(i) = 0,
*          the i-th column of A is a free column.
*          On exit, if JPVT(i) = k, then the i-th column of A*P
*          was the k-th column of A.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(n)
*
*  Each H(i) has the form
*
*     H = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i).
*
*  The matrix P is represented in jpvt as follows: If
*     jpvt(j) = i
*  then the jth column of P is the ith canonical unit vector.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MA, MN, PVT
      DOUBLE PRECISION   AII, TEMP, TEMP2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARF, DLARFG, DORM2R, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DNRM2
      EXTERNAL           IDAMAX, DNRM2
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQPF', -INFO )
         RETURN
      END IF
*
      MN = MIN( M, N )
*
*     Move initial columns up front
*
      ITEMP = 1
      DO 10 I = 1, N
         IF( JPVT( I ).NE.0 ) THEN
            IF( I.NE.ITEMP ) THEN
               CALL DSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            END IF
            ITEMP = ITEMP + 1
         ELSE
            JPVT( I ) = I
         END IF
   10 CONTINUE
      ITEMP = ITEMP - 1
*
*     Compute the QR factorization and update remaining columns
*
      IF( ITEMP.GT.0 ) THEN
         MA = MIN( ITEMP, M )
         CALL DGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         IF( MA.LT.N ) THEN
            CALL DORM2R( 'Left', 'Transpose', M, N-MA, MA, A, LDA, TAU,
     $                   A( 1, MA+1 ), LDA, WORK, INFO )
         END IF
      END IF
*
      IF( ITEMP.LT.MN ) THEN
*
*        Initialize partial column norms. The first n elements of
*        work store the exact column norms.
*
         DO 20 I = ITEMP + 1, N
            WORK( I ) = DNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            WORK( N+I ) = WORK( I )
   20    CONTINUE
*
*        Compute factorization
*
         DO 40 I = ITEMP + 1, MN
*
*           Determine ith pivot column and swap if necessary
*
            PVT = ( I-1 ) + IDAMAX( N-I+1, WORK( I ), 1 )
*
            IF( PVT.NE.I ) THEN
               CALL DSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               WORK( PVT ) = WORK( I )
               WORK( N+PVT ) = WORK( N+I )
            END IF
*
*           Generate elementary reflector H(i)
*
            IF( I.LT.M ) THEN
               CALL DLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) )
            ELSE
               CALL DLARFG( 1, A( M, M ), A( M, M ), 1, TAU( M ) )
            END IF
*
            IF( I.LT.N ) THEN
*
*              Apply H(i) to A(i:m,i+1:n) from the left
*
               AII = A( I, I )
               A( I, I ) = ONE
               CALL DLARF( 'LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                     A( I, I+1 ), LDA, WORK( 2*N+1 ) )
               A( I, I ) = AII
            END IF
*
*           Update partial column norms
*
            DO 30 J = I + 1, N
               IF( WORK( J ).NE.ZERO ) THEN
                  TEMP = ONE - ( ABS( A( I, J ) ) / WORK( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = ONE + 0.05D0*TEMP*
     $                    ( WORK( J ) / WORK( N+J ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     IF( M-I.GT.0 ) THEN
                        WORK( J ) = DNRM2( M-I, A( I+1, J ), 1 )
                        WORK( N+J ) = WORK( J )
                     ELSE
                        WORK( J ) = ZERO
                        WORK( N+J ) = ZERO
                     END IF
                  ELSE
                     WORK( J ) = WORK( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE
*
   40    CONTINUE
      END IF
      RETURN
*
*     End of DGEQPF
*
      END
*CMZ :          02/05/2017  15.25.27  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEQR2 computes a QR factorization of a real m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END
*CMZ :          02/05/2017  15.38.33  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
*CMZ :          02/05/2017  15.35.59  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
*CMZ :          02/05/2017  15.37.15  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
*CMZ :          02/05/2017  15.27.39  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFB applies a real block reflector H or its transpose H' to a
*  real m by n matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'T': apply H' (Transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) DOUBLE PRECISION array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
*          The triangular k by k matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of DLARFB
*
      END
*CMZ :          02/05/2017  15.24.24  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO,
     $                  WORK, 1 )
*
*           C := C - v * w'
*
            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of DLARF
*
      END
*CMZ :          02/05/2017  15.26.35  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of DLARFG
*
      END
*CMZ :          02/05/2017  15.18.34  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFT forms the triangular factor T of a real block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) DOUBLE PRECISION array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
*
                  CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ),
     $                        V( I, 1 ), LDV, V( I, I ), 1, ZERO,
     $                        T( 1, I ), 1 )
               ELSE
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
*
                  CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
*
                     CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ),
     $                           V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO,
     $                           T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
*
                     CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of DLARFT
*
      END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "dmadd.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmadd.F"
*
* $Id: dmadd.F,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dmadd.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmadd.F" 2
          SUBROUTINE          DMADD(M,N,X,X12,X21,Y,Y12,Y21,Z,Z12,Z21)
          DOUBLE PRECISION    X(*), X12(*), X21(*), Y(*), Y12(*), Y21(*)
          DOUBLE PRECISION    Z(*), Z12(*), Z21(*), ADD,  A,      B
          ADD(A,B)  =  A+B
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 16 "dmadd.F" 2

# 1 "dyij.inc" 1
*
* $Id: dyij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyij.inc
*

          IY  =  (LOCF(Y21) - LOCF(Y)) / 2
          JY  =  (LOCF(Y12) - LOCF(Y)) / 2
# 17 "dmadd.F" 2

# 1 "dzij.inc" 1
*
* $Id: dzij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzij.inc
*

          IZ  =  (LOCF(Z21) - LOCF(Z)) / 2
          JZ  =  (LOCF(Z12) - LOCF(Z)) / 2
# 18 "dmadd.F" 2

# 1 "madd.inc" 1
*
* $Id: madd.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: madd.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* madd.inc
*
          MM  =  M
          NN  =  N
          IF(MM .GT. NN) THEN
             MN  =  NN
             NN  =  MM
             MM  =  MN
             IJ  =  JX
             JX  =  IX
             IX  =  IJ
             IJ  =  JY
             JY  =  IY
             IY  =  IJ
             IJ  =  JZ
             JZ  =  IZ
             IZ  =  IJ
          ENDIF
          LXI1  =  1
          LYI1  =  1
          LZI1  =  1
          DO 12     I  =  1, MM
             LXIJ  =  LXI1
             LYIJ  =  LYI1
             LZIJ  =  LZI1
             DO 11  J  =  1, NN
                Z(LZIJ)  =  ADD( X(LXIJ),Y(LYIJ) )
                LXIJ     =  LXIJ + JX
                LYIJ     =  LYIJ + JY
                LZIJ     =  LZIJ + JZ
  11            CONTINUE
             LXI1  =  LXI1 + IX
             LYI1  =  LYI1 + IY
             LZI1  =  LZI1 + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 19 "dmadd.F" 2

# 1 "dmbil.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmbil.F"
*
* $Id: dmbil.F,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dmbil.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmbil.F" 2
          DOUBLE PRECISION FUNCTION DMBIL(N,X,X2,Y,Y12,Y21,Z,Z2)
          DOUBLE PRECISION X(*),X2(*),Y(*),Y12(*),Y21(*),Z(*),Z2(*)
          DOUBLE PRECISION A, B, SUM, ZERO, F, G, SXYZ, SYZ
          F(A,B,SUM)  =  A*B + SUM
          G(A,B,SUM)  =  A*B + SUM
          DATA      ZERO      /  0.D0  /
          SXYZ  =  ZERO
          IF(N .LE. 0)  GOTO 20

# 1 "dxi.inc" 1
*
* $Id: dxi.inc,v 1.1.1.1 1996/02/15 17:48:54 mclareni Exp $
*
* $Log: dxi.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:54  mclareni
* Kernlib
*
*
*
* dxi.inc
*

          IX  =  (LOCF(X2)  - LOCF(X)) / 2
# 19 "dmbil.F" 2

# 1 "dyij.inc" 1
*
* $Id: dyij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyij.inc
*

          IY  =  (LOCF(Y21) - LOCF(Y)) / 2
          JY  =  (LOCF(Y12) - LOCF(Y)) / 2
# 20 "dmbil.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 21 "dmbil.F" 2

# 1 "mbil.inc" 1
*
* $Id: mbil.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mbil.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mbil.inc
*
          LXI  =  1
          LYI  =  1
          DO 12     I  =  1, N
             SYZ   =  ZERO
             LYIJ  =  LYI
             LZJ   =  1
             DO 11  J  =  1, N
                SYZ   =  F(Y(LYIJ),Z(LZJ),SYZ)
                LYIJ  =  LYIJ + JY
                LZJ   =  LZJ + JZ
  11            CONTINUE
             SXYZ  =  G(SYZ,X(LXI),SXYZ)
             LXI   =  LXI + IX
             LYI   =  LYI + IY
  12         CONTINUE
# 22 "dmbil.F" 2
  20      DMBIL  =  SXYZ
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "dmcpy.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmcpy.F"
*
* $Id: dmcpy.F,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dmcpy.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmcpy.F" 2
          SUBROUTINE          DMCPY(M,N,X,X12,X21,Z,Z12,Z21)
          DOUBLE PRECISION    X(*),X12(*),X21(*),Z(*),Z12(*),Z21(*)
          DOUBLE PRECISION    FUNCT, A
          FUNCT(A)  =  A
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 16 "dmcpy.F" 2

# 1 "dzij.inc" 1
*
* $Id: dzij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzij.inc
*

          IZ  =  (LOCF(Z21) - LOCF(Z)) / 2
          JZ  =  (LOCF(Z12) - LOCF(Z)) / 2
# 17 "dmcpy.F" 2

# 1 "mcpy.inc" 1
*
* $Id: mcpy.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mcpy.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mcpy.inc
*
          MM  =  M
          NN  =  N
          IF(MM .GT. NN)  THEN
             MN  =  NN
             NN  =  MM
             MM  =  MN
             IJ  =  JX
             JX  =  IX
             IX  =  IJ
             IJ  =  JZ
             JZ  =  IZ
             IZ  =  IJ
          ENDIF
          LXI1  =  1
          LZI1  =  1
          DO 12     I  =  1, MM
             LXIJ  =  LXI1
             LZIJ  =  LZI1
             DO 11     J  =  1, NN
                Z(LZIJ)  =  FUNCT( X(LXIJ) )
                LXIJ  =  LXIJ + JX
                LZIJ  =  LZIJ + JZ
  11         CONTINUE
             LXI1  =  LXI1 + IX
             LZI1  =  LZI1 + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  15.06.48  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 1 "dmmlt.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "dmmlt.F"
*
* $Id: dmmlt.F,v 1.1.1.1 1996/02/15 17:49:01 mclareni Exp $
*
* $Log: dmmlt.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:01  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "dmmlt.F" 2
          SUBROUTINE        DMMLT(M,N,K,X,X12,X21,Y,Y12,Y21,Z,Z12,Z21,T)
          DOUBLE PRECISION    X(*), X12(*), X21(*), Y(*), Y12(*), Y21(*)
          DOUBLE PRECISION    Z(*), Z12(*), Z21(*), T(*), A, B, CNJF
          DOUBLE PRECISION    ZERO, SUM, DOTF, SQRF
          DOUBLE PRECISION    S11, S21, S22, S31, S41, S51, S52
          DATA      ZERO      /  0.D0  /
          DOTF(A,B,SUM)  =  A*B + SUM
          SQRF(A,SUM)    =  A*A + SUM
          CNJF(A)        =  A

# 1 "dlocf.inc" 1
*
* $Id: dlocf.inc,v 1.1.1.1 1996/02/15 17:49:00 mclareni Exp $
*
* $Log: dlocf.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:00  mclareni
* Kernlib
*
*
*
* dlocf.inc
*

          IF(MIN0(M,N,K) .LE. 0)  RETURN
          LOCX  =  LOCF(X(1))
          LOCY  =  LOCF(Y(1))
          LOCZ  =  LOCF(Z(1))
          IX  =  (LOCF(X21(1)) - LOCX) / 2
          JX  =  (LOCF(X12(1)) - LOCX) / 2
          JY  =  (LOCF(Y21(1)) - LOCY) / 2
          LY  =  (LOCF(Y12(1)) - LOCY) / 2
          IZ  =  (LOCF(Z21(1)) - LOCZ) / 2
          LZ  =  (LOCF(Z12(1)) - LOCZ) / 2
cmsh # 20 "dmmlt.F" 2

# 1 "mmlt.inc" 1
*
* $Id: mmlt.inc,v 1.1.1.1 1996/02/15 17:49:01 mclareni Exp $
*
* $Log: mmlt.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:01  mclareni
* Kernlib
*
*
*
* mmlt.inc
*

# 1 "zisxy.inc" 1
*
* $Id: zisxy.inc,v 1.1.1.1 1996/02/15 17:49:01 mclareni Exp $
*
* $Log: zisxy.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:01  mclareni
* Kernlib
*
*
*
* zisxy.inc
*
          IF(LOCZ .EQ. LOCX)  GOTO 30
          IF(LOCZ .EQ. LOCY)  GOTO 40
          IF(LOCX .EQ. LOCY)  GOTO 20
  10      LY1L  =  1
          LZ1L  =  1
          DO 13     L  =  1, K
             LXI1  =  1
             LZIL  =  LZ1L
             DO 12  I  =  1, M
                S11   =  ZERO
                LXIJ  =  LXI1
                LYJL  =  LY1L
                DO 11  J  =  1, N
                   S11   =  DOTF(X(LXIJ),Y(LYJL),S11)
                   LXIJ  =  LXIJ + JX
                   LYJL  =  LYJL + JY
  11               CONTINUE
                Z(LZIL)  =  S11
                LXI1     =  LXI1 + IX
                LZIL     =  LZIL + IZ
  12            CONTINUE
             LY1L  =  LY1L + LY
             LZ1L  =  LZ1L + LZ
  13         CONTINUE
          RETURN
  20      IF(M .NE. K  .OR.  IX .NE. LY  .OR.  JX .NE. JY)  GOTO 10
          LXI1  =  1
          LZII  =  1
          DO 24     I  =  1, M
             S21   =  ZERO
             LXIJ  =  LXI1
             DO 21  J  =  1, N
                S21   =  SQRF(X(LXIJ),S21)
                LXIJ  =  LXIJ + JX
  21            CONTINUE
             Z(LZII)  =  S21
             IF(I .EQ. M)  GOTO 24
             LXK1  =  LXI1 + IX
             LZIK  =  LZII + LZ
             LZKI  =  LZII + IZ
             DO 23  KDASH  =  I+1, M
                S22   =  ZERO
                LXIJ  =  LXI1
                LXKJ  =  LXK1
                DO 22  J  =  1, N
                   S22   =  DOTF(X(LXIJ),X(LXKJ),S22)
                   LXIJ  =  LXIJ + JX
                   LXKJ  =  LXKJ + JX
  22               CONTINUE
                Z(LZIK)  =  S22
                Z(LZKI)  =  CNJF( Z(LZIK) )
                LXK1  =  LXK1 + IX
                LZIK  =  LZIK + LZ
                LZKI  =  LZKI + IZ
  23            CONTINUE
             LXI1  =  LXI1 + IX
             LZII  =  LZII + IZ + LZ
  24         CONTINUE
          RETURN
  30      IF(LOCX .EQ. LOCY)  GOTO 50
          LXI1  =  1
          DO 34     I  =  1, M
             LY1L  =  1
             LTL   =  1
             DO 32  L  =  1, K
                S31   =  ZERO
                LXIJ  =  LXI1
                LYJL  =  LY1L
                DO 31  J  =  1, N
                   S31   =  DOTF(X(LXIJ),Y(LYJL),S31)
                   LXIJ  =  LXIJ + JX
                   LYJL  =  LYJL + JY
  31               CONTINUE
                T(LTL)  =  S31
                LY1L    =  LY1L + LY
                LTL     =  LTL + 1
  32            CONTINUE
             LXIL  =  LXI1
             LTL   =  1
             DO 33  L  =  1, K
                X(LXIL)  =  T(LTL)
                LXIL     =  LXIL + JX
                LTL      =  LTL + 1
  33            CONTINUE
             LXI1  =  LXI1 + IX
  34         CONTINUE
          RETURN
  40      LY1L  =  1
          DO 44     L  =  1, K
             LXI1  =  1
             LTI   =  1
             DO 42  I  =  1, M
                S41   =  ZERO
                LXIJ  =  LXI1
                LYJL  =  LY1L
                DO 41  J  =  1, N
                   S41   =  DOTF(X(LXIJ),Y(LYJL),S41)
                   LXIJ  =  LXIJ + JX
                   LYJL  =  LYJL + JY
  41               CONTINUE
                T(LTI)  =  S41
                LXI1    =  LXI1 + IX
                LTI     =  LTI + 1
  42            CONTINUE
             LYIL  =  LY1L
             LTI   =  1
             DO 43  I  =  1, M
                Y(LYIL)  =  T(LTI)
                LYIL     =  LYIL + JY
                LTI      =  LTI + 1
  43            CONTINUE
             LY1L  =  LY1L + LY
  44         CONTINUE
          RETURN
# 13 "mmlt.inc" 2

# 1 "xisxxtra.inc" 1
*
* $Id: xisxxtra.inc,v 1.1.1.1 1996/02/15 17:49:01 mclareni Exp $
*
* $Log: xisxxtra.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:01  mclareni
* Kernlib
*
*
*
* xisxxtra.inc
*
  50      LXI1  =  1
          LXII  =  1
          DO 56     I  =  1, M
             S51   =  ZERO
             LXIJ  =  LXI1
             DO 51  J  =  1, N
                S51   =  SQRF(X(LXIJ),S51)
                LXIJ  =  LXIJ + JX
  51            CONTINUE
             T(1)  =  S51
             IF(I .EQ. M)  GOTO 54
             LXK1  =  LXI1 + IX
             LTK  =  2
             DO 53  KDASH  =  I+1, M
                S52   =  ZERO
                LXIJ  =  LXI1
                LXKJ  =  LXK1
                DO 52  J  =  1, N
                   S52   =  DOTF(X(LXIJ),X(LXKJ),S52)
                   LXIJ  =  LXIJ + JX
                   LXKJ  =  LXKJ + JX
  52               CONTINUE
                T(LTK)  =  S52
                LXK1    =  LXK1 + IX
                LTK     =  LTK + 1
  53            CONTINUE
  54         LXIK  =  LXII
             LTK   =  1
             DO 55  KDASH  =  I, M
                X(LXIK)  =  T(LTK)
                LXIK     =  LXIK + JX
                LTK      =  LTK + 1
  55            CONTINUE
             LXI1     =  LXI1 + IX
             LXII     =  LXII + IX + JX
  56         CONTINUE
          IF(M .EQ. 1)  RETURN
          LXII  =  1
          DO 58     I  =  1, M-1
             LXIK  =  LXII + JX
             LXKI  =  LXII + IX
             DO 57  KDASH  =  I+1, M
                X(LXKI)  =  CNJF( X(LXIK) )
                LXIK     =  LXIK + JX
                LXKI     =  LXKI + IX
  57            CONTINUE
             LXII  =  LXII + IX + JX
  58         CONTINUE
          RETURN
cmsh # 14 "mmlt.inc" 2
cmsh # 21 "dmmlt.F" 2
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 18 "dmcpy.F" 2

# 1 "dmmna.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmmna.F"
*
* $Id: dmmna.F,v 1.1.1.1 1996/02/15 17:48:57 mclareni Exp $
*
* $Log: dmmna.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:57  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmmna.F" 2
          SUBROUTINE          DMMNA(M,N,X,X12,X21,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*),X12(*),X21(*),Y(*),Y2(*),Z(*),Z2(*)
          DOUBLE PRECISION    A, B, SUM, F, SIGNF
          F(A,B,SUM)  =  -A*B + SUM
          SIGNF(A)    =  A
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 17 "dmmna.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dmmna.F" 2

# 1 "dzi.inc" 1
*
* $Id: dzi.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzi.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzi.inc
*

          IZ  =  (LOCF(Z2)  - LOCF(Z)) / 2
# 19 "dmmna.F" 2

# 1 "mmpa.inc" 1
*
* $Id: mmpa.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mmpa.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mmpa.inc
*
          LXI1  =  1
          LZI   =  1
          DO 12     I  =  1, M
             LXIJ  =  LXI1
             LYJ   =  1
             SUM   =  SIGNF( Z(LZI) )
             DO 11  J  =  1, N
                SUM  =  F(X(LXIJ),Y(LYJ),SUM)
                LXIJ =  LXIJ + JX
                LYJ  =  LYJ + JY
  11            CONTINUE
             Z(LZI)  =  SUM
             LXI1    =  LXI1 + IX
             LZI     =  LZI + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 20 "dmmna.F" 2

# 1 "dmmns.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmmns.F"
*
* $Id: dmmns.F,v 1.1.1.1 1996/02/15 17:48:57 mclareni Exp $
*
* $Log: dmmns.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:57  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmmns.F" 2
          SUBROUTINE          DMMNS(M,N,X,X12,X21,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*),X12(*),X21(*),Y(*),Y2(*),Z(*),Z2(*)
          DOUBLE PRECISION    A, B, SUM, F, SIGNF
          F(A,B,SUM)  =  -A*B + SUM
          SIGNF(A)    =  -A
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 17 "dmmns.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dmmns.F" 2

# 1 "dzi.inc" 1
*
* $Id: dzi.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzi.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzi.inc
*

          IZ  =  (LOCF(Z2)  - LOCF(Z)) / 2
# 19 "dmmns.F" 2

# 1 "mmpa.inc" 1
*
* $Id: mmpa.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mmpa.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mmpa.inc
*
          LXI1  =  1
          LZI   =  1
          DO 12     I  =  1, M
             LXIJ  =  LXI1
             LYJ   =  1
             SUM   =  SIGNF( Z(LZI) )
             DO 11  J  =  1, N
                SUM  =  F(X(LXIJ),Y(LYJ),SUM)
                LXIJ =  LXIJ + JX
                LYJ  =  LYJ + JY
  11            CONTINUE
             Z(LZI)  =  SUM
             LXI1    =  LXI1 + IX
             LZI     =  LZI + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 20 "dmmns.F" 2

# 1 "dmmpa.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmmpa.F"
*
* $Id: dmmpa.F,v 1.1.1.1 1996/02/15 17:48:57 mclareni Exp $
*
* $Log: dmmpa.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:57  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmmpa.F" 2
          SUBROUTINE          DMMPA(M,N,X,X12,X21,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*),X12(*),X21(*),Y(*),Y2(*),Z(*),Z2(*)
          DOUBLE PRECISION    A, B, SUM, F, SIGNF
          F(A,B,SUM)  =  A*B + SUM
          SIGNF(A)    =  A
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 17 "dmmpa.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dmmpa.F" 2

# 1 "dzi.inc" 1
*
* $Id: dzi.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzi.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzi.inc
*

          IZ  =  (LOCF(Z2)  - LOCF(Z)) / 2
# 19 "dmmpa.F" 2

# 1 "mmpa.inc" 1
*
* $Id: mmpa.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mmpa.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mmpa.inc
*
          LXI1  =  1
          LZI   =  1
          DO 12     I  =  1, M
             LXIJ  =  LXI1
             LYJ   =  1
             SUM   =  SIGNF( Z(LZI) )
             DO 11  J  =  1, N
                SUM  =  F(X(LXIJ),Y(LYJ),SUM)
                LXIJ =  LXIJ + JX
                LYJ  =  LYJ + JY
  11            CONTINUE
             Z(LZI)  =  SUM
             LXI1    =  LXI1 + IX
             LZI     =  LZI + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 20 "dmmpa.F" 2

# 1 "dmmps.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmmps.F"
*
* $Id: dmmps.F,v 1.1.1.1 1996/02/15 17:48:57 mclareni Exp $
*
* $Log: dmmps.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:57  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmmps.F" 2
          SUBROUTINE          DMMPS(M,N,X,X12,X21,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*),X12(*),X21(*),Y(*),Y2(*),Z(*),Z2(*)
          DOUBLE PRECISION    A, B, SUM, F, SIGNF
          F(A,B,SUM)  =  A*B + SUM
          SIGNF(A)    =  -A
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 17 "dmmps.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dmmps.F" 2

# 1 "dzi.inc" 1
*
* $Id: dzi.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzi.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzi.inc
*

          IZ  =  (LOCF(Z2)  - LOCF(Z)) / 2
# 19 "dmmps.F" 2

# 1 "mmpa.inc" 1
*
* $Id: mmpa.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mmpa.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mmpa.inc
*
          LXI1  =  1
          LZI   =  1
          DO 12     I  =  1, M
             LXIJ  =  LXI1
             LYJ   =  1
             SUM   =  SIGNF( Z(LZI) )
             DO 11  J  =  1, N
                SUM  =  F(X(LXIJ),Y(LYJ),SUM)
                LXIJ =  LXIJ + JX
                LYJ  =  LYJ + JY
  11            CONTINUE
             Z(LZI)  =  SUM
             LXI1    =  LXI1 + IX
             LZI     =  LZI + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 20 "dmmps.F" 2

# 1 "dmmpy.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmmpy.F"
*
* $Id: dmmpy.F,v 1.1.1.1 1996/02/15 17:48:58 mclareni Exp $
*
* $Log: dmmpy.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:58  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmmpy.F" 2
          SUBROUTINE          DMMPY(M,N,X,X12,X21,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*),X12(*),X21(*),Y(*),Y2(*),Z(*),Z2(*)
          DOUBLE PRECISION    A, B, SUM, ZERO, F
          F(A,B,SUM)  =  A*B + SUM
          DATA ZERO    / 0.D0 /
          IF(M .LE. 0  .OR.  N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 17 "dmmpy.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dmmpy.F" 2

# 1 "dzi.inc" 1
*
* $Id: dzi.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzi.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzi.inc
*

          IZ  =  (LOCF(Z2)  - LOCF(Z)) / 2
# 19 "dmmpy.F" 2

# 1 "mmpy.inc" 1
*
* $Id: mmpy.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mmpy.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mmpy.inc
*
          LXI1  =  1
          LZI   =  1
          DO 12     I  =  1, M
             LXIJ  =  LXI1
             LYJ   =  1
             SUM   =  ZERO
             DO 11  J  =  1, N
                SUM  =  F(X(LXIJ),Y(LYJ),SUM)
                LXIJ =  LXIJ + JX
                LYJ  =  LYJ + JY
  11            CONTINUE
             Z(LZI)  =  SUM
             LXI1    =  LXI1 + IX
             LZI     =  LZI + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 20 "dmmpy.F" 2

# 1 "dmran.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmran.F"
*
* $Id: dmran.F,v 1.1.1.1 1996/02/15 17:48:58 mclareni Exp $
*
* $Log: dmran.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:58  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmran.F" 2
          SUBROUTINE          DMRAN(M,N,A,B,Z,Z12,Z21)
          DOUBLE PRECISION    A, B, Z(*), Z12(*), Z21(*), C
          DOUBLE PRECISION    DRANF
          IF(M .LE. 0  .OR. N .LE. 0)  RETURN

# 1 "dzij.inc" 1
*
* $Id: dzij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzij.inc
*

          IZ  =  (LOCF(Z21) - LOCF(Z)) / 2
          JZ  =  (LOCF(Z12) - LOCF(Z)) / 2
# 15 "dmran.F" 2
          MM  =  M
          NN  =  N
          IF(MM .GT. NN)  THEN
             MN  =  NN
             NN  =  MM
             MM  =  MN
             IJ  =  JZ
             JZ  =  IZ
             IZ  =  IJ
          ENDIF
          C     =  B - A
          LZI1  =  1
          DO 12     I  =  1, MM
             LZIJ  =  LZI1
             DO 11  J  =  1, NN
                Z(LZIJ)  =  C * DRANF() + A
                LZIJ     =  LZIJ + JZ
  11            CONTINUE
             LZI1  =  LZI1 + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "dmscl.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmscl.F"
*
* $Id: dmscl.F,v 1.1.1.1 1996/02/15 17:48:58 mclareni Exp $
*
* $Log: dmscl.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:58  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmscl.F" 2
          SUBROUTINE          DMSCL(M,N,S,X,X12,X21,Z,Z12,Z21)
          DOUBLE PRECISION    S, X(*),X12(*),X21(*), Z(*),Z12(*),Z21(*)
          DOUBLE PRECISION    FUNCT, A
          FUNCT(A)  =  S*A
          IF(M .LE. 0  .OR. N .LE. 0)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 16 "dmscl.F" 2

# 1 "dzij.inc" 1
*
* $Id: dzij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzij.inc
*

          IZ  =  (LOCF(Z21) - LOCF(Z)) / 2
          JZ  =  (LOCF(Z12) - LOCF(Z)) / 2
# 17 "dmscl.F" 2

# 1 "mcpy.inc" 1
*
* $Id: mcpy.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mcpy.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mcpy.inc
*
          MM  =  M
          NN  =  N
          IF(MM .GT. NN)  THEN
             MN  =  NN
             NN  =  MM
             MM  =  MN
             IJ  =  JX
             JX  =  IX
             IX  =  IJ
             IJ  =  JZ
             JZ  =  IZ
             IZ  =  IJ
          ENDIF
          LXI1  =  1
          LZI1  =  1
          DO 12     I  =  1, MM
             LXIJ  =  LXI1
             LZIJ  =  LZI1
             DO 11     J  =  1, NN
                Z(LZIJ)  =  FUNCT( X(LXIJ) )
                LXIJ  =  LXIJ + JX
                LZIJ  =  LZIJ + JZ
  11         CONTINUE
             LXI1  =  LXI1 + IX
             LZI1  =  LZI1 + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 18 "dmscl.F" 2

# 1 "dmset.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmset.F"
*
* $Id: dmset.F,v 1.1.1.1 1996/02/15 17:48:58 mclareni Exp $
*
* $Log: dmset.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:58  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmset.F" 2
          SUBROUTINE          DMSET(M,N,S,Z,Z12,Z21)
          DOUBLE PRECISION    S, Z(*), Z12(*), Z21(*)
          IF(M .LE. 0  .OR. N .LE. 0)  RETURN

# 1 "dzij.inc" 1
*
* $Id: dzij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzij.inc
*

          IZ  =  (LOCF(Z21) - LOCF(Z)) / 2
          JZ  =  (LOCF(Z12) - LOCF(Z)) / 2
# 14 "dmset.F" 2

# 1 "mset.inc" 1
*
* $Id: mset.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mset.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mset.inc
*
          MM  =  M
          NN  =  N
          IF(MM .GT. NN)  THEN
             MN  =  NN
             NN  =  MM
             MM  =  MN
             IJ  =  JZ
             JZ  =  IZ
             IZ  =  IJ
          ENDIF
          LZI1  =  1
          DO 12     I  =  1, MM
             LZIJ  =  LZI1
             DO 11  J  =  1, NN
                Z(LZIJ)  =  S
                LZIJ     =  LZIJ + JZ
  11            CONTINUE
             LZI1  =  LZI1 + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 15 "dmset.F" 2

# 1 "dmsub.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmsub.F"
*
* $Id: dmsub.F,v 1.1.1.1 1996/02/15 17:48:59 mclareni Exp $
*
* $Log: dmsub.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:59  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmsub.F" 2
          SUBROUTINE          DMSUB(M,N,X,X12,X21,Y,Y12,Y21,Z,Z12,Z21)
          DOUBLE PRECISION    X(*), X12(*), X21(*), Y(*), Y12(*), Y21(*)
          DOUBLE PRECISION    Z(*), Z12(*), Z21(*), ADD,  A,      B
          ADD(A,B)  =  A-B

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 15 "dmsub.F" 2

# 1 "dyij.inc" 1
*
* $Id: dyij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyij.inc
*

          IY  =  (LOCF(Y21) - LOCF(Y)) / 2
          JY  =  (LOCF(Y12) - LOCF(Y)) / 2
# 16 "dmsub.F" 2

# 1 "dzij.inc" 1
*
* $Id: dzij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dzij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dzij.inc
*

          IZ  =  (LOCF(Z21) - LOCF(Z)) / 2
          JZ  =  (LOCF(Z12) - LOCF(Z)) / 2
# 17 "dmsub.F" 2

# 1 "madd.inc" 1
*
* $Id: madd.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: madd.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* madd.inc
*
          MM  =  M
          NN  =  N
          IF(MM .GT. NN) THEN
             MN  =  NN
             NN  =  MM
             MM  =  MN
             IJ  =  JX
             JX  =  IX
             IX  =  IJ
             IJ  =  JY
             JY  =  IY
             IY  =  IJ
             IJ  =  JZ
             JZ  =  IZ
             IZ  =  IJ
          ENDIF
          LXI1  =  1
          LYI1  =  1
          LZI1  =  1
          DO 12     I  =  1, MM
             LXIJ  =  LXI1
             LYIJ  =  LYI1
             LZIJ  =  LZI1
             DO 11  J  =  1, NN
                Z(LZIJ)  =  ADD( X(LXIJ),Y(LYIJ) )
                LXIJ     =  LXIJ + JX
                LYIJ     =  LYIJ + JY
                LZIJ     =  LZIJ + JZ
  11            CONTINUE
             LXI1  =  LXI1 + IX
             LYI1  =  LYI1 + IY
             LZI1  =  LZI1 + IZ
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 18 "dmsub.F" 2

# 1 "dmutl.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dmutl.F"
*
* $Id: dmutl.F,v 1.1.1.1 1996/02/15 17:48:59 mclareni Exp $
*
* $Log: dmutl.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:59  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dmutl.F" 2
          SUBROUTINE          DMUTL(N,X,X12,X21)
          DOUBLE PRECISION    X(*), X12(*), X21(*)
          IF(N .LE. 1)  RETURN

# 1 "dxij.inc" 1
*
* $Id: dxij.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dxij.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dxij.inc
*

          IX  =  (LOCF(X21) - LOCF(X)) / 2
          JX  =  (LOCF(X12) - LOCF(X)) / 2
# 14 "dmutl.F" 2

# 1 "mutl.inc" 1
*
* $Id: mutl.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: mutl.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* mutl.inc
*
          LXII  =  1
          DO 12     IP1  =  2, N
             LXIJ  =  LXII
             LXJI  =  LXII
             DO 11  J  =  IP1, N
                LXIJ  =  LXIJ + JX
                LXJI  =  LXJI + IX
                X(LXJI)  =  X(LXIJ)
  11            CONTINUE
             LXII  =  LXII + IX + JX
  12         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  15.28.59  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
*     ..
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
*CMZ :          02/05/2017  15.23.28  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORM2R overwrites the general real m by n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'T', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'T',
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'T': apply Q' (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORM2R
*
      END
*CMZ :          02/05/2017  10.18.43  by  Michael Scheer
*-- Author :
# 15 "dmutl.F" 2

      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORMQR overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,
     $                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMQR
*
      END
*CMZ :          10/04/2019  12.04.16  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
cmsh # 1 "ranf.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "ranf.F"
*
* $Id: ranf.F,v 1.1.1.1 1996/02/15 17:49:05 mclareni Exp $
*
* $Log: ranf.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:05  mclareni
* Kernlib
*
*

cmsh # 1 "/usr/include/kernnum/pilot.h" 1 3 4
cmsh # 21 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 33 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "ranf.F" 2
cmsh          REAL FUNCTION RANF()
          DOUBLE PRECISION FUNCTION DRANF()
CMSH          DOUBLE PRECISION    DRANF,    G900GT,   G900ST
          DOUBLE PRECISION    G900GT,   G900ST
          DOUBLE PRECISION    DS(2),    DM(2),    DSEED
          DOUBLE PRECISION    DX24,     DX48
          DOUBLE PRECISION    DL,       DC,       DU,       DR
          LOGICAL             SINGLE
          DATA      DS     /  1665 1885.D0, 286 8876.D0  /
          DATA      DM     /  1518 4245.D0, 265 1554.D0  /
          DATA      DX24   /  1677 7216.D0  /
          DATA      DX48   /  281 4749 7671 0656.D0  /
CMSH          SINGLE  =  .TRUE.
CMSH          GOTO 10
CMSH          ENTRY DRANF()
          SINGLE  =  .FALSE.
  10      DL  =  DS(1) * DM(1)
          DC  =  DINT(DL/DX24)
          DL  =  DL - DC*DX24
          DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
          DS(2)  =  DU - DINT(DU/DX24)*DX24
          DS(1)  =  DL
          DR     =  (DS(2)*DX24 + DS(1)) / DX48
CMSH          IF(SINGLE)  THEN
CMSH             RANF  =  SNGL(DR)
CMSH          ELSE
             DRANF  =  DR
CMSH          ENDIF
          RETURN
cmsh          ENTRY G900GT()
cmsh          G900GT  =  DS(2)*DX24 + DS(1)
cmsh          RETURN
cmsh          ENTRY G900ST(DSEED)
cmsh          DS(2)  =  DINT(DSEED/DX24)
cmsh          DS(1)  =  DSEED - DS(2)*DX24
cmsh          G900ST =  DS(1)
cmsh          RETURN
          END
*CMZ :          02/05/2017  15.34.04  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
*CMZ :          02/05/2017  14.53.25  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 15 "dmutl.F" 2

cmsh # 1 "dsinv.F"
cmsh # 1 "<built-in>"
cmsh # 1 "<command-line>"
cmsh # 1 "dsinv.F"
*
* $Id: dsinv.F,v 1.2 1999/09/08 08:05:11 mclareni Exp $
*
* $Log: dsinv.F,v $
* Revision 1.2  1999/09/08 08:05:11  mclareni
* A problem was reported in DSINV which failed on very small numbers, probably
* due to converting to single before a test. The conversion has been removed her
* and also in DSFACT. This resulted in mods to sfact.inc and sfactd.inc which
* meant that some other routines had to be tidied also.
*
* Revision 1.1.1.1  1996/02/15 17:49:05  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 16 "dsinv.F" 2
          SUBROUTINE          DSINV(N,A,IDIM,IFAIL)
          DOUBLE PRECISION    A(IDIM,*),  ZERO,  ONE,  X, Y
          CHARACTER*6         HNAME
          DOUBLE PRECISION    S1, S31, S32, S33,  DOTF
          DOTF(X,Y,S1)  =  X * Y + S1
          DATA      HNAME               /  'DSINV '  /
          DATA      ZERO, ONE           /  0.D0, 1.D0 /
          IF(IDIM .LT. N  .OR.  N .LE. 0)  GOTO 900

# 1 "sfact.inc" 1
*
* $Id: sfact.inc,v 1.2 1999/09/08 08:05:21 mclareni Exp $
*
* $Log: sfact.inc,v $
* Revision 1.2  1999/09/08 08:05:21  mclareni
* A problem was reported in DSINV which failed on very small numbers, probably
* due to converting to single before a test. The conversion has been removed her
* and also in DSFACT. This resulted in mods to sfact.inc and sfactd.inc which
* meant that some other routines had to be tidied also.
*
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*
*
* sfact.inc
*
          IFAIL  =  0
          DO 144    J  =  1, N
             IF((A(J,J)) .LE. ZERO)  GOTO 150
             A(J,J)  =  ONE / A(J,J)
             IF(J .EQ. N)  GOTO 199
 140         JP1  =  J+1
             DO 143   L  =  JP1, N
                A(J,L)  =  A(J,J)*A(L,J)
                S1      =  -A(L,J+1)
                DO 141  I  =  1, J
                   S1  =  DOTF(A(L,I),A(I,J+1),S1)
 141               CONTINUE
                A(L,J+1)  =  -S1
 143            CONTINUE
 144         CONTINUE
 150      IFAIL  =  -1
          RETURN
 199      CONTINUE
cmsh # 25 "dsinv.F" 2

# 1 "sfinv.inc" 1
*
* $Id: sfinv.inc,v 1.1.1.1 1996/02/15 17:49:04 mclareni Exp $
*
* $Log: sfinv.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*
*
* sfinv.inc
*
          IF(N .EQ. 1)  GOTO 399
          A(1,2)  =  -A(1,2)
          A(2,1)  =   A(1,2)*A(2,2)
          IF(N .EQ. 2)  GOTO 320
          DO 314    J  =  3, N
             JM2  =  J - 2
             DO 312 K  =  1, JM2
                S31  =  A(K,J)
                DO 311  I  =  K, JM2
                   S31  =  DOTF(A(K,I+1),A(I+1,J),S31)
 311               CONTINUE
                A(K,J)  =  -S31
                A(J,K)  =  -S31*A(J,J)
 312            CONTINUE
             A(J-1,J)  =  -A(J-1,J)
             A(J,J-1)  =   A(J-1,J)*A(J,J)
 314         CONTINUE
 320      J  =  1
 323         S33  =  A(J,J)
             IF(J .EQ. N)  GOTO 325
             JP1  =  J + 1
             DO 324 I  =  JP1, N
                S33  =  DOTF(A(J,I),A(I,J),S33)
 324            CONTINUE
 325         A(J,J)  =  S33
          JM1  =  J
          J    =  JP1
             DO 328 K  =  1, JM1
                S32  =  ZERO
                DO 327  I  =  J, N
                   S32  =  DOTF(A(K,I),A(I,J),S32)
 327               CONTINUE
                A(K,J)  =  S32
                A(J,K)  =  S32
 328            CONTINUE
          IF(J .LT. N)  GOTO 323
 399      CONTINUE
cmsh # 26 "dsinv.F" 2
          RETURN
 900      CALL TMPRNT(HNAME,N,IDIM,0)
          RETURN
          END
*CMZ :          25/04/2017  14.16.56  by  Michael Scheer
*-- Author :    Michael Scheer   25/04/2017
*
cmsh +PATCH,//MSHCERN/FOR
cmsh +DECK,dsnleq.
cmsh
cmsh Michael Scheer: My changes are marked by cmsh
cmsh
*
* $Id: snleq64.F,v 1.1.1.1 1996/04/01 15:01:52 mclareni Exp $
*
* $Log: snleq64.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
cmsh #include "gen/pilot.h"
cmsh #if defined(CERNLIB_DOUBLE)
      SUBROUTINE DSNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) !cmsh
cmsh #include "gen/imp64.inc"
cmsh #endif
cmsh#if !defined(CERNLIB_DOUBLE)
cmsh      SUBROUTINE  RSNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
cmsh #endif

C     Based on   J.J. More  and  M.Y. Cosnard
C
C       ALGORITHM 554 BRENTM, A Fortran Subroutine for the
C       Numerical Solution of Systems of Nonlinear Equations [C5]
C
C     ACM Trans. Math. Software 6 (1980) 240-251.

      DIMENSION X(N),F(N),W(N,*),MPT(288)
      LOGICAL LCV

      PARAMETER (Z1 = 1, SCALE = 10, P05 = 5*Z1/100)
C**** EPS = SQRT(SMALLEST FP.NUMBER)
C     EPS = 1 / SQRT( 16D0**13 )
cmsh #if !defined(CERNLIB_DOUBLE)
cmsh      PARAMETER (EPS =  0.84293 69702 17878 97282 52636 392E-07)
cmsh #endif
cmsh #if defined(CERNLIB_IBM)
cmsh      PARAMETER (EPS =  0.14901 16119 38476 562D-07)
cmsh #endif
cmsh #if defined(CERNLIB_VAX)
cmsh      PARAMETER (EPS =  0.37252 90298 46191 40625D-08)
cmsh #endif
cmsh #if (defined(CERNLIB_UNIX))&&(defined(CERNLIB_DOUBLE))
cmsh      PARAMETER (EPS =  0.14901 16119 38476 600D-07)
      PARAMETER (EPS =  1.4901161193847656E-008) !cmsh
cmsh #endif
      DATA (MPT(I),I=1,288)
     1/1* 1,1* 2,3* 3,3* 4,4* 5,4* 6,4* 7,4* 8,5* 9,5*10,5*11,5*12,
     2 5*13,5*14,6*15,6*16,5*17,6*18,6*19,6*20,7*21,6*22,6*23,7*24,
     3 6*25,7*26,6*27,7*28,7*29,7*30,7*31,7*32,7*33,7*34,7*35,7*36,
     4 8*37,7*38,7*39,8*40,7*41,8*42,7*43,8*44,8*45,7*46,8*47,8*48/

      INFO=0
      IF(N .LE. 0 .OR. FTOL .LE. 0 .OR. XTOL .LE. 0) RETURN
C
C     Find optimal MOPT for iterative refinement
C
      IF(N .LE. 288) THEN
       MOPT=MPT(N)
      ELSE
       H=0
       DO 1 I = 49,N
       TEMP=LOG(I+Z1)/(N+2*I+1)
       IF(TEMP .LT. H) THEN
        MOPT=I-1
        GO TO 2
       ENDIF
    1  H=TEMP
      ENDIF

    2 IFLAG=0
      NUMF=0
      NFCALL=0

      NIER6=-1
      NIER7=-1
      NIER8=0
      FNORM=0
      DIFIT=0
      XNORM=0
      DO 10 I = 1,N
   10 XNORM=MAX(XNORM,ABS(X(I)))
      DELTA=SCALE*XNORM
      IF(XNORM .EQ. 0) DELTA=SCALE

   20 IF(IPRT .NE. 0) WRITE(6,'(1X,I5,D25.14)') (I,X(I),I=1,N)

      NSING=N
      FNORM1=FNORM
      DIFIT1=DIFIT
      FNORM=0
C
C     Compute step H for the divided difference which approximates
C     the K-th row of the Jacobian matrix
C
      H=EPS*XNORM
      IF(H .EQ. 0) H=EPS
      DO 40 J = 1,N
      DO 30 I = 1,N
   30 W(I,J+3)=0
      W(J,J+3)=H
   40 W(J,2)=X(J)
C
C     Enter a subiteration
C
      DO 150 K = 1,N
      IFLAG=K
      CALL SUB(N,W(1,2),F,IFLAG)
      FKY=F(K)
      NFCALL=NFCALL+1
      NUMF=NFCALL/N
      IF(IFLAG .LT. 0) GO TO 230
      FNORM=MAX(FNORM,ABS(FKY))
C
C     Compute the K-th row of the Jacobian matrix
C
      DO 60 J = K,N
      DO 50 I = 1,N
   50 W(I,3)=W(I,2)+W(I,J+3)
      CALL SUB(N,W(1,3),F,IFLAG)
      FKZ=F(K)
      NFCALL=NFCALL+1
      NUMF=NFCALL/N
      IF(IFLAG .LT. 0) GO TO 230
   60 W(J,1)=FKZ-FKY
      F(K)=FKY
C
C     Compute the Householder transformation to reduce the K-th row
C     of the Jacobian matrix to a multiple of the K-th unit vector
C
      ETA=0
      DO 70 I = K,N
   70 ETA=MAX(ETA,ABS(W(I,1)))
      IF(ETA .EQ. 0) GO TO 150
      NSING=NSING-1
      SKNORM=0
      DO 80 I = K,N
      W(I,1)=W(I,1)/ETA
   80 SKNORM=SKNORM+W(I,1)**2
      SKNORM=SQRT(SKNORM)
      IF(W(K,1) .LT. 0) SKNORM=-SKNORM
      W(K,1)=W(K,1)+SKNORM
C
C     Apply the transformation
C
      DO 90 I = 1,N
   90 W(I,3)=0
      DO 100 J = K,N
      DO 100 I = 1,N
  100 W(I,3)=W(I,3)+W(J,1)*W(I,J+3)
      DO 120 J = K,N
      TEMP=W(J,1)/(SKNORM*W(K,1))
      DO 120 I = 1,N
  120 W(I,J+3)=W(I,J+3)-TEMP*W(I,3)
C
C     Compute the subiterate
C
      W(K,1)=SKNORM*ETA
      TEMP=FKY/W(K,1)
      IF(H*ABS(TEMP) .GT. DELTA) TEMP=SIGN(DELTA/H,TEMP)
      DO 140 I = 1,N
  140 W(I,2)=W(I,2)+TEMP*W(I,K+3)
  150 CONTINUE
C
C     Compute the norms of the iterate and correction vector
C
      XNORM=0
      DIFIT=0
      DO 160 I = 1,N
      XNORM=MAX(XNORM,ABS(W(I,2)))
      DIFIT=MAX(DIFIT,ABS(X(I)-W(I,2)))
  160 X(I)=W(I,2)
C
C     Update the bound on the correction vector
C
      DELTA=MAX(DELTA,SCALE*XNORM)
C
C     Determine the progress of the iteration
C
      LCV=FNORM .LT. FNORM1 .AND. DIFIT .LT. DIFIT1 .AND. NSING .EQ. 0
      NIER6=NIER6+1
      NIER7=NIER7+1
      NIER8=NIER8+1
      IF(LCV) NIER6=0
      IF(FNORM .LT. FNORM1 .OR. DIFIT .LT. DIFIT1) NIER7=0
      IF(DIFIT .GT. EPS*XNORM) NIER8=0
C
C     Tests for convergence
C
      IF(FNORM .LE. FTOL) INFO=1
      IF(DIFIT .LE. XTOL*XNORM .AND. LCV) INFO=2
      IF(FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO=3
      IF(INFO .NE. 0) GO TO 230
C
C     Tests for termination
C
      IF(NUMF .GE. MAXF) INFO=4
      IF(NSING .EQ. N) INFO=5
      IF(NIER6 .EQ. 5) INFO=6
      IF(NIER7 .EQ. 3) INFO=7
      IF(NIER8 .EQ. 4) INFO=8
      IF(INFO .NE. 0) GO TO 230
      IF(.NOT.LCV .OR. DIFIT .GT. P05*XNORM) GO TO 20
C
C     Iterative refinement  (if the iteration is converging)
C
      DO 210 M = 2,MOPT
      FNORM1=FNORM
      FNORM=0
      DO 190 K = 1,N
      IFLAG=K
      CALL SUB(N,W(1,2),F,IFLAG)
      FKY=F(K)
      NFCALL=NFCALL+1
      NUMF=NFCALL/N
      IF(IFLAG .LT. 0) GO TO 230
      FNORM=MAX(FNORM,ABS(FKY))
C
C     Iterative refinement is terminated if it does not give a
C     reduction on residuals
C
      IF(FNORM .GE. FNORM1) THEN
       FNORM=FNORM1
       GO TO 20
      ENDIF
      TEMP=FKY/W(K,1)
      DO 180 I = 1,N
  180 W(I,2)=W(I,2)+TEMP*W(I,K+3)
  190 CONTINUE
C
C     Compute the norms of the iterate and correction vector
C
      XNORM=0
      DIFIT=0
      DO 200 I = 1,N
      XNORM=MAX(XNORM,ABS(W(I,2)))
      DIFIT=MAX(DIFIT,ABS(X(I)-W(I,2)))
  200 X(I)=W(I,2)
C
C     Stopping criteria for iterative refinement
C
      IF(FNORM .LE. FTOL) INFO=1
      IF(DIFIT .LE. XTOL*XNORM) INFO=2
      IF(FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO=3
      IF(NUMF .GE. MAXF .AND. INFO .EQ. 0) INFO=4
      IF(INFO .NE. 0) GO TO 230
  210 CONTINUE
      GO TO 20

  230 IF(IFLAG .LT. 0) INFO=IFLAG
      RETURN
      END
*CMZ :          02/05/2017  14.37.34  by  Michael Scheer
*-- Author :
# 1 "dsumsq.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dsumsq.F"
*
* $Id: dsumsq.F,v 1.1.1.1 1996/04/01 15:02:20 mclareni Exp $
*
* $Log: dsumsq.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*

# 1 "/usr/include/gen/pilot.h" 1 3 4
























# 40 "/usr/include/gen/pilot.h" 3 4

# 57 "/usr/include/gen/pilot.h" 3 4



























































# 10 "dsumsq.F" 2
      SUBROUTINE DSUMSQ(SUB,M,N,NC,A,AL,AU,MODE,EPS,MAXIT,IPRT,
     +                  MFR,IAFR,PHI,DPHI,COV,STD,W,NERROR)


# 1 "/usr/include/gen/imp64.inc" 1 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

# 14 "dsumsq.F" 2

      DIMENSION A(*),AL(*),AU(*),DPHI(*),IAFR(*)
      DIMENSION COV(NC,*),STD(*)
      DIMENSION W(*)
      double precision x(1),y(1),sy(1) !cmsh

      EXTERNAL SUB


***********************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   DSUMSQ IS THE STEERING ROUTINE FOR MINIMIZING A SUM OF SQUARES
*
*   SUBROUTINE CALLED:     D501L1
*
*
*   THE CONSTANTS, VARIABLES AND ARRAYS HAVE THE FOLLOWING MEANING.
*
*   SUB    NAME OF USER-SUPPLIED SUBROUTINE SUBPROGRAM, DECLARED
*          EXTERNAL IN THE CALLING PROGRAM. THIS SUBPROGRAM MUST PROVIDE
*          THE VALUES OF THE FUNCTION AND, IF MODE=1, THE VALUES OF THE
*          DERIVATIVES (SEE EXAMPLE) .
*   M      (INTEGER) NUMBER OF NONLINEAR FUNCTIONS.
*   N      (INTEGER) NUMBER OF UNKNOWN PARAMETERS A.
*   NC     (INTEGER) DECLARED FIRST DIMENSION OF ARRAY  COV  IN THE
*          CALLING PROGRAM, WITH  NC .GE. N .
*   A      (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY OF LENGTH  N .
*          ON ENTRY, A  MUST CONTAIN THE STARTING VALUES OF THE UNKNOWN
*          PARAMETERS FOR THE LEVENBERG-MARQUARDT ALGORITHM.
*          ON EXIT, A  CONTAINS AN APPROXIMATION OF THE MINIMUM POINT.
*   AL     (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY OF LENGTH  N .
*          ON ENTRY, AL  MUST CONTAIN THE LOWER BOUNDS OF  A .
*   AU     (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY OF LENGTH  N .
*          ON ENTRY, AU  MUST CONTAIN THE UPPER BOUNDS OF  A .
*   MODE   (INTEGER)
*          = 0: THE JACOBIAN IS COMPUTED NUMERICALLY
*          = 1: THE JACOBIAN HAS TO BE EVALUATED IN SUBPROGRAM  SUB .
*   EPS    (DOUBLE PRECISION) USER-SUPPLIED TOLERANCE USED TO CONTROL
*          THE TERMINATION CRITERION. EPS SHOULD BE CHOSEN ACCORDING
*          TO THE ACCURACY REQUIRED BY THE UNDERLYING PROBLEM AND TO
*          THE MACHINE ACCURACY ALSO (RECOMMENDED VALUE ON ENTRY:
*          1D-6 ... 1D-12 ).
*   MAXIT  (INTEGER) MAXIMUM PERMITTED NUMBER OF ITERATIONS.
*   IPRT   (INTEGER) PRINTING CONTROL.
*          = 0     : NO PRINTING OF INTERMEDIATE RESULTS
*          = +/- L : PRINTING OF INTERMEDIATE RESULTS AT EVERY ABS(L)-TH
*                    ITERATION; IF  IPRT < 0, PRINTING OF ALL INPUT
*                    PARAMETERS OF DSUMSQ IN ADDITION.
*   MFR    (INTEGER) ON EXIT, MFR CONTAINS THE NUMBER OF FREE VARIABLES
*          AT THE SOLUTION POINT.
*   IAFR   (INTEGER) ONE-DIMENSIONAL ARRAY OF LENGTH  2 * N , USED AS
*          WORKING SPACE. ON EXIT, THE FIRST  MFR  ELEMENTS OF  IAFR
*          CONTAIN THE INDICES OF THE FREE VARIABLES AT THE SOLUTION
*          POINT.
*   PHI    (DOUBLE PRECISION) ON EXIT, PHI  CONTAINS THE VALUE OF THE
*          OBJECTIVE FUNCTION AT THE MINIMUM POINT.
*   DPHI   (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY OF LENGTH  N .
*          ON EXIT, DPHI  CONTAINS THE DERIVATIVES OF THE OBJECTIVE
*          FUNCTION WITH RESPECT TO A (THE GRADIENT) AT THE LAST
*          ITERATION POINT.
*   COV    (DOUBLE PRECISION) TWO-DIMENSIONAL ARRAY OF DIMENSION (NC,N).
*          ON EXIT, COV CONTAINS AN APPROXIMATION TO THE COVARIANCE
*          MATRIX.
*   STD    (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY OF LENGTH  N .
*          ON EXIT, STD  CONTAINS APPROXIMATIONS TO THE STANDARD
*          DEVIATIONS OF THE MODEL PARAMETER ESTIMATORS.
*   W      (DOUBLE PRECISION) ONE-DIMENSIONAL ARRAY OF LENGTH
*          9*N+4*M+2*M*N+3*N*N , USED AS WORKING SPACE.
*   NERROR (INTEGER) ERROR INDICATOR. ON EXIT:
*           = 0: NO ERROR OR WARNING DETECTED.
*           = 1: AT LEAST ONE OF THE CONSTANTS M, N, NC, MAXIT IS
*                ILLEGAL OR AT LEAST FOR ONE J THE RELATION
*                AL(J) .LE. AU(J)  IS NOT TRUE.
*           = 2: THE MAXIMUM NUMBER  MAXIT  OF ITERATIONS HAS BEEN
*                REACHED.
*           = 3: THE OBJECTIVE FUNCTION  PHI  OR ITS DERIVATIVE IS NOT
*                DEFINED FOR THE CURRENT VALUES OF THE UNKNOWN
*                PARAMETER VECTOR  A.
*           = 4: THE ROUTINES  DGEQPF , DORMQR , DTRTRS  OF THE LINEAR
*                ALGEBRA PACKAGE  LAPACK (F001)  WERE UNABLE TO SOLVE
*                THE LINEAR LEAST SQUARES PROBLEMS
*                OR THE ROUTINE  DSINV (F012)  WAS UNABLE TO COMPUTE THE
*                COVARIANCE MATRIX .
*
*************************************************************************
*
*   THE FOLLOWING SUBROUTINE IS A SIMPLE EXAMPLE FOR SUB.
*
*     SUBROUTINE SUB(N,A,M,F,DF,MODE,NERROR)
*     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     PARAMETER (Z0 = 0)
*     DIMENSION A(*),F(*),DF(M,*)
*     NERROR=0
*     F(1)=A(1)-1D6
*     F(2)=A(2)-2D-6
*     F(3)=A(1)*A(2)-2
*     IF(MODE .EQ. 0) RETURN
*     CALL DMSET(M,N,Z0,DF(1,1),DF(1,2),DF(2,1))
*     DF(1,1)=1
*     DF(2,2)=1
*     DF(3,1)=A(2)
*     DF(3,2)=A(1)
*     RETURN
*     END
*
*************************************************************************

      M1=1
      M2=M1+N
      M3=M2+N
      M4=M3+N
      M5=M4+N
      M6=M5+2*M
      M7=M6+3*N
      M8=M7+N
      M9=M8+M+N
      MA=M9+(N+M)*N
      MB=MA+N*N
      MC=MB+N*N
      MD=MC+M

cmsh      CALL D501L1('DSUMSQ',SUB,1,M,X,1,Y,SY,MODE,EPS,MAXIT,
      CALL D501L1('DSUMSQ',SUB,1,M,X(1),1,Y(1),SY(1),MODE,EPS,MAXIT,
     1            IPRT,N,A,AL,AU,PHI,DPHI,IAFR,MFR,COV,NC,STD,
     2            W(M1),W(M2),W(M3),W(M4),W(M5),W(M6),W(M7),W(M8),
     3            W(M9),W(MA),W(MB),W(MC),W(MD),IAFR(N+1),NERROR)

      RETURN

      END
*CMZ :          02/05/2017  15.33.08  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
*CMZ :          02/05/2017  15.31.22  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*A'*B.
*
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END
*CMZ :          02/05/2017  15.34.51  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END
*CMZ :          02/05/2017  15.39.43  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
*CMZ :          02/05/2017  10.18.43  by  Michael Scheer
*-- Author :
      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTRS solves a triangular system of the form
*
*     A * X = B  or  A**T * X = B,
*
*  where A is a triangular matrix of order N, and B is an N-by-NRHS
*  matrix.  A check is made to verify that A is nonsingular.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
*          upper triangular part of the array A contains the upper
*          triangular matrix, and the strictly lower triangular part of
*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
*          triangular part of the array A contains the lower triangular
*          matrix, and the strictly upper triangular part of A is not
*          referenced.  If DIAG = 'U', the diagonal elements of A are
*          also not referenced and are assumed to be 1.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, if INFO = 0, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, the i-th diagonal element of A is zero,
*               indicating that the matrix is singular and the solutions
*               X have not been computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.
     $         LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      END IF
      INFO = 0
*
*     Solve A * x = b  or  A' * x = b.
*
      CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,
     $            LDB )
*
      RETURN
*
*     End of DTRTRS
*
      END
*CMZ :          02/05/2017  14.55.22  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 19 "dvadd.F" 2

cmsh # 1 "dvcpy.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "dvcpy.F"
*
* $Id: dvcpy.F,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: dvcpy.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "dvcpy.F" 2
          SUBROUTINE          DVCPY(N,X,X2,Z,Z2)
          DOUBLE PRECISION    X(*), X2(*), Z(*), Z2(*), FUNCT, A
          FUNCT(A)  =  A
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
cmsh # 15 "dvcpy.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
cmsh # 16 "dvcpy.F" 2

# 1 "vcpy.inc" 1
*
* $Id: vcpy.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: vcpy.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* vcpy.inc
*
          LXJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  FUNCT( X(LXJ) )
             LXJ     =  LXJ + JX
             LZJ     =  LZJ + JZ
  10         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  14.55.22  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 17 "dvcpy.F" 2

# 1 "dvdiv.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvdiv.F"
*
* $Id: dvdiv.F,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: dvdiv.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvdiv.F" 2
          SUBROUTINE          DVDIV(N,X,X2,Y,Y2,Z,Z2,IFAIL)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), Z(*), Z2(*), T
          REALF(T)  =  SNGL(T)
          IFAIL     =  0
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 16 "dvdiv.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 17 "dvdiv.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 18 "dvdiv.F" 2

# 1 "vdiv.inc" 1
*
* $Id: vdiv.inc,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: vdiv.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*
*
* vdiv.inc
*
          LXJ  =  1
          LYJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             IF(REALF(Y(LYJ)) .EQ. 0.)  GOTO 20
             Z(LZJ)  =  X(LXJ) / Y(LYJ)
             LXJ     =  LXJ + JX
             LYJ     =  LYJ + JY
             LZJ     =  LZJ + JZ
  10      CONTINUE
          J  =  0
  20      IFAIL  =  J
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 19 "dvdiv.F" 2

# 1 "dvmpa.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvmpa.F"
*
* $Id: dvmpa.F,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: dvmpa.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvmpa.F" 2
          DOUBLE PRECISION FUNCTION DVMPA(N,X,X2,Y,Y2,S)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), S, A, B
          DOUBLE PRECISION    SUM, MPA
          MPA(A,B,SUM)  =  A*B + SUM
          SUM  =  S
          IF(N .LE. 0)  GOTO 20

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 17 "dvmpa.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dvmpa.F" 2

# 1 "vmpa.inc" 1
*
* $Id: vmpa.inc,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: vmpa.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*
*
* vmpa.inc
*
          LXJ  =  1
          LYJ  =  1
          DO 10     J  =  1, N
             SUM  =  MPA( X(LXJ),Y(LYJ), SUM)
             LXJ  =  LXJ + JX
             LYJ  =  LYJ + JY
  10         CONTINUE
# 19 "dvmpa.F" 2
  20      DVMPA  =  SUM
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "dvmpy.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvmpy.F"
*
* $Id: dvmpy.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvmpy.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvmpy.F" 2
          DOUBLE PRECISION FUNCTION DVMPY(N,X,X2,Y,Y2)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), A, B
          DOUBLE PRECISION    SUM, MPA
          MPA(A,B,SUM)  =  A*B + SUM
          SUM  =  0.D0
          IF(N .LE. 0)  GOTO 20

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 17 "dvmpy.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 18 "dvmpy.F" 2

# 1 "vmpa.inc" 1
*
* $Id: vmpa.inc,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: vmpa.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*
*
* vmpa.inc
*
          LXJ  =  1
          LYJ  =  1
          DO 10     J  =  1, N
             SUM  =  MPA( X(LXJ),Y(LYJ), SUM)
             LXJ  =  LXJ + JX
             LYJ  =  LYJ + JY
  10         CONTINUE
# 19 "dvmpy.F" 2
  20      DVMPY  =  SUM
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 1 "dvmula.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvmula.F"
*
* $Id: dvmula.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvmula.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvmula.F" 2
          SUBROUTINE          DVMULA(N,X,X2,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), Z(*), Z2(*)
          DOUBLE PRECISION    MULA, A, B, C
          MULA(A,B,C)  =  A*B + C
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 16 "dvmula.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 17 "dvmula.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 18 "dvmula.F" 2

# 1 "vmula.inc" 1
*
* $Id: vmula.inc,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: vmula.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*
*
* vmula.inc
*
          LXJ  =  1
          LYJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  MULA( X(LXJ),Y(LYJ),Z(LZJ) )
             LXJ     =  LXJ + JX
             LYJ     =  LYJ + JY
             LZJ     =  LZJ + JZ
  10      CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 19 "dvmula.F" 2

# 1 "dvmul.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvmul.F"
*
* $Id: dvmul.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvmul.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvmul.F" 2
          SUBROUTINE          DVMUL(N,X,X2,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), Z(*), Z2(*)
          DOUBLE PRECISION    ADD, A, B
          ADD(A,B)  =  A*B
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 16 "dvmul.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 17 "dvmul.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 18 "dvmul.F" 2

# 1 "vadd.inc" 1
*
* $Id: vadd.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: vadd.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* vadd.inc
*
          LXJ  =  1
          LYJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  ADD( X(LXJ),Y(LYJ) )
             LXJ     =  LXJ + JX
             LYJ     =  LYJ + JY
             LZJ     =  LZJ + JZ
  10      CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 19 "dvmul.F" 2

# 1 "dvmuna.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvmuna.F"
*
* $Id: dvmuna.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvmuna.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvmuna.F" 2
          SUBROUTINE          DVMUNA(N,X,X2,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), Z(*), Z2(*)
          DOUBLE PRECISION    MULA, A, B, C
          MULA(A,B,C)  =  -A*B + C
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 16 "dvmuna.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
# 17 "dvmuna.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 18 "dvmuna.F" 2

# 1 "vmula.inc" 1
*
* $Id: vmula.inc,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: vmula.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*
*
* vmula.inc
*
          LXJ  =  1
          LYJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  MULA( X(LXJ),Y(LYJ),Z(LZJ) )
             LXJ     =  LXJ + JX
             LYJ     =  LYJ + JY
             LZJ     =  LZJ + JZ
  10      CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  13.26.29  by  Michael Scheer
*-- Author :
# 19 "dvmuna.F" 2

# 1 "dvran.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvran.F"
*
* $Id: dvran.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvran.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvran.F" 2
          SUBROUTINE          DVRAN(N,A,B,Z,Z2)
          DOUBLE PRECISION    A, B, C, Z(*), Z2(*)
          DOUBLE PRECISION    DRANF
          IF(N .LE. 0)  RETURN
          LZJ  =  1

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 16 "dvran.F" 2
          C    =  B - A
          DO 10     J  =  1, N
             Z(LZJ)  =  C*DRANF() + A
             LZJ     =  LZJ + JZ
  10         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  14.52.17  by  Michael Scheer
*-- Author :
cmsh # 1 "dvsca.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "dvsca.F"
*
* $Id: dvsca.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvsca.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "dvsca.F" 2
          SUBROUTINE          DVSCA(N,S,X,X2,Y,Y2,Z,Z2)
          DOUBLE PRECISION    S, X(*), X2(*), Y(*), Y2(*), Z(*), Z2(*)
          DOUBLE PRECISION    ADD, A, B
          ADD(A,B)  =  S*A + B
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
cmsh # 16 "dvsca.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
cmsh # 17 "dvsca.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
cmsh # 18 "dvsca.F" 2

# 1 "vadd.inc" 1
*
* $Id: vadd.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: vadd.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* vadd.inc
*
          LXJ  =  1
          LYJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  ADD( X(LXJ),Y(LYJ) )
             LXJ     =  LXJ + JX
             LYJ     =  LYJ + JY
             LZJ     =  LZJ + JZ
  10      CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  14.52.17  by  Michael Scheer
*-- Author :
cmsh # 19 "dvsca.F" 2

# 1 "dvscl.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvscl.F"
*
* $Id: dvscl.F,v 1.1.1.1 1996/02/15 17:48:52 mclareni Exp $
*
* $Log: dvscl.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:52  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvscl.F" 2
          SUBROUTINE          DVSCL(N,S,X,X2,Z,Z2)
          DOUBLE PRECISION    S, X(*), X2(*), Z(*), Z2(*), FUNCT, A
          FUNCT(A)  =  S*A
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 15 "dvscl.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
# 16 "dvscl.F" 2

# 1 "vcpy.inc" 1
*
* $Id: vcpy.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: vcpy.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* vcpy.inc
*
          LXJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  FUNCT( X(LXJ) )
             LXJ     =  LXJ + JX
             LZJ     =  LZJ + JZ
  10         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  14.51.45  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh# 19 "dvscs.F" 2

cmsh # 1 "dvset.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "dvset.F"
*
* $Id: dvset.F,v 1.1.1.1 1996/02/15 17:48:53 mclareni Exp $
*
* $Log: dvset.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:53  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "dvset.F" 2
          SUBROUTINE          DVSET(N,S,Z,Z2)
          DOUBLE PRECISION    S, Z(*), Z2(*)
          IF(N .LE. 0)  RETURN

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
cmsh # 14 "dvset.F" 2

# 1 "vset.inc" 1
*
* $Id: vset.inc,v 1.1.1.1 1996/02/15 17:48:51 mclareni Exp $
*
* $Log: vset.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:51  mclareni
* Kernlib
*
*
*
* vset.inc
*
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  S
             LZJ     =  LZJ + JZ
  10         CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  14.54.26  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 15 "dvset.F" 2

cmsh # 1 "dvsub.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "dvsub.F"
*
* $Id: dvsub.F,v 1.1.1.1 1996/02/15 17:48:53 mclareni Exp $
*
* $Log: dvsub.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:53  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "dvsub.F" 2
          SUBROUTINE          DVSUB(N,X,X2,Y,Y2,Z,Z2)
          DOUBLE PRECISION    X(*), X2(*), Y(*), Y2(*), Z(*), Z2(*)
          DOUBLE PRECISION    ADD, A, B
          ADD(A,B)  =  A-B
          IF(N .LE. 0)  RETURN

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
cmsh # 16 "dvsub.F" 2

# 1 "dyj.inc" 1
*
* $Id: dyj.inc,v 1.1.1.1 1996/02/15 17:48:55 mclareni Exp $
*
* $Log: dyj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:55  mclareni
* Kernlib
*
*
*
* dyj.inc
*

          JY  =  (LOCF(Y2) - LOCF(Y)) / 2
cmsh # 17 "dvsub.F" 2

# 1 "dzj.inc" 1
*
* $Id: dzj.inc,v 1.1.1.1 1996/02/15 17:48:56 mclareni Exp $
*
* $Log: dzj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:56  mclareni
* Kernlib
*
*
*
* dzj.inc
*

          JZ  =  (LOCF(Z2) - LOCF(Z)) / 2
cmsh # 18 "dvsub.F" 2

# 1 "vadd.inc" 1
*
* $Id: vadd.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: vadd.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* vadd.inc
*
          LXJ  =  1
          LYJ  =  1
          LZJ  =  1
          DO 10     J  =  1, N
             Z(LZJ)  =  ADD( X(LXJ),Y(LYJ) )
             LXJ     =  LXJ + JX
             LYJ     =  LYJ + JY
             LZJ     =  LZJ + JZ
  10      CONTINUE
          RETURN
          END
*CMZ :          02/05/2017  14.54.26  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 19 "dvsub.F" 2

# 1 "dvsum.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dvsum.F"
*
* $Id: dvsum.F,v 1.1.1.1 1996/02/15 17:48:53 mclareni Exp $
*
* $Log: dvsum.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:53  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

# 10 "dvsum.F" 2
          DOUBLE PRECISION FUNCTION DVSUM(N,X,X2)
          DOUBLE PRECISION    X(*), X2(*), SUM
          SUM  =  0.D0
          IF(N .LE. 0)  GOTO 20

# 1 "dxj.inc" 1
*
* $Id: dxj.inc,v 1.1.1.1 1996/02/15 17:48:50 mclareni Exp $
*
* $Log: dxj.inc,v $
* Revision 1.1.1.1  1996/02/15 17:48:50  mclareni
* Kernlib
*
*
*
* dxj.inc
*

          JX  =  (LOCF(X2) - LOCF(X)) / 2
# 15 "dvsum.F" 2
          LXJ  =  1
          DO 10     J  =  1, N
             SUM  =  SUM + X(LXJ)
             LXJ  =  LXJ + JX
  10         CONTINUE
  20      DVSUM  =  SUM
          RETURN
          END
*CMZ :  1.16/04 17/04/2014  11.22.45  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*# 1 "f010pr.F"
*# 1 "<command-line>"
*# 1 "f010pr.F"
*
* $Id: f010pr.F,v 1.1.1.1 1996/02/15 17:48:49 mclareni Exp $
*
* $Log: f010pr.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:49  mclareni
* Kernlib
*
*

*# 1 "/home/scheer/cern/cern/2005/src/packlib/kernlib/kernnum/f011fort/kernnum/pilot.h" 1
*# 21 "/home/scheer/cern/cern/2005/src/packlib/kernlib/kernnum/f011fort/kernnum/pilot.h"

*# 33 "/home/scheer/cern/cern/2005/src/packlib/kernlib/kernnum/f011fort/kernnum/pilot.h"

*# 10 "f010pr.F" 2
      SUBROUTINE F010PR(NAME,N,IDIM,K,KPRNT)
      CHARACTER*6 NAME
      LOGICAL MFLAG,RFLAG
C
C     ******************************************************************
C
C     PRINT ROUTINE FOR PARAMETER ERRORS IN MATRIX SUBROUTINES $EQINV,
C     $EQN, $INV (WHERE $ IS A LETTER SPECIFYING THE ARITHMETIC TYPE).
C
C     NAME         (CHARACTER*6) NAME OF THE CALLING ROUTINE.
C
C     N,IDIM,K     PARAMETERS OF THE CALLING ROUTINE (WITH K=0 IF K IS
C                  NOT TO BE PRINTED).
C
C     KPRNT        PRINT FLAG FOR K (K IS NOT PRINTED IF KPRNT=0).
C
C     ******************************************************************
C
C  START.
      CALL KERMTR('F010.1',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
         IF(LGFILE.EQ.0)  THEN
            IF(KPRNT.EQ.0) WRITE(*,2000) NAME,N,IDIM
            IF(KPRNT.NE.0) WRITE(*,2001) NAME,N,IDIM,K
         ELSE
            IF(KPRNT.EQ.0) WRITE(LGFILE,2000) NAME,N,IDIM
            IF(KPRNT.NE.0) WRITE(LGFILE,2001) NAME,N,IDIM,K
         ENDIF
      ENDIF
      IF(.NOT. RFLAG) CALL ABEND
      RETURN
C
 2000 FORMAT( 7X, 11HSUBROUTINE , A6, 14H ... PARAMETER,
     *        29H ERROR (N.LT.1 OR N.GT.IDIM).,
     *        6X, 3HN =, I4, 6X, 6HIDIM =, I4, 1H. )
 2001 FORMAT( 7X, 11HSUBROUTINE , A6, 14H ... PARAMETER,
     *        39H ERROR (N.LT.1 OR N.GT.IDIM OR K.LT.1).,
     *        6X, 3HN =, I4, 6X, 6HIDIM =, I4, 6X, 3HK =, I4, 1H. )
      END
*CMZ :          07/10/2014  14.31.07  by  Michael Scheer
*-- Author :    Michael Scheer   07/10/2014
*
* $Id: gamma64.F,v 1.1.1.1 1996/04/01 15:01:54 mclareni Exp $
*
* $Log: gamma64.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
* Mathlib gen
*
*

      FUNCTION DGAMMA(X)
C
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C
      CHARACTER*(*) NAME
      PARAMETER(NAME='GAMMA/DGAMMA')

C
      CHARACTER*80 ERRTXT

      DIMENSION C(0:15)

      DATA C( 0) /3.65738 77250 83382 44D0/
      DATA C( 1) /1.95754 34566 61268 27D0/
      DATA C( 2) /0.33829 71138 26160 39D0/
      DATA C( 3) /0.04208 95127 65575 49D0/
      DATA C( 4) /0.00428 76504 82129 09D0/
      DATA C( 5) /0.00036 52121 69294 62D0/
      DATA C( 6) /0.00002 74006 42226 42D0/
      DATA C( 7) /0.00000 18124 02333 65D0/
      DATA C( 8) /0.00000 01096 57758 66D0/
      DATA C( 9) /0.00000 00059 87184 05D0/
      DATA C(10) /0.00000 00003 07690 81D0/
      DATA C(11) /0.00000 00000 14317 93D0/
      DATA C(12) /0.00000 00000 00651 09D0/
      DATA C(13) /0.00000 00000 00025 96D0/
      DATA C(14) /0.00000 00000 00001 11D0/
      DATA C(15) /0.00000 00000 00000 04D0/

      U=X
      IF(U .LE. 0) THEN
       WRITE(ERRTXT,101) U
       CALL MTLPRT(NAME,'C302.1',ERRTXT)
       H=0
       GO TO 9
      ENDIF
    8 F=1
      IF(U .LT. 3) THEN
       DO 1 I = 1,INT(4-U)
       F=F/U
    1  U=U+1
      ELSE
       DO 2 I = 1,INT(U-3)
       U=U-1
    2  F=F*U
      END IF
      H=U+U-7
      ALFA=H+H
      B1=0
      B2=0
      DO 3 I = 15,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    3 B1=B0

    9 DGAMMA=F*(B0-H*B2)

      RETURN
  101 FORMAT('ARGUMENT IS NEGATIVE = ',1P,E15.1)
      END
*CMZ :          21/11/2017  14.38.13  by  Michael Scheer
*-- Author :
# 1 "/opt/cern/pro/src/mathlib/gen/c/gamma64.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "/opt/cern/pro/src/mathlib/gen/c/gamma64.F"
*
* $Id: gamma64.F,v 1.1.1.1 1996/04/01 15:01:54 mclareni Exp $
*
* $Log: gamma64.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:54 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "/opt/cern/pro/src/mathlib/gen/c/gamma64.F" 2
# 20 "/opt/cern/pro/src/mathlib/gen/c/gamma64.F"
      FUNCTION GAMMA(X)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='GAMMA')

C
      CHARACTER*80 ERRTXT

      DIMENSION C(0:15)

      DATA C( 0) /3.65738 77250 83382 44D0/
      DATA C( 1) /1.95754 34566 61268 27D0/
      DATA C( 2) /0.33829 71138 26160 39D0/
      DATA C( 3) /0.04208 95127 65575 49D0/
      DATA C( 4) /0.00428 76504 82129 09D0/
      DATA C( 5) /0.00036 52121 69294 62D0/
      DATA C( 6) /0.00002 74006 42226 42D0/
      DATA C( 7) /0.00000 18124 02333 65D0/
      DATA C( 8) /0.00000 01096 57758 66D0/
      DATA C( 9) /0.00000 00059 87184 05D0/
      DATA C(10) /0.00000 00003 07690 81D0/
      DATA C(11) /0.00000 00000 14317 93D0/
      DATA C(12) /0.00000 00000 00651 09D0/
      DATA C(13) /0.00000 00000 00025 96D0/
      DATA C(14) /0.00000 00000 00001 11D0/
      DATA C(15) /0.00000 00000 00000 04D0/

      U=X
      IF(U .LE. 0) THEN
       WRITE(ERRTXT,101) U
       CALL MTLPRT(NAME,'C302.1',ERRTXT)
       H=0
       GO TO 9
      ENDIF
    8 F=1
      IF(U .LT. 3) THEN
       DO 1 I = 1,INT(4-U)
       F=F/U
    1 U=U+1
      ELSE
       DO 2 I = 1,INT(U-3)
       U=U-1
    2 F=F*U
      END IF
      H=U+U-7
      ALFA=H+H
      B1=0
      B2=0
      DO 3 I = 15,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    3 B1=B0




    9 GAMMA =F*(B0-H*B2)

      RETURN
  101 FORMAT('ARGUMENT IS NEGATIVE = ',1P,E15.1)
      END
*CMZ :          28/08/2014  14.10.36  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX icfnbl.F

# 1 "icfnbl.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "icfnbl.F"
*
* $Id: icfnbl.F,v 1.1.1.1 1996/02/15 17:49:45 mclareni Exp $
*
* $Log: icfnbl.F,v $
* Revision 1.1.1.1 1996/02/15 17:49:45 mclareni
* Kernlib
*
*
# 1 "/usr/include/kerngen/pilot.h" 1 3 4
# 10 "icfnbl.F" 2
      FUNCTION ICFNBL (CHV,JLP,JRP)
C
C CERN PROGLIB# M432 ICFNBL .VERSION KERNFOR 4.21 890323
C ORIG. 04/10/88, JZ
C
C- Find first non-blank character in CHV(JL:JR)

      DIMENSION JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER CHV*(*)

      JJ = JLP(1)
      JR = JRP(1)

   12 IF (JJ.GT.JR) GO TO 19
      IF (CHV(JJ:JJ).EQ.' ') THEN
          JJ = JJ + 1
          GO TO 12
        ENDIF
      NGSLAT = JJ
      ICFNBL = JJ
      RETURN

   19 NGSLAT = 0
      ICFNBL = JJ
      RETURN
      END
*CMZ :          02/05/2017  15.40.53  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
*CMZ :          02/05/2017  15.22.27  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
*CMZ :          02/05/2017  15.21.17  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
*     End of ILAENV
*
      END
*CMZ :  1.16/04 17/04/2014  11.35.08  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*# 1 "kerset.F"
*# 1 "<command-line>"
*# 1 "kerset.F"
*
* $Id: kerset.F,v 1.1.1.1 1996/02/15 17:48:35 mclareni Exp $
*
* $Log: kerset.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:35  mclareni
* Kernlib
*
*

*# 1 "kernnum/pilot.h" 1
*# 21 "kernnum/pilot.h"

*# 33 "kernnum/pilot.h"

*# 10 "kerset.F" 2
          SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
                    PARAMETER(KOUNTE  =  27)
          CHARACTER*6         ERCODE,   CODE(KOUNTE)
          LOGICAL             MFLAG,    RFLAG
          INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)
          DATA      LOGF      /  0  /
          DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 255, 255 /
          DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 255, 255 /
          DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 255, 255 /
          DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 255, 255 /
          DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 255, 255 /
          DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 255, 255 /
          DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 255, 255 /
          DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 255, 255 /
          DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 255, 255 /
          DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 255, 255 /
          DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 255, 255 /
          DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 255, 255 /
          DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 255, 255 /
          DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 255, 255 /
          DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 255, 255 /
          DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 255, 255 /
          DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 255, 255 /
          DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 255, 255 /
          DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 255, 255 /
          DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 255, 255 /
          DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 255, 255 /
          DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 255,   0 /
          DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 255,   0 /
          DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 255,   0 /
          DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 255,   0 /
          DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 255, 255 /
          DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 255, 255 /
          LOGF  =  LGFILE
             L  =  0
          IF(ERCODE .NE. ' ')  THEN
             DO 10  L = 1, 6
                IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
  10            CONTINUE
  12         CONTINUE
          ENDIF
          DO 14     I  =  1, KOUNTE
             IF(L .EQ. 0)  GOTO 13
             IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
  13         IF(LIMITM.GE.0) KNTM(I)  =  LIMITM
             IF(LIMITR.GE.0) KNTR(I)  =  LIMITR
  14         CONTINUE
          RETURN
          ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
          LOG  =  LOGF
          DO 20     I  =  1, KOUNTE
             IF(ERCODE .EQ. CODE(I))  GOTO 21
  20         CONTINUE
          WRITE(*,1000)  ERCODE
          CALL ABEND
          RETURN
  21      RFLAG  =  KNTR(I) .GE. 1
          IF(RFLAG  .AND.  (KNTR(I) .LT. 255))  KNTR(I)  =  KNTR(I) - 1
          MFLAG  =  KNTM(I) .GE. 1
          IF(MFLAG  .AND.  (KNTM(I) .LT. 255))  KNTM(I)  =  KNTM(I) - 1
          IF(.NOT. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1001)  CODE(I)
             ELSE
                WRITE(LOGF,1001)  CODE(I)
             ENDIF
          ENDIF
          IF(MFLAG .AND. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1002)  CODE(I)
             ELSE
                WRITE(LOGF,1002)  CODE(I)
             ENDIF
          ENDIF
          RETURN
1000      FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +           ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +           ' ERROR MONITOR. RUN ABORTED.')
1001      FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +           'CONDITION ',A6)
1002      FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
          END
*CMZ :          02/05/2017  14.47.10  by  Michael Scheer
*-- Author :

*KEEP,cmsh,T=F77.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.

cmsh # 16 "dvxch.F" 2

* $Id: locf.F,v 1.1.1.1 1996/02/15 17:50:37 mclareni Exp $
*
* $Log: locf.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
*KEEP,CMSH.
!
!       Routine were taken from the CERNLIB
!       Changes by Michael Scheer are marked by "cmsh"
!
*KEND.
cmsh #include "kerngen/pilot.h"
cmsh #if defined(CERNLIB_QMMPW)
cmsh #include "mpwgs/locf.F"
cmsh #elif defined(CERNLIB_QMSUN)
cmsh #include "sungs/locf.F"
cmsh #elif defined(CERNLIB_QMVAX)
cmsh #include "vaxgs/locf.F"
cmsh #else
      FUNCTION LOCF (IVAR)
C
C CERN PROGLIBcmsh # N100    LOCF            .VERSION KERNFOR  4.34  930114
C
C-    This is a default which works on several machines
C
cmsh      DIMENSION    IVAR(9)
      double precision ivar(9) !cmsh
cmsh #if defined(CERNLIB_QMLXIA64)
      INTEGER*8    J
cmsh #endif
*    Number of ADdress Units Per Word for Fortran
C                         and its logarithm base 2
      PARAMETER    (NADUPW=4, LADUPW=2)

      J = LOC(IVAR)
cmsh #if defined(CERNLIB_QMLXIA64)
      IF (IAND (J, Z'FFFFFFFF00000000') .NE. 0 ) THEN
        WRITE(*,'(57(1H!))')
        WRITE(*,'(A,Z16,A)') 'LOCF: address #', J,
     &                       '# exceeds the 32 bit space'
        WRITE(*,'(57(1H!))')
       END IF
cmsh #endif
      LOCF = ISHFT (J, -LADUPW)
      END
*CMZ :          02/05/2017  15.07.58  by  Michael Scheer
*-- Author :
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
*CMZ :          29/08/2014  08.41.41  by  Michael Scheer
*-- Author :    Michael Scheer   29/08/2014
      real function msh_hrndm1(id)

      implicit none

      integer id

      print*,'*** Error in msh_hrndm1: Function not yet implemented!'

      msh_hrndm1=0.0*float(id)

      return
      end
*CMZ :          28/08/2014  11.46.46  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
*
* $Id: mtlprt.F,v 1.1.1.1 1996/04/01 15:02:52 mclareni Exp $
*
* $Log: mtlprt.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:52  mclareni
* Mathlib gen
*
*
cmsh #include "gen/pilot.h"
      SUBROUTINE MTLPRT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT
      LOGICAL LMF,LRF

      IF(ERC(5:6).NE.'.0') THEN
        CALL MTLMTR(ERC,MLG,LMF,LRF)
      ELSE
        LMF=.TRUE.
        LRF=.FALSE.
      ENDIF
      IF(LMF) THEN
cmsh        LT=LENOCC(TEXT)
        LT=len_trim(TEXT)
        IF(MLG .LT. 1) WRITE(  *,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
        IF(MLG .GE. 1) WRITE(MLG,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
      ENDIF
      IF(.NOT.LRF) CALL ABEND
      RETURN
100   FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END
*CMZ :          07/10/2014  14.02.53  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
*
* $Id: mtlset.F,v 1.1.1.1 1996/04/01 15:02:53 mclareni Exp $
*
* $Log: mtlset.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
      SUBROUTINE MTLSET(ERC,NLG,MXM,MXR)

      PARAMETER (KTE = 132)
      CHARACTER*6 ERC,CODE(KTE)
      LOGICAL LMF,LRF
      DIMENSION KNTM(KTE),KNTR(KTE)

      DATA ILG /0/

C     renumber the data statements after putting new codes in Unix with:
C     awk -F'[()]' '{ printf"%s(%s)%s(%s)%s(%s)%s\n",$1,NR,$3,NR,$5,NR,$7 }'
C     and modify KTE to the number of lines below

      DATA CODE(1),KNTM(1),KNTR(1) / 'B100.1', 255, 255 /
      DATA CODE(2),KNTM(2),KNTR(2) / 'B300.1', 255, 255 /
      DATA CODE(3),KNTM(3),KNTR(3) / 'B300.2', 255, 255 /
      DATA CODE(4),KNTM(4),KNTR(4) / 'C200.0', 255, 255 /
      DATA CODE(5),KNTM(5),KNTR(5) / 'C200.1', 255, 255 /
      DATA CODE(6),KNTM(6),KNTR(6) / 'C200.2', 255, 255 /
      DATA CODE(7),KNTM(7),KNTR(7) / 'C200.3', 255, 255 /
      DATA CODE(8),KNTM(8),KNTR(8) / 'C201.0', 255, 255 /
      DATA CODE(9),KNTM(9),KNTR(9) / 'C202.0', 255, 255 /
      DATA CODE(10),KNTM(10),KNTR(10) / 'C202.1', 255, 255 /
      DATA CODE(11),KNTM(11),KNTR(11) / 'C202.2', 255, 255 /
      DATA CODE(12),KNTM(12),KNTR(12) / 'C205.1', 255, 255 /
      DATA CODE(13),KNTM(13),KNTR(13) / 'C205.2', 255, 255 /
      DATA CODE(14),KNTM(14),KNTR(14) / 'C207.0', 255, 255 /
      DATA CODE(15),KNTM(15),KNTR(15) / 'C208.0', 255, 255 /
      DATA CODE(16),KNTM(16),KNTR(16) / 'C209.0', 255, 255 /
      DATA CODE(17),KNTM(17),KNTR(17) / 'C209.1', 255, 255 /
      DATA CODE(18),KNTM(18),KNTR(18) / 'C209.2', 255, 255 /
      DATA CODE(19),KNTM(19),KNTR(19) / 'C209.3', 255, 255 /
      DATA CODE(20),KNTM(20),KNTR(20) / 'C210.1', 255, 255 /
      DATA CODE(21),KNTM(21),KNTR(21) / 'C302.1', 255, 255 /
      DATA CODE(22),KNTM(22),KNTR(22) / 'C303.1', 255, 255 /
      DATA CODE(23),KNTM(23),KNTR(23) / 'C304.1', 255, 255 /
      DATA CODE(24),KNTM(24),KNTR(24) / 'C305.1', 255, 255 /
      DATA CODE(25),KNTM(25),KNTR(25) / 'C306.1', 255, 255 /
      DATA CODE(26),KNTM(26),KNTR(26) / 'C307.1', 255, 255 /
      DATA CODE(27),KNTM(27),KNTR(27) / 'C312.1', 255, 255 /
      DATA CODE(28),KNTM(28),KNTR(28) / 'C313.1', 255, 255 /
      DATA CODE(29),KNTM(29),KNTR(29) / 'C315.1', 255, 255 /
      DATA CODE(30),KNTM(30),KNTR(30) / 'C316.1', 255, 255 /
      DATA CODE(31),KNTM(31),KNTR(31) / 'C316.2', 255, 255 /
      DATA CODE(32),KNTM(32),KNTR(32) / 'C320.1', 255, 255 /
      DATA CODE(33),KNTM(33),KNTR(33) / 'C321.1', 255, 255 /
      DATA CODE(34),KNTM(34),KNTR(34) / 'C323.1', 255, 255 /
      DATA CODE(35),KNTM(35),KNTR(35) / 'C327.1', 255, 255 /
      DATA CODE(36),KNTM(36),KNTR(36) / 'C328.1', 255, 255 /
      DATA CODE(37),KNTM(37),KNTR(37) / 'C328.2', 255, 255 /
      DATA CODE(38),KNTM(38),KNTR(38) / 'C328.3', 255, 255 /
      DATA CODE(39),KNTM(39),KNTR(39) / 'C330.1', 255, 255 /
      DATA CODE(40),KNTM(40),KNTR(40) / 'C330.2', 255, 255 /
      DATA CODE(41),KNTM(41),KNTR(41) / 'C330.3', 255, 255 /
      DATA CODE(42),KNTM(42),KNTR(42) / 'C331.1', 255, 255 /
      DATA CODE(43),KNTM(43),KNTR(43) / 'C331.2', 255, 255 /
      DATA CODE(44),KNTM(44),KNTR(44) / 'C334.1', 255, 255 /
      DATA CODE(45),KNTM(45),KNTR(45) / 'C334.2', 255, 255 /
      DATA CODE(46),KNTM(46),KNTR(46) / 'C334.3', 255, 255 /
      DATA CODE(47),KNTM(47),KNTR(47) / 'C334.4', 255, 255 /
      DATA CODE(48),KNTM(48),KNTR(48) / 'C334.5', 255, 255 /
      DATA CODE(49),KNTM(49),KNTR(49) / 'C334.6', 255, 255 /
      DATA CODE(50),KNTM(50),KNTR(50) / 'C336.1', 255, 255 /
      DATA CODE(51),KNTM(51),KNTR(51) / 'C337.1', 255, 255 /
      DATA CODE(52),KNTM(52),KNTR(52) / 'C338.1', 255, 255 /
      DATA CODE(53),KNTM(53),KNTR(53) / 'C340.1', 255, 255 /
      DATA CODE(54),KNTM(54),KNTR(54) / 'C343.1', 255, 255 /
      DATA CODE(55),KNTM(55),KNTR(55) / 'C343.2', 255, 255 /
      DATA CODE(56),KNTM(56),KNTR(56) / 'C343.3', 255, 255 /
      DATA CODE(57),KNTM(57),KNTR(57) / 'C343.4', 255, 255 /
      DATA CODE(58),KNTM(58),KNTR(58) / 'C344.1', 255, 255 /
      DATA CODE(59),KNTM(59),KNTR(59) / 'C344.2', 255, 255 /
      DATA CODE(60),KNTM(60),KNTR(60) / 'C344.3', 255, 255 /
      DATA CODE(61),KNTM(61),KNTR(61) / 'C344.4', 255, 255 /
      DATA CODE(62),KNTM(62),KNTR(62) / 'C345.1', 255, 255 /
      DATA CODE(63),KNTM(63),KNTR(63) / 'C346.1', 255, 255 /
      DATA CODE(64),KNTM(64),KNTR(64) / 'C346.2', 255, 255 /
      DATA CODE(65),KNTM(65),KNTR(65) / 'C346.3', 255, 255 /
      DATA CODE(66),KNTM(66),KNTR(66) / 'C347.1', 255, 255 /
      DATA CODE(67),KNTM(67),KNTR(67) / 'C347.2', 255, 255 /
      DATA CODE(68),KNTM(68),KNTR(68) / 'C347.3', 255, 255 /
      DATA CODE(69),KNTM(69),KNTR(69) / 'C347.4', 255, 255 /
      DATA CODE(70),KNTM(70),KNTR(70) / 'C347.5', 255, 255 /
      DATA CODE(71),KNTM(71),KNTR(71) / 'C347.6', 255, 255 /
      DATA CODE(72),KNTM(72),KNTR(72) / 'C348.1', 255, 255 /
      DATA CODE(73),KNTM(73),KNTR(73) / 'C349.1', 255, 255 /
      DATA CODE(74),KNTM(74),KNTR(74) / 'C349.2', 255, 255 /
      DATA CODE(75),KNTM(75),KNTR(75) / 'C349.3', 255, 255 /
      DATA CODE(76),KNTM(76),KNTR(76) / 'D101.1', 255, 255 /
      DATA CODE(77),KNTM(77),KNTR(77) / 'D103.1', 255, 255 /
      DATA CODE(78),KNTM(78),KNTR(78) / 'D104.1', 255, 255 /
      DATA CODE(79),KNTM(79),KNTR(79) / 'D104.2', 255, 255 /
      DATA CODE(80),KNTM(80),KNTR(80) / 'D105.1', 255, 255 /
      DATA CODE(81),KNTM(81),KNTR(81) / 'D105.2', 255, 255 /
      DATA CODE(82),KNTM(82),KNTR(82) / 'D107.1', 255, 255 /
      DATA CODE(83),KNTM(83),KNTR(83) / 'D110.0', 255, 255 /
      DATA CODE(84),KNTM(84),KNTR(84) / 'D110.1', 255, 255 /
      DATA CODE(85),KNTM(85),KNTR(85) / 'D110.2', 255, 255 /
      DATA CODE(86),KNTM(86),KNTR(86) / 'D110.3', 255, 255 /
      DATA CODE(87),KNTM(87),KNTR(87) / 'D110.4', 255, 255 /
      DATA CODE(88),KNTM(88),KNTR(88) / 'D110.5', 255, 255 /
      DATA CODE(89),KNTM(89),KNTR(89) / 'D110.6', 255, 255 /
      DATA CODE(90),KNTM(90),KNTR(90) / 'D113.1', 255, 255 /
      DATA CODE(91),KNTM(91),KNTR(91) / 'D201.1', 255, 255 /
      DATA CODE(92),KNTM(92),KNTR(92) / 'D202.1', 255, 255 /
      DATA CODE(93),KNTM(93),KNTR(93) / 'D401.1', 255, 255 /
      DATA CODE(94),KNTM(94),KNTR(94) / 'D601.1', 255, 255 /
      DATA CODE(95),KNTM(95),KNTR(95) / 'E210.1', 255, 255 /
      DATA CODE(96),KNTM(96),KNTR(96) / 'E210.2', 255, 255 /
      DATA CODE(97),KNTM(97),KNTR(97) / 'E210.3', 255, 255 /
      DATA CODE(98),KNTM(98),KNTR(98) / 'E210.4', 255, 255 /
      DATA CODE(99),KNTM(99),KNTR(99) / 'E210.5', 255, 255 /
      DATA CODE(100),KNTM(100),KNTR(100) / 'E210.6', 255, 255 /
      DATA CODE(101),KNTM(101),KNTR(101) / 'E210.7', 255, 255 /
      DATA CODE(102),KNTM(102),KNTR(102) / 'E211.0', 255, 255 /
      DATA CODE(103),KNTM(103),KNTR(103) / 'E211.1', 255, 255 /
      DATA CODE(104),KNTM(104),KNTR(104) / 'E211.2', 255, 255 /
      DATA CODE(105),KNTM(105),KNTR(105) / 'E211.3', 255, 255 /
      DATA CODE(106),KNTM(106),KNTR(106) / 'E211.4', 255, 255 /
      DATA CODE(107),KNTM(107),KNTR(107) / 'E406.0', 255, 255 /
      DATA CODE(108),KNTM(108),KNTR(108) / 'E406.1', 255, 255 /
      DATA CODE(109),KNTM(109),KNTR(109) / 'E407.0', 255, 255 /
      DATA CODE(110),KNTM(110),KNTR(110) / 'E408.0', 255, 255 /
      DATA CODE(111),KNTM(111),KNTR(111) / 'E408.1', 255, 255 /
      DATA CODE(112),KNTM(112),KNTR(112) / 'F500.0', 255, 255 /
      DATA CODE(113),KNTM(113),KNTR(113) / 'F500.1', 255, 255 /
      DATA CODE(114),KNTM(114),KNTR(114) / 'F500.2', 255, 255 /
      DATA CODE(115),KNTM(115),KNTR(115) / 'F500.3', 255, 255 /
      DATA CODE(116),KNTM(116),KNTR(116) / 'G100.1', 255, 255 /
      DATA CODE(117),KNTM(117),KNTR(117) / 'G100.2', 255, 255 /
      DATA CODE(118),KNTM(118),KNTR(118) / 'G101.1', 255, 255 /
      DATA CODE(119),KNTM(119),KNTR(119) / 'G101.2', 255, 255 /
      DATA CODE(120),KNTM(120),KNTR(120) / 'G105.1', 255, 255 /
      DATA CODE(121),KNTM(121),KNTR(121) / 'G106.1', 255, 255 /
      DATA CODE(122),KNTM(122),KNTR(122) / 'G106.2', 255, 255 /
      DATA CODE(123),KNTM(123),KNTR(123) / 'G116.1', 255, 255 /
      DATA CODE(124),KNTM(124),KNTR(124) / 'G116.2', 255, 255 /
      DATA CODE(125),KNTM(125),KNTR(125) / 'H101.0', 255, 255 /
      DATA CODE(126),KNTM(126),KNTR(126) / 'H101.1', 255, 255 /
      DATA CODE(127),KNTM(127),KNTR(127) / 'H101.2', 255, 255 /
      DATA CODE(128),KNTM(128),KNTR(128) / 'H301.1', 255, 255 /
      DATA CODE(129),KNTM(129),KNTR(129) / 'U501.1', 255, 255 /
      DATA CODE(130),KNTM(130),KNTR(130) / 'V202.1', 255, 255 /
      DATA CODE(131),KNTM(131),KNTR(131) / 'V202.2', 255, 255 /
      DATA CODE(132),KNTM(132),KNTR(132) / 'V202.3', 255, 255 /

      ILG=NLG
      L=0
      IF(ERC .NE. ' ') THEN
       DO 10 L = 1,6
       IF(ERC(1:L) .EQ. ERC) GOTO 12
   10  CONTINUE
   12  CONTINUE
      ENDIF
      DO 14 I = 1,KTE
      IF(L .EQ. 0 .OR. CODE(I)(1:L) .EQ. ERC(1:L)) THEN
       IF(MXM .GE. 0) KNTM(I)=MXM
       IF(MXR .GE. 0) KNTR(I)=MXR
      ENDIF
   14 CONTINUE
      RETURN

      ENTRY MTLMTR(ERC,MLG,LMF,LRF)

      MLG=ILG
      DO 20 I = 1,KTE
      IF(ERC .EQ. CODE(I))  GOTO 21
   20 CONTINUE
      WRITE(*,100) ERC
      CALL ABEND
      RETURN

   21 LMF=KNTM(I) .GE. 1
      LRF=KNTR(I) .GE. 1
      IF(LMF .AND. KNTM(I) .LT. 255)  KNTM(I)=KNTM(I)-1
      IF(LRF .AND. KNTR(I) .LT. 255)  KNTR(I)=KNTR(I)-1
      IF(.NOT.LRF) THEN
       IF(ILG .LT. 1) WRITE(  *,101) CODE(I)
       IF(ILG .GE. 1) WRITE(ILG,101) CODE(I)
      ENDIF
      RETURN
  100 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR N002: ',
     1'ERROR CODE ',A6,' NOT RECOGNIZED BY ERROR MONITOR. RUN ABORTED.')
  101 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR NOO2.1: ',
     1'RUN TERMINATED BY LIBRARY ERROR CONDITION ',A6)
      END
*CMZ :          28/08/2014  12.54.56  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX r1dp.F

# 1 "r1dp.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "r1dp.F"
*
* $Id: r1dp.F,v 1.1.1.1 1996/04/01 15:01:56 mclareni Exp $
*
* $Log: r1dp.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:56 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "r1dp.F" 2

      FUNCTION C309R1(X,ETA,ZL,PM,EPS,LIMIT,ERR,NPQ,ACC8,ACCH,
     1 LPR,ACCUR,DELL)
C
C (omega) (omega)
C *** Evaluate CF2 = p + PM.q = H (ETA,X)' / H   (ETA,X)
C ZL ZL
C where PM = omega.i
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      LOGICAL LPR
      DOUBLE PRECISION EPS,ERR,ACC8,ACCH,ACCUR,TA,RK
      DOUBLE PRECISION ABSC,HALF

      PARAMETER(HALF = 1D0/2D0)




      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))

      TA=LIMIT+LIMIT
      ETAP=ETA*PM
      XI=1/X
      WI=ETAP+ETAP
      RK=0
      PQ=(1-ETA*XI)*PM
      AA=-(ETA*ETA+ZL*ZL+ZL)+ETAP
      BB=2*(X-ETA+PM)
      RL=XI*PM
      IF(ABSC(BB) .LT. ACCH) THEN
       RL=RL*AA/(AA+RK+WI)
       PQ=PQ+RL*(BB+PM+PM)
       AA=AA+2*(RK+1+WI)
       BB=BB+4*PM
       RK=RK+4
      END IF
      DD=1/BB
      DL=AA*DD*RL
   10 PQ=PQ+DL
      RK=RK+2
      AA=AA+RK+WI
      BB=BB+PM+PM
      DD=1/(AA*DD+BB)
      DL=DL*(BB*DD-1)
      ERR=ABSC(DL)/ABSC(PQ)
      IF(ERR .GE. MAX(EPS,ACC8*RK*HALF) .AND. RK .LE. TA) GO TO 10
C
      NPQ=HALF*RK
      C309R1=PQ+DL
      IF(LPR .AND. NPQ .GE. LIMIT-1 .AND. ERR .GT. ACCUR)
     1 WRITE(6,1000) INT(DIMAG(PM)),NPQ,ERR,ZL+DELL
      RETURN
 1000 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     2 'CF2(',I2,') NOT CONVERGED FULLY IN ',I7,' ITERATIONS'/1X,27X,
     3 'ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL = ',2F8.3)
      END
*CMZ :          28/08/2014  12.03.27  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
*
* $Id: r2dp.F,v 1.1.1.1 1996/04/01 15:01:56 mclareni Exp $
*
* $Log: r2dp.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:56  mclareni
* Mathlib gen
*
*

cmsh #include "gen/pilot.h"
#define CERNLIB_DOUBLE
#define CERNLIB_UNIX

#if defined(CERNLIB_DOUBLE)
      FUNCTION C309R2(X,ETA,ZL,P,EPS,LIMIT,KIND,ERR,NITS,
     1                FPMAX,ACC8,ACC16)
C
C *** evaluate the HYPERGEOMETRIC FUNCTION
C                                        i
C            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i! )
C           1 1              i       i            i
C
C     to accuracy EPS with at most LIMIT terms.
C  If KIND = 0 : using extended precision but real arithmetic only,
C            1 : using normal precision in complex arithmetic,
C   or       2 : using normal complex arithmetic, but with WDIGAM factor
C
C  where AA, BB, and Z are defined below
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX*16 X,ETA,ZL,P,AA,BB,Z,C309R2,WDIGAM
      COMPLEX*16 DD,G,F,AI,BI,T
#if (!defined(CERNLIB_UNIX))&&(!defined(CERNLIB_QMALPH))
      REAL*16 AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN
#endif
#if defined(CERNLIB_UNIX)||defined(CERNLIB_QMALPH)
      DOUBLE PRECISION AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN
#endif

      PARAMETER(TBBB = 3D0/2D0)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif

      ABSC(AA)=ABS(DREAL(AA))+ABS(DIMAG(AA))
      NINTC(AA)=NINT(DREAL(AA))

      C309R2 = 0.

      AA=ZL+1-ETA*P
      BB=2*ZL+2
      Z=2*P*X
      IF(DREAL(BB) .LE. 0 .AND. ABS(BB-NINTC(BB)) .LT.
     1 SQRT(SQRT(ACC8))**3 .AND. DREAL(BB)+LIMIT .GE. TBBB) THEN
       NITS=-1
       RETURN
      END IF
      IF(LIMIT .LE. 0) THEN
       C309R2=0
       ERR=0
       NITS=1
       RETURN
      END IF
      TA=1
      RK=1
      IF(KIND .LE. 0 .AND. ABSC(Z)*ABSC(AA) .GT. ABSC(BB)*1.0) THEN
       DR=1
       DI=0
       GR=1
       GI=0
       AR=DREAL(AA)
       BR=DREAL(BB)
       FI=0
       DO 20 I = 2,LIMIT
       FI1=FI+1
       TR=BR*FI1
       TI=DIMAG(BB)*FI1
       DEN=1/(TR*TR+TI*TI)
       UR=(AR*TR+DIMAG(AA)*TI)*DEN
       UI=(DIMAG(AA)*TR-AR*TI)*DEN
       TR=UR*GR-UI*GI
       TI=UR*GI+UI*GR
       GR=DREAL(Z)*TR-DIMAG(Z)*TI
       GI=DREAL(Z)*TI+DIMAG(Z)*TR
       DR=DR+GR
       DI=DI+GI
       ERR=ABS(GR)+ABS(GI)
       IF(ERR .GT. FPMAX) GO TO 60
       RK=ABS(DR)+ABS(DI)
       TA=MAX(TA,RK)
       IF(ERR .LT. RK*EPS .OR. I .GE. 4 .AND. ERR .LT. ACC16) GO TO 30
       FI=FI1
       AR=AR+1
   20  BR=BR+1
C
   30  C309R2=DR+(0,1)*DI
       ERR=ACC16*TA/RK
      ELSE
C
C*    If REAL*16 arithmetic is not available, (or already using it!),
C*    then use KIND > 0
C
       G=1
       F=1
       IF(KIND .GE. 2) F=WDIGAM(AA)-WDIGAM(BB)-WDIGAM(G)
       DD=F
       DO 40 I = 2,LIMIT
       AI=AA+(I-2)
       BI=BB+(I-2)
       R=I-1
       G=G*Z*AI/(BI*R)
C
C                       multiply by (psi(a+r)-psi(b+r)-psi(1+r))
C
       IF(KIND .EQ. 2) F=F+1/AI-1/BI-1/R
       T=G*F
       DD=DD+T
       ERR=ABSC(T)
       IF(ERR .GT. FPMAX) GO TO 60
       RK=ABSC(DD)
       TA=MAX(TA,RK)
       IF(ERR .LT. RK*EPS .OR. ERR .LT. ACC8 .AND. I .GE. 4) GO TO 50
   40  CONTINUE

   50  ERR=ACC8*TA/RK
       C309R2=DD
      END IF
   60 NITS=I
      RETURN
      END
#endif
*CMZ :          28/08/2014  12.45.09  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX r3dp.F

# 1 "r3dp.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "r3dp.F"
*
* $Id: r3dp.F,v 1.1.1.1 1996/04/01 15:01:56 mclareni Exp $
*
* $Log: r3dp.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:56 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "r3dp.F" 2

      FUNCTION C309R3(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)
C
C evaluate the HYPERGEOMETRIC FUNCTION
C i
C F (AA,BB;;Z) = SUM (AA) (BB) Z / i!
C 2 0 i i i
C
C to accuracy EPS with at most JMAX terms.
C
C if the terms start diverging,
C the corresponding continued fraction is found by RCF
C & evaluated progressively by Steed's method to obtain convergence.
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION X(JMAX,4)
      LOGICAL FINITE
      DOUBLE PRECISION EP,EPS,AT,ATL,ABSC,RE,FPMAX




      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))
      NINTC(W)=NINT(DREAL(W))
C
      RE=0
      X(1,1)=1
      SUM=X(1,1)
      ATL=ABSC(X(1,1))
      F=SUM
      D=1
      DF=SUM
      J=0
      EP=EPS*(10*JMAX)
      MA=-NINTC(AA)
      MB=-NINTC(BB)
      FINITE=ABS(ABS(DREAL(AA))-MA) .LT. EP .AND. ABS(DIMAG(AA)) .LT. EP
     1 .OR. ABS(ABS(DREAL(BB))-MB) .LT. EP .AND. ABS(DIMAG(BB)) .LT. EP
      IMAX=JMAX
      IF(FINITE .AND. MA .GE. 0) IMAX=MIN(MA+1,IMAX)
      IF(FINITE .AND. MB .GE. 0) IMAX=MIN(MB+1,IMAX)
      DO 10 I = 2,IMAX
      X(I,1)=X(I-1,1)*Z*(AA+I-2)*(BB+I-2)/(I-1)
      IF(ABSC(X(I,1)) .GT. FPMAX) THEN
       N=0
       C309R3=SUM
       IF(.NOT.FINITE) RE=AT/ABSC(SUM)
       RETURN
      END IF
      AT=ABSC(X(I,1))
      IF(J .EQ. 0) THEN
       SUM=SUM+X(I,1)
       IF(AT .LT. ABSC(SUM)*EPS) THEN
        N=I
        C309R3=SUM
        IF(.NOT.FINITE) RE=AT/ABSC(SUM)
        RETURN
       END IF
      END IF
      IF(FINITE) GO TO 10
      IF(J .GT. 0 .OR. AT .GT. ATL .OR. I .GE. JMAX-2) J=J+1
      IF(J .EQ. 0) GO TO 10
      CALL C309R7(X(1,1),X(1,2),J,I,X(1,3),EPS)
      IF(I .LT. 0) THEN
       N=0
       C309R3=SUM
       IF(.NOT.FINITE) RE=AT/ABSC(SUM)
       RETURN
      END IF
      DO 50 K = MAX(J,2),I
      D=1/(D*X(K,2)+1)
      DF=DF*D-DF
      F=F+DF
      IF(ABSC(DF) .LT. ABSC(F)*EPS .OR.
     1 DF .EQ. 0 .AND. F .EQ. 0 .AND. I .GE. 4) THEN
       N=K
       C309R3=F
       RE=ABSC(DF)/ABSC(F)
       RETURN
      END IF
   50 CONTINUE
      J=I
   10 ATL=AT
      IF(.NOT.FINITE) I=-JMAX
      N=I
      C309R3=SUM
      IF(.NOT.FINITE) RE=AT/ABSC(SUM)
      RETURN
      END
*CMZ :          28/08/2014  12.58.19  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX r4dp.F

# 1 "r4dp.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "r4dp.F"
*
* $Id: r4dp.F,v 1.1.1.1 1996/04/01 15:01:57 mclareni Exp $
*
* $Log: r4dp.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:57 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "r4dp.F" 2

      FUNCTION C309R4(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     1 FPMIN,FPMAX,LPR)
C
C *** Evaluate CF1 = F'(ZL,ETA,X)/F(ZL,ETA,X)    (REAL)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LPR,ETANE0

      FCL=1
      XI=1/X
      PK=ZL+1
      PX=PK+LIMIT
      EK=ETA/PK
      F=EK+PK*XI
      IF(ABS(F) .LT. FPMIN) F=FPMIN
      D=0
      C=F
      SMALL=SQRT(FPMIN)
      RK2=1+EK*EK
C
C *** begin CF1 loop on PK = k = lambda + 1
C
   10 PK1=PK+1
      TPK1=PK+PK1
      IF(ETANE0) THEN
       EK=ETA/PK
       RK2=1+EK*EK
       TK=TPK1*(XI+EK/PK1)
      ELSE
       TK=TPK1*XI
      END IF
      C=TK-RK2/C
      D=TK-RK2*D
      IF(ABS(C) .LT. FPMIN) C=FPMIN
      IF(ABS(D) .LT. FPMIN) D=FPMIN
      D=1/D
      DF=D*C
      F=F*DF
      FCL=FCL*D*TPK1*XI
      IF(ABS(FCL) .LT. SMALL) FCL=FCL/SMALL
      IF(ABS(FCL) .GT. FPMAX) FCL=FCL*FPMIN
      PK=PK1
      IF(PK .LE. PX) THEN
       IF(ABS(DF-1) .GE. EPS) GO TO 10
       NFP=PK-ZL-1
       ERR=EPS*SQRT(REAL(NFP))
       C309R4=F
      ELSE
       IF(LPR) WRITE (6,1000) LIMIT,ABS(X)
       ERR=2
      END IF
      RETURN
 1000 FORMAT(1X,'***** CERN C309 WCLBES ... CF1 (REAL) HAS FAILED ',
     1'TO CONVERGE AFTER',I10,' ITERATIONS AS ABS(X) =',F15.0)
      END
*CMZ :          28/08/2014  12.57.13  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX r5dp.F

# 1 "r5dp.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "r5dp.F"
*
* $Id: r5dp.F,v 1.1.1.1 1996/04/01 15:01:57 mclareni Exp $
*
* $Log: r5dp.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:57 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "r5dp.F" 2

      FUNCTION C309R5(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     1 FPMIN,FPMAX,LPR)
C
C *** Evaluate CF1 = F'(ZL,ETA,X)/F(ZL,ETA,X)  (COMPLEX)
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      LOGICAL LPR,ETANE0
      DOUBLE PRECISION EPS,ERR,FPMIN,FPMAX,ABSC,SMALL,PX





      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))

      FCL=1
      XI=1/X
      PK=ZL+1
      PX=PK+LIMIT
      EK=ETA/PK
      F=EK+PK*XI
      IF(ABSC(F) .LT. FPMIN) F=FPMIN
      D=0
      C=F
      SMALL=SQRT(FPMIN)
      RK2=1+EK*EK
C
C *** begin CF1 loop on PK = k = lambda + 1
C
   10 PK1=PK+1
      TPK1=PK+PK1
      IF(ETANE0) THEN
       EK=ETA/PK
       RK2=1+EK*EK
       TK=TPK1*(XI+EK/PK1)
      ELSE
       TK=TPK1*XI
      END IF
      C=TK-RK2/C
      D=TK-RK2*D
      IF(ABSC(C) .LT. FPMIN) C=FPMIN
      IF(ABSC(D) .LT. FPMIN) D=FPMIN
      D=1/D
      DF=D*C
      F=F*DF
      FCL=FCL*D*TPK1*XI
      IF(ABSC(FCL) .LT. SMALL) FCL=FCL/SMALL
      IF(ABSC(FCL) .GT. FPMAX) FCL=FCL*FPMIN
      PK=PK1
      IF(DREAL(PK) .LE. PX) THEN
       IF(ABSC(DF-1) .GE. EPS) GO TO 10
       NFP=PK-ZL-1
       ERR=EPS*SQRT(REAL(NFP))
       C309R5=F
      ELSE
       IF(LPR) WRITE (6,1000) LIMIT,ABS(X)
       ERR=2
      END IF
      RETURN
 1000 FORMAT(1X,'***** CERN C309 WCLBES ... CF1 (COMPLEX) HAS FAILED ',
     1'TO CONVERGE AFTER',I10,' ITERATIONS AS ABS(X) =',F15.0)
      END
*CMZ :          28/08/2014  12.59.39  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh: Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX r6dp.F

# 1 "r6dp.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "r6dp.F"
*
* $Id: r6dp.F,v 1.1.1.1 1996/04/01 15:01:57 mclareni Exp $
*
* $Log: r6dp.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:57 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "r6dp.F" 2

      FUNCTION C309R6(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
C
C evaluate the ASYMPTOTIC EXPANSION for the
C LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION
C
C *** CF1A = f = F'(XL,ETA,RHO)/F(XL,ETA,RHO)
C
C that is valid for DREAL(RHO)>0, and best for RHO >> ETA**2, XL,
C and is derived from the 2F0 expansions for H+ and H-
C e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)
C Some lines of this subprogram are for convenience copied from
C Takemasa, Tamura & Wolter CPC 17 (1979) 351.
C
C Evaluate to accuracy EPS with at most NMAX terms.
C
C If the terms start diverging,
C the corresponding continued fraction is found by RCF
C & evaluated progressively by Steed's method to obtain convergence.
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION XX(2,NMAX),G(NMAX),C(NMAX)
      DOUBLE PRECISION RE,EPS,T1,T2,T3,AT,ATL,ABSC,FPMAX
      DOUBLE PRECISION HPI

      PARAMETER(HPI = 1.57079 63267 94896 619D0)




      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))
C
      T1=SIN(DREAL(PSI))
      T2=COS(DREAL(PSI))
      ATL=TANH(DIMAG(PSI))

C GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGN

      COSL=DCMPLX(T2,-T1*ATL)
      TANL=DCMPLX(T1,T2*ATL)/COSL
      RE=0
      XLL1=XL*XL+XL
      ETASQ=ETA*ETA
      SL1=1
      SL=SL1
      SC1=0
      SC=SC1
      TL1=SC
      TL=TL1
      TC1=1-ETA/RHO
      TC=TC1
      FCL=TL+SL*TANL
      G(1)=(TC+SC*TANL)/FCL
      GLAST=G(1)
      ATL=ABSC(GLAST)
      F=GLAST
      D=1
      DF=GLAST
      J=0
      DO 10 N = 2,NMAX
      T1=N-1
      T2=2*T1-1
      T3=T1*T1-T1
      DENOM=2*RHO*T1
      C1=(ETA*T2)/DENOM
      C2=(ETASQ+XLL1-T3)/DENOM
      SL2=C1*SL1-C2*TL1
      TL2=C1*TL1+C2*SL1
      SC2=C1*SC1-C2*TC1-SL2/RHO
      TC2=C1*TC1+C2*SC1-TL2/RHO
      SL=SL+SL2
      TL=TL+TL2
      SC=SC+SC2
      TC=TC+TC2
      SL1=SL2
      TL1=TL2
      SC1=SC2
      TC1=TC2
      FCL=TL+SL*TANL
      IF(ABSC(FCL) .GT. FPMAX .OR. ABSC(FCL) .LT. 1./FPMAX) THEN
       C309R6=G(1)
       FCL=1
       RE=1
       NUSED=0
       RETURN
      END IF
      GSUM=(TC+SC*TANL)/FCL
      G(N)=GSUM-GLAST
      GLAST=GSUM
      AT=ABSC(G(N))
      IF(AT .LT. ABSC(GSUM)*EPS) THEN
       FCL=FCL*COSL
       C309R6=GSUM
       RE=AT/ABSC(GSUM)
       NUSED=N
       RETURN
      END IF
      IF(J .GT. 0 .OR. AT .GT. ATL .OR. N .GE. NMAX-2) J=J+1
      IF(J .EQ. 0) GO TO 10
      CALL C309R7(G,C,J,N,XX,EPS)
      IF(N .LT. 0) THEN
       C309R6=G(1)
       FCL=1
       RE=1
       NUSED=0
       RETURN
      END IF
      DO 60 K = MAX(J,2),N
      D=1/(D*C(K)+1)
      DF=DF*D-DF
      F=F+DF
      IF(ABSC(DF) .LT. ABSC(F)*EPS .OR.
     1 DF .EQ. 0 .AND. F .EQ. 0 .AND. N .GE. 4) THEN
       C309R6=F
       FCL=FCL*COSL
       RE=ABSC(DF)/ABSC(F)
       NUSED=K
       RETURN
      END IF
   60 CONTINUE
      J=N
   10 ATL=AT
      C309R6=F
      FCL=FCL*COSL
      RE=ABSC(DF)/ABSC(F)
      NUSED=-NMAX
      RETURN
      END
*CMZ :          28/08/2014  12.53.19  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX r7dp.F

# 1 "r7dp.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "r7dp.F"
*
* $Id: r7dp.F,v 1.1.1.1 1996/04/01 15:01:57 mclareni Exp $
*
* $Log: r7dp.F,v $
* Revision 1.1.1.1 1996/04/01 15:01:57 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "r7dp.F" 2

      SUBROUTINE C309R7(A,B,IBEG,INUM,XX,EPS)
C
C*******************************************************************
C
C RCF converts polynomial A to the corresponding continued
C fraction, in 'normal' form with coefficients B
C by the 'P algorithm' of Patry & Gupta
C
C A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)
C
C B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)
C
C data:
C A vector A(k), k=1,INUM input
C B vector B(k), k=IBEG,INUM output
C IBEG order of first coef. calc. input
C INUM order of A, even or odd input
C XX auxiliary vector of length .ge. length of vector B
C caller provides space for A,B,XX
C Note that neither of the first two terms A(1) A(2) should be zero
C & the user can start the calculation with any value of
C IBEG provided the c.f. coefs have been already
C calculated up to INUM = IBEG-1
C & the method breaks down as soon as the absolute value
C of a c.f. coef. is less than EPS. At the time of the
C break up INUM has been replaced by minus times the number
C of this coefficient.
C algorithm: J. Patry & S. Gupta, EIR-Bericht 247, November 1973
C Eidg. Institut fur Reaktorforschung
C Wuerenlingen, Switzerland
C see also: Haenggi, Roesel & Trautmann,
C J. Comput. Phys., v. 137, (1980) 242-258
C note: restart procedure modified by I.J.Thompson
C
C*******************************************************************
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION A(100),B(100),XX(2,100)
      LOGICAL EVEN
      DOUBLE PRECISION EPS

      IBN=INUM
      IF(IBEG .GT. 4) GO TO 50
      IF(IBEG .EQ. 4) GO TO 20
      B(1)=A(1)
      IF(IBN .GE. 2) B(2)=-A(2)/A(1)
      IF(IBN .LT. 3) RETURN
      X0=A(3)/A(2)
      XX(2,1)=B(2)
      XX(1,1)=-X0
      XX(1,2)=0
      B(3)=-X0-B(2)
      X0=-B(3)*A(2)
      M=3
      MP12=2
      EVEN=.TRUE.
      IF(IBN .LE. 3) RETURN
   20 IF(ABS(B(3)) .LT. EPS*ABS(X0)) THEN
       INUM=-M
       RETURN
      END IF
      M=4
   30 X1=A(M)
      M2M1=MP12
      MP12=M2M1+1
      IF(EVEN) MP12=M2M1
      DO 40 K = 2,MP12
   40 X1=X1+A(M-K+1)*XX(1,K-1)
      B(M)=-X1/X0
      IF(M .GE. IBN) RETURN
   50 IF(ABS(B(M)) .LT. EPS*ABS(X0)) THEN
       INUM=-M
       RETURN
      END IF
      DO 60 K = M2M1,2,-1
   60 XX(2,K)=XX(1,K)+B(M)*XX(2,K-1)
      XX(2,1)=XX(1,1)+B(M)
      DO 70 K = 1,M2M1
      X0=XX(2,K)
      XX(2,K)=XX(1,K)
   70 XX(1,K)=X0
      X0=X1
      XX(1,M2M1+1)=0
      M=M+1
      EVEN=.NOT.EVEN
      GO TO 30
      END
*CMZ :          28/08/2014  13.48.08  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX radapt.F

# 1 "radapt.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "radapt.F"
*
* $Id: radapt.F,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: radapt.F,v $
* Revision 1.1.1.1 1996/04/01 15:02:13 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "radapt.F" 2
      SUBROUTINE RADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)

C RES = Estimated Integral of F from A to B,
C ERR = Estimated absolute error on RES.
C NSEG specifies how the adaptation is to be done:
C =0 means use previous binning,
C =1 means fully automatic, adapt until tolerance attained.
C =n>1 means first split interval into n equal segments,
C then adapt as necessary to attain tolerance.
C The specified tolerances are:
C relative: RELTOL ; absolute: ABSTOL.
C It stops when one OR the other is satisfied, or number of
C segments exceeds NDIM. Either TOLA or TOLR (but not both!)
C can be set to zero, in which case only the other is used.

      DOUBLE PRECISION TVALS,TERSS
      EXTERNAL F

      PARAMETER (NDIM=100)
      PARAMETER (R1 = 1, HF = R1/2)

      DIMENSION XLO(NDIM),XHI(NDIM),TVAL(NDIM),TERS(NDIM)
      SAVE XLO,XHI,TVAL,TERS,NTER
      DATA NTER /0/

      IF(NSEG .LE. 0) THEN
       IF(NTER .EQ. 0) THEN
        NSEGD=1
        GO TO 2
       ENDIF
       TVALS=0
       TERSS=0
       DO 1 I = 1,NTER
       CALL RGS56P(F,XLO(I),XHI(I),TVAL(I),TE)
       TERS(I)=TE**2
       TVALS=TVALS+TVAL(I)
       TERSS=TERSS+TERS(I)
    1 CONTINUE
       ROOT= SQRT(2*TERSS)
       GO TO 9
      ENDIF
      NSEGD=MIN(NSEG,NDIM)
    2 XHIB=A
      BIN=(B-A)/NSEGD
      DO 3 I = 1,NSEGD
      XLO(I)=XHIB
      XLOB=XLO(I)
      XHI(I)=XHIB+BIN
      IF(I .EQ. NSEGD) XHI(I)=B
      XHIB=XHI(I)
      CALL RGS56P(F,XLOB,XHIB,TVAL(I),TE)
      TERS(I)=TE**2
    3 CONTINUE
      NTER=NSEGD
      DO 4 ITER = 1,NDIM
      TVALS=TVAL(1)
      TERSS=TERS(1)
      DO 5 I = 2,NTER
      TVALS=TVALS+TVAL(I)
      TERSS=TERSS+TERS(I)
    5 CONTINUE
      ROOT= SQRT(2*TERSS)
      IF(ROOT .LE. ABSTOL .OR. ROOT .LE. RELTOL*ABS(TVALS)) GO TO 9
      IF(NTER .EQ. NDIM) GO TO 9
      BIGE=TERS(1)
      IBIG=1
      DO 6 I = 2,NTER
      IF(TERS(I) .GT. BIGE) THEN
       BIGE=TERS(I)
       IBIG=I
      ENDIF
    6 CONTINUE
      NTER=NTER+1
      XHI(NTER)=XHI(IBIG)
      XNEW=HF*(XLO(IBIG)+XHI(IBIG))
      XHI(IBIG)=XNEW
      XLO(NTER)=XNEW
      CALL RGS56P(F,XLO(IBIG),XHI(IBIG),TVAL(IBIG),TE)
      TERS(IBIG)=TE**2
      CALL RGS56P(F,XLO(NTER),XHI(NTER),TVAL(NTER),TE)
      TERS(NTER)=TE**2
    4 CONTINUE
    9 RES=TVALS
      ERR=ROOT
      RETURN
      END
*CMZ :          10/04/2019  11.55.13  by  Michael Scheer
*-- Author :    Michael Scheer   02/05/2017
cmsh # 1 "ranf.F"
# 1 "<built-in>"
# 1 "<command-line>"
cmsh # 1 "ranf.F"
*
* $Id: ranf.F,v 1.1.1.1 1996/02/15 17:49:05 mclareni Exp $
*
* $Log: ranf.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:05  mclareni
* Kernlib
*
*

# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 21 "/usr/include/kernnum/pilot.h" 3 4

# 33 "/usr/include/kernnum/pilot.h" 3 4

# 45 "/usr/include/kernnum/pilot.h" 3 4

cmsh # 10 "ranf.F" 2
          REAL FUNCTION RANF()
cmsh          DOUBLE PRECISION    DRANF,    G900GT,   G900ST
          DOUBLE PRECISION    G900GT,   G900ST
          DOUBLE PRECISION    DS(2),    DM(2),    DSEED
          DOUBLE PRECISION    DX24,     DX48
          DOUBLE PRECISION    DL,       DC,       DU,       DR
          LOGICAL             SINGLE
          DATA      DS     /  1665 1885.D0, 286 8876.D0  /
          DATA      DM     /  1518 4245.D0, 265 1554.D0  /
          DATA      DX24   /  1677 7216.D0  /
          DATA      DX48   /  281 4749 7671 0656.D0  /
          SINGLE  =  .TRUE.
          GOTO 10
cmsh          ENTRY DRANF()
          SINGLE  =  .FALSE.
  10      DL  =  DS(1) * DM(1)
          DC  =  DINT(DL/DX24)
          DL  =  DL - DC*DX24
          DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
          DS(2)  =  DU - DINT(DU/DX24)*DX24
          DS(1)  =  DL
          DR     =  (DS(2)*DX24 + DS(1)) / DX48
          IF(SINGLE)  THEN
             RANF  =  SNGL(DR)
           ELSE
             print*,"*** Error in RANF: Call to DRANF ***"
cmsh             DRANF  =  DR
          ENDIF
          RETURN
          ENTRY G900GT()
          G900GT  =  DS(2)*DX24 + DS(1)
          RETURN
          ENTRY G900ST(DSEED)
          DS(2)  =  DINT(DSEED/DX24)
          DS(1)  =  DSEED - DS(2)*DX24
          G900ST =  DS(1)
          RETURN
          END
*CMZ :          03/12/2024  10.54.19  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX rfft.F

# 1 "rfft.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "rfft.F"
*
* $Id: rfft.F,v 1.1.1.1 1996/02/15 17:48:48 mclareni Exp $
*
* $Log: rfft.F,v $
* Revision 1.1.1.1 1996/02/15 17:48:48 mclareni
* Kernlib
*
*
# 1 "/usr/include/kernnum/pilot.h" 1 3 4
# 10 "rfft.F" 2
      SUBROUTINE RFFT(A,MSIGN)
cmsh      COMPLEX A(1),T1,T2,U,W
      COMPLEX A(2),T1,T2,U,W
      IF(MSIGN.EQ.0) RETURN
      M=IABS(MSIGN)-1
      N=2**M
      U=(0.,1.)
      IF(MSIGN.GT.0) GO TO 2
      CALL CFFT(A,-M)
      F=.25/N
      DO 1 I=1,N
 1    A(I)=F*A(I)
      A(N+1)=A(1)
      U=CONJG(U)
 2    ANGL=3.1415926535898/ISIGN(N,MSIGN)
      W=CMPLX(COS(ANGL),SIN(ANGL))
      N2=N+2
      N1=N2/2
      DO 3 J=1,N1
      K=N2-J
      T1= A(J)+CONJG(A(K))
      T2=(A(J)-CONJG(A(K)))*U
      A(J)= T1+T2
      A(K)=CONJG(T1-T2)
 3    U=U*W
      IF(MSIGN.GT.0) CALL CFFT(A, M)
      RETURN
      END
*CMZ :          25/08/2014  16.15.57  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  14.20.46  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
*
* $Id: rfstft.F,v 1.1 1996/04/17 12:32:04 mclareni Exp $
*
* $Log: rfstft.F,v $
* Revision 1.1  1996/04/17 12:32:04  mclareni
* Add d/rfstft.F (D705) and to Imakefile. cfstft.F becomes D706.
* In tests, add d705m.F for rfstft and d706m.F for cfstft and the corresponding
* additions to main.F and Imakefile.
*
*
      SUBROUTINE RFSTFT(MS,A)

      COMPLEX A(0:*),T,T1,T2,U,W

cmsh      PARAMETER (PI = 3.14159 26535 89793D0)
      PARAMETER (PI = 3.141592653589793)

      IF(MS .EQ. 0) THEN
       A(0)=REAL(A(0))
       RETURN
      ENDIF
      M=ABS(MS)-1
      N=2**M
      U=(0.,1.)
      IF(MS .LT. 0) THEN
       CALL CFSTFT(-M,A)
       F=0.25/N
       DO 1 I = 0,N-1
    1  A(I)=F*A(I)
       A(N)=A(0)
       U=CONJG(U)
      ENDIF
cmsh    2 PHI=PI/SIGN(N,MS)
      PHI=PI/SIGN(N,MS)
      W=CMPLX(COS(PHI),SIN(PHI))
      DO 3 J = 0,N/2
      T=CONJG(A(N-J))
      T1=A(J)+T
      T2=(A(J)-T)*U
      A(J)=T1+T2
      A(N-J)=CONJG(T1-T2)
    3 U=U*W
      IF(MS .GT. 0) CALL CFSTFT(M,A)
      RETURN
      END
*CMZ :          28/08/2014  13.55.30  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX rgs56p.F

# 1 "rgs56p.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "rgs56p.F"
*
* $Id: rgs56p.F,v 1.1.1.1 1996/04/01 15:02:14 mclareni Exp $
*
* $Log: rgs56p.F,v $
* Revision 1.1.1.1 1996/04/01 15:02:14 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "rgs56p.F" 2
      SUBROUTINE RGS56P(F,A,B,RES,ERR)
      DOUBLE PRECISION E5,E6

      PARAMETER (R1 = 1, HF = R1/2)
cmsh      DIMENSION X5(5),W5(5),X6(6),W6(6)
      double precision  X5(5),W5(5),X6(6),W6(6) !msh

      DATA (X5(I),W5(I),I=1,5)
     1/4.6910077030668004D-02, 1.1846344252809454D-01,
     2 2.3076534494715846D-01, 2.3931433524968324D-01,
     3 5.0000000000000000D-01, 2.8444444444444444D-01,
     4 7.6923465505284154D-01, 2.3931433524968324D-01,
     5 9.5308992296933200D-01, 1.1846344252809454D-01/

      DATA (X6(I),W6(I),I=1,6)
     1/3.3765242898423989D-02, 8.5662246189585178D-02,
     2 1.6939530676686775D-01, 1.8038078652406930D-01,
     3 3.8069040695840155D-01, 2.3395696728634552D-01,
     4 6.1930959304159845D-01, 2.3395696728634552D-01,
     5 8.3060469323313225D-01, 1.8038078652406930D-01,
     6 9.6623475710157601D-01, 8.5662246189585178D-02/

      RANG=B-A
      E5=0
      E6=0
      DO 1 I = 1,5
      E5=E5+W5(I)*F(A+RANG*X5(I))
      E6=E6+W6(I)*F(A+RANG*X6(I))
    1 CONTINUE
      E6=E6+W6(6)*F(A+RANG*X6(6))
      RES=HF*(E6+E5)*RANG
      ERR=ABS((E6-E5)*RANG)
      RETURN
      END
*CMZ :          16/04/2025  09.04.52  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2025
# 0 "/home/scheer/cern/github/cernlib_2005/2005/src/mathlib/gen/c/sinint64.F"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "/home/scheer/cern/github/cernlib_2005/2005/src/mathlib/gen/c/sinint64.F"
*
* $Id: sinint64.F,v 1.1.1.1 1996/04/01 15:02:06 mclareni Exp $
*
* $Log: sinint64.F,v $
* Revision 1.1.1.1 1996/04/01 15:02:06 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "/home/scheer/cern/github/cernlib_2005/2005/src/mathlib/gen/c/sinint64.F" 2




      FUNCTION DSININ(X)

# 1 "/usr/include/gen/imp64.inc" 1 3 4

# 1 "/usr/include/gen/imp64.inc" 3 4
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1 1996/04/01 15:02:59 mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
# 17 "/home/scheer/cern/github/cernlib_2005/2005/src/mathlib/gen/c/sinint64.F" 2

# 17 "/home/scheer/cern/github/cernlib_2005/2005/src/mathlib/gen/c/sinint64.F"
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT




      PARAMETER (NAME = 'RCOSIN/DCOSIN')

      DIMENSION S(0:15),C(0:15),P(0:28),Q(0:24)

      PARAMETER (Z1 = 1, R8 = Z1/8, R32 = Z1/32)

      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (CE = 0.57721 56649 01532 86D0)
      PARAMETER (PIH = PI/2)

      DATA S( 0) /+1.95222 09759 53071 08D0/
      DATA S( 1) /-0.68840 42321 25715 44D0/
      DATA S( 2) /+0.45518 55132 25584 84D0/
      DATA S( 3) /-0.18045 71236 83877 85D0/
      DATA S( 4) /+0.04104 22133 75859 24D0/
      DATA S( 5) /-0.00595 86169 55588 85D0/
      DATA S( 6) /+0.00060 01427 41414 43D0/
      DATA S( 7) /-0.00004 44708 32910 75D0/
      DATA S( 8) /+0.00000 25300 78230 75D0/
      DATA S( 9) /-0.00000 01141 30759 30D0/
      DATA S(10) /+0.00000 00041 85783 94D0/
      DATA S(11) /-0.00000 00001 27347 06D0/
      DATA S(12) /+0.00000 00000 03267 36D0/
      DATA S(13) /-0.00000 00000 00071 68D0/
      DATA S(14) /+0.00000 00000 00001 36D0/
      DATA S(15) /-0.00000 00000 00000 02D0/

      DATA C( 0) /+1.94054 91464 83554 93D0/
      DATA C( 1) /+0.94134 09132 86521 34D0/
      DATA C( 2) /-0.57984 50342 92992 76D0/
      DATA C( 3) /+0.30915 72011 15927 13D0/
      DATA C( 4) /-0.09161 01792 20771 34D0/
      DATA C( 5) /+0.01644 37407 51546 25D0/
      DATA C( 6) /-0.00197 13091 95216 41D0/
      DATA C( 7) /+0.00016 92538 85083 50D0/
      DATA C( 8) /-0.00001 09393 29573 11D0/
      DATA C( 9) /+0.00000 05522 38574 84D0/
      DATA C(10) /-0.00000 00223 99493 31D0/
      DATA C(11) /+0.00000 00007 46533 25D0/
      DATA C(12) /-0.00000 00000 20818 33D0/
      DATA C(13) /+0.00000 00000 00493 12D0/
      DATA C(14) /-0.00000 00000 00010 05D0/
      DATA C(15) /+0.00000 00000 00000 18D0/

      DATA P( 0) /+0.96074 78397 52035 96D0/
      DATA P( 1) /-0.03711 38962 12398 06D0/
      DATA P( 2) /+0.00194 14398 88991 90D0/
      DATA P( 3) /-0.00017 16598 84251 47D0/
      DATA P( 4) /+0.00002 11263 77532 31D0/
      DATA P( 5) /-0.00000 32716 32567 12D0/
      DATA P( 6) /+0.00000 06006 92116 15D0/
      DATA P( 7) /-0.00000 01258 67944 03D0/
      DATA P( 8) /+0.00000 00293 25634 58D0/
      DATA P( 9) /-0.00000 00074 56959 21D0/
      DATA P(10) /+0.00000 00020 41054 78D0/
      DATA P(11) /-0.00000 00005 95022 30D0/
      DATA P(12) /+0.00000 00001 83229 67D0/
      DATA P(13) /-0.00000 00000 59205 06D0/
      DATA P(14) /+0.00000 00000 19965 17D0/
      DATA P(15) /-0.00000 00000 06995 11D0/
      DATA P(16) /+0.00000 00000 02536 86D0/
      DATA P(17) /-0.00000 00000 00949 29D0/
      DATA P(18) /+0.00000 00000 00365 52D0/
      DATA P(19) /-0.00000 00000 00144 49D0/
      DATA P(20) /+0.00000 00000 00058 51D0/
      DATA P(21) /-0.00000 00000 00024 23D0/
      DATA P(22) /+0.00000 00000 00010 25D0/
      DATA P(23) /-0.00000 00000 00004 42D0/
      DATA P(24) /+0.00000 00000 00001 94D0/
      DATA P(25) /-0.00000 00000 00000 87D0/
      DATA P(26) /+0.00000 00000 00000 39D0/
      DATA P(27) /-0.00000 00000 00000 18D0/
      DATA P(28) /+0.00000 00000 00000 08D0/

      DATA Q( 0) /+0.98604 06569 62382 60D0/
      DATA Q( 1) /-0.01347 17382 08295 21D0/
      DATA Q( 2) /+0.00045 32928 41165 23D0/
      DATA Q( 3) /-0.00003 06728 86516 55D0/
      DATA Q( 4) /+0.00000 31319 91976 01D0/
      DATA Q( 5) /-0.00000 04211 01964 96D0/
      DATA Q( 6) /+0.00000 00690 72448 30D0/
      DATA Q( 7) /-0.00000 00131 83212 90D0/
      DATA Q( 8) /+0.00000 00028 36974 33D0/
      DATA Q( 9) /-0.00000 00006 73292 34D0/
      DATA Q(10) /+0.00000 00001 73396 87D0/
      DATA Q(11) /-0.00000 00000 47869 39D0/
      DATA Q(12) /+0.00000 00000 14032 35D0/
      DATA Q(13) /-0.00000 00000 04334 96D0/
      DATA Q(14) /+0.00000 00000 01402 73D0/
      DATA Q(15) /-0.00000 00000 00473 06D0/
      DATA Q(16) /+0.00000 00000 00165 58D0/
      DATA Q(17) /-0.00000 00000 00059 94D0/
      DATA Q(18) /+0.00000 00000 00022 37D0/
      DATA Q(19) /-0.00000 00000 00008 59D0/
      DATA Q(20) /+0.00000 00000 00003 38D0/
      DATA Q(21) /-0.00000 00000 00001 36D0/
      DATA Q(22) /+0.00000 00000 00000 56D0/
      DATA Q(23) /-0.00000 00000 00000 24D0/
      DATA Q(24) /+0.00000 00000 00000 10D0/

      IF(ABS(X) .LE. 8) THEN
       Y=R8*X
       H=2*Y**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 15,0,-1
       B0=S(I)+ALFA*B1-B2
       B2=B1
    1 B1=B0
       H=Y*(B0-B2)
      ELSE
       R=1/X
       H=128*R**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 28,0,-1
       B0=P(I)+ALFA*B1-B2
       B2=B1
    2 B1=B0
       PP=B0-H*B2
       B1=0
       B2=0
       DO 3 I = 24,0,-1
       B0=Q(I)+ALFA*B1-B2
       B2=B1
    3 B1=B0
       H=SIGN(PIH,X)-R*(R*PP*SIN(X)+(B0-H*B2)*COS(X))
      END IF
      GO TO 9





      ENTRY DCOSIN(X)


      IF(X .EQ. 0) THEN
       H=0
       CALL MTLPRT(NAME,'C336.1','ARGUMENT X = 0')
      ELSEIF(ABS(X) .LE. 8) THEN
       H=R32*X**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 4 I = 15,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    4 B1=B0
       H=CE+LOG(ABS(X))-B0+H*B2
      ELSE
       R=1/X
       H=128*R**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 5 I = 28,0,-1
       B0=P(I)+ALFA*B1-B2
       B2=B1
    5 B1=B0
       PP=B0-H*B2
       B1=0
       B2=0
       DO 6 I = 24,0,-1
       B0=Q(I)+ALFA*B1-B2
       B2=B1
    6 B1=B0
       H=R*((B0-H*B2)*SIN(X)-R*PP*COS(X))
      END IF




    9 DSININ=H

      RETURN
      END
*CMZ :  1.16/04 17/04/2014  12.59.42  by  Michael Scheer
*-- Author :    Michael Scheer   17/04/2014
*# 1 "tmprnt.F"
*# 1 "<command-line>"
*# 1 "tmprnt.F"
*
* $Id: tmprnt.F,v 1.1.1.1 1996/02/15 17:49:04 mclareni Exp $
*
* $Log: tmprnt.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*

*# 1 "kernnum/pilot.h" 1
*# 21 "kernnum/pilot.h"

*# 33 "kernnum/pilot.h"

*# 10 "tmprnt.F" 2
          SUBROUTINE          TMPRNT(NAME,N,IDIM,K)
          CHARACTER*6         NAME
          LOGICAL             MFLAG,    RFLAG
          IF(NAME(2:2) .EQ. 'S') THEN
             CALL KERMTR('F012.1',LGFILE,MFLAG,RFLAG)
          ELSE
             CALL KERMTR('F011.1',LGFILE,MFLAG,RFLAG)
          ENDIF
cmsh          IF(NAME(3:6) .EQ. 'FEQN') ASSIGN 1002 TO IFMT
cmsh          IF(NAME(3:6) .NE. 'FEQN') ASSIGN 1001 TO IFMT
          IF(MFLAG) THEN
             IF(LGFILE .EQ. 0) THEN
                IF(NAME(3:6) .EQ. 'FEQN') THEN
                   WRITE(*,1002) NAME, N, IDIM, K
                ELSE
                   WRITE(*,1001) NAME, N, IDIM
                ENDIF
             ELSE
                IF(NAME(3:6) .EQ. 'FEQN') THEN
                   WRITE(LGFILE,1002) NAME, N, IDIM, K
                ELSE
                   WRITE(LGFILE,1001) NAME, N, IDIM
                ENDIF
             ENDIF
          ENDIF
          IF(.NOT. RFLAG) CALL ABEND
          RETURN
1001      FORMAT(7X, 31H PARAMETER ERROR IN SUBROUTINE , A6,
     +             27H ... (N.LT.1 OR IDIM.LT.N).,
     +             5X, 3HN =, I4, 5X, 6HIDIM =, I4, 1H. )
1002      FORMAT(7X, 31H PARAMETER ERROR IN SUBROUTINE , A6,
     +             37H ... (N.LT.1 OR IDIM.LT.N OR K.LT.1).,
     +             5X, 3HN =, I4, 5X, 6HIDIM =, I4, 5X, 3HK =, I4,1H.)
          END
*CMZ :          28/08/2014  12.03.44  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014
*
* $Id: wclbes.F,v 1.1.1.1 1996/04/01 15:01:57 mclareni Exp $
*
* $Log: wclbes.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:57  mclareni
* Mathlib gen
*
*

cmsh #include "gen/pilot.h"
#define CERNLIB_DOUBLE

#if defined(CERNLIB_DOUBLE)
      SUBROUTINE WCLBES(ZZ,ETA1,ZLMIN,NL,FC,GC,FCP,GCP,SIG,KFN,MODE1,
     1                  IFAIL,IPR)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C  Original title : COULCC                                             C
C                                                                      C
C  A. R. Barnett           Manchester  March   1981                    C
C  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions C
C                                                                      C
C  The FORM (not the SUBSTANCE) of this program has been modified      C
C   by K.S. KOLBIG (CERN)    December 1987                             C
C                                                                      C
C  original program  RCWFN       in    CPC  8 (1974) 377-395           C
C                 +  RCWFF       in    CPC 11 (1976) 141-142           C
C                 +  COULFG      in    CPC 27 (1982) 147-166           C
C  description of real algorithm in    CPC 21 (1981) 297-314           C
C  description of complex algorithm    JCP 64 (1986) 490-509           C
C  this version written up       in    CPC 36 (1985) 363-372           C
C                                                                      C
C  WCLBES returns F,G,G',G',SIG for complex ETA, ZZ, and ZLMIN,        C
C   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL inclusive,   C
C   thus giving  complex-energy solutions to the Coulomb Schrodinger   C
C   equation,to the Klein-Gordon equation and to suitable forms of     C
C   the Dirac equation ,also spherical & cylindrical Bessel equations  C
C                                                                      C
C  if ABS(MODE1)                                                       C
C            = 1  get F,G,F',G'   for integer-spaced lambda values     C
C            = 2      F,G      unused arrays must be dimensioned in    C
C            = 3      F,  F'          call to at least length (0:1)    C
C            = 4      F                                                C
C            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )    C
C            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in C
C            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC C
C            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )    C
C                                                                      C
C     if MODE1 < 0 then the values returned are scaled by an exponen-  C
C                  tial factor (dependent only on ZZ) to bring nearer  C
C                  unity the functions for large ABS(ZZ), small ETA ,  C
C                  and ABS(ZL) < ABS(ZZ).                              C
C        Define SCALE = (  0        if MODE1 > 0                       C
C                       (  IMAG(ZZ) if MODE1 < 0  &  KFN < 3           C
C                       (  REAL(ZZ) if MODE1 < 0  &  KFN = 3           C
C        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)                 C
C         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )                   C
C               or EXP(SCALE)       * ( H+, H(1), or K)                C
C               or EXP(-SCALE)      * ( H- or H(2) )                   C
C                                                                      C
C  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G    C
C           =  1   spherical Bessel      "      "     "       j & y    C
C           =  2 cylindrical Bessel      "      "     "       J & Y    C
C           =  3 modified cyl. Bessel    "      "     "       I & K    C
C                                                                      C
C          and where Coulomb phase shifts put in SIG if KFN=0 (not -1) C
C                                                                      C
C  The use of MODE and KFN is independent                              C
C    (except that for KFN=3,  H(1) & H(2) are not given)               C
C                                                                      C
C  With negative orders lambda, WCLBES can still be used but with      C
C    reduced accuracy as CF1 becomes unstable. The user is thus        C
C    strongly advised to use reflection formulae based on              C
C    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi) C
C                                                                      C
C  Precision:  results to within 2-3 decimals of 'machine accuracy',   C
C              except in the following cases:                          C
C              (1) if CF1A fails because X too small or ETA too large  C
C               the F solution  is less accurate if it decreases with  C
C               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)   C
C              (2) if ETA is large (e.g. >> 50) and X inside the       C
C                turning point, then progressively less accuracy       C
C              (3) if ZLMIN is around sqrt(ACCUR) distance from an     C
C               integral order and abs(X) << 1, then errors present.   C
C               RERR traces the main roundoff errors.                  C
C                                                                      C
C   WCLBES is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14  C
C          with a section of doubled REAL*16 for less roundoff errors. C
C          (If no doubled precision available, increase JMAX to eg 100)C
C   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32  C
C   For single precision CDC (48 bits) reassign                        C
C        DOUBLE PRECISION=REAL etc.                                    C
C                                                                      C
C   IPR    on input   = 0 : no printing of error messages              C
C                    ne 0 : print error messages on file 6             C
C   IFAIL  in output = -2 : argument out of range                      C
C                    = -1 : one of the continued fractions failed,     C
C                           or arithmetic check before final recursion C
C                    =  0 : All Calculations satisfactory              C
C                    ge 0 : results available for orders up to & at    C
C                             position NL-IFAIL in the output arrays.  C
C                    = -3 : values at ZLMIN not found as over/underflowC
C                    = -4 : roundoff errors make results meaningless   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     Machine dependent constants :                                    C
C                                                                      C
C     ACCU     target bound on relative error (except near 0 crossings)C
C               (ACCUR should be at least 100 * ACC8)                  C
C     ACC8    smallest number with 1+ACC8 .ne.1 in REAL*8  arithmetic  C
C     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic C
C     FPMAX    magnitude of largest floating point number * ACC8       C
C     FPMIN    magnitude of smallest floating point number / ACC8      C
C                                                                      C
C     Parameters determining region of calculations :                  C
C                                                                      C
C        R20      estimate of (2F0 iterations)/(CF2 iterations)        C
C        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily     C
C        XNEAR    minimum ABS(X) for CF2 to converge accurately        C
C        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series  C
C        JMAX     size of work arrays for Pade accelerations           C
C        NDROP    number of successive decrements to define instabilityC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     C309R1 = CF2,   C309R2 = F11,    C309R3 = F20,
C     C309R4 = CF1R,  C309R5 = CF1C,   C309R6 = CF1A,
C     C309R7 = RCF,   C309R8 = CTIDY.
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C
      COMMON /COULC2/ NFP,N11,NPQ1,NPQ2,N20,KAS1,KAS2
      INTEGER NPQ(2),KAS(2)
      EQUIVALENCE (NPQ(1),NPQ1),(NPQ(2),NPQ2)
      EQUIVALENCE (KAS(1),KAS1),(KAS(2),KAS2)
      DOUBLE PRECISION ZERO,ONE,TWO,HALF
      DOUBLE PRECISION R20,ASYM,XNEAR
      DOUBLE PRECISION FPMAX,FPMIN,FPLMAX,FPLMIN
      DOUBLE PRECISION ACCU,ACC8,ACC16
      DOUBLE PRECISION HPI,TLOG
      DOUBLE PRECISION ERR,RERR,ABSC,ACCUR,ACCT,ACCH,ACCB,C309R4
      DOUBLE PRECISION PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,ABSX
      DOUBLE PRECISION DX1,DETA,DZLL

      LOGICAL LPR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2,NPINT

      PARAMETER(ZERO = 0, ONE = 1, TWO = 2, HALF = ONE/TWO)
      PARAMETER(CI = (0,1), CIH = HALF*CI)
      PARAMETER(R20 = 3, ASYM = 3, XNEAR = HALF)
      PARAMETER(LIMIT = 20000, NDROP = 5, JMAX = 50)

CMSH #include "c309prec.inc"
*
* $Id: c309prec.inc,v 1.1.1.1 1996/04/01 15:01:56 mclareni Exp $
*
* $Log: c309prec.inc,v $
* Revision 1.1.1.1  1996/04/01 15:01:56  mclareni
* Mathlib gen
*
*
*
* c309prec.inc
*

#if defined(CERNLIB_IBM)||defined(CERNLIB_IBMAIX)
      PARAMETER(FPMAX = 1D60, FPMIN = 2D-61)
      PARAMETER(ACCU = 1D-14, ACC8 = 2D-16, ACC16 = 3D-33)

#elif defined(CERNLIB_VAX)
      PARAMETER(FPMAX = 1D21, FPMIN = 3D-22)
      PARAMETER(ACCU = 5D-15, ACC8 = 2D-17, ACC16 = 3D-34)

#elif defined(CERNLIB_CDC)
      PARAMETER(FPMAX = 1E250, FPMIN = 1E-250)
      PARAMETER(ACCU = 1D-13, ACC8 = 4D-15, ACC16 = 2D-29)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_SINGLE))
      PARAMETER(FPMAX = 3D292, FPMIN = 2D-293)
      PARAMETER(ACCU = 1D-14, ACC8 = 2D-16, ACC16 = 2D-34)

#elif defined(CERNLIB_NECSX)
      PARAMETER(FPMAX = 4D74 , FPMIN = 5D-79 )
      PARAMETER(ACCU = 1D-14, ACC8 = 2D-16, ACC16 = 2D-16)

#elif (defined(CERNLIB_SINGLE))&&(!defined(CERNLIB_CDC))&&(!defined(CERNLIB_CONVEX))
      PARAMETER(FPMAX = 1E2451, FPMIN = 1E-2451)
      PARAMETER(ACCU = 1D-12, ACC8 = 7D-15, ACC16 = 2D-29)

#elif defined(CERNLIB_QUAD)
      PARAMETER(FPMAX = 3D292, FPMIN = 2D-292)
      PARAMETER(ACCU = 1D-14, ACC8 = 2D-16, ACC16 = 2D-34)

#elif 1
      PARAMETER(FPMAX = 3D292, FPMIN = 2D-292)
      PARAMETER(ACCU = 1D-14, ACC8 = 2D-16, ACC16 = 2D-16)
#endif


      PARAMETER(HPI  = 1.57079 63267 94896 619D0)
      PARAMETER(TLOG = 0.69314 71805 59945 309D0)

      DIMENSION FC(0:*),GC(0:*),FCP(0:*),GCP(0:*),SIG(0:*)
      DIMENSION XRCF(JMAX,4)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
      NINTC(W)=NINT(DREAL(W))
      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))
      NPINT(W,ACCB)=ABSC(NINTC(W)-W) .LT. ACCB .AND. DREAL(W) .LT. HALF
C
      MODE=MOD(ABS(MODE1),10)
      IFCP=MOD(MODE,2) .EQ. 1
      LPR=IPR .NE. 0
      IFAIL=-2
      N11=0
      NFP=0
      KAS(1)=0
      KAS(2)=0
      NPQ(1)=0
      NPQ(2)=0
      N20=0

      ACCUR=MAX(ACCU,50*ACC8)
      ACCT=HALF*ACCUR
      ACCH=SQRT(ACCUR)
      ACCB=SQRT(ACCH)
      RERR=ACCT
      FPLMAX=LOG(FPMAX)
      FPLMIN=LOG(FPMIN)
C
      CIK=ONE
      IF(KFN .GE. 3) CIK=CI*SIGN(ONE,FPMIN-DIMAG(ZZ))
      X=ZZ*CIK
      ETA=ETA1
      IF(KFN .GT. 0) ETA=ZERO
      ETANE0=ABSC(ETA) .GT. ACC8
      ETAI=ETA*CI
      DELL=ZERO
      IF(KFN .GE. 2) DELL=HALF
      ZM1=ZLMIN-DELL
      SCALE=ZERO
      IF(MODE1 .LT. 0) SCALE=DIMAG(X)
C
      M1=1
      L1=M1+NL
      RLEL=ABS(DIMAG(ETA))+ABS(DIMAG(ZM1)) .LT. ACC8
      ABSX=ABS(X)
      AXIAL=RLEL .AND. ABS(DIMAG(X)) .LT. ACC8*ABSX
      IF(MODE .LE. 2 .AND. ABSX .LT. FPMIN) GO TO 310
      XI=ONE/X
      XLOG=LOG(X)

C       log with cut along the negative real axis, see also OMEGA

      ID=1
      DONEM=.FALSE.
      UNSTAB=.FALSE.
      LF=M1
      IFAIL=-1
   10 ZLM=ZM1+(LF-M1)
      ZLL=ZM1+(L1-M1)
C
C ***       ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels
C
      Z11=ZLL
      IF(ID .LT. 0) Z11=ZLM
      P11=CI*SIGN(ONE,ACC8-DIMAG(ETA))
      LAST=L1
C
C ***       Find phase shifts and Gamow factor at lambda = ZLL
C
      PK=ZLL+ONE
      AA=PK-ETAI
      AB=PK+ETAI
      BB=TWO*PK
      ZLNEG=NPINT(BB,ACCB)
      CLGAA=WLOGAM(AA)
      CLGAB=CLGAA
      IF(ETANE0 .AND. .NOT.RLEL) CLGAB=WLOGAM(AB)
      IF(ETANE0 .AND. RLEL) CLGAB=DCONJG(CLGAA)
      SIGMA=(CLGAA-CLGAB)*CIH
      IF(KFN .EQ. 0) SIG(L1)=SIGMA
      IF(.NOT.ZLNEG) CLL=ZLL*TLOG-HPI*ETA-WLOGAM(BB)+(CLGAA+CLGAB)*HALF
      THETA=X-ETA*(XLOG+TLOG)-ZLL*HPI+SIGMA
C
      TA=(DIMAG(AA)**2+DIMAG(AB)**2+ABS(DREAL(AA))+ABS(DREAL(AB)))*HALF
      IF(ID .GT. 0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GO TO 20
C
C ***         use CF1 instead of CF1A, if predicted to converge faster,
C                 (otherwise using CF1A as it treats negative lambda &
C                  recurrence-unstable cases properly)
C
      RK=SIGN(ONE,DREAL(X)+ACC8)
      P=THETA
      IF(RK .LT. 0) P=-X+ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA
      XRCF(1,1)=PK
      F=RK*C309R6(X*RK,ETA*RK,ZLL,P,ACCT,JMAX,NFP,FEST,ERR,FPMAX,XRCF,
     1                                      XRCF(1,3),XRCF(1,4))
      FESL=LOG(FEST)+ABS(DIMAG(X))
      NFP=-NFP
      IF(NFP .LT. 0 .OR. UNSTAB .AND. ERR .LT. ACCB) GO TO 40
      IF(.NOT.ZLNEG .OR. UNSTAB .AND. ERR .GT. ACCB) GO TO 20
      IF(LPR) WRITE(6,1060) '-L',ERR
      IF(ERR.GT.ACCB) GO TO 280
      GO TO 40
C
C ***    evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)
C
   20 IF(AXIAL) THEN
       DX1=X
       DETA=ETA
       DZLL=ZLL
       F=C309R4(DX1,DETA,DZLL,ACC8,SF,RK,ETANE0,LIMIT,ERR,NFP,
     1          FPMIN,FPMAX,LPR)
       FCL=SF
       TPK1=RK
      ELSE
       F=C309R5(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     1          FPMIN,FPMAX,LPR)
      END IF
      IF(ERR .GT. ONE) THEN
       IFAIL=-1
       GO TO 290
      END IF
C
C ***  Make a simple check for CF1 being badly unstable:
C
      IF(ID .GE. 0) THEN
       UNSTAB=DREAL((ONE-ETA*XI)*CI*DIMAG(THETA)/F) .GT. ZERO
     1  .AND. .NOT.AXIAL .AND. ABS(DIMAG(THETA)) .GT. -LOG(ACC8)*HALF
     2  .AND. ABSC(ETA)+ABSC(ZLL) .LT. ABSC(X)
       IF(UNSTAB) THEN
        ID=-1
        LF=L1
        L1=M1
        RERR=ACCT
        GO TO 10
       END IF
      END IF
C
C *** compare accumulated phase FCL with asymptotic phase for G(k+1) :
C     to determine estimate of F(ZLL) (with correct sign) to start recur
C
      W=X*X*(HALF/TPK1+ONE/TPK1**2)+ETA*(ETA-TWO*X)/TPK1
      FESL=(ZLL+ONE)*XLOG+CLL-W-LOG(FCL)
   40 FESL=FESL-ABS(SCALE)
      RK=MAX(DREAL(FESL),FPLMIN*HALF)
      FESL=DCMPLX(MIN(RK,FPLMAX*HALF),DIMAG(FESL))
      FEST=EXP(FESL)
C
      RERR=MAX(RERR,ERR,ACC8*ABS(DREAL(THETA)))
C
      FCL=FEST
      FPL=FCL*F
      IF(IFCP) FCP(L1)=FPL
      FC(L1)=FCL
C
C *** downward recurrence to lambda = ZLM. array GC,if present,stores RL
C
      I=MAX(-ID,0)
      ZL=ZLL+I
      MONO=0
      OFF=ABS(FCL)
      TA=ABSC(SIGMA)

C
C     CORRESPONDS TO   DO 70 L = L1-ID,LF,-ID
C
      L=L1-ID
      LC70=(L1-LF)/ID
   70 IF(LC70 .LE. 0) GO TO 80
      IF(ETANE0) THEN
       IF(RLEL) THEN
        DSIG=ATAN2(DREAL(ETA),DREAL(ZL))
        RL=SQRT(DREAL(ZL)**2+DREAL(ETA)**2)
       ELSE
        AA=ZL-ETAI
        BB=ZL+ETAI
        IF(ABSC(AA) .LT. ACCH .OR. ABSC(BB) .LT. ACCH) GOTO 50
        DSIG=(LOG(AA)-LOG(BB))*CIH
        RL=AA*EXP(CI*DSIG)
       END IF
       IF(ABSC(SIGMA) .LT. TA*HALF) THEN

C               re-calculate SIGMA because of accumulating roundoffs:

        SL=(WLOGAM(ZL+I-ETAI)-WLOGAM(ZL+I+ETAI))*CIH
        RL=(ZL-ETAI)*EXP(CI*ID*(SIGMA-SL))
        SIGMA=SL
        TA=ZERO
       ELSE
        SIGMA=SIGMA-DSIG*ID
       END IF
       TA=MAX(TA,ABSC(SIGMA))
       SL=ETA+ZL*ZL*XI
       PL=ZERO
       IF(ABSC(ZL) .GT. ACCH) PL=(SL*SL-RL*RL)/ZL
       FCL1=(FCL*SL+ID*ZL*FPL)/RL
       SF=ABS(FCL1)
       IF(SF .GT. FPMAX) GO TO 350
       FPL=(FPL*SL+ID*PL*FCL)/RL
       IF(MODE .LE. 1) GCP(L+ID)=PL*ID
      ELSE

C                      ETA = 0, including Bessels.  NB RL==SL

       RL=ZL*XI
       FCL1=FCL*RL+FPL*ID
       SF=ABS(FCL1)
       IF(SF .GT. FPMAX) GO TO 350
       FPL=(FCL1*RL-FCL)*ID
      END IF
      MONO=MONO+1
      IF(SF .GE. OFF) MONO=0
      FCL=FCL1
      OFF=SF
      FC(L)=FCL
      IF(IFCP) FCP(L)=FPL
      IF(KFN .EQ. 0) SIG(L)=SIGMA
      IF(MODE .LE. 2) GC(L+ID)=RL
      ZL=ZL-ID
      IF(MONO .LT. NDROP .OR. AXIAL .OR.
     1             DREAL(ZLM)*ID .GT. -NDROP .AND. .NOT.ETANE0) THEN
       L=L-ID
       LC70=LC70-1
       GO TO 70
      END IF
      UNSTAB=.TRUE.
C
C ***    take action if cannot or should not recur below this ZL:

   50 ZLM=ZL
      LF=L
      IF(ID .LT. 0) GO TO 380
      IF(.NOT.UNSTAB) LF=L+1
      IF(L+MONO .LT. L1-2 .OR. ID .LT. 0 .OR. .NOT.UNSTAB) GO TO 80

C             otherwise, all L values (for stability) should be done
C                        in the reverse direction:

      ID=-1
      LF=L1
      L1=M1
      RERR=ACCT
      GO TO 10

   80 IF(FCL .EQ. ZERO) FCL=ACC8
      F=FPL/FCL
C
C *** Check, if second time around, that the 'f' values agree!
C
      IF(ID .GT. 0) FIRST=F
      IF(DONEM) RERR=MAX(RERR,ABSC(F-FIRST)/ABSC(F))
      IF(DONEM) GO TO 90
C
      NOCF2=.FALSE.
      THETAM=X-ETA*(XLOG+TLOG)-ZLM*HPI+SIGMA
C
C *** on left x-plane, determine OMEGA by requiring cut on -x axis
C     on right x-plane, choose OMEGA (using estimate based on THETAM)
C       so H(omega) is smaller and recurs upwards accurately.
C     (x-plane boundary is shifted to give CF2(LH) a chance to converge)
C
      OMEGA=SIGN(ONE,DIMAG(X)+ACC8)
      IF(DREAL(X) .GE. XNEAR) OMEGA=SIGN(ONE,DIMAG(THETAM)+ACC8)
      SFSH=EXP(OMEGA*SCALE-ABS(SCALE))
      OFF=EXP(MIN(TWO*MAX(ABS(DIMAG(X)),ABS(DIMAG(THETAM)),
     1                         ABS(DIMAG(ZLM))*3),FPLMAX))
      EPS=MAX(ACC8,ACCT*HALF/OFF)
C
C ***    Try first estimated omega, then its opposite,
C        to find the H(omega) linearly independent of F
C        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega)
C
   90 DO 100 L = 1,2
      LH=1
      IF(OMEGA .LT. ZERO) LH=2
      PM=CI*OMEGA
      ETAP=ETA*PM
      IF(DONEM) GO TO 130
      PQ1=ZERO
      PACCQ=ONE
      KASE=0
C
C ***            Check for small X, i.e. whether to avoid CF2 :
C
      IF(MODE .GE. 3 .AND. ABSX .LT. ONE ) GO TO 190
      IF(MODE .LT. 3 .AND. (NOCF2 .OR. ABSX .LT. XNEAR .AND.
     1   ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM) .LT. 4)) THEN
       KASE=5
       GO TO 120
      END IF
C
C ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM
C
      PQ1=C309R1(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH,
     1             LPR,ACCUR,DELL)
      ERR=ERR*MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8))
C
C *** check if impossible to get F-PQ accurately because of cancellation
C                original guess for OMEGA (based on THETAM) was wrong
C                Use KASE 5 or 6 if necessary if Re(X) < XNEAR
      IF(ERR .LT. ACCH) GO TO 110
      NOCF2=DREAL(X) .LT. XNEAR .AND. ABS(DIMAG(X)) .LT. -LOG(ACC8)
  100 OMEGA=-OMEGA
      IF(UNSTAB) GO TO 360
      IF(LPR .AND. DREAL(X) .LT. -XNEAR) WRITE(6,1060) '-X',ERR
  110 RERR=MAX(RERR,ERR)
C
C ***  establish case of calculation required for irregular solution
C
  120 IF(KASE .GE. 5) GO TO 130
      IF(DREAL(X) .GT. XNEAR) THEN

C          estimate errors if KASE 2 or 3 were to be used:

       PACCQ=EPS*OFF*ABSC(PQ1)/MAX(ABS(DIMAG(PQ1)),ACC8)
      END IF
      IF(PACCQ .LT. ACCUR) THEN
       KASE=2
       IF(AXIAL) KASE=3
      ELSE
       KASE=1
       IF(NPQ(1)*R20 .LT. JMAX) KASE=4
C             i.e. change to kase=4 if the 2F0 predicted to converge
      END IF
  130 GO TO (190,140,150,170,190,190), ABS(KASE)
  140 IF(.NOT.DONEM)
C
C ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2)
C
     1  PQ2=C309R1(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH,
     2             LPR,ACCUR,DELL)
      P=(PQ2+PQ1)*HALF
      Q=(PQ2-PQ1)*HALF*PM
      GO TO 160
  150 P=DREAL(PQ1)
      Q=DIMAG(PQ1)
C
C ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*
C
      PQ2=DCONJG(PQ1)
C
C *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL
C
  160 W=(PQ1-F)*(PQ2-F)
      SF=EXP(-ABS(SCALE))
      FCM=SQRT(Q/W)*SF

C                  any SQRT given here is corrected by
C                  using sign for FCM nearest to phase of FCL

      IF(DREAL(FCM/FCL) .LT. ZERO) FCM=-FCM
      GAM=(F-P)/Q
      TA=ABSC(GAM+PM)
      PACCQ=EPS*MAX(TA,ONE/TA)
      HCL=FCM*(GAM+PM)*(SFSH/(SF*SF))

C                            Consider a KASE = 1 Calculation

      IF(PACCQ .GT. ACCUR .AND. KASE .GT. 0) THEN
       F11V=C309R2(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
       IF(ERR .LT. PACCQ) GO TO 200
      END IF
      RERR=MAX(RERR,PACCQ)
      GO TO 230
C
C *** Arrive here if KASE = 4
C     to evaluate the exponentially decreasing H(LH) directly.
C
  170 IF(DONEM) GO TO 180
      AA=ETAP-ZLM
      BB=ETAP+ZLM+ONE
      F20V=C309R3(AA,BB,-HALF*PM*XI,ACCT,JMAX,ERR,FPMAX,N20,XRCF)
      IF(N20 .LE. 0) GO TO 190
      RERR=MAX(RERR,ERR)
      HCL=FPMIN
      IF(ABS(DREAL(PM*THETAM)+OMEGA*SCALE) .GT. FPLMAX) GO TO 330
  180 HCL=F20V*EXP(PM*THETAM+OMEGA*SCALE)
      FCM=SFSH/((F-PQ1)*HCL)
      GO TO 230
C
C *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)
C
C           for small values of X, calculate F(X,SL) directly from 1F1
C               using DREAL*16 arithmetic if possible.
C           where Z11 = ZLL if ID>0, or = ZLM if ID<0
C
  190 F11V=C309R2(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
  200 IF(N11 .LT. 0) THEN

C                               F11 failed from BB = negative integer

       IF(LPR) WRITE(6,1060) '-L',ONE
       IFAIL=-1
       GO TO 290
      END IF

C                      Consider a KASE 2 or 3 calculation :

      IF(ERR .GT. PACCQ .AND. PACCQ .LT. ACCB) THEN
       KASE=-2
       IF(AXIAL) KASE=-3
       GO TO 130
      END IF
      RERR=MAX(RERR,ERR)
      IF(ERR .GT. FPMAX) GO TO 370
      IF(ID .LT. 0) CLL=Z11*TLOG-HPI*ETA-WLOGAM(BB)
     1                        +WLOGAM(Z11+ONE+P11*ETA)-P11*SIGMA
      EK=(Z11+ONE)*XLOG-P11*X+CLL-ABS(SCALE)
      IF(ID .GT. 0) EK=EK-FESL+LOG(FCL)
      IF(DREAL(EK) .GT. FPLMAX) GO TO 350
      IF(DREAL(EK) .LT. FPLMIN) GO TO 340
      FCM=F11V*EXP(EK)
      IF(KASE .GE. 5) THEN
       IF(ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)) .LT. ACCH) KASE=6
C
C ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so
C ***      use an expansion for irregular soln from origin :
C
       SL=ZLM
       ZLNEG=DREAL(ZLM) .LT. -ONE+ACCB
       IF(KASE .EQ. 5 .OR. ZLNEG) SL=-ZLM-ONE
       PK=SL+ONE
       AA=PK-ETAP
       AB=PK+ETAP
       BB=TWO*PK
       CLGAA=WLOGAM(AA)
       CLGAB=CLGAA
       IF(ETANE0) CLGAB=WLOGAM(AB)
       CLGBB=WLOGAM(BB)
       IF(KASE .EQ. 6 .AND. .NOT.ZLNEG) THEN
        IF(NPINT(AA,ACCUR)) CLGAA=CLGAB-TWO*PM*SIGMA
        IF(NPINT(AB,ACCUR)) CLGAB=CLGAA+TWO*PM*SIGMA
       END IF
       CLL=SL*TLOG-HPI*ETA-CLGBB+(CLGAA+CLGAB)*HALF
       DSIG=(CLGAA-CLGAB)*PM*HALF
       IF(KASE .EQ. 6) P11=-PM
       EK=PK*XLOG-P11*X+CLL-ABS(SCALE)
       SF=EXP(-ABS(SCALE))
       CHI=ZERO
       IF(.NOT.(KASE .EQ. 5 .OR. ZLNEG)) GO TO 210
C
C ***  Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)
C
C      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2
C
       CHI=SIGMA-DSIG-(ZLM-SL)*HPI
       F11V=C309R2(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),
     1             FPMAX,ACC8,ACC16)
       RERR=MAX(RERR,ERR)
       IF(KASE .EQ. 6) GO TO 210
       FESL=F11V*EXP(EK)
       FCL1=EXP(PM*CHI)*FCM
       HCL=FCL1-FESL
       RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL))
       HCL=HCL/SIN(CHI)*(SFSH/(SF*SF))
       GO TO 220
C
C *** Use the logarithmic expansion for the irregular solution (KASE 6)
C        for the case that BB is integral so sin(CHI) would be zero.
C
  210  RL=BB-ONE
       N=NINTC(RL)
       ZLOG=XLOG+TLOG-PM*HPI
       CHI=CHI+PM*THETAM+OMEGA*SCALE+AB*ZLOG
       AA=ONE-AA
       IF(NPINT(AA,ACCUR)) THEN
        HCL=ZERO
       ELSE
        IF(ID .GT. 0 .AND. .NOT.ZLNEG) F11V=FCM*EXP(-EK)
        HCL=EXP(CHI-CLGBB-WLOGAM(AA))*(-1)**(N+1)*(F11V*ZLOG+
     1   C309R2(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16))
        RERR=MAX(RERR,ERR)
       END IF
       IF(N .GT. 0) THEN
        EK=CHI+WLOGAM(RL)-CLGAB-RL*ZLOG
        DF=C309R2(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16)
        HCL=HCL+EXP(EK)*DF
       END IF
       RERR=MAX(RERR,TWO*ABS(BB-NINTC(BB)))
  220  PQ1=F-SFSH/(FCM*HCL)
      ELSE
       IF(MODE .LE. 2) HCL=SFSH/((F-PQ1)*FCM)
       KASE=1
      END IF
C
C ***  Now have absolute normalisations for Coulomb Functions
C          FCM & HCL  at lambda = ZLM
C      so determine linear transformations for Functions required :
C
  230 IH=ABS(MODE1)/10
      IF(KFN .EQ. 3) IH=(3-DIMAG(CIK))/2+HALF
      P11=ONE
      IF(IH .EQ. 1) P11=CI
      IF(IH .EQ. 2) P11=-CI
      DF=-PM
      IF(IH .GE. 1) DF=-PM+P11
      IF(ABSC(DF) .LT. ACCH) DF=ZERO
C
C *** Normalisations for spherical or cylindrical Bessel functions
C
      IF(KFN .LE. 0) THEN
       ALFA=ZERO
       BETA=ONE
      ELSE IF(KFN .EQ. 1) THEN
       ALFA=XI
       BETA=XI
      ELSE
       ALFA=XI*HALF
       BETA=SQRT(XI/HPI)
       IF(DREAL(BETA) .LT. ZERO) BETA=-BETA
      END IF
      AA=ONE
      IF(KFN .GT. 0) AA=-P11*BETA

C                Calculate rescaling factors for I & K output

      IF(KFN .GE. 3) THEN
       P=EXP((ZLM+DELL)*HPI*CIK)
       AA=BETA*HPI*P
       BETA=BETA/P
       Q=CIK*ID
      END IF

C                  Calculate rescaling factors for GC output

      IF(IH .EQ. 0) THEN
       TA=ABS(SCALE)+DIMAG(PM)*SCALE
       RK=ZERO
       IF(TA .LT. FPLMAX) RK=EXP(-TA)
      ELSE
       TA=ABS(SCALE)+DIMAG(P11)*SCALE
       IF(ABSC(DF) .GT. ACCH .AND. TA .GT. FPLMAX) GO TO 320
       IF(ABSC(DF) .GT. ACCH) DF=DF*EXP(TA)
       SF=TWO*(LH-IH)*SCALE
       RK=ZERO
       IF(SF .GT. FPLMAX) GO TO 320
       IF(SF .GT. FPLMIN) RK=EXP(SF)
      END IF
      KAS((3-ID)/2)=KASE
      W=FCM/FCL
      IF(LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GO TO 340
      IF(MODE .GE. 3) GO TO 240
      IF(LPR .AND. ABSC(F-PQ1) .LT. ACCH*ABSC(F))
     1                             WRITE(6,1020) LH,ZLM+DELL
      HPL=HCL*PQ1
      IF(ABSC(HPL) .LT. FPMIN .OR. ABSC(HCL) .LT. FPMIN) GO TO 330
C
C *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)
C *** renormalise FC,FCP at each lambda
C ***    ZL   = ZLM - MIN(ID,0) here
C
  240 DO 270 L = LF,L1,ID
      FCL=W*FC(L)
      IF(ABSC(FCL) .LT. FPMIN) GO TO 340
      IF(IFCP) FPL=W*FCP(L)
      FC(L)=BETA*FCL
      IF(IFCP) FCP(L)=BETA*(FPL-ALFA*FCL)*CIK
      FC(L)=C309R8(FC(L),ACCUR)
      IF(IFCP) FCP(L)=C309R8(FCP(L),ACCUR)
      IF(MODE .GE. 3) GO TO 260
      IF(L .EQ. LF) GO TO 250
      ZL=ZL+ID
      ZID=ZL*ID
      RL=GC(L)
      IF(ETANE0) THEN
       SL=ETA+ZL*ZL*XI
       IF(MODE .EQ. 1) THEN
        PL=GCP(L)
       ELSE
        PL=ZERO
        IF(ABSC(ZL) .GT. ACCH) PL=(SL*SL-RL*RL)/ZID
       END IF
       HCL1=(SL*HCL-ZID*HPL)/RL
       HPL=(SL*HPL-PL*HCL)/RL
      ELSE
       HCL1=RL*HCL-HPL*ID
       HPL=(HCL-RL*HCL1)*ID
      END IF
      HCL=HCL1
      IF(ABSC(HCL) .GT. FPMAX) GO TO 320
  250 GC(L)=AA*(RK*HCL+DF*FCL)
      IF(MODE .EQ. 1) GCP(L)=(AA*(RK*HPL+DF*FPL)-ALFA*GC(L))*CIK
      GC(L)=C309R8(GC(L),ACCUR)
      IF(MODE .EQ. 1) GCP(L)=C309R8(GCP(L),ACCUR)
      IF(KFN .GE. 3) AA=AA*Q
  260 IF(KFN .GE. 3) BETA=-BETA*Q
  270 LAST=MIN(LAST,(L1-L)*ID)
      GO TO 280
C
C *** Come here after all soft errors to determine how many L values ok
C
  310 IF(LPR) WRITE(6,1000) ZZ
      GO TO 999
  320 IF(LPR) WRITE(6,1010) ZL+DELL,'IR',HCL,'>',FPMAX
      GO TO 280
  330 IF(LPR) WRITE(6,1010) ZL+DELL,'IR',HCL,'<',FPMIN
      GO TO 280
  340 IF(LPR) WRITE(6,1010) ZL+DELL,'  ',FCL,'<',FPMIN
      GO TO 280
  350 IF(LPR) WRITE(6,1010) ZL+DELL,'  ',FCL,'>',FPMAX
      GO TO 280
  360 IF(LPR) WRITE(6,1030) ZL+DELL
      GO TO 280
  370 IF(LPR) WRITE(6,1040) Z11,I
      IFAIL=-1
      GO TO 290
  380 IF(LPR) WRITE(6,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL
      IFAIL=-1
      GO TO 290
  280 IF(ID .GT. 0 .OR. LAST .EQ. 0) IFAIL=LAST
      IF(ID .LT. 0 .AND. LAST .NE. 0) IFAIL=-3
C
C *** Come here after ALL errors for this L range (ZLM,ZLL)
C
C *** so on first block, 'F' started decreasing monotonically,
C                        or hit bound states for low ZL.
C     thus redo M1 to LF-1 in reverse direction, i.e. do
C      CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)
C
  290 IF(ID .GT. 0 .AND. LF .NE. M1) THEN
       ID=-1
       IF(.NOT.UNSTAB) LF=LF-1
       DONEM=UNSTAB
       LF=MIN(LF,L1)
       L1=M1
       GO TO 10
      END IF
      IF(IFAIL .LT. 0) GO TO 999
      IF(LPR .AND. RERR .GT. ACCB) WRITE(6,1070) RERR
      IF(RERR .GT. 0.1D0) IFAIL=-4
  999 DO 998 L = 1,NL+1
      FC(L-1)=FC(L)
      GC(L-1)=GC(L)
      FCP(L-1)=FCP(L)
      GCP(L-1)=GCP(L)
  998 SIG(L-1)=SIG(L)
      RETURN
C
 1000 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',
     2 1P,2D10.2,' ABS(X) TOO SMALL')
 1010 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P,2E10.1,') ',
     2 A1,E10.1)
 1020 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''H(',I1,
     2 ')'' IS LOST AT ZL = ',2F7.2/1X,'*****',22X,'(EG. COULOMB ',
     3 'EIGENSTATE OR CF1 UNSTABLE)')
 1030 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     2 '(ETA & L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE AT L = ',2F8.2)
 1040 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'OVERFLOW IN 1F1 SERIES AT ZL = ',2F8.3,' AT TERM',I5)
 1050 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR OR MULTIPLE',
     2 ' INSTABILITIES PRESENT'/
     3   1X,'*****',22X,'TRY CALLING TWICE, 1ST FOR ZL FROM',2F8.3,
     4 ' TO',2F8.3,' (INCL)'/1X,'*****',41X,'2ND FOR ZL FROM',2F8.3,
     5 ' TO',2F8.3)
 1060 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'WARNING: AS ''',A2,''' REFLECTION RULES NOT USED ',
     2 'ERRORS CAN BE UP TO',1PD12.2)
 1070 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     1 'WARNING: OVERALL ROUNDOFF ERROR APPROXIMATELY',1PE11.1)
      END
#endif
*CMZ :          02/05/2017  15.09.41  by  Michael Scheer
*-- Author :
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
