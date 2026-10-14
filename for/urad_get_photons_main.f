*CMZ :          14/10/2026  12.20.23  by  Michael Scheer
*CMZ :  4.03/00 07/10/2026  07.28.50  by  Michael Scheer
*CMZ :  4.02/01 02/09/2026  08.56.08  by  Michael Scheer
*CMZ :  4.02/00 16/09/2025  09.13.49  by  Michael Scheer
*CMZ :  4.01/07 18/10/2024  08.57.02  by  Michael Scheer
*CMZ :  4.01/05 16/04/2024  14.41.20  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.32.18  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  10.57.05  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  13.32.32  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.57.49  by  Michael Scheer
*-- Author : Michael Scheer
*KEEP,gplhint.
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
*KEND.

      include 'urad_get_photons.f'

      program urad_get_photons_main

      implicit none

      real*8 ::

     &     Ebeam=1.722        !Beam energy [GeV]
     &  ,  Curr=0.1           !Current [A]
     &  ,  Step=0.2           !Tracking step size [mm]

     &  ,  PinX=10000.0       !X of Pinhole [mm]
     &  ,  PinY=0.0           !Y of Pinhole [mm]
     &  ,  PinZ=0.0           !Z of Pinhole [mm]
     &  ,  PinW=2.0           !Width of Pinhole [mm]
     &  ,  PinH=1.0           !Height of Pinhole [mm]

     &  ,  Perlen=50.0        !Period-length [mm]
     &  ,  Shift=0.0          !Shift [mm] (for spectra only)
     &  ,  Beffv=0.651938438  !B0 vert. [T]
     &  ,  Beffh=0.0          !B0 hori. [T]
     &  ,  Harm=100.          !Harmonic to overwrite Beffv and Beffh, if Nharm is set

     &  ,  EphMin=695.0       !Min. Photon Energy [eV]
     &  ,  EphMax=705.0       !Max. Photon Energy [eV]

     &  ,  Espread=0.001      !Rel. energy spread

     &  ,  xBeta=0.0          !X-reference for optical functions
     &  ,  BetaH=14.0         !Hori. Beta function
     &  ,  AlphaH=0.0         !Hori. Alpha, i.e. -dBetaH/ds/2.0
     &  ,  BetaV=3.4          !Vert. Beta function
     &  ,  AlphaV=0.0         !Vert. Alpha, i.e. -dBetaV/ds/2.0
     &  ,  EmitH=4.4          !Hor. Emit. [nm-rad]
     &  ,  EmitV=0.066        !Ver. Emit. [nm-rad]
     &  ,  Disph=0.0          !Hori. Dispersion [mm]
     &  ,  Dispph=0.0         !Derivative of Hori. Dispersion [mrad]
     &  ,  Dispv=0.0          !Vert. Dispersion [mm]
     &  ,  Disppv=0.0         !Derivative of Vert. Dispersion [mrad]
     &  ,  Pherror=0.0        !Phase errors
     &  ,  PhGshift=9999.0    !Global phase shift. 9999.: Hori. amplitude is zero for (PinX,PinY,PinZ))

      integer ::

     &  Nelec=1000000         !Nelec
     &  ,  Noranone=1         !No random change for first e-

     &  ,  Nharm=1            !Number of harmonic Harm referes to, overwrites Beffv and Beffh

     &  ,  NpinY=51           !Number of vert. points
     &  ,  NpinZ=51           !Number of hori. points
     &  ,  Modepin=1          !Monte-Carlo mode [0,1]

     &  ,  Nper=100           !Number of periods

     &  ,  Nepho=1            !Number of Photon Energies
     &  ,  Modeph=0           !Mode for phase-error
     &  ,  ModeRan=1          !All grid points (zero) or Nelec random grid points (not zero)

      real, dimension (:), allocatable ::  photons

      !photons(1) = Ephoton
      !photons(2) = y
      !photons(3) = z
      !photons(4) = yp
      !photons(5) = zp
      !photons(6) = S0
      !photons(7) = S1
      !photons(8) = S2
      !photons(9) = S3
      !photons(10) = Ephoton
      !...

      integer :: modewave=-1, ndimpho,luno,k, kwrite=1

      real*8 b0eff,park,wlen1,rhv,perl,harm1,gamma

*KEEP,PHYCONparam,T=F77.
c-----------------------------------------------------------------------
c     phyconparam.cmn
c-----------------------------------------------------------------------

      complex*16, parameter :: zone1=(1.0d0,0.0d0), zi1=(0.0d0,1.0d0)

      complex*16, dimension(4,3), parameter ::
     &  vstokes=reshape([
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0, -0.70710678118654746d0),
     &  ( 0.0000000000000000d0, -0.70710678118654746d0),
     &  ( 0.70710678118654746d0, 0.0000000000000000d0),
     &  (-0.70710678118654746d0,-0.70710678118654746d0),
     &  ( 0.70710678118654746d0, 0.0000000000000000d0),
     &  (-0.70710678118654746d0, 0.0000000000000000d0),
     &  (-0.70710678118654746d0, 0.0000000000000000d0)
     &  ],[4,3])

c      vstokes(1,1)=( 0.0d0,        0.0d0)      !horizontal polarization
c      vstokes(1,2)=( 0.0d0,        0.0d0)
c      vstokes(1,3)=(-sqrt(1./2.),       -sqrt(1./2.))
c
c      vstokes(2,1)=( 0.0d0,        0.0d0)      !right handed polarization
c      vstokes(2,2)=( 0.0d0,       -sqrt(1./2.))
c      vstokes(2,3)=(+sqrt(1./2.),        0.0d0)
c
c      vstokes(3,1)=( 0.0d0,        0.0d0)      !left handed polarization
c      vstokes(3,2)=( 0.0d0,       -sqrt(1./2.))
c      vstokes(3,3)=(-sqrt(1./2.),        0.0d0)
c
c      vstokes(4,1)=( 0.0d0,        0.0d0)      !45 degree linear polarization
c      vstokes(4,2)=( sqrt(1./2.),        0.0d0)
c      vstokes(4,3)=(-sqrt(1./2.),        0.0d0)

      double precision, parameter ::
     &  HBAREV1=6.58211889D-16
     &  ,CLIGHT1=2.99792458D8
     &  ,EMASSKG1=9.10938188D-31
     &  ,EMASSE1=0.510998902D6
     &  ,EMASSG1=0.510998902D-3
     &  ,ECHARGE1=1.602176462D-19
     &  ,ERAD1=2.8179380D-15
     &  ,EPS01=8.854187817D-12
     &  ,PI1=3.141592653589793D0
     &  ,rmu04pi1=1.0D-7
     &  ,dnull1=0.0d0
     &  ,done1=1.0d0
     & ,HPLANCK1=6.62606876D-34

      double precision, parameter ::
     & GRARAD1=PI1/180.0d0
     & ,RADGRA1=180.0d0/PI1
     & ,HBAR1=HBAREV1*ECHARGE1
     & ,WTOE1=CLIGHT1*HPLANCK1/ECHARGE1*1.0d9
     & ,CQ1=55.0d0/32.0d0/DSQRT(3.0D0)*HBAR1/EMASSKG1/CLIGHT1
     & ,CGAM1=4.0d0/3.0d0*PI1*ERAD1/EMASSG1**3
     & ,POL1CON1=8.0d0/5.0d0/DSQRT(3.0D0)
     & ,POL2CON1=8.0d0/5.0d0/DSQRT(3.0D0)/2.0d0/PI1/3600.0d0
     &  *EMASSKG1/HBAR1/ERAD1*EMASSG1**5
     & ,TWOPI1=2.0D0*PI1
     & ,HALFPI1=PI1/2.0D0
     & ,sqrttwopi1=sqrt(twopi1)
     & ,rmu01=4.0D0*PI1/1.0D7
     & ,alpha1=echarge1**2/(4.0d0*pi1*eps01*hbar1*clight1)
     & ,gaussn1=1.0d0/sqrt(twopi1)
     & ,cK934=ECHARGE1/(2.0d0*PI1*EMASSKG1*CLIGHT1)/100.0d0
     & ,powcon1=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1
     &  ,gamma1=1.0d0/emassg1
     &  ,emom1=emasse1*dsqrt((gamma1-1.0d0)*(gamma1+1.0d0))
     &  ,rho1=emom1/clight1
     &  ,omegac1=1.5d0*gamma1**3*clight1/rho1
     &  ,ecdipev1=omegac1*hbar1/echarge1
     &  ,ecdipkev1=ecdipev1/1000.0d0

c-----------------------------------------------------------------------
c     end of phyconparam.cmn
c-----------------------------------------------------------------------
*KEND.

      if (nharm.gt.0.and.harm.gt.0.0d0) then

        gamma=ebeam/emassg1

        harm1=harm/nharm
        wlen1=wtoe1/abs(harm/nharm)

        perl=perlen/1000.0d0
        park=2.0d0*(wlen1/(perl*1.0d9/2.0d0/gamma**2)-1.0d0)

        if (park.lt.0.0d0) then
          write(6,*)
     &      '*** Error in urad_get_photons_main:'
          write(6,*)
     &      'Inconsistent values of nharm, harm, and perlen'
          write(6,*)' '
          stop
        endif

        park=sqrt(park)
        b0eff=park/(echarge1*perl/(2.*pi1*emasskg1*clight1))

        if (beffh.eq.0.0d0.and.beffv.ne.0d0) then
          beffv=beffv/abs(beffv)*b0eff
        else if (beffv.eq.0.0d0.and.beffh.ne.0d0) then
          beffh=beffh/abs(beffh)*b0eff
        else
          rhv=beffh/beffv
          beffh=b0eff/sqrt(1.0d0+1.0d0/rhv**2)*beffh/abs(beffh)
          beffv=beffh/rhv
        endif
      endif

      if (modepin.eq.0) then
        ndimpho=9*nepho*nelec*npinz*npiny
        allocate(photons(ndimpho))
      else
        ndimpho=9*nepho*nelec
        allocate(photons(9*nepho*nelec))
      endif

      call urad_get_photons(
     &      nelec,noranone,
     &      perlen,shift,nper,beffv,beffh,
     &      ebeam,curr,step,
     &      pinx,piny,pinz,pinw,pinh,npiny,npinz,modepin,
     &      nepho,ephmin,ephmax,
     &      xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &      disph,dispph,dispv,disppv,
     &      modeph,pherror,phgshift,modewave,ndimpho,photons)

      call util_break

      if (kwrite.ne.0) then

        open(newunit=luno,file='get_photons.pho')

        do k=1,ndimpho,9
          write(luno,*) photons(k),pinx*1000.,photons(k+1:k+8)
        enddo

        close(luno)

      endif

      end

