*CMZ :          13/10/2026  14.13.23  by  Michael Scheer
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

     &     Ebeam=1.722          !Beam energy [GeV]
     &  ,  Curr=0.1          !Current [A]
     &  ,  Step=0.2          !Tracking step size [mm]

     &  ,  PinX=10000.0          !X of Pinhole [mm]
     &  ,  PinY=0.0          !Y of Pinhole [mm]
     &  ,  PinZ=0.0          !Z of Pinhole [mm]
     &  ,  PinW=2.0          !Width of Pinhole [mm]
     &  ,  PinH=1.0          !Height of Pinhole [mm]

     &  ,  Perlen=50.0          !Period-length [mm]
     &  ,  Shift=0.0          !Shift [mm] (for spectra only)
     &  ,  Beffv=0.651938438  !B0 vert. [T]
     &  ,  Beffh=0.0          !B0 hori. [T]

     &  ,  EphMin=695.0          !Min. Photon Energy [eV]
     &  ,  EphMax=705.0          !Max. Photon Energy [eV]

     &  ,  Espread=0.001          !Rel. energy spread

     &  ,  xBeta=0.0           !X-reference for optical functions
     &  ,  BetaH=14.0          !Hori. Beta function
     &  ,  AlphaH=0.0           !Hori. Alpha, i.e. -dBetaH/ds/2.0
     &  ,  BetaV=3.4          !Vert. Beta function
     &  ,  AlphaV=0.0           !Vert. Alpha, i.e. -dBetaV/ds/2.0
     &  ,  EmitH=4.4          !Hor. Emit. [nm-rad]
     &  ,  EmitV=0.066          !Ver. Emit. [nm-rad]
     &  ,  Disph=0.0          !Hori. Dispersion [mm]
     &  ,  Dispph=0.0          !Derivative of Hori. Dispersion [mrad]
     &  ,  Dispv=0.0          !Vert. Dispersion [mm]
     &  ,  Disppv=0.0          !Derivative of Vert. Dispersion [mrad]
     &  ,  Pherror=0.0          !Phase errors
     &  ,  PhGshift=9999.0          !Global phase shift. 9999.: Hori. amplitude is zero for (PinX,PinY,PinZ))

      integer ::

     &  Nelec=1000000          !Nelec
     &,  Noranone=1          !No random change for first e-

     &,  NpinY=51          !Number of vert. points
     &,  NpinZ=51          !Number of hori. points
     &,  Modepin=1          !Monte-Carlo mode [0,1]

     &,  Nper=100          !Number of periods

     &,  Nepho=1          !Number of Photon Energies
     &,  Modeph=0          !Mode for phase-error
     &,  ModeRan=1          !All grid points (zero) or Nelec random grid points (not zero)

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

      integer :: modewave=-1, ndimpho,luno,k, kwrite=0

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

