*CMZ :          03/10/2025  17.10.54  by  Michael Scheer
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
      program urad_phase_main

      use omp_lib
      use uradphasemod
      use wignermod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0),amp(6),ci=(0.0d0,1.0d0),
     &  apoll,apol45,apolh,apolr

      !double precision :: ar(10000),ai(10000)

      double precision, dimension(:), allocatable :: z,y,zm,ym,zprop,yprop,f,ws1,ws2,fg,thez,they,
     &  buffz,buffy,buffe,buffr,buffi,buffu,buffv,buffw,powetot

      double precision, dimension(:,:), allocatable :: s,fzf,stokscr,f2d,fzfprop,f2dprop,
     &  fdspecesour,fdwig,fdwigtzty,
     &  stokesetot,stokespropetot,powe,stosumetot

      double precision, dimension(:,:,:), allocatable :: wigcheck,stosume,
     &  stokese,stokesprope
      double precision, dimension(:,:,:,:), allocatable :: aradf,aradfprop,wig
      double precision, dimension(:,:,:,:,:), allocatable :: wigefold

      complex*16, dimension(:,:), allocatable :: aradscr,esour,esourpin
      complex*16, dimension(:,:,:), allocatable :: esourz,esoury,esourzy2
      complex*16, dimension(:,:,:,:), allocatable :: esourzye,esourzypine
      complex*16, dimension(:,:,:,:,:), allocatable :: esourzy,esourzypin
      complex*16, dimension(:,:,:,:,:,:), allocatable :: esourzypol

      double precision :: banwid=0.001,xbeta=0.0d0,drea,dima,
     &  perlen,shift,ee,ebeam,ebeammean,curr,step,perl,ebeammin,debeam,deltae,g(1000),gsum,
     &  pincen(3),pinw,pinh,park,wlen1,wlen,gamma,
     &  ephmin,ephmax,beffv,beffh,pherror,phgshift,stosum(4),
     &  alphah,alphav,espread,harm1,harm,b0eff,rhv,fsum,ebeamnsig=3.0d0,
     &  betah,betav,pinx,piny,pinz,zz,yy,dzpin,dypin,zmin,ymin,zmax,ymax,
     &  disph,dispph,dispv,disppv,bunchcharge,bunchlen,efi(3),bfi(3),rn(3),
     &  emith,emitv,pinxprop,pinwprop,pinhprop,dy,dz,dthez,dthey,sigz,sigy,
     &  rnsigz=3.0d0,rnsigy=3.0d0,stok(4),stok1,stok2,stok3,stok4,xo,yo,zo,globphase,globphaseprop
     &  ,ajj,gspecnor0,genor0,fdmaxanadist,fdmaxana,fdmax,fluxana,fspecesour,fdwigmax,esourabsmax,
     &  egam,zw,yw,tz,ty,er,ei,specnor_si,gspecnor,genor,dtz,dty,dph,wigsum(4),wigfdmax(4)
     &  ,dum5(5),wigslope(2),wigoffset(2),wigz,wigy,wigzp,wigyp
     &  ,ezr,ezi,eyr,eyi,zob,yob,zp,yp,wigvox,wigflux,wigfluxincut,
     &  wigfluxin,eps=1.0d-12,wigzcut,wigycut,wigmin,wigcenmax,w,fdwigtztymax,w0000,
     &  wigzcutw,wigycutw,depho,wint

      real xran(1),rr(2),
     &  axr,axi,ayr,ayi,azr,azi,
     &  bxr,bxi,byr,byi,bzr,bzi
      real secin,secout

      integer :: idebug=0,noranone,i,ktime=0,iwig,
     &  npiny,npinz,nper,nephogam,nepho,modeold,modeph,modepin,isym,ifixphase,ifold,modesphere,
     &  nharmo,nharm,iy,iz,iobs,kalloerr,nwigpho,
     &  mthreads,nelec,icohere,ihbunch,iepho,iobph,iel,modebunch,ifieldprop,ifieldsym,
     &  modewave=0,isto,nlpoi=0,nobsvprop,npinyprop,npinzprop,iywig,izwig,iypin,izpin,
     &  iwigner,iwignofile,iwigcheck,jz,jy,ic,iobsv,iobfr,nobsv,modepino,
     &  ispline=1,nz,ny,ntz,nty,itz,ity,ianalytic,kpola,kpol,idone,iesourext,
     &  mz,my,ktz,kty,mtz,mty,nzl,nzh,nyl,nyh,kepho,kpola1,kpola2,
     &  ianalytico,ierr,izobs,iyobs,lunwige,lunwig,lunpin,kwigpho,nzyprop,
     &  iwighor,iwigver,nwigcenhit,nwigincut,nosplineefold,iefold,nefold,icbrill,icbrillprop

      namelist/uradphasen/
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,noranone,nobsvprop,npinyprop,npinzprop,iesourext,
     &  pinx,piny,pinz,pinw,pinh,npiny,npinz,modepin,isym,ifixphase,ifold,modesphere,nharm,harm,
     &  nepho,ephmin,ephmax,pherror,phgshift,pinxprop,pinwprop,pinhprop,
     &  mthreads,nelec,icohere,ihbunch,modeph,modebunch,ifieldprop,ifieldsym,iwignofile,
     &  iwigner,
     &  iwigcheck,ianalytic,
     &  betah,betav,alphah,alphav,emith,emitv,espread,wigzp,wigyp,wigz,wigy,
     &  disph,dispph,dispv,disppv,bunchcharge,bunchlen,
     &  nzwig,nywig,thezwig,theywig,nzthewig,nythewig,pinhwig,pinwwig,globphase,globphaseprop,
     &  iwighor,iwigver,
     &  wigzcut,wigycut,wigzcutw,wigycutw,nefold,nosplineefold

      integer, parameter :: nfoldp=16

      integer :: irnsize=64,irnseed(64),ifixseed,irun
      namelist/seedn/irnseed,ifixseed

      integer :: luna,lung,kstat,lunex,lunwigr,lunwigw,lunbun,lunebm,lunflx,lunfld,lunfdp,
     &  lunfldf,lunfdpf,lunwflx,lunwck,lunseed,lunflde,lunfdpe,lunflxe,lunfdpflx,lunfdpflxe

      character(2048) cline

      integer :: ibackspace=8
      character cbs
      equivalence (ibackspace,cbs)

      if (ktime.eq.1) then
        call util_zeit_kommentar_delta(6,'Running urad_phase',1)
      endif

      kalloarad_u=0
      kalloepho_u=0
      kalloobsv_u=0
      kalloaradprop_u=0
      kallostokes_u=0
      kallostokesprop_u=0

      call util_file_delete('urad_phase_prop.flx',kstat)
      call util_file_delete('urad_phase_prop_espread.flx',kstat)
      call util_file_delete('urad_phase.flx',kstat)
      call util_file_delete('urad_phase_espread.flx',kstat)
      call util_file_delete('urad_phase.fld',kstat)
      call util_file_delete('urad_phase_espread.fld',kstat)
      call util_file_delete('urad_phase.fdp',kstat)
      call util_file_delete('urad_phase_espread.fdp',kstat)
      call util_file_delete('urad_phase.wig',kstat)
      call util_file_delete('urad_phase_espread.wig',kstat)

      open(newunit=luna,file='urad_phase.run',status='old',iostat=kstat)
      if (kstat.ne.0) then
        open(newunit=luna,file='urad_phase.run')
        irun=1
      else
        read(luna,*)irun
        rewind(luna)
        irun=irun+1
      endif
      write(luna,*) irun
      flush(luna)
      !close(luna)

      open(newunit=luna,file='urad_phase.nam',status='old',iostat=kstat)
      if (kstat.ne.0) then
        stop "*** Error: Could not open urad_phase.nam"
      endif

      read(luna,uradphasen)
      read(luna,seedn)

      close(luna)

      ianalytic=-ianalytic

      if  (iwigner.eq.0) then
        iwigcheck=0
      else
        if (ianalytic.ne.0) then
          iwigner=-1
        endif
        ifieldprop=1
      endif

      !all util_break

      npiny=(max(1,npiny)/2)*2+1
      npinz=(max(1,npinz)/2)*2+1

      npinyo_u=npiny
      npinzo_u=npinz

      nobsv=npinz*npiny
      nobsv_u=nobsv

      pincen=[pinx,piny,pinz]

      if (npinz.gt.1) then
        zmin=pincen(3)-pinh/2.0d0
        zmax=pincen(3)+pinh/2.0d0
        dzpin=pinh/(npinz-1)
      else
        zmin=pincen(3)
        zmax=pincen(3)
        zmin=pincen(3)
        dzpin=0.0d0
      endif

      if (npiny.gt.1) then
        ymin=pincen(2)-pinh/2.0d0
        ymax=pincen(2)+pinh/2.0d0
        dypin=pinh/(npiny-1)
      else
        ymin=pincen(2)
        ymax=pincen(2)
        ymin=pincen(2)
        dypin=0.0d0
      endif

      if (nepho.le.0) nepho=1
      nepho=(nepho/2)*2+1

      if (nefold.eq.0) nefold=1
      nefold=(nefold/2*2)+1

      if (ianalytic.gt.0 .and. (nefold.gt.1.or.nepho.gt.1)) then
        print*,''
        print*,'*** Warning: Ianalytic > 0 only correct for nefold <= 1 and nepho=1 ***'
        print*,'*** nefold and nepho set one and Epho = (Ephmin+Ephmax)/2.0 ***'
        print*,''
        nefold=1
        nepho=1
        ephmin=(ephmin+ephmax)/2.0d0
        ephmax=ephmin
      endif

      nephogam=nepho*nefold
      ianalytic_u=ianalytic

      if (modeph.ne.0.and.ianalytic.ne.0) then
        write(6,*) '*** Modeph and ianalytic must not be non-zero both ***'
        write(6,*) '*** Modeph overwritten with ianlaytic ***'
      endif

      allocate(esourpin(npinz,npiny),
     &  esourzypine(6,npinz,npiny,nepho),
     &  esourzypin(6,npinz,npiny,nepho,nefold),stosume(4,nepho,nefold),stosumetot(4,nepho))

      esourzypin=(0.0d0,0.0d0)

      allocate(esour(npinzprop,npinyprop),
     &  esourzye(6,npinzprop,npinyprop,nepho),
     &  esourzy(6,npinzprop,npinyprop,nepho,nefold))

      esourzy=(0.0d0,0.0d0)
      esourzye=(0.0d0,0.0d0)

      allocate(
     &  stokese(4,nobsv*nepho,nefold),stokesprope(4,npinzprop*npinyprop*nepho,nefold),
     &  stokesetot(4,nobsv*nepho),stokespropetot(4,npinzprop*npinyprop*nepho))

      allocate(powetot(nobsv),powe(nobsv,nefold))

      if (ianalytic.ne.0) then

        ihbunch=0
        pinxprop=0.0d0

        allocate(arad_u(6,nobsv*nepho),specpow_u(nobsv),stokes_u(4,nobsv*nepho),pow_u(nobsv),
     &    obsv_u(3,nobsv))

        powe=0.0d0
        powetot=0.0d0

        iobsv=0
        yy=ymin-dypin
        do iy=1,npiny
          yy=yy+dypin
          zz=zmin-dzpin
          do iz=1,npinz
            iobsv=iobsv+1
            zz=zz+dzpin
            obsv_u(1,iobsv)=pincen(1)
            obsv_u(2,iobsv)=yy
            obsv_u(3,iobsv)=zz
          enddo
        enddo

        allocate(epho_u(nepho))

        if (nepho.eq.1) then
          epho_u(1)=(ephmin+ephmax)/2.0d0
        else
          depho=(ephmax-ephmin)/dble(nepho-1)
          epho_u(1)=ephmin
          do iepho=2,nepho
            epho_u(iepho)=epho_u(iepho-1)+depho
          enddo
        endif

      else !(ianalytic.eq.0)

        allocate(z(npinz),y(npiny),zm(npinz),ym(npiny))

        y(1)=ymin
        ym(1)=y(1)/1000.0d0
        do iy=2,npiny
          y(iy)=y(iy-1)+dypin
          if(abs(y(iy)).lt.1.0e-9) y(iy)=0.0d0
          ym(iy)=y(iy)/1000.0d0
        enddo

        z(1)=zmin
        zm(1)=z(1)/1000.0d0
        do iz=2,npinz
          z(iz)=z(iz-1)+dzpin
          if(abs(z(iz)).lt.1.0e-9) z(iz)=0.0d0
          zm(iz)=z(iz)/1000.0d0
        enddo

      endif !(ianalytic.eq.0)

      if (nefold.gt.1) then
        ebeammin=ebeam*(1.0d0-ebeamnsig*espread)
        debeam=2.0d0*ebeamnsig*espread/dble(nefold-1)*ebeam
      else
        ebeammin=ebeam
        debeam=0.0d0
      endif

      ebeammean=ebeam
      deltae=espread*ebeam

      icbrill=nobsv/2+1
      icbrillprop=npinzprop*npinyprop/2+1

      gamma=ebeam/emassg1

      if (nharm.gt.0.and.harm.gt.0.0d0) then

        nharmo=nharm
        harm1=harm/nharm
        wlen1=wtoe1/abs(harm/nharm)

        perl=perlen/1000.0d0
        park=2.0d0*(wlen1/(perl*1.0d9/2.0d0/gamma**2)-1.0d0)

        if (park.lt.0.0d0) then
          write(6,*)
     &      '*** Error in urad_phase_main:'
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

      else

        nharmo=0
        perl=perlen/1000.0d0
        b0eff=sqrt(beffh**2+beffv**2)
        park=b0eff*(echarge1*perl/(2.*pi1*emasskg1*clight1))
        wlen1=(1.0d0+park**2/2.0d0)*perl*1.0d9/2.0d0/gamma**2
        harm1=wtoe1/wlen1

        if (nharm.eq.0) then
          nharm=nint(wlen1/(wtoe1/((ephmin+ephmax)/2.0d0)))
        endif

      endif

      defl_u=park

      if (nharmo.gt.0.and.ianalytic.ne.0 .and. (nharm*harm1.lt.ephmin.or.nharm*harm1.gt.ephmax)) then
        write(6,*) '*** Warning: nharm*harm .lt. ephmin .or. nharm*harm.gt.ephmax.'
        write(6,*) '*** This may result in problems for IANALYTIC not zero ***)'
      endif

      nharm_u=nharm

      ebeam=ebeammin
      gamma=ebeam/emassg1

      npinyprop=npinyprop/2*2+1
      npinzprop=npinzprop/2*2+1

      nobsvprop=npinyprop*npinzprop
      dyprop=pinhprop/dble(max(1,npinyprop-1))/1000.0d0
      dzprop=pinwprop/dble(max(1,npinzprop-1))/1000.0d0

      allocate(
     &  obsvzprop_u(npinzprop),obsvyprop_u(npinyprop),obsvprop_u(3,nobsvprop),
     &  aradprop_u(6,nobsvprop*nepho),stokesprop_u(4,nobsvprop*nepho))

      allocate(zprop(npinzprop),yprop(npinyprop),s(npinz,npiny),
     &  f(max(npinz,npiny,npinzprop,npinyprop)),
     &  fg(max(npinz,npiny,npinzprop,npinyprop)),
     &  f2d(npinz,npiny),
     &  f2dprop(npinzprop,npinyprop),
     &  fzf(npinz,npiny),
     &  fzfprop(npinzprop,npinyprop),
     &  aradscr(6,nobsv*nepho),stokscr(4,nobsv*nepho),
     &  aradf(nfoldp,npinz,npiny,nepho),
     &  aradfprop(nfoldp,npinzprop,npinyprop,nepho),
     &  ws1(max(npinz,npiny,npinzprop,npinyprop)),ws2(max(npinz,npiny,npinzprop,npinyprop)))


      obsvyprop_u(1)=-pinhprop/2.0d0/1000.0d0
      yprop(1)=obsvyprop_u(1)*1000.0d0
      do iy=2,npinyprop
        obsvyprop_u(iy)=obsvyprop_u(iy-1)+dyprop
        if (abs(obsvyprop_u(iy)).lt.1.0d-12) obsvyprop_u(iy)=0.0d0
        yprop(iy)=obsvyprop_u(iy)*1000.0d0
      enddo

      obsvzprop_u(1)=-pinwprop/2.0d0/1000.0d0
      zprop(1)=obsvzprop_u(1)*1000.0d0
      do iz=2,npinzprop
        obsvzprop_u(iz)=obsvzprop_u(iz-1)+dzprop
        if (abs(obsvzprop_u(iz)).lt.1.0d-12) obsvzprop_u(iz)=0.0d0
        zprop(iz)=obsvzprop_u(iz)*1000.0d0
      enddo

      iobs=0
      do iy=1,npinyprop
        do iz=1,npinzprop
          iobs=iobs+1
          obsvprop_u(1,iobs)=pinxprop/1000.0d0
          obsvprop_u(2,iobs)=obsvyprop_u(iy)
          obsvprop_u(3,iobs)=obsvzprop_u(iz)
        enddo
      enddo

      nzthewig=(max(1,nzthewig)/2)*2+1
      nythewig=(max(1,nythewig)/2)*2+1

      nz=npinzprop
      ny=npinyprop
      ntz=nzthewig
      nty=nythewig
      nzyprop=nz*ny

      if (iwigner.ne.0) then
        allocate(esourz(nz,ny,nepho),esoury(nz,ny,nepho),esourzy2(2,nz,ny),
     &    wig(nz,ny,ntz,nty))
        allocate(wigefold(nz,ny,ntz,nty,nefold))
c     &      esourzypol(2,nz,ny,nepho,nefold))
      endif

      !all util_break

      allocate(thez(ntz),they(nty),
     &  fdwig(nz,ny),fdwigtzty(ntz,nty),
     &  buffz(max(nz,ntz)),buffy(max(ny,nty)),buffe(nefold),
     &  buffr(max(nz,ntz)),buffi(max(nz,ntz)),
     &  buffw(max(nz,ntz)),
     &  buffu(max(nz,ntz)),buffv(max(nz,ntz)))

      if (ntz.gt.1) then
        dthez=thezwig/dble(max(1,ntz-1))/1000.0d0
        thez(1)=-thezwig/2.0d0/1000.0d0
        do iz=2,ntz
          thez(iz)=thez(iz-1)+dthez
          if (abs(thez(iz)).lt.1.0d-12) thez(iz)=0.0d0
        enddo
      else
        dthez=1.0d0
        thez(1)=0.0d0
      endif

      if (nty.gt.1) then
        dthey=theywig/dble(max(1,nty-1))/1000.0d0
        they(1)=-theywig/2.0d0/1000.0d0
        do iy=2,nythewig
          they(iy)=they(iy-1)+dthey
          if (abs(they(iy)).lt.1.0d-12) they(iy)=0.0d0
        enddo
      else
        dthey=1.0d0
        they(1)=0.0d0
      endif

      wigoffset(1)=wigz
      wigoffset(2)=wigy
      wigslope(1)=wigzp
      wigslope(2)=wigyp

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid !BW
     &  /1.0d6 !m**s -> mm**2

      globphase_u=globphase
      globphaseprop_u=globphaseprop

      rn=0.0d0

      modepino=modepin
      if (ifold.ne.0) modepin=2
      ifixphase_u=ifixphase
      ifold_u=ifold

      pinxprop_u=pinxprop
      pinwprop_u=pinwprop
      pinhprop_u=pinhprop

      npinyprop_u=npinyprop
      npinzprop_u=npinzprop

      nobsvprop_u=nobsvprop

      ifieldsym_u=ifieldsym

      npinz_u=npinz
      npiny_u=npiny
      nobsv_u=nobsv

      pinw_u=pinw
      pinh_u=pinh

      pinxprop_u=pinxprop
      pinwprop_u=pinwprop
      pinhprop_u=pinhprop

      curr_u=curr
      banwid_u=banwid

      ifieldprop_u=ifieldprop
      modepin_u=modepin

      !all util_break

      nywig=npinyprop
      nzwig=npinzprop

      nepho_u=nepho
      nepho_u=max(1,nepho_u)

      emith=emith*1.0d-9
      emitv=emitv*1.0d-9

      bunchlen=bunchlen/1.0d9 !nm->m

      if (ifixseed.ne.0) then
        ifixseed=1
        call util_random_set_seed(irnsize,irnseed)
      endif

      !all util_break

      nelec_u=nelec
      ihbunch_u=ihbunch

      open(newunit=luna,file='urad_phase.nor')
      write(luna,*)curr,banwid
      close(luna)

      open(newunit=lunpin,file='urad_phase.pin')
      write(luna,*)npinz,npiny,pinw,pinh
      write(luna,*)pincen
      write(luna,*)modepino,ifold,ifixphase,ifieldprop,nelec,ihbunch,-ianalytic
      write(luna,*)betah,emith,betav,emitv,espread
      write(luna,*)npinzprop,npinyprop,pinxprop,pinwprop,pinhprop
      write(luna,*)nzthewig,nythewig,thezwig,theywig
      write(luna,*)iwigner,iwignofile,nefold
      close(luna)

      if (mthreads.lt.0) then
        mthreads=OMP_GET_MAX_THREADS()
      else if (mthreads.eq.0) then
        mthreads=1
      endif

      mthreads_u=mthreads

      gsum=0.0d0
      do iefold=1,nefold
        if (deltae.ne.0.0d0) then
          g(iefold)=exp(-((ebeammin+(iefold-1)*debeam-ebeammean)/deltae)**2/2.0d0)/sqrt(twopi1)/deltae
        else
          g(iefold)=1.0d0
        endif
        gsum=gsum+g(iefold)
      enddo

      if (nefold.le.3.or.nosplineefold.ne.0) then
        g=g/gsum
        gsum=1.0d0
      endif

      ebeam=ebeammin-debeam

      if (ianalytic.le.0) then
        open(newunit=lunfld,file='urad_phase.fld')
        open(newunit=lung,file='urad_phase.geo')
        open(newunit=lunflx,file='urad_phase.flx')
      endif

      if (ifieldprop.ne.0.or.ianalytic.gt.0) then
        open(newunit=lunfdp,file='urad_phase.fdp')
        open(newunit=lunfdpflx,file='urad_phase_prop.flx')
      endif

      do iefold=1,nefold

        ebeam=ebeam+debeam
        buffe(iefold)=ebeam
        gamma=ebeam/emassg1

        if (nefold.gt.1) then
          write(cline,*)'Ebeam loop for field calculations:',iefold,' of ',nefold,sngl(ebeam)
          write(6,'(a)',advance='no') cline(1:len_trim(cline))
          do ic=1,len_trim(cline)
            write(6,'(a)',advance='no') cbs
          enddo
        endif

        !allutil_break
        do iepho=1,nepho

          if (ianalytic.gt.0) then

            ifieldprop=0
            pinxprop=0.0d0

            call undulator_source_analytic(curr,banwid,park,gamma,npinzprop,npinyprop,
     &        pinxprop,pinwprop/1000.0d0,pinhprop/1000.0d0,epho_u(iepho)/hbarev1,
     &        perlen/1000.0d0,
     &        nper,0.0d0,wigoffset/1000.0d0,wigslope/1000.0d0,esour,nharm,ajj,
     &        gspecnor,genor,fdmaxana,fluxana)

            esourzy(3,:,:,iepho,iefold)=esour(:,:)

            iobs=0
            do iy=1,npinyprop
              do iz=1,npinzprop
                iobs=iobs+1
                iobph=iobs+nzyprop*(iepho-1)
                aradprop_u(3,iobph)=esourzy(3,iz,iy,iepho,iefold)
                call util_e_to_stokes(aradprop_u(1:3,iobph),specnor_si,stokesprop_u(:,iobph))
                stokesprope(:,iobph,iefold)=stokesprop_u(:,iobph)
              enddo
            enddo

          else if (ianalytic.lt.0) then

            call undulator_source_analytic(curr,banwid,park,gamma,npinz,npiny,
     &        pinx/1000.0d0,pinw/1000.0d0,pinh/1000.0d0,epho_u(iepho)/hbarev1,perlen/1000.0d0,
     &        nper,pinx/1000.0d0,wigoffset/1000.0d0,wigslope/1000.0d0,esourpin,nharm,ajj,
     &        gspecnor,genor,fdmaxana,fluxana)

c              print*,"fdmaxana:",sngl(fdmaxana)

            esourzypin(3,:,:,iepho,iefold)=esourpin(:,:)
            !allutil_break
            iobs=0
            do iy=1,npiny
              do iz=1,npinz
                iobs=iobs+1
                iobph=iobs+nobsv_u*(iepho-1)
                arad_u(3,iobph)=esourpin(iz,iy)
                stokes_u(1,iobph)=abs(esourpin(iz,iy))**2*specnor_si
                stokese(1,iobph,iefold)=stokes_u(1,iobph)
              enddo
            enddo
            !allutil_break

          endif !ianalytic

        enddo !nepho

        if (ianalytic.eq.0) then

          call urad_phase(
     &      mthreads,nelec,noranone,icohere,modebunch,bunchlen,bunchcharge,ihbunch,
     &      perlen,shift,nper,beffv,beffh,
     &      ebeam,curr,step,nlpoi,
     &      pincen,pinw,pinh,npiny,npinz,modepin,modesphere,
     &      nepho,ephmin,ephmax,banwid,
     &      xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &      disph,dispph,dispv,disppv,
     &      modeph,pherror,phgshift,modewave
     &      )

          powe(:,iefold)=pow_u(iefold)

          do iepho=1,nepho
            iobs=0
            do iy=1,npiny
              do iz=1,npinz
                iobs=iobs+1
                iobph=iobs+nzyprop*(iepho-1)
                esourzypin(:,iz,iy,iepho,iefold)=arad_u(:,iobph)
              enddo
            enddo
          enddo

          stokese(:,:,iefold)=stokes_u(:,:)
          !all util_break

        endif !ianalytic.eq.0

        if (ianalytic.le.0) then

c          print*,"S0_max:",sngl(maxval(stokes_u))

          !allutil_break

          do iepho=1,nepho
            do iobs=1,nobsv_u

              iobph=iobs+nobsv_u*(iepho-1)

              if (ianalytic.eq.0) then
                rn(1)=
     &            real(arad_u(2,iobph)*conjg(arad_u(6,iobph))-arad_u(3,iobph)*conjg(arad_u(5,iobph)))
                rn(2)=
     &            real(arad_u(3,iobph)*conjg(arad_u(4,iobph))-arad_u(1,iobph)*conjg(arad_u(6,iobph)))
                rn(3)=
     &            real(arad_u(1,iobph)*conjg(arad_u(5,iobph))-arad_u(2,iobph)*conjg(arad_u(4,iobph)))
                rn=rn/norm2(rn)
              else
                rn=0.0d0
              endif

              axr=real(arad_u(1,iobph))
              axi=imag(arad_u(1,iobph))
              ayr=real(arad_u(2,iobph))
              ayi=imag(arad_u(2,iobph))
              azr=real(arad_u(3,iobph))
              azi=imag(arad_u(3,iobph))

              write(lunfld,'(3(1pe17.8e3),2i10,30(1pe17.8e3))')
     &          obsv_u(1:3,iobs),
     &          iepho,iefold,
     &          epho_u(iepho),buffe(iefold),stokes_u(1:4,iobph),pow_u(iobs),
     &          real(arad_u(1,iobph)),imag(arad_u(1,iobph)),
     &          real(arad_u(2,iobph)),imag(arad_u(2,iobph)),
     &          real(arad_u(3,iobph)),imag(arad_u(3,iobph)),
     &          real(arad_u(4,iobph)),imag(arad_u(4,iobph)),
     &          real(arad_u(5,iobph)),imag(arad_u(5,iobph)),
     &          real(arad_u(6,iobph)),imag(arad_u(6,iobph)),
     &          rn,g(iefold)

            enddo !iobs=1,nobsv_u

          enddo !iepho

        endif !(ianalytic.le.0) then

        !all util_break
        if (ifieldprop.ne.0.and.ianalytic.le.0) then

          call urad_phase_prop(mthreads)

          stokesprope(:,:,iefold)=stokesprop_u(:,:)

          do iepho=1,nepho
            iobs=0
            do iy=1,npinyprop
              do iz=1,npinzprop
                iobs=iobs+1
                iobph=iobs+nobsvprop_u*(iepho-1)
                esourzy(:,iz,iy,iepho,iefold)=aradprop_u(:,iobph)
              enddo
            enddo
          enddo

        endif

        if (ifieldprop.ne.0.or.ianalytic.gt.0) then

          do iobs=1,nobsvprop_u
            do iepho=1,nepho

              iobph=iobs+nobsvprop_u*(iepho-1)

              ! rnx = (eyr+i*eyi)*(bzr-i*bzi) - (ezr+i*ezi)*(byr-i*byi)
              !     = eyr*bzr - i*eyr*bzi + i*eyi*bzr + eyi*bzi
              !     - ezr*byr + i*erz*byi - i*ezi*byr - ezi*byi
              ! real(rnx) = eyr*bzr + eyi*bzi - ezr*byr - ezi*byi

              axr=real(aradprop_u(1,iobph))
              axi=imag(aradprop_u(1,iobph))
              ayr=real(aradprop_u(2,iobph))
              ayi=imag(aradprop_u(2,iobph))
              azr=real(aradprop_u(3,iobph))
              azi=imag(aradprop_u(3,iobph))

              if (ianalytic.eq.0) then
                if (ifieldprop.gt.0) then
                  rn(1)=real(aradprop_u(2,iobph)*conjg(aradprop_u(6,iobph))-aradprop_u(3,iobph)*conjg(aradprop_u(5,iobph)))
                  rn(2)=real(aradprop_u(3,iobph)*conjg(aradprop_u(4,iobph))-aradprop_u(1,iobph)*conjg(aradprop_u(6,iobph)))
                  rn(3)=real(aradprop_u(1,iobph)*conjg(aradprop_u(5,iobph))-aradprop_u(2,iobph)*conjg(aradprop_u(4,iobph)))
                  rn=rn/norm2(rn)
                else if (ifieldprop.eq.-1) then
                  rn(1:3)=aradprop_u(4:6,iobs)
                endif
              endif

              bxr=real(aradprop_u(4,iobph))
              bxi=imag(aradprop_u(4,iobph))
              byr=real(aradprop_u(5,iobph))
              byi=imag(aradprop_u(5,iobph))
              bzr=real(aradprop_u(6,iobph))
              bzi=imag(aradprop_u(6,iobph))

              write(lunfdp,'(3(1pe15.6e3),2i10,30(1pe15.6e3))')
     &          obsvprop_u(1:3,iobs)*1000.0d0,
     &          iepho,iefold,
     &          epho_u(iepho),buffe(iefold),
     &          stokesprop_u(1:4,iobph),
     &          axr,axi,ayr,ayi,azr,azi,
     &          bxr,bxr,byr,byi,bzr,bzi,
     &          rn,g(iefold)

            enddo !nepho

          enddo !iobs=1,nobsvprop_u

          do iepho=1,nepho
            do isto=1,4
              iobs=0
              do iz=1,npinzprop_u
                do iy=1,npinyprop_u
                  iobs=iobs+1
                  iobph=iobs+nzyprop*(iepho-1)
                  s(iz,iy)=stokesprope(isto,iobph,iefold)
                enddo
              enddo
              callutil_break
              call util_spline_integral_2d(npinzprop_u,npinyprop_u,zprop,yprop,s,
     &          stosume(isto,iepho,iefold),kstat)
            enddo !isto

            write(lunfdpflx,*)iepho,iefold,epho_u(iepho),buffe(iefold),stosume(:,iepho,iefold),
     &        g(iefold)
          enddo !nepho

c          print*,"S0_prop_max:",sngl(maxval(stokesprop_u))

        endif !ifieldprop

        !allutil_break
        if (ianalytic.le.0) then

          do iepho=1,nepho

            if (modepin.eq.0.and.npinz_u.ge.3.and.npiny_u.ge.3) then

              do isto=1,4
                iobs=0
                do iz=1,npinz_u
                  do iy=1,npiny_u
                    iobs=iobs+1
                    iobph=iobs+nobsv_u*(iepho-1)
                    s(iz,iy)=stokes_u(isto,iobph)
                  enddo
                enddo
                call util_spline_integral_2d(npinz_u,npiny_u,z,y,s,stosum(isto),kstat)
                stosume(isto,iepho,iefold)=stosum(isto)
              enddo !isto

              write(lunflx,*)iepho,iefold,epho_u(iepho),buffe(iefold),stosume(:,iepho,iefold),
     &          g(iefold)

            else

              iobs=0
              stosum=0.0d0
              do iz=1,npinz_u
                do iy=1,npiny_u
                  iobs=iobs+1
                  iobph=iobs+nobsv_u*(iepho-1)
                  stosum=stosum+stokes_u(1:4,iobph)
                enddo
              enddo
              write(lunflx,*)iepho,iefold,epho_u(iepho),buffe(iefold),stosum/nobsv_u*pinw*pinh
            endif
          enddo

        endif !(ianalytic.le.0) then

      enddo !nefold

      esourzypine=(0.0d0,0.0d0)
      stokesetot=0.0d0

      !all util_break
      if (ianalytic.eq.0) then

        iobs=0
        do iy=1,npiny
          do iz=1,npinz
            iobs=iobs+1
            do iefold=1,nefold
              if (nefold.le.3.or.nosplineefold.ne.0) then
                powetot(iobs)=powetot(iobs)+powe(iobs,iefold)*g(iefold)
              else
                buffr(iefold)=powe(iobs,iefold)*g(iefold)
              endif
            enddo !iefold
            if (nefold.gt.3.and.nosplineefold.eq.0) then
              call util_integral_spline(buffe,buffr,nefold,powetot(iobs))
            endif
          enddo
        enddo
      endif

      !allutil_break

      if (nefold.gt.1) then

        if (ianalytic.le.0) then

          open(newunit=lunflde,file='urad_phase_espread.fld')
          open(newunit=lunflxe,file='urad_phase_espread.flx')

          do iepho=1,nepho

            do iy=1,npiny
              do iz=1,npinz
                do i=1,6
                  do iefold=1,nefold
                    if (nefold.le.3.or.nosplineefold.ne.0) then
                      esourzypine(i,iz,iy,iepho)=esourzypine(i,iz,iy,iepho)+
     &                  esourzypin(i,iz,iy,iepho,iefold)*g(iefold)
                    else
                      buffr(iefold)=dreal(esourzypin(i,iz,iy,iepho,iefold))*g(iefold)
                      buffi(iefold)=dimag(esourzypin(i,iz,iy,iepho,iefold))*g(iefold)
                    endif
                  enddo !iefold
                  if (nefold.gt.3.and.nosplineefold.eq.0) then
                    call util_integral_spline(buffe,buffr,nefold,drea)
                    call util_integral_spline(buffe,buffr,nefold,dima)
                    esourzypine(i,iz,iy,iepho)=dcmplx(drea,dima)
                  endif
                enddo !i
              enddo !iz
            enddo !iy

            iobs=0
            do iy=1,npiny
              do iz=1,npinz
                iobs=iobs+1
                iobph=iobs+nobsv_u*(iepho-1)
                do isto=1,4
                  do iefold=1,nefold
                    if (nefold.le.3.or.nosplineefold.ne.0) then
                      stokesetot(isto,iobph)=stokesetot(isto,iobph)+
     &                  stokese(isto,iobph,iefold)
                    else
                      buffz(iefold)=stokese(isto,iobph,iefold)
                    endif
                  enddo !iefold
                  if (nefold.gt.3.and.nosplineefold.eq.0) then
                    call util_integral_spline(buffe,buffz,nefold,stokesetot(isto,iobph))
                  endif
                enddo !isto
              enddo
            enddo

            do isto=1,4
              do iefold=1,nefold
                if (nefold.le.3.or.nosplineefold.ne.0) then
                  stosumetot(isto,iepho)=stosume(isto,iepho,iefold)+
     &              stosume(isto,iepho,iefold)
                else
                  buffz(iefold)=stosume(isto,iepho,iefold)
                endif
              enddo !iefold
              if (nefold.gt.3.and.nosplineefold.eq.0) then
                call util_integral_spline(buffe,buffz,nefold,stosumetot(isto,iepho))
              endif
            enddo !isto

            iobs=0
            do iy=1,npiny
              do iz=1,npinz

                iobs=iobs+1
                iobph=iobs+nobsv_u*(iepho-1)

                if (ianalytic.eq.0) then
                  rn(1)=
     &              real(esourzypine(2,iz,iy,iepho)*conjg(esourzypine(6,iz,iy,iepho))-esourzypine(3,iz,iy,iepho)*conjg(esourzypine(5,iz,iy,iepho)))
                  rn(2)=
     &              real(esourzypine(3,iz,iy,iepho)*conjg(esourzypine(4,iz,iy,iepho))-esourzypine(1,iz,iy,iepho)*conjg(esourzypine(6,iz,iy,iepho)))
                  rn(3)=
     &              real(esourzypine(1,iz,iy,iepho)*conjg(esourzypine(5,iz,iy,iepho))-esourzypine(2,iz,iy,iepho)*conjg(esourzypine(4,iz,iy,iepho)))
                  rn=rn/norm2(rn)
                else
                  rn=0.0d0
                endif

                axr=real(esourzypine(1,iz,iy,iepho))
                axi=imag(esourzypine(1,iz,iy,iepho))
                ayr=real(esourzypine(2,iz,iy,iepho))
                ayi=imag(esourzypine(2,iz,iy,iepho))
                azr=real(esourzypine(3,iz,iy,iepho))
                azi=imag(esourzypine(3,iz,iy,iepho))

                write(lunflde,'(3(1pe17.8e3),2i10,30(1pe17.8e3))')
     &            obsv_u(1:3,iobs),
     &            iepho,-nefold,
     &            epho_u(iepho),ebeammean,stokesetot(1:4,iobph),powe(iobs,iepho),
     &            real(esourzypine(1,iz,iy,iepho)),imag(esourzypine(1,iz,iy,iepho)),
     &            real(esourzypine(2,iz,iy,iepho)),imag(esourzypine(2,iz,iy,iepho)),
     &            real(esourzypine(3,iz,iy,iepho)),imag(esourzypine(3,iz,iy,iepho)),
     &            real(esourzypine(4,iz,iy,iepho)),imag(esourzypine(4,iz,iy,iepho)),
     &            real(esourzypine(5,iz,iy,iepho)),imag(esourzypine(5,iz,iy,iepho)),
     &            real(esourzypine(6,iz,iy,iepho)),imag(esourzypine(6,iz,iy,iepho)),
     &            rn

              enddo !iz
            enddo !iy

            write(lunflxe,*)iepho,-nefold,epho_u(iepho),ebeammean,stosumetot(:,iepho)

          enddo !iepho
        endif !(ianalytic.le.0)

        if (ifieldprop.gt.0.or.ianalytic.gt.0) then

          open(newunit=lunfdpe,file='urad_phase_espread.fdp')
          open(newunit=lunfdpflxe,file='urad_phase_prop_espread.flx')

          !allutil_break
          do iepho=1,nepho

            iobs=0
            do iy=1,npinyprop
              do iz=1,npinzprop
                iobs=iobs+1
                iobph=iobs+npinyprop*npinzprop*(iepho-1)
                do isto=1,4
                  do iefold=1,nefold
                    if (nefold.le.3.or.nosplineefold.ne.0) then
                      stokespropetot(isto,iobph)=stokespropetot(isto,iobph)+
     &                  stokesprope(isto,iobph,iefold)*g(iefold)*deltae
                    else
                      buffz(iefold)=stokesprope(isto,iobph,iefold)*g(iefold)
                    endif
                  enddo !iefold
                  if (nefold.gt.3.and.nosplineefold.eq.0) then
                    call util_integral_spline(buffe,buffz,nefold,stokespropetot(isto,iobph))
                  endif
                enddo !isto
              enddo
            enddo

            do iy=1,npinyprop
              do iz=1,npinzprop
                do i=1,6
                  do iefold=1,nefold
                    if (nefold.le.3.or.nosplineefold.ne.0) then
                      esourzye(i,iz,iy,iepho)=esourzye(i,iz,iy,iepho)+
     &                  esourzy(i,iz,iy,iepho,iefold)
                    else
                      buffr(iefold)=dreal(esourzy(i,iz,iy,iepho,iefold))*g(iefold)
                      buffi(iefold)=dimag(esourzy(i,iz,iy,iepho,iefold))*g(iefold)
                    endif
                  enddo !iefold
                  if (nefold.gt.3.and.nosplineefold.eq.0) then
                    call util_integral_spline(buffe,buffr,nefold,drea)
                    call util_integral_spline(buffe,buffi,nefold,dima)
                    esourzye(i,iz,iy,iepho)=dcmplx(drea,dima)
                  endif
                enddo !i=1,6
              enddo !iz
            enddo !iy

            !allutil_break
            iobs=0
            do iy=1,npinyprop
              do iz=1,npinzprop

                iobs=iobs+1
                iobph=iobs+nzyprop*(iepho-1)

                if (ianalytic.eq.0) then
                  rn(1)=
     &              real(esourzye(2,iz,iy,iepho)*conjg(esourzye(6,iz,iy,iepho))-esourzye(3,iz,iy,iepho)*conjg(esourzye(5,iz,iy,iepho)))
                  rn(2)=
     &              real(esourzye(3,iz,iy,iepho)*conjg(esourzye(4,iz,iy,iepho))-esourzye(1,iz,iy,iepho)*conjg(esourzye(6,iz,iy,iepho)))
                  rn(3)=
     &              real(esourzye(1,iz,iy,iepho)*conjg(esourzye(5,iz,iy,iepho))-esourzye(2,iz,iy,iepho)*conjg(esourzye(4,iz,iy,iepho)))
                  rn=rn/norm2(rn)
                else
                  rn=0.0d0
                endif

                axr=real(esourzye(1,iz,iy,iepho))
                axi=imag(esourzye(1,iz,iy,iepho))
                ayr=real(esourzye(2,iz,iy,iepho))
                ayi=imag(esourzye(2,iz,iy,iepho))
                azr=real(esourzye(3,iz,iy,iepho))
                azi=imag(esourzye(3,iz,iy,iepho))

                write(lunfdpe,'(3(1pe17.8e3),2i10,30(1pe17.8e3))')
     &            obsvprop_u(1:3,iobs)*1000.0d0,
     &            iepho,-nefold,
     &            epho_u(iepho),ebeammean,stokespropetot(1:4,iobph),
     &            real(esourzye(1,iz,iy,iepho)),imag(esourzye(1,iz,iy,iepho)),
     &            real(esourzye(2,iz,iy,iepho)),imag(esourzye(2,iz,iy,iepho)),
     &            real(esourzye(3,iz,iy,iepho)),imag(esourzye(3,iz,iy,iepho)),
     &            real(esourzye(4,iz,iy,iepho)),imag(esourzye(4,iz,iy,iepho)),
     &            real(esourzye(5,iz,iy,iepho)),imag(esourzye(5,iz,iy,iepho)),
     &            real(esourzye(6,iz,iy,iepho)),imag(esourzye(6,iz,iy,iepho)),
     &            rn

              enddo !iz
            enddo !iy

            do isto=1,4
              iobs=0
              do iz=1,npinzprop_u
                do iy=1,npinyprop_u
                  iobs=iobs+1
                  iobph=iobs+nobsvprop_u*(iepho-1)
                  s(iz,iy)=stokespropetot(isto,iobph)
                enddo
              enddo
              call util_spline_integral_2d(npinzprop,npinyprop,zprop,yprop,s,stosum(isto),kstat)
              stosumetot(isto,iepho)=stosum(isto)
            enddo !isto

            write(lunfdpflxe,*)iepho,iefold,epho_u(iepho),ebeammean,stosumetot(:,iepho)

          enddo !iepho

        endif !ifieldprop

      endif !nefold

      if (iwigner.gt.-5.and.iwigner.ne.0) then

        write(6,*) ''
        write(6,*) 'Starting calculation of Wigner-Distributions'
        write(6,*) ''

        nz=NpinZprop
        ny=NpinYprop
        ntz=nzthewig
        nty=nythewig

        wigvox=dzprop*dyprop*dthey*dthez  !??*4.0d0 ! Factor 4 due to change in integration vari.

        open(newunit=lunwigw,file='urad_phase.wig')

        if (nefold.gt.1) then
          open(newunit=lunwige,file='urad_phase_espread.wig')
        endif

        kpola1=1
        kpola2=4

        if (iwigner.lt.0) then
          kpola1=-iwigner
          kpola2=kpola1
        endif

        do kpola=kpola1,kpola2

          do iepho=1,nepho

            do iefold=1,nefold

              !allutil_break

              if (kpola.eq.1) then
                if (iwigner.eq.-2.or.iwigner.eq.-3.or.iwigner.eq.-4) cycle
                esourzy2(1,:,:)=esourzy(3,:,:,iepho,iefold)
                esourzy2(2,:,:)=conjg(esourzy(3,:,:,iepho,iefold))
              else if (kpola.eq.2) then
                if (iwigner.eq.-1.or.iwigner.eq.-3.or.iwigner.eq.-4) cycle
                esourzy2(1,:,:)=esourzy(2,:,:,iepho,iefold)
                esourzy2(2,:,:)=conjg(esourzy(2,:,:,iepho,iefold))
              else if (kpola.eq.3) then
                if (iwigner.eq.-1.or.iwigner.eq.-2.or.iwigner.eq.-4) cycle
                esourzy2(1,:,:)=esourzy(3,:,:,iepho,iefold)
                esourzy2(2,:,:)=conjg(esourzy(2,:,:,iepho,iefold))
              else if (kpola.eq.4) then
                if (iwigner.eq.-1.or.iwigner.eq.-2.or.iwigner.eq.-3) cycle
                esourzy2(1,:,:)=esourzy(2,:,:,iepho,iefold)
                esourzy2(2,:,:)=conjg(esourzy(3,:,:,iepho,iefold))
              endif

              call undulator_wigner_num(npinzprop,npinyprop,dzprop,dyprop,wtoe1/epho_u(iepho),
     &          esourzy2,ntz,nty,thez,they,wig,curr,banwid)

              wigefold(:,:,:,:,iefold)=wig(:,:,:,:)

              if (iwignofile.eq.0.and.iefold.eq.nefold/2+1) then
                do iy=1,ny
                  do iz=1,nz
                    do itz=1,ntz
                      do ity=1,nty
                        write(lunwigw,'(5i10,5(1pe15.6e3),2i10,10(1pe15.6e3))')
     &                    kpola,iz,iy,itz,ity,
     &                    obsvprop_u(1,1),
     &                    yprop(iy),zprop(iz),they(ity)*1000.0d0,thez(itz)*1000.0d0,
     &                    iepho,iefold,epho_u(iepho),buffe(iefold),
     &                    dreal(esourzy2(1,iz,iy)),dimag(esourzy2(1,iz,iy)),
     &                    dreal(esourzy2(2,iz,iy)),dimag(esourzy2(2,iz,iy)),
     &                    wig(iz,iy,itz,ity)/1.0d12,g(iefold)
                      enddo
                    enddo
                  enddo !nz
                enddo !ny
              endif !iwignofile

            enddo !iefold=1,nefold

            if (nefold.gt.1) then
              do ity=1,nty
                do itz=1,ntz
                  do iy=1,ny
                    do iz=1,nz
                      do iefold=1,nefold
                        buffw(iefold)=wigefold(iz,iy,itz,ity,iefold)*g(iefold)
                      enddo
                      if (nefold.le.3) then
                        wint=sum(buffw)*deltae
                      else
                        call util_integral_spline(buffe,buffw,nefold,wint)
                      endif
                      write(lunwige,'(5i10,5(1pe15.6e3),2i10,11(1pe15.6e3))')
     &                  kpola,iz,iy,itz,ity,
     &                  obsvprop_u(1,1),
     &                  yprop(iy),zprop(iz),they(ity)*1000.0d0,thez(itz)*1000.0d0,
     &                  iepho,-nefold,epho_u(iepho),ebeammean,
     &                  dreal(esourzye(2,iz,iy,iepho)),dimag(esourzye(2,iz,iy,iepho)),
     &                  dreal(esourzye(3,iz,iy,iepho)),dimag(esourzye(3,iz,iy,iepho)),
     &                  wint/1.0d12
                    enddo !iz
                  enddo !iy
                enddo !itz
              enddo !ity
            endif !(nefold.gt.1) then

          enddo !iepho

        enddo !kpola

      endif !iwigner

      print*,' '
      print*,' '

      end
*CMZ :          29/09/2025  12.35.09  by  Michael Scheer
*CMZ :  4.02/00 13/09/2025  10.16.17  by  Michael Scheer
*CMZ :  4.01/07 13/08/2024  10.11.51  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  10.49.56  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase(
     &  mthreads,nelec,noranone,icohere,modebunch,bunchlen,bunchcharge,ihbunch,
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,nlpoi,
     &  pincen,pinw,pinh,npiny,npinz,modepin,modesphere,
     &  nepho,ephmin,ephmax,banwid,
     &  xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  modeph,pherror,phgshift,modewave
     &  )

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      double complex :: rea(3),expsh,esour(max(npinz_u,npinzprop_u),max(npiny_u,npinyprop_u))

      double precision
     &  perlen,shift,ebeam,curr,step,banwid,
     &  pincen(3),pinw,pinh,betah,alphah,betav,alphav,
     &  ephmin,ephmax,beffv,beffh,pherror,phgshift,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,y,z,dy,dz,ymin,zmin,bunchlen,bunchcharge,
     &  xbeta,df,xx,yy,zz,r,xn,yn,zn,h2
     &  ,ajj,enor,fdmax,flux

      integer :: ktime=0,ical=0,
     &  npiny,npinz,nper,nepho,mthreads,nelec,icohere,ihbunch,i,nlpoi,
     &  modeph,modepin,modesphere,modebunch,iy,iz,iobsv,noranone,modewave,
     &  icbrill,iobs,iobfr,ifrq,kalloerr

      save ical

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase',1)

      mthreads_u=mthreads

      nelec_u=nelec
      noranone_u=noranone
      icohere_u=icohere
      modebunch=modebunch_u
      bunchlen_u=bunchlen
      bunchcharge_u=bunchcharge
      ihbunch_u=ihbunch

      perlen_u=perlen/1000.0d0
      shift_u=shift/1000.0d0
      nper_u=nper
      beffv_u=beffv
      beffh_u=beffh

      ebeam_u=ebeam
      gamma_u=ebeam_u/emassg1
      step_u=step/1000.0d0
      nstep_u=max(1,nint(perlen_u/step_u))

      curr_u=curr

      pincen_u=pincen/1000.0d0
      pinw_u=pinw/1000.0d0
      pinh_u=pinh/1000.0d0
      npiny_u=npiny
      npinz_u=npinz
      modepin_u=modepin

      ephmin_u=ephmin_u
      ephmax_u=ephmax_u
      banwid_u=banwid
      nepho_u=nepho

      npiny_u=max(1,npiny_u)
      npinz_u=max(1,npinz_u)

      nlpoi_u=nlpoi

      if (modepin.ne.1) then
        nobsv_u=npiny_u*npinz_u
      else
        npinz_u=1
        npiny_u=1
        nobsv_u=1
      endif

      if (ical.eq.0) then

        allocate(epho_u(nepho),obsv_u(3,nobsv_u),
     &    arad_u(6,nobsv_u*nepho_u),
     &    specpow_u(nobsv_u),
     &    stokes_u(4,nobsv_u*nepho_u),pow_u(nobsv_u))

        if (ihbunch_u.gt.0) then
          allocate(fbunch_u(41,nelec_u/ihbunch_u*nepho_u))
          fbunch_u=0.0d0
        else if (ihbunch_u.lt.0) then
          allocate(fbunch_u(41,nobsv_u*nelec_u/(-ihbunch_u)*nepho_u))
          fbunch_u=0.0d0
        endif
      endif

      stokes_u=0.0d0
      specpow_u=0.0d0
      arad_u=(0.0d0,0.0d0)
      pow_u=0.0d0

      if (npiny_u.eq.1) then
        dy=0.0d0
        ymin=pincen_u(2)
      else
        dy=pinh_u/(npiny_u-1)
        ymin=pincen_u(2)-pinh_u/2.0d0
      endif

      if (npinz_u.eq.1) then
        dz=0.0d0
        zmin=pincen_u(3)
      else
        dz=pinw_u/(npinz_u-1)
        zmin=pincen_u(3)-pinw_u/2.0d0
      endif

      iobsv=0
      y=ymin-dy
      do iy=1,npiny_u
        y=y+dy
        z=zmin-dz
        do iz=1,npinz_u
          iobsv=iobsv+1
          z=z+dz
          obsv_u(1,iobsv)=pincen_u(1)
          obsv_u(2,iobsv)=y
          obsv_u(3,iobsv)=z
          if (modesphere.ne.0) then
            !all util_break
            xx=obsv_u(1,iobsv)
            yy=obsv_u(2,iobsv)
            zz=obsv_u(3,iobsv)
c            r=sqrt(xx*xx+yy*yy+zz*zz)
            h2=(zz**2+yy**2)/xx**2
            if (h2.lt.0.01) then
c              r=xx*(1.0d0+h2/2.0d0-h2**2/8.0d0)
              r=xx*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &      -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
            else
              r=xx*(1.0d0+sqrt(1.0d0+h2))
            endif
            xn=xx/r
            yn=yy/r
            zn=zz/r
            obsv_u(1,iobsv)=xn*xx
            obsv_u(2,iobsv)=yn*xx
            obsv_u(3,iobsv)=zn*xx
          endif
        enddo
      enddo

      nepho_u=max(1,nepho_u)

      if (nepho_u.gt.1) then
        df=(ephmax-ephmin)/(nepho_u-1)
        do i=1,nepho_u
          epho_u(i)=ephmin+(i-1)*df
        enddo
      else
        epho_u(1)=(ephmin+ephmax)/2.0d0
      endif

      xbeta_u=xbeta
      betah_u=betah
      alphah_u=alphah
      betav_u=betav
      alphav_u=alphav
      emith_u=emith
      emitv_u=emitv
      disph_u=disph
      dispph_u=dispph
      dispv_u=dispv
      dispph_u=disppv
      espread_u=espread

      modeph=modeph_u
      pherror_u=pherror
      phgshift_u=phgshift

c      call urad_spline(modewave)
c      stop
c      if (modewave.eq.2) then
c        call urad_nnb(modewave)
c      else if (modewave.eq.3) then
c        call urad_spline(modewave)
c      else
      call urad_amprep(modewave)
c      endif

      stokes_u=stokes_u/1.0d6 ! photons/mm**2

      if (ihbunch.ne.0) then
        fbunch_u(4:14,:)=fbunch_u(4:14,:)*1000.0d0 ! mm
        fbunch_u(17:19,:)=fbunch_u(17:19,:)*1000.0d0 ! mm
        fbunch_u(22:26,:)=fbunch_u(22:26,:)/1.0d6 ! 1/mm**2
      endif

      icbrill=nobsv_u/2+1
      if (abs(phgshift).eq.9999.0d0) then
        do ifrq=1,nepho
          iobfr=icbrill+nobsv_u*(ifrq-1)
          rea(1:2)=(0.0d0,0.0d0)
          rea(3)=arad_u(3,iobfr)
          expsh=rea(3)/abs(rea(3))*1.0d3
          if (phgshift.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
          DO iobs=1,nobsv_u
            iobfr=iobs+nobsv_u*(ifrq-1)
            arad_u(1:6,iobfr)=arad_u(1:6,iobfr)/expsh
          enddo
        enddo
      else if (phgshift.ne.0.0d0) then
        expsh=cdexp(dcmplx(0.0d0,phgshift))*1.0d3
        arad_u=arad_u/expsh
      else
        arad_u=arad_u/1.0d3
      endif

      pincen_u=pincen_u*1000.0d0
      pinw_u=pinw_u*1000.0d0
      pinh_u=pinh_u*1000.0d0
      obsv_u=obsv_u*1000.0d0

c      if (modewave.ne.0) call util_zeit_kommentar(6,'Leaving urad_phase')
      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase',0)

      ical=1
      end
*CMZ :          26/09/2025  11.42.49  by  Michael Scheer
*CMZ :  4.02/00 27/08/2025  14.45.47  by  Michael Scheer
*CMZ :  4.01/07 18/10/2024  09.41.32  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  07.41.13  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  13.39.24  by  Michael Scheer
*CMZ :  4.01/02 14/05/2023  11.47.49  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.34.04  by  Michael Scheer
*CMZ :  4.00/17 05/12/2022  10.30.41  by  Michael Scheer
*CMZ :  4.00/16 17/09/2022  15.46.32  by  Michael Scheer
*CMZ :  4.00/15 02/06/2022  09.45.10  by  Michael Scheer
*CMZ :  4.00/11 28/06/2021  10.33.06  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_amprep(modewave)

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEEP,track.
      include 'track.cmn'
*KEND.
cc+seq,uservar.

      complex*16 :: cde,czero=(0.0d0,0.0d0),ci=(0.0d0,1.0d0),cph00

      double precision :: h2,ddist,wlen,dphi,phase0,cjvsto(4,3)

      double complex , dimension (:,:), allocatable :: aradbuff
      double complex , dimension (:,:,:), allocatable :: arad

      double precision, dimension (:), allocatable :: frq
      double precision, dimension (:,:), allocatable :: wsstokes,pow
      double precision, dimension (:,:,:,:,:), allocatable :: stokesprop
      double precision, dimension (:,:,:), allocatable :: fbunch,stokes

      complex*16, dimension (:), allocatable :: expphiran
      real, dimension (:), allocatable :: pherr,pherrc,phiran
      real, dimension(:,:), allocatable :: pranall,eall

      real eran(6),pran(3),rr(2)

      double complex :: apol,amp0(6),damp(6),amp(6),zexp,
     &  apolh,apolr,apoll,apol45,stokesv(4,3),cero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double precision :: t,udgamtot,upow,vf0,vn,vx0,vx2,vxf0,vxi,vy0,vy2,vyf0,
     &  vyi,vz0,vz2,vzf0,vzi,wlen1,x0,x2,xf0,xi,xlell,y0,y2,yf0,yi,ypi,yy,yyp,
     &  z0,z2,zf0,zi,zpi,zz,zzp,fillb(41),stok1,stok2,stok3,stok4,speknor,
     &  sqnbunch,sqnphsp,specnor,specnor_si,sbnor,rpin,r00(3),xph0,
     &  r(3),r0(3),pw,ph,phsum,pkerr,ppin,parke,pc(3),pcbrill(3),om1,
     &  park,pr,hbarev,obs(3),om,fhigh,flow,gamma,eix,eiy,eiz,emassg,
     &  efx,efy,efz,eharm1,ecdipev,ebeam,dtpho,dt,dtelec,dd0,debeam,
     &  drn0(3),drn00(3),ds,dr0(3),dr00(3),drn(3),dpp,dph,dist,dist0,dobs(3),
     &  bunnor,clight,bunchx,beta,beff,spow,
     &  zp0,yp0,rph,anor,fsum,smax,zob,yob,
     &  xkellip,zampell,yampell,parkv,parkh,zpampell,ypampell,emom,dzpin,dypin,zmin,ymin,phgsh

      double precision xprop,yprop(npinyprop_u),zprop(npinzprop_u),dy,dz,pinwprop,pinhprop
      double complex, dimension(:,:,:,:,:), allocatable :: fprop
      !double complex, dimension(:,:,:), allocatable :: fpriv
      double complex :: fpriv(3,npinzprop_u,npinyprop_u)

      double complex, dimension (:), allocatable ::
     &  uampex,uampey,uampez,uampbx,uampby,uampbz
      double precision, dimension (:,:), allocatable :: utraxyz,ustokes

      double complex :: rea(3),expsh

      integer :: kfreq,iobsv,i,np2,nelec,mbunch,meinbunch,ibu,jbun,
     &  kran=6,icbrill,ilo,kobsv,i1,i2,n,
     &  ifail,ndimu,nstepu,ith,noespread,noemit,jbunch,jubunch,jhbunch,
     &  jcharge=-1,lmodeph,nclo,jeneloss=0,iamppin,
     &  iamppincirc=0,ifrob,iobfr,isub,jvelofield=0,nlbu=0,nepho,ielo,
     &  modewave,iepho,ifieldprop,nzprop,nyprop,im,izm,iym,ifix

      integer, dimension (:), allocatable :: lnbunch

      integer :: idebug=0, lbunch=0, ierr=0, ielec=0
      integer ibunch,ihbunch,mthreads,nobsv,nobsvo,iemit,noranone,iz,iy,ipz,ipy,nobsvz,nobsvy
      integer iobm,iobp,iobfrm,iobfrp
      integer :: ical=0

      save ical
c      integer iuser
c      iuser=user(3)

      nelec_u=max(1,nelec_u)
      mthreads_u=max(1,mthreads_u)
      nelec=nelec_u

      nobsvy=npiny_u
      nobsvz=npinz_u
      nobsvo=npinzo_u*npinyo_u
      dzpin=pinw_u/max(1,npinzo_u-1)
      dypin=pinh_u/max(1,npinyo_u-1)
      zmin=pincen_u(3)-pinw_u/2.0d0
      ymin=pincen_u(2)-pinh_u/2.0d0

      if (nelec_u.gt.1) then
        if (nelec_u.lt.mthreads_u.and.nelec_u.gt.1) then
          mthreads_u=nelec_u
        else
          nelec_u=max(mthreads_u,nelec_u/mthreads_u*mthreads_u)
        endif
      endif

      noranone=noranone_u

      if (nelec_u.eq.1.and.noranone.ne.0) then
        ibunch=0
      else
        ibunch=1
      endif

      ihbunch=ihbunch_u
      mthreads=mthreads_u

      if (modepin_u.eq.1) then
        iamppin=3
        nobsv=1
      else
        iamppin=1
        nobsv=npiny_u*npinz_u
      endif

      icbrill=nobsv/2+1

c      jhbunch=max(0,ihbunch)
      jhbunch=ihbunch
      meinbunch=nelec_u

      lbunch=0

c      if (jhbunch.ne.0) then
c        jhbunch=max(1,jhbunch)
c      endif

      nepho=nepho_u
      if (ifieldprop_u.eq.2) then
        allocate(
c     &    fpriv(3,npinzprop_u,npinyprop_u),
     &    fprop(3,npinzprop_u,npinyprop_u,nepho_u,mthreads),
     &    stokesprop(4,npinzprop_u,npinyprop_u,nepho_u,mthreads))
        stokesprop=0.0d0
      endif

      if (modepin_u.ne.0) then
        if (ical.eq.0) allocate(fieldbunch(7,npinzo_u,npinyo_u,nepho_u),stat=ierr)
        if (ierr.ne.0) then
          print*,""
          print*,"*** Warning in urad_amprep: Could not allocate buffer for beam Ntuple ***"
          print*,""
          return
        endif
        fieldbunch=czero
      endif

      call urad_field_ini(perlen_u,shift_u,beffv_u,beffh_u,modewave)

      if (perlen_u.ne.0.0d0) then
        emom=emasse1*dsqrt((gamma_u-1.0d0)*(gamma_u+1.0d0))
        xkellip=twopi1/perlen_u
        parkh=echarge1*dabs(beffh_u)*perlen_u/(twopi1*emasskg1*clight1)
        parkv=echarge1*dabs(beffv_u)*perlen_u/(twopi1*emasskg1*clight1)

        if (modewave.eq.0) then
c*** OBSOLITE, SEE z0= further down
          zampell=beffv_u*clight1/emom/xkellip**2
          yampell=beffh_u*clight1/emom/xkellip**2
c        zampell=zmx
c        yampell=ymx
c        print*,zpampell
c        ypampell=parkh/gamma_u
c        zpampell=tan(phimx)
c        print*,zpampell
c        stop
        else
          call util_break
          yampell=ymx-ymn
          zampell=zmx-zmn
        endif
      else
        print*,''
        print*,'*** Error in urad_amprep: Zero period-length of undulator ***'
        print*,''
        stop
      endif

      dr00=[1.0d0,0.0d0,0.0d0]
      drn00=dr00/norm2(dr00)
      dr00=drn00*perlen_u
      r00=[0.0d0,0.0d0,0.0d0]

      x0=-perlen_u/2.0d0
      y0=0.0d0
      z0=0.0d0

      beta=dsqrt((1.0d0-1.0d0/gamma_u)*(1.0d0+1.0d0/gamma_u))

      clight=clight1
      hbarev=hbarev1
      ecdipev=ecdipev1
      emassg=emassg1

      if (modewave.eq.0) then
        z0=-zampell*cos(shift_u/2.0d0/perlen_u*twopi1)
        zp0=zpampell*sin(shift_u/2.0d0/perlen_u*twopi1)
        y0=-yampell*cos(shift_u/2.0d0/perlen_u*twopi1)
        yp0=-ypampell*sin(shift_u/2.0d0/perlen_u*twopi1)
      else
        y0=ytrack
        z0=ztrack
        zp0=vztrack/vxtrack
        yp0=vytrack/vxtrack
      endif

      xf0=-x0
      yf0=y0
      zf0=z0

      vn=clight*beta

      vx0=vn/sqrt(1.0d0+(zp0**2+yp0**2))
      vy0=vn*yp0
      vz0=vn*zp0

      vxf0=vx0
      vyf0=vy0
      vzf0=vz0

      vxi=vx0
      vyi=vy0
      vzi=vz0

      r0=r00
      dr0=dr00
      drn0=drn00
      r=r0
      drn=drn0

      vf0=norm2([vxf0,vyf0,vzf0])
      efx=vxf0/vf0
      efy=vyf0/vf0
      efz=vzf0/vf0

      nclo=nint(perlen_u/step_u)+1

      ds=step_u
c      dtim0=ds/beta

      ndimu=nint(nclo*1.1)

      r0=[x0,y0,z0]
      dr0=[xf0-x0,yf0-y0,zf0-z0]
      dr0=[efx,efy,efz]*perlen_u
      r0=r0+dr0/2.0d0

      allocate(frq(nepho_u),
     &  uampex(nepho_u),uampey(nepho_u),uampez(nepho_u),
     &  uampbx(nepho_u),uampby(nepho_u),uampbz(nepho_u),pow(nobsv,mthreads),
     &  utraxyz(14,ndimu),ustokes(4,nepho_u))

      pow=0.0d0
      frq=epho_u

      flow=frq(1)
      fhigh=frq(nepho_u)

      beff=sqrt(beffv_u**2+beffh_u**2)
      park=echarge1*beff*perlen_u/(2.*pi1*emasskg1*clight)
      wlen1=(1+park**2/2.)/2./gamma_u**2*perlen_u*1.0d9

      if (wlen1.ne.0.0) then
        eharm1=wtoe1/wlen1
      else
        eharm1=0.0d0
      endif

      dtpho=perlen_u/clight

      allocate(pherrc(nper_u),pherr(nper_u),arad(6,nepho_u*nobsv,mthreads),
     &  expphiran(max(1,nelec_u)))

      allocate(pranall(2,nelec_u))
      do i=1,nelec_u
        call util_random(3,pran)
        pranall(:,i)=pran(1:2)
        expphiran(i)=exp(dcmplx(0.0d0,twopi1*pran(3)))
      enddo

      if (ibunch.eq.0.or.
     &    emith_u.eq.0.0d0.and.emitv_u.eq.0.0d0.and.espread_u.eq.0.0d0) then
        iemit=0
        ibunch=0
        nelec=1
        nelec_u=1
      else
        iemit=1
      endif

      if (nelec_u.ne.nelec) then
        print*,''
        print*,'--- Warning in urad_amprep: Nelec adjusted to multiple of number of threads:',nelec_u
        print*,''
      endif

      if (iemit.ne.0) then
        allocate(eall(6,nelec_u))
        do i=1,nelec_u
          xi=x0
          if (modepin_u.ne.2) then
            call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &        emith_u,emitv_u,
     &        disph_u,dispph_u,dispv_u,disppv_u,
     &        espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
          else
            ! espread only for folding procedure
            call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &        0.0d0,0.0d0,
     &        disph_u,dispph_u,dispv_u,disppv_u,
     &        espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
          endif
          eall(1,i)=xi-x0
          eall(2,i)=yi
          eall(3,i)=zi
          eall(4,i)=ypi
          eall(5,i)=zpi
          eall(6,i)=dpp
        enddo
        if (noranone.ne.0) eall(:,1)=0.0
      endif

      !allocate(affe(6,nepho_u*nobsv))

      allocate(wsstokes(4,nepho_u*nobsv),stokes(4,nepho_u*nobsv,mthreads))
      stokes=0.0d0
      arad=(0.0d0,0.0d0)

      np2=nper_u/2

      call util_random_gauss_omp(nper_u,pherr,rr)
      pherrc=pherr

      lmodeph=modeph_u

      if (pherror_u.ne.0.0d0.and.(lmodeph.lt.0.or.lmodeph.gt.2)) then
        write(6,*) ""
        write(6,*) "*** Error in urad_amprep: MODEPH must be 0,1, or 2 ***"
        write(6,*) "*** Program aborted ***"
      endif

      if (lmodeph.eq.0.and.eharm1.ne.0.0d0) then
        om1=eharm1/hbarev
        pherr=sngl(pherrc*pherror_u/360.0d0*twopi1/om1)
      else if (lmodeph.eq.1) then
        pherr=sngl(pherr*pherror_u)
      else if (lmodeph.eq.2) then
        pherr(nper_u)=0.0
        phsum=0.0d0
        do i=1,nper_u-1
          pherr(i)=pherr(i)+pherrc(i)
          pherr(i+1)=pherr(i+1)-pherrc(i)
          phsum=phsum+pherr(i)
        enddo
        phsum=phsum+pherr(nper_u)
      else
        pherr=0.0d0
      endif !(lmodeph.eq.0)

      mbunch=max(1,nelec_u)
      nelec=nelec_u

      if (ibunch.ne.0.and.bunchcharge_u.ne.0.0d0) then
        sqnbunch=mbunch
        sqnphsp=sqrt(bunchcharge_u/echarge1)
     &    *meinbunch
     &    /(bunchcharge_u/echarge1)
        bunnor=1.0d0/mbunch
      else
        sqnbunch=mbunch
        sqnphsp=sqrt(dble(nelec_u))
        bunnor=1.0d0/mbunch
      endif

      beff=sqrt(beffv_u**2+beffh_u**2)
      parke=echarge1*beff*perlen_u/(2.*pi1*emasskg1*clight)
      xlell=perlen_u

      ielec=0

      pow=0.0d0
      noemit=0
      noespread=0
      jbunch=ibunch
      jubunch=0
      ebeam=ebeam_u
      debeam=espread_u
      stokesv=vstokes
      specnor=
     &  banwid_u
     &  /(4.0d0*pi1**2*clight*hbarev)
     &  /(4.0d0*pi1*eps01)
     &  *curr_u

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr_u ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid_u !BW
     &  /1.0d6 !m**2 > mm**2

      sbnor=specnor*bunnor
      speknor=specnor

      jeneloss=0
      pw=pinw_u
      ph=pinh_u
      !pr=pinr
      pc=pincen_u
      do iobsv=1,nobsv
        if (abs(obsv_u(2,iobsv)).lt.1.0d-9) obsv_u(2,iobsv)=0.0d0
        if (abs(obsv_u(3,iobsv)).lt.1.0d-9) obsv_u(3,iobsv)=0.0d0
      enddo
      pcbrill=obsv_u(:,icbrill)
      phgsh=phgshift_u

      rea=(0.0d0,0.0d0)
      expsh=(1.0d0,0.0d0)

      if (ifieldprop_u.eq.2) then
        if (npinzprop_u.eq.1) then
          zprop(1)=0.0d0
        else
          dz=pinwprop_u/(npinzprop_u-1)/1000.0d0
          zprop(1)=-pinwprop_u/2.0d0/1000.0d0
          do i=2,npinzprop_u
            zprop(i)=zprop(i-1)+dz
          enddo
        endif

        xprop=pinxprop_u

        if (npinyprop_u.eq.1) then
          yprop(1)=0.0d0
        else
          dy=pinhprop_u/(npinyprop_u-1)/1000.0d0
          yprop(1)=-pinhprop_u/2.0d0/1000.0d0
          do i=2,npinyprop_u
            yprop(i)=yprop(i-1)+dy
          enddo
        endif
      endif !ifieldprop

      ifieldprop=ifieldprop_u
      nyprop=npinyprop_u
      nzprop=npinzprop_u
      cjvsto=dconjg(vstokes)

      pinwprop=pinwprop_u
      pinhprop=pinhprop_u

      ifix=ifixphase_u

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(nepho,nobsvz,nobsvy,nobsv,nelec,frq,nper_u,np2,perlen_u,clight,hbarev,
!$OMP& ifieldprop,xprop,nyprop,nzprop,yprop,zprop,cjvsto,fpriv,pinwprop,pinhprop,
!$OMP& flow,fhigh,czero,cone,rea,expsh,ifix,zob,yob,
!$OMP& x0,y0,z0,xf0,yf0,zf0,vx0,vy0,vz0,vxf0,vyf0,vzf0,gamma_u,sbnor,speknor,
!$OMP& efx,efy,efz,ds,ndimu,curr_u,xlell,parke,amp,amp0,
!$OMP& uampex,uampey,uampez,uampbx,uampby,uampbz,
!$OMP& lmodeph,zp0,yp0,modewave,
!$OMP& jbunch,jubunch,jhbunch,noespread,noemit,ebeam,
!$OMP& stokesv,icbrill,obsv_u,emassg,debeam,dispv_u,disppv_u,
!$OMP& betah_u,alphah_u,betav_u,alphav_u,emith_u,emitv_u,disph_u,dispph_u,
!$OMP& pran,pranall,eall,fillb,r0,dr0,iamppin,iamppincirc,pc,phase0,pr,banwid_u,
!$OMP& pw,ph,idebug,pcbrill,wsstokes,vn,bunchlen_u,modebunch_u,icohere_u)
!$OMP& SHARED(mthreads,stokes,pherr,expphiran,lbunch,lnbunch,modepin_u,fieldbunch,npinzo_u,nobsvo,dzpin,dypin,
!$OMP& fbunch_u,jcharge,jeneloss,jvelofield,iemit,noranone,arad,pow,zmin,ymin,phgsh,fprop,stokesprop)

      jbun=1
      isub=0
      iobsv=0
      ielo=0
      xph0=-perlen_u*dble(nper_u)/2.0d0

!$OMP DO

c      do ilo=1,nelec*nobsv
      do ielec=1,nelec
      do iobsv=1,nobsv

        wsstokes=0.0d0

        !affe=(0.0D0,0.0D0)
        spow=0.0d0

        ith=OMP_GET_THREAD_NUM()+1

c        iobsv=mod(ilo-1,nobsv)+1
c        ibu=(ilo-1)/nobsv+1
        ibu=ielec
        jbun=ibu

        iy=(iobsv-1)/nobsvz+1
        iz=mod(iobsv-1,nobsvz)+1

        !if (iz.gt.nobsvz/2+1) call til_break

c        ielec=ibu

        xi=x0
        yi=y0
        zi=z0

        zpi=vz0/vx0
        ypi=vy0/vx0

        x2=xf0
        y2=yf0
        z2=zf0

        vx2=vxf0
        vy2=vyf0
        vz2=vzf0

        gamma=gamma_u

        dpp=0.0d0

        if (iemit.ne.0) then

          if (noranone.eq.0.or.ielec.ne.1) then

            bunchx=eall(1,ielec)

            xi=xi+bunchx
            yy=eall(2,ielec)
            zz=eall(3,ielec)

            yyp=eall(4,ielec)
            zzp=eall(5,ielec)

            dpp=eall(6,ielec)
            gamma=(1.0d0+dpp)*gamma_u

            ! assume beta(s)=beta0(s)+s**2/beta(0) and alpha0=-s/beta(0)
            ! and a drift transfer-matrix ((1,s),(1,0))

            zi=zz-x0*zzp !inverse transformation
            zpi=zzp

            yi=yy-x0*yyp
            ypi=yyp

            ! simple treatment of closed orbit, assume small angles

            zi=zi+z0
            zpi=zpi+zp0

            yi=yi+y0
            ypi=ypi+yp0

          else

            xi=x0
            yi=y0
            zi=z0
            ypi=yp0
            zpi=zp0

            bunchx=0.0d0

          endif

        else

          xi=x0
          yi=y0
          zi=z0
          ypi=yp0
          zpi=zp0

          bunchx=0.0d0

        endif !iemit

c+self,if=old.
c        zi=zi+dpp*di0
c        zpi=zpi+dpp*dd0
c+self.
        vn=clight*dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))

        vxi=vn/sqrt(1.0d0+ypi**2+zpi**2)
        vyi=vxi*ypi
        vzi=vxi*zpi

        obs=obsv_u(1:3,iobsv)

        if (noranone.eq.0.or.ielec.ne.1.or.iobsv.ne.icbrill) then
          if (iamppin.eq.3) then
            !call util_random(2,pran)
            pran(1:2)=pranall(:,ielec)
            if (iamppincirc.eq.0) then
              obs(2)=pc(2)+(pran(1)-0.5)*ph
              obs(3)=pc(3)+(pran(2)-0.5)*pw
            else
              rpin=(pran(1)-0.5)*pr
              ppin=pran(2)*twopi1
              obs(2)=pc(2)+rpin*cos(ppin)
              obs(3)=pc(3)+rpin*sin(ppin)
            endif
          endif
        endif

        vn=norm2([vxi,vyi,vzi])
        eix=vxi/vn
        eiy=vyi/vn
        eiz=vzi/vn

        h2=((obs(2)-yi)**2+(obs(3)-zi)**2)/(obs(1)-xph0)**2
        if (h2.lt.0.01) then
          rph=abs(obs(1)-xph0)*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &      -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
        else
          rph=sqrt((obs(1)-xph0)**2+((obs(2)-yi)**2+(obs(3)-zi)**2))
        endif

        phase0=(rph-(obsv_u(1,icbrill)-xph0))/clight

        call urad_e_b_field(
     &    jcharge,curr_u,
     &    gamma,udgamtot,
     &    xi,yi,zi,vxi,vyi,vzi,
     &    xf0,yf0,zf0,efx,efy,efz,
     &    x2,y2,z2,vx2,vy2,vz2,dtelec,ds,
     &      0,nstepu,ndimu,utraxyz,phase0,
     &    obs(1),obs(2),obs(3),flow,fhigh,
     &    nepho,frq,uampex,uampey,uampez,uampbx,uampby,uampbz,
     &    ustokes,upow,
     &    jeneloss,jvelofield,ifail,ith,banwid_u,modewave
     &    )

        r0=[xi,yi,zi]
        dr0=[x2-xi,y2-yi,z2-zi]

        drn=dr0/norm2(dr0)
        r0=r0+dr0/2.0d0

        do kfreq=1,nepho

          iobfr=iobsv+nobsv*(kfreq-1)

          om=frq(kfreq)/hbarev

c          if (modewave.eq.0) then
c            amp0=[
c     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
c     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
c     &        ]*1.0d3/sqrt(speknor/curr_u*0.10d0) !urad
c          else
            amp0=[
     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
     &        ]*1.0d3/sqrt(speknor) !urad
c          endif

c          call util_random(1,pran)
c          amp0=amp0*dcmplx(0.0d0,dble(pran(1)*twopi1))

          amp=(0.0d0,0.0d0)
          t=bunchx/vn

          do i=1,nper_u

            r=r0+(i-np2-1)*dr0
            dobs=obs-r
            dist0=norm2(obs-r0)
            dist=norm2(dobs)

            if (kfreq.eq.1) then
              spow=spow+upow*(dist0/dist)**2
              pow(iobsv,ith)=pow(iobsv,ith)+upow*(dist0/dist)**2
            endif

            if (lmodeph.eq.0) then
!!!!!                dt=xlell/clight*((1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
!!!!!     &            (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)
              h2=
     &          (ypi-yp0-dobs(2)/dobs(1))**2 +
     &          (zpi-zp0-dobs(3)/dobs(1))**2
c26.4.2024     &          ((ypi-yp0-dobs(2))/dobs(1))**2 +
c26.4.2024     &          ((zpi-zp0-dobs(3))/dobs(1))**2

c              dt=xlell/clight*
c     &          (
c     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+h2/2.0d0-h2**2/8.0d0
c     &          )

              dph=om*(t+pherr(i))

              dt=xlell/clight*
     &          (
     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
     &          (((((-0.0205078125D0*h2+0.02734375D0)*h2
     &          -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2
     &          )

              t=t+dt
            else if (lmodeph.eq.1.or.lmodeph.eq.2) then
              dph=om*t
              pkerr=parke*(1.0d0+pherr(i))
!!!!!                dt=xlell/clight*((1.0d0+pkerr**2/2.0d0)/2.0d0/gamma**2+
!!!!!     &            (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)
              h2=
     &          (ypi-yp0-dobs(2)/dobs(1))**2 +
     &          (zpi-zp0-dobs(3)/dobs(1))**2
c26.4.2024              h2=((ypi-yp0-dobs(2))**2+(zpi-zp0-dobs(3))**2)/dobs(1)**2
              dt=xlell/clight*
     &          (
c25.4.2024     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
     &          (1.0d0+pkerr**2/2.0d0)/2.0d0/gamma**2+
     &          (((((-0.0205078125D0*h2+0.02734375D0)*h2
     &          -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2
     &          )
              t=t+dt
            endif !lmodeph

            zexp=cdexp(dcmplx(0.0d0,dph))
            damp=amp0*zexp*dist0/dist
            amp=amp+damp

            if (jhbunch.ne.0) then

              if (
     &            ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
     &            mod(ielec,jhbunch).eq.0) .or.
     &            (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

                if (i.eq.1) then
                  fillb(5)=r(1)
                  fillb(6)=r(2)
                  fillb(7)=r(3)
                  fillb(8)=ypi
                  fillb(9)=zpi
                else if (i.eq.nper_u) then

                  if (abs(phgsh).eq.9999.0d0) then
                    rea=amp(1:3)
                    expsh=rea(3)/abs(rea(3))
                    if (phgsh.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
                    amp=amp/expsh
                  else if (phgsh.ne.0.0d0) then
                    expsh=cdexp(dcmplx(0.0d0,phgsh)) !*1.0d3
                    amp=amp/expsh
                  endif

                  if (phgsh.eq.9999.0d0) then
                    rea=amp(1:3)
                    expsh=rea(3)/abs(rea(3))
                    amp=amp/expsh
                  else if (phgsh.eq.-9999.0d0) then
                    rea=amp(1:3)
                    expsh=rea(3)/abs(rea(3))*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
                    amp=amp/expsh
                  else if (phgsh.ne.0.0d0) then
                    expsh=cdexp(dcmplx(0.0d0,phgsh)) !*1.0d3
                    amp=amp/expsh
                  endif

                  fillb(10:12)=r
                  fillb(13)=ypi
                  fillb(14)=zpi
                  fillb(30)=dreal(amp(1))
                  fillb(31)=dimag(amp(1))
                  fillb(32)=dreal(amp(2))
                  fillb(33)=dimag(amp(2))
                  fillb(34)=dreal(amp(3))
                  fillb(35)=dimag(amp(3))
                  fillb(36)=dreal(amp(4))
                  fillb(37)=dimag(amp(4))
                  fillb(38)=dreal(amp(5))
                  fillb(39)=dimag(amp(5))
                  fillb(40)=dreal(amp(6))
                  fillb(41)=dimag(amp(6))
                endif

              endif

            endif

          enddo !nper_u

          if (ifix.ne.0) then
            amp=amp*expphiran(ielec)
          endif

          if (modepin_u.ne.0) then
            iy=int((obs(2)-ymin)/dypin)+1
            iz=int((obs(3)-zmin)/dzpin)+1
            !print*,ilo,ith,obs(3),zmin,dzpin,iz
            fieldbunch(1:6,iz,iy,kfreq)=fieldbunch(1:6,iz,iy,kfreq)+amp(1:6)
            fieldbunch(7,iz,iy,kfreq)=fieldbunch(7,iz,iy,kfreq)+cone
          endif

          apolh=
     &      amp(1)*conjg(stokesv(1,1))
     &      +amp(2)*conjg(stokesv(1,2))
     &      +amp(3)*conjg(stokesv(1,3))

          apolr=
     &      amp(1)*conjg(stokesv(2,1))
     &      +amp(2)*conjg(stokesv(2,2))
     &      +amp(3)*conjg(stokesv(2,3))

          apoll=
     &      amp(1)*conjg(stokesv(3,1))
     &      +amp(2)*conjg(stokesv(3,2))
     &      +amp(3)*conjg(stokesv(3,3))

          apol45=
     &      amp(1)*conjg(stokesv(4,1))
     &      +amp(2)*conjg(stokesv(4,2))
     &      +amp(3)*conjg(stokesv(4,3))

          stok1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
          stok2=dreal(-stok1+2.0d0*apolh*conjg(apolh))
          stok3=dreal(2.0d0*apol45*conjg(apol45)-stok1)
          stok4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

          wsstokes(1,iobfr)=wsstokes(1,iobfr)+stok1*sbnor
          wsstokes(2,iobfr)=wsstokes(2,iobfr)+stok2*sbnor
          wsstokes(3,iobfr)=wsstokes(3,iobfr)+stok3*sbnor
          wsstokes(4,iobfr)=wsstokes(4,iobfr)+stok4*sbnor

          stokes(1:4,iobfr,ith)=stokes(1:4,iobfr,ith)+wsstokes(1:4,iobfr)

          !affe(:,iobfr)=affe(:,iobfr)+amp
          !arad(:,iobfr,ith)=arad(:,iobfr,ith)+affe(:,iobfr)
          arad(:,iobfr,ith)=arad(:,iobfr,ith)+amp

c          if (
c     &        ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
c     &        mod(ielec,jhbunch).eq.0) .or.
c     &        (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then
c            if (ielec.gt.4) then
c              print*,jhbunch,ith,ilo,ielec
c            endif
c          endif

          if (jhbunch.ne.0) then

            if (
     &          ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
     &          mod(ielec,jhbunch).eq.0) .or.
     &          (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

c              print*,jhbunch,ith,ilo,jbun,isub,ibu
              fillb(1)=jbun
              fillb(2)=isub
              fillb(3)=ibu
              fillb(4)=bunchx
              fillb(15)=gamma*emassg
              fillb(16)=udgamtot*emassg
              fillb(17)=obs(1)
              fillb(18)=obs(2)
              fillb(19)=obs(3)
              fillb(20)=kfreq
              fillb(21)=frq(kfreq)

              fillb(22)=wsstokes(1,iobfr)*nelec

              fillb(23)=wsstokes(1,iobfr)*nelec
              fillb(24)=wsstokes(2,iobfr)*nelec
              fillb(25)=wsstokes(3,iobfr)*nelec
              fillb(26)=wsstokes(4,iobfr)*nelec

              fillb(27)=spow
              fillb(28)=1
              fillb(29)=dtelec

              fillb(30)=dreal(amp(1))
              fillb(31)=dimag(amp(1))
              fillb(32)=dreal(amp(2))
              fillb(33)=dimag(amp(2))
              fillb(34)=dreal(amp(3))
              fillb(35)=dimag(amp(3))
              fillb(36)=dreal(amp(4))
              fillb(37)=dimag(amp(4))
              fillb(38)=dreal(amp(5))
              fillb(39)=dimag(amp(5))
              fillb(40)=dreal(amp(6))
              fillb(41)=dimag(amp(6))
              lbunch=lbunch+1
              fbunch_u(:,lbunch)=fillb(:)
            endif !fill

          endif !jhbunch

          if (ifieldprop.eq.2) then
            if (iobsv.eq.1) then
              fprop(1:3,1:nzprop,1:nyprop,kfreq,ith)=(0.0d0,0.0d0)
            endif
            call urad_phase_prop_point(obs,amp(1:3),nzprop,nyprop,
     &        xprop,yprop,zprop,pinwprop,pinhprop,frq(kfreq),fpriv)
            fprop(:,:,:,kfreq,ith)=fprop(:,:,:,kfreq,ith)+fpriv(:,:,:)
            if (iobsv.eq.nobsv) then
              i=0
              do ipy=1,nyprop
                do ipz=1,nzprop
                  i=i+1+nzprop*nyprop*(kfreq-1)
                  apolh=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(1,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(1,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(1,3)

                  apolr=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(2,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(2,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(2,3)

                  apoll=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(3,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(3,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(3,3)

                  apol45=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(4,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(4,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(4,3)

                  stok1=dreal(
     &              apolr*conjg(apolr)+
     &              apoll*conjg(apoll))

                  stok2=-stok1+
     &              dreal(2.*apolh*conjg(apolh))

                  stok3=
     &              dreal(2.*apol45*conjg(apol45))-
     &              stok1

                  stok4=dreal(
     &              apolr*conjg(apolr)-
     &              apoll*conjg(apoll))

                  stokesprop(1,ipz,ipy,kfreq,ith)=stokesprop(1,ipz,ipy,kfreq,ith)+stok1
                  stokesprop(2,ipz,ipy,kfreq,ith)=stokesprop(2,ipz,ipy,kfreq,ith)+stok2
                  stokesprop(3,ipz,ipy,kfreq,ith)=stokesprop(3,ipz,ipy,kfreq,ith)+stok3
                  stokesprop(4,ipz,ipy,kfreq,ith)=stokesprop(4,ipz,ipy,kfreq,ith)+stok4

                enddo
              enddo
            endif
          endif

        enddo !kfreq

      enddo !iobsv
      enddo !nelec
c      enddo !ilo

!$OMP END DO
!$OMP END PARALLEL

      do ith=1,mthreads
        pow_u(:)=pow_u(:)+pow(:,ith)
        arad_u(:,:)=arad_u(:,:)+arad(:,:,ith)
      enddo

      if (globphase_u.eq.9999.0d0) then
        do kfreq=1,nepho
          iobfr=icbrill+nobsv*(kfreq-1)
          cph00=arad_u(3,iobfr)/abs(arad_u(3,iobfr))
          do iobsv=1,nobsv
            iobfr=iobsv+nobsv*(kfreq-1)
            arad_u(:,iobfr)=arad_u(:,iobfr)/cph00
          enddo
        enddo
      else
        arad_u=arad_u*exp(ci*globphase_u)
      endif

      pow_u=pow_u/sqnbunch

      if (ifieldprop_u.eq.2) then
        smax=0.0d0
        do ith=1,mthreads
          do kfreq=1,nepho_u
            i=0
            do iy=1,nyprop
              do iz=1,nzprop
                i=i+1+nobsvprop_u*(kfreq-1)
                stokesprop_u(1:4,i)=stokesprop_u(1:4,i)+stokesprop(1:4,iz,iy,kfreq,ith)
                if(abs(stokesprop_u(1,i)).gt.smax) then
                  smax=abs(stokesprop_u(1,i))
                  izm=iz
                  iym=iy
                  im=i
                endif
              enddo
            enddo
          enddo
        enddo
      endif

      if (icohere_u.eq.0) then

        arad_u=arad_u/sqnbunch

        do ith=1,mthreads
          stokes_u(:,:)=stokes_u(:,:)+stokes(:,:,ith)
        enddo

      else

        do iobsv=1,nobsv
          do kfreq=1,nepho

            iobfr=iobsv+nobsv*(kfreq-1)

            amp(1:3)=arad_u(1:3,iobfr) !/sqnphsp

            apolh=
     &        amp(1)*conjg(stokesv(1,1))
     &        +amp(2)*conjg(stokesv(1,2))
     &        +amp(3)*conjg(stokesv(1,3))

            apolr=
     &        amp(1)*conjg(stokesv(2,1))
     &        +amp(2)*conjg(stokesv(2,2))
     &        +amp(3)*conjg(stokesv(2,3))

            apoll=
     &        amp(1)*conjg(stokesv(3,1))
     &        +amp(2)*conjg(stokesv(3,2))
     &        +amp(3)*conjg(stokesv(3,3))

            apol45=
     &        amp(1)*conjg(stokesv(4,1))
     &        +amp(2)*conjg(stokesv(4,2))
     &        +amp(3)*conjg(stokesv(4,3))

            stok1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
            stok2=dreal(-stok1+2.0d0*apolh*conjg(apolh))
            stok3=dreal(2.0d0*apol45*conjg(apol45)-stok1)
            stok4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

            stokes_u(1,iobfr)=stok1*sbnor
            stokes_u(2,iobfr)=stok2*sbnor
            stokes_u(3,iobfr)=stok3*sbnor
            stokes_u(4,iobfr)=stok4*sbnor

          enddo
        enddo

      endif !icohere_u

c      if (ihbunch.ne.0) then
c        n=0
c        do i=1,nlbu
c          do ith=1,mthreads_u
c            if (fbunch(21,i,ith).ne.0.0d0) then
c              n=n+1
c              fbunch_u(:,n)=fbunch(:,i,ith)
c            endif
c          enddo
c        enddo
c        deallocate(fbunch)
c      endif

      !deallocate(affe)
      deallocate(frq,uampex,uampey,uampez,uampbx,uampby,uampbz,utraxyz,
     &  pherrc,pherr,expphiran,arad,pow,pranall,wsstokes,stokes)

      if (iemit.ne.0) deallocate(eall)

      iobfr=nobsv_u*nepho_u/2+1
      amp(1:3)=arad_u(1:3,iobfr)

      anor=sqrt(stokes_u(1,iobfr)/
     &  (amp(1)*dconjg(amp(1))+amp(2)*dconjg(amp(2))+amp(3)*dconjg(amp(3))))

      arad_u=arad_u*anor/sqrt(specnor_si)

      if (modepin_u.ne.0) then
        do iepho=1,nepho_u
          do iy=1,npinyo_u-1
            do iz=1,npinzo_u-1
              fsum=max(1.0d0,dreal(fieldbunch(7,iz,iy,iepho)))
              fieldbunch(1:6,iz,iy,iepho)=fieldbunch(1:6,iz,iy,iepho)/fsum
c            fieldbunch(1:6,iz,iy,iepho)=fieldbunch(1:6,iz,iy,iepho)
            enddo
          enddo
        enddo
        fieldbunch=fieldbunch*anor/sqrt(specnor_si)
      endif

c      if (ifieldprop_u.eq.2) then
      if (ifieldprop_u.gt.0) then
        stokesprop_u=stokesprop_u*specnor_si
      endif

      ical=1
      return
      end
*CMZ :  4.01/02 09/05/2023  13.11.31  by  Michael Scheer
*CMZ :  4.01/00 10/02/2023  13.52.47  by  Michael Scheer
*-- Author :    Michael Scheer   10/02/2023
      subroutine urad_field_ini(perl,shift,b0v,b0h,modewave)

      implicit none

*KEEP,ampli.
      include 'ampli.cmn'
*KEND.

      double precision shift,perl,b0h,b0v
      integer modewave

      phrshift=shift
      phrperl=perl
      phrb0h=b0h
      phrb0v=b0v

      return
      end
*CMZ :  4.01/07 27/04/2024  09.45.56  by  Michael Scheer
*CMZ :  4.01/05 14/04/2024  07.40.31  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.35.56  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  09.04.01  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  15.28.31  by  Michael Scheer
*CMZ :  4.00/15 28/04/2022  15.32.20  by  Michael Scheer
*CMZ :  4.00/13 16/11/2021  17.32.24  by  Michael Scheer
*CMZ :  4.00/09 15/08/2020  08.51.05  by  Michael Scheer
*CMZ :  3.05/05 10/07/2018  09.19.31  by  Michael Scheer
*CMZ :  3.05/04 05/07/2018  11.10.09  by  Michael Scheer
*CMZ :  3.05/00 27/04/2018  15.22.16  by  Michael Scheer
*CMZ :  3.03/04 13/10/2017  09.16.28  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.32.35  by  Michael Scheer
*CMZ :  3.02/04 13/03/2015  10.38.25  by  Michael Scheer
*CMZ :  2.69/02 02/11/2012  16.40.18  by  Michael Scheer
*CMZ :  2.68/05 04/09/2012  13.30.58  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  11.52.24  by  Michael Scheer
*CMZ :  2.68/03 31/08/2012  09.45.28  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_e_b_field(
     &  icharge,current,
     &  gammai,dgamtot,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efxn,efyn,efzn,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &  nthstep,nstep,ndim,traxyz,phase0,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  nphener,phener,aradex,aradey,aradez,aradbx,aradby,aradbz,stokes,powden,
     &  ieneloss,ivelofield,
     &  istatus,ith,banwid,modewave
     &  )
c123456789123456789_123456789_123456789_123456789_123456789_123456789_12
c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

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

c This subroutine calculates the trajectory and the synchrotron radiation
c of an electron passing a magnetic and electric field. The electric field
c part is very preliminary and not yet fully tested. The FORTRAN random
c generator is used in uradrndm. It should be replaced by a better one
c provided by the user.

c The fields B=(bx,by,bz) and E=(ex,ey,ez) are calculated in the routine
c uradfield_omp(x,y,z,bx,by,bz,ex,ey,ez,gamma,istatus) provided by the user.
c As an example uradbmap.f may be used, which reads a 3D map of the mag. field

c Coordinate system (right handed):
c -----------------------------------

c      x: longitudinal direction
c      y: transversal vertical direction
c      z: transversal horizontal direction

c Units:
c ------

c SI-units: m, Tesla, sec, V etc., but eV for the photon energy
c The flux-density is given in Nph/s/mm**2/0.1%BW, the power-density in W/mm**2

c Input:
c --------

c integer icharge: Particle charge ( +/- 1)
c real*8 gammai: Gamma factor of the e-

c real*8 xelec:  Initial x of e-
c real*8 yelec:  Initial y of e-
c real*8 zelec:  Initial z of e-

c real*8 vxelec:  Initial velocity in x of e-
c real*8 vyelec:  Initial velocity in y of e-
c real*8 vzelec:  Initial velocity in z of e-
c The velocity is internally normalized, so the input norm does not matter

c real*8 xf: x of point in exit plane
c real*8 yf: y of point in exit plane
c real*8 zf: z of point in exit plane

c real*8 [efxn,efyn,efzn]: Normal vector of exit plane

c real*8 vnxex: x component of normal vector of exit plane
c real*8 vnyex: y component of normal vector of exit plane
c real*8 vnzex: z component of normal vector of exit plane

c real*8 ds : step size for tracking

c The tracking stops, when the electron hits the exit plane. The size of the last
c step is corrected such that the plane is hit.

c integer nthstep: If nstep > 0, the trajectory array traxyz is filled,
c                  see below
c integer ndim: Dimension of traxyz, see below

c real*8 phelow: Lowest photon energy / eV
c real*8 phehig: Higest photon energy / eV

c integer nphener: Number of equidistant photon energies

c integer ieneloss:  0: no energy loss due to synchotron radiation
c                    1: continous energy loss due to synchotron radiation
c integer ieneloss: -1: no energy loss due to synchotron radiation with quantum
c                       fluctuations

c integer ivelofield: Contral flag for the calculation of the velocity field:
c                    0: the spectrum includes the velocity field
c                    1: the specrum does not include the velocity field
c                    2: the spectrum includes only the velocity field

c Output:
c -------

c integer: istatus: Status flag:
c  0: no error found
c -1: initial gamma or velocity zero
c -2: dimension ndim of traxyz exceeded
c -3: bad value of ivelofield
c  else: status of uradfield_omp

c real*8 xexit: x of last point of the trajectory
c real*8 yexit: y of last point of the trajectory
c real*8 zexit: z of last point of the trajectory
c real*8 texit: t of last point of the trajectory

c real*8 vnxex: x component of norm. velocity vector of last point
c real*8 vnyex: y component of norm. velocity vector of last point
c real*8 vnzex: z component of norm. velocity vector of last point

c real*8 phener(nphener): Array of equidistant photon energies

c integer nstep: Number of tracking steps done, i.e. used length of traxyz

c real*8 traxyz(1:14,i): Array:
c        traxyz(1,i):  x
c        traxyz(2,i):  y
c        traxyz(3,i):  z
c        traxyz(4,i):  tracking time
c        traxyz(5,i):  x-comp. of norm. velocity vector
c        traxyz(6,i):  y-comp. of norm. velocity vector
c        traxyz(7,i):  z-comp. of norm. velocity vector
c        traxyz(8,i):  x-comp. of mag. field in the center of the step
c        traxyz(9,i):  y-comp. of mag. field in the center of the step
c        traxyz(10,i): z-comp. of mag. field in the center of the step
c        traxyz(11,i): gamma
c        traxyz(12,i): x-comp. of elec. field in the center of the step
c        traxyz(13,i): y-comp. of elec. field in the center of the step
c        traxyz(14,i): z-comp. of elec. field in the center of the step

c real* phase0
c The phase is calculated by phase=phase0+n*dt*dphase,
c where dphase is the phase difference of the nth step dt. Phase0=
c (xobsv-xelec)/clight. The phase factor of the integrand is
c exp(i*omega*phase),where omega referes to the considered photon energy.
c complex*16 aradex(i): x-comp. of amplitude of radiation e-field of phener(i)
c complex*16 aradey(i): y-comp. of amplitude of radiation e-field of phener(i)
c complex*16 aradez(i): z-comp. of amplitude of radiation e-field of phener(i)
c complex*16 aradbx(i): x-comp. of amplitude of radiation b-field of phener(i)
c complex*16 aradby(i): y-comp. of amplitude of radiation b-field of phener(i)
c complex*16 aradbz(i): z-comp. of amplitude of radiation b-field of phener(i)

c real*8 array of Stokes parameters of ith photon energy:
c        stokes(1,i): S0, i.e. total intensity
c        stokes(2,i): S1, Stokes parameter of linear +/- 90 degree polarisation
c        stokes(3,i): S2, Stokes parameter of linear +/- 45 degree polarisation
c        stokes(4,i): S3, Stokes parameter of circular polarisation
c        S0 = sqrt(S1**2+S2**2+S3**2)

c Compilation:
c ------------

c For uradbmap at least F90 is required.
c The line length exceeds 72 characters, please use an appropriate
c compiler option. It is recommended to use compiler options to initialize all
c variables to zero and to treat them as saved''

      implicit none

      complex*16
     &  aradex(nphener),aradey(nphener),aradez(nphener),
     &  aradbx(nphener),aradby(nphener),aradbz(nphener),
     &  ziom,zi,zidom,zone,ziomr1,zicr1,zic,
     &  expom1,expom,dexpomph1,dexpomph,ddexpomph,dexpom,
     &  expomv2,vstokes(4,3),
     &  apolh,apolr,apoll,apol45,dum3,ampex,ampey,ampez,ampbx,ampby,ampbz

      double precision
     &  gammai,dgamtot,dt2,powden,t,phase,phase0,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  phener(nphener),dom2,c,rspn
     &  ,traxyz(14,ndim),stokes(4,nphener),x1,y1,z1,vx1,vy1,
     &  vz1,x2,y2,z2,vx2,vy2,vz2
     &  ,ds,dtim,bshift,x2b,y2b,z2b,bx1,by1,bz1,bx2,by2,bz2
     &  ,dgamsum,gamma,dt
     &  ,vxp,vyp,vzp
     &  ,efx,efy,efz,xf,yf,zf,dist2
     &  ,dtim0,beta,vn,efx2,efy2,efz2,t1,t2,clight,c1,
     &  dgamma,vxsign,bx,by,bz,bpx,bpy,bpz,rarg(5),px,py,pz,
     &  dphase,r,rx,ry,rz
     &  ,dom1,rnbx,rnby,rnbz,dum11,rnr4,br4,b3,rnr2,br2,bet1n,
     &  rnx,rny,rnz,r1,banwid,specnor,pownor,current,
     &  stok1,stok2,stok3,stok4,om,dom,hbarev,echarge,eps0,pi,vsto,dph,
     &  r0,efxn,efyn,efzn

      integer ieneloss,istatus,icharge,nphener,ivelofield,
     &  nthstep,izaehl,nstep,ndim,kstep,lstep,kfreq,isto,ifail,idum,ith,
     &  modewave

      integer :: idebug=0,izlimit

c      integer,save :: ical=0
c+seq,uservar.

      data bshift/0.5d0/
      data clight/2.99792458d8/
      data hbarev/6.58211889D-16/
      data eps0/8.854187817D-12/
      data echarge/1.602176462D-19/
      data pi/3.14159265358979d0/

      data zi/(0.0d0,1.0d0)/
      data zone/(1.0d0,0.0d0)/

      dph=0.0d0
      idum=ith
c      ical=ical+1

      if (nphener.gt.0) phener(1)=phelow
      if (nphener.gt.1) dph=(phehig-phelow)/(nphener-1)

      do kfreq=2,nphener
        phener(kfreq)=phener(kfreq-1)+dph
      enddo

      istatus=0
      ifail=0
      if (icharge.le.0) icharge=-1
      if (icharge.gt.0) icharge=1

      vn=norm2([efxn,efyn,efzn])
      efx=efxn/vn
      efy=efyn/vn
      efz=efzn/vn

      x1=xelec
      y1=yelec
      z1=zelec
      vx1=vxelec
      vy1=vyelec
      vz1=vzelec
      t1=0.0d0

      gamma=gammai
      beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))
      vn=sqrt(vx1*vx1+vy1*vy1+vz1*vz1)
      vx1=vx1/vn*clight*beta
      vy1=vy1/vn*clight*beta
      vz1=vz1/vn*clight*beta
      vn=beta*clight

c vxsign takes care for the direction of flight, since particle must gain
c energy if tracked back

      if (vx1.lt.0) then
        vxsign=-1.0d0
      else
        vxsign=1.0d0
      endif

      dgamsum=0.0d0
      dgamtot=0.0d0
      powden=0.0d0

      aradex=(0.0d0,0.0d0)
      aradey=(0.0d0,0.0d0)
      aradez=(0.0d0,0.0d0)

      aradbx=(0.0d0,0.0d0)
      aradby=(0.0d0,0.0d0)
      aradbz=(0.0d0,0.0d0)

      dtim=ds/vn
      dt=dtim
      dt2=dtim*bshift
      dtim0=dtim

      x2=x1
      y2=y1
      z2=z1
      t2=t1

      vx2=vx1
      vy2=vy1
      vz2=vz1

      x2b=x1+vx1*dt2
      y2b=y1+vy1*dt2
      z2b=z1+vz1*dt2

      call uradfield_omp(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,gamma,istatus,
     &  modewave)
      if (istatus.ne.0) ifail=ifail+abs(istatus)
      istatus=0

      nstep=0
      vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)

      if (gamma.le.0.0d0.or.vn.le.0.0d0) then
        istatus=-1
        return
      endif

      kstep=-1
      nstep=0

      if(nthstep.gt.0) then

        nstep=nstep+1
        if (nstep.gt.ndim) then
          nstep=nstep-1
          istatus=-2
          goto 9000
        endif

        kstep=kstep+1
        if (kstep.eq.nthstep) kstep=0

        traxyz(1,nstep)=x2
        traxyz(2,nstep)=y2
        traxyz(3,nstep)=z2
        traxyz(4,nstep)=t2
        traxyz(5,nstep)=vx2/vn
        traxyz(6,nstep)=vy2/vn
        traxyz(7,nstep)=vz2/vn
        traxyz(8,nstep)=bx2
        traxyz(9,nstep)=by2
        traxyz(10,nstep)=bz2
        traxyz(11,nstep)=gamma
        traxyz(12,nstep)=efx2
        traxyz(13,nstep)=efy2
        traxyz(14,nstep)=efz2

      endif !nthstep.gt.0

      dom=0.0d0
      om=0.0d0
      if (nphener.gt.1) then
        om=phener(1)/hbarev
        dom=(phener(2)-phener(1))/hbarev
      else if (nphener.eq.1) then
        om=phener(1)/hbarev
      endif

      c=clight
      c1=1.0d0/clight

      zidom=zi*dom
      ziom=zi*om
      zic=zi*c

      lstep=0
      t=-dt
      r0=xobsv-xelec
      r=sqrt((xobsv-x1)**2+((yobsv-y1)**2+(zobsv-z1)**2))
c      PHASE=(r-r0)*c1
      phase=phase0
      expom1=zone
      dexpomph1=zone

c--- Loop der Trajektorie

      izaehl=0
      izlimit=(xf-xelec)/ds*2
1000  continue
      !all util_break
      izaehl=izaehl+1
      if (izaehl.gt.izlimit) then
        print*,"*** Warning in urad_e_b_field: Limits of steps exceeded ***"
        istatus=-4
        return
      endif
      if (x2.ne.x2) then
        istatus=-99
        return
      endif

      if (lstep.eq.1) then
        dtim=abs(dist2)/vn
        dt=dtim
        dt2=dtim/2.0d0
      endif

      x1=x2
      y1=y2
      z1=z2

      t1=t2

      vx1=vx2
      vy1=vy2
      vz1=vz2

      bx1=bx2
      by1=by2
      bz1=bz2

      x2b=x1+vx1*dt2
      y2b=y1+vy1*dt2
      z2b=z1+vz1*dt2

      call uradfield_omp(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,gamma,istatus,
     &  modewave)
      if (istatus.ne.0) ifail=ifail+abs(istatus)
      istatus=0

      call uradstep_omp(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,efx2,efy2,efz2,
     &  dtim,x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,ieneloss,
     &  dgamma)

      t2=t1+dtim

      if (ieneloss.ne.0) then
        dgamsum=dgamsum+vxsign*dgamma
        if (abs(dgamsum).gt.gamma*1.0d-8) then
          gamma=gamma+dgamsum
          dgamtot=dgamtot+dgamsum
          dgamsum=0.0d0
        endif
        beta=dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))
        vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
        vx2=vx2/vn*clight*beta
        vy2=vy2/vn*clight*beta
        vz2=vz2/vn*clight*beta
      endif

      BX=VX2*C1
      BY=VY2*C1
      BZ=VZ2*C1

      BPX=VXP*C1
      BPY=VYP*C1
      BPZ=VZP*C1

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

      RX=XOBSV-X2
      RY=YOBSV-Y2
      RZ=ZOBSV-Z2

      R=SQRT(RX*RX+RY*RY+RZ*RZ)
      R1=1.0D0/R
      ZICR1=ZIC*R1

      RNX=RX*R1
      RNY=RY*R1
      RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

      BET1N=(1.0D0-BX*RNX)-BY*RNY-BZ*RNZ

      br2=by**2+bz**2
      rnr2=rny**2+rnz**2
      b3=beta**3
      br4=br2**2
      rnr4=rnr2**2

      if(br2.lt.1.0d-4.and.rnr2.lt.1.0d-4) then
        bet1n=
     &    1.0d0/(1.0d0+beta)/gamma**2
     &    +beta*(rnr2/2.0d0
     &    +rnr4/8.0d0)
     &    +(br2/2.0d0
     &    -br2*rnr2/4.0d0
     &    -br2*rnr4/16.0d0)/beta
     &    +b3*br4*(1.0d0/8.0d0
     &    -rnr2/16.0d0
     &    -rnr4/64.0d0)
     &    -by*rny
     &    -bz*rnz
      endif

      DUM11=1.0D0/BET1N
      DOM1=1.0D0/(R*BET1N*BET1N)

      RNBX=RNX-BX
      RNBY=RNY-BY
      RNBZ=RNZ-BZ

      PX=(RNBY*BPZ-RNBZ*BPY)
      PY=(RNBZ*BPX-RNBX*BPZ)
      PZ=(RNBX*BPY-RNBY*BPX)

      IF (IVELOFIELD.EQ.0.OR.IVELOFIELD.EQ.2) THEN !2 WEGEN POWER
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
        RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
      ELSE IF (IVELOFIELD.EQ.1) THEN
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
        RARG(3)=(RNX*PY-RNY*PX)*DOM1
      ELSE IF (IVELOFIELD.LT.0) THEN
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNX-BX)*DOM2
        RARG(2)=(RNY-BY)*DOM2
        RARG(3)=(RNZ-BZ)*DOM2
      ELSE  !IVELOFIELD
        istatus=-3
        return
      ENDIF !IVELOFIELD

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS      RARG(4)=T+R*C1

      DPHASE=BET1N*DT

      RARG(5)=(RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11

      powden=powden+rarg(5)*dt

C REAL PART OF INTEGRAND }

C COMPLEX PART OF INTEGRAND {

C    ASSUMES phener(I+1)=2*phener(I)   FOR kfreq2P=2
C    OR phener(I+1)=phener(I)+DELTA    FOR kfreq2P>2

C--- LOOP OVER ALL FREQUENCES

      if (nphener.gt.0) then

        kfreq=1

        OM=phener(kfreq)/hbarev
        ZIOM=ZI*OM

        if (izaehl.eq.1) then
          EXPOM1=CDEXP(DCMPLX(0.D0,phase*OM))
        endif

        EXPOM=EXPOM1
        DEXPOMPH1=EXP(ZIOM*DPHASE)
        DEXPOMPH=DEXPOMPH1

        IF(nphener.GT.1) THEN
          DEXPOM=EXP(ZIDOM*PHASE)
          DDEXPOMPH=EXP(ZIDOM*DPHASE)
        ENDIF  !kfreq2P

        IF (IVELOFIELD.NE.2) THEN

          dum3=expom*(zone-dexpomph)/om/bet1n

          ampex=rarg(1)*dum3
          ampey=rarg(2)*dum3
          ampez=rarg(3)*dum3

c          if (idebug.ne.0.and.yobsv.eq.0.0d0.and.zobsv.eq.0.0d0) then
c            write(77,*)izaehl,kfreq,x2,t2,rarg,dreal(ampez),dimag(ampez)
c          endif

          aradex(kfreq)=aradex(kfreq)+ampex
          aradey(kfreq)=aradey(kfreq)+ampey
          aradez(kfreq)=aradez(kfreq)+ampez

c          ampbx=conjg(rny*ampez-rnz*ampey)/clight
c          ampby=conjg(rnz*ampex-rnx*ampez)/clight
c          ampbz=conjg(rnx*ampey-rny*ampex)/clight
          ampbx=(rny*ampez-rnz*ampey)/clight
          ampby=(rnz*ampex-rnx*ampez)/clight
          ampbz=(rnx*ampey-rny*ampex)/clight

          aradbx(kfreq)=aradbx(kfreq)+ampbx
          aradby(kfreq)=aradby(kfreq)+ampby
          aradbz(kfreq)=aradbz(kfreq)+ampbz

        ELSE !IVELOFIELD

          EXPOMV2=R1/BET1N*EXPOM*(ZONE-DEXPOMPH)
          ZIOMR1=ZONE+ZICR1/OM

          ampex=EXPOMV2*(BX-RNX*ZIOMR1)
          ampey=EXPOMV2*(BX-RNX*ZIOMR1)
          ampez=EXPOMV2*(BX-RNX*ZIOMR1)

          aradex(kfreq)=aradex(kfreq)+ampex
          aradey(kfreq)=aradey(kfreq)+ampey
          aradez(kfreq)=aradez(kfreq)+ampez

c          ampbx=conjg(rny*ampez-rnz*ampey)/clight
c          ampby=conjg(rnz*ampex-rnx*ampez)/clight
c          ampbz=conjg(rnx*ampey-rny*ampex)/clight

          ampbx=(rny*ampez-rnz*ampey)/clight
          ampby=(rnz*ampex-rnx*ampez)/clight
          ampbz=(rnx*ampey-rny*ampex)/clight

          aradbx(kfreq)=aradbx(kfreq)+ampbx
          aradby(kfreq)=aradby(kfreq)+ampby
          aradbz(kfreq)=aradbz(kfreq)+ampbz

        ENDIF !IVELOFIELD

        IF (IVELOFIELD.NE.2) THEN

          DO kfreq=2,nphener

            OM=OM+DOM
            EXPOM=EXPOM*DEXPOM
            DEXPOMPH=DEXPOMPH*DDEXPOMPH

            EXPOMV2=1.0D0/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)

            ampex=RARG(1)*EXPOMV2
            ampey=RARG(2)*EXPOMV2
            ampez=RARG(3)*EXPOMV2

            aradex(kfreq)=aradex(kfreq)+ampex
            aradey(kfreq)=aradey(kfreq)+ampey
            aradez(kfreq)=aradez(kfreq)+ampez

c            ampbx=conjg(rny*ampez-rnz*ampey)/clight
c            ampby=conjg(rnz*ampex-rnx*ampez)/clight
c            ampbz=conjg(rnx*ampey-rny*ampex)/clight

            ampbx=(rny*ampez-rnz*ampey)/clight
            ampby=(rnz*ampex-rnx*ampez)/clight
            ampbz=(rnx*ampey-rny*ampex)/clight

            aradbx(kfreq)=aradbx(kfreq)+ampbx
            aradby(kfreq)=aradby(kfreq)+ampby
            aradbz(kfreq)=aradbz(kfreq)+ampbz

c            if (kfreq.eq.nphener/2) then
c              if (user(1).eq.0) then
c                write(56,*)ical,izaehl,x2,dphase,phase,real(aradez(kfreq)),real(RARG(3)*EXPOMV2)
c              else
c                write(57,*)ical,izaehl,x2,dphase,phase,real(aradez(kfreq)),real(RARG(3)*EXPOMV2)
c              endif
c            endif

          ENDDO   !LOOP OVER ALL FREQUENCES

        else

          DO kfreq=2,nphener

            OM=OM+DOM
            EXPOM=EXPOM*DEXPOM
            DEXPOMPH=DEXPOMPH*DDEXPOMPH

            EXPOMV2=R1/BET1N*EXPOM*(ZONE-DEXPOMPH)
            ZIOMR1=ZONE+ZICR1/OM

            ampex=EXPOMV2*(BX-RNX*ZIOMR1)
            ampey=EXPOMV2*(BX-RNX*ZIOMR1)
            ampez=EXPOMV2*(BX-RNX*ZIOMR1)

            aradex(kfreq)=aradex(kfreq)+ampex
            aradey(kfreq)=aradey(kfreq)+ampey
            aradez(kfreq)=aradez(kfreq)+ampez

c            ampbx=conjg(rny*ampez-rnz*ampey)/clight
c            ampby=conjg(rnz*ampex-rnx*ampez)/clight
c            ampbz=conjg(rnx*ampey-rny*ampex)/clight

            ampbx=(rny*ampez-rnz*ampey)/clight
            ampby=(rnz*ampex-rnx*ampez)/clight
            ampbz=(rnx*ampey-rny*ampex)/clight

            aradbx(kfreq)=aradbx(kfreq)+ampbx
            aradby(kfreq)=aradby(kfreq)+ampby
            aradbz(kfreq)=aradbz(kfreq)+ampbz

          ENDDO   !LOOP OVER ALL FREQUENCES

        ENDIF !IVELOFIELD


C COMPLEX PART OF INTEGRAND }

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

        PHASE=PHASE+DPHASE
        EXPOM1=EXPOM1*DEXPOMPH1

      endif !(nphener.gt.0) then

c ef is normal vector of perpendiculare plane at the end of the reference orbit
c dist is distance of electron to this plane
c tracking stops if trajectory hits this plane

      dist2=(x2-xf)*efx+(y2-yf)*efy+(z2-zf)*efz

      if (lstep.eq.0.and.dist2.lt.0.0d0.and.dist2.gt.-2.0d0*ds)  then

        lstep=1
        goto 1000

      else

        kstep=kstep+1
        if (kstep.eq.nthstep) kstep=0

        if(kstep.eq.0) then

          nstep=nstep+1

          if (nstep.gt.ndim) then
            nstep=nstep-1
            istatus=-2
            goto 9000
          endif

          traxyz(1,nstep)=x2
          traxyz(2,nstep)=y2
          traxyz(3,nstep)=z2
          traxyz(4,nstep)=t2
          traxyz(5,nstep)=vx2/vn
          traxyz(6,nstep)=vy2/vn
          traxyz(7,nstep)=vz2/vn
          traxyz(8,nstep)=bx2
          traxyz(9,nstep)=by2
          traxyz(10,nstep)=bz2
          traxyz(11,nstep)=gamma
          traxyz(12,nstep)=efx2
          traxyz(13,nstep)=efy2
          traxyz(14,nstep)=efz2

        endif

        if (lstep.eq.1) goto 9000
        goto 1000

      endif !lstep and dist2

9000  continue

      xexit=x2
      yexit=y2
      zexit=z2

      vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
      vnxex=vx2/vn
      vnyex=vy2/vn
      vnzex=vz2/vn

      texit=t2

      specnor=
     &  banwid
     & /(4.0d0*pi**2*clight*hbarev)
     & /(4.0d0*pi*eps0)
     &  *current/1.0d6 !per mm**2

      pownor=echarge/16.0d0/pi/pi/eps0/clight*current/1.0d6 !W/mm**2

      rspn=sqrt(specnor)

      vstokes(1,1)=( 0.0d0,        0.0d0)      !horizontal polarization
      vstokes(1,2)=( 0.0d0,        0.0d0)
      vstokes(1,3)=(-1.0d0,       -1.0d0)

      vstokes(2,1)=( 0.0d0,        0.0d0)      !right handed polarization
      vstokes(2,2)=( 0.0d0,       -1.0d0)
      vstokes(2,3)=(+1.0d0,        0.0d0)

      vstokes(3,1)=( 0.0d0,        0.0d0)      !left handed polarization
      vstokes(3,2)=( 0.0d0,       -1.0d0)
      vstokes(3,3)=(-1.0d0,        0.0d0)

      vstokes(4,1)=( 0.0d0,        0.0d0)      !45 degree linear polarization
      vstokes(4,2)=( 1.0d0,        0.0d0)
      vstokes(4,3)=(-1.0d0,        0.0d0)

      do isto=1,4
        vsto=dsqrt
     &    (cdabs(vstokes(isto,1))**2
     &    +cdabs(vstokes(isto,2))**2
     &    +cdabs(vstokes(isto,3))**2)
        vstokes(isto,1)=vstokes(isto,1)/vsto
        vstokes(isto,2)=vstokes(isto,2)/vsto
        vstokes(isto,3)=vstokes(isto,3)/vsto

      enddo

      do kfreq=1,nphener

        aradex(kfreq)=aradex(kfreq)*rspn
        aradey(kfreq)=aradey(kfreq)*rspn
        aradez(kfreq)=aradez(kfreq)*rspn

        apolh=
     &    aradex(kfreq)*conjg(vstokes(1,1))
     &    +aradey(kfreq)*conjg(vstokes(1,2))
     &    +aradez(kfreq)*conjg(vstokes(1,3))

        apolr=
     &    aradex(kfreq)*conjg(vstokes(2,1))
     &    +aradey(kfreq)*conjg(vstokes(2,2))
     &    +aradez(kfreq)*conjg(vstokes(2,3))

        apoll=
     &    aradex(kfreq)*conjg(vstokes(3,1))
     &    +aradey(kfreq)*conjg(vstokes(3,2))
     &    +aradez(kfreq)*conjg(vstokes(3,3))

        apol45=
     &    aradex(kfreq)*conjg(vstokes(4,1))
     &    +aradey(kfreq)*conjg(vstokes(4,2))
     &    +aradez(kfreq)*conjg(vstokes(4,3))

        stok1=real(
     &    apolr*conjg(apolr)+
     &    apoll*conjg(apoll))

        stok2=-stok1+
     &    real(2.*apolh*conjg(apolh))

        stok3=
     &    real(2.*apol45*conjg(apol45))-
     &    stok1

        stok4=real(
     &    apolr*conjg(apolr)-
     &    apoll*conjg(apoll))

        stokes(1,kfreq)=stok1
        stokes(2,kfreq)=stok2
        stokes(3,kfreq)=stok3
        stokes(4,kfreq)=stok4

      enddo !nphener

      powden=powden*pownor

      if (istatus.ge.0.and.ifail.ne.0) istatus=ifail

      return
      end
*CMZ :          22/09/2025  10.53.54  by  Michael Scheer
*CMZ :  4.02/00 10/09/2025  13.45.08  by  Michael Scheer
*CMZ :  4.01/07 18/10/2024  14.11.15  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  11.54.00  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      complex*16 :: cph00,ci=(0.0d0,1.0d0)
      real*8 specnor_si

      integer :: mthreads,ktime=0,kfreq,icbrill,iobfr,iobsv

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase_prop',1)

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr_u ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid_u !BW
     &  *1.0d-6 !m**2 -> mm**2

      if (ifieldprop_u.gt.0) then

        if (modepin_u.ne.1) then
          call urad_phase_prop_classic(mthreads)
        else
          call urad_phase_prop_mc(mthreads)
        endif

      else !(ifieldprop_u.gt.0)

        call urad_phase_prop_geo

      endif !(ifieldprop_u.gt.0)

      icbrill=nobsvprop_u/2+1

      if (globphaseprop_u.eq.9999.0d0) then
        do kfreq=1,nepho_u
          iobfr=icbrill+nobsvprop_u*(kfreq-1)
          cph00=aradprop_u(3,iobfr)/abs(aradprop_u(3,iobfr))
          do iobsv=1,nobsvprop_u
            iobfr=iobsv+nobsvprop_u*(kfreq-1)
            aradprop_u(:,iobfr)=aradprop_u(:,iobfr)/cph00
          enddo
        enddo
      else if (globphaseprop_u.ne.0.0d0) then
        aradprop_u=aradprop_u*exp(ci*globphaseprop_u)
      endif

      stokesprop_u=stokesprop_u*specnor_si

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase_prop',0)

      return
      end
*CMZ :          24/09/2025  12.48.34  by  Michael Scheer
*CMZ :  4.02/00 13/09/2025  10.12.52  by  Michael Scheer
*CMZ :  4.01/07 11/08/2024  15.26.17  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  09.37.27  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_classic(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0),a3(3),a3p(3)

      double complex, dimension(:), allocatable :: expom,dexpom,phshift

      double complex :: apolh,apolr,apoll,apol45
      double precision :: dx,dx2,dy,dyph,dzph,dz,y,z,omc,domc,phlowz,phlowy,dzy2,eps(6),
     &  dr,drred,da,x,xobs,yobs,zobs,rlambda1,ans,stok1,stok2,stok3,stok4,stoknor,enor,
     &  eabsmaxprop=-1.0d30

      integer :: ktime=0,i,
     &  mthreads,iy,iz,n,jy,jz,iobs,ieps,ifrq,iobfr,jobs,jobfr

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      nobsvprop_u=npinyprop_u*npinzprop_u

      allocate(expom(nobsv_u*nepho_u),dexpom(nobsv_u*nepho_u),phshift(nobsv_u))

c      aradprop_u=(0.0d0,0.0d0)

      if (npinyprop_u.gt.1) then
        dyph=pinhprop_u/1000.0d0/dble(npinyprop_u-1)
        phlowy=-pinhprop_u/2000.0d0
      else
        dyph=0.0d0
        phlowy=0.0d0
      endif

      if (npinzprop_u.gt.1) then
        dzph=pinwprop_u/1000.0d0/dble(npinzprop_u-1)
        phlowz=-pinwprop_u/2000.0d0
      else
        dzph=0.0d0
        phlowz=0.0d0
      endif

      da=pinw_u/1000.0d0*pinh_u/1000.0d0/dble(max(1,npinz_u-1)*max(1,npiny_u-1))

      n=0

      x=pinxprop_u/1000.0d0
      y=phlowy-dyph
      do iy=1,npinyprop_u
        y=y+dyph
        if (abs(y).lt.1.0d-12) y=0.0d0
        z=phlowz-dzph
        obsvyprop_u(iy)=y
        do iz=1,npinzprop_u
          n=n+1
          z=z+dzph
          if (abs(z).lt.1.0d-12) z=0.0d0
          if (iy.eq.1) obsvzprop_u(iz)=z
          obsvprop_u(1:3,n)=[x,y,z]
        enddo
      enddo

      omc=epho_u(1)/(hbarev1*clight1)

      if(nepho_u.gt.1) then
        domc=(epho_u(2)-epho_u(1))/(hbarev1*clight1)
      endif

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(domc,omc,da,obsvprop_u,obsv_u,nobsv_u,nobsvprop_u,epho_u,nepho_u,aradprop_u,arad_u)

!$OMP DO

      do jobs=1,nobsvprop_u
c        ith=OMP_GET_THREAD_NUM()+1

        x=obsvprop_u(1,jobs)
        y=obsvprop_u(2,jobs)
        z=obsvprop_u(3,jobs)

        DO IOBS=1,NOBSV_u

          XOBS=OBSV_u(1,IOBS)/1000.0d0
          YOBS=OBSV_u(2,IOBS)/1000.0d0
          ZOBS=OBSV_u(3,IOBS)/1000.0d0

          dx=xobs-x
          dx2=dx*dx
          DY=YOBS-y
          DZ=ZOBS-z
          DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

          IF (DZY2.GT.0.01D0*dx2) THEN
            WRITE(6,*)'*** ERROR IN URAD_phase_prop_classic ***'
            WRITE(6,*)'CHECK INPUT FILE AND INCREASE PinX'
            WRITE(6,*)'*** PROGRAM ABORTED ***'
            STOP
          ENDIF

          EPS(1)=DZY2/dx2
          DO IEPS=2,6
            EPS(IEPS)=EPS(IEPS-1)*EPS(1)
          ENDDO !IEPS

c      TAYLOR-EXPANSION DONE WITH REDUCE
c     IN "WTAY1.RED";
c     on rounded;
c     on numval;
c     precision 13;
c     F:=SQRT(1+EPS);
c     DR:=TAY1(F,EPS,6);
c     ON FORT;
c     OUT "RED.FOR";
c     DR;
c     SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

          ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &      -0.0390625D0*eps(4)+
     &      0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

          DR=DABS(dx*(ANS+1.0D0))
          DRRED=-DABS(dx*ANS)

          IF (DR.NE.0.0d0) THEN
            EXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*OMC))/DR
          ELSE
            EXPOM(IOBS)=1.0D0
          ENDIF

          if (nepho_u.gt.1) then
            DEXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*DOMC))
          endif
c            print*,ith,iobs,expom(iobs)
c+seq,dum2.
        ENDDO   !NOBS

        DO ifrq=1,nepho_u

          jOBFR=jOBS+NOBSVprop_u*(ifrq-1)

          RLAMBDA1=epho_u(ifrq)/WTOE1*1.0D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)

          DO IOBS=1,NOBSV_u

            IOBFR=IOBS+NOBSV_u*(ifrq-1)

            IF (ifrq.EQ.1) THEN
              PHSHIFT(IOBS)=EXPOM(IOBS)
            ELSE
              PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
            ENDIF   !(ifrq.EQ.1)

            if (dx.gt.0.0d0) then
              aradprop_u(1:6,jobfr)=aradprop_u(1:6,jobfr)+
     &          arad_u(1:6,iobfr)*PHSHIFT(IOBS)*da*rlambda1
            else
              aradprop_u(1:6,jobfr)=aradprop_u(1:6,jobfr)+
     &          dconjg(arad_u(1:6,iobfr))*PHSHIFT(IOBS)*da*rlambda1
            endif
          ENDDO   !NFREQ

        ENDDO  !NOBSV

      ENDDO !nobsvprop_u

!$OMP END DO

!$OMP END PARALLEL

      eabsmaxprop=-1.0d30

      if (ifieldprop_u.ne.2) then

        do ifrq=1,nepho_u
          do jobs=1,nobsvprop_u

            jobfr=jobs+nobsvprop_u*(ifrq-1)

            apolh=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(1,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(1,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(1,3))

            apolr=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(2,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(2,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(2,3))

            apoll=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(3,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(3,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(3,3))

            apol45=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(4,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(4,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(4,3))

            stok1=dreal(
     &        apolr*conjg(apolr)+
     &        apoll*conjg(apoll))

            stok2=-stok1+
     &        dreal(2.*apolh*conjg(apolh))

            stok3=
     &        dreal(2.*apol45*conjg(apol45))-
     &        stok1

            stok4=dreal(
     &        apolr*conjg(apolr)-
     &        apoll*conjg(apoll))

            stokesprop_u(1,jobfr)=stok1
            stokesprop_u(2,jobfr)=stok2
            stokesprop_u(3,jobfr)=stok3
            stokesprop_u(4,jobfr)=stok4

          enddo
        enddo

      endif !(ifieldprop_u.ne.2)

      IOBFR=nobsv_u/2+1+NOBSV_u*(nepho_u/2)
      jobfr=nobsvprop_u/2+1+nobsvprop_u*(nepho_u/2)

      a3=arad_u(1:3,iobfr)
      a3p=aradprop_u(1:3,jobfr)

c      print*,sqrt(norm2(dreal(a3*dconjg(a3)))),sqrt(norm2(dreal(a3p*dconjg(a3p))))
c      enor=sqrt(norm2(dreal(a3*dconjg(a3)))/norm2(dreal(a3p*dconjg(a3p))))
c      print*,enor
c      enor=norm2(dreal(a3*dconjg(a3)))/norm2(dreal(a3p*dconjg(a3))))
c      print*,enor
c      stoknor=enor**2
c      stop
c      stoknor=enor**2
c      stokesprop_u=stokesprop_u/stoknor
c      aradprop_u=aradprop_u/enor

c      print*,enor,1.0d0/enor,stoknor,1.0d0/stoknor
c      obsvprop_u=obsvprop_u*1000.0d0

      deallocate(expom,dexpom,phshift)

      return
      end
*CMZ :  4.01/07 05/09/2024  16.27.42  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  09.37.27  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_geo

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0),a3(3)

      double complex :: expom,apolh,apolr,apoll,apol45

      double precision dx,dx2,dy,dyph,dzph,dz,y,z,omc,domc,phlowz,phlowy,dzy2,eps(6),
     &  dr,drred,da,x,xobs,yobs,zobs,rlambda1,ans,stok1,stok2,stok3,stok4,rn(3)

      integer :: ifrq,iobph,iobs,ieps,iobfr

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      nobsvprop_u=npinyprop_u*npinzprop_u

      x=pinxprop_u/1000.0d0

      do ifrq=1,nepho_u

        omc=epho_u(ifrq)/(hbarev1*clight1)

        do iobs=1,nobsv_u

          iobph=iobs+nobsv_u*(ifrq-1)

          xobs=obsv_u(1,iobs)/1000.0d0
          yobs=obsv_u(2,iobs)/1000.0d0
          zobs=obsv_u(3,iobs)/1000.0d0

          dx=xobs-x
          dx2=dx*dx

          rn(1)=real(arad_u(2,iobph)*conjg(arad_u(6,iobph))-arad_u(3,iobph)*conjg(arad_u(5,iobph)))
          rn(2)=real(arad_u(3,iobph)*conjg(arad_u(4,iobph))-arad_u(1,iobph)*conjg(arad_u(6,iobph)))
          rn(3)=real(arad_u(1,iobph)*conjg(arad_u(5,iobph))-arad_u(2,iobph)*conjg(arad_u(4,iobph)))
          rn=rn/norm2(rn)

          y=yobs-rn(2)/rn(1)*dx
          z=zobs-rn(3)/rn(1)*dx

          obsvprop_u(1:3,iobph)=[x,y,z]
          obsvprop_u(4:6,iobph)=rn

          dy=yobs-y
          dz=zobs-z
          dzy2=dz*dz+dy*dy
          da=dx2+dzy2 !convert to solid angle

          if (dzy2.le.0.01d0*dx2) THEN

            eps(1)=dzy2/dx2
            do ieps=2,6
              eps(ieps)=eps(ieps-1)*eps(1)
            enddo !ieps

c      TAYLOR-EXPANSION DONE WITH REDUCE
c     IN "WTAY1.RED";
c     on rounded;
c     on numval;
c     precision 13;
c     F:=SQRT(1+EPS);
c     DR:=TAY1(F,EPS,6);
c     ON FORT;
c     OUT "RED.FOR";
c     DR;
c     SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

            ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &        -0.0390625D0*eps(4)+
     &        0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

            dr=dabs(dx*(ans+1.0d0))
            drred=-dabs(dx*ans)

          else
            dr=abs(dx)*sqrt(1.0d0+eps(1))
            drred=dx-dr
          endif

          rlambda1=epho_u(ifrq)/wtoe1*1.0d9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)
          expom=cdexp(dcmplx(0.0d0,drred*omc))/dr

          if (dx.gt.0.0d0) then
            aradprop_u(1:6,iobph)=arad_u(1:6,iobph)*da*expom*rlambda1
          else
            aradprop_u(1:6,iobph)=dconjg(arad_u(1:6,iobph))*da*expom*rlambda1
          endif

        ENDDO   !NFREQ

      ENDDO   !NOBS

      do ifrq=1,nepho_u
        do iobs=1,nobsvprop_u

          iobfr=iobs+nobsvprop_u*(ifrq-1)

          apolh=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(1,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(1,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(1,3))

          apolr=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(2,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(2,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(2,3))

          apoll=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(3,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(3,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(3,3))

          apol45=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(4,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(4,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(4,3))

          stok1=dreal(
     &      apolr*conjg(apolr)+
     &      apoll*conjg(apoll))

          stok2=-stok1+
     &      dreal(2.*apolh*conjg(apolh))

          stok3=
     &      dreal(2.*apol45*conjg(apol45))-
     &      stok1

          stok4=dreal(
     &      apolr*conjg(apolr)-
     &      apoll*conjg(apoll))

          stokesprop_u(1,iobfr)=stok1
          stokesprop_u(2,iobfr)=stok2
          stokesprop_u(3,iobfr)=stok3
          stokesprop_u(4,iobfr)=stok4

        enddo
      enddo

      obsvprop_u=obsvprop_u*1000.0d0

      return
      end
*CMZ :  4.01/07 06/08/2024  14.45.08  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  21.58.40  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_mc(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double complex, dimension(:), allocatable :: expom,dexpom,phshift
      double complex :: apolh,apolr,apoll,apol45

      double precision, dimension(:,:), allocatable :: obsv
      double precision dx,dx2,dy,dyph,dzph,dz,y,z,omc,domc,phlowz,phlowy,dzy2,eps(6),
     &  dr,drred,da,x,xobs,yobs,zobs,rlambda1,ans,stok1,stok2,stok3,stok4,zmin,ymin,fsum

      integer :: mthreads,iy,iz,n,jy,jz,iobs,ieps,ifrq,iobfr,jobs,jobfr,nobsv,iepho

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      nobsv=(npinzo_u-1)*(npinyo_u-1)
      nobsvprop_u=npinyprop_u*npinzprop_u

      allocate(expom(nobsv*nepho_u),dexpom(nobsv*nepho_u),phshift(nobsv),obsv(3,nobsv))

      ymin=pincen_u(2)-pinh_u/2.0d0
      zmin=pincen_u(3)-pinw_u/2.0d0

      dy=pinh_u/max(1,npinyo_u-1)
      dz=pinw_u/max(1,npinzo_u-1)

      iobs=0
      do iy=1,npinyo_u-1
        do iz=1,npinzo_u-1
          iobs=iobs+1
          obsv(1,iobs)=pincen_u(1)
          obsv(2,iobs)=ymin+(dble(iy)-0.5)*dy
          obsv(3,iobs)=zmin+(dble(iz)-0.5)*dz
        enddo
      enddo

      aradprop_u=(0.0d0,0.0d0)

      if (npinyprop_u.gt.1) then
        dyph=pinhprop_u/1000.0d0/dble(npinyprop_u-1)
        phlowy=-pinhprop_u/2000.0d0
      else
        dyph=pinhprop_u/1000.0d0
        phlowy=-pinhprop_u/1000.0d0
      endif

      if (npinzprop_u.gt.1) then
        dzph=pinwprop_u/1000.0d0/dble(npinzprop_u-1)
        phlowz=-pinwprop_u/2000.0d0
      else
        dzph=pinwprop_u/1000.0d0
        phlowz=-pinwprop_u/1000.0d0
      endif

      da=dz*dy

      n=0

      x=pinxprop_u/1000.0d0
      y=phlowy-dyph

      do iy=1,npinyprop_u
        y=y+dyph
        if (abs(y).lt.1.0d-12) y=0.0d0
        z=phlowz-dzph
        do iz=1,npinzprop_u
          n=n+1
          z=z+dzph
          if (abs(z).lt.1.0d-12) z=0.0d0
          obsvprop_u(1:3,n)=[x,y,z]
        enddo
      enddo

      omc=epho_u(1)/(hbarev1*clight1)
      if(nepho_u.gt.1) then
        domc=(epho_u(2)-epho_u(1))/(hbarev1*clight1)
      endif

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(domc,omc,da,obsvprop_u,obsv,npinzo_u,npinyo_u,nobsv,nobsvprop_u,epho_u,nepho_u,aradprop_u,fieldbunch)

!$OMP DO

      do jobs=1,nobsvprop_u
c        ith=OMP_GET_THREAD_NUM()+1

        x=obsvprop_u(1,jobs)
        y=obsvprop_u(2,jobs)
        z=obsvprop_u(3,jobs)

        DO IOBS=1,nobsv

          XOBS=obsv(1,IOBS)
          YOBS=obsv(2,IOBS)
          ZOBS=obsv(3,IOBS)

          dx=xobs-x
          dx2=dx*dx
          DY=YOBS-y
          DZ=ZOBS-z
          DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

          IF (DZY2.GT.0.01D0*dx2) THEN
            WRITE(6,*)'*** ERROR IN URAD_phase_prop_mc ***'
            WRITE(6,*)'CHECK INPUT FILE AND INCREASE PinX'
            WRITE(6,*)'*** PROGRAM ABORTED ***'
            STOP
          ENDIF

          EPS(1)=DZY2/dx2
          DO IEPS=2,6
            EPS(IEPS)=EPS(IEPS-1)*EPS(1)
          ENDDO !IEPS

c      TAYLOR-EXPANSION DONE WITH REDUCE
c     IN "WTAY1.RED";
c     on rounded;
c     on numval;
c     precision 13;
c     F:=SQRT(1+EPS);
c     DR:=TAY1(F,EPS,6);
c     ON FORT;
c     OUT "RED.FOR";
c     DR;
c     SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

          ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &      -0.0390625D0*eps(4)+
     &      0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

          DR=DABS(dx*(ANS+1.0D0))
          DRRED=-DABS(dx*ANS)

          IF (DR.NE.0.0d0) THEN
            EXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*OMC))/DR
          ELSE
            EXPOM(IOBS)=1.0D0
          ENDIF

          if (nepho_u.gt.1) then
            DEXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*DOMC))
          endif
c            print*,ith,iobs,expom(iobs)
c+seq,dum2.
        ENDDO   !NOBS

        DO ifrq=1,nepho_u

          jOBFR=jOBS+NOBSVprop_u*(ifrq-1)

          RLAMBDA1=epho_u(ifrq)/WTOE1*1.0D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)

          iobs=0
          do iy=1,npinyo_u-1
            do iz=1,npinzo_u-1
              iobs=iobs+1

              IOBFR=IOBS+nobsv*(ifrq-1)

              IF (ifrq.EQ.1) THEN
                PHSHIFT(IOBS)=EXPOM(IOBFR)
              ELSE
                PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
              ENDIF   !(ifrq.EQ.1)

              if (dx.gt.0.0d0) then
                aradprop_u(1:6,jobfr)=aradprop_u(1:6,jobfr)+
     &            fieldbunch(1:6,iz,iy,ifrq)*PHSHIFT(IOBS)*da*rlambda1
              else
                aradprop_u(1:6,jobfr)=aradprop_u(1:6,jobfr)+
     &            dconjg(fieldbunch(1:6,iz,iy,ifrq))*PHSHIFT(IOBS)*da*rlambda1
              endif
            ENDDO   !NFREQ

          ENDDO  !NOBSVz
        ENDDO  !NOBSVy

      ENDDO !nobsvprop_u

!$OMP END DO

!$OMP END PARALLEL

      if (ifieldprop_u.ne.2) then

        do ifrq=1,nepho_u
          do jobs=1,nobsvprop_u

            jobfr=jobs+nobsvprop_u*(ifrq-1)

            apolh=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(1,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(1,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(1,3))

            apolr=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(2,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(2,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(2,3))

            apoll=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(3,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(3,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(3,3))

            apol45=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(4,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(4,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(4,3))

            stok1=dreal(
     &        apolr*conjg(apolr)+
     &        apoll*conjg(apoll))

            stok2=-stok1+
     &        dreal(2.*apolh*conjg(apolh))

            stok3=
     &        dreal(2.*apol45*conjg(apol45))-
     &        stok1

            stok4=dreal(
     &        apolr*conjg(apolr)-
     &        apoll*conjg(apoll))

            stokesprop_u(1,jobfr)=stok1
            stokesprop_u(2,jobfr)=stok2
            stokesprop_u(3,jobfr)=stok3
            stokesprop_u(4,jobfr)=stok4

          enddo
        enddo

      endif !(ifieldprop_u.ne.2)

      obsvprop_u=obsvprop_u*1000.0d0

      return
      end
*CMZ :  4.01/07 07/08/2024  08.52.56  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  21.58.40  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_point(sourcepoint,field,
     &  nzprop,nyprop,xprop,yprop,zprop,pinw,pinh,epho,fprop)

      implicit none

      integer nzprop,nyprop

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double complex expom,apolh,apol45,apoll,apolr
      double complex :: field(3),cjfield(3),fprop(3,nzprop,nyprop)

      double precision sourcepoint(3),xprop,yprop(nyprop),zprop(nzprop),pinw,pinh,epho

      double precision dx,dx2,dy,dz,y,z,omc,domc,dzy2,eps(6),
     &  dr,drred,x,xobs,yobs,zobs,darlambda1,ans,stok1,stok2,stok3,stok4,zmin,ymin,fsum

      integer :: iy,iz,n,jy,jz,iobs,ieps,iobfr,jobs,jobfr,nobsv

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      fprop=(0.0d0,0.0d0)

      if (nyprop.gt.1) then
        dy=yprop(2)-yprop(1)
      else
        dy=pinh/1000.0d0
      endif

      if (nzprop.gt.1) then
        dz=zprop(2)-zprop(1)
      else
        dz=pinw/1000.0d0
      endif

      omc=epho/(hbarev1*clight1)
      x=xprop

      xobs=sourcepoint(1)
      yobs=sourcepoint(2)
      zobs=sourcepoint(3)

      daRLAMBDA1=dz*dy*epho/WTOE1*1.0D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)
      cjfield=dconjg(field)

      do iy=1,nyprop

        y=yprop(iy)

        do iz=1,nzprop

          z=zprop(iz)

          dx=xobs-x
          dx2=dx*dx
          DY=YOBS-y
          DZ=ZOBS-z
          DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

          IF (DZY2.GT.0.01D0*dx2) THEN
            WRITE(6,*)'*** ERROR IN URAD_phase_prop_mc ***'
            WRITE(6,*)'CHECK INPUT FILE AND INCREASE PinX'
            WRITE(6,*)'*** PROGRAM ABORTED ***'
            STOP
          ENDIF

          EPS(1)=DZY2/dx2
          DO IEPS=2,6
            EPS(IEPS)=EPS(IEPS-1)*EPS(1)
          ENDDO !IEPS

          ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &      -0.0390625D0*eps(4)+
     &      0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

          DR=DABS(dx*(ANS+1.0D0))
          DRRED=-DABS(dx*ANS)

          IF (DR.NE.0.0d0) THEN
            EXPOM=CDEXP(DCMPLX(0.0d0,DRRED*OMC))/DR
          ELSE
            EXPOM=1.0D0
          ENDIF

          if (dx.gt.0.0d0) then
            fprop(1:3,iz,iy)=fprop(1:3,iz,iy)+field*expom*darlambda1
          else
            fprop(1:3,iz,iy)=dconjg(fprop(1:3,iz,iy))+cjfield*expom*darlambda1
          endif

c          write(77,*)sourcepoint,z,y,dreal(fprop(1:3,iz,iy)),dimag(fprop(1:3,iz,iy))

        ENDDO !nzprop
      enddo !nyprop

      return
      end
*CMZ :  4.01/07 10/08/2024  16.46.49  by  Michael Scheer
*CMZ :  4.01/05 16/04/2024  14.41.20  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.32.18  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  10.57.05  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  13.32.32  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.57.49  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_fold_2d(nx,ny,x,y,fin,sigx,sigy,fold)

      implicit none

      integer ix,iy,nx,ny

      real*8 fin(nx,ny),x(nx),y(ny),sigx,sigy,fold(nx,ny),
     &  fxf(nx,ny),f(max(nx,ny)),fg(max(nx,ny)),ws1(max(nx,ny)),ws2(max(nx,ny))

      do iy=1,ny
        f(1:nx)=fin(1:nx,iy)
        if (sigx.gt.0.0d0) then
          call util_fold_function_gauss_lin(nx,x,f,sigx,3.0d0,fg,ws1,ws2)
          fxf(1:nx,iy)=fg(1:nx)
        else
          fxf(1:nx,iy)=fin(1:nx,iy)
        endif
      enddo

      do ix=1,nx
        if (sigy.gt.0.0d0) then
          f(1:ny)=fxf(ix,1:ny)
          call util_fold_function_gauss_lin(ny,y,f,sigy,3.0d0,fg,ws1,ws2)
          fold(ix,1:ny)=fg(1:ny)
        else
          fold(ix,1:ny)=fxf(ix,1:ny)
        endif
      enddo

      return
      end
*CMZ :  4.02/00 11/09/2025  15.05.09  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2025
      subroutine undulator_wigner_num(nx,ny,dx,dy,wlen,esour,ntx,nty,thex,they,wig,curr,banwid)

      use omp_lib

      implicit none

      include 'phyconparam.cmn'

      integer, intent(in) :: nx,ny,ntx,nty

      complex*16 :: ci=(0.0d0,1.0d0),exptx,expdtx,
     &  cthe,eki,expom,em,ep

      complex*16, intent(in) :: esour(2,nx,ny)
c      complex*16  :: wkern(2*nx,2*ny,nx,ny)
c      complex*16 :: wigc(nx,ny,ntx,nty)

      real*8, intent(in):: wlen,thex(ntx),they(nty),dx,dy
      real*8, intent(out) :: wig(nx,ny,ntx,nty)

      real*8 :: tx,ty,ek,dtx,dty,xm,xp,xpm,ypm,ym,yp,rp,rm,x,y,wlen12,wignor,curr,banwid,
     &  specnor_si

      integer :: ix,iy,itx,ity,kx,ky,jfail,nper,iypm,ixpm,lx,ly,ifound,nmaxth=1

c      real secin,secout
c      print*
c      print*,"     Calculating Wigner Distribution for ",sngl(wlen)," nm"
c      print*

c      secin=secnds(0.0)

      wlen12=1.0d0/(wlen*1.0d-9)**2
      ek=twopi1/abs(wlen*1.0d-9) !1/m
      eki=ci*ek

      wig=0.0d0
c      wigc=(0.0d0,0.0d0)

      if (ntx.gt.1) then
        dtx=thex(2)-thex(1)
      else
        dtx=1.0d0
      endif

      if (nty.gt.1) then
        dty=they(2)-they(1)
      else
        dty=1.0d0
      endif

      nmaxth=OMP_GET_MAX_THREADS()

c      call undulator_wigner_kernel(nx,ny,esour,wkern) bringt's nicht

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid

      wignor=specnor_si*dx*dy*wlen12 !??*4.0d0 ! factor 4 due to change in integration variables!?

c      nmaxth=1
c      print*,"Nur ein Thread!"

c      do iy=1,ny
c        do ix=1,nx
c          write(77,*)ix,iy,
c     &      dreal(esour(1,ix,iy)),dimag(esour(2,ix,iy)),
c     &      dreal(esour(2,ix,iy)),dimag(esour(2,ix,iy))
c        enddo
c      enddo
c      flush(77)
c      close(77)

!$OMP PARALLEL NUM_THREADS(nmaxth) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(nx,ny,ntx,nty,dx,dy,eki,wlen12,they,thex,dtx,dty,wignor)
!$OMP& SHARED(esour,wig)

!$OMP DO

      ! Im Zentrum des Undulators, Gl. 77

      !E-Feld V/m = 1.e-4 / (clight1/1.e8) statvolt/cm) =~ 1.e-4/3.

      do iy=1,ny

        do ix=1,nx

          do iypm=-ny+1,ny-1

            ypm=dy*iypm

            ky=iy-iypm/2
            ly=iy+iypm/2

            if (ky.lt.1.or.ky.gt.ny) cycle
            if (ly.gt.ny.or.ly.lt.1) cycle

            do ity=1,nty
              ty=they(ity)

              do ixpm=-nx+1,nx-1

                xpm=dx*ixpm

                kx=ix-ixpm/2
                lx=ix+ixpm/2

                if (kx.lt.1.or.kx.gt.nx) cycle
                if (lx.gt.nx.or.lx.lt.1) cycle

                do itx=1,ntx

                  em=esour(1,kx,ky)
                  ep=esour(2,lx,ly)

                  if (wlen.gt.0.0d0) then
                    if (itx.eq.1) then
                      tx=thex(itx)
                      expom=exp(-eki*(xpm*tx+ypm*ty))
                      expdtx=exp(-eki*xpm*dtx)
                    else
                      expom=expom*expdtx
                    endif

                  else
                    if (itx.eq.1) then
                      tx=thex(itx)
                      expom=exp(eki*(xpm*tx+ypm*ty))
                      expdtx=exp(eki*xpm*dtx)
                    else
                      expom=expom*expdtx
                    endif
                  endif

C                  expom=exp(-eki*(xpm*tx+ypm*ty))

                  wig(ix,iy,itx,ity)=wig(ix,iy,itx,ity)+
     &              dreal(em*ep*expom)*wignor
c                  wigc(ix,iy,itx,ity)=wigc(ix,iy,itx,ity)+
c     &              em*ep*expom*wignor

                enddo
              enddo

            enddo
          enddo

        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c      wig=dreal(wig)
c      secout=secnds(0.0)

c      print*,"     Seconds used:",secout-secin

      end
*CMZ :          28/09/2025  14.46.51  by  Michael Scheer
*CMZ :  4.02/00 11/09/2025  13.56.35  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2025
      subroutine undulator_source_analytic(curr,banwid,defl,gamma,nx,ny,pinx,pinw,pinh,om,
     &  perlen,nper,dist,offset,slope,esour,nord,ajj,specnor,enor,fdmax,flux)

      implicit none

      include 'phyconparam.cmn'

      integer nx,ny,nord

      complex*16 :: ci=(0.0d0,1.0d0),esour(nx,ny),cthe,expos

      real*8 defl,gamma,pinx,pinw,pinh,om,perlen,devlen,dist,Com,eharm,omharm,wlharm,specnor
      real*8 dsinin,sinc,echarge_CGS,clight_CGS,hbar_CGS,fdmax,flux,curr,banwid
     &  ,specnor_SI,enorbase_CGS,enorbase_SI,enorbase,specnor_CGS,offset(2),slope(2),
     &  besin,unduk

      real*8 :: r2,dx,dy,ajj,enor,wlc,besn1,besn,x,y,the2,arg,dom,eps=1.0d-12
      integer :: ix,iy,jfail,nper,ioffs,islope

      echarge_CGS=4.8032d-10 !statcoulomb
      clight_CGS=clight1*100.0d0
      hbar_CGS=hbar1*1.0d7 !erg-sec

      ! Gl. 28 und Gl. 83
      !E-Feld V/m = 1.e-4 / (clight1/1.e8) statvolt/cm) =~ 1.e-4/3.
      SPECNOR_CGS=
     &  banwid
     &  *curr/echarge1 ! Strom bzw. Zahl der Elektronen pro Zeit
     &  *clight_CGS
     &  /hbar_CGS/4.0d0/pi1**2
     &  /1.0d6 ! rad**2->mrad**2

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid !BW
     &  /1.0d6 ! m**2->mm**2

c      besin=defl**2*dble(nord)/(4.0d0+2.0d0*defl**2)
      unduk=twopi1/perlen
      besin=om*defl**2/(8.0d0*unduk*clight1*gamma**2)

      call util_bessel((nord-1)/2,besin,besn1,jfail)
      call util_bessel((nord+1)/2,besin,besn,jfail)

      ajj=besn1-besn

      devlen=perlen*dble(nper)

      if (
     &    slope(1).ne.0.0d0 .or.slope(2).ne.0.0d0
     &    ) then
        islope=1
      else
        islope=0
      endif

      if (
     &    offset(1).ne.0.0d0.or.offset(2).ne.0.0d0.or.islope.ne.0
     &    ) then
        ioffs=1
      else
        ioffs=0
      endif

      wlharm=(1.0d0+defl**2/2.0d0)*perlen/2.0d0/gamma**2*1.0d9/dble(nord)
      eharm=wtoe1/wlharm
      omharm=eharm/hbarev1
      dom=om-omharm

      if (nx.gt.1) then
        dx=pinw/(nx-1)
      else
        dx=0.0d0
      endif

      if (ny.gt.1) then
        dy=pinh/(ny-1)
        y=-pinh/2.0d0
      else
        y=0.0d0
        dy=0.0d0
      endif

      if (dist.ne.0.0d0) then
        specnor_CGS=specnor_CGS*(pinx*100.0d0)**2
        enorbase_CGS=-defl*echarge_CGS/(2.0d0*clight_CGS**2*gamma)*ajj * devlen/dist ! Gl. (77)
        specnor_SI=specnor_SI*pinx**2
        enorbase=sqrt(specnor_CGS*enorbase_CGS**2/specnor_SI)
      else
        enorbase_CGS=defl*echarge_CGS/(2.0d0*clight_CGS**2*gamma)*ajj ! Gl. (75,77)
        enorbase=sqrt(specnor_CGS*10000.0d0*enorbase_CGS**2/specnor_SI) !10000: cm**2 -> m**2
      endif

c      enorbase_SI=echarge1/(4.0d0*pi1*eps01)/clight1**2 ! Gl. 1, korrigiert: c -> c**2
        specnor=specnor_SI

      if (dist.eq.0.0d0) then

        ! Im Zentrum des Undulators, Gl. 77

        !E-Feld V/m = 1.e-4 / (clight1/1.e8) statvolt/cm) =~ 1.e-4/3.
c        enor=defl*om*echarge_CGS/(2.0d0*clight_CGS**2*gamma)*ajj ! Gl. (75,77)
        enor=enorbase*om

        wlc=om/devlen/clight1 !1/m**2

        expos=(1.0d0,0.0d0)

        do iy=1,ny
          if (nx.gt.1) then
            x=-pinw/2.0d0
          else
            x=0.0d0
          endif
          do ix=1,nx
            if (islope.eq.1) then
              expos=exp(ci*om/clight1*(slope(1)*(x-offset(1))+slope(2)*(y-offset(2))))
            endif
            r2=(x-offset(1))**2+(y-offset(2))**2 !m**2
            if (r2.lt.eps) r2=0.0d0
            esour(ix,iy)=ci*enor*expos*(pi1-2.0d0*dsinin(wlc*r2)) !Gl. (77), virtual source E(0,r)
            x=x+dx
          enddo
          y=y+dy
        enddo

      else !dist.eq.0

        if (slope(1).ne.0.0d0.or.slope(2).ne.0.0d0.or.offset(1).ne.0.0d0.or.offset(2).ne.0.0d0) then
          print*,"*** Warning in undulator_source_analytic: Offset and slope not yet implemented"
        endif

        Com=twopi1/perlen*dom/om*dble(nord)

c       enor=  defl*om*echarge_CGS/(2.0d0*clight_CGS**2*gamma)*ajj ! Gl. (75) für dist=0
c        enor=-defl*om* echarge_CGS/(2.0d0*clight_CGS**2*gamma)*ajj * devlen/dist ! Gl. (77)
        enor=enorbase*om

        cthe=ci*om*dist/2.0d0/clight1

        do iy=1,ny
          if (nx.gt.1) then
            x=-pinw/2.0d0
          else
            x=0.0d0
          endif
          do ix=1,nx
            the2=(x**2+y**2)/dist**2
            wlc=devlen/2.0d0*(Com+om*the2/clight1/2.0d0) ! Gl. 164 und 175 im Anhang A
            ! Gl. 175: E(z0,the)
            if (wlc.ne.0.0d0) then
              esour(ix,iy)=enor*exp(cthe*the2)*sin(wlc)/wlc
c              write(66,*)ix,iy,x,y,wlc,sin(wlc)/wlc,dreal(esour(ix,iy)),dimag(esour(ix,iy))
            else
              esour(ix,iy)=enor
            endif
            x=x+dx
            if (abs(x).lt.eps) x=0.0d0
          enddo
          y=y+dy
          if (abs(y).lt.eps) y=0.0d0
        enddo

      endif !dist

      fdmax=curr/echarge1*alpha1*(defl*ajj*devlen/(wlharm*1.0d-9*2.0d0*gamma))**2    !Gl. (84)

      flux=curr/echarge1*pi1*alpha1*defl**2*ajj**2*dble(nper)/(2.0d0*(1.0d0+defl**2/2.0d0)) !Gl. (86)
      flux=flux*dble(nord) !Gl. (86), meine Korrektur: Mit Ordnung multipliziert

      fdmax=fdmax*banwid*1.0d-6 ! per mrad**2
      flux=flux*banwid

      end
*CMZ :  4.01/04 17/11/2023  11.59.33  by  Michael Scheer
*CMZ :  4.00/17 04/10/2022  08.10.22  by  Michael Scheer
*CMZ :  4.00/11 27/05/2021  09.41.25  by  Michael Scheer
*CMZ :  3.02/04 03/12/2014  15.11.16  by  Michael Scheer
*-- Author :    Michael Scheer   03/12/2014
      subroutine util_break
*KEEP,debugwave.
      include 'debugwave.cmn'
*KEND.
      return
      end
*CMZ :  4.00/07 06/04/2020  08.51.14  by  Michael Scheer
*CMZ :  3.06/00 18/02/2019  19.28.56  by  Michael Scheer
*CMZ :  2.70/02 14/12/2012  10.35.01  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE util_fold_function_gauss_lin(NF,XF,F,SIGMA,DNSIGMA,FG,WS1,WS2)
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

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG(X)=INT{F(XF)*G(XF-X),DXF}

C--   INPUT:

C-       NF:   NUMBER OF XF,F-VALUES
C-       XF:   ARRAY OF X-VALUES (MUST BE IN ASCENDING ORDER)
C-       F: ARRAY OF FUNCTION-VALUES
C-       SIGMA:  SIGMA OF GAUSSIAN
C-       DNSIGMA: NUMBER OF SIGMAS TO BE CONSIDERED

C--   OUTPUT:

C-       FG:   FG(X0) IS CALCULATED

      IMPLICIT NONE

      EXTERNAL FUNCTION DERF

      INTEGER NF,IH,IL,I

      REAL*8 XF(NF),F(NF),SIGMA,X0,FG(NF),DNSIGMA
      REAL*8 WS1(NF),WS2(NF)

      REAL*8 XL,XH,YL,YH,H,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,EPS

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

C- CHECK ASCENDING ORDER

      DO I=2,NF
        IF (XF(I).LE.XF(I-1))
     &    STOP '*** ERROR SR UTIL_FOLD_FUNCTION_GAUSS_LIN:
     &    ARRAY XF NOT IN ASCENDING ORDER ***'
      ENDDO

      EPS=(XF(NF)-XF(1))*1.0D-10

      DO IL=1,NF-1

        IH=IL+1

        XL=XF(IL)
        XH=XF(IH)
        YL=F(IL)
        YH=F(IH)

        XHXL=XH*XL
        XH2=XH*XH
        XL2=XL*XL

        H=XH-XL

        IF (H.LE.0.0D0) THEN
          PRINT*,
     &      '*** ERROR SR UTIL_FOLD_GAUSS_LIN:'
          PRINT*,
     &      '*** ARRAY XF NOT IN ASCENDING ORDER'
          STOP
        ENDIF

        WS1(IL)=(XH*YL-XL*YH)/h
        WS2(IL)=(YH-YL)/h

        FG(IL)=0.0D0

      ENDDO !NF-1

      FG(NF)=0.0D0

      SN=DNSIGMA*SIGMA
      S2=SIGMA*SIGMA
      S22=2.0D0*S2
      SNR21=1.0D0/(ROOT2*SIGMA)
      SR2PI1=R2PI1/SIGMA

      DO I=1,NF


        X0=XF(I)

        X02=X0*X0
        X02S2=S2+X02
        X023S2=S22+X02S2

        IF (X0-SN.GE.XF(1)-EPS.AND.X0+SN.LE.XF(NF)+EPS) THEN

C UPPER BRANCH

          DO IL=I,NF-1

            IH=IL+1

            XL=XF(IL)
            XH=XF(IH)

            IF (XL-X0.LE.SN) THEN

              IF (XH-X0.GT.SN) XH=X0+SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(-EXP(-DX**2/S22)*S2*WS2(IL)
     &          +SQPI2*SIGMA*(WS1(IL)+X0*WS2(IL))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(-EXP(-DX**2/S22)*S2*WS2(IL)
     &          +SQPI2*SIGMA*(WS1(IL)+X0*WS2(IL))*
     &          DERF(DX*SNR21))

              FG(I)=FG(I)+FGH-FGL

            ELSE
              GOTO 81
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IL

 81       CONTINUE

C LOWER BRANCH

          DO IH=I,2,-1

            IL=IH-1

            XL=XF(IL)
            XH=XF(IH)

            IF (X0-XH.LE.SN) THEN

              IF (X0-XL.GT.SN) XL=X0-SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)
     &          +X0*(WS2(IL)))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)
     &          +X0*(WS2(IL)))*
     &          DERF(DX*SNR21))

              FG(I)=FG(I)+FGH-FGL

            ELSE
              GOTO 82
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IH

 82       CONTINUE

        ELSE IF (X0+SN.GT.XF(NF)) THEN

          GOTO 88

        ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

      ENDDO !NF

 88   CONTINUE

      RETURN
      END
*CMZ :  3.06/00 18/02/2019  18.55.48  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  12.31.20  by  Michael Scheer
*CMZ :  2.63/03 02/05/2008  14.41.00  by  Michael Scheer
*CMZ :  2.52/04 12/07/2004  16.15.58  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.24.55  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_FOLD_FUNCTION_GAUSS(NF,XF,F,SIGMA,RNSIGMA,FG,
     &  COEF,WS1,WS2,WS3,WS4)
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

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG(X)=INT{F(XF)*G(XF-X),DXF}

C--   INPUT:

C-       NF:   NUMBER OF XF,F-VALUES
C-       XF:   ARRAY OF X-VALUES (MUST BE IN ASCENDING ORDER)
C-       F: ARRAY OF FUNCTION-VALUES
C-       SIGMA:  SIGMA OF GAUSSIAN
C-       RNSIGMA: NUMBER OF SIGMAS TO BE CONSIDERED

C--   OUTPUT:

C-       FG:   FG(X0) IS CALCULATED

      IMPLICIT NONE

      EXTERNAL FUNCTION DERF

      INTEGER NF,IH,IL,I

      REAL*8 XF(NF),F(NF),SIGMA,X0,FG(NF),RNSIGMA
      REAL*8 COEF(NF)
      REAL*8 WS1(NF),WS2(NF),WS3(NF),WS4(NF)

      REAL*8 CH,CL,CH2,CL2,CHCL,CH2CL,CHCL2,XL,XH,YL,YH,H,H61,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,EPS

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

C- CHECK ASCENDING ORDER

      DO I=2,NF
        IF (XF(I).LE.XF(I-1))
     &    STOP '*** ERROR SR UTIL_FOLD_FUNCTION_GAUSS:
     &    ARRAY XF NOT IN ASCENDING ORDER ***'
      ENDDO

      EPS=(XF(NF)-XF(1))*1.0D-10

C- SPLINES OF FUNCTION F

      CALL UTIL_SPLINE_COEF(XF,F,NF,-9999.0d0,-9999.0d0,COEF,WS1,WS2,WS3,WS4)

      DO IL=1,NF-1

        IH=IL+1

        XL=XF(IL)
        XH=XF(IH)
        YL=F(IL)
        YH=F(IH)
        CL=COEF(IL)
        CH=COEF(IH)

        CL2=2.0D0*CL
        CH2=2.0D0*CH
        CHCL=CH-CL
        CHCL2=CH+CL2
        CH2CL=CH2+CL

        XHXL=XH*XL
        XH2=XH*XH
        XL2=XL*XL

        H=XH-XL

        IF (H.LE.0.0D0) THEN
          PRINT*,
     &      '*** ERROR SR UTIL_FOLD_GAUSS:'
          PRINT*,
     &      '*** ARRAY XF NOT IN ASCENDING ORDER'
          STOP
        ENDIF

        H61=1.0D0/(6.0D0*H)

        WS1(IL)=((CHCL2*XH-CH2CL*XL)*XHXL+6.0D0*(XH*YL-XL*YH))*H61
        WS2(IL)=((CH2-CL2)*XHXL-CHCL2*XH2+CH2CL*XL2+6.0D0*(YH-YL))*H61
        WS3(IL)=(-CH*XL+CL*XH)/(2.0D0*H)
        WS4(IL)=CHCL*H61

        FG(IL)=0.0D0

      ENDDO !NF-1

      FG(NF)=0.0D0

      SN=RNSIGMA*SIGMA
      S2=SIGMA*SIGMA
      S22=2.0D0*S2
      SNR21=1.0D0/(ROOT2*SIGMA)
      SR2PI1=R2PI1/SIGMA

      DO I=1,NF

        X0=XF(I)

        X02=X0*X0
        X02S2=S2+X02
        X023S2=S22+X02S2

        IF (X0-SN.GE.XF(1)-EPS.AND.X0+SN.LE.XF(NF)+EPS) THEN

C UPPER BRANCH

          DO IL=I,NF-1

            IH=IL+1

            XL=XF(IL)
            XH=XF(IH)

            IF (XL-X0.LE.SN) THEN

              IF (XH-X0.GT.SN) XH=X0+SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XH+X0)+WS4(IL)*(S22+XH**2+XH*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XL+X0)+WS4(IL)*(S22+XL**2+XL*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              FG(I)=FG(I)+FGH-FGL

            ELSE
              GOTO 81
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IL

 81       CONTINUE

C LOWER BRANCH

          DO IH=I,2,-1

            IL=IH-1

            XL=XF(IL)
            XH=XF(IH)

            IF (X0-XH.LE.SN) THEN

              IF (X0-XL.GT.SN) XL=X0-SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XH+X0)+WS4(IL)*(S22+XH**2+XH*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XL+X0)+WS4(IL)*(S22+XL**2+XL*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              FG(I)=FG(I)+FGH-FGL

            ELSE
              GOTO 82
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IH

 82       CONTINUE

        ELSE IF (X0+SN.GT.XF(NF)) THEN

          GOTO 88

        ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

      ENDDO !NF

 88   CONTINUE

      RETURN
      END
*CMZ :  4.01/07 16/08/2024  09.21.28  by  Michael Scheer
*-- Author :    Michael Scheer   15/08/2024
*CMZ :          15/08/2024  11.02.29  by  Michael Scheer
      subroutine util_fold_gauss_lin_2d(nx,ny,x,y,fin,rnsigx,sigx,rnsigy,sigy,fold)

      implicit none

      integer ix,iy,nx,ny

      real*8 fin(nx,ny),x(nx),y(ny),sigx,sigy,fold(nx,ny),rnsigx,rnsigy,
     &  fxf(nx,ny),f(max(nx,ny)),fg(max(nx,ny)),ws1(max(nx,ny)),ws2(max(nx,ny))

      do iy=1,ny
        f(1:nx)=fin(1:nx,iy)
        if (sigx.gt.0.0d0) then
          call util_fold_function_gauss_lin(nx,x,f,sigx,rnsigx,fg,ws1,ws2)
          fxf(1:nx,iy)=fg(1:nx)
        else
          fxf(1:nx,iy)=fin(1:nx,iy)
        endif
      enddo

      do ix=1,nx
        if (sigy.gt.0.0d0) then
          f(1:ny)=fxf(ix,1:ny)
          call util_fold_function_gauss_lin(ny,y,f,sigy,rnsigy,fg,ws1,ws2)
          fold(ix,1:ny)=fg(1:ny)
        else
          fold(ix,1:ny)=fxf(ix,1:ny)
        endif
      enddo

      return
      end
*CMZ :  4.02/00 06/03/2025  18.47.11  by  Michael Scheer
*CMZ :  4.01/07 23/08/2024  14.52.21  by  Michael Scheer
*CMZ : 00.00/15 07/12/2012  20.04.19  by  Michael Scheer
*-- Author :    Michael Scheer   06/12/2012
      subroutine util_fold_gauss_2d(nx,ny,x,y,f,sigx,rnsigx,sigy,rnsigy,fg,ispline,istat)

c Folding of f(x(ix),y(iy)) with a 2D Gaussian.
c The Gaussian is considered from -rnsig*sig -> +rnsig*sig
C IT'S ONLY CORRECT FOR X,Y FAR ENOUGH FROM THE EDGES!!

c Dimensions f(nx,ny), fg(nx,ny)

      implicit none

      double precision, dimension(:), allocatable :: wf,wfg,w1,w2,w3,w4,coef

      double precision
     &  sigx,rnsigx,sigy,rnsigy,x(nx),y(ny),f(nx,ny),fg(nx,ny)

      integer nx,ny,istat,ix,iy,ispline
      integer :: nallox=0,nalloy=0

      save

      if (2.0d0*rnsigx*sigx.ge.x(nx)-x(1).or.2.0d0*rnsigy*sigy.ge.y(ny)-y(1)) then
        istat=-1
        fg=0.0d0
        return
      endif

      if (ispline.eq.0) then
        call util_fold_gauss_lin_2d(nx,ny,x,y,f,rnsigx,sigx,rnsigy,sigy,fg)
        istat=0
        return
      endif

      if (nx.gt.0.and.istat.lt.0) deallocate(wf,wfg,w1,w2,w3,w4,coef)

      istat=0
      fg=0.0d0

      if (nx.lt.3.or.ny.lt.3) then
        istat=-1
        return
      endif

      if (nx.gt.nallox.or.ny.gt.nalloy) then
        if (nx.eq.0) deallocate(wf,wfg,w1,w2,w3,w4,coef)
        allocate(wf(max(nx,ny)))
        allocate(wfg(max(nx,ny)))
        allocate(w1(max(nx,ny)))
        allocate(w2(max(nx,ny)))
        allocate(w3(max(nx,ny)))
        allocate(w4(max(nx,ny)))
        allocate(coef(max(nx,ny)))
        nallox=nx
        nalloy=ny
      endif

      do iy=1,ny
        wf=f(1:nx,iy)
        call util_fold_function_gauss(nx,x,wf,sigx,rnsigx,wfg,coef,w1,w2,w3,w4)
        fg(1:nx,iy)=wfg(1:nx)
      enddo !iy

      do ix=1,nx
        wf=fg(ix,1:ny)
        call util_fold_function_gauss(ny,y,wf,sigy,rnsigy,wfg,coef,w1,w2,w3,w4)
        fg(ix,1:ny)=wfg(1:ny)
      enddo !iy

      return
      end
*CMZ :  4.01/05 05/01/2024  11.27.46  by  Michael Scheer
*CMZ :  4.01/04 17/12/2023  11.45.19  by  Michael Scheer
*CMZ :  4.01/02 08/05/2023  13.06.52  by  Michael Scheer
*CMZ :  4.01/00 10/02/2023  13.27.16  by  Michael Scheer
*CMZ :  4.00/15 28/04/2022  15.32.20  by  Michael Scheer
*CMZ :  4.00/13 16/11/2021  17.32.24  by  Michael Scheer
*CMZ :  4.00/09 15/08/2020  08.51.05  by  Michael Scheer
*CMZ :  3.05/05 10/07/2018  09.19.31  by  Michael Scheer
*CMZ :  3.05/04 05/07/2018  11.10.09  by  Michael Scheer
*CMZ :  3.05/00 27/04/2018  15.22.16  by  Michael Scheer
*CMZ :  3.03/04 13/10/2017  09.16.28  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.32.35  by  Michael Scheer
*CMZ :  3.02/04 13/03/2015  10.38.25  by  Michael Scheer
*CMZ :  2.69/02 02/11/2012  16.40.18  by  Michael Scheer
*CMZ :  2.68/05 04/09/2012  13.30.58  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  11.52.24  by  Michael Scheer
*CMZ :  2.68/03 31/08/2012  09.45.28  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_omp(icharge,
     &  gammai,dgamtot,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efxn,efyn,efzn,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &  nthstep,nstep,ndim,traxyz,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  nphener,phener,aradx,arady,aradz,stokes,powden,
     &  ieneloss,ivelofield
     &  ,istatus,ith,modewave)
c123456789123456789_123456789_123456789_123456789_123456789_123456789_12
c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

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

c This subroutine calculates the trajectory and the synchrotron radiation
c of an electron passing a magnetic and electric field. The electric field
c part is very preliminary and not yet fully tested. The FORTRAN random
c generator is used in uradrndm. It should be replaced by a better one
c provided by the user.

c The fields B=(bx,by,bz) and E=(ex,ey,ez) are calculated in the routine
c uradfield_omp(x,y,z,bx,by,bz,ex,ey,ez,gamma,istatus) provided by the user.
c As an example uradbmap.f may be used, which reads a 3D map of the mag. field

c Coordinate system (right handed):
c -----------------------------------

c      x: longitudinal direction
c      y: transversal vertical direction
c      z: transversal horizontal direction

c Units:
c ------

c SI-units: m, Tesla, sec, V etc., but eV for the photon energy
c The flux-density is given in Nph/s/mm**2/0.1%BW, the power-density in W/mm**2

c Input:
c --------

c integer icharge: Particle charge ( +/- 1)
c real*8 gammai: Gamma factor of the e-

c real*8 xelec:  Initial x of e-
c real*8 yelec:  Initial y of e-
c real*8 zelec:  Initial z of e-

c real*8 vxelec:  Initial velocity in x of e-
c real*8 vyelec:  Initial velocity in y of e-
c real*8 vzelec:  Initial velocity in z of e-
c The velocity is internally normalized, so the input norm does not matter

c real*8 xf: x of point in exit plane
c real*8 yf: y of point in exit plane
c real*8 zf: z of point in exit plane

c real*8 [efxn,efyn,efzn]: Normal vector of exit plane

c real*8 vnxex: x component of normal vector of exit plane
c real*8 vnyex: y component of normal vector of exit plane
c real*8 vnzex: z component of normal vector of exit plane

c real*8 ds : step size for tracking

c The tracking stops, when the electron hits the exit plane. The size of the last
c step is corrected such that the plane is hit.

c integer nthstep: If nstep > 0, the trajectory array traxyz is filled,
c                  see below
c integer ndim: Dimension of traxyz, see below

c real*8 phelow: Lowest photon energy / eV
c real*8 phehig: Higest photon energy / eV

c integer nphener: Number of equidistant photon energies

c integer ieneloss:  0: no energy loss due to synchotron radiation
c                    1: continous energy loss due to synchotron radiation
c integer ieneloss: -1: no energy loss due to synchotron radiation with quantum
c                       fluctuations

c integer ivelofield: Contral flag for the calculation of the velocity field:
c                    0: the spectrum includes the velocity field
c                    1: the specrum does not include the velocity field
c                    2: the spectrum includes only the velocity field

c Output:
c -------

c integer: istatus: Status flag:
c  0: no error found
c -1: initial gamma or velocity zero
c -2: dimension ndim of traxyz exceeded
c -3: bad value of ivelofield
c  else: status of uradfield_omp

c real*8 xexit: x of last point of the trajectory
c real*8 yexit: y of last point of the trajectory
c real*8 zexit: z of last point of the trajectory
c real*8 texit: t of last point of the trajectory

c real*8 vnxex: x component of norm. velocity vector of last point
c real*8 vnyex: y component of norm. velocity vector of last point
c real*8 vnzex: z component of norm. velocity vector of last point

c real*8 phener(nphener): Array of equidistant photon energies

c integer nstep: Number of tracking steps done, i.e. used length of traxyz

c real*8 traxyz(1:14,i): Array:
c        traxyz(1,i):  x
c        traxyz(2,i):  y
c        traxyz(3,i):  z
c        traxyz(4,i):  tracking time
c        traxyz(5,i):  x-comp. of norm. velocity vector
c        traxyz(6,i):  y-comp. of norm. velocity vector
c        traxyz(7,i):  z-comp. of norm. velocity vector
c        traxyz(8,i):  x-comp. of mag. field in the center of the step
c        traxyz(9,i):  y-comp. of mag. field in the center of the step
c        traxyz(10,i): z-comp. of mag. field in the center of the step
c        traxyz(11,i): gamma
c        traxyz(12,i): x-comp. of elec. field in the center of the step
c        traxyz(13,i): y-comp. of elec. field in the center of the step
c        traxyz(14,i): z-comp. of elec. field in the center of the step

c The phase is calculated by phase=phase0+n*dt*dphase,
c where dphase is the phase difference of the nth step dt. Phase0=
c (xobsv-xelec)/clight. The phase factor of the integrand is
c exp(i*omega*phase),where omega referes to the considered photon energy.
c complex*16 aradx(i): x-comp. of amplitude of radiation field of phener(i)
c complex*16 arady(i): y-comp. of amplitude of radiation field of phener(i)
c complex*16 aradz(i): z-comp. of amplitude of radiation field of phener(i)

c real*8 array of Stokes parameters of ith photon energy:
c        stokes(1,i): S0, i.e. total intensity
c        stokes(2,i): S1, Stokes parameter of linear +/- 90 degree polarisation
c        stokes(3,i): S2, Stokes parameter of linear +/- 45 degree polarisation
c        stokes(4,i): S3, Stokes parameter of circular polarisation
c        S0 = sqrt(S1**2+S2**2+S3**2)

c Compilation:
c ------------

c For uradbmap at least F90 is required.
c The line length exceeds 72 characters, please use an appropriate
c compiler option. It is recommended to use compiler options to initialize all
c variables to zero and to treat them as ,,saved''

      implicit none

      complex*16 aradx(nphener),arady(nphener),aradz(nphener),
     &  ziom,zi,zidom,zone,ziomr1,zicr1,zic,
     &  expom1,expom,dexpomph1,dexpomph,ddexpomph,dexpom,
     &  expomv2,vstokes(4,3),
     &  apolh,apolr,apoll,apol45,dum3

      double precision
     &  gammai,dgamtot,dt2,powden,t,phase,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  phener(nphener),dom2,c,rspn
     &  ,traxyz(14,ndim),stokes(4,nphener),x1,y1,z1,vx1,vy1,
     &  vz1,x2,y2,z2,vx2,vy2,vz2
     &  ,ds,dtim,bshift,x2b,y2b,z2b,bx1,by1,bz1,bx2,by2,bz2
     &  ,dgamsum,gamma,dt
     &  ,x2int,y2int,z2int,ddt,dddt,ddt2,dddt2
     &  ,vx2int,vy2int,vz2int,vxpint,vypint,vzpint
     &  ,vxp,vyp,vzp
     &  ,x3int,y3int,z3int,vx3int,vy3int,vz3int,ddddt,ddddt2
     &  ,efx,efy,efz,xf,yf,zf,dist1,dist2,disti
     &  ,dtim0,beta,vn,efx2,efy2,efz2,t1,t2,clight,c1,
     &  dgamma,vxsign,bx,by,bz,bpx,bpy,bpz,rarg(5),px,py,pz,
     &  dphase,r,rx,ry,rz
     &  ,dom1,rnbx,rnby,rnbz,dum11,rnr4,br4,b3,rnr2,br2,bet1n,
     &  rnx,rny,rnz,r1,banwid,specnor,pownor,current,
     &  stok1,stok2,stok3,stok4,om,dom,hbarev,echarge,eps0,pi,vsto,dph,
     &  r0,efxn,efyn,efzn

      integer ieneloss,istatus,icharge,nphener,ivelofield,
     &  nthstep,izaehl,nstep,ndim,kstep,lstep,ifreq,isto,ifail,ith
      integer :: kcount=1,modewave

c      integer,save :: ical=0
c+seq,uservar.


      data bshift/0.5d0/
      data clight/2.99792458d8/
      data hbarev/6.58211889D-16/
      data banwid/1.0d-3/
      data current/0.10d0/ ! Be care, so also ucur in URADBUE etc.
      data eps0/8.854187817D-12/
      data echarge/1.602176462D-19/
      data pi/3.14159265358979d0/

      data zi/(0.0d0,1.0d0)/
      data zone/(1.0d0,0.0d0)/

      dph=0.0d0
c      ical=ical+1

      if (nphener.gt.0) phener(1)=phelow
      if (nphener.gt.1) dph=(phehig-phelow)/(nphener-1)

      do ifreq=2,nphener
        phener(ifreq)=phener(ifreq-1)+dph
      enddo

      istatus=0
      ifail=0
      if (icharge.le.0) icharge=-1
      if (icharge.gt.0) icharge=1

      vn=norm2([efxn,efyn,efzn])
      efx=efxn/vn
      efy=efyn/vn
      efz=efzn/vn

      x1=xelec
      y1=yelec
      z1=zelec
      vx1=vxelec
      vy1=vyelec
      vz1=vzelec
      t1=0.0d0

      gamma=gammai
      beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))
      vn=sqrt(vx1*vx1+vy1*vy1+vz1*vz1)
      vx1=vx1/vn*clight*beta
      vy1=vy1/vn*clight*beta
      vz1=vz1/vn*clight*beta
      vn=beta*clight

c vxsign takes care for the direction of flight, since particle must gain
c energy if tracked back

      if (vx1.lt.0) then
        vxsign=-1.0d0
      else
        vxsign=1.0d0
      endif

      dgamsum=0.0d0
      dgamtot=0.0d0
      powden=0.0d0
      aradx=(0.0d0,0.0d0)
      arady=(0.0d0,0.0d0)
      aradz=(0.0d0,0.0d0)

      dtim=ds/vn
      dt=dtim
      dt2=dtim*bshift
      dtim0=dtim

      x2=x1
      y2=y1
      z2=z1
      t2=t1

      vx2=vx1
      vy2=vy1
      vz2=vz1

      x2b=x1+vx1*dt2
      y2b=y1+vy1*dt2
      z2b=z1+vz1*dt2

      call uradfield_omp(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,gamma,istatus,
     &  modewave)
      if (istatus.ne.0) ifail=ifail+abs(istatus)
      istatus=0

      nstep=0
      vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)

      if (gamma.le.0.0d0.or.vn.le.0.0d0) then
        istatus=-1
        return
      endif

      kstep=-1
      nstep=0

      if(nthstep.gt.0) then

        nstep=nstep+1
        if (nstep.gt.ndim) then
          nstep=nstep-1
          istatus=-2
          goto 9000
        endif

        kstep=kstep+1
        if (kstep.eq.nthstep) kstep=0

        traxyz(1,nstep)=x2
        traxyz(2,nstep)=y2
        traxyz(3,nstep)=z2
        traxyz(4,nstep)=t2
        traxyz(5,nstep)=vx2/vn
        traxyz(6,nstep)=vy2/vn
        traxyz(7,nstep)=vz2/vn
        traxyz(8,nstep)=bx2
        traxyz(9,nstep)=by2
        traxyz(10,nstep)=bz2
        traxyz(11,nstep)=gamma
        traxyz(12,nstep)=efx2
        traxyz(13,nstep)=efy2
        traxyz(14,nstep)=efz2

      endif !nthstep.gt.0

      dom=0.0d0
      om=0.0d0
      if (nphener.gt.1) then
        om=phener(1)/hbarev
        dom=(phener(2)-phener(1))/hbarev
      else if (nphener.eq.1) then
        om=phener(1)/hbarev
      endif

      c=clight
      c1=1.0d0/clight

      zidom=zi*dom
      ziom=zi*om
      zic=zi*c

      lstep=0
      t=-dt
      r0=xobsv-xelec
      r=sqrt((xobsv-x1)**2+((yobsv-y1)**2+(zobsv-z1)**2))
      PHASE=(r-r0)*c1
      expom1=zone
      dexpomph1=zone

c--- Loop der Trajektorie

      izaehl=0
1000  continue

      izaehl=izaehl+1
c      print*,ith,izaehl,x2
      if (x2.ne.x2) then
        istatus=-99
        return
      endif

      if (lstep.eq.1) then
        dtim=abs(dist2)/vn
        dt=dtim
        dt2=dtim/2.0d0
      endif

      x1=x2
      y1=y2
      z1=z2

      t1=t2

      vx1=vx2
      vy1=vy2
      vz1=vz2

      bx1=bx2
      by1=by2
      bz1=bz2

      x2b=x1+vx1*dt2
      y2b=y1+vy1*dt2
      z2b=z1+vz1*dt2

      call uradfield_omp(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,gamma,istatus,
     &  modewave)
      if (istatus.ne.0) ifail=ifail+abs(istatus)
      istatus=0

      call uradstep_omp(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,efx2,efy2,efz2,
     &  dtim,x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,ieneloss,
     &  dgamma)

      if (ieneloss.ne.0) then
        dgamsum=dgamsum+vxsign*dgamma
        if (abs(dgamsum).gt.gamma*1.0d-8) then
          gamma=gamma+dgamsum
          dgamtot=dgamtot+dgamsum
          dgamsum=0.0d0
        endif
        beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))
        vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
        vx2=vx2/vn*clight*beta
        vy2=vy2/vn*clight*beta
        vz2=vz2/vn*clight*beta
      endif

      BX=VX2*C1
      BY=VY2*C1
      BZ=VZ2*C1

      BPX=VXP*C1
      BPY=VYP*C1
      BPZ=VZP*C1

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

      RX=XOBSV-X2
      RY=YOBSV-Y2
      RZ=ZOBSV-Z2

      R=SQRT(RX*RX+RY*RY+RZ*RZ)
      R1=1.D0/R
      ZICR1=ZIC*R1

      RNX=RX*R1
      RNY=RY*R1
      RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

      BET1N=(1.0D0-BX*RNX)-BY*RNY-BZ*RNZ

      br2=by**2+bz**2
      rnr2=rny**2+rnz**2
      b3=beta**3
      br4=br2**2
      rnr4=rnr2**2

      if(br2.lt.1.0d-4.and.rnr2.lt.1.0d-4) then
        bet1n=
     &    1.0d0/(1.0d0+beta)/gamma**2
     &    +beta*(rnr2/2.0d0
     &    +rnr4/8.0d0)
     &    +(br2/2.0d0
     &    -br2*rnr2/4.0d0
     &    -br2*rnr4/16.0d0)/beta
     &    +b3*br4*(1.0d0/8.0d0
     &    -rnr2/16.0d0
     &    -rnr4/64.0d0)
     &    -by*rny
     &    -bz*rnz
      endif

      DUM11=1.D0/BET1N
      DOM1=1.D0/(R*BET1N*BET1N)

      RNBX=RNX-BX
      RNBY=RNY-BY
      RNBZ=RNZ-BZ

      PX=(RNBY*BPZ-RNBZ*BPY)
      PY=(RNBZ*BPX-RNBX*BPZ)
      PZ=(RNBX*BPY-RNBY*BPX)

      IF (IVELOFIELD.EQ.0.OR.IVELOFIELD.EQ.2) THEN !2 WEGEN POWER
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
        RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
      ELSE IF (IVELOFIELD.EQ.1) THEN
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
        RARG(3)=(RNX*PY-RNY*PX)*DOM1
      ELSE IF (IVELOFIELD.LT.0) THEN
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNX-BX)*DOM2
        RARG(2)=(RNY-BY)*DOM2
        RARG(3)=(RNZ-BZ)*DOM2
      ELSE  !IVELOFIELD
        istatus=-3
        return
      ENDIF !IVELOFIELD

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS      RARG(4)=T+R*C1

      DPHASE=BET1N*DT

      RARG(5)=(RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11

      powden=powden+rarg(5)*dt

C REAL PART OF INTEGRAND }

C COMPLEX PART OF INTEGRAND {

C    ASSUMES phener(I+1)=2*phener(I)   FOR IFREQ2P=2
C    OR phener(I+1)=phener(I)+DELTA    FOR IFREQ2P>2

C--- LOOP OVER ALL FREQUENCES

      if (nphener.gt.0) then

        IFREQ=1

        OM=phener(IFREQ)/hbarev
        ZIOM=ZI*OM

        EXPOM=EXPOM1
        DEXPOMPH1=EXP(ZIOM*DPHASE)
        DEXPOMPH=DEXPOMPH1

        IF(nphener.GT.1) THEN
          DEXPOM=EXP(ZIDOM*PHASE)
          DDEXPOMPH=EXP(ZIDOM*DPHASE)
        ENDIF  !IFREQ2P

        IF (IVELOFIELD.NE.2) THEN

          dum3=expom*(zone-dexpomph)/om/bet1n

          aradx(ifreq)=aradx(ifreq)+rarg(1)*dum3
          arady(ifreq)=arady(ifreq)+rarg(2)*dum3
          aradz(ifreq)=aradz(ifreq)+rarg(3)*dum3

        ELSE !IVELOFIELD

          EXPOMV2=R1/BET1N*EXPOM*(ZONE-DEXPOMPH)
          ZIOMR1=ZONE+ZICR1/OM

          aradx(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
          arady(ifreq)=arady(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
          aradz(ifreq)=aradz(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)

        ENDIF !IVELOFIELD

        IF (IVELOFIELD.NE.2) THEN

          DO IFREQ=2,nphener

            OM=OM+DOM
            EXPOM=EXPOM*DEXPOM
            DEXPOMPH=DEXPOMPH*DDEXPOMPH

            EXPOMV2=1.0D0/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)

            aradx(ifreq)=aradx(ifreq)+RARG(1)*EXPOMV2
            arady(ifreq)=arady(ifreq)+RARG(2)*EXPOMV2
            aradz(ifreq)=aradz(ifreq)+RARG(3)*EXPOMV2

          ENDDO   !LOOP OVER ALL FREQUENCES

        else

          DO IFREQ=2,nphener

            OM=OM+DOM
            EXPOM=EXPOM*DEXPOM
            DEXPOMPH=DEXPOMPH*DDEXPOMPH

            EXPOMV2=R1/BET1N*EXPOM*(ZONE-DEXPOMPH)
            ZIOMR1=ZONE+ZICR1/OM

            aradx(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
            arady(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
            aradz(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)

          ENDDO   !LOOP OVER ALL FREQUENCES

        ENDIF !IVELOFIELD


C COMPLEX PART OF INTEGRAND }

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

        PHASE=PHASE+DPHASE
        EXPOM1=EXPOM1*DEXPOMPH1

      endif !(nphener.gt.0) then

      t2=t1+dtim

c ef is normal vector of perpendiculare plane at the end of the reference orbit
c dist is distance of electron to this plane
c tracking stops if trajectory hits this plane

      dist2=(x2-xf)*efx+(y2-yf)*efy+(z2-zf)*efz

      if (lstep.eq.0.and.dist2.lt.0.0d0.and.dist2.gt.-2.0d0*ds)  then

        lstep=1
        goto 1000

      else

        kstep=kstep+1
        if (kstep.eq.nthstep) kstep=0

        if(kstep.eq.0) then

          nstep=nstep+1

          if (nstep.gt.ndim) then
            nstep=nstep-1
            istatus=-2
            goto 9000
          endif

          traxyz(1,nstep)=x2
          traxyz(2,nstep)=y2
          traxyz(3,nstep)=z2
          traxyz(4,nstep)=t2
          traxyz(5,nstep)=vx2/vn
          traxyz(6,nstep)=vy2/vn
          traxyz(7,nstep)=vz2/vn
          traxyz(8,nstep)=bx2
          traxyz(9,nstep)=by2
          traxyz(10,nstep)=bz2
          traxyz(11,nstep)=gamma
          traxyz(12,nstep)=efx2
          traxyz(13,nstep)=efy2
          traxyz(14,nstep)=efz2

        endif

        if (lstep.eq.1) goto 9000
        goto 1000

      endif !lstep and dist2

9000  continue

      xexit=x2
      yexit=y2
      zexit=z2

      vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
      vnxex=vx2/vn
      vnyex=vy2/vn
      vnzex=vz2/vn

      texit=t2

      specnor=
     &  banwid
     & /(4.0d0*pi**2*clight*hbarev)
     & /(4.0d0*pi*eps0)
     &  *current/1.0d6 !per mm**2

      pownor=echarge/16.0d0/pi/pi/eps0/clight*current/1.0d6 !W/mm**2

      rspn=sqrt(specnor)

      vstokes(1,1)=( 0.0d0,        0.0d0)      !horizontal polarization
      vstokes(1,2)=( 0.0d0,        0.0d0)
      vstokes(1,3)=(-1.0d0,       -1.0d0)

      vstokes(2,1)=( 0.0d0,        0.0d0)      !right handed polarization
      vstokes(2,2)=( 0.0d0,       -1.0d0)
      vstokes(2,3)=(+1.0d0,        0.0d0)

      vstokes(3,1)=( 0.0d0,        0.0d0)      !left handed polarization
      vstokes(3,2)=( 0.0d0,       -1.0d0)
      vstokes(3,3)=(-1.0d0,        0.0d0)

      vstokes(4,1)=( 0.0d0,        0.0d0)      !45 degree linear polarization
      vstokes(4,2)=( 1.0d0,        0.0d0)
      vstokes(4,3)=(-1.0d0,        0.0d0)

      do isto=1,4
        vsto=dsqrt
     &    (cdabs(vstokes(isto,1))**2
     &    +cdabs(vstokes(isto,2))**2
     &    +cdabs(vstokes(isto,3))**2)
        vstokes(isto,1)=vstokes(isto,1)/vsto
        vstokes(isto,2)=vstokes(isto,2)/vsto
        vstokes(isto,3)=vstokes(isto,3)/vsto

      enddo

      do ifreq=1,nphener

        aradx(ifreq)=aradx(ifreq)*rspn
        arady(ifreq)=arady(ifreq)*rspn
        aradz(ifreq)=aradz(ifreq)*rspn

        apolh=
     &    aradx(ifreq)*conjg(vstokes(1,1))
     &    +arady(ifreq)*conjg(vstokes(1,2))
     &    +aradz(ifreq)*conjg(vstokes(1,3))

        apolr=
     &    aradx(ifreq)*conjg(vstokes(2,1))
     &    +arady(ifreq)*conjg(vstokes(2,2))
     &    +aradz(ifreq)*conjg(vstokes(2,3))

        apoll=
     &    aradx(ifreq)*conjg(vstokes(3,1))
     &    +arady(ifreq)*conjg(vstokes(3,2))
     &    +aradz(ifreq)*conjg(vstokes(3,3))

        apol45=
     &    aradx(ifreq)*conjg(vstokes(4,1))
     &    +arady(ifreq)*conjg(vstokes(4,2))
     &    +aradz(ifreq)*conjg(vstokes(4,3))

        stok1=
     &    apolr*conjg(apolr)+
     &    apoll*conjg(apoll)

        stok2=-stok1+
     &    2.*apolh*conjg(apolh)

        stok3=
     &    2.*apol45*conjg(apol45)-
     &    stok1

        stok4=
     &    apolr*conjg(apolr)-
     &    apoll*conjg(apoll)

        stokes(1,ifreq)=stok1
        stokes(2,ifreq)=stok2
        stokes(3,ifreq)=stok3
        stokes(4,ifreq)=stok4

      enddo !nphener

      powden=powden*pownor

      if (istatus.ge.0.and.ifail.ne.0) istatus=ifail

      return
      end
*CMZ :  4.01/04 25/11/2023  13.39.02  by  Michael Scheer
*CMZ :  4.01/02 09/05/2023  13.15.30  by  Michael Scheer
*CMZ :  4.00/15 28/04/2022  11.45.17  by  Michael Scheer
*CMZ :  4.00/13 16/11/2021  17.18.53  by  Michael Scheer
*CMZ :  3.05/05 09/07/2018  15.22.23  by  Michael Scheer
*CMZ :  3.05/04 05/07/2018  08.56.42  by  Michael Scheer
*CMZ :  3.04/00 23/01/2018  17.17.28  by  Michael Scheer
*CMZ :  3.03/04 04/12/2017  15.56.53  by  Michael Scheer
*CMZ :  3.03/02 18/11/2015  12.56.22  by  Michael Scheer
*CMZ :  3.02/04 13/03/2015  10.36.11  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/04 04/09/2012  09.38.42  by  Michael Scheer
*CMZ :  2.68/03 29/08/2012  12.25.27  by  Michael Scheer
*-- Author :    Michael Scheer   22/08/2012
      subroutine uradfield_omp(x,y,z,bxout,byout,bzout,ex,ey,ez,gamma,istatus,
     &  modewave)

      implicit none

*KEEP,ampli.
      include 'ampli.cmn'
*KEND.

      double precision :: x,y,z,bx,by,bz,ex,ey,ez,
     &  bxout,byout,bzout,gamma,axout,ayout,azout

      integer :: istatus,modewave

      ex=0.0d0
      ey=0.0d0
      ez=0.0d0

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0

      if (phrb0v.ne.0.0d0) then

        call bhalba_omp(phrb0v,phrperl,x+phrshift/2.0d0,y,z,bx,by,bz)

        bxout=bxout+bx
        byout=byout+by
        bzout=bzout+bz

      endif

      if (phrb0h.ne.0.0d0) then

        call bhalba_omp(phrb0h,phrperl,x-phrshift/2.0d0,y,z,bx,by,bz)

        bxout=bxout+bx
        byout=byout+bz
        bzout=bzout-by

      endif

      istatus=0

      return
      end
*CMZ :  4.00/15 27/04/2022  11.54.41  by  Michael Scheer
*CMZ :  3.05/20 31/10/2018  10.51.37  by  Michael Scheer
*CMZ :  3.05/04 04/07/2018  14.02.25  by  Michael Scheer
*CMZ :  3.05/00 02/05/2018  11.33.05  by  Michael Scheer
*CMZ :  3.04/00 01/03/2018  15.00.40  by  Michael Scheer
*CMZ :  3.02/04 11/03/2015  13.19.06  by  Michael Scheer
*CMZ :  3.01/04 20/05/2014  12.49.15  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/05 06/09/2012  15.53.22  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  09.07.11  by  Michael Scheer
*CMZ :  2.68/03 27/08/2012  09.35.31  by  Michael Scheer
*-- Author :    Michael Scheer
c*******************************************************************************
      subroutine uradstep_omp(xin,yin,zin,vx1,vy1,vz1,bxin,byin,bzin,
     &  efieldx,efieldy,efieldz,
     &  dtim,
     &  x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,
     &  ieneloss,dgamma)
c*******************************************************************************

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

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

      implicit none

      double precision bbet1,vsenk1,vsenk2,tz,tz2,tz3,tz4,tz5,
     &  vsenkzyk,acc,
     &  efieldx,efieldy,efieldz,
     &  dk,z12,y12,vz12,vy12,vx12,xin,yin,zin

      double precision bx2,by2,bz2,bsq,bbet,bux,buy,buz,v0sq,
     &  vx1,vy1,vz1,
     &  v0bet,vpar,vparx,vpary,vparz,vparsq,vsenk,dtim,bxin,byin,bzin
      double precision x1n,y1n,z1n,x2n,y2n,z2n

      double precision x2,y2,z2,vx2,vy2,vz2,x1,y1,z1,vxp,vyp,vzp,zyk,
     &  gamma,sz,cz,zp,yp,dy,dz,brho
     &  ,dgamma

      double precision bmovecut,dgammao,v0,v12n(3),emomgev,pel(3),
     &  dpphoton(3)

      double precision erad1,cgam1,pdum,echarge1,emasskg1,emassg1,
     &  clight1,pi1

c      save

      integer icharge,ieneloss,ical,icorrz,icorry,isbig

      data ical/0/

      data pi1/3.14159265358979d0/
      data emasskg1/9.10938188D-31/
      data emassg1/0.510998902D-3/
      data clight1/2.99792458d8/
      data echarge1/1.602176462D-19/

      data bmovecut/1.0d-6/

c      if (ical.eq.0) then
        erad1=2.8179380D-15
        cgam1=4.0d0/3.0d0*pi1*erad1/emassg1**3
        pdum=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1
c        ical=1
c      endif

      dgamma=0.0d0

      isbig=0
      x1=0.0d0
      y1=0.0d0
      z1=0.0d0

      if (icharge.gt.0) then
        bx2=-bxin
        by2=-byin
        bz2=-bzin
      else
        bx2=bxin
        by2=byin
        bz2=bzin
      endif

      if(dabs(bx2).lt.bmovecut.and.dabs(by2).lt.bmovecut
     &    .and.dabs(bz2).lt.bmovecut) then

        vxp=0.d0
        vyp=0.d0
        vzp=0.d0

        x2=x1+vx1*dtim
        y2=y1+vy1*dtim
        z2=z1+vz1*dtim

        vx2=vx1
        vy2=vy1
        vz2=vz1

        goto 999

      endif !b-cut

      bsq=bx2*bx2+by2*by2+bz2*bz2
      bbet=dsqrt(bsq)
      bbet1=1.0d0/bbet

      bux=bx2*bbet1
      buy=by2*bbet1
      buz=bz2*bbet1

c
c   Betrag von v0 paralell und senkrecht
c
      v0sq=vx1*vx1+vy1*vy1+vz1*vz1
      v0bet=dsqrt(v0sq)

c
c  vpar
c
      vpar=vx1*bux+vy1*buy+vz1*buz
      vparsq=vpar*vpar
      vparx=vpar*bux
      vpary=vpar*buy
      vparz=vpar*buz

c      vsenk2=(v0sq-vparsq)
      vsenk2=(v0bet-vpar)*(v0bet+vpar)

      if (vsenk2.le.0.0d0) then

c
c  Zeitableitung der Geschwindigkeit
c
        vxp=0.d0
        vyp=0.d0
        vzp=0.d0

c
c v(dtim),x(dtim) berechnen
c

        x2=x1+vx1*dtim
        y2=y1+vy1*dtim
        z2=z1+vz1*dtim

        vx2=vx1
        vy2=vy1
        vz2=vz1

        goto 999

      else    !(vsenk2.lt.0.0)

        vsenk=dsqrt(vsenk2)
        vsenk1=1.0d0/vsenk

c
c   vektor n1 berechnen
c

        x1n=(vx1-vpar*bux)*vsenk1
        y1n=(vy1-vpar*buy)*vsenk1
        z1n=(vz1-vpar*buz)*vsenk1

c
c  Vektor n2=(bux,buy,buz) kreuz n1
c

        x2n = buy*z1n - buz*y1n
        y2n = buz*x1n - bux*z1n
        z2n = bux*y1n - buy*x1n

c
c Zyklotronfrequenz
c
        zyk=(echarge1/(gamma*emasskg1))*bbet
c
c

        tz=zyk*dtim

        if (tz.le.0.03d0) then
          tz2=tz*tz
          tz3=tz2*tz
          tz4=tz3*tz
          tz5=tz4*tz
          cz=1.0d0-tz2/2.0d0+tz4/24.0d0
          sz=tz-tz3/6.0d0+tz5/120.0d0
        else
          cz=cos(tz)
          sz=sin(tz)
          isbig=1
        endif

c
c  Zeitableitung der Geschwindigkeit
c

cerror 4.7.2018        vxp=vsenk*zyk*x2n
cerror 4.7.2018        vyp=vsenk*zyk*y2n
cerror 4.7.2018        vzp=vsenk*zyk*z2n

c
c v(dtim),x(dtim) berechnen
c

        vx2=vparx + vsenk*(x1n*cz+x2n*sz)
        vy2=vpary + vsenk*(y1n*cz+y2n*sz)
        vz2=vparz + vsenk*(z1n*cz+z2n*sz)

c
c x(dtim) berechnen
c

        vsenkzyk=vsenk/zyk

        x2=x1+vsenkzyk*(x2n+x1n*sz-x2n*cz)+vparx*dtim
        y2=y1+vsenkzyk*(y2n+y1n*sz-y2n*cz)+vpary*dtim
        z2=z1+vsenkzyk*(z2n+z1n*sz-z2n*cz)+vparz*dtim

        emomgev=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
        brho=emomgev*1.0d9/clight1

        if (ieneloss.ne.0) then
          dgamma=-pdum*bsq*vsenk2/v0sq*gamma**2*dtim
          if (ieneloss.eq.-1) then
            v0=sqrt(v0sq)
            v12n(1)=vx2/v0
            v12n(2)=vy2/v0
            v12n(3)=vz2/v0
            call uradphoton(v12n,gamma,bx2,by2,bz2,
     &        dgamma,dtim,dpphoton) !dgamma will be overwritten
            if (dgamma.ne.0.0d0) then
              emomgev=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !Gev
              pel=emomgev*v12n+dpphoton
              dgamma=
     &          sqrt(1.0d0+(pel(1)**2+pel(2)**2+pel(3)**2)/emassg1**2)-
     &          gamma
              emomgev=sqrt(pel(1)**2+pel(2)**2+pel(3)**2)
              vx2=pel(1)/((gamma+dgamma)*emassg1)*clight1
              vy2=pel(2)/((gamma+dgamma)*emassg1)*clight1
              vz2=pel(3)/((gamma+dgamma)*emassg1)*clight1
            endif !dgamma
          endif !ieneloss .eq. -1
        endif !ieneloss

      endif   !(vsenk.lt.0.0)

      acc=icharge*echarge1/(gamma*emasskg1)
cerror 4.7.2018, due to error above, now here
      vx12=(vx1+vx2)/2.0d0
      vy12=(vy1+vy2)/2.0d0
      vz12=(vz1+vz2)/2.0d0

      vxp=acc*(vy12*bzin-vz12*byin)
      vyp=acc*(vz12*bxin-vx12*bzin)
      vzp=acc*(vx12*byin-vy12*bxin)

      if (isbig.eq.0) then

        dk=(clight1*dtim)**2/brho

        dy=abs(bz2*dk)
        y12=abs(y1-y2)
        icorry=0
        if (dy.ge.2.0d0*y12
     &      .or.gamma.gt.1.0d10
     &      ) then
          y2=y1+vy12*dtim
          icorry=1
        endif

        dz=abs(by2*dk)
        z12=abs(z1-z2)
        icorrz=0
        if (dz.ge.2.0d0*z12
     &      .or.gamma.gt.1.0d10
     &      ) then
          z2=z1+vz12*dtim
          icorrz=1
        endif

        if (icorry.eq.1.and.icorrz.eq.1) then
          x2=x1+vx1*dtim
        endif

      endif !(isbig.eq.0) then

999   continue

      dgammao=dgamma
      if (efieldx.ne.0.0d0.or.efieldy.ne.0.0d0.or.efieldz.ne.0.0d0) then
        call uradestep(icharge,x2,y2,z2,vx2,vy2,vz2,
     &    efieldx,efieldy,efieldz,dtim,gamma,dgamma)
        dgamma=dgammao+dgamma
        vxp=(vx2-vx1)/dtim
        vyp=(vy2-vy1)/dtim
        vzp=(vz2-vz1)/dtim
      endif

      x2=x2+xin
      y2=y2+yin
      z2=z2+zin

      return
      end
*CMZ :  4.00/15 28/04/2022  11.44.06  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.64/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.36/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/10 21/08/96  12.30.33  by  Michael Scheer
*CMZ : 00.01/04 30/11/94  14.09.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.18  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.57  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine bhalba_omp(b0halba,perlen,XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT)

      implicit none

      double precision, parameter :: twopi1=6.2831853071795862d0
      double precision zkz,yky,dnszkz,dshyky,dsnzkz,dcszkz,dchyky,
     &  b0halba,perlen,xin,yin,zin,bxout,byout,bzout,halk

      halk=twopi1/perlen

      yky=halk*yin
      zkz=halk*xin

      dshyky=dsinh(yky)
      dchyky=dsqrt(1.0d0+dshyky*dshyky)
      dsnzkz=dsin(zkz)
      dcszkz=dcos(zkz)

      bxout=-b0halba*dshyky*dsnzkz
      byout= b0halba*dchyky*dcszkz
      bzout=0.0d0

      return
      end
*CMZ :  4.00/13 27/10/2021  13.40.54  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.20.35  by  Michael Scheer
*CMZ :  3.05/04 28/06/2018  09.35.09  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.46.53  by  Michael Scheer
*CMZ :  3.02/05 25/03/2015  09.51.10  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  15.45.26  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/03 01/09/2012  16.13.42  by  Michael Scheer
*CMZ :  2.68/02 08/06/2012  09.54.11  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine uradphoton(veln,gamma,bx,by,bz,dgamma,dtim,dpphoton)

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

c NO WARRANTY

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      integer, parameter :: nbing1=1000
      integer :: ical=0,i

      double precision veln(3),bmag(3),bx,by,bz,ebeam,elmom,gamma,
     &  bparn,bper(3),bpern,epho,eec,bpervn(3),
     &  dgamma,b2per,dpphoton(3),
     &  dtim,ec,photons,de,deecg1,eecg1,g1,yrnint10

      real rnrn
c      double precision, dimension (:), allocatable ::
c     &  xrn,yrn,yrnint,coef,work1,work2,work3,work4

      double precision ::
     &  xrn(nbing1)=0.0d0,yrn(nbing1)=0.0d0,yrnint(nbing1)=0.0d0,
     &  coef(nbing1)=0.0d0,
     &  work1(nbing1)=0.0d0,work2(nbing1)=0.0d0,
     &  work3(nbing1)=0.0d0,work4(nbing1)=0.0d0

      save ical,xrn,yrn,yrnint,coef

      double precision :: eecmaxg1=5.0d0

      if (ical.eq.0) then

        deecg1=eecmaxg1/(nbing1-1)
        eecg1=0.0d0

        do i=1,10
          eecg1=eecg1+deecg1/10.0d0
          call util_g1_static(eecg1,g1)
          xrn(i)=eecg1
          yrn(i)=g1/eecg1
        enddo

        yrnint(1)=0.0d0
        do i=2,10
          yrnint(i)=yrnint(i-1)
     &      +(yrn(i)+yrn(i-1))/2.0d0*(xrn(i)-xrn(i-1))
        enddo
        yrnint10=yrnint(10)

        eecg1=0.0d0
        do i=10,nbing1
          eecg1=eecg1+deecg1
          call util_g1_static(eecg1,g1)
          xrn(i)=eecg1
          yrn(i)=g1/eecg1
        enddo

        call util_spline_running_integral(
     &    xrn(10:nbing1),yrn(10:nbing1),nbing1-10+1,yrnint(10:nbing1),
     &    coef,work1,work2,work3,work4)

        yrnint(10)=yrnint10
        yrnint(11:nbing1)=yrnint(11:nbing1)+yrnint10
        yrnint=yrnint/yrnint(nbing1)

        do i=2,nbing1
          if (yrnint(i).le.yrnint(i-1)) then
            stop '*** Error in photon: Bad integration of G1 ***'
          endif
        enddo

        call util_spline_coef(
     &    yrnint,xrn,nbing1,0.0d0,0.0d0,
     &    coef,work1,work2,work3,work4)

        ical=1
      endif !ical

      bmag(1)=bx
      bmag(2)=by
      bmag(3)=bz

      elmom=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
      ebeam=emassg1*gamma !GeV

      bparn=(bmag(1)*veln(1)+bmag(2)*veln(2)+bmag(3)*veln(3))
      bper=bmag-bparn*veln
      b2per=bper(1)**2+bper(2)**2+bper(3)**2
      bpern=sqrt(b2per)
      bpervn=bper/bpern

      ec=ecdipkev1*bpern*ebeam**2*1.0d-6 !GeV

      call uradrndm(rnrn)  !S. 39

      !dgamma = pdum * gamma**2 * b2per * dtim
      !dN = 15*sqrt(3)/8 * dE/Ec = 3.2476 * de/ec

      de=powcon1*b2per*gamma*ebeam*dtim !GeV

      if (ec.ne.0.0d0) then
        photons=3.2476d0*de/ec !number of photons
      else
        photons=0.0d0
      endif

      call uradrndm(rnrn)

      if(rnrn.le.photons) then

        call uradrndm(rnrn)  !s. 39
        call util_spline_inter(yrnint,xrn,coef,nbing1,
     &    dble(rnrn),eec,-1)
        if (eec.lt.0.0d0) then
          print*,
     &      '*** Warning in PHOTON: Negative photon energy occured ***'
          print*,'rnrn:',rnrn
          print*,'setting Epho/Ec = 1.e-6'
          eec=1.0d-6
        endif

        epho=eec*ec
        dgamma=-epho/ebeam*gamma
        dpphoton=-veln*epho

      else

        dpphoton=0.0d0
        dgamma=0.0d0

      endif !(rnrn.le.wrad)

      return
      end
*CMZ :  3.05/20 31/10/2018  10.51.37  by  Michael Scheer
*CMZ :  3.05/04 04/07/2018  14.02.25  by  Michael Scheer
*CMZ :  3.05/00 02/05/2018  11.33.05  by  Michael Scheer
*CMZ :  3.04/00 01/03/2018  15.00.40  by  Michael Scheer
*CMZ :  3.02/04 11/03/2015  13.19.06  by  Michael Scheer
*CMZ :  3.01/04 20/05/2014  12.49.15  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/05 06/09/2012  15.53.22  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  09.07.11  by  Michael Scheer
*CMZ :  2.68/03 27/08/2012  09.35.31  by  Michael Scheer
*-- Author :    Michael Scheer
c*******************************************************************************
      subroutine uradstep(xin,yin,zin,vx1,vy1,vz1,bxin,byin,bzin,
     &  efieldx,efieldy,efieldz,
     &  dtim,
     &  x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,
     &  ieneloss,dgamma)
c*******************************************************************************

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

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

      implicit none

      double precision bbet1,vsenk1,vsenk2,tz,tz2,tz3,tz4,tz5,
     &  vsenkzyk,acc,
     &  efieldx,efieldy,efieldz,
     &  dk,z12,y12,vz12,vy12,vx12,xin,yin,zin

      double precision bx2,by2,bz2,bsq,bbet,bux,buy,buz,v0sq,
     &  vx1,vy1,vz1,
     &  v0bet,vpar,vparx,vpary,vparz,vparsq,vsenk,dtim,bxin,byin,bzin
      double precision x1n,y1n,z1n,x2n,y2n,z2n

      double precision x2,y2,z2,vx2,vy2,vz2,x1,y1,z1,vxp,vyp,vzp,zyk,
     &  gamma,sz,cz,zp,yp,dy,dz,brho
     &  ,dgamma

      double precision bmovecut,dgammao,v0,v12n(3),emomgev,pel(3),
     &  dpphoton(3)

      double precision erad1,cgam1,pdum,echarge1,emasskg1,emassg1,
     &  clight1,pi1

      save

      integer icharge,ieneloss,ical,icorrz,icorry,isbig

      data ical/0/

      data pi1/3.14159265358979d0/
      data emasskg1/9.10938188D-31/
      data emassg1/0.510998902D-3/
      data clight1/2.99792458d8/
      data echarge1/1.602176462D-19/

      data bmovecut/1.0d-6/

      if (ical.eq.0) then
        erad1=2.8179380D-15
        cgam1=4.0d0/3.0d0*pi1*erad1/emassg1**3
        pdum=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1
        ical=1
      endif

      dgamma=0.0d0

      isbig=0
      x1=0.0d0
      y1=0.0d0
      z1=0.0d0

      if (icharge.gt.0) then
        bx2=-bxin
        by2=-byin
        bz2=-bzin
      else
        bx2=bxin
        by2=byin
        bz2=bzin
      endif

      if(dabs(bx2).lt.bmovecut.and.dabs(by2).lt.bmovecut
     &    .and.dabs(bz2).lt.bmovecut) then

        vxp=0.d0
        vyp=0.d0
        vzp=0.d0

        x2=x1+vx1*dtim
        y2=y1+vy1*dtim
        z2=z1+vz1*dtim

        vx2=vx1
        vy2=vy1
        vz2=vz1

        goto 999

      endif !b-cut

      bsq=bx2*bx2+by2*by2+bz2*bz2
      bbet=dsqrt(bsq)
      bbet1=1.0d0/bbet

      bux=bx2*bbet1
      buy=by2*bbet1
      buz=bz2*bbet1

c
c   Betrag von v0 paralell und senkrecht
c
      v0sq=vx1*vx1+vy1*vy1+vz1*vz1
      v0bet=dsqrt(v0sq)

c
c  vpar
c
      vpar=vx1*bux+vy1*buy+vz1*buz
      vparsq=vpar*vpar
      vparx=vpar*bux
      vpary=vpar*buy
      vparz=vpar*buz

c      vsenk2=(v0sq-vparsq)
      vsenk2=(v0bet-vpar)*(v0bet+vpar)

      if (vsenk2.le.0.0d0) then

c
c  Zeitableitung der Geschwindigkeit
c
        vxp=0.d0
        vyp=0.d0
        vzp=0.d0

c
c v(dtim),x(dtim) berechnen
c

        x2=x1+vx1*dtim
        y2=y1+vy1*dtim
        z2=z1+vz1*dtim

        vx2=vx1
        vy2=vy1
        vz2=vz1

        goto 999

      else    !(vsenk2.lt.0.0)

        vsenk=dsqrt(vsenk2)
        vsenk1=1.0d0/vsenk

c
c   vektor n1 berechnen
c

        x1n=(vx1-vpar*bux)*vsenk1
        y1n=(vy1-vpar*buy)*vsenk1
        z1n=(vz1-vpar*buz)*vsenk1

c
c  Vektor n2=(bux,buy,buz) kreuz n1
c

        x2n = buy*z1n - buz*y1n
        y2n = buz*x1n - bux*z1n
        z2n = bux*y1n - buy*x1n

c
c Zyklotronfrequenz
c
        zyk=(echarge1/(gamma*emasskg1))*bbet
c
c

        tz=zyk*dtim

        if (tz.le.0.03d0) then
          tz2=tz*tz
          tz3=tz2*tz
          tz4=tz3*tz
          tz5=tz4*tz
          cz=1.0d0-tz2/2.0d0+tz4/24.0d0
          sz=tz-tz3/6.0d0+tz5/120.0d0
        else
          cz=cos(tz)
          sz=sin(tz)
          isbig=1
        endif

c
c  Zeitableitung der Geschwindigkeit
c

cerror 4.7.2018        vxp=vsenk*zyk*x2n
cerror 4.7.2018        vyp=vsenk*zyk*y2n
cerror 4.7.2018        vzp=vsenk*zyk*z2n

c
c v(dtim),x(dtim) berechnen
c

        vx2=vparx + vsenk*(x1n*cz+x2n*sz)
        vy2=vpary + vsenk*(y1n*cz+y2n*sz)
        vz2=vparz + vsenk*(z1n*cz+z2n*sz)

c
c x(dtim) berechnen
c

        vsenkzyk=vsenk/zyk

        x2=x1+vsenkzyk*(x2n+x1n*sz-x2n*cz)+vparx*dtim
        y2=y1+vsenkzyk*(y2n+y1n*sz-y2n*cz)+vpary*dtim
        z2=z1+vsenkzyk*(z2n+z1n*sz-z2n*cz)+vparz*dtim

        emomgev=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
        brho=emomgev*1.0d9/clight1

        if (ieneloss.ne.0) then
          dgamma=-pdum*bsq*vsenk2/v0sq*gamma**2*dtim
          if (ieneloss.eq.-1) then
            v0=sqrt(v0sq)
            v12n(1)=vx2/v0
            v12n(2)=vy2/v0
            v12n(3)=vz2/v0
            call uradphoton(v12n,gamma,bx2,by2,bz2,
     &        dgamma,dtim,dpphoton) !dgamma will be overwritten
            if (dgamma.ne.0.0d0) then
              emomgev=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !Gev
              pel=emomgev*v12n+dpphoton
              dgamma=
     &          sqrt(1.0d0+(pel(1)**2+pel(2)**2+pel(3)**2)/emassg1**2)-
     &          gamma
              emomgev=sqrt(pel(1)**2+pel(2)**2+pel(3)**2)
              vx2=pel(1)/((gamma+dgamma)*emassg1)*clight1
              vy2=pel(2)/((gamma+dgamma)*emassg1)*clight1
              vz2=pel(3)/((gamma+dgamma)*emassg1)*clight1
            endif !dgamma
          endif !ieneloss .eq. -1
        endif !ieneloss

      endif   !(vsenk.lt.0.0)

      acc=icharge*echarge1/(gamma*emasskg1)
cerror 4.7.2018, due to error above, now here
      vx12=(vx1+vx2)/2.0d0
      vy12=(vy1+vy2)/2.0d0
      vz12=(vz1+vz2)/2.0d0

      vxp=acc*(vy12*bzin-vz12*byin)
      vyp=acc*(vz12*bxin-vx12*bzin)
      vzp=acc*(vx12*byin-vy12*bxin)

      if (isbig.eq.0) then

        dk=(clight1*dtim)**2/brho

        dy=abs(bz2*dk)
        y12=abs(y1-y2)
        icorry=0
        if (dy.ge.2.0d0*y12
     &      .or.gamma.gt.1.0d10
     &      ) then
          y2=y1+vy12*dtim
          icorry=1
        endif

        dz=abs(by2*dk)
        z12=abs(z1-z2)
        icorrz=0
        if (dz.ge.2.0d0*z12
     &      .or.gamma.gt.1.0d10
     &      ) then
          z2=z1+vz12*dtim
          icorrz=1
        endif

        if (icorry.eq.1.and.icorrz.eq.1) then
          x2=x1+vx1*dtim
        endif

      endif !(isbig.eq.0) then

999   continue

      dgammao=dgamma
      if (efieldx.ne.0.0d0.or.efieldy.ne.0.0d0.or.efieldz.ne.0.0d0) then
        call uradestep(icharge,x2,y2,z2,vx2,vy2,vz2,
     &    efieldx,efieldy,efieldz,dtim,gamma,dgamma)
        dgamma=dgammao+dgamma
        vxp=(vx2-vx1)/dtim
        vyp=(vy2-vy1)/dtim
        vzp=(vz2-vz1)/dtim
      endif

      x2=x2+xin
      y2=y2+yin
      z2=z2+zin

      return
      end
*CMZ :  3.03/02 19/11/2015  13.37.52  by  Michael Scheer
*CMZ :  3.02/04 11/03/2015  13.19.06  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  13.24.02  by  Michael Scheer
*CMZ :  2.68/03 25/08/2012  14.52.34  by  Michael Scheer
*CMZ :  2.66/21 22/11/2011  13.51.05  by  Michael Scheer
*-- Author :    Michael Scheer
      subroutine uradestep(icharge,x,y,z,vx,vy,vz,
     &  efieldx,efieldy,efieldz,
     &  dt,gamma,dgamma)

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

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

c Simple approach to electrical fields into account
c The energy and gamma of the particle is not changed, but returned in dgamma

      double precision x,y,z,vx,vy,vz,dt,gamma,dgamma
      double precision px,py,pz,qdt,eg,vn,pn,
     &  px0,py0,pz0,pn0,pn2,pn02,de,dpx,dpy,dpz,dtm,
     &  efieldx,efieldy,efieldz,
     &  echarge1,emasskg1,emassg1

      integer icharge

      data echarge1/1.602176462D-19/
      data emasskg1/9.10938188D-31/
      data emassg1/0.510998902D-3/

      if (efieldx.eq.0.0d0.and.efieldy.eq.0.0d0.and.efieldz.eq.0.0d0)
     &    then
        dgamma=0.0d0
        return
      endif

      qdt=icharge*echarge1*dt
      eg=emasskg1*gamma

      px0=vx*eg !SI-units
      py0=vy*eg
      pz0=vz*eg
      pn02=px0*px0+py0*py0+pz0*pz0
      pn0=sqrt(pn02)

      dpx=efieldx*qdt !SI-units
      dpy=efieldy*qdt
      dpz=efieldz*qdt

      px=px0+dpx !SI-units
      py=py0+dpy
      pz=pz0+dpz

      vn=sqrt(vx*vx+vy*vy+vz*vz)
      pn2=px*px+py*py+pz*pz
      pn=sqrt(pn2)

c total momentum and energy are kept!!
      vx=px/pn*vn
      vy=py/pn*vn
      vz=pz/pn*vn

      dtm=0.5d0*dt/(emasskg1*gamma) ! F=dp/dt, m=m_e*gamma, a=F/m*dp/dt/m_e/gamma

      x=x+dtm*dpx
      y=y+dtm*dpy
      z=z+dtm*dpz

      de=(vx*dpx+dpy*vy+dpz*vz)/echarge1/1.0d9 !GeV
      dgamma=de/emassg1

      return
      end
*CMZ :  3.02/04 13/03/2015  10.24.40  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  15.45.49  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/03 29/08/2012  09.57.11  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2012
      subroutine uradrndm(rn)

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

c NO WARRANTY

      implicit none
      real rn

      save

      call random_number(rn)

      return
      end
