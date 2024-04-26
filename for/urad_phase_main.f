*CMZ :          16/04/2024  14.41.20  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.32.18  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  10.57.05  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  13.32.32  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.57.49  by  Michael Scheer
*-- Author : Michael Scheer
*KEEP,gplhint.
*KEND.
      program urad_phase_main

      use omp_lib
      use uradphasemod
      use wignermod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double precision :: ar(10000),ai(10000)

      double precision, dimension(:), allocatable :: z,y
      double precision, dimension(:,:), allocatable :: s

      double precision :: banwid=0.001,xbeta=0.0d0,
     &  perlen,shift,ebeam,curr,step,perl,
     &  pincen(3),pinw,pinh,park,wlen1,gamma,
     &  ephmin,ephmax,beffv,beffh,pherror,stosum(4),
     &  alphah,alphav,espread,harm,b0eff,rhv,fsum,
     &  betah,betav,pinx,piny,pinz,zz,yy,dzpin,dypin,zmin,ymin,
     &  disph,dispph,dispv,disppv,bunchcharge,bunchlen,efi(3),bfi(3),rn(3),
     &  emith,emitv,pinxprop,pinwprop,pinhprop,dy,dz,dthez,dthey,thezwig,theywig

      real xran(1),rr(2),axr,axi,ayr,ayi,azr,azi

      integer :: idebug=0,noranone,i,ktime=1,
     &  npiny,npinz,nper,nepho,modeph,modepin,modesphere,nharm,iy,iz,iobs,
     &  mthreads,nelec,icohere,ihbunch,iepho,iobph,iel,modebunch,ifieldprop,
     &  modewave=0,isto,nlpoi=0,nobsvprop,npinyprop,npinzprop,
     &  iwigner,jz,jy,
     &  ispline=0

      namelist/uradphasen/
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,noranone,nobsvprop,npinyprop,npinzprop,
     &  pinx,piny,pinz,pinw,pinh,npiny,npinz,modepin,modesphere,nharm,harm,
     &  nepho,ephmin,ephmax,pherror,pinxprop,pinwprop,pinhprop,
     &  mthreads,nelec,icohere,ihbunch,modeph,modebunch,ifieldprop,iwigner,
     &  betah,betav,alphah,alphav,emith,emitv,espread,
     &  disph,dispph,dispv,disppv,bunchcharge,bunchlen,
     &  nzwig,nywig,thezwig,theywig,nzthewig,nythewig,pinhwig,pinwwig

      integer :: irnsize=64,irnseed(64),ifixseed
      namelist/seedn/irnseed,ifixseed

      integer :: luna,istat,kalloc=1

c      if (ktime.eq.1) then
c        call util_zeit_kommentar_delta(6,'Running urad_phase',1)
c      endif

      open(newunit=luna,file='urad_phase.nam',status='old',iostat=istat)
      if (istat.ne.0) then
        stop "*** Error: Could not open urad_phase.nam"
      endif

      read(luna,uradphasen)
      read(luna,seedn)

      close(luna)

      if(nelec.eq.1.and.noranone.eq.0) then
        noranone=1
        print*
        print*,'*** Changed NORANONE=0 to NORANONE=1, since NELEC=1'
        print*
      endif

      pincen=[pinx,piny,pinz]
      emith=emith*1.0d-9
      emitv=emitv*1.0d-9

      bunchlen=bunchlen/1.0d9 !nm->m

      !print*,"sigz, sizp:",sigz/1000.0d0,sigzp/1000.0d0
      !print*,"sigy, siyp:",sigy/1000.0d0,sigyp/1000.0d0

      if (ifixseed.ne.0) then
        ifixseed=1
        call util_random_set_seed(irnsize,irnseed)
      endif

      if (nharm.gt.0.and.harm.gt.0.0d0) then

        gamma=ebeam/emassg1
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

      endif

      npiny=max(1,npiny)
      npinz=max(1,npinz)

      npinyo_u=npiny
      npinzo_u=npinz

      open(newunit=luna,file='urad_phase.pin')
      write(luna,*)npinz,npiny,pinw,pinh
      write(luna,*)pincen
      close(luna)

      ymin=pincen(2)-pinh/2.0d0
      zmin=pincen(3)-pinw/2.0d0
      dypin=pinh/max(1,npiny-1)
      dzpin=pinw/max(1,npinz-1)

      if (iwigner.ne.0) then

        if (ifieldprop.eq.0) then
          pinxprop=pincen(1)
          pinwprop=pinw
          pinhprop=pinh
          npinzprop=npinz
          npinyprop=npiny
          dyprop=pinh/dble(max(1,npiny-1))/1000.0d0 !meter
          dzprop=pinw/dble(max(1,npinz-1))/1000.0d0 !meter
        else
          dyprop=pinhprop/dble(max(1,npinyprop-1))/1000.0d0
          dzprop=pinwprop/dble(max(1,npinzprop-1))/1000.0d0
        endif

        allocate(wigz(nzwig),wigy(nywig),wigthez(nzthewig),wigthey(nythewig))

        if (theywig.le.0.0d0) theywig=pinhprop/pinxprop*1000.0d0

        if (nywig.gt.1) then
          dywig=pinhwig/dble(nywig-1)
          wigy(1)=-pinhwig/2.0d0
          do iy=2,nywig
            wigy(iy)=wigy(iy-1)+dywig
            if (abs(wigy(iy)).lt.1.0d-12) wigy(iy)=0.0d0
          enddo
        else
          dywig=pinhwig
          wigy(1)=0.0d0
        endif

        if (nythewig.gt.1) then
          dthey=theywig/dble(max(1,nythewig-1))
          wigthey(1)=-theywig/2.0d0
          do iy=2,nythewig
            wigthey(iy)=wigthey(iy-1)+dthey
            if (abs(wigthey(iy)).lt.1.0d-12) wigthey(iy)=0.0d0
          enddo
        else
          dthey=theywig
          wigthey(1)=0.0d0
        endif

        if (thezwig.le.0.0d0) thezwig=pinwprop/pinxprop*1000.0d0

        if (nzwig.gt.1) then
          dzwig=pinwwig/dble(nzwig-1)/1000.0d0
          wigz(1)=-pinwwig/2.0d0/1000.0d0 !meter
          do iz=2,nzwig
            wigz(iz)=wigz(iz-1)+dzwig
            if (abs(wigz(iz)).lt.1.0d-12) wigz(iz)=0.0d0
          enddo
        else
          dzwig=pinhwig
          wigz(1)=0.0d0
        endif

        if (nzthewig.gt.1) then
          dthez=thezwig/dble(max(1,nzthewig-1))/1000.0d0
          wigthez(1)=-thezwig/2.0d0/1000.0d0
          do iz=2,nzthewig
            wigthez(iz)=wigthez(iz-1)+dthez
            if (abs(wigthez(iz)).lt.1.0d-12) wigthez(iz)=0.0d0
          enddo
        else
          dthez=thezwig
          wigthez(1)=0.0d0
        endif

      endif

      if (mthreads.lt.0) then
        mthreads=OMP_GET_MAX_THREADS()
      else if (mthreads.eq.0) then
        mthreads=1
      endif

      mthreads_u=mthreads

      call urad_phase(
     &  mthreads,nelec,noranone,icohere,modebunch,bunchlen,bunchcharge,ihbunch,
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,nlpoi,
     &  pincen,pinw,pinh,npiny,npinz,modepin,modesphere,
     &  nepho,ephmin,ephmax,banwid,
     &  xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  modeph,pherror,modewave
     &  )

      allocate(z(nobsv_u),y(nobsv_u),s(npinz_u,npiny_u))

      if (ifieldprop.ne.0.or.iwigner.ne.0) then

        pinxprop_u=pinxprop
        pinwprop_u=pinwprop
        pinhprop_u=pinhprop

        npinyprop_u=npinyprop
        npinzprop_u=npinzprop

        if (ifieldprop.ne.0) then
          call urad_phase_prop(mthreads)
        else
          nobsvprop_u=nobsv_u
          obsvprop_u=obsv_u
          stokesprop_u=stokes_u
          aradprop_u=arad_u
        endif


        open(newunit=luna,file='urad_phase.fdp')
        do iobs=1,nobsvprop_u
          do iepho=1,nepho_u
            iobph=iobs+nobsvprop_u*(iepho-1)
            ! rnx = (eyr+i*eyi)*(bzr-i*bzi) - (ezr+i*ezi)*(byr-i*byi)
            !     = eyr*bzr - i*eyr*bzi + i*eyi*bzr + eyi*bzi
            !     - ezr*byr + i*erz*byi - i*ezi*byr - ezi*byi
            ! real(rnx) = eyr*bzr + eyi*bzi - ezr*byr - ezi*byi
            rn(1)=real(aradprop_u(2,iobph)*conjg(aradprop_u(6,iobph))-aradprop_u(3,iobph)*conjg(aradprop_u(5,iobph)))
            rn(2)=real(aradprop_u(3,iobph)*conjg(aradprop_u(4,iobph))-aradprop_u(1,iobph)*conjg(aradprop_u(6,iobph)))
            rn(3)=real(aradprop_u(1,iobph)*conjg(aradprop_u(5,iobph))-aradprop_u(2,iobph)*conjg(aradprop_u(4,iobph)))
            rn=rn/norm2(rn)
            axr=real(aradprop_u(1,iobph))
            axi=imag(aradprop_u(1,iobph))
            ayr=real(aradprop_u(2,iobph))
            ayi=imag(aradprop_u(2,iobph))
            azr=real(aradprop_u(3,iobph))
            azi=imag(aradprop_u(3,iobph))
            write(luna,'(3(1pe15.6e3),i10,20(1pe15.6e3))')
     &        obsvprop_u(1:3,iobs),iepho,epho_u(iepho),stokesprop_u(1:4,iobph),
     &        real(aradprop_u(1,iobph)),imag(aradprop_u(1,iobph)),
     &        real(aradprop_u(2,iobph)),imag(aradprop_u(2,iobph)),
     &        real(aradprop_u(3,iobph)),imag(aradprop_u(3,iobph)),
     &        real(aradprop_u(4,iobph)),imag(aradprop_u(4,iobph)),
     &        real(aradprop_u(5,iobph)),imag(aradprop_u(5,iobph)),
     &        real(aradprop_u(6,iobph)),imag(aradprop_u(6,iobph)),
     &        rn
          enddo
        enddo
        close(luna)

        if (iwigner.ne.0) then
          call urad_phase_wigner(ispline)
        endif

      endif

      if (ktime.eq.1) then
        call util_zeit_kommentar_delta(6,'Writing results to files',1)
      endif

      open(newunit=luna,file='urad_phase.fld')

      do iobs=1,nobsv_u
        do iepho=1,nepho_u
          iobph=iobs+nobsv_u*(iepho-1)

          rn(1)=real(arad_u(2,iobph)*conjg(arad_u(6,iobph))-arad_u(3,iobph)*conjg(arad_u(5,iobph)))
          rn(2)=real(arad_u(3,iobph)*conjg(arad_u(4,iobph))-arad_u(1,iobph)*conjg(arad_u(6,iobph)))
          rn(3)=real(arad_u(1,iobph)*conjg(arad_u(5,iobph))-arad_u(2,iobph)*conjg(arad_u(4,iobph)))
          rn=rn/norm2(rn)

          axr=real(arad_u(1,iobph))
          axi=imag(arad_u(1,iobph))
          ayr=real(arad_u(2,iobph))
          ayi=imag(arad_u(2,iobph))
          azr=real(arad_u(3,iobph))
          azi=imag(arad_u(3,iobph))

          write(luna,'(3(1pe15.6e3),i10,21(1pe15.6e3))')
     &      obsv_u(1:3,iobs),iepho,epho_u(iepho),stokes_u(1:4,iobph),pow_u(iobs),
     &      real(arad_u(1,iobph)),imag(arad_u(1,iobph)),
     &      real(arad_u(2,iobph)),imag(arad_u(2,iobph)),
     &      real(arad_u(3,iobph)),imag(arad_u(3,iobph)),
     &      real(arad_u(4,iobph)),imag(arad_u(4,iobph)),
     &      real(arad_u(5,iobph)),imag(arad_u(5,iobph)),
     &      real(arad_u(6,iobph)),imag(arad_u(6,iobph)),
     &      rn

        enddo
      enddo
      close(luna)

      if (modepin_u.ne.0) then
        open(newunit=luna,file='urad_phase.ebm')
        do iepho=1,nepho_u
          do iy=1,npinyo_u-1
            yy=ymin+(dble(iy)-0.5)*dypin
            if(abs(yy).lt.1.0e-9) yy=0.0d0
            do iz=1,npinzo_u-1
              zz=zmin+(dble(iz)-0.5)*dzpin
              if(abs(zz).lt.1.0e-9) zz=0.0d0
c              fsum=max(1.0d0,dreal(fieldbunch(7,iz,iy,iepho)))
              write(luna,*)iepho,iz,iy,epho_u(iepho),zz,yy,
     &          dreal(fieldbunch(1,iz,iy,iepho)),
     &          dimag(fieldbunch(1,iz,iy,iepho)),
     &          dreal(fieldbunch(2,iz,iy,iepho)),
     &          dimag(fieldbunch(2,iz,iy,iepho)),
     &          dreal(fieldbunch(3,iz,iy,iepho)),
     &          dimag(fieldbunch(3,iz,iy,iepho)),
     &          dreal(fieldbunch(4,iz,iy,iepho)),
     &          dimag(fieldbunch(4,iz,iy,iepho)),
     &          dreal(fieldbunch(5,iz,iy,iepho)),
     &          dimag(fieldbunch(5,iz,iy,iepho)),
     &          dreal(fieldbunch(6,iz,iy,iepho)),
     &          dimag(fieldbunch(6,iz,iy,iepho)),
     &          dreal(fieldbunch(7,iz,iy,iepho))
            enddo
          enddo
        enddo
        close(luna)
      endif

      open(newunit=luna,file='urad_phase.flx')

      z(1:npinz_u)=obsv_u(3,1:npinz_u)

      do iy=1,npiny_u
        iobs=npinz_u*(iy-1)+1
        y(iy)=obsv_u(2,iobs)
      enddo

      do iepho=1,nepho_u
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
            call util_spline_integral_2d(npinz_u,npiny_u,z,y,s,stosum(isto),
     &        istat,kalloc)
            kalloc=0
          enddo !isto
          write(luna,*)iepho,epho_u(iepho),stosum
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
          write(luna,*)iepho,epho_u(iepho),stosum/nobsv_u*pinw*pinh
        endif
      enddo
      close(luna)

      if (ihbunch.gt.0) then
        open(newunit=luna,file='urad_phase.bun')
        do iel=1,nelec_u/iabs(ihbunch_u)*nepho_u
          if(fbunch_u(21,iel).ne.0.0d0) then
            write(luna,*)fbunch_u(:,iel)
          endif
        enddo
        close(luna)
      else if (ihbunch.lt.0) then
        open(newunit=luna,file='urad_phase.bun')
        do iel=1,nobsv_u*nelec_u/iabs(ihbunch_u)*nepho_u
          if(fbunch_u(21,iel).ne.0.0d0) then
            write(luna,*)fbunch_u(:,iel)
          endif
        enddo
        close(luna)
      endif

      call  util_random_get_seed(irnsize,irnseed)

      open(newunit=luna,file='urad_phase.seeds',status='unknown')
      write(luna,*)irnsize
      do i=1,irnsize
        write(luna,*)i,irnseed(i)
      enddo
      flush(luna)
      close(luna)

      if (iwigner.ne.0) then
        open(newunit=luna,file='urad_phase.wig')
        do iepho=1,nepho
          do iz=1,nzwig
            do iy=1,nywig
              do jz=1,nzthewig
                do jy=1,nythewig
                  write(luna,'(4i10,5(1pe15.6e3),i10,9(1pe15.6e3))')
     &              iz,iy,jz,jy,pinXprop,
     &              wigy(iy)*1000.0d0,wigz(iz)*1000.0d0,wigthey(jy)*1000.0d0,wigthez(jz)*1000.0d0,
     &              iepho,epho_u(iepho),
     &              wigr(iz,iy,jz,jy,iepho,1)*1000.0d0,wigi(iz,iy,jz,jy,iepho,1)*1000.0d0,
     &              wigr(iz,iy,jz,jy,iepho,2)*1000.0d0,wigi(iz,iy,jz,jy,iepho,2)*1000.0d0,
     &              wigr(iz,iy,jz,jy,iepho,3)*1000.0d0,wigi(iz,iy,jz,jy,iepho,3)*1000.0d0,
     &              wigr(iz,iy,jz,jy,iepho,4)*1000.0d0,wigi(iz,iy,jz,jy,iepho,4)*1000.0d0
                enddo
              enddo
            enddo
          enddo
        enddo
        close(luna)
      endif

      if (ktime.eq.1) then
        call util_zeit_kommentar_delta(6,'Done',0)
      endif

      end
*CMZ :          26/04/2024  10.49.56  by  Michael Scheer
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
     &  modeph,pherror,modewave
     &  )

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      double precision
     &  perlen,shift,ebeam,curr,step,banwid,
     &  pincen(3),pinw,pinh,betah,alphah,betav,alphav,
     &  ephmin,ephmax,beffv,beffh,pherror,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,y,z,dy,dz,ymin,zmin,bunchlen,bunchcharge,
     &  xbeta,df,xx,yy,zz,r,xn,yn,zn,h2

      integer
     &  npiny,npinz,nper,nepho,mthreads,nelec,icohere,ihbunch,i,nlpoi,
     &  modeph,modepin,modesphere,modebunch,iy,iz,iobsv,noranone,modewave

c      if (modewave.ne.0) call util_zeit_kommentar(6,'Entered urad_phase')
      call util_zeit_kommentar_delta(6,'Entered urad_phase',1)

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

      if (modepin.eq.0) then
        nobsv_u=npiny_u*npinz_u
      else
        npinz_u=1
        npiny_u=1
        nobsv_u=1
      endif

      allocate(epho_u(nepho),obsv_u(3,nobsv_u),
     &  arad_u(6,nobsv_u*nepho_u),
     &  specpow_u(nobsv_u),
     &  stokes_u(4,nobsv_u*nepho_u),pow_u(nobsv_u)
     &  )

      if (ihbunch_u.gt.0) then
        allocate(fbunch_u(41,nelec_u/ihbunch_u*nepho_u))
        fbunch_u=0.0d0
      else if (ihbunch_u.lt.0) then
        allocate(fbunch_u(41,nobsv_u*nelec_u/(-ihbunch_u)*nepho_u))
        fbunch_u=0.0d0
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
      arad_u=arad_u/1.0d3

      obsv_u=obsv_u*1000.0d0

c      if (modewave.ne.0) call util_zeit_kommentar(6,'Leaving urad_phase')
      call util_zeit_kommentar_delta(6,'Leaving urad_phase',0)

      end
*CMZ :          26/04/2024  07.41.13  by  Michael Scheer
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

      complex*16 :: cde,czero=(0.0d0,0.0d0)
      double precision :: h2,ddist,wlen,dphi,phase0

      double complex , dimension (:,:), allocatable :: aradbuff
      double complex , dimension (:,:,:), allocatable :: arad

      double precision, dimension (:), allocatable :: frq
      double precision, dimension (:,:), allocatable :: wsstokes,pow
      double precision, dimension (:,:,:), allocatable :: fbunch,stokes

      real, dimension (:), allocatable :: pherr,pherrc,phiran
      real, dimension(:,:), allocatable :: pranall,eall

      real eran(6),pran(2),rr(2)

      double complex :: apol,amp0(6),damp(6),amp(6),zexp,
     &  apolh,apolr,apoll,apol45,stokesv(4,3),cero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double precision :: t,udgamtot,upow,vf0,vn,vx0,vx2,vxf0,vxi,vy0,vy2,vyf0,
     &  vyi,vz0,vz2,vzf0,vzi,wlen1,x0,x2,xf0,xi,xlell,y0,y2,yf0,yi,ypi,yy,yyp,
     &  z0,z2,zf0,zi,zpi,zz,zzp,fillb(41),stok1,stok2,stok3,stok4,speknor,
     &  sqnbunch,sqnphsp,specnor,sbnor,rpin,r00(3),xph0,
     &  r(3),r0(3),pw,ph,phsum,pkerr,pherror,ppin,parke,pc(3),pcbrill(3),om1,
     &  park,pr,hbarev,obs(3),om,fhigh,flow,gamma,eix,eiy,eiz,emassg,
     &  efx,efy,efz,eharm1,ecdipev,ebeam,dtpho,dt,dtelec,dd0,debeam,
     &  drn0(3),drn00(3),ds,dr0(3),dr00(3),drn(3),dpp,dph,dist,dist0,dobs(3),
     &  bunnor,clight,bunchx,beta,beff,spow,
     &  zp0,yp0,rph,anor,fsum,
     &  xkellip,zampell,yampell,parkv,parkh,zpampell,ypampell,emom,dzpin,dypin,zmin,ymin

      double complex, dimension (:), allocatable ::
     &  uampex,uampey,uampez,uampbx,uampby,uampbz
      double precision, dimension (:,:), allocatable :: utraxyz,ustokes

      integer :: kfreq,iobsv,i,np2,nelec,mbunch,meinbunch,ibu,jbun,
     &  kran=6,icbrill,ilo,kobsv,i1,i2,n,
     &  ifail,ndimu,nstepu,ith,noespread,noemit,jbunch,jubunch,jhbunch,
     &  jcharge=-1,lmodeph,nclo,jeneloss=0,iamppin,
     &  iamppincirc=0,ifrob,iobfr,isub,jvelofield=0,nlbu=0,nepho,ielo,
     &  modewave,iepho

      integer, dimension (:), allocatable :: lnbunch

      integer :: idebug=0, lbunch=0, ierr=0, ielec=0
      integer ibunch,ihbunch,mthreads,nobsv,nobsvo,iemit,noranone,iz,iy,nobsvz,nobsvy
      integer iobm,iobp,iobfrm,iobfrp

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

      if (nelec_u.ne.nelec) then
        print*,''
        print*,'--- Warning in urad_amprep: Nelec adjusted to multiple of number of threads:',nelec_u
        print*,''
      endif

      if (nelec_u.eq.1) then
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

      if (modepin_u.ne.0) then
        allocate(fieldbunch(7,npinzo_u,npinyo_u,nepho_u),stat=ierr)
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
c*** OBSOLITE, SEE z0= further down
        xkellip=twopi1/perlen_u
        zampell=beffv_u*clight1/emom/xkellip**2
        yampell=beffh_u*clight1/emom/xkellip**2
c        zampell=zmx
c        yampell=ymx
        parkh=echarge1*dabs(beffh_u)*perlen_u/(twopi1*emasskg1*clight1)
        parkv=echarge1*dabs(beffv_u)*perlen_u/(twopi1*emasskg1*clight1)
        zpampell=parkv/gamma_u
c        print*,zpampell
c        ypampell=parkh/gamma_u
c        zpampell=tan(phimx)
c        print*,zpampell
c        stop
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
     &  phiran(max(1,nelec_u)))

      allocate(pranall(2,nelec_u))
      do i=1,nelec_u
        call util_random(2,pran)
        pranall(:,i)=pran
      enddo

      if (ibunch.eq.0.or.
     &    emith_u.eq.0.0d0.and.emitv_u.eq.0.0d0.and.espread_u.eq.0.0d0) then
        iemit=0
      else
        iemit=1
      endif

      if (iemit.ne.0) then
        allocate(eall(6,nelec_u))
        do i=1,nelec_u
          xi=x0
          call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &      emith_u,emitv_u,
     &      disph_u,dispph_u,dispv_u,disppv_u,
     &      espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
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

      if (pherror.ne.0.0d0.and.(lmodeph.lt.0.or.lmodeph.gt.2)) then
        write(6,*) ""
        write(6,*) "*** Error in urad_amprep: MODEPH must be 0,1, or 2 ***"
        write(6,*) "*** Program aborted ***"
      endif

      if (lmodeph.eq.0.and.eharm1.ne.0.0d0) then
        om1=eharm1/hbarev
        pherr=sngl(pherrc*pherror/360.0d0*twopi1/om1)
      else if (lmodeph.eq.1) then
        pherr=sngl(pherr*pherror)
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

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(nepho,nobsvz,nobsvy,nobsv,nelec,frq,nper_u,np2,perlen_u,clight,hbarev,
!$OMP& flow,fhigh,czero,cone,
!$OMP& x0,y0,z0,xf0,yf0,zf0,vx0,vy0,vz0,vxf0,vyf0,vzf0,gamma_u,sbnor,speknor,
!$OMP& efx,efy,efz,ds,ndimu,curr_u,xlell,parke,amp,amp0,
!$OMP& uampex,uampey,uampez,uampbx,uampby,uampbz,
!$OMP& lmodeph,zp0,yp0,modewave,
!$OMP& jbunch,jubunch,jhbunch,noespread,noemit,ebeam,
!$OMP& stokesv,icbrill,obsv_u,emassg,debeam,dispv_u,disppv_u,
!$OMP& betah_u,alphah_u,betav_u,alphav_u,emith_u,emitv_u,disph_u,dispph_u,
!$OMP& pran,pranall,eall,fillb,r0,dr0,iamppin,iamppincirc,pc,phase0,pr,banwid_u,
!$OMP& pw,ph,idebug,pcbrill,wsstokes,vn,bunchlen_u,modebunch_u,icohere_u)
!$OMP& SHARED(mthreads,stokes,pherr,lbunch,lnbunch,modepin_u,fieldbunch,npinzo_u,nobsvo,dzpin,dypin,
!$OMP& fbunch_u,jcharge,jeneloss,jvelofield,iemit,noranone,arad,pow,zmin,ymin)

      jbun=1
      isub=0
      iobsv=0
      ielo=0
      xph0=-perlen_u*dble(nper_u)/2.0d0

!$OMP DO

      do ilo=1,nelec*nobsv

        wsstokes=0.0d0

        !affe=(0.0D0,0.0D0)
        spow=0.0d0

        ith=OMP_GET_THREAD_NUM()+1

        iobsv=mod(ilo-1,nobsv)+1
        ibu=(ilo-1)/nobsv+1
        jbun=ibu

        iy=(iobsv-1)/nobsvz+1
        iz=mod(iobsv-1,nobsvz)+1

        !if (iz.gt.nobsvz/2+1) call til_break

        ielec=ibu

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

        t=bunchx/vn

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
            pran(:)=pranall(:,ielec)
            if (iamppincirc.eq.0) then
              obs(2)=pc(2)+(pran(1)-0.5)*pw
              obs(3)=pc(3)+(pran(2)-0.5)*ph
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
        r0=r0+dr0/2.0d0

        do kfreq=1,nepho

          iobfr=iobsv+nobsv*(kfreq-1)

          om=frq(kfreq)/hbarev

          if (modewave.eq.0) then
            amp0=[
     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
     &        ]*1.0d3/sqrt(speknor/curr_u*0.10d0) !urad
          else
            amp0=[
     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
     &        ]*1.0d3/sqrt(speknor) !urad
          endif

c          call util_random(1,pran)
c          amp0=amp0*dcmplx(0.0d0,dble(pran(1)*twopi1))

          amp=(0.0d0,0.0d0)

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

        enddo !kfreq

      enddo !nbunch

!$OMP END DO
!$OMP END PARALLEL

      do ith=1,mthreads
        pow_u(:)=pow_u(:)+pow(:,ith)
        arad_u(:,:)=arad_u(:,:)+arad(:,:,ith)
      enddo
      pow_u=pow_u/sqnbunch

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
     &  pherrc,pherr,phiran,arad,pow,pranall,wsstokes,stokes)

      if (iemit.ne.0) deallocate(eall)

      iobfr=nobsv_u*nepho_u/2+1
      amp(1:3)=arad_u(1:3,iobfr)

      anor=sqrt(stokes_u(1,iobfr)/
     &  (amp(1)*dconjg(amp(1))+amp(2)*dconjg(amp(2))+amp(3)*dconjg(amp(3))))
      arad_u=arad_u*anor

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
      endif

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
*CMZ :          14/04/2024  07.40.31  by  Michael Scheer
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
     &  nthstep,izaehl,nstep,ndim,kstep,lstep,kfreq,isto,ifail,ith,
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
      ith=ith
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
*CMZ :          15/04/2024  11.54.00  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      integer :: mthreads,ktime=1

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase_prop',1)
      if (modepin_u.eq.0) then
        call urad_phase_prop_classic(mthreads)
      else
        call urad_phase_prop_mc(mthreads)
      endif
      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase_prop',0)

      return
      end
*CMZ :          15/04/2024  09.37.27  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_classic(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double complex, dimension(:), allocatable :: expom,dexpom,phshift

      double complex :: apolh,apolr,apoll,apol45

      double precision dx,dx2,dy,dyph,dzph,dz,y,z,omc,domc,phlowz,phlowy,dzy2,eps(6),
     &  dr,drred,da,x,xobs,yobs,zobs,rlambda1,ans,stok1,stok2,stok3,stok4

      integer :: ktime=1,
     &  mthreads,iy,iz,n,jy,jz,iobs,ieps,ifrq,iobfr,jobs,jobfr

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      nobsvprop_u=npinyprop_u*npinzprop_u

      allocate(
     &  obsvprop_u(3,nobsvprop_u),stokesprop_u(4,nobsvprop_u*nepho_u),
     &  aradprop_u(6,nobsvprop_u*nepho_u),
     &  expom(nobsv_u*nepho_u),dexpom(nobsv_u*nepho_u),
     &  phshift(nobsv_u))

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

      da=pinw_u*pinh_u/dble(max(1,npinz_u-1)*max(1,npiny_u-1))

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
              PHSHIFT(IOBS)=EXPOM(IOBFR)
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

      do ifrq=1,nepho_u
        do jobs=1,nobsvprop_u

          jobfr=jobs+nobsvprop_u*(ifrq-1)

          apolh=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(1,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(1,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(1,3))

          apolr=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(2,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(2,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(2,3))

          apoll=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(3,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(3,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(3,3))

          apol45=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(4,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(4,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(4,3))

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

          stokesprop_u(1,jobfr)=stok1
          stokesprop_u(2,jobfr)=stok2
          stokesprop_u(3,jobfr)=stok3
          stokesprop_u(4,jobfr)=stok4

        enddo
      enddo

      obsvprop_u=obsvprop_u*1000.0d0

      return
      end
*CMZ :          15/04/2024  21.58.40  by  Michael Scheer
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

      allocate(
     &  obsvprop_u(3,nobsvprop_u),stokesprop_u(4,nobsvprop_u*nepho_u),
     &  aradprop_u(6,nobsvprop_u*nepho_u),
     &  expom(nobsv*nepho_u),dexpom(nobsv*nepho_u),
     &  phshift(nobsv),obsv(3,nobsv))

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

      do ifrq=1,nepho_u
        do jobs=1,nobsvprop_u

          jobfr=jobs+nobsvprop_u*(ifrq-1)

          apolh=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(1,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(1,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(1,3))

          apolr=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(2,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(2,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(2,3))

          apoll=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(3,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(3,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(3,3))

          apol45=
     &      aradprop_u(1,jobfr)*dconjg(vstokes(4,1))+
     &      aradprop_u(2,jobfr)*dconjg(vstokes(4,2))+
     &      aradprop_u(3,jobfr)*dconjg(vstokes(4,3))

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

          stokesprop_u(1,jobfr)=stok1
          stokesprop_u(2,jobfr)=stok2
          stokesprop_u(3,jobfr)=stok3
          stokesprop_u(4,jobfr)=stok4

        enddo
      enddo

      obsvprop_u=obsvprop_u*1000.0d0

      return
      end
*CMZ :          15/04/2024  11.54.00  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_wigner(ispline)

      use omp_lib
      use uradphasemod
      use wignermod

      implicit none

      double complex, dimension(:), allocatable :: cw
      double precision, dimension(:), allocatable :: expomz

      double complex :: czero=(0.0d0,0.0d0),ci=(0.0d0,1.0d0),omc2i,dcw,ap,am,
     &  expomy,dexpomy,expthey,dexpomthey,dexpomz,expthez,dexpomthez,edumy,edumz

      double precision :: dwz,dwy,thez,they,dzy,zp,yp,za,ya,yw,zw,ymin,ymax,zmin,zmax,dz,dy,
     &  dzzp,dzzm,dyyp,dyym,zzm,zzp,yym,yyp,pinhpropm,pinwpropm,hbar,c

      integer :: ktime=1,ispline,iepho,ixy,iz,iy,jz,jy,kx,ky,lzm,lzp,lym,lyp,lobsm,lobsp,
     &  lmph,lpph,izth,iyth,icwz,icwy,ith

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
!+seq,uservar.

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase_wigner',1)

      allocate(
     &  wigr(nzwig,nywig,nzthewig,nythewig,nepho_u,4),
     &  wigi(nzwig,nywig,nzthewig,nythewig,nepho_u,4),
     &  cw(mthreads_u),expomz(npinzprop_u)
     &  )

      wigr=0.0d0
      wigi=0.0d0

      pinhpropm=pinhprop_u/1000.0d0
      pinwpropm=pinwprop_u/1000.0d0

      dzy=dzprop*dyprop

      ymin=-pinhpropm/2.0d0
      zmin=-pinwpropm/2.0d0

      if (nywig.gt.1) then
        dwy=wigy(2)-wigy(1)
      else
        dwy=0.0d0
      endif

      if (nzwig.gt.1) then
        dwz=wigz(2)-wigz(1)
      else
        dwz=0.0d0
      endif

      c=clight1
      hbar=hbarev1

      do ixy=1,4

        kx=3
        ky=3
        if (ixy.eq.1) then
          kx=3
          ky=3
        else if (ixy.eq.2) then
          kx=3
          ky=2
        else if (ixy.eq.3) then
          kx=2
          ky=3
        else if (ixy.eq.4) then
          kx=2
          ky=2
        endif

        do iepho=1,nepho_u

          do iyth=1,nythewig

            they=wigthey(iyth)

              do izth=1,nzthewig

                thez=wigthez(izth)

                do iy=1,nywig

                  yw=wigy(iy)
                  omc2i=dcmplx(0.0d0,2.0d0*epho_u(iepho)/hbar/c)

                  do iz=1,nzwig

                    zw=wigz(iz)

                    cw=czero
                    icwy=0

                    expomy=cdexp(omc2i*they*(ymin-dyprop))
                    dexpomy=cdexp(omc2i*they*dyprop)

                    expomz(1)=cdexp(omc2i*thez*(zmin-dzprop))
                    dexpomz=cdexp(omc2i*thez*dzprop)

                    do jz=2,npinzprop_u
                      expomz(jz)=expomz(jz-1)*dexpomz
                    enddo

                    do jy=1,npinyprop_u

                      expomy=expomy*dexpomy

                      yp=ymin+(jy-1)*dyprop
                      lym=int((yw-yp-ymin)/dyprop)+1
                      if (lym.lt.1.or.lym.gt.npinyprop_u) cycle
                      lyp=int((yw+yp-ymin)/dyprop)+1
                      if (lyp.lt.1.or.lyp.gt.npinyprop_u) cycle

                      icwz=0

!$omp PARALLEL NUM_THREADS(mthreads_u) DEFAULT(PRIVATE)
!$omp& SHARED(nepho_u,nythewig,nzthewig,nzwig,nywig,wigz,wigy,wigthey,wigthez,epho_u,hbar,c,iepho,
!$omp& dyprop,dzprop,czero,npinyprop_u,npinzprop_u,zmin,ymin,aradprop_u,dzy,wigr,wigi,
!$omp& lyp,lym,ixy,kx,ky,expomz,expomy,zw,cw)

      ith=OMP_GET_THREAD_NUM()+1

!$omp DO

                      do jz=1,npinzprop_u

                        !expomz=expomz*dexpomz
                        zp=zmin+(jz-1)*dzprop
                        zzm=zw-zp
                        zzp=zw+zp
                        lzm=int((zzm-zmin)/dzprop)+1

                        if (lzm.lt.1.or.lzm.gt.npinzprop_u) cycle
                        lzp=int((zzp-zmin)/dzprop)+1
                        if (lzp.lt.1.or.lzp.gt.npinzprop_u) cycle
                        lobsm=lzm+(lym-1)*npinzprop_u
                        lobsp=lzp+(lyp-1)*npinzprop_u
                        lmph=lobsm+nobsvprop_u*(iepho-1)
                        lpph=lobsp+nobsvprop_u*(iepho-1)
c                        if(lpph.lt.1) then
c                          print*,kx,ky,jz,zw,zp
c                        endif
                        dcw=dconjg(aradprop_u(kx,lpph))*aradprop_u(ky,lmph)*expomz(jz)*expomy
c                        if (icwz.eq.0) dcw=dcw/2.0d0 !Simpson
c                        if (icwy.eq.0) dcw=dcw/2.0d0 !Simpson
                        cw(ith)=cw(ith)+dcw
                        icwz=1
                      enddo !jz
!$omp END DO
!$omp END PARALLEL
c                      cw=cw-dcw/2.0d0 !Simpson
                      expomy=expomy*dexpomy
                      icwy=1

                    enddo !jy

c                    cw=cw-dcw*3.0d0/4.0d0 !Simpson

                    do ith=1,mthreads_u
                      wigr(iz,iy,izth,iyth,iepho,ixy)=
     &                  wigr(iz,iy,izth,iyth,iepho,ixy)+dreal(cw(ith))*4.0d0*dzy
                      wigi(iz,iy,izth,iyth,iepho,ixy)=
     &                  wigi(iz,iy,izth,iyth,iepho,ixy)+dimag(cw(ith))*4.0d0*dzy
                    enddo
                  enddo !iy

                enddo !ix

              enddo !thez
            enddo !they

        enddo !nepho

      enddo !ixy

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase_wigner',0)

      return
      end
*CMZ :          05/01/2024  11.27.46  by  Michael Scheer
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
*KEND.

c NO WARRANTY

      implicit none
      real rn

      save

      call random_number(rn)

      return
      end
