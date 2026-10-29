*CMZ :          29/10/2026  09.03.40  by  Michael Scheer
*CMZ :  4.03/00 08/09/2026  12.41.50  by  Michael Scheer
*CMZ :  4.02/01 31/08/2026  07.37.54  by  Michael Scheer
*CMZ :  4.02/00 13/09/2025  10.16.17  by  Michael Scheer
*CMZ :  4.01/07 13/08/2024  10.11.51  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  10.49.56  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
*KEEP,WAVE_SEEDS.
      module waveseeds

      integer ::
     &  irnsize_wave=64, ifixseed_wave=-9999
     &  ,irnseed_wave(64)=[
     &  376577121
     & ,52147852
     & ,-1273034815
     & ,-1963249100
     & ,1195262240
     & ,-1718716574
     & ,-224354675
     & ,432587481
     & ,1692325775
     & ,1934175653
     & ,-107106772
     & ,648589804
     & ,-1919014861
     & ,-1763988460
     & ,-1039845022
     & ,1414926465
     & ,-1214705659
     & ,560082688
     & ,527470902
     & ,-1078636718
     & ,272932485
     & ,-356992740
     & ,-2013991490
     & ,-588501795
     & ,-1010120436
     & ,-1558306344
     & ,-1116776222
     & ,794926823
     & ,-1157173406
     & ,63711032
     & ,-1870802148
     & ,-674825931
     & ,-690546468
     & ,1671737514
     & ,-224394481
     & ,2026233226
     & ,-1141752469
     & ,2061158685
     & ,-1225625467
     & ,-147464566
     & ,1692325775
     & ,1934175653
     & ,-107106772
     & ,648589804
     & ,-1919014861
     & ,-1763988460
     & ,-1039845022
     & ,1414926465
     & ,-1214705659
     & ,560082688
     & ,527470902
     & ,-1078636718
     & ,272932485
     & ,-356992740
     & ,-2013991490
     & ,-588501795
     & ,-1010120436
     & ,-1558306344
     & ,-1116776222
     & ,794926823
     & ,-1157173406
     & ,63711032
     & ,-1870802148
     &  ,-674825931
     &  ]

      end module waveseeds
*KEEP,URADPHASEMOD.
      module uradphasemod

      double complex, dimension(:,:,:,:), allocatable :: fieldbunch
      double complex, dimension(:,:), allocatable :: arad_u,aradprop_u

      double precision, dimension(:,:), allocatable :: obsv_u, stokes_u,
     &  track_u, fbunch_u, obsvprop_u,stokesprop_u
      double precision, dimension(:), allocatable ::  epho_u,specpow_u,pow_u,
     &  obsvzprop_u,obsvyprop_u

      integer*8, dimension(:), allocatable :: nrad_u

      double precision
     &  ebeam_u,gamma_u,curr_u,banwid_u,
     &  xi_u,xe_u,yi_u,ye_u,zi_u,ze_u,step_u,
     &  pincen_u(3),pinw_u,pinh_u,
     &  ephmin_u,ephmax_u,emith_u,emitv_u,
     &  perlen_u,shift_u,beffv_u,beffh_u,pherror_u,phgshift_u,
     &  xbeta_u,betah_u,alphah_u,betav_u,alphav_u,espread_u,
     &  disph_u,dispph_u,dispv_u,disppv_u,bunchlen_u,bunchcharge_u,
     &  pinxprop_u,pinwprop_u,pinhprop_u,globphase_u,globphaseprop_u,
     &  wavzstart,wavystart,wavzamp,wavyamp,defl_u,harm_u

      real, dimension(:), allocatable :: photons_u,electrons_u,propweights_u

      integer nstep_u,nepho_u,nobsv_u,nbunch_u,npiny_u,npinz_u,ianalytic_u,nharm_u,
     &  nper_u,modeph_u,modepin_u,modesphere_u,noranone_u,nlpoi_u,ifixphase_u,ifold_u,
     &  nobsvprop_u,npinyprop_u,npinzprop_u,npinzo_u,npinyo_u,ifieldprop_u,ifieldsym_u

      integer nelecampgenpho_u,igenpho_u,ngenpho_u

      integer
     &  ibunch_u,ihbunch_u,mthreadso_u,mthreads_u,nelec_u,icohere_u,modebunch_u,modewave_u,
     &  kalloarad_u,kalloepho_u,kalloobsv_u,kalloaradprop_u,kallostokes_u,kallostokesprop_u

      end module uradphasemod
*KEND.

      subroutine urad_get_photons(
     &  nelec,noranone,
     &  perlen,shift,nper,beffv,beffh,
     &  nharm,harm,
     &  ebeam,curr,step,
     &  pinx,piny,pinz,pinw,pinh,npiny,npinz,modepin,
     &  nepho,ephmin,ephmax,
     &  xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  modeph,pherror,phgshift,modewave,
     &  ifieldprop,npinzprop,npinyprop,
     &  pinxprop,pinwprop,pinhprop,
     &  ndimpho,npho,photons,propweights,s0thresh)

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,PHYCONPARAM.
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

      complex*16 rea(3),expsh

      real*8
     &  perlen,shift,
     &  beffv,beffh,harm,
     &  ebeam,curr,step,
     &  pinx,piny,pinz,pinw,pinh,
     &  ephmin,ephmax,
     &  xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  pherror,phgshift,
     &  pinxprop,pinwprop,pinhprop,
     &  df,z,zmin,dz,y,ymin,dy,dyprop,dzprop,
     &  wlen1,rhv,perl,park,harm1,gamma,b0eff,s0thresh

      integer nelec,noranone,nper,npiny,npinz,modepin,nepho,modeph,ifieldprop,iepho,modewave,nharm,
     &  i,iz,iy,k,ierr,iobs,icbrill,ifrq,iobfr,ndimpho,nobsvprop,npinzprop,npinyprop,l,npho,
     &  iobrf

      integer :: igenpho=1

      integer ModeRan,ihbunch !All grid points (zero) or Nelec random grid points (not zero)

      real photons(9*npho),propweights(npho)

      real s0ran(3)

      real*8 :: SPECNOR_SI = 5.0007272079722076d+39

      moderan=1

      igenpho_u=igenpho

      ebeam_u=ebeam
      gamma_u=ebeam_u/emassg1

      if (modewave.eq.-1) mthreads_u=OMP_GET_MAX_THREADS()
      !mthreads_u=1
      !print*,"mthreads!!"

      npinzo_u=npinz
      npinyo_u=npiny

      nelec_u=nelec
      noranone_u=noranone
      icohere_u=0
      bunchlen_u=0.0d0
      bunchcharge_u=0.0d0

      if (modewave.eq.-1) ihbunch_u=0
      ihbunch=ihbunch_u

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

      perlen_u=perlen/1000.0d0
      shift_u=shift/1000.0d0
      nper_u=nper
      beffv_u=beffv
      beffh_u=beffh

      step_u=step/1000.0d0
      nstep_u=max(1,nint(perlen_u/step_u))

      curr_u=curr

      pincen_u=[pinx,piny,pinz]/1000.d0
      pinw_u=pinw/1000.0d0
      pinh_u=pinh/1000.0d0
      npiny_u=npiny
      npinz_u=npinz
      modepin_u=modepin

      ephmin_u=ephmin
      ephmax_u=ephmax
      banwid_u=0.001d0
      nepho_u=nepho

      npiny_u=max(1,npiny_u)
      npinz_u=max(1,npinz_u)

      phgshift_u=0.0d0

      nobsv_u=npiny_u*npinz_u
      nobsv_u=npiny_u*npinz_u

      deallocate(arad_u,stat=ierr)
      deallocate(epho_u,stat=ierr)
      deallocate(obsv_u,stat=ierr)
      deallocate(specpow_u,stat=ierr)
      deallocate(stokes_u,stat=ierr)
      deallocate(pow_u,stat=ierr)
      deallocate(nrad_u,stat=ierr)
      deallocate(photons_u,stat=ierr)
      deallocate(electrons_u,stat=ierr)

      allocate(arad_u(6,nobsv_u*nepho_u))
      allocate(epho_u(nepho_u))
      allocate(obsv_u(3,nobsv_u))
      allocate(specpow_u(nobsv_u))
      allocate(stokes_u(4,nobsv_u*nepho_u))
      allocate(pow_u(nobsv_u))
      allocate(nrad_u(nobsv_u+1))
      allocate(photons_u(9*npho))
      allocate(electrons_u(5*nelec_u))

      photons_u=0.0
      arad_u=(0.0d,0.0d0)
      specpow_u=0.0d0
      stokes_u=0.0d0
      pow_u=0.0d0
      nrad_u=0

      !allutil_break

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

      iobs=0
      y=ymin-dy
      do iy=1,npiny_u
        y=y+dy
        z=zmin-dz
        if (abs(y).le.1.0d-15) y=0.0d0
        do iz=1,npinz_u
          iobs=iobs+1
          z=z+dz
          if (abs(z).le.1.0d-15) z=0.0d0
          obsv_u(1,iobs)=pincen_u(1)
          obsv_u(2,iobs)=y
          obsv_u(3,iobs)=z
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
      emith_u=emith/1.0d9
      emitv_u=emitv/1.0d9
      disph_u=disph
      dispph_u=dispph
      dispv_u=dispv
      dispph_u=disppv
      espread_u=espread

      modeph_u=modeph
      pherror_u=pherror

      if (ifieldprop.ne.0) then
        nobsvprop=npinzprop*npinyprop
        npinzprop_u=npinzprop
        npinyprop_u=npinyprop
        nobsvprop_u=nobsvprop
        pinxprop_u=pinxprop
        pinwprop_u=pinwprop
        pinhprop_u=pinhprop
        deallocate(obsvzprop_u,obsvyprop_u,obsvprop_u,aradprop_u,stokesprop_u,stat=ierr)
        allocate(
     &    obsvzprop_u(npinzprop),obsvyprop_u(npinyprop),obsvprop_u(3,nobsvprop),
     &    aradprop_u(6,nobsvprop*nepho),stokesprop_u(4,nobsvprop*nepho))
        aradprop_u=(0.0d0,0.0d0)
        stokesprop_u=0.0d0
        !all util_break
        dyprop=pinhprop/dble(max(1,npinyprop-1))/1000.0d0
        dzprop=pinwprop/dble(max(1,npinzprop-1))/1000.0d0
        obsvyprop_u(1)=-pinhprop/2.0d0/1000.0d0
        do iy=2,npinyprop
          obsvyprop_u(iy)=obsvyprop_u(iy-1)+dyprop
          if (abs(obsvyprop_u(iy)).lt.1.0d-12) obsvyprop_u(iy)=0.0d0
        enddo
        obsvzprop_u(1)=-pinwprop/2.0d0/1000.0d0
        do iz=2,npinzprop
          obsvzprop_u(iz)=obsvzprop_u(iz-1)+dzprop
          if (abs(obsvzprop_u(iz)).lt.1.0d-12) obsvzprop_u(iz)=0.0d0
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
        !all util_break
        if (modepin.ne.0) then
          modepin_u=0
          nelec_u=1
          emith_u=0.0d0
          emitv_u=0.0d0
          espread_u=0.0d0
          Ifieldprop_u = 0
          ihbunch_u=0
          call urad_amprep(modewave)
          call urad_phase_prop_classic(mthreads_u)
          stokesprop_u=stokesprop_u*specnor_si/1.0d6
          modepin_u=1
          nelec_u=nelec
          emith_u=emith/1.0d9
          emitv_u=emitv/1.0d9
          espread_u=espread
        endif
      endif

      ihbunch_u=ihbunch

      call urad_amprep(modewave)

      stokes_u=stokes_u/1.0d6

      nelec=nelec_u

      call urad_source_folded_photons(nelec,ifieldprop,nepho,ndimpho,npho,photons,propweights,
     &  s0thresh)

      ihbunch_u=ihbunch

      end
*CMZ :          28/10/2026  09.50.30  by  Michael Scheer
*CMZ :  4.03/00 08/09/2026  12.41.50  by  Michael Scheer
*CMZ :  4.02/01 31/08/2026  07.37.54  by  Michael Scheer
*CMZ :  4.02/00 13/09/2025  10.16.17  by  Michael Scheer
*CMZ :  4.01/07 13/08/2024  10.11.51  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  10.49.56  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_source_folded_photons(nelec,ifieldprop,nepho,ndimpho,nphoin,
     &  photons,propweights,s0thresh)

      use uradphasemod

      implicit none

      real*8 depho,s0sum(nepho),s0meanprop(nepho),s0maxprop(nepho),p,dy,dz,s0max(nepho),s0thresh,
     &  s0fluxprop,s0fluxpin

      real :: photons(9*nphoin),propweights(nphoin)
      real s0ran(3)

      integer :: nelec,ifieldprop,nepho,ndimpho,npho,nphoin,nphothresh,
     &  k,l,iy,iz,iepho,iobfr,iobs,kcount,maxcount=100000,ks0max

      if (nepho.gt.1) then
        depho=epho_u(2)-epho_u(1)
      else
        depho=0.0d0
      endif

      npho=0
      propweights=0.0

      if (ifieldprop.ne.0) then

        s0meanprop=0.0d0
        s0sum=0.0d0
        s0maxprop=0.0d0
        s0max=0.0d0

c        do iobs=1,nobsv_u
c          do iepho=1,nepho_u
c            iobfr=iobs+(iepho-1)*nobsv_u
c            if(stokes_u(1,iobfr).gt.s0max(iepho)) s0max(iepho)=stokes_u(1,iobfr)
c          enddo
c        enddo

        do l=1,ndimpho,9
          if (nepho.gt.1) then
            iepho=nint(photons_u(l)/depho)+1
          else
            iepho=1
          endif
          if (photons_u(l+5).gt.s0max(iepho)) then
            ks0max=l
            s0max(iepho)=photons_u(l+5)
          endif
        enddo

        nphothresh=0
        do l=1,ndimpho,9
          if (nepho.gt.1) then
            iepho=nint(photons_u(l)/depho)+1
          else
            iepho=1
          endif
          if (photons_u(l+5).ge.s0max(iepho)*s0thresh) then
            nphothresh=nphothresh+1
            s0sum(iepho)=s0sum(iepho)+photons_u(l+5)
          endif
        enddo

c        print*,ks0max,photons_u(ks0max+5)

        do iobs=1,nobsv_u
          iepho=nepho/2+1
          iobfr=iobs+(iepho-1)*nobsv_u
          if (nrad_u(iobs).ne.0) then
            s0fluxpin=s0fluxpin+stokes_u(1,iobfr)/nrad_u(iobs)
          endif
        enddo

        do iobs=1,nobsvprop_u
          do iepho=1,nepho_u
            iobfr=iobs+(iepho-1)*nobsvprop_u
            if(stokesprop_u(1,iobfr).gt.s0maxprop(iepho)) s0maxprop(iepho)=stokesprop_u(1,iobfr)
            s0meanprop(iepho)=s0meanprop(iepho)+stokesprop_u(1,iobfr)
          enddo
        enddo

        s0meanprop=s0meanprop/nobsvprop_u

        s0fluxpin=s0fluxpin/nobsv_u*pinw_u*pinh_u
        s0fluxprop=s0meanprop(nepho/2+1)*pinwprop_u/1000.0d0*pinhprop_u/1000.0d0

c        stokesprop_u=stokesprop_u/s0fluxprop*s0fluxpin

        do iobs=1,nobsvprop_u
          do iepho=1,nepho_u
            iobfr=iobs+(iepho-1)*nobsvprop_u
            if(stokesprop_u(1,iobfr).gt.s0maxprop(iepho)) s0maxprop(iepho)=stokesprop_u(1,iobfr)
            s0meanprop(iepho)=s0meanprop(iepho)+stokesprop_u(1,iobfr)
          enddo
        enddo

        s0meanprop=s0meanprop/nobsvprop_u

c        if (ifieldprop.ne.0) then
c          print*,s0fluxpin,s0fluxprop,s0fluxpin/s0fluxprop
c        endif

c        print*,s0max,s0maxprop
c        print*,s0sum,s0meanprop

        k=1
        npho=0
        do l=1,ndimpho,9
          npho=npho+1
          k=l
          if (nepho.gt.1) then
            iepho=nint(photons_u(l)/depho)+1
          else
            iepho=1
          endif
          if (photons_u(l+5).lt.s0max(iepho)*s0thresh) then
            photons(k)=photons_u(k)
            photons(k+1:k+4)=photons_u(k+1:k+4)*1.0e3
            photons(k+5:k+8)=photons_u(k+5:k+8)/1.0e6
     &        *s0fluxpin/s0fluxprop
            cycle
          endif
c          kcount=0
c          do while (.true.)
c            kcount=kcount+1
c            if(kcount.gt.maxcount) then
c              print*,"*** Warning in urad_source_folded_photons: More than 100000 tries to randomize photon exceeded"
c              print*,"*** Skipping photon"
c              exit
c            endif
c          call util_random(3,s0ran)
            call util_random(2,s0ran)
            iz=int(s0ran(1)*npinzprop_u)+1
            iy=int(s0ran(2)*npinyprop_u)+1
            iobs=(iy-1)*npinzprop_u+iz
            iobfr=iobs+(iepho-1)*nobsvprop_u
            dy=obsvyprop_u(iy)
            dz=obsvyprop_u(iz)
            p=stokesprop_u(1,iobfr)
c           if (s0ran(3)*s0maxprop(iepho).le.p) then
              photons(k)=photons_u(l)
              photons(k+1)=(photons_u(l+1)+dy)*1.0e3
              photons(k+2)=(photons_u(l+2)+dz)*1.0e3
              photons(k+3:k+4)=photons_u(l+3:l+4)*1.0e3
              photons(k+5:k+8)=photons_u(l+5:l+8)/1.0e6
              propweights(npho)=p/s0meanprop(iepho)
c              exit
c            endif
c          enddo !while
        enddo
      else
        propweights_u=1.0
        do k=1,ndimpho,9
          npho=npho+1
          photons(k)=photons_u(k)
          photons(k+1:k+4)=photons_u(k+1:k+4)*1.0e3
          photons(k+5:k+8)=photons_u(k+5:k+8)/1.0e6
     &      *s0fluxpin/s0fluxprop
          propweights(npho)=1.0
        enddo
      endif
c      print*,"npho:",npho
      !allutil_break

      end
*CMZ :          24/10/2026  15.05.14  by  Michael Scheer
*CMZ :  4.03/00 07/10/2026  07.27.51  by  Michael Scheer
*CMZ :  4.02/01 02/09/2026  09.04.08  by  Michael Scheer
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

c     Calculates stokes_u, arad_u, and fbunch on the meter scale

      use omp_lib
      use uradphasemod

      implicit none

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
*KEEP,TRACK.
      INTEGER NCO,MCO,IBYONLY,ITRACK

      DOUBLE PRECISION PHIX,PHIMX,PHIMN
     &  ,XMX,XMN,YMX,YMN,ZMX,ZMN
     &  ,BXMX,BXMN,BYMX,BYMN,BZMX,BZMN
     &  ,XMXAE,XMNAE,YMXAE,YMNAE,ZMXAE,ZMNAE
     &  ,BXMXAE,BXMNAE,BYMXAE,BYMNAE,BZMXAE,BZMNAE
     &  ,BINT0,PINT,tint
     &  ,BINT1X,BINT1Y,BINT1Z,BINT2X,BINT2Y,BINT2Z,BINT3Y,BINT3YA
     &  ,DTIM0,DS0
     &  ,B2INTY
     &  ,WTRA2I,TTRA2I,HTRA2I,WTRA2IC,
     &  XTRACK,YTRACK,ZTRACK,
     &  VXTRACK,VYTRACK,VZTRACK

      COMMON/MINMAXC/
     &  PHIX,PHIMX,PHIMN,XMX,XMN,YMX,YMN,ZMX,ZMN
     &  ,BXMX,BXMN,BYMX,BYMN,BZMX,BZMN
     &  ,XMXAE,XMNAE,YMXAE,YMNAE,ZMXAE,ZMNAE
     &  ,BXMXAE,BXMNAE,BYMXAE,BYMNAE,BZMXAE,BZMNAE
     &  ,BINT0,PINT,tint
     &  ,BINT1X,BINT1Y,BINT1Z,BINT2X,BINT2Y,BINT2Z,BINT3Y,BINT3YA
     &  ,DTIM0,DS0
     &  ,WTRA2I,TTRA2I,HTRA2I,WTRA2IC
     &  ,B2INTY
     &  ,XTRACK,YTRACK,ZTRACK
     &  ,VXTRACK,VYTRACK,VZTRACK
     &  ,NCO,MCO,IBYONLY,ITRACK
*KEND.
cc+seq,uservar.

      complex*16 :: cde,czero=(0.0d0,0.0d0),ci=(0.0d0,1.0d0),cph00
     &  ,cjvsto(4,3)

      double precision :: h2,ddist,wlen,dphi,phase0

      double complex , dimension (:,:), allocatable :: aradbuff
      double complex , dimension (:,:,:), allocatable :: arad
      double precision, dimension (:), allocatable :: frq
      double precision, dimension (:,:), allocatable :: wsstokes,pow
      double precision, dimension (:,:,:,:,:), allocatable :: stokesprop
      double precision, dimension (:,:,:), allocatable :: fbunch,stokes

      complex*16, dimension (:), allocatable :: expphiran
      real, dimension (:), allocatable :: pherr,pherrc,phiran
      real, dimension(:,:), allocatable :: pranall,eall,photons

      real eran(6),pran(3),rr(2)

      double complex :: apol,amp0(6),damp(6),amp(6),ampn(6),zexp,
     &  apolh,apolr,apoll,apol45,stokesv(4,3),cero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double precision :: t,udgamtot,upow,vf0,vn,vx0,vx2,vxf0,vxi,vy0,vy2,vyf0,
     &  vyi,vz0,vz2,vzf0,vzi,wlen1,x0,x2,xf0,xi,xlell,y0,y2,yf0,yi,ypi,yy,yyp,
     &  z0,z2,zf0,zi,zpi,zz,zzp,fillb(41),stok1,stok2,stok3,stok4,speknor,
     &  sqnbunch,sqnphsp,specnor,specnor_si,root_specnor_si,sbnor,rpin,r00(3),xph0,
     &  r(3),r0(3),pw,ph,phsum,pkerr,ppin,parke,pc(3),pcbrill(3),om1,
     &  park,pr,hbarev,obs(3),om,fhigh,flow,gamma,eix,eiy,eiz,emassg,
     &  efx,efy,efz,eharm1,ecdipev,ebeam,dtpho,dt,dtelec,dd0,debeam,
     &  drn0(3),drn00(3),ds,dr0(3),dr00(3),drn(3),dpp,dph,dist,dist0,dobs(3),
     &  bunnor,clight,bunchx,beta,beff,spow,
     &  zp0,yp0,rph,a2,anor,fsum,smax,zob,yob,
     &  xkellip,zampell,yampell,parkv,parkh,zpampell,ypampell,emom,dzpin,dypin,zmin,ymin,phgsh,
     &  emitho,emitvo,rn(3)

      double precision xprop,yprop(npinyprop_u),zprop(npinzprop_u),dy,dz,pinwprop,pinhprop
      double complex, dimension(:,:,:,:,:), allocatable :: fprop
      !double complex, dimension(:,:,:), allocatable :: fpriv
      double complex :: fpriv(3,npinzprop_u,npinyprop_u)

      double complex, dimension (:), allocatable ::
     &  uampex,uampey,uampez,uampbx,uampby,uampbz
      double precision, dimension (:,:), allocatable :: utraxyz,ustokes

      double complex :: rea(3),expsh

      integer :: l,kfreq,iobsv,i,np2,nelec,mbunch,meinbunch,ibu,jbun,kmod,
     &  kran=6,icbrill,ilo,kobsv,i1,i2,n,
     &  ifail,ndimu,nstepu,ith,noespread,noemit,jbunch,jubunch,jhbunch,
     &  jcharge=-1,lmodeph,nclo,jeneloss=0,iamppin,
     &  iamppincirc=0,ifrob,iobfr,isub,jvelofield=0,nlbu=0,nepho,ielo,
     &  modewave,iepho,ifieldprop,nzprop,nyprop,im,izm,iym,ifix,lz,ly,
     &  lunbun,npho(1000),lpho,mpho

      integer, dimension (:), allocatable :: lnbunch

      integer :: idebug=0, lbunch=0, ierr=0, ielec=0,lele
      integer ibunch,ihbunch,mthreads,nobsv,nobsvo,iemit,noranone,iz,iy,ipz,ipy,nobsvz,nobsvy
      integer iobm,iobp,iobfrm,iobfrp,neleco,noranoneo
      integer :: ical=0

      save ical
c      integer iuser
c      iuser=user(3)

      emitho=emith_u
      emitvo=emitv_u

      npho=1

c      if (nelecampgenpho_u.ne.0) then
c        emitho=0.0d0
c        emitvo=0.0d0
c      endif

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

      if (abs(modepin_u).eq.1) then
        if (modepin_u.eq.1) then
          iamppin=3
        else
          iamppin=-3
        endif
        nobsv=1
      else
        iamppin=1
        nobsv=npiny_u*npinz_u
      endif

      icbrill=nobsvo/2+1

c      jhbunch=max(0,ihbunch)
      jhbunch=ihbunch
      meinbunch=nelec_u

      lbunch=0

c      if (jhbunch.ne.0) then
c        jhbunch=max(1,jhbunch)
c      endif
      !print*,size(photons_u)
      !allu

      if (nelec_u.lt.mthreads) mthreads=1

      if (igenpho_u.gt.0) then
        deallocate(photons,stat=ierr)
        mpho=size(photons_u)/mthreads
        allocate(photons(mpho,mthreads))
        photons=0.0
      endif

      nepho=nepho_u
      if (ifieldprop_u.eq.2) then
        allocate(
c     &    fpriv(3,npinzprop_u,npinyprop_u),
     &    fprop(3,npinzprop_u,npinyprop_u,nepho_u,mthreads),
     &    stokesprop(4,npinzprop_u,npinyprop_u,nepho_u,mthreads))
        stokesprop=0.0d0
      endif

c Array fieldbunch taken out. Therefor, modepin=1 and ifieldprop=1 not allowed anymore.
c To reinvent it, check +self,if=fieldbunch in CMZ-file wave.cmz

      call urad_field_ini(perlen_u,shift_u,beffv_u,beffh_u,modewave)

      if (perlen_u.ne.0.0d0) then

        emom=emasse1*dsqrt((gamma_u-1.0d0)*(gamma_u+1.0d0))

        xkellip=twopi1/perlen_u
        parkh=echarge1*dabs(beffh_u)*perlen_u/(twopi1*emasskg1*clight1)
        parkv=echarge1*dabs(beffv_u)*perlen_u/(twopi1*emasskg1*clight1)

        if (modewave.le.0) then
c*** OBSOLeTE, SEE z0= further down
          zampell=beffv_u*clight1/emom/xkellip**2
          yampell=beffh_u*clight1/emom/xkellip**2
c        zampell=zmx
c        yampell=ymx
c        print*,zpampell
          zpampell=parkh/gamma_u
          ypampell=parkv/gamma_u
c        zpampell=tan(phimx)
c        print*,zpampell
c        stop
        else
c          !all util_break
          yampell=(ymx-ymn)/2.0d0
          zampell=(zmx-zmn)/2.0d0
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

      if (modewave.le.0) then
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

      zi_u=z0
      yi_u=y0

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
c      dr0=[xf0-x0,yf0-y0,zf0-z0]
      dr0=[efx,efy,efz]*perlen_u
      r0=r0+dr0/2.0d0

      allocate(frq(nepho_u),
     &  uampex(nepho_u),uampey(nepho_u),uampez(nepho_u),
     &  uampbx(nepho_u),uampby(nepho_u),uampbz(nepho_u),
     &  utraxyz(14,ndimu),ustokes(4,nepho_u))

      allocate(pow(size(pow_u),mthreads))

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

      allocate(pherrc(nper_u),pherr(nper_u),
     &  arad(6,nepho_u*nobsv_u,mthreads),
     &  expphiran(max(1,nelec_u)))

      allocate(pranall(2,nelec_u))
      do i=1,nelec_u
        call util_random(3,pran)
        pranall(:,i)=pran(1:2)
        expphiran(i)=exp(dcmplx(0.0d0,twopi1*pran(3)))
      enddo

      !allutil_break
      if (nelec.gt.1.and.(ibunch.eq.0.or.
     &  emith_u.eq.0.0d0.and.emitv_u.eq.0.0d0.and.espread_u.eq.0.0d0.and.modewave.ge.0)
c     &    .or.
c     &    nelecampgenpho_u.ne.0
     &    ) then
        iemit=0
        ibunch=0
        nelec=1
        nelec_u=1
        print*,"*** Warning: Due to zero emittance and e-spread, Nelec set 1! ***"
      else
        iemit=1
      endif

      if (nelec_u.ne.nelec) then
        print*,''
        print*,'--- Warning in urad_amprep: Nelec adjusted to multiple of number of threads:',nelec_u
        print*,''
      endif

      if (modewave.eq.-1) then
        electrons_u=-9999.0
        electrons_u(1)=yi
        electrons_u(2)=zi
        electrons_u(3)=vyi/vxi
        electrons_u(4)=vzi/vxi
        electrons_u(5)=0.0
      endif

      lele=0
      if (iemit.ne.0) then

        allocate(eall(6,nelec_u))

        if (modebunch_u.eq.-1) then
          open(newunit=lunbun,file='urad_phase.ele',status='old')
          do i=1,nelec_u
            read(lunbun,*) eall(:,i)
          enddo
          close(lunbun)
        else
          do i=1,nelec_u
            xi=x0
            if (modepin_u.ne.2) then
              call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &          emith_u,emitv_u,
     &          disph_u,dispph_u,dispv_u,disppv_u,
     &          espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
            else
              ! espread only for folding procedure
              call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &          0.0d0,0.0d0,
     &          disph_u,dispph_u,dispv_u,disppv_u,
     &          espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
            endif
            eall(1,i)=xi-x0
            eall(2,i)=yi
            eall(3,i)=zi
            eall(4,i)=ypi
            eall(5,i)=zpi
            eall(6,i)=dpp
            if (modewave.ne.-1) then
              kmod=0
              if (i.eq.1) then
                kmod=1
              else
                if (jhbunch.ne.0) then
                  if(mod(i,iabs(jhbunch)).eq.0) then
                    kmod=1
                  endif
                endif
              endif
              if (kmod.ne.0) then
                electrons_u(lele+1:lele+1)=yi-ypi*x0
                electrons_u(lele+2:lele+2)=zi-zpi*x0
                electrons_u(lele+3:lele+3)=ypi
                electrons_u(lele+4:lele+4)=zpi
                electrons_u(lele+5:lele+5)=dpp
                lele=lele+5
              endif
            endif
          enddo
        endif

        if (noranone.ne.0) eall(:,1)=0.0

      endif

      !allocate(affe(6,nepho_u*nobsv))

      allocate(wsstokes(4,nepho_u*nobsvo),stokes(4,nepho_u*nobsvo,mthreads))
      stokes=0.0d0
      arad=(0.0d0,0.0d0)

      np2=nper_u/2

      call util_random_gauss_omp(nper_u,pherr,rr)
      pherrc=pherr

      lmodeph=modeph_u

      if (pherror_u.ne.0.0d0.and.(lmodeph.lt.0.or.lmodeph.gt.2)) then
        print*, ""
        print*, "*** Error in urad_amprep: MODEPH must be 0,1, or 2 ***"
        print*, "*** Program aborted ***"
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

c      if (ibunch.ne.0.and.bunchcharge_u.ne.0.0d0) then
c        sqnbunch=mbunch
c        sqnphsp=sqrt(bunchcharge_u/echarge1)
c     &    *meinbunch
c     &    /(bunchcharge_u/echarge1)
c        bunnor=1.0d0/mbunch
c      else
c        sqnbunch=mbunch
c        sqnphsp=sqrt(dble(nelec_u))
c        bunnor=1.0d0/mbunch
c      endif

      sqnbunch=1.0d0
      sqnphsp=1.0d0
      bunnor=1.0d0

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

c      specnor=
c     &  banwid_u
c     &  /(4.0d0*pi1**2*clight*hbarev)
c     &  /(4.0d0*pi1*eps01)
c     &  *curr_u
      specnor=1.0d0

      !Note: Amp0/sqrt(specnor_si) cancels the mm-scale,
      !but amp0*1000 converts everything to then m-scale
      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr_u ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid_u !BW
     &  /1.0d6 !m**2 > mm**2

c      sbnor=specnor*bunnor
c      speknor=specnor
c      sbnor=bunnor
      sbnor=bunnor*specnor_si
      root_specnor_si=sqrt(specnor_si)

c      anor=sqrt(sbnor/specnor_si)
c      anor=sqrt(1.0d0/specnor_si)
      anor=1.0d0

      jeneloss=0
      pw=pinw_u
      ph=pinh_u
      !pr=pinr
      pc=pincen_u
      do iobsv=1,nobsvo
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

      !allutil_break

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(ly,lz,kobsv,nepho,nobsvz,nobsvy,nobsv,nelec,frq,nper_u,np2,perlen_u,clight,hbarev,
!$OMP& ifieldprop,xprop,nyprop,nzprop,yprop,zprop,cjvsto,fpriv,pinwprop,pinhprop,
!$OMP& flow,fhigh,czero,cone,rea,expsh,ifix,zob,yob,
!$OMP& x0,y0,z0,xf0,yf0,zf0,vx0,vy0,vz0,vxf0,vyf0,vzf0,gamma_u,sbnor,speknor,bunnor,
!$OMP& efx,efy,efz,ds,ndimu,curr_u,xlell,parke,amp,ampn,amp0,
!$OMP& uampex,uampey,uampez,uampbx,uampby,uampbz,a2,anor,specnor_si,root_specnor_si,
!$OMP& lmodeph,zp0,yp0,modewave,
!$OMP& jbunch,jubunch,jhbunch,noespread,noemit,ebeam,
!$OMP& stokesv,icbrill,obsv_u,emassg,debeam,dispv_u,disppv_u,
!$OMP& betah_u,alphah_u,betav_u,alphav_u,emith_u,emitv_u,disph_u,dispph_u,
!$OMP& pran,pranall,eall,fillb,r0,dr0,iamppin,iamppincirc,pc,phase0,pr,banwid_u,
!$OMP& pw,ph,idebug,pcbrill,wsstokes,vn,bunchlen_u,modebunch_u,icohere_u,igenpho_u)
!$OMP& SHARED(mthreads,stokes,pherr,expphiran,lbunch,lnbunch,modepin_u,fieldbunch,npinzo_u,nobsvo,dzpin,dypin,
!$OMP& fbunch_u,jcharge,jeneloss,jvelofield,iemit,noranone,arad,pow,photons,
!$OMP& nrad_u,zmin,ymin,phgsh,fprop,stokesprop,npho,mpho)

      jbun=1
      isub=0
      iobsv=0
      ielo=0
      xph0=-perlen_u*dble(nper_u)/2.0d0

!$OMP DO

c      do ilo=1,nelec*nobsv
      do ielec=1,nelec

        !allutil_break

        do iobsv=1,nobsv

          wsstokes=0.0d0

          !affe=(0.0D0,0.0D0)
          spow=0.0d0

          ith=OMP_GET_THREAD_NUM()+1

          ! iobsv=mod(ilo-1,nobsv)+1
          ! ibu=(ilo-1)/nobsv+1
          ibu=ielec
          jbun=ibu

          iy=(iobsv-1)/nobsvz+1
          iz=mod(iobsv-1,nobsvz)+1

          !if (iz.gt.nobsvz/2+1) !allutil_break
          !allutil_break

          !ielec=ibu

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

          if (noranone.ne.0.and.ielec.eq.1.and.abs(iamppin).eq.3) then
            kobsv=icbrill
          else
            kobsv=iobsv
          endif

          !allutil_break
          obs=obsv_u(1:3,kobsv)

          if (noranone.eq.0.or.ielec.ne.1.or.kobsv.ne.icbrill) then

            if (abs(iamppin).eq.3) then
              !call util_random(2,pran)
              pran(1:2)=pranall(:,ielec)
              if (iamppincirc.eq.0) then
                obs(2)=pc(2)+(pran(1)-0.5)*ph
                obs(3)=pc(3)+(pran(2)-0.5)*pw
                iy=nint((obs(2)-ymin)/dypin)+1
                iz=nint((obs(3)-zmin)/dzpin)+1
                kobsv=iz+nobsvz*(iy-1)
                obs=obsv_u(:,kobsv)
              else
                rpin=(pran(1)-0.5)*pr
                ppin=pran(2)*twopi1
                obs(2)=pc(2)+rpin*cos(ppin)
                obs(3)=pc(3)+rpin*sin(ppin)
                iy=nint((obs(2)-ymin)/dypin)+1
                iz=nint((obs(3)-zmin)/dzpin)+1
                kobsv=iz+nobsvz*(iy-1)
                obs=obsv_u(:,kobsv)
              endif
            endif

            if (iamppin.eq.-3) then

              ly=1+nint((obs(2)-(pc(2)-ph/2.0d0))/dypin)
              lz=1+nint((obs(3)-(pc(3)-pw/2.0d0))/dzpin)
              kobsv=lz+(ly-1)*npinzo_u

              obs(2)=ymin+(ly-1)*dypin
              obs(3)=zmin+(lz-1)*dzpin

              if (abs(obs(2)).lt.1.0d-12) obs(2)=0.0d0
              if (abs(obs(3)).lt.1.0d-12) obs(3)=0.0d0

            endif

          endif

          vn=norm2([vxi,vyi,vzi])
          eix=vxi/vn
          eiy=vyi/vn
          eiz=vzi/vn

          h2=((obs(2)-yi)**2+(obs(3)-zi)**2)/(obs(1)-xph0)**2
          if (h2.lt.0.01) then
            rph=abs(obs(1)-xph0)*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &        -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
          else
            rph=sqrt((obs(1)-xph0)**2+((obs(2)-yi)**2+(obs(3)-zi)**2))
          endif

          phase0=(rph-(obsv_u(1,icbrill)-xph0))/clight

          call urad_e_b_field(
     &      jcharge,curr_u,
     &      gamma,udgamtot,
     &      xi,yi,zi,vxi,vyi,vzi,
     &      xf0,yf0,zf0,efx,efy,efz,
     &      x2,y2,z2,vx2,vy2,vz2,dtelec,ds,
     &      0,nstepu,ndimu,utraxyz,phase0,
     &      obs(1),obs(2),obs(3),flow,fhigh,
     &      nepho,frq,uampex,uampey,uampez,uampbx,uampby,uampbz,
     &      ustokes,upow,
     &      jeneloss,jvelofield,ifail,ith,banwid_u,modewave
     &      )

          r0=[xi,yi,zi]
          dr0=[x2-xi,y2-yi,z2-zi]

          drn=dr0/norm2(dr0)
          r0=r0+dr0/2.0d0

          do kfreq=1,nepho

            !allutil_break
            iobfr=kobsv+nobsvo*(kfreq-1)

            om=frq(kfreq)/hbarev

c          if (modewave.eq.0) then
c            amp0=[
c     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
c     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
c     &        ]*1.0d3/sqrt(speknor/curr_u*0.10d0) !urad
c          else
            !Note: (uampex**2 = S0)
            !Amp0/sqrt(specnor_si) * 1000 converts everything to m-scale
            amp0=[
     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
c     &        ]*1.0d3/sqrt(speknor) !urad
     &        ] / root_specnor_si *1000.0d0
c          endif

c          call util_random(1,pran)
c          amp0=amp0*dcmplx(0.0d0,dble(pran(1)*twopi1))

            amp=(0.0d0,0.0d0)
            t=bunchx/vn

            !allutil_break
            do i=1,nper_u

              !if (i.eq.1.or.i.eq.nper_u) !allutil_break

              r=r0+(i-np2-1)*dr0
              dobs=obs-r
              dist0=norm2(obs-r0)
              dist=norm2(dobs)

              if (kfreq.eq.1) then
                !allu
                spow=spow+upow*(dist0/dist)**2
                pow(kobsv,ith)=pow(kobsv,ith)+upow*(dist0/dist)**2
              endif

              if (lmodeph.eq.0) then
!!!!!                dt=xlell/clight*((1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
!!!!!     &            (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)
                h2=
     &            (ypi-yp0-dobs(2)/dobs(1))**2 +
     &            (zpi-zp0-dobs(3)/dobs(1))**2
c26.4.2024     &          ((ypi-yp0-dobs(2))/dobs(1))**2 +
c26.4.2024     &          ((zpi-zp0-dobs(3))/dobs(1))**2

c              dt=xlell/clight*
c     &          (
c     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+h2/2.0d0-h2**2/8.0d0
c     &          )

                dph=om*(t+pherr(i))

                dt=xlell/clight*
     &            (
     &            (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
     &            (((((-0.0205078125D0*h2+0.02734375D0)*h2
     &            -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2
     &            )

                t=t+dt
              else if (lmodeph.eq.1.or.lmodeph.eq.2) then
                dph=om*t
                pkerr=parke*(1.0d0+pherr(i))
!!!!!                dt=xlell/clight*((1.0d0+pkerr**2/2.0d0)/2.0d0/gamma**2+
!!!!!     &            (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)
                h2=
     &            (ypi-yp0-dobs(2)/dobs(1))**2 +
     &            (zpi-zp0-dobs(3)/dobs(1))**2
c26.4.2024              h2=((ypi-yp0-dobs(2))**2+(zpi-zp0-dobs(3))**2)/dobs(1)**2
                dt=xlell/clight*
     &            (
c25.4.2024     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
     &            (1.0d0+pkerr**2/2.0d0)/2.0d0/gamma**2+
     &            (((((-0.0205078125D0*h2+0.02734375D0)*h2
     &            -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2
     &            )
                t=t+dt
              endif !lmodeph

              zexp=cdexp(dcmplx(0.0d0,dph))
              damp=amp0*zexp*dist0/dist
              amp=amp+damp

              if (i.eq.nper_u) then
                if (ifix.ne.0) then
                  amp=amp*expphiran(ielec)
                endif
              endif

              if (jhbunch.ne.0) then

c4.8.2026              if (
c4.8.2026     &            ((iamppin.eq.3.and.kobsv.eq.icbrill).and.jhbunch.gt.0.and.
c4.8.2026     &            mod(ielec,jhbunch).eq.0) .or.
c4.8.2026     &            (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

                if (
     &              ((iamppin.eq.3.or.kobsv.eq.icbrill).and.jhbunch.gt.0.and.
     &              mod(ielec,jhbunch).eq.0) .or.
     &              (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

                  if (i.eq.1) then
                    !allutil_break
                    fillb(5)=r(1)
                    fillb(6)=r(2)
                    fillb(7)=r(3)
                    fillb(8)=ypi
                    fillb(9)=zpi
                  else if (i.eq.nper_u) then

C Bug, if nobsv.gt.1, since phase not globally shifted
c                  if (abs(phgsh).eq.9999.0d0) then
c                    rea=amp(1:3)
c                    expsh=rea(3)/abs(rea(3))
c                    if (phgsh.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
c                    amp=amp/expsh
c                  else if (phgsh.ne.0.0d0) then
c                    expsh=cdexp(dcmplx(0.0d0,phgsh)) !*1.0d3
c                    amp=amp/expsh
c                  endif
c
c                  if (phgsh.eq.9999.0d0) then
c                    rea=amp(1:3)
c                    expsh=rea(3)/abs(rea(3))
c                    amp=amp/expsh
c                  else if (phgsh.eq.-9999.0d0) then
c                    rea=amp(1:3)
c                    expsh=rea(3)/abs(rea(3))*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
c                    amp=amp/expsh
c                  else if (phgsh.ne.0.0d0) then
c                    expsh=cdexp(dcmplx(0.0d0,phgsh)) !*1.0d3
c                    amp=amp/expsh
c                  endif

c                  a2=sum(abs(amp)**2)
c                  if (a2.gt.0.0d0) then
c                    stok1=a2*sbnor
c                    print*,"S0:",stok1
c                    anor=sqrt(stok1/a2/specnor_si)
c                    anor=sqrt(sbnor/specnor_si)
c                  else
c                    anor=1.0d0
c                  endif

c                  ampn=amp*anor/1000.0d0
c                  ampn=amp/1000.0d0
                    ampn=amp
c                  print*,"Probe:",
c     &              sum(abs(ampn)**2)*specnor_si,
c     &              ampn(3)

                    fillb(10:12)=r+dr0
                    fillb(13)=ypi
                    fillb(14)=zpi
                    fillb(30)=dreal(ampn(1))
                    fillb(31)=dimag(ampn(1))
                    fillb(32)=dreal(ampn(2))
                    fillb(33)=dimag(ampn(2))
                    fillb(34)=dreal(ampn(3))
                    fillb(35)=dimag(ampn(3))
                    fillb(36)=dreal(ampn(4))
                    fillb(37)=dimag(ampn(4))
                    fillb(38)=dreal(ampn(5))
                    fillb(39)=dimag(ampn(5))
                    fillb(40)=dreal(ampn(6))
                    fillb(41)=dimag(ampn(6))

                  endif

                endif

              endif !jhbunch.ne.0

            enddo !nper_u

            !allutil_break

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

            !if (kobsv.eq.icbrill) !all util_break
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

            !allu

            if (sum(abs(amp)).gt.0.0d0) amp=amp/sum(abs(amp))

            rn(1)=
     &        real(amp(2)*conjg(amp(6))-amp(3)*conjg(amp(5)))
            rn(2)=
     &            real(amp(3)*conjg(amp(4))-amp(1)*conjg(amp(6)))
            rn(3)=
     &        real(amp(1)*conjg(amp(5))-amp(2)*conjg(amp(4)))

            rn(2)=rn(2)/rn(1)
            rn(3)=rn(3)/rn(1)
            rn(1)=sqrt(1.0d0-(rn(2)**2+rn(3)**2))

            !allu

            if (modewave.lt.0.and.ifieldprop.eq.0.or.igenpho_u.gt.0) then

              kmod=0

              if (ielec.eq.1.or.modewave.lt.0) then
                kmod=1
              else
                if (jhbunch.ne.0) then
                  if(mod(ielec,iabs(jhbunch)).eq.0) then
                    kmod=1
                  endif
                endif
              endif

              if (kmod.ne.0) then
                if (npho(ith).le.mpho-8) then
                  photons(npho(ith),ith)=frq(kfreq)
                  photons(npho(ith)+1,ith)=obs(2)
                  photons(npho(ith)+2,ith)=obs(3)
                  photons(npho(ith)+3,ith)=rn(2)
                  photons(npho(ith)+4,ith)=rn(3)
                  photons(npho(ith)+5,ith)=wsstokes(1,iobfr)
                  photons(npho(ith)+6,ith)=wsstokes(2,iobfr)
                  photons(npho(ith)+7,ith)=wsstokes(3,iobfr)
                  photons(npho(ith)+8,ith)=wsstokes(4,iobfr)
                  npho(ith)=npho(ith)+9
                endif
              endif

            endif

            if (kfreq.eq.1) then
              nrad_u(kobsv)=nrad_u(kobsv)+1
              !allu
c              print*,ielec,kobsv,iobfr
              if (iobfr.eq.19) then
                !print*,ielec,kobsv, wsstokes(1,iobfr),nradth(kobsv,ith),stokes(1,iobfr,ith)
                !allutil_break
              endif
            endif

c          if (
c     &        ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
c     &        mod(ielec,jhbunch).eq.0) .or.
c     &        (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then
c            if (ielec.gt.4) then
c              print*,jhbunch,ith,ilo,ielec
c            endif
c          endif

            if (jhbunch.ne.0.and.modewave.eq.0) then
              if (
     &            ((iamppin.eq.3.or.kobsv.eq.icbrill).and.jhbunch.gt.0.and.
     &            mod(ielec,jhbunch).eq.0) .or.
     &            (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

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

                ! 29.8.2026 fillb(22)=wsstokes(1,iobfr)*nelec

                ! 29.8.2026 fillb(23)=wsstokes(1,iobfr)*nelec
                ! 29.8.2026 fillb(24)=wsstokes(2,iobfr)*nelec
                ! 29.8.2026 fillb(25)=wsstokes(3,iobfr)*nelec
                ! 29.8.2026 fillb(26)=wsstokes(4,iobfr)*nelec

                fillb(22)=wsstokes(1,iobfr)

                fillb(23)=wsstokes(1,iobfr)
                fillb(24)=wsstokes(2,iobfr)
                fillb(25)=wsstokes(3,iobfr)
                fillb(26)=wsstokes(4,iobfr)

                fillb(27)=spow
                fillb(28)=1
                fillb(29)=dtelec

c              anor=sqrt(sbnor/specnor_si)
c              ampn=amp*anor/1000.0d0
c              ampn=amp/1000.0d0
                ampn=amp

                fillb(30)=dreal(ampn(1))
                fillb(31)=dimag(ampn(1))
                fillb(32)=dreal(ampn(2))
                fillb(33)=dimag(ampn(2))
                fillb(34)=dreal(ampn(3))
                fillb(35)=dimag(ampn(3))
                fillb(36)=dreal(ampn(4))
                fillb(37)=dimag(ampn(4))
                fillb(38)=dreal(ampn(5))
                fillb(39)=dimag(ampn(5))
                fillb(40)=dreal(ampn(6))
                fillb(41)=dimag(ampn(6))
                lbunch=lbunch+1
                fbunch_u(:,lbunch)=fillb(:)
              endif !fill

            endif !jhbunch

            if (ifieldprop.eq.2) then
              if (kobsv.eq.1) then
                fprop(1:3,1:nzprop,1:nyprop,kfreq,ith)=(0.0d0,0.0d0)
              endif
              call urad_phase_prop_point(obs,amp(1:3),nzprop,nyprop,
     &          xprop,yprop,zprop,pinwprop,pinhprop,frq(kfreq),fpriv)
              fprop(:,:,:,kfreq,ith)=fprop(:,:,:,kfreq,ith)+fpriv(:,:,:)
              if (kobsv.eq.nobsv) then
                i=0
                do ipy=1,nyprop
                  do ipz=1,nzprop
                    i=i+1+nzprop*nyprop*(kfreq-1)
                    apolh=
     &                fprop(1,ipz,ipy,kfreq,ith)*cjvsto(1,1)+
     &                fprop(2,ipz,ipy,kfreq,ith)*cjvsto(1,2)+
     &                fprop(3,ipz,ipy,kfreq,ith)*cjvsto(1,3)

                    apolr=
     &                fprop(1,ipz,ipy,kfreq,ith)*cjvsto(2,1)+
     &                fprop(2,ipz,ipy,kfreq,ith)*cjvsto(2,2)+
     &                fprop(3,ipz,ipy,kfreq,ith)*cjvsto(2,3)

                    apoll=
     &                fprop(1,ipz,ipy,kfreq,ith)*cjvsto(3,1)+
     &                fprop(2,ipz,ipy,kfreq,ith)*cjvsto(3,2)+
     &                fprop(3,ipz,ipy,kfreq,ith)*cjvsto(3,3)

                    apol45=
     &                fprop(1,ipz,ipy,kfreq,ith)*cjvsto(4,1)+
     &                fprop(2,ipz,ipy,kfreq,ith)*cjvsto(4,2)+
     &                fprop(3,ipz,ipy,kfreq,ith)*cjvsto(4,3)

                    stok1=dreal(
     &                apolr*conjg(apolr)+
     &                apoll*conjg(apoll))

                    stok2=-stok1+
     &                dreal(2.*apolh*conjg(apolh))

                    stok3=
     &                dreal(2.*apol45*conjg(apol45))-
     &                stok1

                    stok4=dreal(
     &                apolr*conjg(apolr)-
     &                apoll*conjg(apoll))

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

      !allu
      lpho=0
      do ith=1,mthreads
        if (modewave.lt.0.or.igenpho_u.gt.0) then
          mpho=npho(ith)-1
          photons_u(lpho+1:lpho+mpho)=photons(1:mpho,ith)
          lpho=lpho+mpho
        endif
        pow_u(:)=pow_u(:)+pow(:,ith)
        arad_u(:,:)=arad_u(:,:)+arad(:,:,ith)
      enddo


      if (globphase_u.eq.9999.0d0) then
        do kfreq=1,nepho
          iobfr=icbrill+nobsvo*(kfreq-1)
          cph00=arad_u(3,iobfr)/abs(arad_u(3,iobfr))
          do iobsv=1,nobsv
            iobfr=iobsv+nobsvo*(kfreq-1)
            arad_u(:,iobfr)=arad_u(:,iobfr)/cph00
          enddo
        enddo
      else if (globphase_u.ne.0.0d0) then
        arad_u=arad_u*exp(ci*globphase_u)
      endif

      if (sqnbunch.ne.1.0d0) pow_u=pow_u/sqnbunch

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

      !allutil_break

      if (icohere_u.eq.0) then

        if (sqnbunch.ne.1.0d0) arad_u=arad_u/sqnbunch
        do ith=1,mthreads
          stokes_u(:,:)=stokes_u(:,:)+stokes(:,:,ith)
          !print*,stokes_u(1,19)
        enddo

      else

        do iobsv=1,nobsvo
          do kfreq=1,nepho

            iobfr=iobsv+nobsvo*(kfreq-1)

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
     &  pherrc,pherr,expphiran,arad,pow,pranall,wsstokes,stokes,
     &  )
      if (iemit.ne.0) deallocate(eall)

c      iobfr=nobsv_u*nepho_u/2+1
c
c      !all util_break
c      if (nrad(iobfr).eq.0) then
c        do iobsv=1,nobsvo
c          do kfreq=1,nepho
c            iobfr=iobsv+(kfreq-1)*nobsvo
c            if (nrad(iobfr).ne.0) then
c              iobfr=-iobfr
c              exit
c            endif
c          enddo
c          if (iobfr.lt.0) then
c            iobfr=-iobfr
c            exit
c          endif
c        enddo
c      endif

c      amp(1:3)=arad_u(1:3,iobfr)
c      anor=sqrt(stokes_u(1,iobfr)/
c     &  sum(abs(amp)**2)/specnor_si)
c      anor=sqrt(sbnor/specnor_si)

c      arad_u=arad_u*anor

c      if (ifieldprop_u.eq.2) then
      if (ifieldprop_u.gt.0) then
        stokesprop_u=stokesprop_u*specnor_si
      endif

c      print*,"Nrad_Max, iz,iy,iobsv:",nradmax,izradmax,iyradmax,iobsvradmax

c      print*,specnor_si
c      print'(5(1pe12.3))',stokes_u(1,1),abs(arad_u(3,1))**2*specnor_si

9999  ical=1

c      if (nelecampgenpho_u.ne.0) then
c        emitho=emith_u
c        emitvo=emitv_u
c      endif

      end
*CMZ :  4.01/02 09/05/2023  13.11.31  by  Michael Scheer
*CMZ :  4.01/00 10/02/2023  13.52.47  by  Michael Scheer
*-- Author :    Michael Scheer   10/02/2023
      subroutine urad_field_ini(perl,shift,b0v,b0h,modewave)

      implicit none

*KEEP,AMPLI.

      INTEGER IAMPDIMP
      PARAMETER (IAMPDIMP=10000)

      INTEGER IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD,modeph
     &  ,IAMPOBSV,IAMPCOMP,IAMPREP,IAMPSUP,IAMPSEED,iamppin,iamppincirc
     &  ,iampcoh,iampincoh,nampbunchharm,noemitph,noespreadph,mampthreads

      DOUBLE PRECISION AMPFREQ,ampr2corr,AMPPHI(IAMPDIMP),AMPSHIFT(IAMPDIMP)
     &  ,AMPSCALE(IAMPDIMP),AMPRAN,FACAMPLI,phrERROR,
     &  phrshift,phrb0h,phrb0v,phrperl,phdx,phdy,phdz,phrxbeta,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv,
     &  ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen

      COMMON/AMPLIC/
     &  ampfreq,ampr2corr,AMPPHI,AMPSHIFT,AMPSCALE,AMPRAN,FACAMPLI,phrERROR
     &  ,ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56,phrxbeta,
     &  phrshift,phrb0h,phrb0v,phrperl,phdx,phdy,phdz,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen,
     &  IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD,modeph
     &  ,IAMPOBSV,IAMPCOMP,IAMPREP,IAMPSUP,IAMPSEED,iamppin,iamppincirc
     &  ,iampcoh,iampincoh,nampbunchharm,noemitph,noespreadph,mampthreads

      NAMELIST/AMPLIN/
     &  IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD
     & ,IAMPCOMP,IAMPREP
     & ,AMPSHIFT,AMPSCALE,ampfreq,ampr2corr,AMPPHI,IAMPSUP,AMPRAN,IAMPSEED
     &  ,ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56
     &  ,iampcoh,iampincoh,nampbunchharm

      NAMELIST/PHASEREPN/phrERROR,noemitph,noespreadph,modeph,phrxbeta,
     &  phrshift,phrperl,phrb0h,phrb0v,phdx,phdy,phdz,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv
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

*KEEP,GPLHINT.
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
*CMZ :  4.01/04 17/11/2023  11.59.33  by  Michael Scheer
*CMZ :  4.00/17 04/10/2022  08.10.22  by  Michael Scheer
*CMZ :  4.00/11 27/05/2021  09.41.25  by  Michael Scheer
*CMZ :  3.02/04 03/12/2014  15.11.16  by  Michael Scheer
*-- Author :    Michael Scheer   03/12/2014
      subroutine util_break
*KEEP,DEBUGWAVE.
      double precision x_debug,y_debug,z_debug,a_debug(100)
      integer i_debug,k_debug

      common/c_debug/x_debug,y_debug,z_debug,a_debug,i_debug,k_debug
*KEND.
      return
      end
*CMZ :  4.01/05 15/04/2024  15.07.21  by  Michael Scheer
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
*CMZ :  3.03/02 29/02/2016  16.24.34  by  Michael Scheer
*CMZ :  3.02/03 05/09/2014  12.29.51  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random(n,r)
*KEEP,RANDOM.
      integer*8 irancalls
      integer irnseed(64),irnmode,irnsize,irnseedi(64)
      common /randomc/ irancalls,irnseed,irnmode,irnsize,irnseedi
      namelist /randomn/ irnmode,irnseed
*KEND.

      real r(n)

      call random_number(r)

      irancalls=irancalls+n

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
*CMZ :  4.02/01 12/05/2026  08.28.47  by  Michael Scheer
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
          WRITE(6,*)'X:',x
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
*CMZ :  4.01/07 29/11/2024  18.43.16  by  Michael Scheer
*CMZ :  4.01/05 31/10/2022  17.05.45  by  Michael Scheer
*CMZ : 00.00/16 19/03/2014  12.30.26  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.27.13  by  Michael Scheer
*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_PARABEL(Xin,Yin,A,YP,XOPT,yopt,IFAIL)

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
*CMZ :  4.00/11 28/05/2021  09.19.32  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.18.57  by  Michael Scheer
*CMZ :  2.68/05 03/09/2012  09.27.27  by  Michael Scheer
*CMZ : 00.00/11 11/02/2011  15.34.09  by  Michael Scheer
*-- Author :    Michael Scheer   11/02/2011
      SUBROUTINE util_spline_running_integral(X,Y,N,RESULT
     &                                 ,COEF,WORK1,WORK2,WORK3,WORK4)
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

*KEEP,AMPLI.

      INTEGER IAMPDIMP
      PARAMETER (IAMPDIMP=10000)

      INTEGER IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD,modeph
     &  ,IAMPOBSV,IAMPCOMP,IAMPREP,IAMPSUP,IAMPSEED,iamppin,iamppincirc
     &  ,iampcoh,iampincoh,nampbunchharm,noemitph,noespreadph,mampthreads

      DOUBLE PRECISION AMPFREQ,ampr2corr,AMPPHI(IAMPDIMP),AMPSHIFT(IAMPDIMP)
     &  ,AMPSCALE(IAMPDIMP),AMPRAN,FACAMPLI,phrERROR,
     &  phrshift,phrb0h,phrb0v,phrperl,phdx,phdy,phdz,phrxbeta,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv,
     &  ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen

      COMMON/AMPLIC/
     &  ampfreq,ampr2corr,AMPPHI,AMPSHIFT,AMPSCALE,AMPRAN,FACAMPLI,phrERROR
     &  ,ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56,phrxbeta,
     &  phrshift,phrb0h,phrb0v,phrperl,phdx,phdy,phdz,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen,
     &  IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD,modeph
     &  ,IAMPOBSV,IAMPCOMP,IAMPREP,IAMPSUP,IAMPSEED,iamppin,iamppincirc
     &  ,iampcoh,iampincoh,nampbunchharm,noemitph,noespreadph,mampthreads

      NAMELIST/AMPLIN/
     &  IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD
     & ,IAMPCOMP,IAMPREP
     & ,AMPSHIFT,AMPSCALE,ampfreq,ampr2corr,AMPPHI,IAMPSUP,AMPRAN,IAMPSEED
     &  ,ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56
     &  ,iampcoh,iampincoh,nampbunchharm

      NAMELIST/PHASEREPN/phrERROR,noemitph,noespreadph,modeph,phrxbeta,
     &  phrshift,phrperl,phrb0h,phrb0v,phdx,phdy,phdz,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv
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
*CMZ :  4.03/00 13/09/2026  08.38.55  by  Michael Scheer
*CMZ :  4.02/01 24/09/2025  12.48.34  by  Michael Scheer
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

c13.9.2026      da=pinw_u/1000.0d0*pinh_u/1000.0d0/dble(max(1,npinz_u-1)*max(1,npiny_u-1))
      da=pinw_u*pinh_u/dble(max(1,npinz_u-1)*max(1,npiny_u-1))

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

c13.9.2026          XOBS=OBSV_u(1,IOBS)/1000.0d0
c13.9.2026          YOBS=OBSV_u(2,IOBS)/1000.0d0
c13.9.2026          ZOBS=OBSV_u(3,IOBS)/1000.0d0

          XOBS=OBSV_u(1,IOBS)
          YOBS=OBSV_u(2,IOBS)
          ZOBS=OBSV_u(3,IOBS)

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
