*CMZ :          02/04/2024  12.41.56  by  Michael Scheer
*CMZ :  4.01/00 14/01/2023  21.16.25  by  Michael Scheer
*-- Author :    Michael Scheer   14/01/2023
*KEEP,URADPHASEMOD.
      module uradphasemod

      double complex, dimension(:,:,:,:), allocatable :: fieldbunch
      double complex, dimension(:,:), allocatable :: arad_u,aradprop_u

      double precision, dimension(:,:), allocatable :: obsv_u, stokes_u,
     &  track_u, fbunch_u, obsvprop_u,stokesprop_u
      double precision, dimension(:), allocatable ::  epho_u,specpow_u,pow_u

      double precision
     &  ebeam_u,gamma_u,curr_u,banwid_u,
     &  xi_u,xe_u,yi_u,ye_u,zi_u,ze_u,step_u,
     &  pincen_u(3),pinw_u,pinh_u,
     &  ephmin_u,ephmax_u,emith_u,emitv_u,
     &  perlen_u,shift_u,beffv_u,beffh_u,pherror_u,
     &  xbeta_u,betah_u,alphah_u,betav_u,alphav_u,espread_u,
     &  disph_u,dispph_u,dispv_u,disppv_u,bunchlen_u,bunchcharge_u,
     &  pinxprop_u,pinwprop_u,pinhprop_u

      integer nstep_u,nepho_u,nobsv_u,nbunch_u,npiny_u,npinz_u,
     &  nper_u,modeph_u,modepin_u,modesphere_u,noranone_u,nlpoi_u,
     &  nobsvprop_u,npinyprop_u,npinzprop_u,npinzo_u,npinyo_u

      integer
     &  ibunch_u,ihbunch_u,mthreads_u,nelec_u,icohere_u,modebunch_u

      end module uradphasemod
*KEEP,WIGNERMOD.
      module wignermod

      double precision, dimension(:,:,:,:,:,:), allocatable :: wigr,wigi
      double precision, dimension(:), allocatable :: wigthez,wigthey,wigz,wigy
      double precision :: dzprop,dyprop,dzwig,dywig,pinwwig,pinhwig

      integer nywig,nzwig,nzthewig,nythewig,nzfringe,nyfringe

      end module wignermod
