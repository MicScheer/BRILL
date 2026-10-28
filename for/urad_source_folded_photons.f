*CMZ :          25/10/2026  12.48.34  by  Michael Scheer
*CMZ :  4.03/00 08/09/2026  12.41.50  by  Michael Scheer
*CMZ :  4.02/01 31/08/2026  07.37.54  by  Michael Scheer
*CMZ :  4.02/00 13/09/2025  10.16.17  by  Michael Scheer
*CMZ :  4.01/07 13/08/2024  10.11.51  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  10.49.56  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_source_folded_photons(nelec,ifieldprop,nepho,ndimpho,npho,
     &  photons,propweights,s0thrash)

      use uradphasemod

      implicit none

      real*8 depho,s0sum(npho),s0maxprop(npho),p,dy,dz,s0max(npho),s0thrash
      real :: photons(9*npho),propweights(npho)
      real s0ran(3)

      integer :: nelec,ifieldprop,nepho,ndimpho,npho,
     &  k,l,iy,iz,iepho,iobfr,iobs,kcount,maxcount=100000

      if (nepho.gt.1) then
        depho=epho_u(2)-epho_u(1)
      else
        depho=0.0d0
      endif

      npho=0
      propweights=0.0

      if (ifieldprop.ne.0) then

        s0sum=0.0d0
        s0maxprop=0.0d0
        s0max=0.0d0

        do iobs=1,nobsv_u
          do iepho=1,nepho_u
            iobfr=iobs+(iepho-1)*nobsv_u
            if(stokes_u(1,iobfr).gt.s0max(iepho)) s0max(iepho)=stokes_u(1,iobfr)
          enddo
        enddo

        do iobs=1,nobsvprop_u
          do iepho=1,nepho_u
            iobfr=iobs+(iepho-1)*nobsvprop_u
            if(stokesprop_u(1,iobfr).gt.s0maxprop(iepho)) s0maxprop(iepho)=stokesprop_u(1,iobfr)
            s0sum(iepho)=s0sum(iepho)+stokesprop_u(1,iobfr)
          enddo
        enddo

        k=1
        do l=1,ndimpho,9
          if (nepho.gt.1) then
            iepho=nint(photons_u(l)/depho)+1
          else
            iepho=1
          endif
          if (photons_u(l+5).lt.s0max(iepho)*s0thrash) then
            npho=npho+1
            k=l
            photons(k)=photons_u(k)
            photons(k+1:k+4)=photons_u(k+1:k+4)*1.0e3
            photons(k+5:k+8)=photons_u(k+5:k+8)/1.0e6
            cycle
          endif
          kcount=0
          do while (.true.)
            kcount=kcount+1
            if(kcount.gt.maxcount) then
              print*,"*** Warning in urad_source_folded_photons: More than 100000 tries to randomize photon exceeded"
              print*,"*** Skipping photon"
              exit
            endif
            call util_random(3,s0ran)
            iz=int(s0ran(1)*npinzprop_u)+1
            iy=int(s0ran(2)*npinyprop_u)+1
            iobs=(iy-1)*npinzprop_u+iz
            iobfr=iobs+(iepho-1)*nobsvprop_u
            dy=obsvyprop_u(iy)
            dz=obsvyprop_u(iz)
            p=stokesprop_u(1,iobfr)
            if (s0ran(3)*s0maxprop(iepho).le.p) then
              k=l
              photons(k)=photons_u(l)
              photons(k+1)=(photons_u(l+1)+dy)*1.0e3
              photons(k+2)=(photons_u(l+2)+dz)*1.0e3
              photons(k+3:k+4)=photons_u(l+3:l+4)*1.0e3
              photons(k+5:k+8)=photons_u(l+5:l+8)/1.0e6
              npho=npho+1
              propweights(npho)=1.0
              exit
            endif
          enddo !while
        enddo
      else
        propweights_u=1.0
        do k=1,ndimpho,9
          npho=npho+1
          photons(k)=photons_u(k)
          photons(k+1:k+4)=photons_u(k+1:k+4)*1.0e3
          photons(k+5:k+8)=photons_u(k+5:k+8)/1.0e6
        enddo
      endif
c      print*,"npho:",npho
      !allutil_break

      end
