*CMZ :  2.00/00 06/12/2017  16.27.22  by  Michael Scheer
*CMZ :  1.02/01 08/04/2015  15.32.14  by  Michael Scheer
*CMZ :  1.02/00 07/10/2013  15.31.10  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.06.36  by  Michael Scheer
*CMZ :  1.01/00 11/09/2013  14.04.22  by  Michael Scheer
*CMZ :  1.00/03 06/09/2013  11.51.22  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.18.24  by  Michael Scheer
*CMZ :  1.00/00 04/09/2013  15.16.56  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  13.16.25  by  Michael Scheer
*CMZ :  0.00/02 27/06/2011  16.48.57  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.27.25  by  Michael Scheer
*-- Author :    Andreas Gaupp
	program brilliance
c	plot of central brilliance for undulators, wigglers and bends

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEEP,datetime.
      character(10) dtday,dttime,dtzone
      integer idatetime(8)
      common/datec/idatetime,dtday,dttime,dtzone
*KEND.
        integer nmenuep,imenue
        parameter (nmenuep=1000)
        integer i,imesh,iflag,irun,icolor(nmenuep)
        character(80) scratch_filename, host
        character(256) comment(nmenuep)
        character ctab
        integer ictab
        equivalence(ctab,ictab)

        data scratch_filename/' '/
        data ictab/9/

        open(unit=99,file='brill.run',status='old')
        read(99,*)irun
        irun=irun+1
        rewind(99)
        write(99,*)irun
        call hostnm(host)
        write(99,*)host(1:len_trim(host))
        call date_and_time(dtday,dttime,dtzone,idatetime)
        write(99,*)'     ',dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(1:4)
        close(99)

        open(unit=99,file='brill_files.lis',status='unknown')
        write(99,'(a)')'brill.run'
        write(99,'(a)')'brill.in'
        write(99,'(a)')'brfl.dat'

	open (unit=10, file='brfl.dat',form='FORMATTED', status='NEW')
	open (unit= 1, file='pb.com',form='FORMATTED', status='NEW')
	open (unit= 2, file='pf.com',form='FORMATTED', status='NEW')

	pi  = 4.*atan(1.)
	pi1 = pi
	ndatfil=0
c	default fuer papier plot
c	Papier ist im Hochformat DINA4
	xpap = 20.5
	ypap = 29.5
c	Abstand des Ursprungs vom unteren Papierrand
	yorig = 3.5
c	Abstand des Ursprungs vom linken Papierrand
	xorig = 4.
c	Abstand des linken  Bildrandes vom linken Papierrand
	xr1   = 1.
c	Abstand des rechten Bildrandes vom linken Papierrand
	xr2   = 19.
c	Abstand des unteren Bildrandes vom unteren Papierrand
	yr1   = 2.
c	Abstand des oberen Bildrandes vom unteren Papierrand
	yr2   = 25.5
	nlab  = 0
	do 89 i=1,50
89	labsize(i) = 40

c	write date and time onto BRFL.DAT
cmsh	call date (datstring)
cmsh	call time (timstring)
cmsh	write (10,90) datstring, timstring
cmsh 90	format (//,5x,' RUN OF PROGRAM BRILLIANCE' /
cmsh	1          5x,' -------------------------',//
cmsh	2          7x,' date: ',9a1,'   time:',8a1//)

        write(10,*)' '
        write(10,*)'     RUN ',irun,' OF PROGRAM BRILLIANCE on ',
     &    host(1:len_trim(host))
        write(10,*)'     -------------------------'
        write(10,*)'     ',dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(1:4)
        write(10,*)' '

c	set BESSY - II data as default
c	June 6, 1986
	gamma = 3327.
	epsx  =  6.1e-9
        epsy  =  6.1e-11
	curr  = 0.1
	rho   = 4.36
	alength = 4.71
	betax   = 16.
	betay   = 3.1
	depthoffield = 0.
	coherent = 0.
	sourceprop = 0.
	ncurve  = 1
        icurve  = 1
        imenue=0

c	main menue
        write(6,*)' '
        write(6,*)'--- The formulas for the wiggler and bending magnet should'
        write(6,*)'--- probably be updated according to G. Geloni.'
1       write(6,*)' '
        print*,' main menue'
        print*,' 1) storage ring parameters'
        print*,' 2) insertion parameters'
	print*,' 3) undulator (according to Walker)'
        print*,' 4) wiggler (old mode, see also item 13)'
	print*,' 5) bending magnet (according to Kim)'
	print*,' 6) exit'
	print*,' 7) label'
	print*,' 8) layout'
	print*,' 9) depth of field '
	print*,' 10) coherent flux '
        print*,' 11) files with source properties '
	print*,' 12) elliptical undulator (according to Walker)'
        print*,' 13) wiggler (N/2 poles are considered, recommended mode)'
        print*,' 14) undulator (acc. to Kim)'
        print*,' 15) elliptical undulator (according to Kim)'
c        print*,'    next number of curves (<1000) : ',i3
        print*,' Enter your choice:   ',icurve

778     read(5,'(a)')cline
        if(cline(1:1).eq.'*'.or.cline(1:1).eq.'#'.or.cline.eq.''
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          print*,cline(1:len_trim(cline))
          read (cline,*)  iflag
        endif
c        write(6,*)' '
c        write(6,*)cline(1:len_trim(cline))
c        write(6,*)' '
        if (iflag.eq.3.or.iflag.eq.4.or.iflag.eq.5.or.
     &      iflag.eq.12.or.iflag.eq.13.or.iflag.eq.14.or.iflag.eq.15) then
          imenue=imenue+1
          if (imenue.gt.ncurve) then
            write(6,*)'*** Maximum number of menue items reached, exiting!'
            iflag=6
          else
            read (cline,*,err=97,end=97)  icolor(imenue),icolor(imenue)
            goto 98
97          continue
            icolor=1
98          continue
            write(comment(imenue),*)'* run ',irun
            do ic=1,len_trim(cline)
              if (cline(ic:ic).eq.'!') then
                comment(imenue)=cline(ic+1:len_trim(cline))
                goto 1234
              endif
            enddo
          endif
        endif
1234    if (iflag.le.0 .or. iflag.gt.15) goto 1

	write (6,'(" iflag = ",i4)') iflag

        goto (11,12,13,14,15,16,17,18,19,20,21,22,113,114,115), iflag
11	call ring
	goto 1
12	call insert
	goto 1
13	ndatfil = ndatfil+1
	call undula_walker
	goto 1
14	ndatfil = ndatfil+1
	call wiggler
	goto 1
15	ndatfil = ndatfil+1
	call bend
	goto 1
17	continue
	call label (1)
	goto 1
18      continue
	call layout
	goto 1
19	continue
	read (5,*) depthoffield
	write (6 ,191) depthoffield
	write (10,191) depthoffield
191	format (
	1 //,'##############switch on depth of field#######',/
	2    ' depth of field = ',f2.0,/)
	goto 1
20	continue
	read(5,*) coherent
	write (6,201) coherent						
	write (10,201) coherent
201	format (//' coherent flux is plotted, coherent = ',f5.1,/
	1         ' photon/sec/.1%BW',  //)
	goto 1
21	continue
	read (5,*) sourceprop
	write (6 ,211) ifix(sourceprop)
	write (10,211) ifix(sourceprop)
211	format
	1(/,' generate files  source.',i3,'ff  to store',/,
	2' effective cross section and divergence',/)
	goto 1
22	continue
	ndatfil = ndatfil+1
	call ellipticalundu_walker
	goto 1
113     ndatfil = ndatfil+1
	call wiggler_msh
	goto 1
114     ndatfil = ndatfil+1
	call undula
	goto 1
115   continue
	ndatfil = ndatfil+1
	call ellipticalundu
	goto 1

16	continue
	ncurve = ncurve - 1

c	plot flux as stored in arrays  curvex and curvey
cmsh	call plotb
cmsh	call plotf
cmsh	call label (2)
cmsh	write (1,999)
cmsh	write (2,999)
cmsh 999	format (
cmsh	1 'plend',/
cmsh	2 'exit')

	close (unit= 1)
	close (unit= 2)
	close (unit=10)

c	Schreibe Kurven einzeln auf SCRATCH f"ur PAW
	do icurve = 1, ncurve
cmsh        encode (3,101,number) icurve
        write(number,fmt=101) icurve
101     format (i3)
        scratch_filename(1:5) = 'brfl.'
        if (icurve.lt.10) then
          scratch_filename(6:6) = number(3:3)
        else if (icurve.lt.100) then
          scratch_filename(6:7) = number(2:3)
        else
          scratch_filename(6:8) = number(3:3)
        endif
        open (unit=51,file=scratch_filename,status='unknown')
        write(cline,*)'* ',dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(1:4)
        write(51,'(a)')cline(2:len_trim(cline))
        write(cline,*)'* ',irun,' ',host(1:len_trim(host)),' ! run and host'
        write(51,'(a)')cline(2:len_trim(cline))
        write(cline,*)'* ',icolor(icurve),' ! color'
        write(51,'(a)')cline(2:len_trim(cline))
        last=len_trim(comment(icurve))
        do i=last,1,-1
          if (comment(icurve)(i:i).ne.ctab) goto 4488
        enddo
4488    write(cline,'(a)')'* '//comment(icurve)(1:i)
        write(51,'(a)')cline(1:len_trim(cline))
        write(51,'(a)')
     &    '* Eph/eV                  Brill           Flux(no ESpread)      Bright'
	do imesh=1,mesh(icurve)-1
          if (brill(imesh,icurve) .le. 0.) brill(imesh,icurve)=1.
          if (flux (imesh,icurve) .le. 0.) flux (imesh,icurve)=1.
          write (51,*)
     &      hanue(imesh,icurve), brill(imesh,icurve), flux(imesh,icurve)
     &      , fluxdens(imesh,icurve)
        enddo
	close (unit=51)
        write(99,'(a)')scratch_filename(1:len_trim(scratch_filename))
	enddo	! icurve

        write(6,*)' --- All done! ---'
cmsh	write (6,161)
cmsh 161	format (' All done!',/' Now start PAW for graphics'
cmsh	3 ' to plot using appropriate .kumac files.'
cmsh	2' and/or "print BRFL.dat" for results on printer',/)

        close(99)
	stop
	end
*CMZ :  2.00/00 06/12/2017  17.06.58  by  Michael Scheer
*CMZ :  1.02/01 21/11/2017  14.02.35  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.13.29  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  09.21.00  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.42.03  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.30.20  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	subroutine bend
	
c	************************************************************
c			   version from Jan 1994
c				B E N D N E W
c	************************************************************

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.
        external function BSKR3
        external function akanue
        real hanuec,hanuemin,hanuemax,hanue1,y1,df2,psimax,dpsi,d1f,psi,x,
     &    arg,dd1F,df,sigrp,alamda,sigr,sigx,sigy,sigxp,sigyp,anorm,b1,eps,
     &    h2,bskr3,akanue				
        integer iend,imesh,ierror,npsi,ipsi
        character a

1	continue	
	write (6,101)
101     format ('$ enter betax, betay (m/rad): ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) betax,betay
        endif
cmsh	read(5,*) betax, betay
2	continue
c	critical energy according to Winick
c	note: rho in meter, hanuec in eV, wavelength in m
	hanuec = 3.*gamma**3/(806500.*4.*pi*rho)
	write (6,91) rho, betax, betay, hanuec
91	format (/,' bending magnet',/
	1 ' radius (meter)          ',f10.4,/
	2 ' betax (m/rad)           ',f10.4,/
	3 ' betay (m/rad)           ',f10.4,/
	4 ' hanuec (eV)             ',f10.4,//
	9 '$correct (y/n)? ')
779       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          a=cline(1:1)
        endif
cmsh	read(5,92) a
92	format (a1)
	if (a.ne.'y') goto 1
	write (10,110) betax, betay, rho, hanuec
110	format (//' spectrum of bending magnet',/
	1         ' --------------------------',/
	1 5x,' betax (m/rad)                    ',f10.4,/
	2 5x,' betay (m/rad)                    ',f10.4,//
	3 5x,' rho (meter)                      ',f10.4,/
	4 5x,' critical photon energy (eV)      ',f10.4,//
	9 5x,' h*nue (eV)      Brilliance           Flux',/
	1 5x '                   photon             photon',/
	2 5x '            sec*(mm*mrad)**2*0.1% B  sec*mrad*0.1% BW',/)
	write (6,102)
102	format ('$ enter photon energy range in eV, hanuemin, hanuemax: ')
710     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) hanuemin, hanuemax
        endif
cmsh	read(5,*) hanuemin, hanuemax
	write (6,103)
103	format ('$ enter number of points in this range, mesh: ')
712     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 712
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256
	iend = mesh(icurve)

c	find number of points
	do 199 imesh = 1, mesh(icurve)
	hanue1 = hanuemin
	1	+ (hanuemax-hanuemin)*(float(imesh-1)/float(mesh(icurve)))**2
	hanue(imesh,icurve) = hanue1
	y1 = hanue1/hanuec
	if (y1.gt.5. . and. iend.ge.mesh(icurve) ) iend = imesh
199 	continue

	mesh(icurve) = iend

	if (sourceprop.gt.0.) then
	call fileopen
	write (11,*) mesh(icurve)-1
	write (12,*) mesh(icurve)-1
	write (13,*) mesh(icurve)-1
	write (14,*) mesh(icurve)-1
	endif

c	main loop
	do 200 imesh = 1, mesh(icurve)
	hanue1 = hanuemin
	1	+ (hanuemax-hanuemin)*(float(imesh-1)/float(mesh(icurve)))**2
	hanue(imesh,icurve) = hanue1
	y1 = hanue1/hanuec

c	use Kwang-Je's Formula
c	units are photons/(sec mrad**2 0.1%BW)
	H2  = (y1*akanue(2./3., y1/2.,ierror))**2
c	this is sometimes called brightness
	dF2 = 3.461e6 * curr * gamma**2 * H2

c	vertically integrated intensity
c	units are photons/(sec mrad 0.1%BW)
c	G1 = y1 * G0(y1)
c	dF = 1.256e10 * curr * gamma * G1

c	numerical integration in vertical
c	for details see WIGGNEW

        psimax = 20.*2./(gamma*sqrt(y1))
        npsi = 400
        dpsi = psimax/npsi
        d1F = 0.
        do ipsi = 1,npsi
          psi = (psimax/npsi)*(ipsi-.5)
          x = gamma * psi
          arg = y1 * (1.+X**2)**(3./2.) / 2.

          if (arg.gt.40.) exit
          bskr31=bskr3(arg,1)
          bskr32=bskr3(arg,2)
          if (bskr31.lt.1.0e-20) bskr31=0.0
          if (bskr32.lt.1.0e-20) bskr32=0.0
          dd1F = 3.461 e6 * gamma**2 * curr * (y1 * (1.+X**2))**2*
     &      (BSKR32**2 + X**2/(1.+X**2) * BSKR31**2)

c        print*,"BENDNEW, BSKR3:",arg,bskr3(arg,2),bskr3(arg,1)
c	        dd1F = 3.461e6*gamma**2*curr * (y1*(1.+X**2))**2*
c	1 	(BSKR3(arg,2)**2 + x**2/(1.+x**2)*BSKR3(arg,1)**2)

          d1F = d1F + dd1F*dpsi*1000.
        enddo

        d1F = d1F*2.
        dF = d1F

c	Divergence of radiation in rad
	sigrp = (1./(1000.*sqrt(2.*pi))) * dF/dF2

c	wavelength in meter
	alamda  = 1./(hanue1*806500.)

c	diffraction limited source size in m
	sigr  = alamda/(4.*pi*sigrp)

c	e beam cross section and divergence
c	this assumes that the tangent point is a beam waist
c	(highly unlikely) and ignores contributions to beam size due
c	to energy spread
	sigx     = sqrt(epsx*betax)
	sigy     = sqrt(epsy*betay)
	sigxp    = sqrt(epsx/betax)
	sigyp    = sqrt(epsy/betay)

	if (sourceprop.gt.0.) then
	write (11,666)
	1 hanue1, sqrt(sigx**2+sigr**2), sigx, sigr
666	format(
	1 2x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6 )
	write (12,667)
	1 hanue1, sqrt(sigxp**2+sigrp**2), sigxp, sigrp
667	format (
	1 2x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6 )

	write (13,668)
	1 hanue1, sqrt(sigy**2+sigr**2), sigy, sigr
668	format(
	1 2x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6 )
	write (14,669)
	1 hanue1, sqrt(sigyp**2+sigrp**2), sigyp, sigrp
669	format (
	1 2x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6 )

	endif

c	Normalisation according to
c	K.J.Kim in NIM A261, 44 (1987) (Novosibirsk)
c	units are  mm**2
c	anorm = 2.*pi*sqrt(
c	1 (sigx**2+sigr**2)*(sigy**2*(1.+(sigyp/sigrp)**2)+sigr**2))
c	make expression more symmetric in x and y. AG 94/01/05
	anorm = 2.*pi*sqrt(
	1 (sigx**2+sigr**2)*(sigy**2+sigr**2)*(1.+(sigyp/sigrp)**2))

c	units are photons/(sec (mm*mrad)**2 0.1%BW)
	b1 = dF2/(anorm*1.e6)
	brill(imesh,icurve) = b1
	
c	vertically integrated flux for 1 mrad horizontal angle
c	at curvex in curvey
	flux (imesh, icurve) = dF

c	#######################
c	coherent flux
c	#######################
	if (coherent.gt.0.) then
c	formula due to K.J.Kim, LBL 2236
c	see comment in undula.for   Nov. 30, 1992   AG

c	photon energy in keV
	eps = hanue(imesh,icurve)/1000.
	flux(imesh,icurve) = .385*b1/(hanue(imesh,icurve)**2)
	endif

c	write (6,201) hanue(imesh,icurve), brill(imesh,icurve),
c	1             flux(imesh,icurve)
c201	format (' hanue=',e10.4,'  b = ',e10.4,' flux=',e10.4)

c	hanl1 = alog10(hanue1)
c	hanl2 = alog10(hanue2)
c	blog1 = alog10(b1)
c	blog2 = alog10(b2)
c	if (hanl1.lt.hnminl .or. hanl2.gt.hnmaxl) goto 200
c	if (blog1.gt.bmaxl   .or. blog2.gt.bmaxl  ) goto 200
c	if (blog1.lt.bminl   .or. blog2.lt.bminl  ) goto 200
c	y1 = yscale*(blog1-bminl)
c	x1 = xscale*(hanl1-hnminl)
c	y2 = yscale*(blog2-bminl)
c	x2 = xscale*(hanl2-hnminl)
c	if (y1.lt.0. .or. y2.lt.0.) goto 200
c	call vector (x1,y1, x2,y2)


c	write (10,203) hanue(imesh,icurve), brill(imesh,icurve),
c	1              flux(imesh,icurve)
	write (10,203) hanue(imesh,icurve), brill(imesh,icurve),
	1              flux(imesh,icurve)
203	format (2x,e12.6,10x,e12.6,10x,e12.6)
200	continue

	if (sourceprop.gt.0.) then

	write (11,211) iend-1
211	format (' hanue(eV), SIGMAX, sigmax, sigmar',/
	1' number of points ',i5)
	write (12,212) iend-1
212	format (' hanue(eV), SIGMAX", sigmax", sigmar',/
	1' number of points ',i5)
	write (13,213) iend-1
213	format (' hanue(eV), SIGMAY,  sigmay, sigmar',/
	1' number of points ',i5)
	write (14,214) iend-1
214	format (' hanue(eV), SIGMAY", sigmay", sigmar',/
	1' number of points ',i5)

	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif

 	icurve = icurve+1
	ncurve = ncurve+1
	return
	end
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.13.04  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.30.45  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	****************************************************************
c				P L O T B
c	****************************************************************
	subroutine plotb
c	plots central brillance

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	iunit = 1
c	set all default
c	Hochformat auf DINA 4
c	---------------------
	xpap = 20.5
	ypap = 28.5

c	Extremwerte

	hmn  = hanue(1,1)
	hmx  = hmn
	bmn  = brill(1,1)
	bmx  = bmn

	do 80 icurve = 1, ncurve
	do 80 imesh  = 1, mesh(icurve)
	if (hmn.gt.hanue(imesh,icurve)) hmn = hanue(imesh,icurve)
	if (hmx.lt.hanue(imesh,icurve)) hmx = hanue(imesh,icurve)
	if (bmn.gt.brill(imesh,icurve)) bmn = brill(imesh,icurve)
	if (bmx.lt.brill(imesh,icurve)) bmx = brill(imesh,icurve)
80	continue

	hnmin = hmn* .98
	hnmax = hmx*1.02
	bmin  = bmn* .98
	bmax  = bmx*1.02
	write (6,101) hmn, hmx, bmn, bmx
101	format (' Plot of Brilliance',/
	1' Result of min/max search',/
	1' hanue min (eV)   : ',e10.4,/
	2' hanue max (eV)   : ',e10.4,/
	3' brill min        : ',e10.4,/
	4' brill max        : ',e10.4,/
	5' enter new h nue min and h nue max in eV: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) s1, s2
        endif
	if (s1.gt.0.) hnmin = s1
	if (s2.gt.0.) hnmax = s2

	write (6,102)
102	format (
	1' enter new brill min and brill max  ph/[sec(mm*mrad)**2*0.1%BW] :')
779      read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) s1, s2
        endif
	if (s1.gt.0.) bmin = s1
	if (s2.gt.0.) bmax = s2

	hnminl = alog10(hnmin)
	hnmaxl = alog10(hnmax)
	bminl  = alog10(bmin)
	bmaxl  = alog10(bmax)

c	Masstab cm/Zehnerpotenz
	xscale = (xr2-xorig)/(hnmaxl-hnminl)
	ybscale= (yr2-yorig)/(bmaxl -bminl )

c	retten der Parameter auf tape 10 und auf Schirm
	write (10,2) hnmin,hnmax,bmin,bmax,hnminl,hnmaxl,bminl,bmaxl,
	1            xscale,ybscale,xr1,xr2,yr1,yr2,xorig,yorig,
	2	     (xlabb(i), ylabb(i),(lab(j,i),j=1,8),i=1,nlab)
	write ( 6,2) hnmin,hnmax,bmin,bmax,hnminl,hnmaxl,bminl,bmaxl,
	1            xscale,ybscale,xr1,xr2,yr1,yr2,xorig,yorig,
	2	     (xlabb(i), ylabb(i),(lab(j,i),j=1,8),i=1,nlab)
2	format (/,'  Parameter fuer Plot der axialen Brillanz',/
	1 ' hnmin, hnmax (eV)            ',f10.4,2x,f13.4,/
	3 ' bmin , bmax                  ',e10.4,2x,e10.4,/
	5 ' hnminl hnmaxl                ',f10.4,2x,f10.4,/
	6 ' bminl, bmaxl                 ',f10.4,2x,f10.4,/
	7 ' xscale, ybscale (cm/potenz)  ',f10.4,2x,f10.4/
	9 ' xr1, xr2 (cm)                ',f10.4,2x,f10.4,/
	1 ' yr1, yr2 (cm)                ',f10.4,2x,f10.4,/
	2 ' xorig, yorig (cm)            ',f10.4,2x,f10.4,//
	3 ' labels:',/
	4  (' (',f5.1,',',f5.1,') : ',8a1,'   (',f5.1,',',f5.1,') : ',8a1))

c	schreibe NEMO Kommandos auf File PB.COM
	write (iunit,1)
1	format (
	1 '$set term/nowrap',/
	2 '$set noverify',/
c	2 '$run prg:nemo',/
	2 '$run prg:nemo94',/
c	2 '$run prg:nemo94.exe_alpha',/
	3 'input',/
	4 'clear',/
	5 'plinit        51'/
	6 'setsym   1.     ')
c	Achsenkreuz
	x1 = 0.
	x2 = xscale*(hnmaxl-hnminl)
	y1 = 0.
	y2 = ybscale*(bmaxl - bminl )
	call vector (x1,y1, x1,y2)
	call vector (x1,y2, x2,y2)
	call vector (x2,y2, x2,y1)
	call vector (x2,y1, x1,y1)

	call xachse
	call yachse

	do 500 icurve = 1,ncurve
	do 501 imesh  = 1,mesh(icurve)-1
	x1 = hanue(imesh  ,icurve)
	x2 = hanue(imesh+1,icurve)
	y1 = brill(imesh  ,icurve)
	y2 = brill(imesh+1,icurve)
	if (x1.le.0. .or. x2.le.0. .or. y1.le.0. .or. y2.le.0.) goto 501
	if (x1.lt.hnmin .or. x2.lt.hnmin .or. y1.lt.bmin .or. y2.lt.bmin
	1  .or.
	2   x1.gt.hnmax .or. x2.gt.hnmax .or. y1.gt.bmax .or. y2.gt.bmax)
	3      goto 501
	x1 = xscale*(alog10(x1) - hnminl)
	x2 = xscale*(alog10(x2) - hnminl)
	y1 =ybscale*(alog10(y1) - bminl )
	y2 =ybscale*(alog10(y2) - bminl )
	call vector (x1, y1, x2, y2)
501	continue
500	continue

	return
	end

c	****************************************************************
c				P L O T F
c	****************************************************************
	subroutine plotf
c	plot of total flux

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	iunit = 2
c	set all default
c	Hochformat auf DINA 4
c	---------------------
	xpap = 20.5
	ypap = 28.5

c	Extremwerte
	hmn  = hanue(1,1)
	hmx  = hmn
	fmn  = flux(1,1)
	fmx  = fmn

	do 80 icurve = 1, ncurve
	do 80 imesh  = 1, mesh(icurve)
	if (hmn.gt.hanue(imesh,icurve)) hmn = hanue(imesh,icurve)
	if (hmx.lt.hanue(imesh,icurve)) hmx = hanue(imesh,icurve)
	if (fmn.gt.flux (imesh,icurve)) fmn = flux (imesh,icurve)
	if (fmx.lt.flux (imesh,icurve)) fmx = flux (imesh,icurve)
80	continue

	hnmin = hmn* .98
	hnmax = hmx*1.02
	fmin  = fmn* .98
	fmax  = fmx*1.02
	write (6,101) hmn, hmx, fmn, fmx
101	format (' Result of min/max search',/
	1' hanue min (eV)   : ',e10.4,/
	2' hanue max (eV)   : ',e10.4,/
	3' flux  min        : ',e10.4,/
	4' flux  max        : ',e10.4,/
	5' enter new h nue min and h nue max in eV: ')
710      read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) s1, s2
        endif
	if (s1.gt.0.) hnmin = s1
	if (s2.gt.0.) hnmax = s2

	write (6,102)
102	format (
	1' enter new flux min and flux max in ph/[sec*0.1%BW] :')
711      read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          read(cline,*) s1, s2
        endif
	if (s1.gt.0.) fmin = s1
	if (s2.gt.0.) fmax = s2

	hnminl = alog10(hnmin)
	hnmaxl = alog10(hnmax)
	fminl  = alog10(fmin)
	fmaxl  = alog10(fmax)

c	Masstab cm/Zehnerpotenz
	xscale = (xr2-xorig)/(hnmaxl-hnminl)
	yfscale= (yr2-yorig)/(fmaxl-fminl)

c	retten der Parameter auf tape 10 und auf Schirm
	write (10,2) hnmin,hnmax,fmin,fmax,hnminl,hnmaxl,fminl,fmaxl,
	1            xscale,yfscale,xr1,xr2,yr1,yr2,xorig,yorig,
	2	     (xlabf(i), ylabf(i),(lab(j,i),j=1,8),i=1,nlab)
	write ( 6,2) hnmin,hnmax,fmin,fmax,hnminl,hnmaxl,fminl,fmaxl,
	1            xscale,yfscale,xr1,xr2,yr1,yr2,xorig,yorig,
	2	     (xlabf(i), ylabf(i),(lab(j,i),j=1,8),i=1,nlab)
2	format (/,'  Parameter fuer Plot des Flusses',/
	1 ' hnmin, hnmax (eV)            ',f10.2,2x,f13.2,/
	3 ' fmin , fmax                  ',e10.3,2x,e10.3,/
	5 ' hnminl hnmaxl                ',f10.4,2x,f10.4,/
	6 ' fminl, fmaxl                 ',f10.4,2x,f10.4,/
	7 ' xscale, yfscale (cm/potenz)  ',f10.4,2x,f10.4/
	9 ' Rand xr1, xr2 (cm)           ',f10.4,2x,f10.4,/
	1 '      yr1, yr2 (cm)           ',f10.4,2x,f10.4,/
	2 ' xorig, yorig (cm)            ',f10.4,2x,f10.4,//
	3 ' labels:',/
	4  (' (',f6.1,',',f6.1,') : ',8a1,'   (',f6.1,',',f6.1,') : ',8a1))

c	schreibe NEMO Kommandos auf File PF.COM
	write (iunit,1)
1	format (
	1 '$set term/nowrap',/
	2 '$set noverify',/
	2 '$run prg:nemo94',/
	3 'input',/
	4 'clear',/
	5 'plinit        51'/
	6 'setsym   1.     ')

c	Achsenkreuz
	x1 = 0.
	x2 = xscale*(hnmaxl-hnminl)
	y1 = 0.
	y2 = yfscale*(fmaxl-fminl)
	call vector (x1,y1, x1,y2)
	call vector (x1,y2, x2,y2)
	call vector (x2,y2, x2,y1)
	call vector (x2,y1, x1,y1)

	call xachse
	call yachse

	do 500 icurve = 1,ncurve
	do 501 imesh  = 1,mesh(icurve)-1
	x1 = hanue(imesh  ,icurve)
	x2 = hanue(imesh+1,icurve)
	y1 = flux (imesh  ,icurve)
	y2 = flux (imesh+1,icurve)
	if (x1.le.0. .or. x2.le.0. .or. y1.le.0. .or. y2.le.0.) goto 501
	if (x1.lt.hnmin .or. x2.lt.hnmin .or. y1.lt.fmin .or. y2.lt.fmin
	1  .or.
	2   x1.gt.hnmax .or. x2.gt.hnmax .or. y1.gt.fmax .or. y2.gt.fmax)
	3      goto 501
	x1 = xscale*(alog10(x1) - hnminl)
	x2 = xscale*(alog10(x2) - hnminl)
	y1 =yfscale*(alog10(y1) - fminl )
	y2 =yfscale*(alog10(y2) - fminl )
	call vector (x1, y1, x2, y2)
501	continue
500	continue

	return
	end

c	**************************************************************
c				X A C H S E
c	**************************************************************
	subroutine xachse
c	beschrifte x-Achse
c	------------------

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	xpap = 20.5
	ypap = 28.5
	n = int(hnminl)
c	if (abs(n-hnminl).gt. .01) n = n+1
11	continue
c	mache Strich bei 10**n
	x1 = xscale*(n-hnminl)
	y1 = 0.
	x2 = x1
	y2 = y1+.2

	if (x1.lt.0. .or. x1.gt.xscale*(hnmaxl-hnminl)) goto 114
	call vector (x1,y1,x2,y2)
c	number of units
	nx1=int(100.*(y1-1.0+yorig))
	ny1=int(100.*(xpap-(x1+xorig-0.35)))
	nheight=40
	richt=-90.
	write (iunit,112) nx1,ny1,nheight,richt
112	format ('label   ',3i4,f4.0,'10')
c	nx1=int(100.*(y1+yorig-0.3))
c	nx1=int(100.*(y1+yorig-0.7))
	nx1=int(100.*(y1+yorig-0.5))
c	ny1=int(100.*(xpap-(x1+xorig+0.1)))
	ny1=int(100.*(xpap-(x1+xorig+0.5)))
	nheight=25
	nexp=n
	write (iunit,113) nx1,ny1,nheight,richt,nexp
113	format ('label   ',3i4,f4.0,i2)

114	y2 = y1 + 0.1
	do 111 i=2,9
	x1 = xscale * (n+alog10(float(i))-hnminl)
	if (x1.le.0. .or. x1.ge.(xscale*(hnmaxl-hnminl))) goto 111
	call vector (x1,y1,x1,y2)
111	continue
	n = n+1
	if (float(n).le.(hnmaxl+1.)) goto 11

c	Schreibe Text und Einheiten
c	---------------------------
	nx1=int(100.*(yorig-2.0))
	ny1=int(100.*(xpap-(xorig+4.8)))
	nheight=40
	richt=-90.
	write (iunit,115) nx1,ny1,nheight,richt
115	format('label   ',3i4,f4.0,'photon energy [eV]')

	return
	end
c	************************************************************
c			   Y A C H S E
c	************************************************************
	subroutine yachse
c	beschriften der y-Achse
c	-----------------------

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	xpap = 20.5
	ypap = 29.5
	if (iunit.ge.2) goto 500

	n = int(bminl)
c	mache Strich bei 10**n
12	continue
	x1 = 0.
       	y1 = ybscale*(n-bminl)
	x2 = x1+.2
	y2 = y1
	if (y1.lt.0. .or. y1.gt. ybscale*(bmaxl-bminl)) goto 114
	call vector (x1,y1,x2,y2)
c	number of units
	nx1=int(100.*( y1+yorig-0.2))
	ny1=int(100.*(xpap-(x1+xorig-1.3)))
	nheight=40
	richt=-90.
	write (iunit,112) nx1,ny1,nheight,richt
112	format ('label   ',3i4,f4.0,'10')
	nx1=int(100.*(y1+yorig+0.2))
	ny1=int(100.*(xpap-(x1+xorig-.5)))
	nheight=25
	nexp=n
	write (iunit,113) nx1,ny1,nheight,richt,nexp
113	format ('label   ',3i4,f4.0,i2)

114	x2 = x1 + 0.1
	do 121 i=2,9
	y1 = ybscale*(n-bminl+alog10(float(i)))
	if (y1.ge.0. .and. y1.le.ybscale*(bmaxl-bminl))
	1	call vector(x1,y1,x2,y1)
121	continue
	n = n+1
	if (float(n).le.bmaxl+1.) goto 12

	nx1=int(100.*(yorig+6.7))
	ny1=int(100.*(xpap-(-1.8+xorig)))
	nheight=40
	richt=0.
	write (iunit,214) nx1,ny1,nheight,richt
214	format('label   ',3i4,f4.0,'brilliance [')
	nx1 = nx1 + 520
	ny1 = ny1 + 40
	nheight = 25
	write (iunit,215) nx1,ny1,nheight,richt
215	format('label   ',3i4,f4.0,'photons')
	nx1 = nx1-200
	ny1 = ny1-58
	write (iunit,216) nx1,ny1,nheight,richt
216	format ('label   ',3i4,f4.0,'sec (mm*mrad)**2 0.1% BW')
	nx1 = nx1 + 490
	ny1 = ny1 + 58 - 40
	nheight = 40
	write (iunit,217) nx1,ny1,nheight,richt
217	format ('label   ',3i4,f4.0,' ]')
	nx1 = nx1 - 500
	ny1 = ny1 +  20
	nx2 = nx1 + 505
	ny2 = ny1
	write (iunit,218) nx1,ny1,nx2,ny2
218	format ('vector  ',4i4)	
	return

500	continue
	n = int(fminl)
c	mache Strich bei 10**n
512	continue
	x1 = 0.
       	y1 = yfscale*(n-fminl)
	x2 = x1 + .2
	y2 = y1
	if (y1.lt.0. .or. y1.gt. yfscale*(fmaxl-fminl)) goto 614
	call vector (x1,y1,x2,y2)
c	number of units
	nx1=int(100.*( y1+yorig-0.2))
	ny1=int(100.*(xpap-(x1+xorig-1.3)))
	nheight=40
	richt=-90.
	write (iunit,612) nx1,ny1,nheight,richt
612	format ('label   ',3i4,f4.0,'10')
	nx1=int(100.*(y1+yorig+0.2))
	ny1=int(100.*(xpap-(x1+xorig-.6)))
	nheight=25
	nexp=n
	write (iunit,613) nx1,ny1,nheight,richt,nexp
613	format ('label   ',3i4,f4.0,i2)

614	x1 = 0.
	x2 = x1 + 0.1
	do 621 i=2,9
	y1 = yfscale*(n-fminl+alog10(float(i)))
	if (y1.ge.0. .and. y1.le.yfscale*(fmaxl-fminl))
	1	call vector(x1,y1,x2,y1)
621	continue
	n = n+1
	if (float(n).le.fmaxl+1.) goto 512

c	Beschriftung fuer Fluss-Plot
	nx1=int(100.*(yorig+4.2))
	ny1=int(100.*(xpap-(-1.8+xorig)))
	nheight=40
	richt=0.
c	########################
c	coherent power
c	########################
	if (coherent.gt.0.)then
	write (iunit,5141) nx1,ny1,nheight,richt
5141	format('label   ',3i4,f4.0,'coherent flux [')
	nx1 = nx1 + 500
	ny1 = ny1 + 40
	nheight = 25
	write (iunit,5151) nx1,ny1,nheight,richt
5151	format('label   ',3i4,f4.0,'photons')
	nx1 = nx1-40
	ny1 = ny1-58
	write (iunit,5161) nx1,ny1,nheight,richt
5161	format ('label   ',3i4,f4.0,'sec 0.1%BW')
	nx1 = nx1 + 180
	ny1 = ny1 - 40 + 58
	nheight = 40
	write (iunit,5171) nx1,ny1,nheight,richt
5171	format ('label   ',3i4,f4.0,' ]')
	nx1 = nx1 - 180
	ny1 = ny1 + 20
	nx2 = nx1 + 200
	ny2 = ny1
	write (iunit,5181) nx1,ny1,nx2,ny2
5181	format ('vector  ',4i4)	


	else
	write (iunit,514) nx1,ny1,nheight,richt
514	format('label   ',3i4,f4.0,'flux [')
	nx1 = nx1 + 310
	ny1 = ny1 + 40
	
	nheight = 25
	write (iunit,515) nx1,ny1,nheight,richt
515	format('label   ',3i4,f4.0,'photons')
	nx1 = nx1-100
	ny1 = ny1-58
	write (iunit,516) nx1,ny1,nheight,richt
516	format ('label   ',3i4,f4.0,'sec * 0.1% BW')
	nx1 = nx1 + 280
	ny1 = ny1 - 40 + 58
	nheight = 40
	write (iunit,517) nx1,ny1,nheight,richt
517	format ('label   ',3i4,f4.0,' ]')
	nx1 = nx1 - 290
	ny1 = ny1 + 20
	nx2 = nx1 + 290
	ny2 = ny1
	write (iunit,518) nx1,ny1,nx2,ny2
518	format ('vector  ',4i4)	
	nx1 = int(100.*(yorig+12.0))
	ny1 = int(100.*(xpap-(-1.8+xorig)))
	nheight = 30
	write (iunit,519) nx1,ny1,nheight,richt
519	format ('label   ',3i4,f4.0,' (wiggler & dipole: 1 mrad)')
c	write (iunit,519) nx1,ny1,nheight,richt
c519	format ('label   ',3i4,f4.0,' (horiz. acceptance: 20 mrad)')
	endif

	return
	end
c	******************************************************************
c				V E C T O R
c	******************************************************************
	subroutine vector(x1,y1,x2,y2)

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.
	
	xpap= 20.5
c	does a line, coordinates are in cm with respect to the
c	origin
C	NEMO dreht das Bild auf dem Papier
	scalex = 1.
	scaley = 1.
	nx1 = int(scaley*100.*(y1+yorig))
	ny1 = int(scalex*100.*(xpap-(x1+xorig)))
	nx2 = int(scaley*100.*(y2+yorig))
	ny2 = int(scalex*100.*(xpap-(x2+xorig)))
	write (iunit,1) nx1,ny1,nx2,ny2
1	format('vector  ',4i4)
	return
	end
*CMZ :  1.02/01 21/11/2017  14.03.10  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.13.51  by  Michael Scheer
*CMZ :  1.01/00 12/09/2013  15.01.41  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.12.29  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.32.20  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	****************************************************************
c			  ELLIPTICAL UNDULATOR
c	****************************************************************
	subroutine ellipticalundu
c 	expanded for elliptical undulator 04/12 AG
c	normalization changed according to K.J.Kim  April 11, 1988

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

	write (6,998)
998	format (' Elliptical undulator under development Dec 2004')

1	write (6,20)
20	format (' enter periode length in m: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alamd0
        endif
	write (6,10)
10	format (' enter number of undulator periodes: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) an
        endif
	if ( an*alamd0.ge.alength) write (6,21) alamd0, an, alength
21	format (' **** warning ****',/
	1 ' lamda0  = ',f10.2,' m',/
	2 '      N  = ',f10.0,/
	3 ' alength = ',f10.2,' m '/
	4 ' *** not compatible ***'/)
	
	write (6,30)
30	format(' enter range in effective wiggler strength
     1	 KeffMIN and KeffMAX: ')
710      read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) akmin, akmax
        endif
	if (akmin.ge.0. .and. akmax.ge.akymin) goto 32
	write (6,31) akmin, akmax
31	format (' **** warning ****',/
	1 ' AKMIN   = ',f10.3,'   AKMAX  = ',f10.3,' *** inconsistent'/)
32	continue
	write (6,35)
35	format (' enter ratio of horizontal to vertical wiggler '
     1  ' strength akxky = q =  Kx/Ky.'/,
     1  ' P**2 I in APPLE II is optimized for ',/
     2  ' first harmonic    q = 0.9999',/
     3  ' third harmonic    q = 0.42  ',/
     4  ' fifth harmonic    q = 0.32  ',/
     5  ' seventh harmonic  q = 0.27  ',/
     6  ' nineth harmonic   q = 0.24  ',/
     7  ' eleventh harmonic q = 0.22  ',/
     9  /)
711     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          read(cline,*) akxky
        endif
	write (6,'(" Kx/Ky = ",e12.6)') akxky

	if (akxky.gt.1.) then
	write (6,'(//," ?????????  in Ellipticalundu Kx/Ky  must be <=1 ",//)')
	write (10,'(//," ????????? in Ellipticalundu Kx/Ky  must be <=1 ",//)')
	endif

	akymin=akmin/(sqrt(1.+akxky**2))
	akymax=akmax/(sqrt(1.+akxky**2))

	write (6,40)
40	format (' enter order of undulator spectrum (1, 3, 5,..): ')
713     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          read(cline,*) smallk
        endif
	write (6,'(" smallk = ",f10.3)') smallk

	if (smallk.lt.1. .or.
	1   abs(smallk-float(int(smallk))).gt.1.e-4) write (6,41) smallk
41	format (' **** warning ****',/
	1 ' bad value for order of spectrum. SMALLK = ',f10.3)
2	continue

c	get minimum gap for hybrid-REC structure using Halbach's formula
c	for REC without 10% reduction
	gapmin = alamd0*(1.519-sqrt(2.309+.5556*alog(.003223*akymax/alamd0)))	
	write (6,100) an, alamd0, akymin, akymax, akxky, smallk, gapmin
100	format (/,' undulator parameters',/
	1 ' number of periodes             AN      = ',f10.0,/
	2 ' periode length (m)             alamd0  = ',f10.5,/
	3 ' minimum K                      KeffMIN = ',f10.2,/
	4 ' maximum K                      KeffMAX = ',f10.2,/
	4 ' ratio Kx/Ky                    akxky   = ',f10.2,/
	5 ' order of spectrum              smallk  = ',f10.0,/
	6 ' minimum gap (hybrid REC) (m)   gapmin  = ',f10.5,//
	9 '$ correct? (y/n)  ')
714       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 714
        else
          a=cline(1:1)
        endif
cmsh	read(5,101) a
101	format (a1)
	if (a.eq.'y') goto 109
	goto 1
109	continue

c	plot undulator spectrum
	write (6,110)
110	format (' enter photon energy range (eV) of this curve',/
	1 '$ hanuemin and hanuemax:  ')
715       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 715
        else
          read(cline,*) hanuemin, hanuemax
        endif
	if (hanuemax.le.0.) then
		hanuemax = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*
	2		(1.+0.5*aKymin**2*(1.+akxky**2))))
	endif
	if (hanuemin.le.0.) then
		hanuemin = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*
	2		(1.+0.5*aKymax**2*(1.+akxky**2))))
	endif
	write (6,120)
120	format ('Enter number of points on this curve, mesh:  ')
716     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 716
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256

c	calculate normalisation of brilliance for K = 0
c	only for initial printout, not for final curve
	theta0 = 1./(gamma*sqrt(an*smallk))
c	note alength is in m
	sigmax = sqrt(epsx*betax)
  	sigmay = sqrt(epsy*betay)
	sigmaxp= sqrt(theta0**2+epsx/betax)
	sigmayp= sqrt(theta0**2+epsy/betay)
	write (6,121) theta0, sigmax, sigmay, sigmaxp, sigmayp, dgamma
121	format (/
	1 ' theta0 (rad, K=0) = ',e10.4,/
	2 ' sigmax (m)        = ',e10.4,/
	3 ' sigmay (m)        = ',e10.4,/
	4 ' sigmaxp (rad) K=0 = ',e10.4,/
	5 ' sigmayp (rad) K=0 = ',e10.4,/
	6 ' sigma E           = ',e10.4,//)
	
c	get maximum power according to K.J.Kim
	call power (an, aKymax, alamd0, gamma, curr, Ptot, d2P)


	write (10,122) alamd0, an, akymin, akymax, akxky, smallk, gapmin,
	1              theta0,sigmax, sigmay,sigmaxp,sigmayp,dgamma
	write (10,123) Ptot, d2P

122	format (/
	1 ' properties of undulator',/
	1 ' -----------------------',/
	2 5x,' periode length (m)             	',f10.4,/
	3 5x,' number of periodes              	',f10.4,/
	4 5x,' minimum Ky                      	',f10.4,/
	5 5x,' maximum Ky                       ',f10.4,/
	5 5x,' ratio Kx/Ky                      ',f10.4,/
	6 5x,' order of spectrum k              ',f10.4,/
	7 5x,' minimum gap (REC hybrid, m)      ',f10.4,//
	1 ' properties of radiation',/
	2 5x,' theta0 (rad) K = 0               ',e10.4,/
	3 5x,' sigmax (m)                       ',e10.4,/
	4 5x,' sigmay (m)                       ',e10.4,/
	5 5x,' sigmaxp (rad) K = 0              ',e10.4,/
	6 5x,' sigmayp (rad) K = 0              ',e10.4,/
	7 5x,' dgamma                           ',e10.4)
123	format (
	7 5x,' Ptot (Watt)                      ',e10.4,/
	7 5x,' dP/dOmega (Watt/mrad^2)          ',e10.4,/
	9 5x,' h*nue (eV)         Brilliance              Flux         K
	1       Ptot    gap',/
	2 5x '                       photon              photon
	3     (Watt)     (m)',/
	4 5x '          sec*(mm*mrad)**2*0.1% BW      sec*0.1% BW  ',/)

	if (sourceprop.gt.0.) then
	call fileopen
	write (11,*) mesh(icurve)-1
	write (12,*) mesh(icurve)-1
	write (13,*) mesh(icurve)-1
	write (14,*) mesh(icurve)-1

c	e beam cross section and divergence
	sigmax = sqrt(epsx*betax)
	sigmay = sqrt(epsy*betay)
	sigxpr = sqrt(epsx/betax)
	sigypr = sqrt(epsy/betay)

	endif

	do 200 imesh=1,mesh(icurve)
	hanue1 = hanuemin+(hanuemax-hanuemin)*(imesh-1)/float(mesh(icurve)-1)
c	8065 cm**-1  = 1 eV
	ak1 = 2.*(smallk*2.*gamma**2/(hanue1*alamd0*806500.)-1.)
	if (ak1.le.0.) goto 200

	aky1 = sqrt(ak1/(1.+akxky**2))
	akx1 = aky1*akxky
	if (aky1.lt.akymin .or. aky1.gt.akymax) goto 200

c	formula due to Krinsky in Koch's book, eq.281, p.152
c	checked on April 9, 1985
c	units are  photons/[sec*mrad**2*0.1%BW*Amp])
c	our efka is the same as F_n(K) in xray data booklet
c	b1 = 4.56E7*(an*gamma)**2*curr*efka(smallk,ak1)
	b1= 4.56E7*(an*gamma)**2*curr*efkaelli(aky1,akx1)
c	this is the spectral brightness, which is reduced by a
c	finite energy spread
c	here we take the spectral distribution of undulator radiation
c	to be a normal with width sigma = 1./(smallk*2.6*N)
cmsh	to be a normal with width sigma = 1./(smallk*2.*sqrt(2.)*N)
c	The energy spread of the e beam is multiplied by 2 because
c	gamma**2 enters the energy formula
c	This is necessary as discovered by Howard Padmore at ALS in
c	winter 1993/94.                         AG 1994/01/16

cmsh        b1= b1/sqrt(1.+(2.*dgamma*smallk*2.6*aN)**2)
        b1= b1/sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2)
        dgamred=1./sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2) !msh

cmsh{

c	save on axis spectral brightness into array fluxdens
c	need to correct for electron beam divergence
c	e beam divergence
	sigxpr = sqrt(epsx/betax)
	sigypr = sqrt(epsy/betay)
c	natural divergence of radiation
        alamda = (alamd0/(smallk*2.*gamma**2))*(1.+(aky1**2+akx1**2)/2.)
	aL = alamd0*an
	sigmarp = sqrt(alamda/(aL))

        print*,'*** Fluxdens unclear in ellipticalundu! ***'
	fluxdens(imesh,icurve) = b1/
	1 sqrt(1.+(sigxpr**2+sigypr**2)/sigmarp**2)

cmsh}

c	flux(imesh,icurve) is the total flux at energy curvex
c	here we make a very crude approximation,
c	since we calculate the flux by multiplying the angular brightness
c	by 2*pi*theta0**2
c	our formula is the same as eq.17 in xray data booklet
c	the formula for sigmar is consistent with Walker
c	ESRP - IRM - 54/84
c	***** sigmar = sqrt (lambda / L) *****

c	J.B. 30.5.89
	sigmar = 1./gamma*sqrt((1.+(akx1**2+aky1**2)/2.)/(smallk*2.*an))
	flux(imesh,icurve) = 2.*pi*b1*(sigmar*1000.)**2
	hanue (imesh, icurve) = hanue1

	call normuelli (anorm, aky1,akx1)

	if (imesh.gt.256 .or. icurve.gt.100 .or. anorm.le.0.)
     &    write (6,*) ' *** imesh=',imesh,'   icurve=',icurve,
     &    ' *** anorm=',anorm
        b1= flux(imesh,icurve)/anorm
        brill(imesh,icurve)= b1

	if (sourceprop.gt.0.) then
c	natural divergence of radiation
          alamda = (alamd0/(smallk*2.*gamma**2))*(1.+(aky1**2+akx1**2)/2.)
          aL = alamd0*an
          sigmarp = sqrt(alamda/(aL))
c	radiation cross section
          sigmar = sqrt(alamda*al)/(4.*pi)
c	MSH behauptet basierend auf einem Vergleich der zweiten Momente
c	(2000/06):
c	sigmar=(1/2 pi) sqrt(2 lambda L)
c	sigmar * sigmarp = (1/2 pi) lambda (*0.82)

          size_xt = sqrt(sigmax**2+sigmar**2)
          dive_xt = sqrt(sigxpr**2 + sigmarp**2)
          size_yt = sqrt(sigmay**2+sigmar**2)
          dive_yt = sqrt(sigypr**2 + sigmarp**2)

          write (11,*) hanue1, size_xt, sigmax, sigmar
          write (12,*) hanue1, dive_xt, sigxpr, sigmarp
          write (13,*) hanue1, size_yt, sigmay, sigmar
          write (14,*) hanue1, dive_yt, sigypr, sigmarp
        endif

c	#########################
c	coherent flux
c	#########################
	if (coherent.gt.0.) then
c	use formula due to K.J.Kim, LBL-22236 (march 1987)
c	published in PAC1987, p 194
c	converting into units photons/sec/0.1%BW I get
c	coherent Flux = .385 Photons/sec/0.1%BW
c	                * B(0,0)/(phot/sec/mm**2 mrad**2 0.1%BW)
c	                * (1/(hnue/eV)**2)

          flux(imesh,icurve)= .385*b1/(hanue1**2)
	endif

c Die Eneriebreite des Strahles hat kaum Einfluss auf den Fluss, deshalb
        flux(imesh,icurve)=flux(imesh,icurve)/dgamred

c	write (6,201)
c	1 hanue(imesh,icurve), ak1, smallk, brill(imesh,icurve)
c201	format (' hanue/eV =',f9.3,' K=',f9.3,' k=',f2.0,' B=',e10.4,
c	1 ' ph/(sec (mm*mrad)**2 .1%)')
c
c	total power in Watt according to K.J.Kim
	ptot = 1.774e-8*gamma**2*(aky1**2+akx1**2)*an*curr/alamd0

	dummy = 2.309+.5556*alog(.003223*aky1/alamd0)
	gap = 0.
	if (dummy.gt.0.) gap = alamd0*(1.519-sqrt(dummy))	

	write (10,203) hanue(imesh,icurve), brill(imesh,icurve),
	1   flux(imesh,icurve), ak1, ptot, gap
203	format (3x,e12.4,8x,e12.4,8x,e12.6,3x,f7.3,3x,e10.3,2x,f8.6)

200	continue
	icurve = icurve+1
	ncurve = ncurve+1

	if (sourceprop.gt.0.) then
	write (11, 211)
211	format (//,' hanue,   SIGMAX,    sigmax,    sigmar ')
	write (12, 212)
212	format (//,' hanue,  SIGMAX",   sigmax",    sigmar" ')
	write (13, 213)
213	format (//,' hanue,   SIGMAY,    sigmay,    sigmar ')
	write (14, 214)
214	format (//,' hanue,  SIGMAY",   sigmay",    sigmar" ')
	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif

	write (6,301)
301	format ('Enter another harmonic: ')
717     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 717
        else
          read(cline,*) smallk
        endif
	if (smallk.le.0.) return
	goto 109
	end
c	************************************************************
c				N O R M U E L L I
c	************************************************************
	subroutine normuelli (anorm, aky1,akx1)

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

c	calculate normalisation of brilliance
c	changed on march 28, 1988: phase space area is 2*pi * sigma * sigmap
c	units of anormx and anormy are  mm * mrad

c	##########################################################
c	Nov. 26, 1992:
c	Zur Normierung: Wir (JB und AG) beschliessen, ab sofort nur
c	noch Kurven zu erzeugen mit
c	sigmar  = (1/4pi) sqrt (lambda*L)
c	sigmarp =         sqrt (lambda/L)
c
c	Dies ist Konvention, da die transversalen Profile eindeutig
c	keine Gauss Kurven sind (siehe USEM Studien vom 11.11.92)
c	aber es gilt:
c	2 pi * sigmar * sigmarp = 1/2 * lambda
c	##########################################################
c
c	Normx = 2*pi * SIGMAX * SIGMAXP
c
c	e beam cross section
	sigmax = sqrt(epsx*betax)
c	e beam divergence
	sigxpr = sqrt(epsx/betax)

c	natural divergence of radiation
	alamda = (alamd0/(smallk*2.*gamma**2))*(1.+(aky1**2+akx1**2)/2.)
	aL = alamd0*an
	sigmarp = sqrt(alamda/(aL))

c	radiation cross section
c	This is FAITH!!!!!!!!!!!!!!!!!!!!!!!!
c	here we take K.J.Kim's convention
	sigmar = sqrt(alamda*al)/(4.*pi)

	anormx = 1.e6*2.*pi*sqrt((sigmax**2 + sigmar**2) *
	1                        (sigxpr**2 + sigmarp**2))

c	e beam cross section
	sigmay = sqrt(epsy*betay)
c	e beam divergence
	sigypr = sqrt(epsy/betay)

	anormy = 1.e6*2.*pi*sqrt((sigmay**2 + sigmar**2)*
	1	                 (sigypr**2 + sigmarp**2))

	anorm = anormx*anormy

	if (depthoffield.eq.1. ) then
c 	section added October 1992
c	This section is to include depth of field affects according to
c	K.J.Kim, LBL-22236, March 1987, published in PAC 1987, p.194
c	The philosophy is as follows:
c	the on-axis brightness is written as an integral along the
c	undulator axis. The contribution from a path element dz
c	is assumed to be a gaussian times a form factor g(z).
c	The free parameters are choosen to satisfy several criteria
c	which are:
c	- single electron flux is completely coherent
c	- spatial flux density is the brilliance integrated over angle
c	- angular flux density is the brilliance integrated over source area
c	- Total flux is double integral over brilliance

c	the last condition is met within 7%

c	2004/03/17: Comparing to WAVE. In one case numerical agreement
c	is obtained between WAVE (author: M.Scheer, BESSY) and BRILL
c	with depth of field effect. The depth of field effect reduces
c	the brilliance by a factor of 5 in that case.

	sigmar = sqrt(2.*alamda*al)/(4.*pi)
c----------------- corrected 18.3.2004
c	sigmarp = sqrt(alamda/2.*al)
	sigmarp = sqrt(alamda/(2.*al))
c----------------- 18.3.2004

c	this re-definition of sigmar and sigmarp has
c	consequences only for the depth of field effect calculation.
c	AG 2004/03/17

	anorm = 0.
c	4 points per periode
	nn = 4.*an
	dz = aL/nn
	do i=1,nn	
	z = -aL/2. + (i-0.5)*dz
c	g(z) = [ 1 + (2pi*z/L)**2] / 4L
	g = (1.+(2.*pi*z/aL)**2)/(4.*aL)
	Omegax = 2.*pi
	1	*sqrt((sigmax**2+sigmar**2)*(sigxpr**2+sigmarp**2) +
	2	(z*sigxpr*sigmarp)**2)
	Omegay = 2.*pi
	1	*sqrt((sigmay**2+sigmar**2)*(sigypr**2+sigmarp**2) +
	2	(z*sigypr*sigmarp)**2)
	anorm = anorm + dz * g / (Omegax * Omegay)

c	anorm = anorm + (1.+(2.*pi*z/aL)**2) /
c	1 (sqrt((sigmax**2+sigmar**2)*(sigxpr**2+sigmarp**2) +
c	2	(z*sigxpr*sigmarp)**2)
c	3 *sqrt((sigmay**2+sigmar**2)*(sigypr**2+sigmarp**2) +
c	4	(z*sigypr*sigmarp)**2))
	enddo
c	anorm = (2.*anorm*dz)/((2.*pi)**2*4.*aL)
	anorm = 1.e12/anorm
	endif			! of depth of field = 1

	if (depthoffield.eq.2.) then
c	this section is to treat the depth of field by superimposing
c	intensities from different path elements without any interference
c	not yet implemented 11/92 AG
	endif

	return
	end
c	*************************************************************
c				E F K A E L L I
c	*************************************************************
	real function efkaelli(aky,akx)
c	changed dimension to (0:11), 8.12.2004, JB
c	dimension B(0:9)
	dimension B(0:11)

c	the argument smallk is transmitted via COMMON
c	this is the function defined e.g. by Krinsky in Koch's book
c	eq. 276, p.151

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

c	cf. Elleaume in Undulator book, p.98, eq. 82
c	identical to Walker's Y except for Abs
	arg = (0.25*smallk*Abs(aKy**2-aKx**2))/(1.+(aKy**2+aKx**2)/2.)

c	use simplier CERN library routine  26.Oct 1995

c	increase dimension,(0:9) to (0:11), 8.12.2004, JB
	if (smallk.gt.11.) stop 'increase dimension B(0:11) in efkaelli'

	orderm = (smallk-1.)/2.
	orderp = (smallk+1.)/2.
c	order must be integer
	if (abs(orderm-int(orderm)) .gt. .001)
	1 write (6,*) ' in EFKA orderm not integer! orderm=',orderm
	if (abs(orderp-int(orderp)) .gt. .001)
	1 write (6,*) ' in EFKA orderp not integer! orderp=',orderp
c	we always need integer order of Bessel function
c	0. .le. A .lt. 1. is the order of Bessel function in the SR call
	A = 0.
	Nm = int(orderm)
	Np = int(orderp)
	if (Np.gt.11) stop 'increase dimension B(0:11) in efkaelli'
c	call BSJA (arg, A, Np, 6, B)
c	we don't have to increase 5->6 ???  04/12/14 AG
        print*,"BSJA:",arg,A,np,B
	call BSJA (arg, A, Np, 5, B)

c	?????????????????????????????
c	In case of an APPLE undulator we know that the relative phase
c	between Kx and Ky is always pi/2.
c	Pascal Elleaume is calculating the flux for a given polarisation
c	direction.
c	When I calculate the flux with the corresponding circular
c	polarisation I obtain a simple sum: (Ky+Kx)**2
c	thus this is the flux of circularly polarised radiation and not
c	elliptically polarised!! This difference may not be very
c	important as long as we are close to the maximum of the degree of
c	circular polarisation.
c	But users may want to use the code in strange situations
c	the other orthogonal polarisation, i.e. circular polarisation
c	with opposite helicity to the magnetic field yields (Ky-Kx)**2
c	the total flux irrespective of polarisation, i.e. S0 is proportional to
c	(Ky+Kx)**2 + (Ky-Kx)**2 = 2*(Ky**2+Kx**2)

c	Walker has   Ax = Ky (B(Np)-B(Nm))
c	and          Ay = i Kx (B(Np)+B(Nm))
c	(eq. A4ff)

c	Elleaume has Fn ~ (B(Np)-B(Nm))**2 for both polarisations
c	(eq. 86)

c	?????????????????????????????
c
cd	the unit polarisation vector is (ex+iey)/sqrt(2). This explains
c	the factor 1/2 in (aky+akx)**2/2
c	Elleaume
c	efkaelli = (smallk/(1.+(aKy**2+aKx**2)/2.))**2*
c	1	(aKy+aKx)**2/2.*(B(Np)-B(Nm))**2

c	there remains a difference with respect to Walker for Kx != Ky
c	on top Walker has (B(Np)+B(Nm)) which is difficult to understand.
c	a numerical study by WAVE by MSH indicates better agreement with
c	Walker's expression than with Elleaume's expression. Both
c	formulas agree for Kx/Ky=0 and Kx/Ky=1.
c	(note: for  Kx/Ky=1 only B(Np) contributes because the argument
c	of the Bessel functions vanishes)
c	The expressions differ by 7% for Kx/Ky=.5 and Walker agrees with WAVE
c	JB 2004/12/20
c	Note: The comparision is done for the elliptically polarised
c	intensity and not for the circularly polarised intensity!
c	Walker appendix
c	good only for Kx/Ky <= 1.

	Ax = aKy*(B(Np)-B(Nm))
	Ay = aKx*(B(Np)+B(Nm))
	efkaelli = (smallk/(1.+(aKy**2+aKx**2)/2.))**2*
	1	(Ay**2+Ax**2)

c	write (6,777)   Nm, Np, B(Nm), B(Np), efkaelli
c777	format ( ' in SR EFKA:',/,
c	1 ' Nm=',I4,' Np=',I4,' J(Nm)=',e10.4,' J(Np)=',e10.4,
c	2 ' efkaelli = ',e10.4)

	return
	end
	
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.15.03  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.32.49  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	******************************************************************
c				L A B E L
c	******************************************************************
	subroutine label (iflag)
c	male den Label fuer eine der Kurven

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	xpap = 20.5
	ypap = 29.5
	if (iflag.gt.1) goto 500
c	Einlesen von neuen Labels
	nlab = nlab+1
	write (6,101)
101	format (' enter label, ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,105) llab(nlab)
        endif
105	format (a8)
	write (6,102)
102	format ('  and coordinates in cm ',/
	2' (first (x,y) for brilliance plot, then (x,y) for flux plot): ')
779      read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) xlabb(nlab), ylabb(nlab), xlabf(nlab), ylabf(nlab)
        endif
150	continue
	if (nlab.gt.100) nlab = 100
	return

500 	continue
	do 599 ilab=1,nlab

	nheight = labsize(ilab)
	richt   = -90.
	nx1     = int(100.*(      yorig+ylabb(ilab)))
	ny1     = int(100.*(xpap-(xorig+xlabb(ilab))))
	write (1,501) nx1,ny1,nheight,richt,(lab(j,ilab),j=1,8)
501	format ('label   ',3i4,f4.0,8a1)
	nx1     = int(100.*(      yorig+ylabf(ilab)))
	ny1     = int(100.*(xpap-(xorig+xlabf(ilab))))
	write (2,501) nx1,ny1,nheight,richt,(lab(j,ilab),j=1,8)
599	continue
	return
	end
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.11.21  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.33.28  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	****************************************************************
c				L A Y O U T
c  	****************************************************************
	subroutine layout

c	establishes details of layout

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	size = 0.4
c	size = 0.5
100	continue
	write (6,101) xpap, ypap, xorig, yorig,
	1	 xr1, xr2, yr1, yr2, size
101	format (/
	1' 0) exit                        ',/
	2' 1) Papiergroesse xpap (cm)     ',f10.2,/
	2' 2)               ypap (cm)     ',f10.2,/
	3' 3) Ursprung      xorigin (cm)  ',f10.2,/
	4' 4)               yorigin (cm)  ',f10.2,/
	5' 5) Rand  links von links xr1   ',f10.2,/
	6' 6) Rand rechts von links xr2   ',f10.2,/
	7' 7) Rand unten von unten  yr1   ',f10.2,/
	8' 8) Rand  oben von unten  yr2   ',f10.2,/
	9' 9) Groesse der Label     size  ',f10.2,/
	1/'$Waehle Zahl und neuen Wert: ')

778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) iflag
        endif
	if (iflag.le.0) return
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) dummy
        endif
	goto (1,2,3,4,5,6,7,8,9) iflag
1	xpap    = dummy
	goto 100
2	ypap    = dummy
	goto 100
3	xorig   = dummy
	goto 100
4	yorig   = dummy
	goto 100
5	xr1     = dummy
	goto 100
6	xr2     = dummy
	goto 100
7	yr1     = dummy
	goto 100
8	yr2     = dummy
	goto 100
9	size    = dummy
	do 91 ilabel = 1, nlab
91	labsize (ilabel) = int (100.*size)
	goto 100
	end
*CMZ :  0.01/00 03/09/2013  13.20.32  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.34.01  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	************************************************************
c	                       P O W E R
c	************************************************************

	subroutine power (N, K, lambda0, gamma, I, Ptot, d2P)

	real*4 N, K, lambda0, I
c	N	number of periodes
c	K	wiggler strength
c	lambda0	wiggler period
c	gamma	electron energy
c	I	electron current
c	Ptot	total power
c	d2P	power density

	pi = 4.*atan(1.)
	Z_0 = 377.
	e   = 1.6e-19
	c   = 3.0e8

	Ptot = N/6. * Z_0 * I * 2.*pi*e*c/lambda0 * (gamma*K)**2
	GK   = K*(K**6+24./7.*K**4+4.*K**2+16./7.)/(1.+K**2)**3.5
c	only on axis power density
	d2P  = Ptot*21.*gamma**2/(16.*pi*K)*GK

	write (6,101) N, K, lambda0, gamma, I, Ptot, d2P
101	format (' call to subroutine power with parameters',/
	1' N                    ',e12.6,/
	2' K                    ',e12.6,/
	3' lambda0              ',e12.6,/
	4' gamma                ',e12.6,/
	5' I                    ',e12.6,/
	6' Ptot                 ',e12.6,/
	7' d2P                  ',e12.6,//
	8' ',//)
cmsh	8' No match to subroutine power !!!!!!!!',//)

	return
	end
*CMZ :  1.02/01 16/07/2015  09.01.06  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.09  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.11.05  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.34.31  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	subroutine ring
c	*************************************************************
c                          R I N G
c	*************************************************************

        character a

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	goto 2
1	write (6,10)
10	format ('Enter electron energy in units of rest mass: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) gamma
          open(unit=999,file='brill.ebeam',status='unknown')
          write(999,*)gamma*0.000511
          close(999)
        endif
	write (6,12)
12	format (' enter fractional energy spread: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) dgamma
        endif
	write (6,20)
20	format (' enter horizontal and vertical electron beam emittance',/
	1 '$ in units of  pi m rad: ')
710       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) epsx, epsy
        endif
	write (6,30)
30	format ('Enter average stored electron current in Amps: ... ')
711     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          read(cline,*) curr
          open(unit=999,file='brill.curr',status='unknown')
          write(999,*)curr
          close(999)
        endif
	write (6,40)
40	format ('Enter bending radius in meter: ...................')
712     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 712
        else
          read(cline,*) rho
        endif
2	write (6,100) gamma, dgamma, rho, epsx, epsy, curr
100	format (//' storage ring parameters',/
	1          ' -----------------------',/
	1 ' electron energy in units of rest mass   gamma = ',f10.1,/
	2 ' fractional energy spread               dgamma = ',e10.4,/
	2 ' bending radius (meter)                    rho = ',f10.5,/
	3 ' horizontal emittance ( Pi m rad )        epsx = ',e10.4,/
	4 ' vertical   emittance                     epsy = ',e10.4,/
	5 ' stored current (Ampere)                  curr = ',f10.3,/
	9 // '$correct? (y/n)    ')
713       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          a=cline(1:1)
        endif
cmsh	read(5,101) a
101	format (a1)
        if (a.eq.'y') goto 90
	goto 1
90	continue
	write (10,91) gamma, dgamma, rho, epsx, epsy, curr
91	format (/
	1 ' Ring Parameters',/
	1 ' ---------------',/
	2 5x,' gamma			',f10.0,/
	2 5x,' dgamma			',e10.4,/
	3 5x,' rho (meter)		',f10.4,/
	4 5x,' epsx ( pi m rad)    	',e10.4,/
	5 5x,' epsy ( pi m rad) 	',e10.4,/
	6 5x,' current (amp)		',f10.4,/)
	return
	end
*CMZ :  1.02/01 21/11/2017  14.03.10  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.09  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.14.24  by  Michael Scheer
*CMZ :  1.01/00 12/09/2013  15.01.41  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.10.31  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.17  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	****************************************************************
c			      version Jan 1994
c				U N D U N E W
c	****************************************************************
c	subroutine undula (ndatfil)
	subroutine undula


c	2008/09/19
c	the following explanation is included to help comparing to
c	other expressions
c
c	outline:
c	0) the spectral angular brightness of an undulator is obtained
c		for zero emittance and zero energy spread
c	1) the spectral angular brightness is corrected for
c		ebeam energy spread
c	2) the spectral flux is obtained by multiplying the spectral
c		brightness with the opening angle of undulator radiation
c	3) the spectral brilliance is obtained by dividing the spectral flux
c		by the effective transverse transverse phase source
c		space volume
c		

c	We have (this defines quantities)
c	lamda_n = lamda_0/(n * 2*gamma^2) * (1+(1/2) K^2)
c	(below we label the undulator harmonics by smallk)
c	We deal with a planar undulator
c	K ( = ak1) is the ordinary wiggler parameter

c	the spectral brightness is ( xray data booklet eq. 13)
c	B1 = alpha N^2 gamma^2 (Delta omega/omega) (I/e) F_n(K)

c	= 1.744*10^14 N^2 (E/GeV)^2 (I/Amp) F_n(K) photons/(sec mrad^2 0.1%BW)

c	We take spectral distribution of undulator radiation to be a
c	normal distribution with fractional width
c	sigrph = 1 / (2.6 N n)
c	sigrph = 1 / (2 * sqrt(2)) N n) !msh
c
c	(this is an attempt to match a (sin x/x)^2 spectral distribution
c	curve to a normal distribution)

c	the ebeam energy distribution is taken to be normal with
c	fractional width
c	(delta gamma) = (Delta gamma) / gamma

c	the spectral brightness corrected for ebeam energy spread is

c	B1(ebeam) = B1 * sigrph/ sqrt(sigrph^2 + (2*delta gamma)^2)
c	= B1/sqrt(1. + (2* delta gamma * n * 2.6 * N)^2 ) !msh
c	= B1/sqrt(1. + (2* delta gamma * n * 2.0*sqrt(2.0) * N)^2 )

c	the factor 2 in front of delta gamma is due to the fact that the
c	ebeam energy enters the photon energy as gamma^2

c	the spectral flux is
c	flux = 2 pi B1 * sigrp^2
c	with sigrp = sqrt(lamda/L)
c	Note: sigrp is the angular width of the undulator radiation
c	      sigrph is the spectral width of the undulator radiation

c	the spectral brilliance is
c	B = flux / Norm
c	where Norm = Normx * Normy
c	with
c	Normx = sqrt((sigmax^2+sigmar^2)*(sigmarp^2+sigmaxp^2))
c	Normy = sqrt((sigmay^2+sigmar^2)*(sigmarp^2+sigmayp^2))
c	here
c	sigmax  = sqrt(epsilonx*betax)
c	sigmaxp = sqrt(epsilonx/betax)
c	sigmar  = (1/4 pi) sqrt(lamda*L)
c	sigmarp = sqrt(lamda/L)
c	note: sigmar*sigmarp = lamda/(4pi)

c	2008/09/19 Andreas Gaupp



c	normalization changed according to K.J.Kim  April 11, 1988

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

1	write (6,20)
20	format ('$ enter periode length in m: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alamd0
        endif
	write (6,10)
10	format ('$ enter number of undulator periodes: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) an
        endif
	if ( an*alamd0.ge.alength) write (6,21) alamd0, an, alength
21	format (' **** warning ****',/
	1 ' lamda0  = ',f10.2,' m',/
	2 '      N  = ',f10.0,/
	3 ' alength = ',f10.2,' m '/
	4 ' *** not compatible ***'/)
	
	write (6,30)
30	format('$ enter range in wiggler strength, KMIN and KMAX: ')
710     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) akmin, akmax
        endif
	if (akmin.ge.0. .and. akmax.ge.akmin) goto 32
	write (6,31) akmin, akmax
31	format (' **** warning ****',/
	1 ' AKMIN   = ',f10.3,'   AKMAX  = ',f10.3,' *** inconsistent'/)
32	write (6,40)
40	format ('$ enter order of undulator spectrum (1, 3, 5,..): ')
711     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          read(cline,*) smallk
        endif
	if (smallk.lt.1. .or.
	1   abs(smallk-float(int(smallk))).gt.1.e-4) write (6,41) smallk
41	format (' **** warning ****',/
	1 ' bad value for order of spectrum. SMALLK = ',f10.3)
2	continue

c	get minimum gap for hybrid-REC structure using Halbach's formula
c	for REC without 10% reduction
	gapmin = alamd0*(1.519-sqrt(2.309+.5556*alog(.003223*akmax/alamd0)))	
	write (6,100) an, alamd0, akmin, akmax, smallk, gapmin
100	format (/,' undulator parameters',/
	1 ' number of periodes             AN     = ',f10.0,/
	2 ' periode length (m)             alamd0 = ',f10.5,/
	3 ' minimum K                      KMIN   = ',f10.2,/
	4 ' maximum K                      KMAX   = ',f10.2,/
	5 ' order of spectrum              smallk = ',f10.0,/
	6 ' minimum gap (hybrid REC) (m)   gapmin = ',f10.5,//
	9 '$ correct? (y/n)  ')
712       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 712
        else
          a=cline(1:1)
        endif
cmsh	read(5,101) a
101	format (a1)
	if (a.eq.'y') goto 109
	goto 1
109	continue

c	plot undulator spectrum
	write (6,110)
110	format (' enter photon energy range (eV) of this curve',/
	1 '$ hanuemin and hanuemax:  ')
713       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          read(cline,*) hanuemin, hanuemax
        endif
	if (hanuemax.le.0.) then
		hanuemax = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*(1.+0.5*aKmin**2)))
	endif
	if (hanuemin.le.0.) then
		hanuemin = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*(1.+0.5*aKmax**2)))
	endif
	write (6,120)
120	format ('Enter number of points on this curve, mesh:  ')
714     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 714
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256

c	calculate normalisation of brilliance for K = 0
c	only for initial printout, not for final curve
	theta0 = 1./(gamma*sqrt(an*smallk))
c	note alength is in m
	sigmax = sqrt(epsx*betax)
  	sigmay = sqrt(epsy*betay)
	sigmaxp= sqrt(theta0**2+epsx/betax)
	sigmayp= sqrt(theta0**2+epsy/betay)
	write (6,121) theta0, sigmax, sigmay, sigmaxp, sigmayp, dgamma
121	format (/
	1 ' theta0 (rad, K=0) = ',e10.4,/
	2 ' sigmax (m)        = ',e10.4,/
	3 ' sigmay (m)        = ',e10.4,/
	4 ' sigmaxp (rad) K=0 = ',e10.4,/
	5 ' sigmayp (rad) K=0 = ',e10.4,/
	6 ' sigma E           = ',e10.4,//)
	
c	get maximum power according to K.J.Kim
	call power (an, aKmax, alamd0, gamma, curr, Ptot, d2P)

	write (10,122) alamd0, an, akmin, akmax, smallk, gapmin,
	1              theta0,sigmax, sigmay,sigmaxp,sigmayp,dgamma
	write (10,123) Ptot, d2P

122	format (/
	1 ' properties of undulator',/
	1 ' -----------------------',/
	2 5x,' periode length (m)             	',f10.4,/
	3 5x,' number of periodes              	',f10.4,/
	4 5x,' minimum K                       	',f10.4,/
	5 5x,' maximum K                        ',f10.4,/
	6 5x,' order of spectrum k              ',f10.4,/
	7 5x,' minimum gap (REC hybrid, m)      ',f10.4,//
	1 ' properties of radiation',/
	2 5x,' theta0 (rad) K = 0               ',e10.4,/
	3 5x,' sigmax (m)                       ',e10.4,/
	4 5x,' sigmay (m)                       ',e10.4,/
	5 5x,' sigmaxp (rad) K = 0              ',e10.4,/
	6 5x,' sigmayp (rad) K = 0              ',e10.4,/
	7 5x,' dgamma                           ',e10.4)
123	format (
	7 5x,' Ptot (Watt)                      ',e10.4,/
	7 5x,' dP/dOmega (Watt/mrad^2)          ',e10.4,/
	9 5x,' h*nue (eV)    Brilliance         Flux   Flux density    K
	1       Ptot    gap',/
	2 5x '                  photon         photon      photon
	3     (Watt)     (m)',/
	4 5x '     sec*(mm*mrad)**2*0.1% BW  sec*0.1%BW  sec*0.1%BW*mrad^2'
	5 ,/)

	if (sourceprop.gt.0.) then
	call fileopen
	write (11,*) mesh(icurve)-1
	write (12,*) mesh(icurve)-1
	write (13,*) mesh(icurve)-1
	write (14,*) mesh(icurve)-1

c	e beam cross section and divergence
	sigmax = sqrt(epsx*betax)
	sigmay = sqrt(epsy*betay)
	sigxpr = sqrt(epsx/betax)
	sigypr = sqrt(epsy/betay)

	endif

	do 200 imesh=1,mesh(icurve)
	hanue1 = hanuemin+(hanuemax-hanuemin)*(imesh-1)/float(mesh(icurve)-1)
c	8065 cm**-1  = 1 eV
	ak1 = 2.*(smallk*2.*gamma**2/(hanue1*alamd0*806500.)-1.)
	if (ak1.le.0.) goto 200
	ak1 = sqrt(ak1)
C	this ak1 is wiggler constant from
c	lamda_n = lamda_0/(2.*gamma**2)*(1+(1./2.)*K**2)
c	i.e. planar undulator


	if (ak1.lt.akmin .or. ak1.gt.akmax) goto 200

c	formula due to Krinsky in Koch's book, eq.281, p.152
c	checked on April 9, 1985
c	units are  photons/[sec*mrad**2*0.1%BW*Amp])
c	our efka is the same as F_n(K) in xray data booklet
c	this is the same as xray data booklet eq. 13
c	b1 is on axis spectral brightness
c	units: curr in Amp
c	       b1 in photons/(sec mrad^2 0.1% )

c	b1 = 4.56E7*(an*gamma)**2*curr*efka(smallk,ak1)
	b1= 4.56E7*(an*gamma)**2*curr*efka(ak1)

c	this is the spectral brightness, which is reduced by a
c	finite energy spread
c	 here we take the spectral distribution of undulator radiation
c	to be a normal with fractional width sigma = 1./(smallk*2.6*N)
c	to be a normal with fractional width sigma = 1./(smallk*2.0*sqrt(2.0)*N) !MSH
c	(not 1/(2*N) !!!!!!!!!!!!!!

c	the distribution of ebeam energy is taken to be a normal
c	distribution with fractional width dgamma = Delta gamma / gamma
c	The energy spread of the e beam is multiplied by 2 because
c	gamma**2 enters the energy formula
c	This is necessary as discovered by Howard Padmore at ALS in
c	winter 1993/94.                         AG 1994/01/16

cmsh	b1= b1/sqrt(1.+(2.*dgamma*smallk*2.6*aN)**2)
	b1= b1/sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2)
        dgamred=1./sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2) !msh
c	save on axis spectral brightness into array fluxdens
c	need to correct for electron beam divergence
c	e beam divergence
	sigxpr = sqrt(epsx/betax)
	sigypr = sqrt(epsy/betay)
c	natural divergence of radiation
	alamda = (alamd0/(smallk*2.*gamma**2))*(1.+ak1**2/2.)
	aL = alamd0*an
	sigmarp = sqrt(alamda/(aL))

        print*,'*** Fluxdens unclear in undunew! ***'
	fluxdens(imesh,icurve) = b1/
	1 sqrt(1.+(sigxpr**2+sigypr**2)/sigmarp**2)

c	flux(imesh,icurve) is the total flux at energy curvex
c	here we make a very crude approximation,
c	since we calculate the flux by multiplying the angular brightness
c	by 2*pi*theta0**2
c	our formula is the same as eq.17 in xray data booklet
c	the formula for sigmar is consistent with Walker
c	ESRP - IRM - 54/84
c	***** sigmar = sqrt (lambda / L) *****
c	Note: the quantity called sigmar here should better be called
c	sigmarp


c	J.B. 30.5.89
c	this is sigmar = sqrt (n lamda / L)
	sigmar = 1./gamma*sqrt((1.+ak1**2/2.)/(smallk*2.*an))
	flux(imesh,icurve) = 2.*pi*b1*(sigmar*1000.)**2
	hanue (imesh, icurve) = hanue1
	call normu (anorm, ak1)
	if (imesh.gt.256 .or. icurve.gt.100 .or. anorm.le.0.)
	1	write (6,*) ' *** imesh=',imesh,'   icurve=',icurve,
	2                   ' *** anorm=',anorm
	b1= flux(imesh,icurve)/anorm
	brill (imesh,icurve) = b1

	if (sourceprop.gt.0.) then
c	natural divergence of radiation
	alamda = (alamd0/(smallk*2.*gamma**2))*(1.+ak1**2/2.)
	aL = alamd0*an
	sigmarp = sqrt(alamda/(aL))
c	radiation cross section
	sigmar = sqrt(alamda*al)/(4.*pi)
c	MSH behauptet basierend auf einem Vergleich der zweiten Momente
c	(2000/06):
c	sigmar=(1/2 pi) sqrt(2 lambda L)
c	sigmar * sigmarp = (1/2 pi) lambda (*0.82)

	size_xt = sqrt(sigmax**2+sigmar**2)
	dive_xt = sqrt(sigxpr**2 + sigmarp**2)
	size_yt = sqrt(sigmay**2+sigmar**2)
	dive_yt = sqrt(sigypr**2 + sigmarp**2)

	write (11,*) hanue1, size_xt, sigmax, sigmar
	write (12,*) hanue1, dive_xt, sigxpr, sigmarp
	write (13,*) hanue1, size_yt, sigmay, sigmar
	write (14,*) hanue1, dive_yt, sigypr, sigmarp
	endif



c	#########################
c	coherent flux
c	#########################
	if (coherent.gt.0.) then
c	use formula due to K.J.Kim, LBL-22236 (march 1987)
c	published in PAC1987, p 194
c	converting into units photons/sec/0.1%BW I get
c	coherent Flux = .385 Photons/sec/0.1%BW
c	                * B(0,0)/(phot/sec/mm**2 mrad**2 0.1%BW)
c	                * (1/(hnue/eV)**2)
				
c	attention!!!!!!!!!!!!!!!!!!!!!!!!!1
	flux(imesh,icurve) = .385*b1/(hanue1**2)
	endif

c Die Eneriebreite des Strahles hat kaum Einfluss auf den Fluss, deshalb
        flux(imesh,icurve)=flux(imesh,icurve)/dgamred

c	write (6,201)
c	1 hanue(imesh,icurve), ak1, smallk, brill(imesh,icurve)
c201	format (' hanue/eV =',f9.3,' K=',f9.3,' k=',f2.0,' B=',e10.4,
c	1 ' ph/(sec (mm*mrad)**2 .1%)')
c
c	total power in Watt according to K.J.Kim
	ptot = 1.774e-8*(gamma*ak1)**2*an*curr/alamd0

	dummy = 2.309+.5556*alog(.003223*ak1/alamd0)
	gap = 0.
	if (dummy.gt.0.) gap = alamd0*(1.519-sqrt(dummy))	

	write (10,203) hanue(imesh,icurve), brill(imesh,icurve),
	1   flux(imesh,icurve), fluxdens(imesh,icurve), ak1, ptot, gap
203	format (3x,e12.4, 3x,e12.4, 3x,e12.6, 3x,e12.6,
	1 3x,f7.3, 3x,e10.3, 2x,f8.6)


c	hanl1 = alog10(hanue1)
c	hanl2 = alog10(hanue2)
c	blog1 = alog10(b1)
c	blog2 = alog10(b2)
c	if (hanl1.lt.hnminl  .or. hanl2.gt.hnmaxl ) goto 200
c	if (blog1.gt.bmaxl   .or. blog2.gt.bmaxl  ) goto 200
c	if (blog1.lt.bminl   .or. blog2.lt.bminl  ) goto 200
c	y1 = yscale*(blog1-bminl)
c	x1 = xscale*(hanl1-hnminl)
c	y2 = yscale*(blog2-bminl)
c	x2 = xscale*(hanl2-hnminl)
c	if (y1.lt.0. .or. y2.lt.0.) goto 200
c	call vector (x1,y1, x2,y2)

200	continue
	icurve = icurve+1
	ncurve = ncurve+1

	if (sourceprop.gt.0.) then
	write (11, 211)
211	format (//,' hanue,   SIGMAX,    sigmax,    sigmar ')
	write (12, 212)
212	format (//,' hanue,  SIGMAX",   sigmax",    sigmar" ')
	write (13, 213)
213	format (//,' hanue,   SIGMAY,    sigmay,    sigmar ')
	write (14, 214)
214	format (//,' hanue,  SIGMAY",   sigmay",    sigmar" ')
	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif

	write (6,301)
301	format ('Enter another harmonic: ')
715     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 715
        else
          read(cline,*) smallk
        endif
	if (smallk.le.0.) return
	goto 109
	end
c	************************************************************
c				N O R M U
c	************************************************************
	subroutine normu (anorm, ak1)

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

c	calculate normalisation of brilliance
c	changed on march 28, 1988: phase space area is 2*pi * sigma * sigmap
c	units of anormx and anormy are  mm * mrad

c	##########################################################
c	Nov. 26, 1992:
c	Zur Normierung: Wir (JB und AG) beschliessen, ab sofort nur
c	noch Kurven zu erzeugen mit
c	sigmar  = (1/4pi) sqrt (lambda*L)
c	sigmarp =         sqrt (lambda/L)
c
c	Dies ist Konvention, da die transversalen Profile eindeutig
c	keine Gauss Kurven sind (siehe USEM Studien vom 11.11.92)
c	aber es gilt:
c	2 pi * sigmar * sigmarp = 1/2 * lambda
c	##########################################################
c
c	Normx = 2*pi * SIGMAX * SIGMAXP
c
c	e beam cross section
	sigmax = sqrt(epsx*betax)
c	e beam divergence
	sigxpr = sqrt(epsx/betax)

c	natural divergence of radiation
	alamda = (alamd0/(smallk*2.*gamma**2))*(1.+ak1**2/2.)
	aL = alamd0*an
	sigmarp = sqrt(alamda/(aL))

c	radiation cross section
c	This is FAITH!!!!!!!!!!!!!!!!!!!!!!!!
c	here we take K.J.Kim's convention
	sigmar = sqrt(alamda*al)/(4.*pi)

	anormx = 1.e6*2.*pi*sqrt((sigmax**2 + sigmar**2) *
	1                        (sigxpr**2 + sigmarp**2))

c	e beam cross section
	sigmay = sqrt(epsy*betay)
c	e beam divergence
	sigypr = sqrt(epsy/betay)

	anormy = 1.e6*2.*pi*sqrt((sigmay**2 + sigmar**2)*
	1	                 (sigypr**2 + sigmarp**2))

	anorm = anormx*anormy

	if (depthoffield.eq.1. ) then
c 	section added October 1992
c	This section is to include depth of field affects according to
c	K.J.Kim, LBL-22236, March 1987, published in PAC 1987, p.194
c	The philosophy is as follows:
c	the on axis brightness is written as an integral along the
c	undulator axis. The contribution from a path element dz
c	is assumed to be a gaussian times a form factor g(z).
c	The free parameters are choosen to satisfy several criteria
c	which are:
c	- single electron flux is completely coherent
c	- spatial flux density is the brilliance integrated over angle
c	- angular flux density is the brilliance integrated over source area
c	- Total flux is double integral over brilliance

c	the last condition is met within 7%

c	2004/03/17: Comparing to WAVE. In one case numerical agreement
c	is obtained between WAVE (author: M.Scheer, BESSY) and BRILL
c	with depth of field effect. The depth of field effect reduces
c	the brilliance by a factor of 5 in that case.

	sigmar = sqrt(2.*alamda*al)/(4.*pi)
c----------------- corrected 18.3.2004
c	sigmarp = sqrt(alamda/2.*al)
	sigmarp = sqrt(alamda/(2.*al))
c----------------- 18.3.2004

c	this re-definition of sigmar and sigmarp has
c	consequences only for the depth of field effect calculation.
c	AG 2004/03/17

	anorm = 0.
c	4 points per periode
	nn = 4.*an
	dz = aL/nn
	do i=1,nn	
	z = -aL/2. + (i-0.5)*dz
c	g(z) = [ 1 + (2pi*z/L)**2] / 4L
	g = (1.+(2.*pi*z/aL)**2)/(4.*aL)
	Omegax = 2.*pi
	1	*sqrt((sigmax**2+sigmar**2)*(sigxpr**2+sigmarp**2) +
	2	(z*sigxpr*sigmarp)**2)
	Omegay = 2.*pi
	1	*sqrt((sigmay**2+sigmar**2)*(sigypr**2+sigmarp**2) +
	2	(z*sigypr*sigmarp)**2)
	anorm = anorm + dz * g / (Omegax * Omegay)

c	anorm = anorm + (1.+(2.*pi*z/aL)**2) /
c	1 (sqrt((sigmax**2+sigmar**2)*(sigxpr**2+sigmarp**2) +
c	2	(z*sigxpr*sigmarp)**2)
c	3 *sqrt((sigmay**2+sigmar**2)*(sigypr**2+sigmarp**2) +
c	4	(z*sigypr*sigmarp)**2))
	enddo
c	anorm = (2.*anorm*dz)/((2.*pi)**2*4.*aL)
	anorm = 1.e12/anorm
	endif			! of depth of field = 1

	if (depthoffield.eq.2.) then
c	this section is to treat the depth of field by superimposing
c	intensities from different path elements without any interference
c	not yet implemented 11/92 AG
	endif
	return
	end
c	*************************************************************
c				E F K A
c	*************************************************************
c	real function efka(smallk, ak)
	real function efka(ak)
c	changed dimension to (0:11), 8.12.2004, JB
c	dimension B(0:9)
	dimension B(0:11)

c	the argument smallk is transmitted via COMMON
c	this is the function defined e.g. by Krinsky in Koch's book
c	eq. 276, p.151

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	arg = (smallk*(ak/2.)**2)/(1.+ak**2/2.)
c	use CERN routine BSJA
c	tested 30.10.95 AG
	goto 999

c	ord1 = (smallk-1.)/2.
c	ord2 = (smallk+1.)/2.
c	if (abs(ord1-float(int(ord1))).ge.1.e-4  .or.
c	1   abs(ord2-float(int(ord2))).ge.1.e-4)
c	1                      write (6,10) smallk, ord1,ord2
c10	format (/,' **** warning ****',/
c	1 ' smallk  = ',f10.2,/
c	2 ' ord1    = ',f10.2,/
c	3 ' ord2    = ',f10.2,/
c	4 ' *** not ok ***')
c	n = int (max(ord1,ord2))
c	x    = arg
c	y    = 0.
c	alpha= 0.
c	beta = 0.
c	call combes (x, y, alpha, beta, n, bjre, bjim, yre, yim, ny)
c	n1 = ord1+1.
c	n2 = ord2+1.
c	FK = smallk*ak/(1.+ak**2/2.)*(bjre(n1)-bjre(n2))
c	FK = FK**2
cc	write (6,101) smallk, ak, fk
cc101	format (' k = ',f3.0,' K = ',f6.2,' FK = ',e10.4)
c	efka = fk
c	return

c	####################################################
999	continue
c	use simplier CERN library routine  26.Oct 1995

c	increased dimension,(0:9) to (0:11), 8.12.2004, JB
c	if (smallk.gt.9.) stop 'increase dimension B(0:9) in efka'
	if (smallk.gt.11.) stop 'increase dimension B(0:11) in efka'

	orderm = (smallk-1.)/2.
	orderp = (smallk+1.)/2.
c	order must be integer
	if (abs(orderm-int(orderm)) .gt. .001)
	1 write (6,*) ' in EFKA orderm not integer! orderm=',orderm
	if (abs(orderp-int(orderp)) .gt. .001)
	1 write (6,*) ' in EFKA orderp not integer! orderp=',orderp
c	we always need integer order of Bessel function
c	0. .le. A .lt. 1. is the order of Bessel function in the SR call
	A = 0.
	Nm = int(orderm)
	Np = int(orderp)
	if (Np.gt.11) stop 'increase dimension B(0:11) in efka'
        print*,"BSJA:",arg,A,np,B
	call BSJA (arg, A, Np, 5, B)

	efka = ( (smallk*ak/(1.+ak**2/2.)) * (B(Nm)-B(Np)) )**2

c	write (6,777)   Nm, Np, B(Nm), B(Np), efka
c777	format ( ' in SR EFKA:',/,
c	1 ' Nm=',I4,' Np=',I4,' J(Nm)=',e10.4,' J(Np)=',e10.4,
c	2 ' efka = ',e10.4)

	return
	end
*CMZ :  2.00/00 06/12/2017  17.02.44  by  Michael Scheer
*CMZ :  1.02/01 21/11/2017  14.02.35  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.09  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.09.35  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.56  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	subroutine wiggler

c	********************************************************************
c			    version from Jan 1994
c			        W I G G N E W
c	********************************************************************
c	this is reconstructed 12/95
c	changed to avaoid K_(cline/3)

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

1	write (6,20)
20	format ('$ enter length of periode (m) ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alamd0
        endif
	write (6,10)
10	format ('$ enter number of periodes: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) an
        endif
	if(an*alamd0.gt.alength) write (6,21) an,alamd0,alength
21	format (' **** warning ****',/
	1 ' number of periodes                AN = ',f10.3,/
	2 ' length of periode (m)        alamd0  = ',f10.3,/
	3 ' length of long straight (m)  alength = ',f10.3,/
	4 ' *** inconsisten ***',/)
	write (6,30)
30	format ('$ enter wiggler strength AK:  ')
710     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) AK
        endif
				
2	continue
c	characteristic wavelength in m
c	checked on April 4, 1985
c	numerically the same as  hanuec=1.73e-8(eV/Gauss)*gamma**2/B
	alamdc = (2./3.)*alamd0/(gamma**2*ak)
	hanuec = 1./(alamdc*806500.)
	B_0  = AK*2.*pi*.911e-27*9.e20/(4.8e-10*alamd0*100.)

	write (6,100)  alamd0, an, ak, hanuec, B_0
100	format (/,' wiggler parameter ',/
	2 ' periode length (m)          alamd0 = ',f10.3,/
	1 ' number of periodes              an = ',f10.3,/
	3 ' wiggler strength                AK = ',f10.3,/
	4 ' critical photon energy (eV) hanuec = ',f10.3,/
	4 ' Magnetic field amplitude (Gauss)   = ',f10.3,//
	5 '$correct ? (y/n)  ')
711       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          a=cline(1:1)
        endif
cmsh	read(5,101)  a
101	format (a1)
	if (a.ne.'y' .and. a.ne.'Y') goto 1

c	calculate wiggler spectrum
	write (6,110)
110	format ('Enter photon energy range (in eV), hanuemin, hanuemax:  ')
712     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 712
        else
          read(cline,*) hanuemin, hanuemax
        endif
	write (6,120)
120	format ('$ enter number of points in this range: ')
713     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256

c	total power in Watt according to K.J.Kim, NIM A246, 67 (1986)
	call power (aN, aK, alamd0, gamma, curr, Ptot, dPdOmega)
c	Ptot = (an/6.)*377.*curr*(2.*pi*1.6e-19*3e8)/alamd0
c	1	*gamma**2*aK**2

c	Ptot = 1.774e-8*(gamma*ak)**2*an*curr/alamd0

c	peak power angular density according to K.J.Kim
c	GK = aK*(aK**6+24./7.*aK**4+4.*aK**2+16./7.)/
c	1	(1.+aK**2)**3.5
c	dPdOmega = Ptot*21.*gamma**2/(16.*pi*aK)*GK

	write (10,121) alamd0, an, ak, hanuec, alamdc*1.e9, Ptot,
	1	dPdOmega, B_0
121	format (/,' Properties of wiggler',/
	1         ' ---------------------',/
	1 5x,' alamd0 (m)                               ',f10.4,/
	2 5x,' number of periodes                       ',f10.4,/
	3 5x,' wiggler strength                         ',f10.4,//
	4 5x,' critical photon energy (eV)              ',f10.4,/
	5 5x,' critical wavelength (nm)                 ',f10.4,/
	6 5x,' total power (Watt)                       ',f10.3,/
	6 5x,' power density (Watt/rad^2)               ',e10.4,/
	7 5x,' magn. field amplitude (Gauss)            ',f10.4,//
	9 5x,' h*nue (eV)         Brilliance           Flux',/
	1 5x '                     photon             photon',/
	2 5x '          sec*(mm*mrad)**2*0.1% BW   sec*mrad*0.1% BW')

	if (sourceprop.gt.0.) then
	call fileopen
	write (11,*) mesh(icurve)-1
	write (12,*) mesh(icurve)-1
	write (13,*) mesh(icurve)-1
	write (14,*) mesh(icurve)-1
	endif

c	find number of points
	iend = mesh(icurve)
	do 199 imesh=1,mesh(icurve)
	hanue1 = hanuemin
	1 +(hanuemax-hanuemin)*((imesh-1)/float(mesh(icurve)))**2
	y1 = hanue1/hanuec
	if (y1.gt.6. .and. iend.ge.mesh(icurve)) iend = imesh
199	continue	

	mesh(icurve) = iend
	write (6,*) ' In WIGGNEW MESH(icurve):',mesh(icurve)

c	*********************************************************
	do 200 imesh=1,mesh(icurve)		! start main loop
	hanue1 = hanuemin
	1 +(hanuemax-hanuemin)*((imesh-1)/float(mesh(icurve)))**2
	hanue(imesh,icurve) = hanue1
c	wavelength in m
	alamda1 = 1./(hanue1*806500.)
c	this is y as defined in x ray data booklet
	y1 = hanue1/hanuec

c	if (y1.gt.6.) then
c	if (iend.ge.mesh(icurve)) iend = imesh
c		goto 200
c	endif

c	this is H_2(y1) as used by Kwang-Je Kim in Xray data booklet
c	akanue is modified Bessel function of second kind K_2/3
	H2 = (y1*akanue(2./3., y1/2.,ierror))**2

c	This is d2 F / (d theta d psi) of single pole as quoted by
c	Kwang Je Kim which is sometimes called spectral brightness
c	units are phot/(sec mrad**2 0.1BW)
	d2F = 3.461e06 * curr * gamma**2 * H2
c	write (6,*) ' H2(y1=',y1,')=',H2, ' d2F=',d2F

c	here we try to integrate the vertical distribution of SR
c	this is an alternative to using G_0 which requires K_5/3
c
c	we have
c	
c	    d^2F  	           photons                I
c	----------- = 3.461 e12 --------------- gamma^2 (---) *
c	dtheta dpsi             sec mr^2 0.1%BW          Amp
c
c	                                X^2
c	y1^2 (1+X^2)^2 [K_2/3^2(ksi) + ----- K_1/3^2(ksi)]
c	                               1+X^2
c
c	integrate over ten times the diffraction limit
	psimax = 10.*2./(gamma*sqrt(y1))
	npsi = 100
	dpsi = psimax/float(npsi)
	d1F = 0.
c	write (6,*) ' psimax:', psimax,'rad  dpsi:',dpsi

	do ipsi = 1, npsi
c	integrate 0 ... psimax
          psi = dpsi * (float(ipsi) - .5)
c	ksi = y1(1+X^2)^3/2 / 2

c	calculate angular distribution of radiation
c	units are photon/(sec mr^2 Amp 0.1%BW)
          X = gamma * psi
          arg = y1 * (1.+X**2)**(3./2.) / 2.
c	d2F = 3.461 e12 * gamma**2 * curr * y1**2 * (1.+X**2)**2*
c	1 (BSKR3(arg,2)**2 + X**2/(1.+X**2) * BSKR3(arg,1)**2)
          if (arg.gt.40.) exit
          bskr31=bskr3(arg,1)
          bskr32=bskr3(arg,2)
c          print*,"WIGGNEW, BSKR3:",arg,bskr32,bskr31
          dd1F = 3.461 e6 * gamma**2 * curr * (y1 * (1.+X**2))**2*
     &      (BSKR32**2 + X**2/(1.+X**2) * BSKR31**2)
c	write (6,*) ' arg:',arg,' dd1f:',dd1f

c	factor 1000 to convert to 1/mrad
          d1F = d1F + dd1F*dpsi*1000.
        enddo	 	!   vertical integration

c	this is to allow for above and below the ring plane
	d1F = d1F * 2.
	dF = d1F
555	continue

c	write (66,*) alog10(y1), alog10(G1), y1, G1
c	total flux is just flux of one pole * 2 N
	flux (imesh, icurve) = 2.*an*dF

c	get sigrp from on-axis brightness and vertically integrated flux
c	We get 0.64/gamma at hnuec
c	units are rad and m
	sigrp = (1./(1000.*sqrt(2.*pi))) * dF/d2F

c	get sigr from completely coherent radiation requirement
	sigr  = alamda1/(4.*pi*sigrp)	

c	write (6,666) y1, d2F, dF, sigrp, sigr
c666	format (
c	1' y1=',f6.4,' d2F=',e10.4,' dF=',e10.4,
c	2' sigrp(rad)=',e10.4,' sigr(m)=',e10.4,/)

c	do normalisation for brilliance
	call normw (anorm, y1, AK, sigrp, sigr)
c	write (6,*) ' after call to NORMW anorm:',anorm
	brill(imesh,icurve) = d2F * anorm

	if (sourceprop.gt.0.) then
c	this section tries to estimate some effective cross section
c	and divergence as follows

c 	e beam cross section
c	this assumes an e beam waist at the center of the wiggler
c	units are m and rad
	sigx  = sqrt(epsx*betax)
	sigy  = sqrt(epsy*betay)
c	e beam divergence
	sigxp = sqrt(epsx/betax)
	sigyp = sqrt(epsy/betay)
c	e beam excursion
	x0    = AK*alamd0/(gamma*2.*pi)

	sumx = 0.
	sumy = 0.
	n = an
	do i = 1, n
c	distance of source point to reference plane
	zn1 = alamd0*(i-.75) - an*alamd0/2.
	zn2 = alamd0*(i-.25) - an*alamd0/2.
	sumx = sumx +
	1	exp(-.5*x0**2/(sigx**2+sigr**2+(zn1*sigxp)**2))
	2	/(sqrt(sigx**2+sigr**2+(zn1*sigxp)**2))
	sumx = sumx +
	3       exp(-.5*x0**2/(sigx**2+sigr**2+(zn2*sigxp)**2))
	4	/(sqrt(sigx**2+sigr**2+(zn2*sigxp)**2))
	sumy = sumy
	1 + 1./sqrt(sigy**2+sigr**2+(zn1*sigyp)**2)
	2 + 1./sqrt(sigy**2+sigr**2+(zn2*sigyp)**2)
	enddo		! loop over Lichterkette
	sigxeff = 1./sumx
	sigyeff = 1./sumy
	sigxpeff = sqrt(sigxp**2+sigrp**2)
	sigypeff = sqrt(sigyp**2+sigrp**2)
c	2011/01/01
c	naeheres Nachdenken und nervigen Fragen von MSH:
c	diese Addition von sigrp muss in der Lichterkette erfolgen

	write (11,*) hanue1, sigxeff, sigx, sigr, x0
	write (12,*) hanue1, sigxpeff, sigxp, sigrp
	write (13,*) hanue1, sigyeff, sigy, sigr
	write (14,*) hanue1, sigypeff, sigyp, sigrp

	endif		! sourceprop

c	############################
c	coherent flux
c	############################
	if (coherent.gt.0.) then
c	formula due to K.J.Kim, LBL-22236 (march 1987)
c	see comment in undula.for   AG   Nov.30, 1992

c	flux(imesh,icurve)=.385*b1/(hanue(imesh,icurve)**2)
	flux(imesh,icurve)=.385*brill(imesh,icurve)
	1	/(hanue(imesh,icurve)**2)

	endif		! coherent
c	############################

	write (6,220) hanue(imesh,icurve), brill(imesh,icurve),
	1    flux(imesh,icurve)
220	format (' h nue = ',f8.2,' b = ',e10.4,
	1 ' ph/[sec (mm*mrad)**2 0.1%BW]   flux='e10.4)

	write (10,222) hanue(imesh,icurve), brill(imesh,icurve),
	1    flux(imesh,icurve)
222	format (2x,e12.4,10x,e12.6,10x,e12.6)
200	continue		! end main loop over hanue
c	**************************************************

	if (sourceprop.gt.0.) then

	write (11,211) iend-1
211	format (' hanue,  SIGMAXeff, sigmax, sigmar, x0',/
	1' number of points is ',i4)
	write (12,212) iend-1
212	format (' hanue,  SIGMAXPeff, sigmaxp, sigmarp',/
	1' number of points is ',i4)
	write (13,213) iend-1
213	format (' hanue,  SIGMAYeff, sigmay, sigmar',/
	1' number of points is ',i4)
	write (14,214) iend-1
214	format (' hanue,  SIGMAYPeff, sigmayp, sigmarp',/
	1' number of points is ',i4)

	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif		! source properties

	icurve = icurve+1
	ncurve = ncurve+1
	return
	end
*CMZ :  0.01/00 03/09/2013  08.41.47  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.38.28  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.30.20  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
c	************************************************************
c			   F I L E O P E N
c	************************************************************
	subroutine fileopen
c	include 'userdisk_5:[gaupp.brill]brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.
	
c	generates file name on scratch with number ndatfil as extension
c	and opens these files
c	these files are to be used to store the effective source dimensions

c	write (6,99) ndatfil, sourceprop-1
c99	format ('  anfang in  generate_name: ndatfile= ',i5,
c	1 ' sourceprop-1= ',e10.4)
cmsh	encode (3,101,number) ifix(ndatfil+sourceprop-1)
	write(number,fmt=101) ifix(ndatfil+sourceprop-1)
101	format (i3)
c	write (6,102) number
c102	format (' in generate_name: number=',a3)

	if (ifix(ndatfil+sourceprop-1).lt.1000) then
c	write (6,*) ' filenumber .lt.1000'
	filename_xs = 'scratch:size_xt.nnn'
	filename_xd = 'scratch:dive_xt.nnn'
	filename_ys = 'scratch:size_yt.nnn'
	filename_yd = 'scratch:dive_yt.nnn'
c	write (6,105) filename_s, filename_d
	filename_xs(17:19) = number(1:3)
	filename_xd(17:19) = number(1:3)
	filename_ys(17:19) = number(1:3)
	filename_yd(17:19) = number(1:3)
c	write (6,105) filename_s, filename_d
	endif
	if (ifix(ndatfil+sourceprop-1).lt.100) then
c	write (6,*) ' filenumber .lt.100'
	filename_xs = 'scratch:size_xt.nnn'
	filename_xd = 'scratch:dive_xt.nnn'
	filename_ys = 'scratch:size_yt.nnn'
	filename_yd = 'scratch:dive_yt.nnn'
c	write (6,105)
c	1 filename_xs, filename_xd, filename_ys, filename_yd
	filename_xs(18:19) = number(2:3)
	filename_xd(18:19) = number(2:3)
	filename_ys(18:19) = number(2:3)
	filename_yd(18:19) = number(2:3)
c	write (6,105)
c	1 filename_xs, filename_xd, filename_ys, filename_yd
	filename_xs(17:17) = '0'
	filename_xd(17:17) = '0'
	filename_ys(17:17) = '0'
	filename_yd(17:17) = '0'
c	write (6,105)
c	1 filename_xs, filename_xd, filename_ys, filename_yd
	endif
	if (ifix(ndatfil+sourceprop-1).lt.10) then
c	write (6,*) ' filenumber .lt.10'
	filename_xs = 'scratch:size_xt.nnn'
	filename_xd = 'scratch:dive_xt.nnn'
	filename_ys = 'scratch:size_yt.nnn'
	filename_yd = 'scratch:dive_yt.nnn'
c	write (6,105) filename_s, filename_d
	filename_xs(19:19) = number(3:3)
	filename_xd(19:19) = number(3:3)
	filename_ys(19:19) = number(3:3)
	filename_yd(19:19) = number(3:3)
c	write (6,105) filename_s, filename_d
	filename_xs(17:18) = '00'
	filename_xd(17:18) = '00'
	filename_ys(17:18) = '00'
	filename_yd(17:18) = '00'
c	write (6,105) filename_s, filename_d
	endif

	write (6,105)
	1 filename_xs, filename_xd, filename_ys, filename_yd
105	format (' in FILEOPEN:',/
	1' filename_xs = ',a80,/
	2' filename_xd = ',a80,/
	1' filename_ys = ',a80,/
	2' filename_yd = ',a80)
	open (unit=11,file=filename_xs,status='unknown')
	open (unit=12,file=filename_xd,status='unknown')
	open (unit=13,file=filename_ys,status='unknown')
	open (unit=14,file=filename_yd,status='unknown')
 	return
	end
*CMZ :  1.02/00 20/09/2013  12.06.09  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.08.55  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  09.08.32  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.34.31  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	subroutine insert

c	*************************************************************
c				 I N S E R T
c	*************************************************************

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

	goto 2
1	write (6,10)
10	format ('$ length of straight section in m: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alength
        endif
	write (6,20)
20	format (' horizontal and vertical betatron function in center',/
	1 '$ of straight section (m/rad): ')
779       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) betax, betay
        endif
2	write (6,100) alength, betax, betay
100	format(///' properties of insertion',/
	1 ' length of straight section (m)         ALENGTH = ',f10.1,/
	2 ' horizontal betatron function (m/rad)     BETAX = ',f10.1,/
	3 ' vertical   betatron function (m/rad)     BETAY = ',f10.1,//
	9 '$correct? (y/n)    ')
710       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          a=cline(1:1)
        endif
cmsh	read(5,101) a
101	format (a1)
	if (a.eq.'y') goto 90
	goto 1
90	continue
	write (10,91) alength, betax, betay
91	format (/
	1 ' insertion parameters',/
	1 ' --------------------',/
	2 5x,' length of section (m)	',f10.4,/
	3 5x,' betax (m/rad)		',f10.4,/
	4 5x,' betay (m/rad)		',f10.4,/)
	return
	end
*CMZ :  0.01/00 03/09/2013  09.15.01  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.56  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	subroutine normw (anorm, y, AK, sigrp, sigr)
c	*****************************************************************
c				N O R M W
c	*****************************************************************

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

c	units are m and rad
c 	e beam cross section
	sigx  = sqrt(epsx*betax)
	sigy  = sqrt(epsy*betay)
c	e beam divergence
	sigxp = sqrt(epsx/betax)
	sigyp = sqrt(epsy/betay)
c	e beam excursion
	x0    = AK*alamd0/(gamma*2.*pi)
c	write (6,667) sigx,sigxp, sigy, sigyp, x0, sigr, sigrp
c667	format ('  sigx=',e10.4,' sigxp = ',e10.4,/
c	1'  sigy=',e10.4,' sigyp = ',e10.4,
c	2' x0 = ',e10.4,/
c	3' sigr = ',e10.4,'m sigrp = ',e10.4,'rad')

c	reduction of brilliance due to depth of field according to
c	K.J.Kim, LBL 22317 and NIM A261, 44 (1987) (SRI Novosibirsk 1986)
c	Kwang-Je Kim gives the following expression for the brightness of a
c	wiggler.
c	He points out that the term z**2 * sigrp**2 should not be there.
c
c	This is the spatial and angular distribution of brightness
c
c	B(x,y,phi,psi) = (d2 F / dphi dpsi)  1/2pi  *
c	
c	SUM(j) exp{-(1/2) [ (x-xj+phi*zj)**2/(sigx**2+sigr**2+zj**2*sigxp**2)
c		      + (zj*psi+y)**2/Sigma1**2 + psi**2/Sigma2**2
c		      + y**2/Sigma3**2 ] }
c				      /
c	        (sqrt{sigx**2+sigr**2+zj**2*sigxp**2}*Sigma1)
c
c	where
c	Sigma1**2 = (sigy**2+sigr**2)*(1 + sigyp**2/sigrp**2 ) + zj**2*sigyp**2
c
c	Sigma2**2 = sigrp**2 + sigyp**2 + (zj*siyp*sigrp)**2/(sigy**2+sigr**2)
c
c	Sigma3**2 = z**2*sigrp**2 + (sigr**2+sigy**2)*(1+sirp**2/siyp**2)
c
c			Nov. 15, 1993
c	Note: There is an inconsistency for x and y in Sigma1 and Sigma2
c	numerically the difference is insignificant as long as
c	lambda >> emittance.     AG 94/01/05

	n = an
	sum = 0.
	do 100 i = 1,n
	zn1 = alamd0*(i-.75) - an*alamd0/2.
	zn2 = alamd0*(i-.25) - an*alamd0/2.
c	type *, ' zn1:',zn1,'  zn2:',zn2

	sigma1=sqrt(sigy**2+sigr**2+(zn1*sigyp)**2+
	1 (epsy**2+(sigyp*sigr)**2)/sigrp**2)
c	write (6,*) ' sigma1:',sigma1

c	sigma1=sqrt((sigy**2+sigr**2)*(1.+(sigyp/sigrp)**2)+(zn1*sigyp)**2)
c	sigma1=sqrt((sigy**2+sigr**2+(zn1*sigyp)**2)*(1.+(sigyp/sigrp)**2))
	sum = sum +
	1	exp(-.5*x0**2/(sigx**2+sigr**2+(zn1*sigxp)**2))
	2	/(sqrt(sigx**2+sigr**2+(zn1*sigxp)**2)*sigma1)
c	arg = -.5*x0**2/(sigx**2+sigr**2+(zn1*sigxp)**2)
c	write (6,*) ' arg:',arg, '  sum:', sum

	sigma1=sqrt(sigy**2+sigr**2+(zn2*sigyp)**2+
	1 (epsy**2+(sigyp*sigr)**2)/sigrp**2)
c	write (6,*) ' sigma1:',sigma1
c	sigma1=sqrt((sigy**2+sigr**2)*(1.+(sigyp/sigrp)**2)+(zn2*sigyp)**2)
c	sigma1=sqrt((sigy**2+sigr**2+(zn2*sigyp)**2)*(1.+(sigyp/sigrp)**2))
	sum = sum +
	3       exp(-.5*x0**2/(sigx**2+sigr**2+(zn2*sigxp)**2))
	4	/(sqrt(sigx**2+sigr**2+(zn2*sigxp)**2)*sigma1)
c	arg = -.5*x0**2/(sigx**2+sigr**2+(zn2*sigxp)**2)
c	write (6,*) ' arg:',arg, '  sum:', sum

100	continue
c	type*, 'sum/(2.*pi):',sum/(2.*pi)
c	units are mm**2
	anorm = 1.e-6*sum/(2.*pi)
	return
	end
*CMZ :  2.00/00 06/12/2017  17.09.20  by  Michael Scheer
*CMZ :  1.02/01 21/11/2017  14.02.35  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  09.16.27  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.56  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	real function akanue (anue,arg,ierror)
c	****************************************************************
c				A K A N U E
c	****************************************************************
	common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi
	real imjpl, imjmn, imkanue
	real Jp53, Jm53, Ip53, Im53
	dimension B(0:9)

c	calculate modified bessel function of the second kind K_2/3 (arg)
c	and K_5/3(arg)
c	according to Jackson eq. 3.101, p.75, and eq.3.86, p.71

c	argument and function value of K_2/3 are real numbers
c	use real subroutine BSKR3
	if (abs(anue-2./3.) .lt. .01) goto 333
	if (abs(anue-5./3.) .lt. .01) goto 555
	write (6,*) ' Call to AKANUE with anue = ',anue

cc	the following code is not used
cc	However COMBES is a complex function
c
c	ierror = 0
c	if (abs(2.*anue-float(int(2.*anue))).ge.1.e-4) goto 2
c	write (6,1) anue, arg
c1	format (/' **** warning ****   call of modified bessl function with'
c	1 ' half integer order',/
c	2 ' nue = ',e10.4,'    arg = ',e10.4/)
c2	continue
c	if (anue.gt.0.) goto 4
c	write (6,3) anue
c3	format (/' **** warning ****    call to modified bessel function with'
c	1 ' negativ index',/
c	2 ' nue = ',e10.4,/)
c4	continue
cc	need to have two calls to combes for positive and negative index
c	x     = 0.
c	y     = arg
c	alpha = anue - int(anue)
c	beta  = 0.
c	n1    = int(anue)
c	ny    = 0
c	call combes (x, y, alpha, beta, n1, bjre, bjim, yre, yim, ny)
c	rejpl = bjre (n1+1)
c	imjpl = bjim (n1+1)
cc	write (6,51) x, y, alpha, beta, n1, rejpl,imjpl
cc51	format (' call to COMBES with x=',e10.4,' y=',e10.4,
cc	1 ' alpha=',e10.4,' beta=',e10.4,' n1=',i2,
cc	2 ' result: rejpl=',e10.4,' imjpl=',e10.4)
c	n2    = -int(anue)-1
c	alpha = -anue-(int(-anue)-1.)
c	call combes (x, y, alpha, beta, n2, bjre, bjim, yre, yim, ny)
c	rejmn = bjre (-n2+1)
c	imjmn = bjim (-n2+1)
cc	write (6,52) x, y, alpha, beta, n2, rejmn, imjmn
cc52	format (' call to COMBES with x=',e10.4,' y=',e10.4,
cc	1 ' alpha=',e10.4,' beta=',e10.4,' n2=',i2,
cc	2 ' result: rejmn=',e10.4,' imjmn=',e10.4)
cc	write (6,101) anue, x, y, rejpl, imjpl, anue, x, y, rejmn, imjmn
cc101	format (/' bessel function of order  anue',/
cc	1 ' j(+',f10.4,', ',e10.4,' +i*',e10.4,') = ',e10.4,' + i*',e10.4,/
cc	2 ' j(-',f10.4,', ',e10.4,' +i*',e10.4,') = ',e10.4,' + i*',e10.4,/)
c	cospi  = cos (pi*anue)
c	sinpi  = sin (pi*anue)
c	cospi2 = cos ((anue+1.)*pi/2.)
c	sinpi2 = sin ((anue+1.)*pi/2.)
c	expr1  = rejpl-(imjpl*cospi-imjmn)/sinpi
c	expr2  = imjpl+(rejpl*cospi-rejmn)/sinpi
c	rekanue = (pi/2.)*(cospi2*expr1 - sinpi2*expr2)
c	imkanue = (pi/2.)*(sinpi2*expr1 + cospi2*expr2)
cc	if (abs(imkanue/rekanue).le.3.e-4) goto 99
cc	
cc	write (6,102) anue,arg,rekanue, imkanue
cc102	format (' **** warning ****',/
cc	1 ' real (ka(nue=',e10.4,',arg=',e10.4,')) = ',e10.4,/
cc	2 ' imag (ka(nue',11x,',arg',11x,')) = ',e10.4,/)
c	if (abs(imkanue/rekanue).gt. .01) ierror = 1
c99	akanue = rekanue
c	return

c	###################################################
333	continue
c	nue = 2/3
c	On Oct 26, 1995 I had the problem that GENLIB was not available
c	due to hard ware problems. Now I try to call
c	BSKR3(x,2) (CERNlib C340) which is simply K_2/3
	if (abs(anue-2./3.) .gt. .001)
	1 write (6,*) ' nue =',anue,' .e. 2/3 in AKANUE'
	nue = int(anue*3. + .0001)
c        print*,"AKANUE, BSKR3:",arg,nue,bskr3(arg,nue)
	akanue = BSKR3 (arg,nue)

c	write (6,334) anue, arg, akanue
c334	format (' call to SR EFKA with anue=',e10.4,' arg=',e10.4,
c	1' akanue=',e10.4)
	return

555	continue
c	nue = 5/3. This part is only used in Function G0
c	Bessel functions of order 5/3 is only included in CERN MATHLIB
c	 must use formula by jackson to derive K_5/3 from J_5/3 and J_-5/3
	if (abs(anue - 5./3.) .gt. .001)
	1 write (6,*) ' nue = ',anue,' .ne. 5/3 in AKANUE'

c	2003/07/09
c	JB finds that akanue(5/3,arg) is an oscillating function of arg.
c	this is not to be expected!
c	The conversion from the Bessel function with an imaginary
c	argument requires the calculation of the Bessel function I(5/3)
c	as in Jackson formula 3.100, p.75. I.e. the call to the CERN
c	routine should be BSIA (instead of BSJA).
c	Further more the argument should not be doubled.
c	There is no problem for BRILL since there is no  call to akanue
c	with nue=5/3. The calculation is done by numerically integrating
c	the vertical distribution.

c	I change the code but keep the numerical integration for flux

c	call BSJA(2.*arg, anue-1., 1, 5, B)
c	this is J_(+5/3) (2 arg)
c	Jp53 = B(1)
c	call BSJA(2.*arg, 1./3., -2, 5, B)
c	this is J_(-5/3) (2 arg)
c	Jm53 = B(2)

	call BSIA(arg, anue-1., 1, 5, B)
c	this is I_(+5/3) (arg)
	Ip53 = B(1)
	call BSIA(arg, 1./3., -2, 5, B)
c	this is I_(-5/3) (arg)
	Im53 = B(2)

c	2003/07/09: this expression is wrong!
c	akanue = pi/2. *
c	1 (cos(4.*pi/3.)*Jp53 - sin(4.*pi/3.)/sin(5.*pi/3.)*
c	2	(Jp53*cos(5.*pi/3.) - Jm53) )
c	the calculation in the book as of 2003/07/09 uses:
c	XRAY DATA BOOKLET K.J.Kim eq. 5
c	Jackson eq. 3.101
c	Jackson eq. 3.86
c	Jackson eq. 3.85
c	Jackson eq. 3.100

c	we find
	akanue = pi/Sqrt(3.) * (Ip53 - Im53)

c	this expression needs to be integrated as in XRAY DATA BOOKLET
c	eq. 5 to yield the universal function G_1(y):
c	G_1(y) = y * Integral(_y, ^Infinity) K_5/3(x) dx
c	not yet verified!	

	write (6,*) ' Call to BSIA with arg=',arg,' AKANUE=',akanue

	return
	end
*CMZ :  0.01/00 03/09/2013  09.15.28  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.56  by  Andreas Gaupp
*-- Author :    Andreas Gaupp
	real function G0(eta)
c	************************************************************
c				G 0
c	************************************************************
c	should never be called again!     Dec 1995
c	calculates integral of vertically integrated flux formula
c	eta is lower limit of integral
c	upper limit is at infinity
c	this is the function defined by Green as quoted by H.Winick

	anue    = 5./3.
	sum     = 0.
c	initial step size
	du0     = eta/10.
	if (du0.le.0.) then
		write (6,7) du0
7		format (' initial step size in G0=',e10.4)
		stop
		endif
	dsum0   = akanue(anue,eta,ierror)*du0
	u       = eta
	du      = du0
	sum     = sum + dsum0

10	u    =  u + du
	dsum = akanue(anue,u,ierror)*du
	sum  = sum + dsum
	if (abs(dsum).lt.abs(dsum0)*.01) goto 90
	if (ierror.ge.1)                 goto 90
	du   = du*1.1
	goto 10
90	continue

c	write (6,91) eta, dsum0, du, dsum, sum
c91	format (' eta=',e10.4,' dsum0=',e10.4,' du=',
c	1	e10.4,' dsum=',e10.4,' G0=',e10.4)
	G0   = sum
	return
	end
*CMZ :  2.00/00 06/12/2017  17.02.44  by  Michael Scheer
*CMZ :  1.02/01 21/11/2017  14.02.35  by  Michael Scheer
*CMZ :  1.02/00 08/10/2013  11.20.58  by  Michael Scheer
*CMZ :  1.01/00 11/09/2013  13.18.56  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.09.35  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.56  by  Andreas Gaupp
*-- Author : Michael Scheer
	subroutine wiggler_msh

c Wie wiggnew, aber ohne die exponentiellen Abhängigkeiten, dafür mit Faktor
c 4pi statt 2pi. Siehe auch brill-wiggler-dipole.pdf
c Die alte Formel ist unklar, siehe Brill-Ordner. Sie ist in
c brill-wiggler.kumac als Modus implementiert, stimmt aber nicht mit den
c Resulaten dieser Brill-Version überein!? (Sep 2013)


c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

1	write (6,20)
20	format ('$ enter length of periode (m) ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alamd0
        endif
	write (6,10)
10	format ('$ enter number of periodes: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) an
        endif
	if(an*alamd0.gt.alength) write (6,21) an,alamd0,alength
21	format (' **** warning ****',/
	1 ' number of periodes                AN = ',f10.3,/
	2 ' length of periode (m)        alamd0  = ',f10.3,/
	3 ' length of long straight (m)  alength = ',f10.3,/
	4 ' *** inconsisten ***',/)
	write (6,30)
30	format ('$ enter wiggler strength AK:  ')
710     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) AK
        endif
				
2	continue
c	characteristic wavelength in m
c	checked on April 4, 1985
c	numerically the same as  hanuec=1.73e-8(eV/Gauss)*gamma**2/B
	alamdc = (2./3.)*alamd0/(gamma**2*ak)
	hanuec = 1./(alamdc*806500.)
	B_0  = AK*2.*pi*.911e-27*9.e20/(4.8e-10*alamd0*100.)

	write (6,100)  alamd0, an, ak, hanuec, B_0
100	format (/,' wiggler parameter ',/
	2 ' periode length (m)          alamd0 = ',f10.3,/
	1 ' number of periodes              an = ',f10.3,/
	3 ' wiggler strength                AK = ',f10.3,/
	4 ' critical photon energy (eV) hanuec = ',f10.3,/
	4 ' Magnetic field amplitude (Gauss)   = ',f10.3,//
	5 '$correct ? (y/n)  ')
711       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          a=cline(1:1)
        endif
cmsh	read(5,101)  a
101	format (a1)
	if (a.ne.'y' .and. a.ne.'Y') goto 1

c	calculate wiggler spectrum
	write (6,110)
110	format ('Enter photon energy range (in eV), hanuemin, hanuemax:  ')
712     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 712
        else
          read(cline,*) hanuemin, hanuemax
        endif
	write (6,120)
120	format ('$ enter number of points in this range: ')
713     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256

c	total power in Watt according to K.J.Kim, NIM A246, 67 (1986)
	call power (aN, aK, alamd0, gamma, curr, Ptot, dPdOmega)
c	Ptot = (an/6.)*377.*curr*(2.*pi*1.6e-19*3e8)/alamd0
c	1	*gamma**2*aK**2

c	Ptot = 1.774e-8*(gamma*ak)**2*an*curr/alamd0

c	peak power angular density according to K.J.Kim
c	GK = aK*(aK**6+24./7.*aK**4+4.*aK**2+16./7.)/
c	1	(1.+aK**2)**3.5
c	dPdOmega = Ptot*21.*gamma**2/(16.*pi*aK)*GK

	write (10,121) alamd0, an, ak, hanuec, alamdc*1.e9, Ptot,
	1	dPdOmega, B_0
121	format (/,' Properties of wiggler',/
	1         ' ---------------------',/
	1 5x,' alamd0 (m)                               ',f10.4,/
	2 5x,' number of periodes                       ',f10.4,/
	3 5x,' wiggler strength                         ',f10.4,//
	4 5x,' critical photon energy (eV)              ',f10.4,/
	5 5x,' critical wavelength (nm)                 ',f10.4,/
	6 5x,' total power (Watt)                       ',f10.3,/
	6 5x,' power density (Watt/rad^2)               ',e10.4,/
	7 5x,' magn. field amplitude (Gauss)            ',f10.4,//
	9 5x,' h*nue (eV)         Brilliance           Flux',/
	1 5x '                     photon             photon',/
	2 5x '          sec*(mm*mrad)**2*0.1% BW   sec*mrad*0.1% BW')

	if (sourceprop.gt.0.) then
	call fileopen
	write (11,*) mesh(icurve)-1
	write (12,*) mesh(icurve)-1
	write (13,*) mesh(icurve)-1
	write (14,*) mesh(icurve)-1
	endif

c	find number of points
	iend = mesh(icurve)
	do 199 imesh=1,mesh(icurve)
	hanue1 = hanuemin
	1 +(hanuemax-hanuemin)*((imesh-1)/float(mesh(icurve)))**2
	y1 = hanue1/hanuec
	if (y1.gt.6. .and. iend.ge.mesh(icurve)) iend = imesh
199	continue	

	mesh(icurve) = iend
	write (6,*) ' In WIGGNEW MESH(icurve):',mesh(icurve)

c	*********************************************************
	do imesh=1,mesh(icurve)		! start main loop
	hanue1 = hanuemin
	1 +(hanuemax-hanuemin)*((imesh-1)/float(mesh(icurve)))**2
	hanue(imesh,icurve) = hanue1
c	wavelength in m
	alamda1 = 1./(hanue1*806500.)
c	this is y as defined in x ray data booklet
	y1 = hanue1/hanuec

c	if (y1.gt.6.) then
c	if (iend.ge.mesh(icurve)) iend = imesh
c		goto 200
c	endif

c	this is H_2(y1) as used by Kwang-Je Kim in Xray data booklet
c	akanue is modified Bessel function of second kind K_2/3
	H2 = (y1*akanue(2./3., y1/2.,ierror))**2

c	This is d2 F / (d theta d psi) of single pole as quoted by
c	Kwang Je Kim which is sometimes called spectral brightness
c	units are phot/(sec mrad**2 0.1BW)
	d2F = 3.461e06 * curr * gamma**2 * H2
c	write (6,*) ' H2(y1=',y1,')=',H2, ' d2F=',d2F

c	here we try to integrate the vertical distribution of SR
c	this is an alternative to using G_0 which requires K_5/3
c
c	we have
c	
c	    d^2F  	           photons                I
c	----------- = 3.461 e12 --------------- gamma^2 (---) *
c	dtheta dpsi             sec mr^2 0.1%BW          Amp
c
c	                                X^2
c	y1^2 (1+X^2)^2 [K_2/3^2(ksi) + ----- K_1/3^2(ksi)]
c	                               1+X^2
c
c	integrate over ten times the diffraction limit
	psimax = 10.*2./(gamma*sqrt(y1))
	npsi = 100
	dpsi = psimax/float(npsi)
	d1F = 0.
c	write (6,*) ' psimax:', psimax,'rad  dpsi:',dpsi

	do ipsi = 1, npsi
c	integrate 0 ... psimax
          psi = dpsi * (float(ipsi) - .5)
c	ksi = y1(1+X^2)^3/2 / 2

c	calculate angular distribution of radiation
c	units are photon/(sec mr^2 Amp 0.1%BW)
          X = gamma * psi
          arg = y1 * (1.+X**2)**(3./2.) / 2.
c	d2F = 3.461 e12 * gamma**2 * curr * y1**2 * (1.+X**2)**2*
c	1 (BSKR3(arg,2)**2 + X**2/(1.+X**2) * BSKR3(arg,1)**2)
          if (arg.gt.40.) exit
          bskr31=bskr3(arg,1)
          bskr32=bskr3(arg,2)
c          print*,"WIGGLER_MSH, BSKR3:",arg,bskr32,bskr31
          if (bskr31.lt.1.0e-20) bskr31=0.0
          if (bskr32.lt.1.0e-20) bskr32=0.0
          dd1F = 3.461 e6 * gamma**2 * curr * (y1 * (1.+X**2))**2*
     &      (BSKR32**2 + X**2/(1.+X**2) * BSKR31**2)
c	write (6,*) ' arg:',arg,' dd1f:',dd1f

c	factor 1000 to convert to 1/mrad
          d1F = d1F + dd1F*dpsi*1000.

c        print*,"ipsi:",ipsi
        enddo	 	!   vertical integration

c	this is to allow for above and below the ring plane
	d1F = d1F * 2.
	dF = d1F
555	continue

c	write (66,*) alog10(y1), alog10(G1), y1, G1
c	total flux is just flux of one pole * 2 N
	flux (imesh, icurve) = 2.*an*dF

c	get sigrp from on-axis brightness and vertically integrated flux
c	We get 0.64/gamma at hnuec
c	units are rad and m
	sigrp = (1./(1000.*sqrt(2.*pi))) * dF/d2F

c	get sigr from completely coherent radiation requirement
	sigr  = alamda1/(4.*pi*sigrp)	

c	write (6,666) y1, d2F, dF, sigrp, sigr
c666	format (
c	1' y1=',f6.4,' d2F=',e10.4,' dF=',e10.4,
c	2' sigrp(rad)=',e10.4,' sigr(m)=',e10.4,/)

c	do normalisation for brilliance
	call normw_msh (anorm, y1, AK, sigrp, sigr)
c	write (6,*) ' after call to NORMW anorm:',anorm
        brill(imesh,icurve) = d2F * anorm
        fluxdens(imesh,icurve)=d2F*2.0*an !msh 8.10.2013

	if (sourceprop.gt.0.) then
c	this section tries to estimate some effective cross section
c	and divergence as follows

c 	e beam cross section
c	this assumes an e beam waist at the center of the wiggler
c	units are m and rad
	sigx  = sqrt(epsx*betax)
	sigy  = sqrt(epsy*betay)
c	e beam divergence
	sigxp = sqrt(epsx/betax)
	sigyp = sqrt(epsy/betay)
c	e beam excursion
	x0    = AK*alamd0/(gamma*2.*pi)

	sumx = 0.
	sumy = 0.
	n = an
	do i = 1, n
c	distance of source point to reference plane
	zn1 = alamd0*(i-.75) - an*alamd0/2.
	zn2 = alamd0*(i-.25) - an*alamd0/2.
	sumx = sumx +
cmsh 	1	exp(-.5*x0**2/(sigx**2+sigr**2+(zn1*sigxp)**2))
     &    0.5 ! nur eine Lichterkette eingeziehen.
	2	/(sqrt(sigx**2+sigr**2+(zn1*sigxp)**2))
	sumx = sumx +
cmsh	3       exp(-.5*x0**2/(sigx**2+sigr**2+(zn2*sigxp)**2))
     &    0.5 ! nur eine Lichterkette eingeziehen.
	4	/(sqrt(sigx**2+sigr**2+(zn2*sigxp)**2))
	sumy = sumy
	1 + 1./sqrt(sigy**2+sigr**2+(zn1*sigyp)**2)
	2 + 1./sqrt(sigy**2+sigr**2+(zn2*sigyp)**2)
	enddo		! loop over Lichterkette
	sigxeff = 1./sumx
	sigyeff = 1./sumy
	sigxpeff = sqrt(sigxp**2+sigrp**2)
	sigypeff = sqrt(sigyp**2+sigrp**2)
c	2011/01/01
c	naeheres Nachdenken und nervigen Fragen von MSH:
c	diese Addition von sigrp muss in der Lichterkette erfolgen

	write (11,*) hanue1, sigxeff, sigx, sigr, x0
	write (12,*) hanue1, sigxpeff, sigxp, sigrp
	write (13,*) hanue1, sigyeff, sigy, sigr
	write (14,*) hanue1, sigypeff, sigyp, sigrp

	endif		! sourceprop

c	############################
c	coherent flux
c	############################
	if (coherent.gt.0.) then
c	formula due to K.J.Kim, LBL-22236 (march 1987)
c	see comment in undula.for   AG   Nov.30, 1992

c	flux(imesh,icurve)=.385*b1/(hanue(imesh,icurve)**2)
	flux(imesh,icurve)=.385*brill(imesh,icurve)
	1	/(hanue(imesh,icurve)**2)

	endif		! coherent
c	############################

	write (6,220) hanue(imesh,icurve), brill(imesh,icurve),
	1    flux(imesh,icurve)
220	format (' h nue = ',f8.2,' b = ',e10.4,
	1 ' ph/[sec (mm*mrad)**2 0.1%BW]   flux='e10.4)

	write (10,222) hanue(imesh,icurve), brill(imesh,icurve),
	1    flux(imesh,icurve)
222	format (2x,e12.4,10x,e12.6,10x,e12.6)
200     continue		! end main loop over hanue
        enddo
c	**************************************************

	if (sourceprop.gt.0.) then

	write (11,211) iend-1
211	format (' hanue,  SIGMAXeff, sigmax, sigmar, x0',/
	1' number of points is ',i4)
	write (12,212) iend-1
212	format (' hanue,  SIGMAXPeff, sigmaxp, sigmarp',/
	1' number of points is ',i4)
	write (13,213) iend-1
213	format (' hanue,  SIGMAYeff, sigmay, sigmar',/
	1' number of points is ',i4)
	write (14,214) iend-1
214	format (' hanue,  SIGMAYPeff, sigmayp, sigmarp',/
	1' number of points is ',i4)

	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif		! source properties

	icurve = icurve+1
	ncurve = ncurve+1
	return
	end
*CMZ :  1.01/00 11/09/2013  13.29.09  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  09.15.01  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.56  by  Andreas Gaupp
*-- Author :    Michael Scheer
	subroutine normw_msh (anorm, y, AK, sigrp, sigr)
c	*****************************************************************
c				N O R M W
c	*****************************************************************

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

c	units are m and rad
c 	e beam cross section
	sigx  = sqrt(epsx*betax)
	sigy  = sqrt(epsy*betay)
c	e beam divergence
	sigxp = sqrt(epsx/betax)
	sigyp = sqrt(epsy/betay)
c	e beam excursion
	x0    = AK*alamd0/(gamma*2.*pi)
c	write (6,667) sigx,sigxp, sigy, sigyp, x0, sigr, sigrp
c667	format ('  sigx=',e10.4,' sigxp = ',e10.4,/
c	1'  sigy=',e10.4,' sigyp = ',e10.4,
c	2' x0 = ',e10.4,/
c	3' sigr = ',e10.4,'m sigrp = ',e10.4,'rad')

c	reduction of brilliance due to depth of field according to
c	K.J.Kim, LBL 22317 and NIM A261, 44 (1987) (SRI Novosibirsk 1986)
c	Kwang-Je Kim gives the following expression for the brightness of a
c	wiggler.
c	He points out that the term z**2 * sigrp**2 should not be there.
c
c	This is the spatial and angular distribution of brightness
c
c	B(x,y,phi,psi) = (d2 F / dphi dpsi)  1/2pi  *
c	
c	SUM(j) exp{-(1/2) [ (x-xj+phi*zj)**2/(sigx**2+sigr**2+zj**2*sigxp**2)
c		      + (zj*psi+y)**2/Sigma1**2 + psi**2/Sigma2**2
c		      + y**2/Sigma3**2 ] }
c				      /
c	        (sqrt{sigx**2+sigr**2+zj**2*sigxp**2}*Sigma1)
c
c	where
c	Sigma1**2 = (sigy**2+sigr**2)*(1 + sigyp**2/sigrp**2 ) + zj**2*sigyp**2
c
c	Sigma2**2 = sigrp**2 + sigyp**2 + (zj*siyp*sigrp)**2/(sigy**2+sigr**2)
c
c	Sigma3**2 = z**2*sigrp**2 + (sigr**2+sigy**2)*(1+sirp**2/siyp**2)
c
c			Nov. 15, 1993
c	Note: There is an inconsistency for x and y in Sigma1 and Sigma2
c	numerically the difference is insignificant as long as
c	lambda >> emittance.     AG 94/01/05

	n = an
	sum = 0.
	do 100 i = 1,n
	zn1 = alamd0*(i-.75) - an*alamd0/2.
	zn2 = alamd0*(i-.25) - an*alamd0/2.
c	type *, ' zn1:',zn1,'  zn2:',zn2

	sigma1=sqrt(sigy**2+sigr**2+(zn1*sigyp)**2+
	1 (epsy**2+(sigyp*sigr)**2)/sigrp**2)
c	write (6,*) ' sigma1:',sigma1

c	sigma1=sqrt((sigy**2+sigr**2)*(1.+(sigyp/sigrp)**2)+(zn1*sigyp)**2)
c	sigma1=sqrt((sigy**2+sigr**2+(zn1*sigyp)**2)*(1.+(sigyp/sigrp)**2))
	sum = sum +
cmsh	1	exp(-.5*x0**2/(sigx**2+sigr**2+(zn1*sigxp)**2))
     &    0.5 ! nur eine Lichterkette
	2	/(sqrt(sigx**2+sigr**2+(zn1*sigxp)**2)*sigma1)
c	arg = -.5*x0**2/(sigx**2+sigr**2+(zn1*sigxp)**2)
c	write (6,*) ' arg:',arg, '  sum:', sum

	sigma1=sqrt(sigy**2+sigr**2+(zn2*sigyp)**2+
	1 (epsy**2+(sigyp*sigr)**2)/sigrp**2)
c	write (6,*) ' sigma1:',sigma1
c	sigma1=sqrt((sigy**2+sigr**2)*(1.+(sigyp/sigrp)**2)+(zn2*sigyp)**2)
c	sigma1=sqrt((sigy**2+sigr**2+(zn2*sigyp)**2)*(1.+(sigyp/sigrp)**2))
	sum = sum +
cmsh	3       exp(-.5*x0**2/(sigx**2+sigr**2+(zn2*sigxp)**2))
     &    0.5 ! nur eine Lichterkette
	4	/(sqrt(sigx**2+sigr**2+(zn2*sigxp)**2)*sigma1)
c	arg = -.5*x0**2/(sigx**2+sigr**2+(zn2*sigxp)**2)
c	write (6,*) ' arg:',arg, '  sum:', sum

100	continue
c	type*, 'sum/(2.*pi):',sum/(2.*pi)
c	units are mm**2
	anorm = 1.e-6*sum/(2.*pi)
	return
	end
*CMZ :  1.02/01 16/07/2015  09.19.08  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.09  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.14.24  by  Michael Scheer
*CMZ :  1.01/00 12/09/2013  15.01.41  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.10.31  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.17  by  Andreas Gaupp
*-- Author :    Michael Scheer
c	****************************************************************
c	Based on undunew of Andreas Gaupp, see comments there
c	****************************************************************
	subroutine undula_walker

c Stand 8.3.2015

c	sigmar  = sqrt(2*lamda*L)/(2*pi)
c	sigmarp = sqrt(lamda/2.*L)
c	note: sigmar*sigmarp = lamda/(2pi)

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

1	write (6,20)
20	format ('$ enter periode length in m: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alamd0
        endif
	write (6,10)
10	format ('$ enter number of undulator periodes: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) an
        endif
	if ( an*alamd0.ge.alength) write (6,21) alamd0, an, alength
21	format (' **** warning ****',/
	1 ' lamda0  = ',f10.2,' m',/
	2 '      N  = ',f10.0,/
	3 ' alength = ',f10.2,' m '/
	4 ' *** not compatible ***'/)
	
	write (6,30)
30	format('$ enter range in wiggler strength, KMIN and KMAX: ')
710     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) akmin, akmax
        endif
	if (akmin.ge.0. .and. akmax.ge.akmin) goto 32
	write (6,31) akmin, akmax
31	format (' **** warning ****',/
	1 ' AKMIN   = ',f10.3,'   AKMAX  = ',f10.3,' *** inconsistent'/)
32	write (6,40)
40	format ('$ enter order of undulator spectrum (1, 3, 5,..): ')
711     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          read(cline,*) smallk
        endif
	if (smallk.lt.1. .or.
	1   abs(smallk-float(int(smallk))).gt.1.e-4) write (6,41) smallk
41	format (' **** warning ****',/
	1 ' bad value for order of spectrum. SMALLK = ',f10.3)
2	continue

c	get minimum gap for hybrid-REC structure using Halbach's formula
c	for REC without 10% reduction
	gapmin = alamd0*(1.519-sqrt(2.309+.5556*alog(.003223*akmax/alamd0)))	
	write (6,100) an, alamd0, akmin, akmax, smallk, gapmin
100	format (/,' undulator parameters',/
	1 ' number of periodes             AN     = ',f10.0,/
	2 ' periode length (m)             alamd0 = ',f10.5,/
	3 ' minimum K                      KMIN   = ',f10.2,/
	4 ' maximum K                      KMAX   = ',f10.2,/
	5 ' order of spectrum              smallk = ',f10.0,/
	6 ' minimum gap (hybrid REC) (m)   gapmin = ',f10.5,//
	9 '$ correct? (y/n)  ')
712       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 712
        else
          a=cline(1:1)
        endif
cmsh	read(5,101) a
101	format (a1)
	if (a.eq.'y') goto 109
	goto 1
109	continue

c	plot undulator spectrum
	write (6,110)
110	format (' enter photon energy range (eV) of this curve',/
	1 '$ hanuemin and hanuemax:  ')
713       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          read(cline,*) hanuemin, hanuemax
        endif
	if (hanuemax.le.0.) then
		hanuemax = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*(1.+0.5*aKmin**2)))
	endif
	if (hanuemin.le.0.) then
		hanuemin = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*(1.+0.5*aKmax**2)))
	endif
	write (6,120)
120	format ('Enter number of points on this curve, mesh:  ')
714     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 714
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256

c	calculate normalisation of brilliance for K = 0
c	only for initial printout, not for final curve
	theta0 = 1./(gamma*sqrt(an*smallk))
c	note alength is in m
	sigmax = sqrt(epsx*betax)
  	sigmay = sqrt(epsy*betay)
	sigmaxp= sqrt(theta0**2+epsx/betax)
	sigmayp= sqrt(theta0**2+epsy/betay)
	write (6,121) theta0, sigmax, sigmay, sigmaxp, sigmayp, dgamma
121	format (/
	1 ' theta0 (rad, K=0) = ',e10.4,/
	2 ' sigmax (m)        = ',e10.4,/
	3 ' sigmay (m)        = ',e10.4,/
	4 ' sigmaxp (rad) K=0 = ',e10.4,/
	5 ' sigmayp (rad) K=0 = ',e10.4,/
	6 ' sigma E           = ',e10.4,//)
	
c	get maximum power according to K.J.Kim
	call power (an, aKmax, alamd0, gamma, curr, Ptot, d2P)

	write (10,122) alamd0, an, akmin, akmax, smallk, gapmin,
	1              theta0,sigmax, sigmay,sigmaxp,sigmayp,dgamma
	write (10,123) Ptot, d2P

122	format (/
	1 ' properties of undulator',/
	1 ' -----------------------',/
	2 5x,' periode length (m)             	',f10.4,/
	3 5x,' number of periodes              	',f10.4,/
	4 5x,' minimum K                       	',f10.4,/
	5 5x,' maximum K                        ',f10.4,/
	6 5x,' order of spectrum k              ',f10.4,/
	7 5x,' minimum gap (REC hybrid, m)      ',f10.4,//
	1 ' properties of radiation',/
	2 5x,' theta0 (rad) K = 0               ',e10.4,/
	3 5x,' sigmax (m)                       ',e10.4,/
	4 5x,' sigmay (m)                       ',e10.4,/
	5 5x,' sigmaxp (rad) K = 0              ',e10.4,/
	6 5x,' sigmayp (rad) K = 0              ',e10.4,/
	7 5x,' dgamma                           ',e10.4)
123	format (
	7 5x,' Ptot (Watt)                      ',e10.4,/
	7 5x,' dP/dOmega (Watt/mrad^2)          ',e10.4,/
	9 5x,' h*nue (eV)    Brilliance         Flux   Flux density    K
	1       Ptot    gap',/
	2 5x '                  photon         photon      photon
	3     (Watt)     (m)',/
	4 5x '     sec*(mm*mrad)**2*0.1% BW  sec*0.1%BW  sec*0.1%BW*mrad^2'
	5 ,/)

	if (sourceprop.gt.0.) then
	call fileopen
	write (11,*) mesh(icurve)-1
	write (12,*) mesh(icurve)-1
	write (13,*) mesh(icurve)-1
	write (14,*) mesh(icurve)-1

c	e beam cross section and divergence
	sigmax = sqrt(epsx*betax)
	sigmay = sqrt(epsy*betay)
	sigxpr = sqrt(epsx/betax)
	sigypr = sqrt(epsy/betay)

	endif

	do 200 imesh=1,mesh(icurve)
	hanue1 = hanuemin+(hanuemax-hanuemin)*(imesh-1)/float(mesh(icurve)-1)
c	8065 cm**-1  = 1 eV
	ak1 = 2.*(smallk*2.*gamma**2/(hanue1*alamd0*806500.)-1.)
	if (ak1.le.0.) goto 200
	ak1 = sqrt(ak1)
C	this ak1 is wiggler constant from
c	lamda_n = lamda_0/(2.*gamma**2)*(1+(1./2.)*K**2)
c	i.e. planar undulator


	if (ak1.lt.akmin .or. ak1.gt.akmax) goto 200

c	formula due to Krinsky in Koch's book, eq.281, p.152
c	checked on April 9, 1985
c	units are  photons/[sec*mrad**2*0.1%BW*Amp])
c	our efka is the same as F_n(K) in xray data booklet
c	this is the same as xray data booklet eq. 13
c	b1 is on axis spectral brightness
c	units: curr in Amp
c	       b1 in photons/(sec mrad^2 0.1% )

c	b1 = 4.56E7*(an*gamma)**2*curr*efka(smallk,ak1)
	b1= 4.56E7*(an*gamma)**2*curr*efka(ak1)

c	this is the spectral brightness, which is reduced by a
c	finite energy spread
c	 here we take the spectral distribution of undulator radiation
c	to be a normal with fractional width sigma = 1./(smallk*2.6*N)
c	to be a normal with fractional width sigma = 1./(smallk*2.0d0*sqrt(2.0d0)*N) !msh
c	(not 1/(2*N) !!!!!!!!!!!!!!

c	the distribution of ebeam energy is taken to be a normal
c	distribution with fractional width dgamma = Delta gamma / gamma
c	The energy spread of the e beam is multiplied by 2 because
c	gamma**2 enters the energy formula
c	This is necessary as discovered by Howard Padmore at ALS in
c	winter 1993/94.                         AG 1994/01/16

cmsh	b1= b1/sqrt(1.+(2.*dgamma*smallk*2.6*aN)**2)
	b1= b1/sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2)
        dgamred=1./sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2) !msh

c	save on axis spectral brightness into array fluxdens
c	need to correct for electron beam divergence
c	e beam divergence

	sigxpr = sqrt(epsx/betax)
	sigypr = sqrt(epsy/betay)
c	natural divergence of radiation
	alamda = (alamd0/(smallk*2.*gamma**2))*(1.+ak1**2/2.)
	aL = alamd0*an

        sigmarp = sqrt(alamda/(2.0d0*aL))

        print*,'*** Fluxdens unclear in undula_walker! ***'
	fluxdens(imesh,icurve) = b1/
     &    sqrt(1.+(sigxpr**2+sigypr**2)/sigmarp**2)


c	J.B. 30.5.89
c	this is sigmarp = sqrt (n lamda / L)
        sigmarp = 1./gamma*sqrt((1.+ak1**2/2.)/(smallk*2.*an)) / sqrt(2.0d0)

	flux(imesh,icurve) = 2.*pi*b1*(sigmarp*1000.)**2
	hanue (imesh, icurve) = hanue1

        call normu_walker (anorm, ak1)

        if (imesh.gt.256 .or. icurve.gt.100 .or. anorm.le.0.)
	1	write (6,*) ' *** imesh=',imesh,'   icurve=',icurve,
	2                   ' *** anorm=',anorm
	b1= flux(imesh,icurve)/anorm
	brill (imesh,icurve) = b1

        if (sourceprop.gt.0.) then
c	natural divergence of radiation
          alamda = (alamd0/(smallk*2.*gamma**2))*(1.+ak1**2/2.)
          aL = alamd0*an
          sigmarp=sqrt(alamda/(2.0d0*aL))
c	radiation cross section
          sigmar=sqrt(2.0d0*alamda*al)/(2.0d0*pi)
c	MSH behauptet basierend auf einem Vergleich der zweiten Momente
c	(2000/06):
c	sigmar=(1/2 pi) sqrt(2 lambda L)
c	sigmar * sigmarp = (1/2 pi) lambda (*0.82)

          size_xt = sqrt(sigmax**2+sigmar**2)
          dive_xt = sqrt(sigxpr**2 + sigmarp**2)
          size_yt = sqrt(sigmay**2+sigmar**2)
          dive_yt = sqrt(sigypr**2 + sigmarp**2)

          write (11,*) hanue1, size_xt, sigmax, sigmar
          write (12,*) hanue1, dive_xt, sigxpr, sigmarp
          write (13,*) hanue1, size_yt, sigmay, sigmar
          write (14,*) hanue1, dive_yt, sigypr, sigmarp
        endif



c	#########################
c	coherent flux
c	#########################
	if (coherent.gt.0.) then
c	use formula due to K.J.Kim, LBL-22236 (march 1987)
c	published in PAC1987, p 194
c	converting into units photons/sec/0.1%BW I get
c	coherent Flux = .385 Photons/sec/0.1%BW
c	                * B(0,0)/(phot/sec/mm**2 mrad**2 0.1%BW)
c	                * (1/(hnue/eV)**2)
				
c	attention!!!!!!!!!!!!!!!!!!!!!!!!!1
	flux(imesh,icurve) = .385*b1/(hanue1**2)
	endif

c Die Eneriebreite des Strahles hat kaum Einfluss auf den Fluss, deshalb
        flux(imesh,icurve)=flux(imesh,icurve)/dgamred

c	write (6,201)
c	1 hanue(imesh,icurve), ak1, smallk, brill(imesh,icurve)
c201	format (' hanue/eV =',f9.3,' K=',f9.3,' k=',f2.0,' B=',e10.4,
c	1 ' ph/(sec (mm*mrad)**2 .1%)')
c
c	total power in Watt according to K.J.Kim
	ptot = 1.774e-8*(gamma*ak1)**2*an*curr/alamd0

	dummy = 2.309+.5556*alog(.003223*ak1/alamd0)
	gap = 0.
	if (dummy.gt.0.) gap = alamd0*(1.519-sqrt(dummy))	

	write (10,203) hanue(imesh,icurve), brill(imesh,icurve),
	1   flux(imesh,icurve), fluxdens(imesh,icurve), ak1, ptot, gap
203	format (3x,e12.4, 3x,e12.4, 3x,e12.6, 3x,e12.6,
	1 3x,f7.3, 3x,e10.3, 2x,f8.6)


c	hanl1 = alog10(hanue1)
c	hanl2 = alog10(hanue2)
c	blog1 = alog10(b1)
c	blog2 = alog10(b2)
c	if (hanl1.lt.hnminl  .or. hanl2.gt.hnmaxl ) goto 200
c	if (blog1.gt.bmaxl   .or. blog2.gt.bmaxl  ) goto 200
c	if (blog1.lt.bminl   .or. blog2.lt.bminl  ) goto 200
c	y1 = yscale*(blog1-bminl)
c	x1 = xscale*(hanl1-hnminl)
c	y2 = yscale*(blog2-bminl)
c	x2 = xscale*(hanl2-hnminl)
c	if (y1.lt.0. .or. y2.lt.0.) goto 200
c	call vector (x1,y1, x2,y2)

200	continue
	icurve = icurve+1
	ncurve = ncurve+1

	if (sourceprop.gt.0.) then
	write (11, 211)
211	format (//,' hanue,   SIGMAX,    sigmax,    sigmar ')
	write (12, 212)
212	format (//,' hanue,  SIGMAX",   sigmax",    sigmar" ')
	write (13, 213)
213	format (//,' hanue,   SIGMAY,    sigmay,    sigmar ')
	write (14, 214)
214	format (//,' hanue,  SIGMAY",   sigmay",    sigmar" ')
	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif

	write (6,301)
301	format ('Enter another harmonic: ')
715     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 715
        else
          read(cline,*) smallk
        endif
	if (smallk.le.0.) return
	goto 109
	end
*CMZ :  1.02/01 16/04/2015  15.44.19  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.13.51  by  Michael Scheer
*CMZ :  1.01/00 12/09/2013  15.01.41  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.12.29  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.32.20  by  Andreas Gaupp
*-- Author :    Michael Scheer
c	****************************************************************
c			  ELLIPTICAL UNDULATOR
c       Based on ellipticalundu by Andreas Gaupp
c	****************************************************************
	subroutine ellipticalundu_walker
c 	expanded for elliptical undulator 04/12 AG
c	normalization changed according to K.J.Kim  April 11, 1988

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

        character a

	write (6,998)
998	format (' Elliptical undulator under development Dec 2004')

1	write (6,20)
20	format (' enter periode length in m: ')
778     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 778
        else
          read(cline,*) alamd0
        endif
	write (6,10)
10	format (' enter number of undulator periodes: ')
779     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 779
        else
          read(cline,*) an
        endif
	if ( an*alamd0.ge.alength) write (6,21) alamd0, an, alength
21	format (' **** warning ****',/
	1 ' lamda0  = ',f10.2,' m',/
	2 '      N  = ',f10.0,/
	3 ' alength = ',f10.2,' m '/
	4 ' *** not compatible ***'/)
	
	write (6,30)
30	format(' enter range in effective wiggler strength
     1	 KeffMIN and KeffMAX: ')
710      read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 710
        else
          read(cline,*) akmin, akmax
        endif
	if (akmin.ge.0. .and. akmax.ge.akymin) goto 32
	write (6,31) akmin, akmax
31	format (' **** warning ****',/
	1 ' AKMIN   = ',f10.3,'   AKMAX  = ',f10.3,' *** inconsistent'/)
32	continue
	write (6,35)
35	format (' enter ratio of horizontal to vertical wiggler '
     1  ' strength akxky = q =  Kx/Ky.'/,
     1  ' P**2 I in APPLE II is optimized for ',/
     2  ' first harmonic    q = 0.9999',/
     3  ' third harmonic    q = 0.42  ',/
     4  ' fifth harmonic    q = 0.32  ',/
     5  ' seventh harmonic  q = 0.27  ',/
     6  ' nineth harmonic   q = 0.24  ',/
     7  ' eleventh harmonic q = 0.22  ',/
     9  /)
711     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 711
        else
          read(cline,*) akxky
        endif
	write (6,'(" Kx/Ky = ",e12.6)') akxky

	if (akxky.gt.1.) then
	write (6,'(//," ?????????  in Ellipticalundu Kx/Ky  must be <=1 ",//)')
	write (10,'(//," ????????? in Ellipticalundu Kx/Ky  must be <=1 ",//)')
	endif

	akymin=akmin/(sqrt(1.+akxky**2))
	akymax=akmax/(sqrt(1.+akxky**2))

	write (6,40)
40	format (' enter order of undulator spectrum (1, 3, 5,..): ')
713     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 713
        else
          read(cline,*) smallk
        endif
	write (6,'(" smallk = ",f10.3)') smallk

	if (smallk.lt.1. .or.
	1   abs(smallk-float(int(smallk))).gt.1.e-4) write (6,41) smallk
41	format (' **** warning ****',/
	1 ' bad value for order of spectrum. SMALLK = ',f10.3)
2	continue

c	get minimum gap for hybrid-REC structure using Halbach's formula
c	for REC without 10% reduction
	gapmin = alamd0*(1.519-sqrt(2.309+.5556*alog(.003223*akymax/alamd0)))	
	write (6,100) an, alamd0, akymin, akymax, akxky, smallk, gapmin
100	format (/,' undulator parameters',/
	1 ' number of periodes             AN      = ',f10.0,/
	2 ' periode length (m)             alamd0  = ',f10.5,/
	3 ' minimum K                      KeffMIN = ',f10.2,/
	4 ' maximum K                      KeffMAX = ',f10.2,/
	4 ' ratio Kx/Ky                    akxky   = ',f10.2,/
	5 ' order of spectrum              smallk  = ',f10.0,/
	6 ' minimum gap (hybrid REC) (m)   gapmin  = ',f10.5,//
	9 '$ correct? (y/n)  ')
714       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 714
        else
          a=cline(1:1)
        endif
cmsh	read(5,101) a
101	format (a1)
	if (a.eq.'y') goto 109
	goto 1
109	continue

c	plot undulator spectrum
	write (6,110)
110	format (' enter photon energy range (eV) of this curve',/
	1 '$ hanuemin and hanuemax:  ')
715       read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 715
        else
          read(cline,*) hanuemin, hanuemax
        endif
	if (hanuemax.le.0.) then
		hanuemax = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*
	2		(1.+0.5*aKymin**2*(1.+akxky**2))))
	endif
	if (hanuemin.le.0.) then
		hanuemin = 1240.e-9 * smallk *
	1		((2.*gamma**2)/(alamd0*
	2		(1.+0.5*aKymax**2*(1.+akxky**2))))
	endif
	write (6,120)
120	format ('Enter number of points on this curve, mesh:  ')
716     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 716
        else
          read(cline,*) mesh(icurve)
        endif
	if (mesh(icurve).gt.256) mesh(icurve) = 256

c	calculate normalisation of brilliance for K = 0
c	only for initial printout, not for final curve
	theta0 = 1./(gamma*sqrt(an*smallk))
c	note alength is in m
	sigmax = sqrt(epsx*betax)
  	sigmay = sqrt(epsy*betay)
	sigmaxp= sqrt(theta0**2+epsx/betax)
	sigmayp= sqrt(theta0**2+epsy/betay)
	write (6,121) theta0, sigmax, sigmay, sigmaxp, sigmayp, dgamma
121	format (/
	1 ' theta0 (rad, K=0) = ',e10.4,/
	2 ' sigmax (m)        = ',e10.4,/
	3 ' sigmay (m)        = ',e10.4,/
	4 ' sigmaxp (rad) K=0 = ',e10.4,/
	5 ' sigmayp (rad) K=0 = ',e10.4,/
	6 ' sigma E           = ',e10.4,//)
	
c	get maximum power according to K.J.Kim
	call power (an, aKymax, alamd0, gamma, curr, Ptot, d2P)


	write (10,122) alamd0, an, akymin, akymax, akxky, smallk, gapmin,
	1              theta0,sigmax, sigmay,sigmaxp,sigmayp,dgamma
	write (10,123) Ptot, d2P

122	format (/
	1 ' properties of undulator',/
	1 ' -----------------------',/
	2 5x,' periode length (m)             	',f10.4,/
	3 5x,' number of periodes              	',f10.4,/
	4 5x,' minimum Ky                      	',f10.4,/
	5 5x,' maximum Ky                       ',f10.4,/
	5 5x,' ratio Kx/Ky                      ',f10.4,/
	6 5x,' order of spectrum k              ',f10.4,/
	7 5x,' minimum gap (REC hybrid, m)      ',f10.4,//
	1 ' properties of radiation',/
	2 5x,' theta0 (rad) K = 0               ',e10.4,/
	3 5x,' sigmax (m)                       ',e10.4,/
	4 5x,' sigmay (m)                       ',e10.4,/
	5 5x,' sigmaxp (rad) K = 0              ',e10.4,/
	6 5x,' sigmayp (rad) K = 0              ',e10.4,/
	7 5x,' dgamma                           ',e10.4)
123	format (
	7 5x,' Ptot (Watt)                      ',e10.4,/
	7 5x,' dP/dOmega (Watt/mrad^2)          ',e10.4,/
	9 5x,' h*nue (eV)         Brilliance              Flux         K
	1       Ptot    gap',/
	2 5x '                       photon              photon
	3     (Watt)     (m)',/
	4 5x '          sec*(mm*mrad)**2*0.1% BW      sec*0.1% BW  ',/)

      if (sourceprop.gt.0.) then
        call fileopen
        write (11,*) mesh(icurve)-1
        write (12,*) mesh(icurve)-1
        write (13,*) mesh(icurve)-1
        write (14,*) mesh(icurve)-1

c	e beam cross section and divergence
        sigmax = sqrt(epsx*betax)
        sigmay = sqrt(epsy*betay)
        sigxpr = sqrt(epsx/betax)
        sigypr = sqrt(epsy/betay)

      endif

      do 200 imesh=1,mesh(icurve)
        hanue1 = hanuemin+(hanuemax-hanuemin)*(imesh-1)/float(mesh(icurve)-1)
c	8065 cm**-1  = 1 eV
        ak1 = 2.*(smallk*2.*gamma**2/(hanue1*alamd0*806500.)-1.)
        if (ak1.le.0.) goto 200

        aky1 = sqrt(ak1/(1.+akxky**2))
        akx1 = aky1*akxky
	if (aky1.lt.akymin .or. aky1.gt.akymax) goto 200

c	formula due to Krinsky in Koch's book, eq.281, p.152
c	checked on April 9, 1985
c	units are  photons/[sec*mrad**2*0.1%BW*Amp])
c	our efka is the same as F_n(K) in xray data booklet
c	b1 = 4.56E7*(an*gamma)**2*curr*efka(smallk,ak1)
        b1= 4.56E7*(an*gamma)**2*curr*efkaelli(aky1,akx1)
c	this is the spectral brightness, which is reduced by a
c	finite energy spread
c	here we take the spectral distribution of undulator radiation
c	to be a normal with width sigma = 1./(smallk*2.6*N)
c	to be a normal with width sigma = 1./(smallk*2.0d0*sqrt(2.0d0)*N) !msh
c	The energy spread of the e beam is multiplied by 2 because
c	gamma**2 enters the energy formula
c	This is necessary as discovered by Howard Padmore at ALS in
c	winter 1993/94.                         AG 1994/01/16

cmsh        b1= b1/sqrt(1.+(2.*dgamma*smallk*2.6*aN)**2)
        b1= b1/sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2)
        dgamred=1./sqrt(1.+(2.*dgamma*smallk*2.0d0*sqrt(2.0d0)*aN)**2) !msh

cmsh{

c	save on axis spectral brightness into array fluxdens
c	need to correct for electron beam divergence
c	e beam divergence
        sigxpr = sqrt(epsx/betax)
        sigypr = sqrt(epsy/betay)
c	natural divergence of radiation
        alamda = (alamd0/(smallk*2.*gamma**2))*(1.+(aky1**2+akx1**2)/2.)
        aL = alamd0*an

        sigmarp=sqrt(alamda/(2.0d0*aL))

        print*,'*** Fluxdens unclear in ellipticalundu_walker! ***'
        fluxdens(imesh,icurve) = b1/
     &    sqrt(1.+(sigxpr**2+sigypr**2)/sigmarp**2)

cmsh}

c	flux(imesh,icurve) is the total flux at energy curvex
c	here we make a very crude approximation,
c	since we calculate the flux by multiplying the angular brightness
c	by 2*pi*theta0**2
c	our formula is the same as eq.17 in xray data booklet
c	the formula for sigmar is consistent with Walker
c	ESRP - IRM - 54/84
c	***** sigmar = sqrt (lambda / L) *****

c	J.B. 30.5.89
        sigmarp=1./gamma*sqrt((1.+(akx1**2+aky1**2)/2.)/(smallk*2.*an))/sqrt(2.0d0)
        flux(imesh,icurve) = 2.*pi*b1*(sigmarp*1000.)**2
        hanue (imesh, icurve) = hanue1

        call normuelli_walker(anorm, aky1,akx1)

	if (imesh.gt.256 .or. icurve.gt.100 .or. anorm.le.0.)
     &    write (6,*) ' *** imesh=',imesh,'   icurve=',icurve,
     &    ' *** anorm=',anorm
        b1= flux(imesh,icurve)/anorm
        brill(imesh,icurve)= b1

	if (sourceprop.gt.0.) then
c	natural divergence of radiation
          alamda = (alamd0/(smallk*2.*gamma**2))*(1.+(aky1**2+akx1**2)/2.)
          aL = alamd0*an
          sigmarp = sqrt(alamda/(2.0d0*aL))
c	radiation cross section
          sigmar = sqrt(2.0d0*alamda*al)/(2.0d0*pi)
c	MSH behauptet basierend auf einem Vergleich der zweiten Momente
c	(2000/06):
c	sigmar=(1/2 pi) sqrt(2 lambda L)
c	sigmar * sigmarp = (1/2 pi) lambda (*0.82)

          size_xt = sqrt(sigmax**2+sigmar**2)
          dive_xt = sqrt(sigxpr**2 + sigmarp**2)
          size_yt = sqrt(sigmay**2+sigmar**2)
          dive_yt = sqrt(sigypr**2 + sigmarp**2)

          write (11,*) hanue1, size_xt, sigmax, sigmar
          write (12,*) hanue1, dive_xt, sigxpr, sigmarp
          write (13,*) hanue1, size_yt, sigmay, sigmar
          write (14,*) hanue1, dive_yt, sigypr, sigmarp
        endif

c	#########################
c	coherent flux
c	#########################
	if (coherent.gt.0.) then
c	use formula due to K.J.Kim, LBL-22236 (march 1987)
c	published in PAC1987, p 194
c	converting into units photons/sec/0.1%BW I get
c	coherent Flux = .385 Photons/sec/0.1%BW
c	                * B(0,0)/(phot/sec/mm**2 mrad**2 0.1%BW)
c	                * (1/(hnue/eV)**2)

          flux(imesh,icurve)= .385*b1/(hanue1**2)
	endif

c Die Eneriebreite des Strahles hat kaum Einfluss auf den Fluss, deshalb
        flux(imesh,icurve)=flux(imesh,icurve)/dgamred

c	write (6,201)
c	1 hanue(imesh,icurve), ak1, smallk, brill(imesh,icurve)
c201	format (' hanue/eV =',f9.3,' K=',f9.3,' k=',f2.0,' B=',e10.4,
c	1 ' ph/(sec (mm*mrad)**2 .1%)')
c
c	total power in Watt according to K.J.Kim
	ptot = 1.774e-8*gamma**2*(aky1**2+akx1**2)*an*curr/alamd0

	dummy = 2.309+.5556*alog(.003223*aky1/alamd0)
	gap = 0.
	if (dummy.gt.0.) gap = alamd0*(1.519-sqrt(dummy))	

	write (10,203) hanue(imesh,icurve), brill(imesh,icurve),
	1   flux(imesh,icurve), ak1, ptot, gap
203	format (3x,e12.4,8x,e12.4,8x,e12.6,3x,f7.3,3x,e10.3,2x,f8.6)

200	continue
	icurve = icurve+1
	ncurve = ncurve+1

	if (sourceprop.gt.0.) then
	write (11, 211)
211	format (//,' hanue,   SIGMAX,    sigmax,    sigmar ')
	write (12, 212)
212	format (//,' hanue,  SIGMAX",   sigmax",    sigmar" ')
	write (13, 213)
213	format (//,' hanue,   SIGMAY,    sigmay,    sigmar ')
	write (14, 214)
214	format (//,' hanue,  SIGMAY",   sigmay",    sigmar" ')
	close (unit=11)
	close (unit=12)
	close (unit=13)
	close (unit=14)
	endif

	write (6,301)
301	format ('Enter another harmonic: ')
717     read(5,'(a)')cline
        if(cline(1:1).eq.''.or.cline(1:1).eq.'*'.or.cline(1:1).eq.'#'
     &      .or.cline(1:1).eq.'!'.or.cline(1:1).eq.'%') then
          goto 717
        else
          read(cline,*) smallk
        endif
	if (smallk.le.0.) return
	goto 109
	end
*CMZ :  1.02/01 08/04/2015  15.55.31  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.09  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.14.24  by  Michael Scheer
*CMZ :  1.01/00 12/09/2013  15.01.41  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.10.31  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.35.17  by  Andreas Gaupp
*-- Author :    Michael Scheer
c	************************************************************
c				N O R M U
c	************************************************************
	subroutine normu_walker(anorm, ak1)

c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

	sigmax = sqrt(epsx*betax)
	sigxpr = sqrt(epsx/betax)

c	natural divergence of radiation
	alamda = (alamd0/(smallk*2.*gamma**2))*(1.+ak1**2/2.)
        aL = alamd0*an

        sigmarp = sqrt(alamda/(2.0d0*aL))
        sigmar = sqrt(2.0d0*alamda*al)/(2.0d0*pi)

	anormx = 1.e6*2.*pi*sqrt((sigmax**2 + sigmar**2) *
     &    (sigxpr**2 + sigmarp**2))

	sigmay = sqrt(epsy*betay)
	sigypr = sqrt(epsy/betay)

	anormy = 1.e6*2.*pi*sqrt((sigmay**2 + sigmar**2)*
     &    (sigypr**2 + sigmarp**2))

	anorm = anormx*anormy

        if (depthoffield.ne.0. ) then
          print*,'*** Warning Depth of Field not implemented in normu_walker ***'
          stop '*** brill aborted ***'
c 	section added October 1992
c	This section is to include depth of field affects according to
c	K.J.Kim, LBL-22236, March 1987, published in PAC 1987, p.194
c	The philosophy is as follows:
c	the on axis brightness is written as an integral along the
c	undulator axis. The contribution from a path element dz
c	is assumed to be a gaussian times a form factor g(z).
c	The free parameters are choosen to satisfy several criteria
c	which are:
c	- single electron flux is completely coherent
c	- spatial flux density is the brilliance integrated over angle
c	- angular flux density is the brilliance integrated over source area
c	- Total flux is double integral over brilliance

c	the last condition is met within 7%

c	2004/03/17: Comparing to WAVE. In one case numerical agreement
c	is obtained between WAVE (author: M.Scheer, BESSY) and BRILL
c	with depth of field effect. The depth of field effect reduces
c	the brilliance by a factor of 5 in that case.
	endif			! of depth of field

	return
	end
*CMZ :  1.02/01 08/04/2015  15.59.15  by  Michael Scheer
*CMZ :  1.02/00 20/09/2013  12.06.08  by  Michael Scheer
*CMZ :  1.01/01 12/09/2013  15.13.51  by  Michael Scheer
*CMZ :  1.01/00 12/09/2013  15.01.41  by  Michael Scheer
*CMZ :  1.00/02 06/09/2013  09.12.29  by  Michael Scheer
*CMZ :  0.01/00 03/09/2013  12.44.33  by  Michael Scheer
*CMZ :  0.00/02 28/06/2011  09.28.02  by  Michael Scheer
*CMZ :  0.00/01 27/06/2011  16.32.20  by  Andreas Gaupp
*-- Author :    Michael Scheer
c	************************************************************
c				N O R M U E L L I
c	************************************************************
	subroutine normuelli_walker(anorm, aky1,akx1)

c Replaces normuelli according to  Walker's formulars, 8.3.2015
c	include 'brill.cmn'
*KEEP,brill.
c	all common for code BRILL
      real*8       llab(8)

      real brill,flux,hanue,fluxdens
      integer icurve,ncurve,mesh

      common brill (256,100), flux (256,100), hanue (256,100),
     &  fluxdens(256,100), icurve, ncurve, mesh(100)
	
      real pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      common /main/ pi,
     &  gamma   , dgamma  , epsx    , epsy    , curr    , rho  ,
     &  alength , betax   , betay   ,
     &  an      , alamd0  , akmin   , akmax   ,
     &  akymin  , akymax  , akxmin  , akxmax  , akxky   ,
     &  smallk , gapmin ,
     &  depthoffield, coherent, sourceprop
	
      real bjre, bjim, yre, yim, pi1
      common /bessel/ bjre(128), bjim(128), yre(128), yim(128), pi1

      integer ndatfil
      common/undul/ndatfil
	
      real xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb, ylabb, xlabf, ylabf
      integer
     &  labsize, nlab,
     &  iunit		

      common/plt/ xorig, yorig, xr1, yr1, xr2, yr2, xpap, ypap,
     &  xscale, ybscale, yfscale,
     &  hnminl, hnmaxl, bminl, bmaxl, fminl, fmaxl,dummy,
     &  xlabb( 100), ylabb( 100), xlabf( 100), ylabf( 100),
     &  lab( 8, 100), labsize(100), nlab,
     &  iunit		

c	real*8       llab(8)
      byte lab
      byte datstring (9), timstring(8)
      character*80 filename_xs, filename_xd, filename_ys, filename_yd
      character*3 number

      equivalence (lab(1,1),llab(1))

      character(256)cline
      common/clinec/cline
*KEND.

c	sigmar  = (1/2pi) sqrt (2*lambda*L)
c	sigmarp =         sqrt (lambda/(2L))
c
	sigmax = sqrt(epsx*betax)
        sigxpr = sqrt(epsx/betax)

c	natural divergence of radiation
        alamda = (alamd0/(smallk*2.*gamma**2))*(1.+(aky1**2+akx1**2)/2.)
        aL = alamd0*an

        sigmarp = sqrt(alamda/(2.0d0*aL))
        sigmar = sqrt(2.0d0*alamda*al)/(2.0d0*pi)

        anormx = 1.e6*2.*pi*sqrt((sigmax**2 + sigmar**2) *
     &    (sigxpr**2 + sigmarp**2))

	sigmay = sqrt(epsy*betay)
	sigypr = sqrt(epsy/betay)

	anormy = 1.e6*2.*pi*sqrt((sigmay**2 + sigmar**2)*
     &    (sigypr**2 + sigmarp**2))

        anorm = anormx*anormy

        if (depthoffield.ne.0.0) then
          stop '*** Error in normuelli_walker: Depth of field not implemented ***'
	endif

	return
	end
