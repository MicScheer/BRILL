
##############################################################################
#
#      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
#      Hahn-Meitner-Platz 1
#      D-14109 Berlin
#      Germany
#
#      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
#
# -----------------------------------------------------------------------
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy (wave_gpl.txt) of the GNU General Public
#    License along with this program.
#    If not, see <http://www.gnu.org/licenses/>.
#
#    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
#    der GNU General Public License, wie von der Free Software Foundation,
#    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
#    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
#
#    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
#    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
#    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
#    Siehe die GNU General Public License fuer weitere Details.
#
#    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
#    zusammen mit diesem Programm erhalten haben. Wenn nicht,
#    siehe <http://www.gnu.org/licenses/>.
#
##############################################################################

#+PATCH,//BRILL/PYTHON
#+DECK,pyDipole,T=PYTHON.

import os,sys,platform,shutil,time,re

import tkinter as tk
from tkinter import *

import numpy as np
import pandas as pd

from scipy import special

import matplotlib as mpl
import matplotlib.pyplot as plt


global \
clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
ecdipev,ecdipkev,g1const,g1max,h2const,h2max

g1max=0.9212
g1const=2.457e13
h2max=1.474
h2const=1.327e13

hbarev1=6.58211889e-16
clight1=2.99792458e8
emasskg1=9.10938188e-31
emasse1=0.510998902e6
emassg1=0.510998902e-3
echarge1=1.602176462e-19
erad1=2.8179380e-15
eps01=8.854187817e-12
pi1=3.141592653589793e0
grarad1=pi1/180.e0
radgra1=180.e0/pi1
hplanck1=6.62606876e-34
hbar1=hbarev1*echarge1
wtoe1=clight1*hplanck1/echarge1*1.e9
cq1=55.e0/32.e0/(3.0e0)**0.5*hbar1/emasskg1/clight1
cgam1=4.e0/3.e0*pi1*erad1/emassg1**3
pol1con1=8.e0/5.e0/(3.0e0)**0.5
pol2con1=8.e0/5.e0/(3.0e0)**0.5/2.e0/pi1/3600.e0*emasskg1/hbar1/erad1*emassg1**5

twopi1=2.0e0*pi1
halfpi1=pi1/2.0e0
sqrttwopi1=(twopi1)**0.5
dnull1=0.0e0
done1=1.0e0
rmu01=4.0e0*pi1/1.0e7
rmu04pi1=1.0e-7
alpha1=echarge1**2/(4.0e0*pi1*eps01*hbar1*clight1)
gaussn1=1.0e0/(twopi1)**0.5

ck934=echarge1/(2.0e0*pi1*emasskg1*clight1)/100.0e0


# +PATCH,//WAVES/PYTHON
# +KEEP,mshutil,T=PYTHON.


global Klold, Nold, Xa1old, Xanold

Klold = 1
Nold = -99
Xa1old = -9999.
Xanold = -9999.

import numpy as np

def pe2(x):
  try: return '{:.2e}'.format(float(x))
  except: return x
def pe3(x):
  try: return '{:.3e}'.format(float(x))
  except: return x
def pe4(x):
  try: return '{:.4e}'.format(float(x))
  except: return x
def pe5(x):
  try: return '{:.5e}'.format(float(x))
  except: return x
def pe6(x):
  try: return '{:.6e}'.format(float(x))
  except: return x
def pe7(x):
  try: return '{:.7e}'.format(float(x))
  except: return x
def pe8(x):
  try: return '{:.8e}'.format(float(x))
  except: return x
def pe9(x):
  try: return '{:.9e}'.format(float(x))
  except: return x
def pe10(x):
  try: return '{:.10e}'.format(float(x))
  except: return x

def pg2(x):
  try: return '{:.2g}'.format(float(x))
  except: return x
def pg3(x):
  try: return '{:.3g}'.format(float(x))
  except: return x
def pg4(x):
  try: return '{:.4g}'.format(float(x))
  except: return x
def pg5(x):
  return '{:.5g}'.format(float(x))
#  try: return '{:.5g}'.format(float(x))
#  except: return x
def pg6(x):
  try: return '{:.6g}'.format(float(x))
  except: return x
def pg7(x):
  try: return '{:.7g}'.format(float(x))
  except: return x
def pg8(x):
  try: return '{:.8g}'.format(float(x))
  except: return x
def pg9(x):
  try: return '{:.9g}'.format(float(x))
  except: return x
def pg10(x):
  try: return '{:.10g}'.format(float(x))
  except: return x

def readfloat(s,default=-9999):
  print(s)
  ans = input()
  if ans == '': return default
  else: return float(ans)
#enddef

def readint(s,default=-9999):
  print(s)
  ans = input()
  if ans == '': return default
  else: return int(float(ans))
#enddef

def printnl(line):
  print("\n",line,"\n")
#enddef printnl()

def set_console_title(console='Python'):
#+seq,mshimportsind.
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  # Histograms and Ntuples
  global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
  Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
  H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
  Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
  Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
  H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
  FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman
#+KEEP,plotglobind,T=PYTHON.
#*CMZ :          28/09/2019  14.39.13  by  Michael Scheer
  global MPLmain, MPLmaster, Nfigs,Figgeom, Figgeom2, FiggeomR, FiggeomL, XtermGeo, Figs,Fig,Ax,\
  Fig1,Ax1,Fig6,Ax6,Fig2,Ax2,Fig7,Ax7,Fig3,Ax3,Fig8,Ax8, Figgeoms, \
  Fig4,Ax4,Fig9,Ax9,Fig5,Ax5,Fig10,Ax10,\
  Screewidth, Screenheight, ScaleSizeX, ScaleSizeY, \
  FirstConsole, Console, Igetconsole,Klegend, Fwidth, Fheight, Fxoff, Fyoff, \
  Kfig, Kax, Ihist,Iprof, Imarker, Ierr, Isurf, Iinter, Isame, Itight, IsameGlobal, Iline, CMap, Cmap, Tcmap, Surfcolor, Cmaps, \
  Iplotopt, Ispline, Kecho, Kdump,Kpdf, Ndump,Npdf, Legend, \
  Kplots,Nwins, Zones, Kzone, Nxzone, Nyzone, Zone, Axes, Icmap, \
  Mode3d,Mode3D, Mode2d,Mode2D, CanButId, CanButIds, \
  MarkerSize, MarkerType, MarkerColor, \
  Markersize, Markertype, Markercolor, \
  Fillstyle, FillStyle, \
  Textcolor, WaveFilePrefix,WaveDump, \
  LineStyle, LineWidth, LineColor, \
  Linestyle, Linewidth, Linecolor, \
  Author, \
  Tightpad, Xtightpad,Ytightpad, ColorbarPad,\
  LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace, \
  Histcolor, Histedgecolor, Histbarwidth, Kdate, Kfit, Kstat, YTitle, \
  Icont3d, Iboxes, Inoempty, Iclosed,Itrisurf, Iscatter, Iscat3d, Ifill1d, TitPad, Xtitle, Ytitle, \
  Gtit,Xtit,Ytit,Ztit,Ttit,Ptit,Colors, Surfcolors,Linestyles, Markertypes, \
  LexpX,LexpY,LexpRot,LexpPow,\
  GtitFontSize,Titfontsize,Atitfontsize,Axislabelsize,Textfontsize,Datefontsize,\
  Statfontsize, Axislabeldist, Axislabeldist3d, Axisdist, Axisdist3d, \
  XFit, YFit, Xfit, Yfit,Ystat, YStat, \
  GtitFontSize,TitFontSize,AtitFontSize,AxisLabelSize,TextFontSize,DateFontSize,\
  StatFontSize, AxisLabelDist, AxisLabelDist3d, AxisTitleDist, AxisTitleDist3d, \
  AtitFontSize3d, Atitfontsize3d, NXtick,NXtick3d, Nxtick,Nxtick3d, Ktitles,  Dummy,\
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,ZoomZmin,ZoomZmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull2D,Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux,Ishow

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ,Tnpa,Tnone

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  Console = console

  if platform.system() == 'Linux':
    sys.stdout.write("\x1b]2;" + console + "\x07")
  elif platform.system() == 'Windows':
    #ctypes.windll.kernel32.SetConsoleTitleW(console)
    system("title "+console)
  #endif

#enddef set_console_title()

def util_spline_coef_periodic(x,y):

  """
c--- calculates spline coefficients for periodic function
c--- the interval must be closed, i.e. x(n)-x(1)=periodlength and y(n)=y(1)

c--   input:

c-       n: number of x,y-values, must be at least five
c-       x: array of x-values
c-       y: array of y-values

c--   ouput:

c-       ypp:   spline-coefficients
c-     ifail:   error status

c--   workingspace: aa(n),bb(n),cc(n),c(n),cn(n)
  """
  n = len(x)

  ypp =  np.zeros_like(x)

  aa = np.zeros_like(x)
  bb = np.zeros_like(x)
  cc = np.zeros_like(x)
  c  = np.zeros_like(x)
  cn = np.zeros_like(x)

  if n < 5:
    ifail=-1
    return ifail, ypp
  #endif

  ymax=0.0
  ya = abs(y)
  ymax = ya.max()

  n -= 1 # now n is last index

  if abs(y[n]-y[0])/ymax > 1.0e-9:
    ifail=-2
    return ifail, yp
  else:
    y[0]=(y[n]+y[0])/2.0
    y[n]=y[0]
  #endif

  #cn=0.0 #letzte Spalte der Matrix

  n1=n-1
  n2=n-2

  aa[0]=(x[n]-x[n1])/6.e0
  bb[0]=(x[1]-x[0]+(x[n]-x[n1]))/3.e0
  cc[0]=(x[1]-x[0])/6.e0
  c[0]=(y[1]-y[0])/(x[1]-x[0])-(y[n]-y[n1])/(x[n]-x[n1])

  for j in range(1,n):
    #do j=2,n1
    aa[j]=(x[j]-x[j-1])/6.0
    bb[j]=(x[j+1]-x[j-1])/3.0
    cc[j]=(x[j+1]-x[j])/6.0
    c[j]=(y[j+1]-y[j])/(x[j+1]-x[j])-(y[j]-y[j-1])/(x[j]-x[j-1])
    #print(j,aa[j],bb[j],cc[j],c[j])
    #print(j,x[j],y[j])
  #enddo #j

  # Auf Dreiecksmatrix bringen

  # Oberste Zeile

  cc[0]=cc[0]/bb[0]
  aa[0]=aa[0]/bb[0]
  c[0] = c[0]/bb[0]
  bb[0]=1.0

  # 2. Zeile, d.h. die erste regulaere, cn[j] ist die letzte Spalte der Matrix

  bb[1]=bb[1]/aa[1]
  cc[1]=cc[1]/aa[1]
  c[1] = c[1]/aa[1]
  aa[1]=1.0

  aa[1]=0.0
  bb[1]=bb[1]-cc[0]
  cn[1]=-aa[0]
  c[1] = c[1]-c[0]

  cc[1]=cc[1]/bb[1]
  cn[1]=cn[1]/bb[1]
  c[1] = c[1]/bb[1]

  bb[1]=1.0

  # Nun die hoeheren bis n-3

  for j in range(2,n-3):
    #do j=3,n-3

    bb[j]=bb[j]/aa[j]
    cc[j]=cc[j]/aa[j]
    c[j] = c[j]/aa[j]

    aa[j]=0.0
    bb[j]=bb[j]-cc[j-1]
    cn[j]=-cn[j-1]
    c[j] = c[j]-c[j-1]

    cc[j]=cc[j]/bb[j]
    cn[j]=cn[j]/bb[j]
    c[j] = c[j]/bb[j]

    bb[j]=1.0

  #enddo

  # vorletzte zeile

  bb[n2]=bb[n2]/aa[n2]
  cc[n2]=cc[n2]/aa[n2]
  c[n2] = c[n2]/aa[n2]
  aa[n2]=1.0

  aa[n2]=0.0
  bb[n2]=bb[n2]-cc[n2-1]
  cc[n2]=cc[n2]-cn[n2-1]
  c[n2] = c[n2]-c[n2-1]

  cc[n2]=cc[n2]/bb[n2]
  c[n2] = c[n2]/bb[n2]

  bb[n2]=1.0

  # Letzte Zeile

  ypp[n2]=aa[n1]/cc[n1]
  ypp[n1]=bb[n1]/cc[n1]
  c[n1]=c[n1]/cc[n1]
  cc[n1]=1.0
  ypp[0]=cc[n1]

  # Oberste Zeile abziehen

  ypp[0]=0.0
  ypp[1]=-cc[0]
  ypp[n1]=ypp[n1]-aa[0]
  c[n1]=c[n1]-c[0]

  for j in range(1,n1):
    #do j=2,n2

    c[n1]=c[n1]/ypp[j]
    ypp=ypp/ypp[j]

    ypp[j]=ypp[j]-bb[j]
    ypp[j+1]=ypp[j+1]-cc[j]
    ypp[n1]=ypp[n1]-cn[j]
    c[n1]=c[n1]-c[j]

  #enddo

  c[n1]=c[n1]/ypp[n1]

  # Ernten

  ypp[n1]=c[n1]

  # Vorletzte Zeile

  bb[n2]=bb[n2]/cc[n2]
  c[n2]=c[n2]/cc[n2]
  cc[n2]=1.0

  c[n2]=c[n2]-c[n1]
  cc[n2]=0.0

  c[n2]=c[n2]/bb[n2]
  bb[n2]=1.0

  ypp[n2]=c[n2]

  # Letzte Spale nullen und normieren, letzte und vorletzte sind bereits fertig

  bb[0]=bb[0]/aa[0]
  cc[0]=cc[0]/aa[0]
  c[0]=c[0]/aa[0]
  aa[0]=1.0

  c[0]=c[0]-c[n1]
  aa[0]=0.0

  cc[0]=cc[0]/bb[0]
  c[0]=c[0]/bb[0]
  bb[0]=1.0


  # Regulaere Zeilen
  for j in range(1,n-3):
    #do j=2,n-3

    bb[j]=bb[j]/cn[j]
    cc[j]=cc[j]/cn[j]
    c[j]=c[j]/cn[j]
    cn[j]=1.0

    c[j]=c[j]-c[n1]
    cn[j]=0.0

    cc[j]=cc[j]/bb[j]
    c[j]=c[j]/bb[j]
    bb[j]=1.0

  #enddo

  j=n-2
  while j > 1:
    j -= 1
    #do j=n-3,2,-1
    ypp[j]=c[j]-cc[j]*ypp[j+1]
  #enddo

  ypp[0]=(c[0]-cc[0]*ypp[1]-aa[0]*ypp[n1])/bb[0]
  ypp[n]=ypp[0]

  ifail=0

  return ifail, ypp
#enddef util_spline_coef_periodic(x,y)

def util_spline_inter(xa,ya,y2a,x,mode):


  global Klold, Nold, Xa1old, Xanold


  """
C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       YA2:  ARRAY SPLINE COEFFICIENTS
C-       X: Y(X) IS CALCULATED
C-       MODE: CONTROL FLAG:
C-             MODE.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             MODE.LT.0: NEW INITIALIZATION

C--   OUTPUT: Y(X), DY/DX(X), D2Y/DX2(X)

  """
  y = None
  yp = None
  ypp = None

  n = len(xa) - 1

  eps=abs(xa[n]-xa[0])/1.0e10
  xx=x

  if xa[0] > xa[n]:
    print('*** Error in util_spline_inter: x-values must be in ascending order ***')
    Quit()
  #endif xa[0] > xa[n]:

  if xx < xa[0] and xx > xa[0]-eps:
    xx=xa[0]
  elif xx > xa[n] and xx < xa[n]+eps:
    xx=xa[n]
  #endif

  if xx < xa[0] or xx > xa[n]:
    print('xa[0], xa[n]:',xa[0], xa[n])
    print('x:',x)
    print('*** Error in util_spline_inter: x out of range ***')
    return y,yp,ypp
  #endif

  if mode < 0 or Klold >= n:
    klo=0
  elif Nold == n and  xa[0] == Xa1old and  xa[n] == Xanold and  xx > xa[Klold]:
    klo=Klold
  else:
    klo=0
  #endif

  if xx < xa[klo+1]:
    khi=klo+1
  else:
    khi=n
    while (khi-klo) > 1:
      k=int((khi+klo)/2)
      if xa[k] > xx:
        khi=k
      else:
        klo=k
      #endif
    #endwhile
  #endif

  h=xa[khi]-xa[klo]

  if h <= 0.0:
    print('*** error in util_spline_inter: bad input ***')
    return y,yp,ypp
  #endif

  a=(xa[khi]-xx)/h
  b=(xx-xa[klo])/h

  yl = ya[klo]
  yh = ya[khi]
  y2l = y2a[klo]
  y2h = y2a[khi]

  y = a*yl+b*yh+(a*(a+1.)*(a-1.)*y2l+b*(b+1.)*(b-1.)*y2h)*(h**2)/6.
  yp = (yh-yl)/h+((3.0*b*b-1.0)*y2h-(3.0*a*a-1.0)*y2l)/6.0*h
  ypp = y2l + (y2h-y2l)/h*(xx-xa[klo])

  Klold=klo
  Nold=n
  Xa1old=xa[0]
  Xanold=xa[n]

  return y,yp,ypp

#enddef util_spline_inter(xa,ya,y2a,x,y,mode)

def util_spline_coef(x,y,yp1=9999.,ypn=9999.):

#+seq,mshimportsind.
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  # Histograms and Ntuples
  global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
  Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
  H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
  Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
  Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
  H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
  FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman
#+KEEP,plotglobind,T=PYTHON.
#*CMZ :          28/09/2019  14.39.13  by  Michael Scheer
  global MPLmain, MPLmaster, Nfigs,Figgeom, Figgeom2, FiggeomR, FiggeomL, XtermGeo, Figs,Fig,Ax,\
  Fig1,Ax1,Fig6,Ax6,Fig2,Ax2,Fig7,Ax7,Fig3,Ax3,Fig8,Ax8, Figgeoms, \
  Fig4,Ax4,Fig9,Ax9,Fig5,Ax5,Fig10,Ax10,\
  Screewidth, Screenheight, ScaleSizeX, ScaleSizeY, \
  FirstConsole, Console, Igetconsole,Klegend, Fwidth, Fheight, Fxoff, Fyoff, \
  Kfig, Kax, Ihist,Iprof, Imarker, Ierr, Isurf, Iinter, Isame, Itight, IsameGlobal, Iline, CMap, Cmap, Tcmap, Surfcolor, Cmaps, \
  Iplotopt, Ispline, Kecho, Kdump,Kpdf, Ndump,Npdf, Legend, \
  Kplots,Nwins, Zones, Kzone, Nxzone, Nyzone, Zone, Axes, Icmap, \
  Mode3d,Mode3D, Mode2d,Mode2D, CanButId, CanButIds, \
  MarkerSize, MarkerType, MarkerColor, \
  Markersize, Markertype, Markercolor, \
  Fillstyle, FillStyle, \
  Textcolor, WaveFilePrefix,WaveDump, \
  LineStyle, LineWidth, LineColor, \
  Linestyle, Linewidth, Linecolor, \
  Author, \
  Tightpad, Xtightpad,Ytightpad, ColorbarPad,\
  LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace, \
  Histcolor, Histedgecolor, Histbarwidth, Kdate, Kfit, Kstat, YTitle, \
  Icont3d, Iboxes, Inoempty, Iclosed,Itrisurf, Iscatter, Iscat3d, Ifill1d, TitPad, Xtitle, Ytitle, \
  Gtit,Xtit,Ytit,Ztit,Ttit,Ptit,Colors, Surfcolors,Linestyles, Markertypes, \
  LexpX,LexpY,LexpRot,LexpPow,\
  GtitFontSize,Titfontsize,Atitfontsize,Axislabelsize,Textfontsize,Datefontsize,\
  Statfontsize, Axislabeldist, Axislabeldist3d, Axisdist, Axisdist3d, \
  XFit, YFit, Xfit, Yfit,Ystat, YStat, \
  GtitFontSize,TitFontSize,AtitFontSize,AxisLabelSize,TextFontSize,DateFontSize,\
  StatFontSize, AxisLabelDist, AxisLabelDist3d, AxisTitleDist, AxisTitleDist3d, \
  AtitFontSize3d, Atitfontsize3d, NXtick,NXtick3d, Nxtick,Nxtick3d, Ktitles,  Dummy,\
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,ZoomZmin,ZoomZmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull2D,Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux,Ishow

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ,Tnpa,Tnone

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  """
C--- calculates spline coefficients

C--   input:

C-       x: array of x-values
C-       y: array of y-values
C-       yp1:  second derivative at first x-value
C-       ypn:  second derivative at last x-value

C--   ouput:

C-       y2:   spline-coefficients

C--   workingspace: aa(n),bb(n),cc(n),c(n)
  """
  n = len(x)

  y2 = np.zeros_like(x)

  aa = np.zeros_like(x)
  bb = np.zeros_like(x)
  cc = np.zeros_like(x)
  c  = np.zeros_like(x)

  ifail = 0

  if n < 3 :
      if abs(yp1) == 9999.0 :
          y2[0]=0.0
      else:
          y2[0]=yp1
      #endif
      if abs(ypn) == 9999.0 :
          y2[n-1]=0.0
      else:
          y2[n-1]=ypn
      #endif
      return ifail, y2
  #endif

  if abs(yp1) == 9999.0 :
      xx=x[0:3]
      yy=y[0:3]
      a, yp, opt, kfail = util_parabel(xx,yy)
      if kfail == 0 or kfail == 2:
          y2[0]=2.0*a[2]
      else:
          y2[0]=0.0
      #endif
  else:
      y2[0]=yp1
  #endif

  if abs(ypn) == 9999.0 :
      xx=x[-3:]
      yy=y[-3:]
      a, yp, opt, kfail = util_parabel(xx,yy)
      if kfail == 0 or kfail == 2:
          y2[n-1]=2.0*a[2]
      else:
          y2[n-1]=0.0
      #endif
  else:
      y2[n-1]=ypn
  #endif

  c[0]=y2[0]
  c[n-1]=y2[n-1]

  bb[0]=1.0
  cc[0]=0.0
  cc[n-1]=1.0

  j=1
  while j < n-1:
    if x[j+1] == x[j] :
      print('*** error in util_spline_coef: intervall of zero length')
      print('j, x[j], x[j+1]:',j,x[j],x[j+1])
      return -1
    #endif
    aa[j]=(x[j]-x[j-1])/6.0
    bb[j]=(x[j+1]-x[j-1])/3.0
    cc[j]=(x[j+1]-x[j])/6.0
    c[j]=(y[j+1]-y[j])/(x[j+1]-x[j])-(y[j]-y[j-1])/(x[j]-x[j-1])
    j += 1
  #enddo !j

  j = 0
  while j < n - 2:
  #do j=2,n-1

    j += 1

    bb[j]=bb[j]-aa[j]*cc[j-1]
    c[j]= c[j]-aa[j]* c[j-1]
    #aa[j]=aa[j]-aa[j]*bb[j-1]

    cc[j]=cc[j]/bb[j]
    c[j]= c[j]/bb[j]
    bb[j]=1.0

  #enddo !j

  j = n - 1
  while j > 1:
  #do j=n-1,2,-1
    j -= 1
    y2[j]=c[j]-cc[j]*y2[j+1]
    if abs(y2[j]) < 1.0e-15: y2[j]=0.0
  #enddo
  #endwhile j > 2

  return ifail, y2

#enddef util_spline_coeff(x,y,yp1=9999.,ypn=9999.)

def Quit(*args, delay=0):
  #reakpoint()
  nargs =  len(args)

  text = ''
  for i in range(nargs):
    text += str(args[i]) + " "
  #endif

  if delay > 0:

    if len(text):
      print("\n",text, "\nWaiting",delay," seconds before kill")
      #time.sleep(delay)
    else:
      print("\nWaiting",delay," seconds before kill")
      #time.sleep(delay)
    #endif len(text):

    set_console_title(os.getcwd())

    if platform.system() == 'Windows':
      stat = os.system("sleep " + str(delay) + " && taskkill /F /PID " + str(os.getpid()) + " &")
    else:
      stat = os.system("sleep " + str(delay) + " && kill " + str(os.getpid()) + " &")
    #endif platform.system() == 'Windows'

  elif delay < 0:
    return
  else:
    print("\n",text)
    set_console_title(os.getcwd())
    if platform.system() == 'Windows':
      stat = os.system("taskkill /F /PID " + str(os.getpid()))
    else:
      stat = os.system("kill " + str(os.getpid()))
    #endif platform.system() == 'Windows'

#enddef Quit(text = '', delay=0)

def exit(text = ''): Quit(text)

def Wexit(ew=''):
  if type(ew) == str:
    Quit(ew)
  else:
    Quit()

def qwait(delay=3):
  Quit("",delay)

def sleep(isec):
  time.sleep(isec)

def wait(isec=1000,text='waiting'):
  if text: print("--- ",text," for ",isec," seconds ---")
  time.sleep(isec)

def util_parabel(xin,yin):

  """
  C--- CALCULATES A(1),A(2),A(3), THE DERIVATIVES YP(X(1)),YP(X(2)),YP(X(3)),
  C    AND THE EXTREMUM (XOPT,A(XOPT)) OF PARABOLA Y=A1+A2*X+A3*X**2
  C    FROM COORDINATES OF THE THREE POINTS (X(1),Y(1)),(X(2),Y(2)),(X(3),Y(3))
  C
  """

  #reakpoint()
  a = [-9999.,-9999.,-9999.]
  yp = [-9999.,-9999.,-9999.]
  opt = [-9999.,-9999.]
  ifail = 0

# calculate f=a0+a1*(x-x0)+a2*(x-x0)**2
#  = a0 + a1*x - a0*x0 + a2*x**2 - 2*a2*x*x0 + a2*x0**2
#  = a0 + (a2*x0 -a0)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

# change system: (x0,s0)->(0,0), i.e.
# calculate f=a1*dx+a2*dx**2
#  df/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

  #print("util_parabel:",xin,"\n",yin)
  xy = pd.DataFrame(columns=['x','y'])
  xy.x = xin
  xy.y = yin
  xy = xy.sort_values(by='x')

  x = xy.x
  y = xy.y

  x.index = range(len(x))
  y.index = range(len(y))

  x0=x[1]
  f0=y[1]

  fm=y[0]-f0
  fp=y[2]-f0

  dxm=x[0]-x0
  dxp=x[2]-x0

  # fm=a1*dxm+a2*dxm**2
  # fp=a1*dxp+a2*dxp**2

  # (dxm dxm2) (a1) = (y[0])
  # (dxp dxp2) (a2) = (y[2])

  dxm2=dxm*dxm
  dxp2=dxp*dxp

  det=dxm*dxp2-dxp*dxm2

  if det != 0.0:
    a1=(fm*dxp2-fp*dxm2)/det
    a2=(fp*dxm-fm*dxp)/det
  else:
    ifail=1
    return a,yp,opt,ifail
  #endif

  if a2 != 0.0:
    dxmax=-a1/(2.0*a2)
    dymax=(a1+a2*dxmax)*dxmax
    opt[0]=x0+dxmax
    opt[1]=f0+dymax
  #endif

  # calculate f=f0+a1*dx+a2*dx**2
  # = a1*x - a1*x0 + a2*x**2 + a2*x0**2 - 2*a2*x*x0
  #  f = f0 + (a2*x0 -a1)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

  a22=2.0*a2

  a[0]=f0 + (a2*x0 -a1)*x0
  a[1]=a1 - a22*x0
  a[2]=a2

  # calculate yp=a1+2*a2*dx

  yp[0]=a1+a22*dxm
  yp[1]=a1
  yp[2]=a1+a22*dxp

  if opt[0] <= min(xin) or opt[0] >= max(xin): ifail = 2

  return a,yp,opt,ifail
#enddef util_parabel(xin,yin)
def util_eqn(a,b):

  x = []
  istat = 0

  npar = len(b)
  nrow = len(a)

  if nrow < 1:
    istat = -1
    return istat,x
  #endif nrow < 1

  ncolumn = len(a[0])
  if nrow != ncolumn:
    istat = -2
    return istat,x
  #endif nrow < 1

  if nrow != npar:
    istat = -3
    return istat,x
  #endif nrow < 1

  x = np.linalg.solve(a,b)
  return istat,x

#enddef util_eqn(a,b)

def util_linear_fit_data(datapoints,fitfun):

  istat = 0
  param = np.array(1)

  nffun = len(fitfun[0][0])
  nfpoi = len(fitfun)
  npar = len(fitfun[0])
  nfun = len(datapoints[0])
  npoi = len(datapoints)

 # print(npar,nfun,npoi)
 # print(fitfun)
 # print(len(fitfun))

  a=np.zeros([npar,npar], dtype = 'float')
  param=np.zeros([npar], dtype = 'float')

  if len(fitfun) != npoi:
    print("*** Error in util_linear_fit_data: Number of functions differ for datapoints and fit-functions")
    istat = 1
    return istat, param
  #endif

  if nffun != nfun:
    print("*** Error in util_linear_fit_data: Number of functions differ for datapoints and fit-functions")
    istat = 1
    return istat, param
  #endif

#  print(datapoints)

  for ipar in range(npar):
    for ifun in range(nfun):
      for ipoi in range(npoi):
        param[ipar] += datapoints[ipoi][ifun] * fitfun[ipoi][ipar][ifun]
        #print(ipoi,ipar,ifun, datapoints[ipoi][ifun],fitfun[ipoi][ipar][ifun])
      #endfor ipoi in range(npoi)
    #endfor ifun in range(nfun)
  #endfor ipar in range(npar)

  for jpar in range(npar):
    for ipar in range(npar):
      for ifun in range(nfun):
        for ipoi in range(npoi):
          a[ipar][jpar] += fitfun[ipoi][ipar][ifun] * fitfun[ipoi][jpar][ifun]
        #endfor ipoi in range(npoi)
      #endfor ifun in range(nfun)
    #endfor ipar in range(npar)
  #endfor jpar in range(npar)

  param = np.linalg.solve(a,param)

  return istat,param

#endif

def util_vnorm(v):
  vn = 0
  for i in range(len(v)): vn += v[i]*v[i]
  return np.sqrt(vn)
#enddef util_vnorm

def util_rotate(cen,vrot,phi,vin,eps=1.0e-10):

      istat=0
      vlen=util_vnorm(vrot)

      if vlen == 0.0:
        vout=vin
        istat=1
        return istat, vout
      #endif

      o = vrot/vlen

      s = np.sin(phi)
      c = np.cos(phi)

      c1 = 1.0 - c

      rm11 = o[0] * o[0] * c1 + c
      rm22 = o[1] * o[1] * c1 + c
      rm33 = o[2] * o[2] * c1 + c

      rm12 = o[0] * o[1] * c1 - o[2] * s
      rm13 = o[0] * o[2] * c1 + o[1] * s

      rm21 = o[0] * o[1] * c1 + o[2] * s
      rm23 = o[1] * o[2] * c1 - o[0] * s

      rm31 = o[0] * o[2] * c1 - o[1] * s
      rm32 = o[1] * o[2] * c1 + o[0] * s

      r = np.array(vin) - np.array(cen)

      vout1 = rm11 * r[0] + rm12 * r[1] + rm13 * r[2] + cen[0]
      vout2 = rm21 * r[0] + rm22 * r[1] + rm23 * r[2] + cen[1]
      vout3 = rm31 * r[0] + rm32 * r[1] + rm33 * r[2] + cen[2]

      rm = [[rm11,rm12,rm13],[rm21,rm22,rm23],[rm31,rm32,rm33]]
      vout = [vout1,vout2,vout3]

      if eps:
        for i in range(3):
          if abs(vout[i]) < eps: vout[i] = 0.0
          for j in range(3):
            if abs(rm[i][j]) < eps: rm[i][j] = 0.0
          #endfor
        #endfor
      #endif

      return istat,vout,rm
#def util_rotate(cen,vrot,phi,vin,vout,istat)

def evnm(eg): return wtoe1/eg
def nmev(wl): return wtoe1/wl

def fsize(cfile):
  return os.stat(cfile).st_size
#endif

def fexist(f):
  import os
  if os.path.exists(f): return True
  return False
#enddef

def fwrite(F,*args):

  nargs =  len(args)

  text = ''
  for i in range(nargs-1):
    text += str(args[i]) + " "
  #endif
  text += str(args[nargs-1]) + "\n"
  F.write(text)
#enddef

def util_determinante(a):
  return np.linalg.det(a)
#enddef

def util_solve(a,x):
  return np.linalg.solve(a,x)
#enddef

import os,sys,platform,shutil,time,re

import tkinter as tk
from tkinter import *

import numpy as np
import pandas as pd

from scipy import special

import matplotlib as mpl
import matplotlib.pyplot as plt

global ebeam,curr,b,elow,ehigh,nener,bw,dist,psimax,psimin,npsi,g,ec,psi,ener,\
eec,dip0,dip,absmuen_default,absmufe_default,abupb_default,\
thickbe,thickfe,thickpb,absmuen,absmube,absmufe,absmupl,powsum

absmuen = []
absmube = []
absmufe = []
absmupb = []

thickbe= 0.0
thickfe= 0.0
thickpb= 0.0

densebe = 1845.0 # kg/m**3
densefe = 7870.  # kg/m**3
densepb = 11350. # kg/m**3

#SIGMA(tot_a) Z=8; NUCL. DATA TABLES A7,576 (1970); DENSITY OF WATER USED !!!
#1000.	!UNITS: KG AND METER	; DENSITY OF WATER USED !!!
absmuen_default = [\
  [1000.000,455.4440], \
  [1500.000,153.5712], \
  [2000.000,68.88120], \
  [3000.000,21.53008], \
  [4000.000,9.221800], \
  [5000.000,4.705000], \
  [6000.000,2.706316], \
  [8000.000,1.091560], \
  [10000.00,0.5344880], \
  [15000.00,0.1486780], \
  [20000.00,0.5947120E-01], \
  [30000.00,0.1674980E-01], \
  [40000.00,0.7302160E-02], \
  [50000.00,0.4290960E-02], \
  [60000.00,0.3131648E-02], \
  [80000.00,0.2416488E-02], \
  [100000.0,0.2314860E-02], \
  [150000.0,0.2484240E-02], \
  [200000.0,0.2657384E-02], \
  [300000.0,0.2871932E-02]]

#SIGMA(tot) Z=4; NUCL. DATA TABLES A7,575 (1970)
#1845.   !UNITS: KG AND METER
#  28
absmube_default = [\
[200.0000,4959.528], \
[300.0000,1792.058], \
[400.0000,841.5970], \
[500.0000,450.0480], \
[600.0000,268.5272], \
[700.0000,171.6888], \
[800.0000,116.0137], \
[900.0000,81.99165], \
[1000.000,59.61236], \
[1500.000,17.84361], \
[2000.000,7.351300], \
[3000.000,2.018266], \
[4000.000,0.8086430], \
[5000.000,0.4023166], \
[6000.000,0.2339050], \
[8000.000,0.1055914], \
[10000.00,0.6148360E-01], \
[15000.00,0.3007350E-01], \
[20000.00,0.2232122E-01], \
[30000.00,0.1791044E-01], \
[40000.00,0.1644018E-01], \
[50000.00,0.1557139E-01], \
[60000.00,0.1496992E-01], \
[80000.00,0.1403430E-01], \
[100000.0,0.1329917E-01], \
[150000.0,0.1189574E-01], \
[200000.0,0.1089329E-01], \
[300000.0,0.9423030E-02]]

#  TOT-H of FE; NUCL. DATA TABL. A7, 586 (1970)
#  7870.	!UNITS: KG AND METER
#  22
absmufe_default = [ \
[1000.000,887.1940], \
[1500.000,329.8680], \
[2000.000,158.4660], \
[3000.000,54.97800], \
[4000.000,24.79400], \
[5000.000,13.47500], \
[6000.000,8.095779], \
[7112.000,5.098940], \
[7112.000,41.17960], \
[8000.000,30.39960], \
[10000.00,16.92460], \
[15000.00,5.594820], \
[20000.00,2.490180], \
[30000.00, 0.7675360], \
[40000.00, 0.3287900], \
[50000.00, 0.1681680], \
[60000.00, 0.9820579E-01], \
[80000.00, 0.4204200E-01], \
[100000.0, 0.2209900E-01], \
[150000.0, 0.8031099E-02], \
[200000.0, 0.4861780E-02], \
[300000.0, 0.3384920E-02]]

#  TOT of PB; NUCL. DATA TABL. A7, 624 (1970)
#  11350.	!UNITS: KG AND METER (columns: eV and m**2/kg)
#  38
absmupb_default = [
[1000.000,511.6320], \
[1500.000,234.0135], \
[2000.000,127.0359], \
[2484.000,79.36110], \
[2484.000,210.4668], \
[2586.000,190.9899], \
[2586.000,277.3278], \
[3000.000,196.2225], \
[3066.000,184.3038], \
[3066.000,213.3738], \
[3554.000,148.8384], \
[3554.000,157.8501], \
[3851.000,130.2336], \
[3851.000,136.0476], \
[4000.000,124.1289], \
[5000.000,72.09360], \
[6000.000,46.22130], \
[8000.000,22.41297], \
[10000.00,12.76173], \
[13035.00,6.598890], \
[13035.00,16.27920], \
[15000.00,11.04660], \
[15200.00,10.75590], \
[15200.00,14.88384], \
[15861.00,13.43034], \
[15861.00,15.49431], \
[20000.00,8.517510], \
[30000.00,2.994210], \
[40000.00,1.418616], \
[50000.00, 0.7877970], \
[60000.00, 0.4941900], \
[80000.00, 0.2369205], \
[88004.00, 0.1875015], \
[88004.00, 0.7558200], \
[100000.0, 0.5523300], \
[150000.0, 0.2011644], \
[200000.0, 0.9941940E-01], \
[300000.0 ,0.3982590E-01]]

def H2(y): return (y * special.kv(2./3.,complex(y/2.)).real)**2

def G1a(y):
    # from Shaukart Khan
    if y < 4.:
        spectrum = 391.8 * y**0.333 * np.exp(-y*0.8307) \
        -192.0 * y**0.500 * np.exp(-y*0.7880)
    else:
        spectrum = 164.0 * y**0.500 * np.exp(-y)
    #endif
    return spectrum*0.007565
#enddef

def dfdtdp(y,psi,gamma,cur,banwid):

# Output:
# Hori., vert. component of the spectral flux-density of a dipole
# in dN/dt/BW/dA in 1/s/banwid/m**2, and the power density in W/m**2
# YES, m**2!

# Input:
# y: Egam / Echar with the photon energy Egam in eV and the characteristic
# energy of the dipole
# psi: The vertical observation angle in rad
# gamma: Rel. gamma factor of the beam electron
# cur: The beam current in Ampere
# banwid: BW = dEgam/Egam * banwid

    const=3.0/4.0*alpha1/pi1**2*gamma**2*banwid*cur/echarge1

    x=gamma*psi
    xx=x*x
    xx1=xx+1.0
    xi=y*np.sqrt(xx1)**3/2.0

    bk13=special.kv(1./3.,complex(xi)).real
    bk23=special.kv(2./3.,complex(xi)).real

    par=const*(y*xx1)**2*bk23**2
    per=const*y*y*xx*xx1*bk13**2

    dfdtdp=par+per

    powr=7.0/16.0*echarge1*echarge1/4.0/pi1/eps01 \
      /(np.sqrt(psi*psi+1.0/gamma/gamma))**5 \
      *(1.0+5.0/7.0*xx/xx1)*cur/echarge1

    return par,per,powr
#enddef

def dip_readin():

    global ebeam,curr,b,elow,ehigh,nener,bw,dist,psimax,psimin,npsi,g,ec,psi,\
    ener,eec,dip,dip0,thickbe,thickfe,thickpb,idose

    try:

        Fin = open('pyDipole.in','r')
        cin = Fin.readlines()
        Fin.close()

        ebeam = float(cin[0].split(':')[1])
        curr = float(cin[1].split(':')[1])
        b = float(cin[2].split(':')[1])
        elow = float(cin[3].split(':')[1])
        ehigh = float(cin[4].split(':')[1])
        nener = int(cin[5].split(':')[1])
        bw = float(cin[6].split(':')[1])
        dist = float(cin[7].split(':')[1])
        psimin = float(cin[8].split(':')[1])
        psimax = float(cin[9].split(':')[1])
        npsi = int(cin[10].split(':')[1])

        thickbe = float(cin[11].split(':')[1])
        thickfe = float(cin[12].split(':')[1])
        thickpb = float(cin[13].split(':')[1])

        idose = int(cin[14].split(':')[1])

    except:
    #else:

        ebeam = readfloat("\nBeam energy [GeV]:")
        curr = readfloat("Beam current [Amp]:")

        b = readfloat("Field strength [T]:")

        elow = readfloat("Lowest photon energy [keV]:")
        ehigh = readfloat("Highest photon energy [keV]:")
        nener = readint("Number of energies:")

        bw = readfloat("Bandwidth:")

        dist = readfloat("Distance to source [m]:")

        psimin = readfloat("Min vertical angle [mrad]:")
        psimax = readfloat("Max vertical angle [mrad]:")
        npsi = readint("Number of angles:")

        thickbe = readint("Thickness of Be filter [mm]:")
        thickfe = readint("Thickness of Fe filter [mm]:")
        thickpb = readint("Thickness of Pb filter [mm]:")

        idose = 0

    #endtry

    if ebeam <= 0: ebeam = 1.7
    g = ebeam / emassg1

    if curr <= 0: curr = 0.3
    if dist <= 0.0: dist = 10.

    if b <= 0.0: b = 1.3

    if nener <= 0: nener = 100
    if elow <= 0.0: elow = 0.1
    if ehigh < elow: ehigh = 100. * elow
    if bw <= 0: bw = 0.001

    npsi = int(npsi / 2) * 2 + 1
    if npsi <= 0: npsi = 101
    if psimax <= 0.: psimax = 5. / g * 1000. # mrad

    ec = 0.665 * b * ebeam**2

    ener = np.linspace(elow,ehigh,nener)
    psi = np.linspace(psimin,psimax,npsi)
    eec = ener / ec

#enddef dip_readin()

def dip_printin():

    global ebeam,curr,b,elow,ehigh,nener,bw,dist,psimax,psimin,npsi,g,ec,psi,\
    ener,eec,dip,dip0,thickbe,thickfe,thickpb,idose

    Fout = open("pyDipole.in",'w')

    print("\n\nBeam energy [GeV]:",ebeam)
    print("Beam current [Amp]:",curr)

    print("Field strength [T]:",b)

    print("Characteristic energy [keV]:",ec)

    print("Gamma:",g)

    print("Lowest photon energy [keV]:",elow)
    print("Highest photon energy [keV]:",ehigh)
    print("Number of energies:",nener)
    print("BW:",bw)

    print("Distance to source [m]:",dist)

    print("Min. vertical angle [mrad]:",psimin)
    print("Max. vertical angle [mrad]:",psimax)
    print("Number of angles:",npsi)

    print("Thickness of Be filter [mm]:",thickbe)
    print("Thickness of Fe filter [mm]:",thickfe)
    print("Thickness of Pb filter [mm]:",thickpb)

    Fout.write("Beam energy [GeV]: " + str(ebeam) + '\n')
    Fout.write("Beam current [Amp]: " + str(curr) + '\n')

    Fout.write("Field strength [T]: " + str(b) + '\n')

    Fout.write("Lowest photon energy [keV]: " + str(elow) + '\n')
    Fout.write("Highest photon energy [keV]: " + str(ehigh) + '\n')
    Fout.write("Number of energies: " + str(nener) + '\n')
    Fout.write("BW: " + str(bw) + '\n')

    Fout.write("Distance to source [m]: " + str(dist) + '\n')

    Fout.write("Min. vertical angle [mrad]: " + str(psimin) + '\n')
    Fout.write("Max. vertical angle [mrad]: " + str(psimax) + '\n')
    Fout.write("Number of angles: " + str(npsi) + '\n')

    Fout.write("Thickness of Be filter [mm]: " + str(thickbe) + '\n')
    Fout.write("Thickness of Fe filter [mm]: " + str(thickfe) + '\n')
    Fout.write("Thickness of Pb filter [mm]: " + str(thickpb) + '\n')

    Fout.write("DoCalcDose: " + str(idose) + '\n')

    Fout.close()
    #reakpoint()
#endef dip_printin()

def dip_calc_and_writeout():

    global ebeam,curr,b,elow,ehigh,nener,bw,dist,psimax,psimin,npsi,g,ec,psi,\
    ener,eec,dip,dip0,powsum

    Fspec = open("pyDipole_unfiltered.out","w")
    Fspec0 = open("pyDipole_0_unfiltered.out","w")
    Fspec0.write("* E=" + str(ebeam) + " GeV, Current=" + str(curr) + "  A , Ec=" + str(ec) + " keV, BW=" + str(bw) + " \n")
    Fspec.write("* E=" + str(ebeam) + " GeV, Current=" + str(curr) + "  A , Ec=" + str(ec) + " keV, BW=" + str(bw) + " \n")
    Fspec.write("* y [mm] Egamma [keV] Flux-density [dN/s/BW/mm**2]  Power-density [W/mm**2] \n")
    Fspec0.write("* Egamma [keV] Flux-density [dN/s/BW/mm**2] Flux [dN/s/mrad/BW] H2 G1 Power-density [W/mm**2] FluxDen/Flux [mrad]\n")

    dist2 = dist**2

    powsum = 0.0

    ia = 0
    for ang in psi:
        ia += 1
        for ee in ener:
          y = ee / ec
          par,per,powr = dfdtdp(y,ang/1000.,g,curr,bw)
          if ia > 1 and ia <= npsi:
            powsum += powr
          else:
            powsum += powr/2.0
          #endif
          spec = (par + per) / 1.e6 / dist2
          if ang == 0:
            h2 = H2(y)
            g1 = G1a(y)
            flux = 2.457e13*ebeam*curr*g1
            Fspec0.write( \
            '{:.4g}'.format(ee) + " " \
            + '{:.4g}'.format(spec) + " " \
            + '{:.4g}'.format(flux) + " " \
            + '{:.4g}'.format(h2) + " " \
            + '{:.4g}'.format(g1) + "  " \
            + '{:.4g}'.format(powr/1.e6/dist2) + " " \
            + '{:.4g}'.format(flux/par*1.e6) + "\n" \
            )
          #endif
          Fspec.write('{:.4g}'.format(ang*dist) + " " + '{:.4g}'.format(ee) \
          + " " + " " + '{:.4g}'.format(spec) + " " \
          + '{:.4g}'.format(powr/1.e6/dist2) + " " + "\n")
        #endfor ener
    #endfor ang

    Fspec.close()
    Fspec0.close()

    powsum = powsum/1.0e6*(psi[1]-psi[0])/1000.
    print("\nUnfiltered total power:",'{:.4g}'.format(powsum),' [W/mrad]\n')

#enddef dip_writeout()

def dip_writeout():
    #reakpoint()
    global ebeam,curr,b,elow,ehigh,nener,bw,dist,psimax,psimin,npsi,g,ec,psi,\
    ener,eec,dip,dip0

    Fspec0 = open("pyDipole_0.out","w")
    Fspec0.write("* E=" + str(ebeam) + " GeV, Current=" + str(curr) + "  A , Ec=" + str(ec) + " keV, BW=" + str(bw) + " \n")
    Fspec0.write("* Egamma [keV] Flux-density [dN/s/BW/mm**2] Flux [dN/s/mrad/BW] H2 G1 Power-density [W/mm**2] FluxDen/Flux [mrad]\n")
    Fspec0.close()

    dip0.to_csv( \
    'pyDipole_0.out',header=False,index=False,sep=' ',decimal='.',mode='a')

    Fspec = open("pyDipole.out","w")
    Fspec.write("* E=" + str(ebeam) + " GeV, Current=" + str(curr) + "  A , Ec=" + str(ec) + " keV, BW=" + str(bw) + " \n")
    Fspec.write("* y [mm] Egamma [keV] Flux-density [dN/s/BW/mm**2]  Power-density [W/mm**2] \n")
    Fspec.close()

    dip.to_csv( \
    'pyDipole.out',header=False,index=False,sep=' ',decimal='.',mode='a')

#enddef dip_writeout()

def dip_readout():

    global dip0,dip,nener,ener

    #reakpoint()
    fname = 'pyDipole_0_unfiltered.out'
    try:
        dip0 = pd.read_csv(fname,header=None,delim_whitespace=True,comment='*')
        pass
    except:
        print("\n*** Error nread(...): Failed to read " + fname)
    #endif

    dip0.columns = ['Egam','FluxDen','Flux','H2','G1','PowDen','vSiz']

    fname = 'pyDipole_unfiltered.out'
    try:
        dip = pd.read_csv(fname,delim_whitespace=True,header=None, comment='*')
    except:
        print("\n*** Error nread(...): Failed to read " + fname)
    #endif
    dip.columns = ['hPos','Egam','FluxDen','PowDen']

    ener = np.array(dip0.Egam)
    nener = len(ener)

#enddef dip_readout()

def dip_plot():

    global dip0,dip,fddose,flxdose,thickfe,thickpb
    global absmuen,absmube,absmufe,absmupb,\
    absmuen_default,absmube_default,absmufe_default,absmupb_default

    #reakpoint()

    dip_readout()

    emin = dip0.Egam.min()
    emax = dip0.Egam.max()
    nener = len(dip0.Egam)
    dener =  (emax - emin) / max(nener-1,1)

    thresh = 1.0e-20

    Fig = plt.gcf()
    Ax = plt.gca()

    #d0pl = dip0.query('FluxDen > ' + str(thresh))
    #plt.plot(d0pl.Egam,d0pl.FluxDen)
    plt.plot(dip0.Egam,dip0.FluxDen)

    #Ax.set_yscale('log')

    Ax.set_title('Flux-Density')
    Ax.set_xlabel('E$_{\gamma}$ [keV]')
    Ax.set_ylabel('flux-density [N$_\gamma$/mm$^{2}$/s/BW]')
    titx = 0.7
    tity = 0.85

    if fddose:
        Ax.text(titx,tity,'Dose Rate: \n'+'{:.3g}'.format(fddose)+' [mGy/h]',transform=Ax.transAxes)

    plt.savefig("dipole_flux_density.pdf")
    print("dipole_flux_density.pdf written")
    plt.show(block=False)

    plt.figure()
    Fig = plt.gcf()
    Ax = plt.gca()

    #d0pl = dip0.query('Flux > ' + str(thresh))
    #plt.plot(d0pl.Egam,d0pl.Flux)
    plt.plot(dip0.Egam,dip0.Flux)

    #Ax.set_yscale('log')

    Ax.set_title('Flux')
    Ax.set_xlabel('E$_{\gamma}$ [keV]')
    Ax.set_ylabel('flux [N$_\gamma$/s/BW]')
    if fluxdose:
        Ax.text(titx,tity,'Mean Dose Rate: \n'+'{:.3g}'.format(fluxdose)+' [mGy/h]',transform=Ax.transAxes)
        Ax.text(titx,tity,'Mean Dose Rate: \n'+'{:.3g}'.format(fluxdose)+' [mGy/h]',transform=Ax.transAxes)
    #endif
    plt.savefig("dipole_flux.pdf")
    print("dipole_flux.pdf written")
    plt.show(block=False)

    dpl = dip.query('Egam == ' + str(emin))
    if dpl.FluxDen.max() > thresh:

      plt.figure()
      Fig = plt.gcf()
      Ax = plt.gca()

      plt.plot(dpl.hPos,dpl.FluxDen)
      Ax.set_title('Vert. Flux-Density Distribution for E = ' + str(emin) + '[keV]')
      Ax.set_xlabel('vert. pos. [mm]')
      Ax.set_ylabel('flux-density [N$_\gamma$/mm$^{2}$/s/BW]')

      y = np.array(dpl.hPos)
      fd = np.array(dpl.FluxDen)

      w = 0.0
      ysum2 = 0.0
      wsum = 0.0
      for i in range(len(dpl.hPos)):
        w = fd[i]
        ysum2 += y[i]**2 * w
        wsum += w
      #endfor

      if wsum >0: rms = np.sqrt(ysum2/wsum)
      else: rms = 0

      titx = 0.7

      Ax.text(titx,tity,'RMS: '+'{:.3g}'.format(rms)+' [mm]',transform=Ax.transAxes)

      plt.savefig("dipole_flux_density_vert_dist_Emin.pdf")
      print("dipole_flux_density_vert_dist_Emin.pdf written")
      plt.show(block=False)

    #endif

    #reakpoint()
    dpl = dip.query('Egam == ' + str(emax))

    if dpl.FluxDen.max() > thresh:

      plt.figure()
      Fig = plt.gcf()
      Ax = plt.gca()

      plt.plot(dpl.hPos,dpl.FluxDen)
      Ax.set_title('Vert. Flux-Density Distribution for E = ' + str(emax) + '[keV]')
      Ax.set_xlabel('vert. pos. [mm]')
      Ax.set_ylabel('flux-density [N$_\gamma$/mm$^{2}$/s/BW]')

      y = np.array(dpl.hPos)
      fd = np.array(dpl.FluxDen)

      w = 0.0
      ysum2 = 0.0
      wsum = 0.0
      for i in range(len(dpl.hPos)):
        w = fd[i]
        ysum2 += y[i]**2 * w
        wsum += w
      #endfor

      if wsum > 0: rms = np.sqrt(ysum2/wsum)
      else : rms = 0

      Ax.text(titx,tity,'RMS: '+'{:.3g}'.format(rms)+' [mm]',transform=Ax.transAxes)

      plt.savefig("dipole_flux_density_vert_dist_Emax.pdf")
      print("dipole_flux_density_vert_dist_Emax.pdf written")

      plt.show(block=False)
    #endif

    plt.figure()
    Fig = plt.gcf()
    Ax = plt.gca()

    plt.plot(dip0.Egam,dip0.vSiz)

    Ax.set_title('Vertical Beam Size')
    Ax.set_xlabel('E$_\gamma$ [keV]')
    Ax.set_ylabel('Flux-density / Flux [mrad]')

    plt.savefig("dipole_vert_beam_size.pdf")
    print("dipole_vert_beam_size.pdf written")

    plt.show(block=False)

#enddef dip_plot()

fddose = 0.0
fluxdose = 0.0

def dip_dose(fabs = 'ABSORPDOSE.RP'):

  global dip,dip0, bw, psimax, psimin, dist, fddose,fluxdose

  #reakpoint()

  fdd = np.array(dip0.FluxDen)
  flxdd = np.array(dip0.Flux)

  for iene in range(nener):
    absc = dip_absmuint(ener[iene],'Tissue')
    if absc == 0.0:
      print('*** Error in dip_dose(): Egam ',ener[iene],'keV is out of table ***')
      print("*** Contribution will be ignored !! ***")
    #endif
    fdd[iene] = fdd[iene] * absc
    flxdd[iene] = flxdd[iene] * absc
  #endfor

  for iene in range(nener-1):
    de = (ener[iene+1]-ener[iene]) * 1000.0
    fddose += (fdd[iene]+fdd[iene+1])/2. * de
    fluxdose += (flxdd[iene]+flxdd[iene+1])/2. * de
  #endfor

  fddose = fddose * 1.0e6 / bw * echarge1 * 3600. * 1000.
  #fluxdose = fluxdose * 1.0e6 / bw * echarge1 / dist / dist * 3600.0 * 1000.
  fluxdose = fluxdose * 1.0e6 / bw * echarge1 * 3600.0 * 1000.
  #Vorsicht, dieser Dosisbegriff ist sinnlos, da er vom hori. Winkelintervall
  #abhängt.

  vsiz_eff = fluxdose/(fddose*dist**2)*dist

  print("Flux-density Dose Rate [mGy/h]:",'{:.3g}'.format(fddose))
  print("Effective vert. aperture [mm]:",'{:.3g}'.format(vsiz_eff),'\n')
  #print("Flux Dose [mGy/h/mrad]:",'{:.3g}'.format(fluxdose))

#enddef dip_dose()

def dip_filter(fil='Fe'):

  global dip,dip0, npsi,thickfe,thickbe,thickpb,densefe,densebe,densepb

  ener = np.array(dip0.Egam)
  nener = len(ener)

  if fil == 'Be':
    thick = thickbe/1000.
    dense = densefe
  elif fil == 'Fe':
    thick = thickfe/1000.
    dense = densefe
  elif fil == 'Pb':
    thick = thickpb/1000.
    dense = densepb
  else:
    print('*** Error in dip_filter: Unknown material: ',fil)
    return
  #endif

  #reakpoint()
  #print(dip.FluxDen.max())

  for iene in range(nener):
    try:
        absc = dip_absmuint(ener[iene],fil,thick)
    except:
        print("*** Error in dip_absmuint for material ",fil,' and Egam =',ener[iene])
    #endtry
    absor =  np.exp(-thick*dense*absc)
    dip0.FluxDen[iene] = dip0.FluxDen[iene] * absor
    dip0.Flux[iene] = dip0.Flux[iene] * absor
    for ipsi in range(npsi):
        k = iene + nener*ipsi
        dip.FluxDen[k] = dip.FluxDen[k] * absor
    #endfor
  #endfor
  #print(dip.FluxDen.max())

#enddef dip_filter()

global klold,egold
klold = 0
egold = -1.0

def dip_absmuint(eg, fil='Fe', thick=0.0):

  global nener,ener,absmuen,absmube,absmufe,absmupb,klold, \
  absmuen_default,absmube_default,absmufe_default,absmupb_default

  global klold,egold

  if fil != 'Tissue' and thick <= 0.0: return 1.0

  if len(absmuen) == 0:
    try:
        absmuen = pd.read_csv('ABSORPDOSE.RP',header=None,delim_whitespace=True,comment='*',skiprows=3)
    except:
        absmuen = pd.DataFrame(absmuen_default,columns=['egam','mue'])
    #endtry
    absmuen.columns = ['egam','mue']
    absmuen.egam /= 1000.0
  #endif

  if len(absmube) == 0:
    try:
        absmube = pd.read_csv('sigma_tot.be',header=None,delim_whitespace=True,comment='*',skiprows=3)
    except:
        absmube = pd.DataFrame(absmube_default,columns=['egam','mue'])
    #endtry
    absmube.columns = ['egam','mu']
    absmube.egam /= 1000.0
  #endif

  if len(absmufe) == 0:
    try:
        absmufe = pd.read_csv('sigma_toth.fe',header=None,delim_whitespace=True,comment='*',skiprows=3)
    except:
        absmufe = pd.DataFrame(absmufe_default,columns=['egam','mue'])
    #endtry
    absmufe.columns = ['egam','mu']
    absmufe.egam /= 1000.0
  #endif

  if len(absmupb) == 0:
    try:
        absmupb = pd.read_csv('sigma_tot.pb',header=None,delim_whitespace=True,comment='*',skiprows=3)
    except:
        absmupb = pd.DataFrame(absmupb_default,columns=['egam','mue'])
    #endtry
    absmupb.columns = ['egam','mu']
    absmupb.egam /= 1000.0
  #endif

  #reakpoint()
  if fil == 'Be':
    egam = absmube.egam
    a = absmube.mu
  elif fil == 'Fe':
    egam = absmufe.egam
    a = absmufe.mu
  elif fil == 'Pb':
    egam = absmupb.egam
    a = absmupb.mu
  elif fil == 'Tissue':
    egam = absmuen.egam
    a = absmuen.mue
  else:
   print('*** Error in dip_absmuint: Unknown matrial',fil)
   return 0.0
  #endif

  if eg < egam[0] or eg > egam[len(egam)-1]:
    print('*** Error in dip_absmuint: Egam out of table material: ',fil)
    return 0.0
  #endif

  if eg >= egold:klo = klold
  else: klo = 0

  if eg >= egam[klo] and eg < egam[klo+1]:
    khi = klo +1
  else:
    khi = len(egam) - 1
    while khi-klo > 1:
      k=int((khi+klo)/2)
      if egam[k] > eg:
          khi=k
      else:
          klo=k
      #endif
    #end while
  #endif

  h=egam[khi]-egam[klo]

# interpolation of y by y=aa*x**bb

  if h != 0.0:
      bb=np.log(a[khi]/a[klo]) / np.log(egam[khi]/egam[klo])
      res=a[klo]*np.exp(bb*(np.log(eg/egam[klo])))
  else:
      res=a[klo]
  #endif

  klold = klo
  egold = eg

  return res

#enddef absmuint():

#reakpoint()

dip_readin()
dip_printin()
dip_calc_and_writeout()

try:
    dip_readout()
    try:
      if thickbe > 0: dip_filter('Be')
      if thickfe > 0: dip_filter('Fe')
      if thickpb > 0: dip_filter('Pb')
    except:
      print('*** Errors occured in dip_filter() ***')
    #endtry
    if idose:
        try:
            dip_dose()
        except:
            print('*** Errors occured in dip_dose() ***')
        #endtry
        try:
            dip_plot()
        except:
            print('*** Errors occured in dip_plot() ***')
        #endtry
    #endif
except:
    print('*** Errors occured in dip_readout() ***')
#endtry

#reakpoint()

dip_writeout()
