
global WavesMode
WavesMode = 'M_Hbook'

#!/usr/bin/env python

# +PATCH,//WAVES/PYTHON
# +KEEP,imports,T=PYTHON.

import subprocess, glob, os, sys, time, re, fileinput, shutil, platform, copy

import matplotlib as mpl
mpl.use('TkAgg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm #color maps

from tkinter import *
from tkinter.font import Font
import tkinter as tk
from tkinter import ttk

from math import *

import numpy as np
import pandas as pd

from scipy import stats
from scipy import special
from scipy import interpolate
from scipy.stats import *
from scipy.interpolate import interp1d
from scipy.spatial import ConvexHull, convex_hull_plot_2d

from copy import *
###################################################

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

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

  ypp =  zeros_like(x)

  aa = zeros_like(x)
  bb = zeros_like(x)
  cc = zeros_like(x)
  c  = zeros_like(x)
  cn = zeros_like(x)

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

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

  y2 = zeros_like(x)

  aa = zeros_like(x)
  bb = zeros_like(x)
  cc = zeros_like(x)
  c  = zeros_like(x)

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

#+seq,waveglobal.
#+seq,hisglobal.
# {
# +PATCH,//WAVES/PYTHON
# +KEEP,mhbglobal,T=PYTHON.

#+KEEP,vecglobal,T=PYTHON.
#*CMZ :          28/09/2019  14.39.13  by  Michael Scheer

global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
VxyzX,VxyzY,VxyzZ

VxyzX=None
VxyzY=None
VxyzZ=None

VsplX=None
VsplY=None
Vspl1=None
Vspl2=None
VsplI=None
VsplI=None
VsplCoef=None
Vxint=None
Vyint=None
SplineMode = 'new'

Nspline = pd.DataFrame(columns=['x','y','yp','ypp','yint'])
Ninter = pd.DataFrame(columns=['x','y','yp','ypp','yint'])
Nfitxy = pd.DataFrame(columns=['x','y','ey','fit'])
Nfitint = pd.DataFrame(columns=['x','fit'])
global NL, BL
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobal,T=PYTHON.
global Istatus, WarningText, ErrorText, Gdebug, Platform, System, Uname

Istatus = 0
WarningText = ''
ErrorText = ''

Platform = platform.platform().upper()
System = platform.system().upper()

if System != 'WINDOWS':
    os.system('uname > .uname')
    fun = open('.uname','r')
    Uname = fun.readline().strip()
    System = Uname[0:5].upper()
#if System != 'WINDOWS':

Gdebug = 0
# Histograms and Ntuples
global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman

Tdf = type(pd.DataFrame())

Tfig = None
Figman = None
Tax2d = None
Tax3d = None

H1h = None
H1 = []
H1head = []
Hdir = []
Ndir = 0
Kdir = 0
Cdir = ''

H1HLast = None

H2h = None
H2 = []
H2head = []

Ntup = []
Nctup = []
Ncind = []
H1ind = []
H2ind = []
Nind = []
Nhead = []

Nh1 = 0
Nh2 = 0

Nntup = 0
Nnctup = 0

H1Last = None
H2Last = None

H1h = None
H2h = None
N = None

Ind = -1
IndLast = -1

H1I = -1
H2I = -1
H1ILast = -1
H2ILast = -1

N = None
Nsel = None
NLast = None

Ni = -1
NiLast = -1

Nlines = 0
Ncolon = 0

#+seq,plotglobal.
#+PATCH,//WAVES/PYTHON
#+KEEP,nxyzglobal,T=PYTHON.
global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv,Nx,Nxy,Nxyz
N1 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N2 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N3 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N4 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N5 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N6 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N7 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N8 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
N9 = pd.DataFrame(columns=['x','y','z','s','t','bx','by','bz','b'])
Nv = pd.DataFrame(columns=['x'])
Nx = pd.DataFrame(columns=['x'])
Nxy = pd.DataFrame(columns=['x','y'])
Nxy = pd.DataFrame(columns=['x','y'])
Nxyz = pd.DataFrame(columns=['x','y','z'])

#CanButIds = []
#CanButId = 0
#Legend = []

Mode3ds = ['none','boxes','cont3d','hist','inter','surf','trisurf']
# +PATCH,//WAVES/PYTHON
# +KEEP,mhbglobal,T=PYTHON.
# }
#+seq,nxyzglobal.
#+seq,wglobal.
#+seq,statusglobal.


#+PATCH,//WAVES/PYTHON
#+KEEP,plotglobal,T=PYTHON.

import matplotlib as mpl

global MPLmain, MPLmaster, Nfigs,Figgeom, Figgeom2, FiggeomR, FiggeomL, \
XtermGeo, Figs,Fig,Ax, Figgeoms, \
Screewidth, Screenheight, ScaleSizeX, ScaleSizeY, \
Fig1,Ax1,Fig6,Ax6,Fig2,Ax2,Fig7,Ax7,Fig3,Ax3,Fig8,Ax8,\
Fig4,Ax4,Fig9,Ax9,Fig5,Ax5,Fig10,Ax10,\
FirstConsole, Console, Igetconsole,Klegend, Fwidth, Fheight, Fxoff, Fyoff, \
Kfig, Kax, Ihist,Iprof, Imarker, Ierr, Isurf, Iinter, Isame, Itight, IsameGlobal, Iline, CMap, Cmap, Tcmap, Surfcolor, \
Kplots, Nwins, Zones, Kzone, Nxzone, Nyzone, Zone, Axes, Legend, CanButId, CanButIds, \
Iplotopt, Ispline, Kecho, Kdump,Kpdf, Ndump,Npdf, Icmap, \
MarkerSize, MarkerType, MarkerColor, \
Markersize, Markertype, Markercolor, \
Fillstyle, FillStyle, WaveFilePrefix, \
Mode3d,Mode3D, Mode2d,Mode2D, \
Textcolor,Linestyle, Linewidth, Linecolor, Author, \
Histcolor, Histedgecolor, HistEdgeColor, Histbarwidth, Kdate, Kstat, Kfit, \
Icont3d, Iboxes, Inoempty, Iclosed, Itrisurf, Iscatter, Iscat3d, Ifill1d, TitPad, Xtitle, Ytitle, \
YTitle, Tightpad, Xtightpad,Ytightpad, ColorbarPad,\
LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace, \
Gtit,Xtit,Ytit,Ztit,Ttit,Ptit,Surfcolors, Cmaps, Colors, Linestyles, Markertypes, \
LexpX,LexpY,LexpRot,LexpPow, Xstat, XStat, XFit, YFit, Xfit, Yfit, Ystat, YStat, \
GtitFontSize,Titfontsize,Atitfontsize,Axislabelsize,Textfontsize,Datefontsize,\
Statfontsize, Axislabeldist, Axislabeldist3d, Axisdist, Axisdist3d,\
GtitFontSize,TitFontSize,AtitFontSize,AxisLabelSize,TextFontSize,DateFontSize,\
StatFontSize, AxisLabelDist, AxisLabelDist3d, AxisTitleDist, AxisTitleDist3d,\
AtitFontSize3d, Atitfontsize3d, NXtick,NXtick3d, Nxtick,Nxtick3d, Ktitles, Dummy,\
ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax, Tdate, TdateOv, Trun, TrunOv, \
LogX, LogY, LogZ, NxbBinMax, Khdeleted,WisLinux, Waveplot, \
Mrun, Mcomment, Mdate, ROFx, ROFy, Hull3D, Kgrid,KyAxis,KxAxis,KzAxis,Kbox, \
FillColor


FillColor = 'none'

NxBinMax = 0
Khdeleted =0

Tdate = None
TdateOv = None
Trun = None
TrunOv = None

WaveFilePrefix = 'wave_plot_'

ZoomXmin = -1.e30
ZoomXmax = 1.e30
ZoomYmin = -1.e30
ZoomYmax = 1.e30
LogX = 0
LogY = 0
LogZ = 0

Console = 'wavesPython'
FirstConsole = 1

Figgeom = ""
Figgeom2 = ""
FiggeomR = ""
FiggeomL = ""
XtermGeo = ""
Figgeoms = []

Fwidth= 500
Fheight = 300
Fxoff = 10
Fyoff = 10

Kdate = False
Kstat = False

Ktitles = 1
Kecho = 1
Kdump = True
Kpdf = True
Ndump = 0
Npdf = 0

Author = ''

Waveplot = None
Mrun = None
Mcomment = None
Mdate = None
ROFx = 0.03
ROFy = 0.95

fcfg = ''
if WavesMode == 'WAVES' or WavesMode == 'WPLOT': fcfg = 'waveplot.cfg'
elif WavesMode == 'UNDUMAG': fcfg = 'undugui.cfg'
else:
  if fexist('waveplot.cfg'):
    WavesMode = 'WAVES'
    fcfg = 'waveplot.cfg'
  else:
    fcfg = 'ntupplot.cfg'
  #endif
#endif

YFit = 0.8
Yfit = 0.8
Xfit = 0.75
XFit = 0.75

YStat = 0.8
Ystat = 0.8
Xstat = 0.05
XStat = 0.05

if fcfg != '' and os.path.exists(fcfg):

  wcon = open(fcfg,"r")

  while True:
    line = wcon.readline()
    if not len(line): break
    line = line.split('!')
    line = line[0].strip()
    words = line.split(':')
    if words[0] == 'Kstat':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kstat = True
      else: Kstat = False
    if words[0] == 'Kfit':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kstat = True
      else: KFit = False
    elif words[0] == 'Kbeam':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kbeam = True
      else: Kbeam = False
    elif words[0] == 'Kcurr':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kcurr = True
      else: Kcurr = False
    elif words[0] == 'Kcode':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kcode = True
      else: Kcode = False
    elif words[0] == 'Kecho':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kecho = True
      else: Kecho = False
    elif words[0] == 'Kpreload':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kpreload = True
      else: Kpreload = False
    elif words[0] == 'Koverview':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Ioverview = True
      else: Ioverview = False
    elif words[0] == 'Kgrid':
      if words[1].upper().strip() == 'FALSE' or words[1].upper().strip() == 'NO' or words[1].strip() == '0': Kgrid = False
      else: Kgrid = True
    elif words[0] == 'Kdate':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kdate = True
      else: Kdate = False
    elif words[0] == 'Kdump':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kdate = True
      else: Kdump = False
    elif words[0] == 'Kpdf':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Kdate = True
      else: Kpdf = False
    elif words[0] == 'xFit': XFit = float(words[1])
    elif words[0] == 'yFit': YFit = float(words[1])
    elif words[0] == 'xStat': XStat = float(words[1])
    elif words[0] == 'yStat': YStat = float(words[1])
    elif words[0] == 'ROFx': ROFx = float(words[1])
    elif words[0] == 'ROFy': ROFy = float(words[1])
    elif words[0] == 'Kconsole':
      if words[1].upper().strip() == 'TRUE' or words[1].upper().strip() == 'YES' or words[1].strip() == '1': Igetconsole = True
      else: Igetconsole = False
    elif words[0] == 'Author':
      if len(words) > 1: Author = words[1].strip()
  #endwhile True

  wcon.close()

else:
  Ioverview = 1
  Author = 'USER'
  Kfit = False
  Kstat = False
  Kecho = False
  Kbeam = True
  Kcurr = True
  Kdate = True
  Kpdf = True
  Kdump = True
  #Krun = True

#endif os.path.exists(".waveplot_config")

if Author.lower().strip() == 'user':

  try: Author = os.environ['USER']
  except:
    try: Author = os.environ['USERNAME']
    except:
      try: Author = os.getcwd()
      except: pass

  if Author != '':
    Al = list(Author)
    Al[0] = Al[0].upper()
    Author = ''
    for c in Al: Author += c
#endif Author.lower().strip() == 'user'

Kplots = list(np.linspace(1,1000,1000)*0)
Nwins = 0
Fig = None
Ax = None
Figs = []
Nfigs = 0
IsameGlobal = 0
Ifill1d = 0

TitPad = 12

GtitFontSize = 'large'
TitFontSize = 'large'
AtitFontSize = 11
AtitFontSize3d = AtitFontSize
AxisLabelSize = 10
TextFontSize = 11
DateFontSize = 9
StatFontSize = -9
AxisTitleDist = 6. # distance of title from axis
AxisTitleDist3d = 5. # distance of title from axis
Axistitledist3d = 5. # distance of title from axis
AxisLabelDist3d = 1. # distance of tick-labels from axis
AxisLabelDist = 5. # distance of tick-labels from axis

NXtick = -9
NXtick3d = 9

Xtitle = 0.5
Ytitle = 1.03
YTitle = 1.03
Itight = 0
Zones = []
Nxzone = 0
Nyzone = 0
Kzone = 0
Axes = []

Icmap = 1
CMap=None
CMap='winter'
CMap='hot'
CMap='autumn'
CMap='summer'
CMap='magma'
CMap='plasma'
CMap='RdPu'
CMap = 'copper'
CMap = 'rainbow'
CMap='jet'
Cmaps = ['none','autumn','copper','jet','hot','magma','plasma','rainbow','RdPu','summer','viridis','winter']
Cmap = CMap
Tcmap = mpl.colorbar.Colorbar

Mode3D = 'surf'
Mode2D = '!'

Surfcolor='salmon'
Surfcolor = 'red'
Surfcolor = 'navajowhite'
Surfcolor = 'mediumspringgreen'

Ispline = 0

Linecolor = 'red'
Linewidth = 1.
Markersize = 4.
Markercolor = 'red'
Linestyle = 'solid'
Markertype = 'o'
Fillstyle = 'full'
FillStyle = 'full'

LineColor = 'red'
LineWidth = 1.
MarkerSize = 6.
MarkerColor = 'red'
LineStyle = 'solid'
MarkerType = 'o'

Textcolor = 'black'

Histbarwidth = 1.
Histcolor = 'yellow'
HistEdgeColor = 'blue'
Histedgecolor = 'blue'

Klegend = 0
Igetconsole = True

Imarker = 0
Ihist = 0
Iprof = 0
Ispline = 0
Isame = 0
Iscatter = 0
Iscat3d = 0
Icont3d = 0
Itrisurf = 0
Iclosed = 0
Iboxes = 0
Inoempty = 0
Ierr = 0
Isurf = 0
Iline = 0
Iinter = 0
Ifill1d = 0
Gtit = ''
Colors = ['black','red','blue','green','cyan','magenta','gray','yellow', 'white']
Surfcolors = ['none','blue','cyan','gray','green','navajowhite','magenta','mediumspringgreen','red','salmon','yellow']
Linestyles = ['none','dashed','dotted','dashdot','solid']
Markertypes = ['none','.','o','v','s','*','^','<','>','+','P','x','X']
Tightpad = 2.5
Xtightpad = 0.5
Ytightpad = 1.0
ColorbarPad = '!'
LeftMargin = 0.12
RightMargin = 0.92
TopMargin = 0.86
BottomMargin = 0.13
#ScaleSizeX = 1.4143 # for savefig in printfig to adjust aspect ratio
ScaleSizeX = 1. # for savefig in printfig to adjust aspect ratio
ScaleSizeY = 1. # for savefig in printfig to adjust aspect ratio
Xspace = 0.4
Yspace = 0.5
LexpX = 0.
LexpY = 0.05
LexpRot = 'horizontal'
LexpPow = 2
Kgrid = True
Kbox = True
KxAxis = True
KyAxis = True
KzAxis = True

CanButIds = []
CanButId = 0


# +PATCH,//WAVES/PYTHON
# +KEEP,m_hbook,T=PYTHON.

from scipy import *
#+PATCH,//WAVES/PYTHON
#+KEEP,mshimports,T=PYTHON.

# import mshutil as u
# from mshutil import *

##############################################################

# begin of sequence mshimports

import tkinter as tk
from tkinter import *
from tkinter.font import Font
from tkinter import ttk
from tkinter import filedialog

import sys
import os
import time
import platform
import ctypes
import re
import fileinput
import shutil
from shutil import *
import subprocess
import glob

import numpy as np
import pandas as pd

from math import *

from scipy import *

import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm #color maps


global \
clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
ecdipev,ecdipkev,fwhmgauss1,fwhmsinxx21,rmssinxx21

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
hplanck1=6.626176e-34
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

fwhmgauss1=np.sqrt(2.0*np.log(2))*2.0
fwhmsinxx21=2.783115
rmssinxx21=1.05244

global Ftyp,Ftype
Ff = open("ftypedum","w")
Ftyp = type(Ff)
Ftype =type(Ff)
Ff.close()
os.remove("ftypedum")

global Legend
Legend = []
# end of sequence mshimports

#+seq,plotglobal.
#CMap = 'rainbow'
#HistEdgeColor = 'blue'
#AxisLabelSize = 10
#Nxzone = 1
#Nyzone = 1

def debug(kmenu=None,kitem=None):
  global \
  Wave, Root, WaveOut, Editor, WinGeo, \
  Debug,FWAVEIN,FWVS,Wavein,WaveinO,Wmenu,Mapping,MenuVeto,MenuAllVeto,NMenuAllVeto,\
  WAVECom, ROOTCom, EWOUTCom, EDICom, MenuActive, Nmenuactive, \
  MappingVeto,Veto,InvVeto,Variables,Arrays,Specvar,Trigger,Calc,Help,NHelp,Nmenu,\
  Nmenuveto,NmenuOld,Mother,MenuMother,MenuOld,Daughter,IMother,Nveto,\
  Ninvveto,Nmap,Nvar,Iarr,Narr,Ntrigger,Ncalc,NULL,ONE,MONE,SNULL,SONE,SMONE,Lastvar,\
  Nwavein,kWaveinRead, KWAVES,MMitem,Nmitem,Kmitem,Imenu,Ipmenu,\
  Kmitemold, Iback, Istak, Tcolor, Vetocolor, SFrame, SComment, SMitem, \
  VarToWaveIn, PosX, PosY, WinPos, Nsitem,Pmenu,PadX,PadY,SMexist, \
  I,ZONE,ZNULL,Ical,Lmitem,FIOitem,PMenuGeo, Nfocus, MyWavesFont,Ifocus, \
  ScreenW,ScreeH,WinX,WinY,CanW,CanH
  pass
#  print("\n\n debug::kmenu,kitem",kmenu,kitem,SMitem[kmenu])
#  print("\n\n debug::kmenu,kitem",kmenu,kitem)
#enddef debug(kmenu,kitem)

global Narg,Argv
Narg = len(sys.argv)
Argv = sys.argv

global Icallfromoverview, Krun
global Foverview
Foverview = None

Icallfromoverview = None
Krun = None
Waveplot = 0

global VarHead
VarHead = ''

global Vlocate, CanVlocate, \
VlocX,VlocY, VlocDistX,VlocDistY, VlocDist, \
VlocXO,VlocYO, VlocDistXO,VlocDistYO, VlocDistO

SplineMode = 'new'
ColorbarPad = 0.1
Vlocate = []
VlocXO = 0
VlocYO = 0
VlocX = 0
VlocY = 0
VlocDistX = 0
VlocDistY = 0
VlocDist = 0

global Aspect
Aspect = "auto"

def set_aspect(asp='!'):
  global Aspect

  if type(asp) == str:
    if asp == '!': return
    if asp == '0': asp = 'auto'
  else:
    if asp == 0: asp = 'auto'
  #endif

  Aspect = asp
#enddef

NL = "\n"

def sind(phi): return sin(phi/180.*np.pi)
def cosd(phi): return cos(phi/180.*np.pi)
def tand(phi): return tan(phi/180.*np.pi)
def asind(x): return asin(x)/np.pi*180.
def acosd(x): return acos(x)/np.pi*180.
def atand(x): return atan(x)/np.pi*180.
def cotd(phi): return cot(phi/180.*np.pi)

def set_y_title(y='!'):
  global YTitle, Ytitle, Nxzone, Nyzone
  if y == '!':
    y = YTitle
  if Nyzone > 1:
    y += (Nyzone-1)*0.05
  Ytitle = y
def get_y_title():
  return Ytitle

def set_x_fit(x='!'):
  global XFit, Xfit, Nxzone, Nyzone
  if x == '!': x = XFit
  Xfit = x
#enddef

def get_x_fit(): return Xfit

def set_y_fit(y='!'):
  global YFit, Yfit, Nxzone, Nyzone
  if y == '!': y = YFit
  if Nyzone > 1:y -= (Nyzone-1)*0.15
  Yfit = y
  YFit = y
#enddef

def set_x_stat(x='!'):
  global XStat, Xstat, Nxzone, Nyzone
  if x == '!':
    x = XStat
  elif x == '+':
    x = XStat + 0.2
  #endif
  if Nxzone > 1: x = (Nxzone-1)*0.15
  Xstat = x
  XStat = x
#enddef

def get_x_stat():
  global Xstat
  return Xstat
#endif

def set_y_stat(y='!'):
  global YStat, Ystat, Nxzone, Nyzone
  if y == '!': y = YStat
  if Nyzone > 1:y = 0.8 - (Nyzone-1)*0.15
  Ystat = y
  YStat = y
#enddef

def get_y_stat():
  global Ystat
  return Ystat
#endif

def ellipse(x0,y0,a,b,alpha=0.0,n=1000):

  cosa = cos(alpha)
  sina = sin(alpha)

  phi =vcre(n,0.,2.0*pi)

  x = x0 + cosa*a*cos(phi) - sina*b*sin(phi)
  y = y0 + sina*a*cos(phi) + cosa*b*sin(phi)

  istat,area,pathlen = polygon2d(x,y)

  return x,y,phi,area[:-1],pathlen[:-1]

#enddef ellipse

def nellipse(nt,x0=0.0,y0=0.0,a=1.0,b=1.0,alpha=0.0,n=1000):
  nt = ncre(nt,nt,"x:y:phi:area:path")
  nt.x,nt.y,nt.phi,nt.area,nt.path = ellipse(x0,y0,a,b,alpha,n)
  nupdate_header(nt)
  return nt
#enddef

def polygon2d(x="?",y="?"):

  if type(x) == str:
    print("polygon2d(x,y): Returns istat, the area[0:n+1] and pathlen[0:n+1];")
    print("pathlen[n+1] referes to the closed curve.")
    return
  #endif

  n = len(x)

  if n != len(y):
    istat = -1
    pathlen = -9999.
    area = -9999.
    print("*** Error in polygon2d: X and y must have the same length ***")
    return istat,area,pathlen
  #endif

  if n < 2:
    istat = 0
    pathlen = 0.0
    area = 0.0
  elif n == 2:
    istat = 0
    pathlen = 2.0*sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)
    area = 0.0
  #endif

  istat = 0
  area = vcre(n+1,0.0)
  pathlen = vcre(n+1,0.0)

  n1 = n -1
  for i in range(n1):
#    dx = x[i+1]-x[i]
#    dy = y[i+1]-y[i]
#    dl =  sqrt(dx**2+dy**2)
#    print(i,dx,dy,dl)
#    pathlen[i+1] = pathlen[i] + dl
    pathlen[i+1] = pathlen[i] + sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
    area[i+1] = area[i] - (x[i+1]-x[i]) * (y[i+1]+y[i])/2.0
  #endfor

  n2 = n -2

  pathlen[n] = pathlen[n1] + sqrt((x[0]-x[n1])**2 + (y[0]-y[n1])**2)
  area[n] = area[n1] - (x[0]-x[n1]) * (y[0]+y[n1])/2.0

  return istat,area,pathlen
#enddef

def get_master():
  wmain = plt.gcf()
  return wmain.canvas.toolbar.master
#enddef

def get_isame():
  global Isame
  return Isame
#enddef

def get_kpdf():
  global Kpdf
  return Kpdf
#enddef

def optdump(kdump='True'):
  global Kdump
  Kdump = kdump
#enddef

def optpdf(kpdf='True'):
  global Kpdf
  Kpdf = kpdf
#enddef

def store_kdump_kpdf():
  global KdumpOld,KpdfOld,Kdump,Kpdf
  KdumpOld = Kdump
  KpdfOld = Kpdf
#enddef

def restore_kdump_kpdf():
  global KdumpOld,KpdfOld,Kdump,Kpdf
  Kdump = KdumpOld
  Kpdf = KpdfOld
#enddef

def set_kpdf(k):
  global Kpdf
  Kpdf = k
#enddef

def get_kecho():
  global Kecho
  return Kecho
#enddef

def set_kecho(k):
  global Kecho
  Kecho = k
#enddef

def get_kdump():
  global Kdump
  return Kdump
#enddef

def set_kdump(k):
  global Kdump
  Kdump = k
#enddef

def latex(xNDC,yNDC,text,color='black',fontsize=12,halign='left',valign='top', \
bbox='none',bbstyle='round',fc='white',ec='black',pad=0.2):

  mpl.rcParams['text.usetex'] = True
  ax = plt.gca()

  x,y = NDC_to_WC(xNDC,yNDC)

  if bbox != 'none':
    ax.text(x,y,text,color=color, fontsize=fontsize, \
    horizontalalignment=halign, verticalalignment=valign, \
    bbox=dict(boxstyle=bbstyle, fc=fc, ec=ec, pad=pad))
  else:
    ax.text(x,y,text,color=color, fontsize=fontsize, \
    horizontalalignment=halign, verticalalignment=valign)
  #endif

#enddef latex()

global MShWelcome
MShWelcome = False

def get_mshwelcome():
  global MShWelcome
  return MShWelcome
#enddef

def mshwelcome(program='Wave-Plot',year='2021', getconsole=False):

  global Kdate, Kbox, KxAxis, KyAxis, MShWelcome

  KdateO, KboxO, KxAxisO, KyAxisO = Kdate, Kbox, KxAxis, KyAxis

  MShWelcome = True

  optdate(False)

  window(program,getconsole=getconsole)

  null(0,10,0,10)
  optxaxis(False); optyaxis(False); optbox(False)

  text = ""
  text += "Welcome to " + program + "\n\n"
  text += "by Michael Scheer \n Helmholtz-Zentrum Berlin\n\n"
  textndc(0.5,0.8,text,fontsize=15,color='magenta')

  text = ""
  text += "Copyright " + str(year) + " Helmholtz-Zentrum Berlin (HZB)"
  text += "\n      Hahn-Meitner-Platz 1"
  text += "\n      D-14109 Berlin"
  text += "\n      Germany"
  text += "\n"
  text += "\n      Michael.Scheer@Helmholtz-Berlin.de"
  text += "\n"
  text += "\n -----------------------------------------------------------------------"
  text += "\n"
  text += "\n    This program is free software: you can redistribute it and/or modify"
  text += "\n    it under the terms of the GNU General Public License as published by"
  text += "\n    the Free Software Foundation, either version 3 of the License, or"
  text += "\n    (at your option) any later version."
  text += "\n"
  text += "\n    This program is distributed in the hope that it will be useful,"
  text += "\n    but WITHOUT ANY WARRANTY; without even the implied warranty of"
  text += "\n    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
  text += "\n    GNU General Public License for more details."

  textndc(0.5,0.4,text,fontsize=10,color='magenta')

  Kdate, Kbox, KxAxis, KyAxis = KdateO, KboxO, KxAxisO, KyAxisO
#enddef mshwelcome()

def macro(fmacro=''):

  if not len(fmacro): return

  if os.path.exists(fmacro):
    exec(compile(open(fmacro).read(), fmacro, 'exec'))
  else:
    if os.path.exists(fmacro + ".py"):
      fmacro += '.py'
      exec(compile(open(fmacro).read(), fmacro, 'exec'))
    elif os.path.exists("../py/" + fmacro):
      fmacro = "../py/" + fmacro
      exec(compile(open(fmacro).read(), fmacro, 'exec'))
    elif os.path.exists("../py/" + fmacro + ".py"):
      fmacro = "../py/" + fmacro + ".py"
      exec(compile(open(fmacro).read(), fmacro, 'exec'))
    else:
      print("*** " + fmacro + " not found")
  #endif
#enddef macro(fmacro='')

def rotate2d(x,y,phi=0.0, x0=0.0, y0=0.0, mode='grad'):
  x = np.array(x)
  y = np.array(y)

  if mode.lower() == 'grad': phi = phi / 180. * pi

  cosphi = cos(phi)
  sinphi = sin(phi)

  xx0 = x - x0
  yy0 = y - y0

  xr = x0 + xx0 * cosphi - yy0 * sinphi
  yr = y0 + xx0 * sinphi + yy0 * cosphi

  return xr, yr
#enddef rotate2d()

def draw_circle(x=0.0,y=0.0,r=1.,fill=False,lwid=1):
    ax = plt.gca()
    circ = plt.Circle((x,y),radius=r,fill=fill,linewidth=lwid)
    ax.add_patch(circ)
    showplot()
#enddef

def setmode2d(m='!'):
  global Mode2d, Mode2D
  if m == '!': Mode2d = Mode2D
  else: Mode2d = m
#enddef setmode2d(m='!')

def getmode2d():
  global Mode2D
  return Mode2D
#enddef getmode2d()

def setmode3d(m='surf'):
  global Mode3D; Mode3D = m
  _setmode3d()
#enddef setmode3d(m='surf')

def getmode3d():
  global Mode3D
  return Mode3D
#enddef getmode3d()

def set_frame_geom(wf=0,hf=0,window='!'):

    global Wmaster

    if window == '!': window = Wmaster

    swid,shei,w,ho,canw,canho,x,y,bbox,axw,axho = get_geo_all()

    if wf <= 0: wf = w
    if hf <= 0: hf = ho

    r = wf / axw
    wf = int(wf * r)
    hf = int(hf * r)

    sgeo =  str(wf) + 'x' + str(hf) + '+' + str(x) + '+' + str(y)
    Wmaster.geometry(sgeo)
    Wmaster.update()

    swid,shei,w,ho,canw,canho,x,y,bbox,axw,axho = get_geo_all()

    dh = axw-axho
    hnew = int(ho+dh)
    sgeo =  str(int(w)) + 'x' + str(hnew) + '+' + str(x) + '+' + str(y)
    Wmaster.geometry(sgeo)
    Wmaster.update()

    swid,shei,w,h,canw,canh,x,y,bbox,axw,axh = get_geo_all()

    # axh = axh0 + da/dh|0 * (hw-h0) =! axw*hf/wf
    hnew = int((axw*hf/wf - axho) * (h-ho)/(axh-axho) + ho)

    sgeo =  str(int(w)) + 'x' + str(hnew) + '+' + str(x) + '+' + str(y)
    Wmaster.geometry(sgeo)
    Wmaster.update()

    #print(get_geo_all())

#enddef square_geom()

def set_frame_square(wf=0,window='!'):
  set_frame_geom(wf,wf,window='!')
#enddef

def get_geo_all():
    global Wmaster,WinX,WinY,ScreenW,ScreenH,CanW,CanH

    wid,h,WinX,WinY = getgeo()
    ScreenW = Wmaster.winfo_screenwidth()
    ScreenH = Wmaster.winfo_screenheight()

    fig = plt.gcf()
    CanW,CanH = fig.canvas.get_width_height()

    ax = plt.gca()

    bbox = ax.get_window_extent()
    axw = bbox.width
    axh = bbox.height

    return ScreenW,ScreenH,wid,h,CanW,CanH,WinX,WinY,bbox,axw,axh
#enddef get_geo_all()

#Elast

#Internal functions {
def _clearCanvas():
  window_clear()
  showplot(False)
#enddef _clearCanvas()

def _delPlot():
    fig = plt.gcf()
    plot = plt.gca()
    fig.delaxes(plot)
    shpl()
#enddef _delPlot()

def _zones():
    global Wmain, Wmaster, Winz, Erows, Ecols, Ekzon, Myfont,Nyzone,Nxzone

    Winz = Toplevel()
    Winz.attributes('-topmost', 1)

    fr = Frame(Winz)

    lrow = Label(fr,text='number of rows',font=Myfont)
    lrow.pack(side=LEFT)

    Erows = Entry(fr,width=3,font=Myfont)
    Erows.insert(0,Nyzone)
    Erows.pack(fill=X,side=RIGHT)

    fc = Frame(Winz)
    lcol = Label(fc,text='number of columns',font=Myfont)
    lcol.pack(side=LEFT)
    Ecols = Entry(fc,width=3,font=Myfont)
    Ecols.insert(1,Nxzone)
    Ecols.pack(fill=X,side=RIGHT)

    fk = Frame(Winz)
    lk = Label(fk,text='selection',font=Myfont)
    lk.pack(side=LEFT)

    Ekzon = Entry(fk,width=3,font=Myfont)
    Ekzon.insert(2,Kzone)
    Ekzon.pack(fill=X,side=RIGHT)

    fr.pack(fill=BOTH)
    fc.pack(fill=BOTH)
    fk.pack(fill=BOTH)

    cbClear= Checkbutton(Winz,text="Same canvas",  onvalue=1,
                         offvalue=0, variable=IsameCanvas).pack()

    bClose = Button(Winz,text='Ok',command=_setzones).pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-170)
    Winz.geometry(sgeo)

    Wmaster.wait_window(Winz)
#enddef _zones():

def _lines():
    global Winl, Myfont, Ewid, combocol, combosty, Wmaster

    Winl = Toplevel()
    Winl.attributes('-topmost', 1)

    flc = Frame(Winl)
    llc = Label(flc,text='color',width=5,font=Myfont)
    llc.pack(side=LEFT)
    combocol = ttk.Combobox(flc,values=Colors)
    idx = getcolorindex(Linecolor)
    if idx >= 0 and idx < len(Colors): combocol.current(idx)
    combocol.bind("<<ComboboxSelected>>",_combocol)
    combocol.pack(side=RIGHT)
    flc.pack()

    fls = Frame(Winl)
    lls = Label(fls,text='style',width=5,font=Myfont)
    lls.pack(side=LEFT)
    combosty = ttk.Combobox(fls,values=Linestyles)
    idx = getlinestyleindex(Linestyle)
    if idx >= 0 and idx < len(Linestyles): combosty.current(idx)
    combosty.bind("<<ComboboxSelected>>",_combosty)
    combosty.pack(side=RIGHT)
    fls.pack()

    fw = Frame(Winl)
    lwid = Label(fw,text='width',font=Myfont)
    lwid.pack(side=LEFT)
    Ewid = Entry(fw,width=3,font=Myfont)
    Ewid.insert(1,Linewidth)
    Ewid.pack(fill=X,side=RIGHT)
    fw.pack(fill=X,padx=5)

    bClose = Button(Winl,text='Ok',command=_closewinl).pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-140)
    Winl.geometry(sgeo)

    Wmaster.wait_window(Winl)
#enddef _lines()

def _marker():
    global Winm, Myfont, Esiz, combocol, combosty, Markersize

    Winm = Toplevel()
    Winm.attributes('-topmost', 1)

    fmc = Frame(Winm)
    lmc = Label(fmc,text='color',width=5,font=Myfont)
    lmc.pack(side=LEFT)
    combomcol = ttk.Combobox(fmc,values=Colors)
    idx = getcolorindex(Linecolor)
    if idx >= 0 and idx < len(Colors): combomcol.current(idx)
    combomcol.bind("<<ComboboxSelected>>",_combomcol)
    combomcol.pack(side=RIGHT)
    fmc.pack()

    fmt = Frame(Winm)
    lmt = Label(fmt,text='type',width=5,font=Myfont)
    lmt.pack(side=LEFT)
    combosty = ttk.Combobox(fmt,values=Markertypes)
    idx = getmarkertypeindex(Markertype)
    if idx >= 0 and idx < len(Markertypes): combosty.current(idx)
    combosty.bind("<<ComboboxSelected>>",_combosty)
    combosty.pack(side=RIGHT)
    fmt.pack()

    fsiz = Frame(Winm)
    lwid = Label(fsiz,text='size',font=Myfont)
    lwid.pack(side=LEFT)
    Esiz = Entry(fsiz,width=3,font=Myfont)
    Esiz.insert(1,Markersize)
    Esiz.pack(fill=X,side=RIGHT)
    fsiz.pack(fill=X,padx=5)

    bClose = Button(Winm,text='Ok',command=_closewinm).pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-140)
    Winm.geometry(sgeo)

    global Wmaster
    Wmaster.wait_window(Winm)
#enddef _marker()

def _options3d():
    global Wmain, Wmaster, Win3d, Myfont, comboscol, combocmap, combomod3d, Cmaps

    Win3d = Toplevel()
    Win3d.attributes('-topmost', 1)

    fmo = Frame(Win3d)
    lcm = Label(fmo,text='Mode',width=10,font=Myfont)
    lcm.pack(side=LEFT)
    combomod3d = ttk.Combobox(fmo,values=Mode3ds)
    idx = getmode3dindex(Mode3d)
    if idx >= 0 and idx < len(Mode3ds): combomod3d.current(idx)
    combomod3d.bind("<<ComboboxSelected>>",_combomod3d)
    combomod3d.pack(side=RIGHT)
    fmo.pack()

    fcm = Frame(Win3d)
    lcm = Label(fcm,text='color map',width=10,font=Myfont)
    lcm.pack(side=LEFT)
    combocmap = ttk.Combobox(fcm,values=Cmaps)
    idx = getcmapindex(Cmap)
    if idx >= 0 and idx < len(Cmaps): combocmap.current(idx)
    combocmap.bind("<<ComboboxSelected>>",_combocmap)
    combocmap.pack(side=RIGHT)
    fcm.pack()

    flc = Frame(Win3d)
    llc = Label(flc,text='surface color',width=10,font=Myfont)
    llc.pack(side=LEFT)
    comboscol = ttk.Combobox(flc,values=Surfcolors)
    idx = getsurfcolorindex(Surfcolor)
    if idx >= 0 and idx < len(Surfcolors): comboscol.current(idx)
    comboscol.bind("<<ComboboxSelected>>",_comboscol)
    comboscol.pack(side=RIGHT)
    flc.pack()

    cbCmap= Checkbutton(Win3d,text="Use color map",  onvalue=1, offvalue=0, variable=Icmap).pack()

    bClose = Button(Win3d,text='Ok',command=_setopt3d).pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-170)
    Win3d.geometry(sgeo)

    Wmaster.wait_window(Win3d)
#enddef _options3d()

def _optdate():
  global Kdate
  if Kdate == False: Kdate = True
  else: Kdate = False
  date_on_figure()
#enddef _optdate()

def _optgrid():
  global Kgrid
  if Kgrid == False: Kgrid = True
  else: Kgrid = False
  plt.grid()
  showplot()
#enddef _optgrid()

def _optxaxis():
  global KxAxis
  if KxAxis == False: KxAxis = True
  else: KxAxis = False
  ax = plt.gca()
  ax.xaxis.set_visible(KxAxis)
  showplot()
#enddef _optxaxis()

def _optbox():
  global Kbox
  if Kbox == False: Kbox = True
  else: Kbox = False
  plt.box(Kbox)
  showplot()
#enddef _optbox()

def _optyaxis():
  global KyAxis
  if KyAxis == False: KyAxis = True
  else: KyAxis = False
  ax = plt.gca()
  ax.yaxis.set_visible(KyAxis)
  showplot()
#enddef _optyaxis()

def _togglelinlogx():
    global LogX
    if LogX: LogX = 0
    else: LogX = 1
#enddef

def _togglelinlogy():
    global LogY
    if LogY: LogY = 0
    else: LogY = 1
#enddef

def _togglelinlogz():
    global LogZ
    if LogZ: LogZ = 0
    else: LogZ = 1
#enddef

def _setlinlog():

    global WinLinLog

    try: Mmenu.unpost()
    except: pass

    WinLinLog = Toplevel()
    WinLinLog.attributes('-topmost', 1)

    bllx = Button(WinLinLog,text='Toggle Lin/Log x',command=_togglelinlogx)
    bllx.pack()

    blly = Button(WinLinLog,text='Toggle Lin/Log y',command=_togglelinlogy)
    blly.pack()

    bllz = Button(WinLinLog,text='Toggle Lin/Log z',command=_togglelinlogz)
    bllz.pack()

    bClose = Button(WinLinLog,text='Ok',command=_closewinlinlog)
    bClose.pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-180)
    WinLinLog.geometry(sgeo)

    Wmaster.wait_window(WinLinLog)

#enddef _setlinlog()

def _setzoom():

    global WinZoom, Wmaster, \
    EZoomXmin, EZoomXmax, EZoomYmin, EZoomYmax,\
    ZoomXmin, ZoomXmax, ZoomYmin, ZoomYmax

    try: Mmenu.unpost()
    except: pass

    WinZoom = Toplevel()
    WinZoom.attributes('-topmost', 1)

    fxmin = Frame(WinZoom)
    lxmin = Label(fxmin,text='ZoomXmin',font=Myfont)
    lxmin.pack(side=LEFT)
    EZoomXmin = Entry(fxmin,width=10,font=Myfont)
    EZoomXmin.insert(1," {:.4g}".format(ZoomXmin))
    EZoomXmin.pack(fill=X,side=RIGHT)
    fxmin.pack(fill=X,padx=5)

    fxmax = Frame(WinZoom)
    lxmax = Label(fxmax,text='ZoomXmax',font=Myfont)
    lxmax.pack(side=LEFT)
    EZoomXmax = Entry(fxmax,width=10,font=Myfont)
    EZoomXmax.insert(1," {:.4g}".format(ZoomXmax))
    EZoomXmax.pack(fill=X,side=RIGHT)
    fxmax.pack(fill=X,padx=5)

    fymin = Frame(WinZoom)
    lymin = Label(fymin,text='ZoomYmin',font=Myfont)
    lymin.pack(side=LEFT)
    EZoomYmin = Entry(fymin,width=10,font=Myfont)
    EZoomYmin.insert(1," {:.4g}".format(ZoomYmin))
    EZoomYmin.pack(fill=X,side=RIGHT)
    fymin.pack(fill=X,padx=5)

    fymax = Frame(WinZoom)
    lymax = Label(fymax,text='ZoomYmax',font=Myfont)
    lymax.pack(side=LEFT)
    EZoomYmax = Entry(fymax,width=10,font=Myfont)
    EZoomYmax.insert(1," {:.4g}".format(ZoomYmax))
    EZoomYmax.pack(fill=X,side=RIGHT)
    fymax.pack(fill=X,padx=5)

    bAbort = Button(WinZoom,text='Cancel',command=_abortwinzoom)
    bAbort.pack()

    bClose = Button(WinZoom,text='Ok',command=_closewinzoom)
    bClose.pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-180)
    WinZoom.geometry(sgeo)

    Wmaster.wait_window(WinZoom)

#enddef _setzoom()

def _setuser():

    global WinUser, Author, EUs, Wmaster

    try: Mmenu.unpost()
    except: pass

    WinUser = Toplevel()
    WinUser.attributes('-topmost', 1)

    fus = Frame(WinUser)
    lwid = Label(fus,text='User name',font=Myfont)
    lwid.pack(side=LEFT)
    EUs = Entry(fus,width=8,font=Myfont)
    EUs.insert(1,Author)
    EUs.pack(fill=X,side=RIGHT)
    fus.pack(fill=X,padx=5)

    bClose = Button(WinUser,text='Ok',command=_closewinuser)
    bClose.pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-160)
    WinUser.geometry(sgeo)

    Wmaster.wait_window(WinUser)

#enddef _setuser()

def _setmode2d(m='!'):
  global Debug
  #if Debug: Quit("Ende in set_plot_params")
  setmode2d(m)
#def _setmode2d(m='!')

def _setmode3d(m='!'):
  global Mode3d, Mode3D
  if m == '!': Mode3d = Mode3D
  else: Mode3d = m
#enddef _setmode3d(m='!')

def _optstat():
  global Kstat
  if Kstat == False: Kstat = True
  else: Kstat = False

def _optzaxis():
  global KzAxis
  if KzAxis == False: KzAxis = True
  else: KzAxis = False
  ax = plt.gca()
  ax.zaxis.set_visible(KzAxis)
  showplot()
#enddef _optzaxis()

def _optrun():
  global Krun,Kruns
  if Krun == False: Krun = True
  else: Krun = False
  run_on_figure()
#enddef

def _settitles():
  global Wint, Eglob, Eplot, Extit, Eytit, Eztit

  gtit = Eglob.get()
  pltit = Eplot.get()
  xtit = Extit.get()
  ytit = Eytit.get()

  set_global_title(gtit)
  set_title(pltit)
  set_x_title(xtit)
  set_y_title(ytit)

  if hasattr(Ax,'zaxis'):
    ztit = Eztit.get()
    set_z_title(ztit)
    # not yet: set_t_title(Ttit)

  Wint.destroy()
#enddef _settitles()

def _titles():
    global Wint, Eglob, Eplot, Extit, Eytit, Eztit, Myfont, Ptit

    Wint = Toplevel()
    Wint.attributes('-topmost', 1)

    fg = Frame(Wint)
    lglob = Label(fg,text='Global Title',font=Myfont)
    lglob.pack(side=LEFT)
    Eglob = Entry(fg,width=45,font=Myfont)
    Eglob.insert(0,Gtit)
    Eglob.pack(fill=X,side=RIGHT)

    fpl = Frame(Wint)
    lplot = Label(fpl,text='Plot Title',font=Myfont)
    lplot.pack(side=LEFT)
    Eplot = Entry(fpl,width=45,font=Myfont)
    Eplot.insert(0,Ptit)
    Eplot.pack(fill=X,side=RIGHT)

    fx = Frame(Wint)
    lxtit = Label(fx,text='x Title',font=Myfont)
    lxtit.pack(side=LEFT)
    Extit = Entry(fx,width=45,font=Myfont)
    Extit.insert(0,Xtit)
    Extit.pack(fill=X,side=RIGHT)

    fy = Frame(Wint)
    lytit = Label(fy,text='y Title',font=Myfont)
    lytit.pack(side=LEFT)
    Eytit = Entry(fy,width=45,font=Myfont)
    Eytit.insert(0,Ytit)
    Eytit.pack(fill=X,side=RIGHT)

    fg.pack(fill=BOTH)
    fpl.pack(fill=BOTH)
    fx.pack(fill=BOTH)
    fy.pack(fill=BOTH)

    if hasattr(Ax,'zaxis'):
      fz = Frame(Wint)
      lztit = Label(fz,text='z Title',font=Myfont)
      lztit.pack(side=LEFT)
      Eztit = Entry(fz,width=45,font=Myfont)
      Eztit.insert(0,Ztit)
      Eztit.pack(fill=X,side=RIGHT)
      fz.pack(fill=BOTH)

    bClose = Button(Wint,text='Ok',command=_settitles).pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-170)
    Wint.geometry(sgeo)

    global Wmaster
    Wmaster.wait_window(Wint)
#enddef _titles()

def _options3d():
    global Wmain, Wmaster, Win3d, Myfont, comboscol, combocmap, combomod3d, Cmaps

    Win3d = Toplevel()
    Win3d.attributes('-topmost', 1)

    fmo = Frame(Win3d)
    lcm = Label(fmo,text='Mode',width=10,font=Myfont)
    lcm.pack(side=LEFT)
    combomod3d = ttk.Combobox(fmo,values=Mode3ds)
    idx = getmode3dindex(Mode3d)
    if idx >= 0 and idx < len(Mode3ds): combomod3d.current(idx)
    combomod3d.bind("<<ComboboxSelected>>",_combomod3d)
    combomod3d.pack(side=RIGHT)
    fmo.pack()

    fcm = Frame(Win3d)
    lcm = Label(fcm,text='color map',width=10,font=Myfont)
    lcm.pack(side=LEFT)
    combocmap = ttk.Combobox(fcm,values=Cmaps)
    idx = getcmapindex(Cmap)
    if idx >= 0 and idx < len(Cmaps): combocmap.current(idx)
    combocmap.bind("<<ComboboxSelected>>",_combocmap)
    combocmap.pack(side=RIGHT)
    fcm.pack()

    flc = Frame(Win3d)
    llc = Label(flc,text='surface color',width=10,font=Myfont)
    llc.pack(side=LEFT)
    comboscol = ttk.Combobox(flc,values=Surfcolors)
    idx = getsurfcolorindex(Surfcolor)
    if idx >= 0 and idx < len(Surfcolors): comboscol.current(idx)
    comboscol.bind("<<ComboboxSelected>>",_comboscol)
    comboscol.pack(side=RIGHT)
    flc.pack()

    cbCmap= Checkbutton(Win3d,text="Use color map",  onvalue=1, offvalue=0, variable=Icmap).pack()

    bClose = Button(Win3d,text='Ok',command=_setopt3d).pack()

    wid,h,x,y = getgeo()
    sgeo = '+' + str(x+wid//3) + '+' + str(y+h-170)
    Win3d.geometry(sgeo)

    Wmaster.wait_window(Win3d)
#enddef _options3d()

def _nextzone():
    global Wins, Erows, Ecols, Ekzon, Winz, Nyzone, Nxzone

    ny = Nyzone
    nx = Nxzone
    k = Kzone + 1

    if k > ny * nx:
      k = 1
      zone(nx,ny,k)
    else:
      zone(nx,ny,k,'s')
#enddef _nextzone()

def _setzones():
    global Wins, Erows, Ecols, Ekzon, Winz, Nyzone, Nxzone

    ny=int(Erows.get())
    nx=int(Ecols.get())
    k=int(Ekzon.get())

    if not ny: ny = Nyzone
    if not nx: nx = Nxzone
    if not k: k = Kzone

    #print("\n\n_setzones:",IsameCanvas.get(),"\n\n")

    if IsameCanvas.get(): zone(nx,ny,k,'s')
    else: zone(nx,ny,k)

    Winz.destroy()
#    Ax = plt.gca()
#enddef _setzones()

def _setopt3d():
    global Win3d

    Win3d.destroy()
#def _setopt3d()
def _vlocate(ev):

  global Vlocate, CanVlocate, Fig, \
  VlocX,VlocY, VlocDistX,VlocDistY, VlocDist, \
  VlocXO,VlocYO, VlocDistXO,VlocDistYO, VlocDistO

  but = ev.button
  dclick = ev.dblclick

  if dclick:
    Fig.canvas.mpl_disconnect(CanVlocate)
    print("\n Vlocate:\n")
    for elem in Vlocate:
      print(elem)
    #endfor elem in Vlocate
    return
  #endif dclick

  if but == 1:
    VlocXO = VlocX
    VlocYO = VlocY
    VlocX = ev.xdata
    VlocY = ev.ydata
    VlocDistX = VlocX - VlocXO
    VlocDistY = VlocY - VlocYO
    VlocDist = sqrt(VlocDistX*VlocDistX + VlocDistY*VlocDistY)
    Vlocate.append([VlocX,VlocY])
    print("Point:",VlocX,VlocY)
    if len(Vlocate) > 1: print("Distances:",VlocDistX,VlocDistY,VlocDist)
  elif but == 2:
    print(Vlocate.pop(-1)," deleted")
  #endif but == 1:

#def _vlocate(ev)

def _mhclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
#cid = fig.canvas.mpl_connect('button_press_event', onclick)
#def _mhclick(event)

def _sethistcolor(hc='!'):
  global Histedgecolor, HistEdgeColor
  if hc == '!': hc = HistEdgeColor
  Histedgecolor = hc
#enddef

def _setaxislabelsize(size='!'):
  global AxisLabelSize,Axislabelsize,Nxzone,Nyzone
  if size == '!': size = AxisLabelSize
  if Nxzone*Nyzone > 1:
    size *= 0.8
  Axislabelsize = size
#enddef

def _setaxistitlesize3d(size='!'):
  global AtitFontSize3d,Nxzone,Nyzone,Atitfontsize3d
  if size == '!': size = AtitFontSize3d
  if Nxzone*Nyzone > 1:
    size *= 0.9
  Atitfontsize3d = size
  plt.rcParams['axes.labelsize'] = size
#enddef

def _setaxistitlesize(size='!'):
  global AtitFontSize,Nxzone,Nyzone,Atitfontsize
  if size == '!': size = AtitFontSize
  if Nxzone*Nyzone > 1:
    size *= 0.9
  Atitfontsize = size
  plt.rcParams['axes.labelsize'] = size
def getaxistitlesize(): return Atitfontsize

def _setaxistitledist(dist='!'):
  global AxisTitleDist ,Nxzone,Nyzone, Axistitledist
  if dist == '!': dist = AxisTitleDist
  Axistitledist = dist
  if Nxzone*Nyzone > 1:
    dist *= 0.5
  mpl.rcParams['axes.labelpad'] = Axistitledist
#enddef

def _setaxistitledist3d(dist='!'):
  global AxisTitleDist3d ,Nxzone,Nyzone, Axistitledist3d
  if dist == '!': dist = AxisTitleDist3d
  if Nxzone*Nyzone > 1:
    dist -= 4.
  Axistitledist3d = dist
  mpl.rcParams['axes.labelpad'] = Axistitledist3d
#enddef

def _setaxislabeldist(dist='!'):
  global AxisLabelDist, Nxzone, Nyzone, Axislabeldist
  if dist == '!': dist = AxisLabelDist
  if Nxzone*Nyzone > 1:
    dist *= 0.5
  Axislabeldist = dist
#enddef

def _setaxislabeldist3d(dist='!'):
  global AxisLabelDist3d, Nxzone, Nyzone, Axislabeldist3d
  if dist == '!': dist = AxisLabelDist3d
  if Nxzone*Nyzone > 1:
    dist *= 0.5
  Axislabeldist3d = dist
#enddef

def _set_number_of_ticks_3d(n=-9): #, axis=''):
  global Nxtick3d, NXtick3d, Nxzone, Nyzone
  if NXtick3d <= 0: return # automatically
  if n <= 0: Nxtick3d = NXtick3d
  if Nxzone*Nyzone > 1: Nxtick3d -= 3
  #plt.locator_params(axis='x',nbins=Nxtick3d)
  plt.locator_params(tight=True,nbins=Nxtick3d)
#enddef

def _set_number_of_ticks(n=-9): #, axis=''):
  global Nxtick, NXtick, Nxzone, Nyzone
  if NXtick <= 0: return # automatically
  if n <= 0: Nxtick = NXtick
  if Nxzone*Nyzone > 1: Nxtick -= 3
  #plt.locator_params(axis=axis,nbins=Nxtick)
  plt.locator_params(tight=True,nbins=Nxtick)
#enddef

def _setcolormap(cmap='!'):
  global Cmap, CMap
  if cmap == '!': cmap = CMap
  Cmap = cmap
#enddef

def _setcolorbarpad(pad='!'):
  global ColorbarPad, Colorbarpad
  if pad == '!': pad = ColorbarPad
  Colorbarpad = pad
#enddef

#}Internal functions

def optbox(k=True):
  global Kbox
  Kbox = k
  plt.box(Kbox)
  showplot()
#enddef optbox()

def optxaxis(k=True):
  global KxAxis
  KxAxis = k
  ax = plt.gca()
  ax.xaxis.set_visible(KxAxis)
  showplot()
#enddef optxaxis()

def optyaxis(k=True):
  global KyAxis
  KyAxis = k
  ax = plt.gca()
  ax.yaxis.set_visible(KyAxis)
  showplot()
#enddef optyaxis()

def optzaxis(k=True):
  global KzAxis
  KzAxis = k
  ax = plt.gca()
  ax.zaxis.set_visible(KzAxis)
  showplot()
#enddef optzaxis()

def yesno(sin):
  sin = str(sin).upper()
  if sin == '0' or sin[0] == 'N' or sin == 'F': sout = 'no'
  elif sin == '1' or sin[0] == 'Y' or sin == 'T' or sin[0] == 'J': sout = 'yes'
  else:
    print("Bad input for " + sin)
    return 'BAD'
  #endif
  return sout
#enddef yesno(sin)

def yesno_01(sin):
  sin = str(sin).upper()
  if sin == '0' or sin[0] == 'N' or sin == 'F': iout = 0
  elif sin == '1' or sin[0] == 'Y' or sin == 'T' or sin[0] == 'J': iout = 1
  else:
    print("Bad input for " + sin)
    return -9999
  #endif
  return iout
#enddef yesno_01(sin)

####################################################################

def vcre(nx=10,xmin=1,xmax=10):

  if type(nx) == int:
    if nx <= 0:
      print("*** Error in vre: Negative or zero dimension ***")
      return -1
    #endif nx <= 0
    return vec(xmin,xmax,nx)
  else:
    if xmax <= 0:
      print("*** Error in vre: Negative or zero dimension ***")
      return -1
    #endif nx <= 0
    return vec(nx,xmin,xmax)
  #endif

#def vcre(n=10,xmin=1,xmax=10)

def getplotsize():
  w,h = plt.gcf().get_size_inches()
  return w*2.54,h*2.54
#enddef getplotsize()

def plothull3d(isame=0,facecolor='blue',alpha=0.5,edgecolor='black', ishow=1, mode='face'):
  global Hull3D, Ax, Isame

  Isame = isame
  getzone('3d')

  if not isame:
    xmina = []
    xmaxa = []
    ymina = []
    ymaxa = []
    zmina = []
    zmaxa = []
    for p in Hull3D:
      xmin = amin(tuple(zip(p[0],p[1],p[2]))[0])
      xmax = amax(tuple(zip(p[0],p[1],p[2]))[0])
      ymin = amin(tuple(zip(p[0],p[1],p[2]))[1])
      ymax = amax(tuple(zip(p[0],p[1],p[2]))[1])
      zmin = amin(tuple(zip(p[0],p[1],p[2]))[2])
      zmax = amax(tuple(zip(p[0],p[1],p[2]))[2])
      xmina.append(xmin)
      xmaxa.append(xmax)
      zmina.append(zmin)
      zmaxa.append(zmax)
      ymina.append(ymin)
      ymaxa.append(ymax)
    #endfor p  in Hull3D
    xmin = amin(xmina)
    xmax = amax(xmaxa)
    ymin = amin(ymina)
    ymax = amax(ymaxa)
    zmin = amin(zmina)
    zmax = amax(zmaxa)
  #endif not isame:

  if mode == 'face':

    face = mplot3d.art3d.Poly3DCollection(Hull3D)
    face.set_color(facecolor)
    face.set_edgecolor(edgecolor)
    face.set_alpha(alpha)

    Ax.add_collection3d(face)

  else:

    for closed in Hull3D:
      Ax.plot(*zip(*closed),edgecolor)
    #endfor closed in Hull3D

  #endif mode == 'face'

  if not isame:
    dx = xmax - xmin
    Ax.set_xlim3d(xmin-dx/10.,xmax+dx/10.)
    dy = ymax - ymin
    Ax.set_ylim3d(ymin-dy/10.,ymax+dy/10.)
    dz = zmax - zmin
    Ax.set_zlim3d(zmin-dz/10.,zmax+dz/10.)
  #endif not isame:

  if ishow: showplot()
#enddef plothull3d()

def nmerge(n1,n2,n12,vars1='',vars2='',vars12=''):

  idx1 = GetIndexN(n1)

  if idx1 == -1:
    if type(n1) != str:
      print("*** Error in nmerge: n1 must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #endeif idx1 == -1

  idx2 = GetIndexN(n2)

  if idx2 == -1:
    if type(n2) != str:
      print("*** Error in nmerge: n2 must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #endeif idx2 == -1

  n1 = Ntup[idx1]
  n2 = Ntup[idx2]

  if vars1 =='':
    vars1 = ""
    nc = n1.columns
    for iv in range(len(nc)-1):
      vars1 += nc[iv] + ':'
    #endfor iv in range(len(n1.columns)):
    vars1 += nc[len(nc)-1]
  #endif vars1 =='':

  if vars2 =='':
    vars2 = ""
    nc = n2.columns
    for iv in range(len(nc)-1):
      vars2 += nc[iv] + ':'
    #endfor iv in range(len(n2.columns)):
    vars2 += nc[len(nc)-1]
  #endif vars2 =='':

  if vars12 == '': vars12 = vars1 + ":" + vars2

  n12 = ncre(n12,n12,vars12)

  varl1 = nlistcolon(vars1)
  varl2 = nlistcolon(vars2)
  varl12 = nlistcolon(vars12)

  for iv in range(len(varl1)):
    v = varl1[iv]
    w = varl12[iv]
    n12[w] = n1[v]
  #endfor iv in range(len(varl1))

  iw = len(varl1)

  for iv in range(len(varl2)):
    v = varl2[iv]
    w = varl12[iv+iw]
    n12[w] = n2[v]
  #endfor iv in range(len(varl1))

  nupdate_header(n12)

  return n12
#enddef nmerge(...)

def bundu(gl,a=3.598,b=-3.840,c=0.631):
  # a=3.598,b=-3.840,c=0.631) is Beff for 77K
  # a=3.694,b=-5.068,c=1.520) is Bmax !! for 300K
  return a*exp(b*gl+c*gl*gl)
#enddef bundu()

def mn(): setmarkersize(6.)
def mm(): setmarkersize(4.)
def ms(): setmarkersize(2.)

def ncoph():
  print("\n *** Use nproj1 or nproj2 *** \n")
#enddef ncoph()

def legend(lege='',loc='best',x='',y=''):

  global Legend

  if len(lege):
    Legend.append(lege)
    return
  #endif

  if type(lege) == str:
    lege = Legend
  #endif type(lege) == str

  if type(lege) != list:
    print("\n*** Error in legend(lege): Legend must be a list ***\n")
    return
  #endif type(lege) != list
  #reakpoint()
  if type(x) != str and type(y) != str:
    plt.legend(lege,bbox_to_anchor=(x,y))
  elif type(x) != str:
    plt.legend(lege,bbox_to_anchor=(x,0.95))
  elif type(y) != str:
    plt.legend(lege,bbox_to_anchor=(0.95,y))
  else:
    plt.legend(lege,loc=loc)
  #endif

  showplot()
#def legend(lege)

def wans(text="Q or q to quit:"):
  if os.path.exists(".mhbook_quit"): Quit("Terminated, found file .mhbook_quit!")
  print("\n",text)
  ans = input()
  if ans.lower() == 'q': Quit()
#enddef wans(text="Hit return to continue:")

def nfun(ntup='?',fun='x', xvals='default', ioverwrite=1):

  if type(ntup) == str and ntup == '?' or type(fun) != str or not fun:
    print("\nUsage: nfun(ntupo, fun='x', xvals='default',ioverwrite=1) creates new Ntuple with function.\n")
    return
  #endif type(ntup) == str and ntup == '?'

  if type(ntup) != str:
    print("*** Error in nfun: Ntuple must be of type string ***")
    return -1
  #endif type(ntup) != str

  if type(fun) != str:
    print("*** Error in nfun: Function must be of type string ***")
    return -1
  #endif type(fun) != str

  nt = ncre(ntup,fun,'x:y')

  xv = pd.DataFrame(xvals,columns=['x'])
  nt.x = xv.x
  x = nt.x

  exec('nt.y = ' + fun)

  nupdate_header(ntup)

  return nt

#enddef nfun(ntup='?',fun='x', xvals='default', ioverwrite=1)

def hdump(hist='?',filh='hdump.dat'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(hist) == str and hist == '?':
    dump("\nUsage: hdump(Histo): dumps histogram Histo.\n")
    return
  #endif type(h) == str and h == '?'

  idx = GetIndexH1(hist)

  if idx < 0 or Khdeleted:
    idx = GetIndexH2(hist)
    if idx >= 0 or Khdeleted:
      print("*** hdump for 2D-Histos not yet implemented ***")
      return
    #endif
    print(" Error in hdump: Non-existing histogram",hist)
    return
  #endif idx < 0 or Khdelete

  Fh = open(filh,"w")
  fwrite(Fh,"* x       y        ave       ey       y2        n")
  for i in range(len(H1h)):
    fwrite(Fh,H1h.x[i],H1h.y[i],H1h.ave[i],H1h.ey[i],H1h.y2[i],H1h.n[i])
  #endfor
  Fh.close()

#enddef hdump(hist='?'):

def hprint(hist='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(hist) == str and hist == '?':
    print("\nUsage: hprint(Histo): Prints histogram Histo.\n")
    return
  #endif type(h) == str and h == '?'

  idx = GetIndexH1(hist)

  if idx < 0 or Khdeleted:
    print(" Error in hprint: Non-existing histogram",hist)
    return
  #endif idx < 0 or Khdelete

  print(H1h)
#enddef hprint(hist='?'):

def hfun(hist='?',fun='x', nx=101, xmin=-0.5, xmax=100.5):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(hist) == str and hist == '?' or type(fun) != str or not fun:
    print("\nUsage: hfun(Histo, fun='x', nx=101, xmin=-0.5, xmax=100.5) fills new histogram Histo with function or replaces content of an existing histogram by the function\n")
    return
  #endif type(h) == str and h == '?'

  if Kecho:
    if type(hist) == str:
      print("hfun(hist='" + hist + "',fun='", fun, "', nx=", nx, ", xmin=",xmin, ", xmax=", xmax, ")")
    else:
      print("hfun(hist=", hist, ",fun='", fun, "', nx=", nx, ", xmin=",xmin, ", xmax=", xmax, ")")
  #endif Kecho

  idx = GetIndexH1(hist)

  if idx < 0 or Khdelete:
    dx = (xmax-xmin) / max(1,nx-1)
    hret = hbook1(hist,fun,nx,xmin-dx/2.,xmax+dx/2.,True)
  #endif idx < 0 or Khdelete

  h = hget(hist)

  x = h.x

  exec('h.y = ' + fun)

  h.y2 = h.y * h.y
  h.ave = h.y
  h.ey = 0.0
  h.n = 1.0

  h1header_update(hist)

  return
#enddef hfun(...)

def h2header_update(hist='?'):
  print('*** h2header_update not yet implemented ***')
#enddef h2header_update(hist='?')

def weednan(h):
  hwn = deepcopy(h)
  hwn[np.isnan(hwn)] = 0.0
  return hwn
#enddef

def h1header_update(hist='?'):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(hist) == str and hist == '?':
    print("\nUsage: h1header_upate(Histo): Updates header of histogram Histo.\n")
    return
  #endif type(h) == str and h == '?'

  idx = GetIndexH1(hist)

  if idx == -1:
    print("*** Non-existing histogram ***")
  #endif idx == -1

  h = H1h

  nx = len(h)

  xmin = h.x.min()
  xmax = h.x.max()
  dx = (xmax - xmin) / max(nx-1,1)
  xmin -= dx / 2.
  xmax += dx / 2.

  h.y = weednan(h.y)
  h.y2 = weednan(h.y2)
  h.n = weednan(h.n)

  h['ave'] = h.y/h.n
  #h.ave[np.isnan(h.ave)] = 0.0
  h['ave'] = weednan(h['ave'])

  h['ey'] = (h.y2/h.n-h.ave**2)**0.5
  #h.ey[np.isnan(h.ey)] = 0.0
  h['ey'] = weednan(h['ey'])

  H1h = h
  H1[idx] = H1h

  head1 = H1head[idx]

  head1[2] = nx
  head1[3] = xmin
  head1[4] = xmax

  head1[7] = min(h.y)
  head1[8] = max(h.y)
  head1[9] = h.n.sum()

  sumy = h.y.sum()
  ya = abs(h.y)
  sumya = ya.sum()

  xmean = 0.
  xrms = 0.
  xsqsum = 0.

  if sumya:
    xsqsum = ((h.x)**2*ya).sum()
    xmean = (h.x*h.y).sum()/sumya
    xrms = np.sqrt(max(0.0,xsqsum/sumya-xmean**2))
  #endif

  head1[10] = sumy

  head1[11] = xsqsum
  head1[12] = xmean
  head1[13] = xrms

  quer = 'x < ' + str(xmin)
  hund = h.query(quer)
  nunder = hund.n.sum()
  under = hund.y.sum()

  quer = 'x >= ' + str(xmax)
  hov = h.query(quer)
  nover = hov.n.sum()
  over = hov.y.sum()

  head1[14] = nunder
  head1[16] = under
  head1[15] = nover
  head1[17] = over

  head1[18] = h.x.sum()
  head1[19] = ((h.x)**2).sum()

#enddef h1header_update(hist='?')

def hreset(h):

  idx = GetIndexH1(h)

  if idx == -1:
    idx = GetIndexH2(h)
    if idx == -1:
      print("*** hrest: Non-existing histogram ***")
      return
    else:
      h2reset(h)
      return
    #endif
  else:
    h1reset(h)
    return
  #endif

  print("*** hrest: Non-existing histogram ***")

#enddef hreset(h)

def h1reset(h):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  idh = GetIndexH1(h)

  if idh < 0:
      print("*** Error in hreset(Histo): Non-existing histogram")
      return
  #endif

  H1h.y = H1h.x*0.0
  H1h.y2 = H1h.x*0.0
  H1h.n = H1h.x*0.0
  H1h.ave = H1h.x*0.0
  H1h.ey = H1h.x*0.0

  h1header_update(h)

#enddef h1reset(h='?')

def hdelete(h='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(h) == str and h == '?':
    print("\nUsage: hdelete(Histo)\nThe histogram is not really deleted but marked as.")
    return
  #endif type(h) == str and h == '?'

  idh = GetIndexH1(h)

  if idh < 0:
    idh = GetIndexH2(h)
    if idh < 0:
      print("*** Error in hdelete(Histo): Non-existing histogram")
      return -1
    else:
      H2hh[38] = 1
    #endif idh < 0
  else:
    H1hh[20] = 1
  #if idh < 0

  Khdeleted = 1

#enddef hdelete(h='?')

def hmin(h='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(h) == str and h == '?':
    print("\nUsage: hmin(Histo)")
    return
  #endif type(h) == str and h == '?'

  idh = GetIndexH1(h)

  if idh < 0:
    idh = GetIndexH2(h)
    if idh < 0:
      print("*** Error in hmin(Histo): Non-existing histogram")
      return -1
    else:
      return H2hh[11]
    #endif idh < 0
  else:
    return H1hh[7]
  #if idh < 0

#enddef hmin(h='?')

def hmax(h='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(h) == str and h == '?':
    print("\nUsage: hmax(Histo)")
    return
  #endif type(h) == str and h == '?'

  idh = GetIndexH1(h)

  if idh < 0:
    idh = GetIndexH2(h)
    if idh < 0:
      print("*** Error in hmax(Histo): Non-existing histogram")
      return -1
    else:
      return H2hh[12]
    #endif idh < 0
  else:
    return H1hh[8]
  #if idh < 0

#enddef hmax(h='?')

def printplopt():
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  if type(Klegend) == int or type(Klegend) == bool:
    iledg = Klegend
    Klegend = IntVar(get_master())
    Klegend.set(iledg)
  #if type(Klegend) == int or type(Klegend) == bool

  if type(Icmap) == int:
    icmap = Icmap
    Icmap = IntVar(get_master())
    Icmap.set(icmap)
  #endif type(Icmap) == int

  if type(IsameGlobal) == int:
    isame = IsameGlobal
    IsameGlobal = IntVar(get_master())
    IsameGlobal.set(isame)
  #endif type(IsameGlobal) == int

  print("boxes; Iboxes:",Iboxes)
  print("cont3d; Icont3d:",Icont3d)
  print("err, E; Ierr:",Ierr)
  print("fill1d, F; Ifill1d:",Ifill1d)
  print("hist, H; Ihist:",Ihist)
  print("inter; Iinter:",Iinter)
  print("line, L; Iline:",Iline)
  print("marker, M; Imarker:",Imarker)
  print("noempty, N; Inoempty:",Inoempty)
  print("prof, P; Iprof:",Ifill1d)
  print("same, S; Isame, IsameGlobal:",Isame,IsameGlobal.get())
  print("scat3d; Iscat3d:",Iscat3d)
  print("scatter; Iscatterf:",Iscatter)
  print("spline, C; Ispline:",Ispline)
  print("surf; Isurf:",Isurf)
  print("trisurf; Itrisurf:",Itrisurf)
  print("closed; Iclosed:",Iclosed)
#enddef printplopt()

def plotoptions(plopt=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  Iplotopt = 0

  if type(plopt) != str: return

  if plopt.lower() == '2d' or plopt == '': plopt = Mode2d
  elif plopt.lower() == '3d' or plopt == '': plopt = Mode3d
  #if re.search('2d',plopt.lower()): plopt = re.sub('2[dD]',Mode2d,plopt)
  #elif re.search('3d',plopt.lower()): plopt = re.sub('3[dD]',Mode3d,plopt)

  Ihist = 0
  Iprof = 0
  Ispline = 0
  Iboxes = 0
  Isame = 0
  Iscatter = 0
  Iscat3d = 0
  Icont3d = 0
  Itrisurf = 0
  Ierr = 0
  Isurf = 0
  Iline = 0
  Iinter = 0
  Ifill1d = 0
  Imarker = 0
  Inoempty = 0
  Iclosed = 0

  if not Fig: window()

  if type(Klegend) == int or type(Klegend) == bool:
    iledg = Klegend
    Klegend = IntVar(get_master())
    Klegend.set(iledg)
  #endif type(Klegend) == int or type(Klegend) == bool

  if type(Icmap) == int:
    icmap = Icmap
    Icmap = IntVar(get_master())
    Icmap.set(icmap)
  #endif type(Icmap) == int

  if type(IsameGlobal) == int:
    isame = IsameGlobal
    IsameGlobal = IntVar(get_master())
    IsameGlobal.set(isame)
  #endif type(IsameGlobal) == int

  if not Kplots[Kzone-1]:
    if re.search('same',plopt): plopt = re.sub("same","",plopt)
    if re.search('S',plopt): plopt = re.sub("S","",plopt)
  #endif

  if re.search('marker',plopt) or re.search('M',plopt) or re.search('P',plopt):
    Iplotopt = 1;  Imarker = 1
  if re.search('spline',plopt) or re.search('C',plopt):
    Iplotopt = 1;  Ispline = 1; plopt = re.sub("spline","",plopt)
  if re.search('same',plopt) or re.search('S',plopt) or IsameGlobal.get() == 1:
    Iplotopt = 1;  Isame = 1
  if re.search('scatter',plopt):
    Iplotopt = 1;  Iscatter = 1

  if re.search('boxes',plopt): Iplotopt = 1;  Iboxes = 1

  if not Iboxes and re.search('scat3d',plopt): Iplotopt = 1;  Iscat3d = 1
  if not Iboxes and re.search('cont3d',plopt): Iplotopt = 1;  Icont3d = 1
  if not Iboxes and re.search('trisurf',plopt): Iplotopt = 1;  Itrisurf = 1
  if not Iboxes and re.search('inter',plopt): Iplotopt = 1;  Iinter = 1
  if not Iboxes and re.search('hist',plopt) or re.search('H',plopt): Iplotopt = 1;  Ihist = 1

  if re.search('err',plopt) or re.search('E',plopt): Iplotopt = 1;  Ierr = 1
  if Itrisurf == 0 and re.search('surf',plopt): Iplotopt = 1;  Isurf = 1
  if re.search('closed',plopt): Iplotopt = 1;  Iclosed = 1
  if re.search('line',plopt) or re.search('L',plopt): Iplotopt = 1;  Iline = 1
  if re.search('fill1d',plopt) or re.search('F',plopt): Iplotopt = 1;  Ifill1d = 1
  if re.search('prof',plopt) or re.search('P',plopt): Iplotopt = 1;  Iprof = 1
  if re.search('sprof',plopt) or re.search('profs',plopt) or re.search('spread',plopt): Iplotopt = 1;  Iprof = 1; Ierr = 1
  if re.search('noempty',plopt) or re.search('N',plopt): Iplotopt = 1;  Inoempty = 1

  if not Linestyle or Linestyle == 'none' or Linewidth <= 0.: Iline = 0
  if not Markertype or Markertype == 'none' or Markersize <= 0.: Imarker = 0

#enddef plotoptions(plopt=''):

def help(key=''):

  keys = []
  keys.append("hplot1d") #1
  keys.append("plopt") #2

  ikey = -1

  if type(key) == int:
    ikey = key
  else:
    for k in range(len(keys)):
      if key == keys[k]:
        ikey = k + 1
        break
    #for k in range(len(keys))
  #endif type(key) == int

  if ikey < 0:
    print("\nProvide number of keyword for help; keys are:\n")
    for k in range(len(keys)):
      print(k+1,":",keys[k])
    #endfor k in range(len(keys))
    print("\nNot found what you are looking for? Try the desired command with no argument...")
  elif ikey == 1:
    hplot1d()
    print("See also help('plopt')")
  elif ikey == 2:
    printplopt()
    print("\nMode2d and Mode3d contain the default options for the plot commands.")
    print("\nThey can be set or got by \n setmode2d(m='!') and getmode2d() etc..")
  #endif ikey < 0

#def help(key='')

def setnticksx(nx=0):
  global NxBinMax
  NxBinMax = nx
#enddef setnticksx(nx=0)

def vfft(vx,vy):
  om = np.fft.fftfreq(len(vx),vx[1]-vx[0])
  fft = np.fft.fft(vy)
  return om, fft
#def vfft(vx,vy)

def reset_zoom():
  global ZoomXmin,ZoomXmax,ZoomYmin,ZoomYmax,NxBinMax
  ZoomXmin = -1.e30
  ZoomXmax=1.e30
  ZoomYmin = -1.e30
  ZoomYmax=1.e30
#def reset_zoom()

def mhb_mkdir(chdir='!'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  #if Kecho: print("\nmhb_mkdir::chdir, Kdir, Ndir:",chdir,Kdir,Ndir)

  if chdir == '!': chdir = "Hdir" + str(Ndir-1)

  ifound = 0
  for i in range(Ndir):
    #if Kecho: print("\nmhb_mkdir::i, Hdir[i][0]:",i,Hdir[i][0])
    if Hdir[i][0] == chdir:
      ifound = i
      chdir += '_' + str(Ndir)
      break
    #endif
  #endfor

  Kdir = Ndir
  Ndir += 1

  Nh1 = 0
  Nh2 = 0
  Nntup =0
  Nnctup = 0

  #Cdir = Hdir[Kdir][0]

  #Nh1 = Hdir[Kdir][2]
  #Nh2 = Hdir[Kdir][3]
  #Nntup = Hdir[Kdir][4]
  #Nnctup = Hdir[Kdir][5]

  #H1head = Hdir[Kdir][6]
  #H1 = Hdir[Kdir][7]
  #H2head = Hdir[Kdir][8]
  #H2 = Hdir[Kdir][9]

  #Nhead = Hdir[Kdir][10]
  #Ntup = Hdir[Kdir][11]
  #Nctup = Hdir[Kdir][12]
  #Nind = Hdir[Kdir][13]

  #             0     1    2   3   4      5    6  7  8  9  10 11 12 13 14 15 16
  Hdir.append([chdir,Ndir,Nh1,Nh2,Nntup,Nnctup,[],[],[],[],[],[],[],[],[],[],[]])
#enddef mhb_mkdir(chdir='!')

def mhb_ldir():
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  l = 0
  for k in range(Ndir):
    if len(Wfiles[k]) > l:
      l = len(Wfiles[k])
    #endif
  #endfor

  l += 10
  sf = "{0:<" + str(l) + "} {1:<20}"
  for k in range(Ndir): print(k," ",sf.format(Wfiles[k],Hdir[k][0]))

#enddef mhb_ldir

def mhb_pwd(isilent=0,iretval=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if not isilent: print(Wfiles[Kdir],NL,Hdir[Kdir][0])
  if iretval: return Wfiles[Kdir],Hdir[Kdir][0]

#enddef mhb_pwd

def mhb_cd_up(): mhb_cd('+')
def mhb_cd_down(): mhb_cd('-')

def mhb_cd(cdir='!'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  #if Kecho:
  #  print("mhb_cd::cdir:",cdir)
  #  print("mhb_cd::Kdir,Ndir:",Kdir,Ndir)
  #endif

  kdir = -1

  #print("\nCreating internal N-tuple directory: "+cdir)
  if type(cdir) == str:
    if cdir == '!':
      Kdir = Ndir - 1
      kdir = Kdir
    elif cdir == '-':
      Kdir -= 1
      if Kdir < 0: Kdir = Ndir - 1
      kdir = Kdir
    elif cdir == '+':
      Kdir += 1
      if Kdir >= Ndir: Kdir = 0
      kdir = Kdir
    else:
      for i in range(Ndir):
        #print("mhb_cd::i,Hdir[i][0]:",Hdir[i][0])
        if Hdir[i][0] == cdir:
          #if Kecho: print("mhb_cd::i,Hdir[i][0]:",Hdir[i][0])
          kdir = i
          break
  elif type(cdir) == int and cdir >= 0 and cdir <= Ndir:
    kdir = cdir
  #endif

  if kdir == -1:
    print('\n*** Warning in mhb_cd: Unknown Ntuple directory ***\n')
  else:

    Kdir = kdir
    Cdir = Hdir[Kdir][0]
    Fdir = Wfiles[Kdir]

    Nh1 = Hdir[Kdir][2]
    Nh2 = Hdir[Kdir][3]
    Nntup = Hdir[Kdir][4]
    Nnctup = Hdir[Kdir][5]

    H1head = Hdir[Kdir][6]
    H1 = Hdir[Kdir][7]
    H2head = Hdir[Kdir][8]
    H2 = Hdir[Kdir][9]

    Nhead = Hdir[Kdir][10]
    Ntup = Hdir[Kdir][11]
    Nctup = Hdir[Kdir][12]
    Nind = Hdir[Kdir][13]
    Ncind = Hdir[Kdir][14]

    H1ind = Hdir[Kdir][15]
    H2ind = Hdir[Kdir][16]

    global Wrun, Krun
    Wrun = int(H1[0].y[1])

    print("\nCWD:", Kdir,"   ",Fdir,"    ",Cdir,"    Run:", Wrun)

    if Krun and Trun:
      Trun.set_text("Run "+str(Wrun))
      Trun.set_visible(True)
    #endif

#enddef mhb_cd(cdir='!')

def zoom(xmin,xmax,ymin=-1.2345e30,ymax=1.2345e30):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  yn, yx = Ax.get_ylim()
  if ymin == -1.2345e30: ymn=yn
  else: ymn = ymin
  yn, yx = Ax.get_ylim()
  if ymax == 1.2345e30: ymx=yx
  else: ymx = ymax
  Ax.axis([xmin,xmax,ymn,ymx])
  ZoomXmin=xmin
  ZoomXmax=xmax
  ZoomYmin=ymin
  ZoomYmax=ymax
  showplot()

def setplotsize(fig='',w='A4',h='landscape',pad=2.):

  if type(fig) == str: fig = plt.gcf()

  if type(w) == str and w.upper() == 'A4':
    if type(h) == str and ( h.lower() == 'l' or h.lower() == 'landscape'):
      w = 29.7-2.*pad
      h = 21.0-2.*pad
    else:
      h = 29.7-2.*pad
      w = 21.0-2.*pad
    #endif type(h) == str and ( h.lower() == 'l' or h.lower() == 'landscape'
  #endif type(w) == str and w.upper() == 'A4'

  fig.set_figwidth(w/2.54)
  fig.set_figheight(h/2.54)

#enddef setplotsize():

def pplot(pname="WavePlot.pdf",w=0,h=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  Fig = plt.gcf()
  wo,ho = getplotsize()

  if type(w) ==  float and type(h) == float and w*h != 0.0 or type(w)== str:
    setplotsize(Fig,w,h)
  #endif type(w) ==  float and type(h) == float and w*h != 0.0 or type(w)== str

  try:
    Fig.savefig(pname)
    print("\nFigure written to ",pname)
  except:
    print("\nCould not write PDF-Dokument, try another format... ")
  #endtry

  setplotsize(Fig,wo,ho)

#enddef pplot("pname")

def h1pack(idh='?', data=None):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == str and idh == '?':
    print("\nUsage: h1pack(idh, data) To be tested!")
    return

  idx = GetIndexH1(idh)

  if idx == -1:
    print("*** Error in h1pack: Non-existing histogram")
    return -1

  if type(data) == list: data = pd.DataFrame(data,columns=['y'])
  else: data.columns = ['y']

  npack = ncre("npack","npack","x:y",ioverwrite=1)
  npack = nfill("npack",data['y'])

  nproj1("npack","x:y", idh=idx)

#enddef h1pack(idh='?', )

def hcopn(idh='?', nt='', varlis='x:y:ey', ntit='!',kweedzero=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == str and idh == '?':
    print("\nUsage: nt = hcopn(idh, nt,varlis, ntit='!',kweedzero=1)")
    return
  #endif

  is2d = 0
  idx = GetIndexH1(idh)

  if idx == -1:
    idx = GetIndexH1(idh)
    is2d = 1
  #endif

  if idx == -1:
    print("*** Error in hcopn: Non-existing histogram")
    return -1
  #endif

  varl = nlistcolon(varlis)

  if not is2d:

    if ntit == '!': ntit = nt

    nt = ncre(ntit,ntit,varlis,1)
    nt[varl[0]] = H1[idx]['x']

    if len(varl) > 2:
      if varl[1]:
        nt[varl[1]] = H1[idx]['y']
      else:
        nt[varl[1]] = H1[idx]['x']*0.0
      #endif
      if varl[2]:
        nt[varl[2]] = H1[idx]['ey']
      else:
        nt[varl[2]] = H1[idx]['x']*0.0
      #endif
    elif len(varl) > 1:
      if varl[1]:
        nt[varl[1]] = H1[idx]['y']
      else:
        nt[varl[1]] = H1[idx]['x']*0.0
      #endif
    #endif

  else:

    if ntit == '!':  ntit = nt
    nt = ncre(ntit,ntit,varlis,1)
    nt[varl[0]] = H2[idx]['x']
    if varl[1]:
      nt[varl[1]] = H2[idx]['y']
    else:
      nt[varl[1]] = H2[idx]['x']*0.0
    if varl[2]:
      nt[varl[2]] = H2[idx]['z']
    else:
      nt[varl[3]] = H2[idx]['x']*0.0
    if varl[3]:
      nt[varl[3]] = H2[idx]['ez']
    else:
      nt[varl[3]] = H2[idx]['x']*0.0

  #endif not is2d

  nupdate_header(nt)

  if kweedzero:
    nt = ncopn(nt,ntit,varlis='',select='y != 0 or ey != 0',ioverwrite=1)

  return nt

#enddef hcopn(idh='?', nt='', varlis='')

def nrandom(nt='?',varlis='', n=100, mode='u', iplot=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: nt = nrandom(nt,varlis, n=100, mode='u', iplot=0)")
    return

  varl = nlistcolon(varlis)
  dim = Ncolon + 1

  if mode == 'u':
    points = np.random.rand(n,dim)
    name = nt + "_" + str(n) + "_uniform_" + str(dim)
  else:
    points = np.random.randn(n,dim)
    name = nt + "_" + str(n) + "_normal_" + str(dim)
  #endif modu == 'u'

  nt = ncre(nt,name,varlis)
  nt = nfill(nt,points)

  if iplot and dim < 4:
    nplot(nt,varlis)
  #endif iplot

  return nt
#enddef nrandom(nt='?',varlis='', mode='u', dim=1, iplot=0)

def nhull2d(nt='?',varlis='',select='', iplot=1, iretval=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: hull = nhull2d(nt,varlis,select,iplot=1)")
    return
  #endif

  idn = GetIndexN(nt)

  if idn == -1:
    if type(nt) != str:
      print("*** Error in nhull2d: Nt must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #if type(nt) != str

  varl = nlistcolon(varlis)
  if Ncolon != 1:
    print("Bad list of variables: Must contain exactly tow, but is",varlis)
    return -1
  #endif

  nt = Ntup[idn]

  if len(select):
    N = nt.query(select)
    Nsel = N
  else:
    Nsel = N
  #endif len(select)

  vx, vy = ncopv(Nsel,varlis,select)

  vxmin = vx.min()
  vxmax = vx.max()
  vymin = vy.min()
  vymax = vy.max()

  dvx = vxmax - vxmin
  if dvx:
    wx = (vx - vxmin) / dvx
  else:
    wx = vx
  #endif

  dvy = vymax - vymin
  if dvy:
    wy = (vy - vymin) / dvy
  else:
    wy = vy
  #endif

  points = np.array([wx,wy]).T

  nhull = ncre("Nhull2d","Nhull2d","ipoi:iedge:x:y",ioverwrite=1)
  fill = []
  nedge = 0

  try:

    hull = ConvexHull(points)
    ms = len(wx)

    s = []

    for m in range(ms):
      s.append(list(hull.simplices[m]))
    #endfor

    sc = []
    sc.append(s[0])

    kl = s[0][0]
    il = s[0][1]

    for i in range(1,ms):
      for k in range(1,ms):
        i1 = s[k][0]
        i2 = s[k][1]
        if i1 == il:
          sc.append(s[k])
          il = i2
        elif i2 == il:
          sw = s[k]
          sc.append([sw[1],sw[0]])
          il = i1
        #endif
      #endfor
      if il == kl: break
    #endfor

    for ed in sc:
      nedge += 1
      i1 = ed[0]
      i2 = ed[1]
      #x1 = points[i1][0]
      #y1 = points[i1][1]
      #x2 = points[i2][0]
      #y2 = points[i2][1]
      x1 = vx[i1]
      y1 = vy[i1]
      x2 = vx[i2]
      y2 = vy[i2]
      fill.append([i1,nedge,x1,y1])
      fill.append([i2,nedge,x2,y2])
    #endfor

    fill = np.array(fill)
    nhull = nfill(nhull,fill)

    if iplot:
      if Isame:
        nplot(nhull,"x:y",plopt='sameline')
      else:
        nplot(nhull,"x:y",plopt='line')
      #endif
    #endif iplot

  except:

    npoi = len(wx)

    if dvx == 0 and dvy == 0:

      fill.append([0,0,vx[0],vy[0]])
      fill.append([0,0,vx[0],vy[0]])
      fill = np.array(fill)
      nhull = nfill(nhull,fill)

      if iplot:
        if Isame:
          nplot(nhull,"x:y",plopt='samemarker')
        else:
          nplot(nhull,"x:y",plopt='marker')
        #endif
      #endif iplot

    else:

      if dvx < dvy:
        ww =wy
      else:
        ww =wx
      #endif
      vmn = ww.min()
      for imin in range(npoi):
        if ww[imin] == vmn:
          break
        #endif
      #endfor
      vmx = ww.max()
      for imax in range(npoi):
        if ww[imax] == vmx:
          break
        #endif
      #endfor

    #endif

      fill.append([0,1,vx[imin],vy[imin]])
      fill.append([0,1,vx[imax],vy[imax]])
      fill = np.array(fill)
      nhull = nfill(nhull,fill)

      if iplot:
        if Isame:
          nplot(nhull,"x:y",plopt='samemarker')
        else:
          nplot(nhull,"x:y",plopt='marker')
        #endif
      #endif iplot

  #endtry

  if iretval:
    return nhull
  else:
    return
  #endif

#enddef nhull2d(nt='?')

def vhull2d(vx,vy,varlis='',iplot=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz



  vxmin = vx.min()
  vxmax = vx.max()
  vymin = vy.min()
  vymax = vy.max()

  dvx = vxmax - vxmin
  wx = (vx - vxmin) / dvx
  dvy = vymax - vymin
  wy = (vy - vymin) / dvy

  points = np.array([wx,wy]).T

  hull = ConvexHull(points)

  n=-1
  wx = []
  wy = []

  for i in hull.vertices:
    n += 1
    wx.append(points[n][0])
    wy.append(points[n][1])
    wx.append(vx[n])
    wy.append(vy[n])
  #endfor s in hull.simplices

  if iplot:
    vplxy(wx,vy)
  #endif iplot

  return wx,wy

#enddef vhull2d(nt='?')

def hull3d(points):
#+seq,mhbglobind.
  global Hull3D

  xyz = pd.DataFrame(points)
  xyz.columns = ['x','y','z']

  vx = xyz.x
  vy = xyz.y
  vz = xyz.z

  vxmin = vx.min()
  vxmax = vx.max()
  vymin = vy.min()
  vymax = vy.max()
  vzmin = vz.min()
  vzmax = vz.max()

  dvx = vxmax - vxmin
  wx = (vx - vxmin) / dvx
  dvy = vymax - vymin
  wy = (vy - vymin) / dvy
  dvz = vzmax - vzmin
  wz = (vz - vzmin) / dvz

  points = np.array([wx,wy,wz]).T

  hull = ConvexHull(points)

  # Find planes

  ntria = len(hull.simplices)
  done = np.zeros([ntria],dtype=int)

  nplan = 0

  iplanes = []
  planes = []

  for i in range(ntria):

    if done[i]: continue

    iplan = []
    plan = []

    s1 = hull.simplices[i]

    i1 = s1[0]
    i2 = s1[1]
    i3 = s1[2]

    x1 = wx[i1]
    x2 = wx[i2]
    x3 = wx[i3]
    y1 = wy[i1]
    y2 = wy[i2]
    y3 = wy[i3]
    z1 = wz[i1]
    z2 = wz[i2]
    z3 = wz[i3]

    iplan.append([i1,i2,i3])
    plan.append([vx[i3],vy[i3],vz[i3]])
    plan.append([vx[i1],vy[i1],vz[i1]])
    plan.append([vx[i2],vy[i2],vz[i2]])
    plan.append([vx[i3],vy[i3],vz[i3]])

    v = [x2-x1,y2-y1,z2-z1]
    w = [x3-x1,y3-y1,z3-z1]
    vn = [ v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0] ]
    vn = vn / (vn[0]*vn[0] + vn[1]*vn[1] + vn[2]*vn[2])**0.5

    done[i] = 1

    j=i+1

    while j < ntria:

      s2 = hull.simplices[j]

      ihit = 0
      for k in [0,1,2]:

        l = s2[k]

        xx = wx[l]
        yy = wy[l]
        zz = wz[l]

        dr = [xx-x1,yy-y1,zz-z1]
        dist = dr[0]*vn[0] + dr[1]*vn[1] + dr[2]*vn[2]

        if abs(dist) < 1.e-9: ihit += 1

      #endfor k in [0,1,2]

      if ihit == 3:

        done[j] = 1

        i1 = s2[0]
        i2 = s2[1]
        i3 = s2[2]

        iplan.append([i1,i2,i3])
        plan.append([vx[i1],vy[i1],vz[i1]])
        plan.append([vx[i2],vy[i2],vz[i2]])
        plan.append([vx[i3],vy[i3],vz[i3]])

      #endif ihit == 3:

      j += 1
    #endwhile j <= ntria:

    nplan += 1

    vnmax = abs(vn).max()
    for iax in [0,1,2]:
      if abs(vn[iax]) == vnmax:
        break
    #endfor iax in [0,1,2]

    plan2d = []
    if iax == 2:
      for ip in range(len(plan)):
        plan2d.append([plan[ip][0],plan[ip][1]])
    elif iax == 1:
      for ip in range(len(plan)):
        plan2d.append([plan[ip][0],plan[ip][2]])
    elif iax == 0:
      for ip in range(len(plan)):
        plan2d.append([plan[ip][1],plan[ip][2]])
    #endif i == 0

    hull2d = ConvexHull(plan2d)

    ip2d = []
    p2d = []
    for iv in hull2d.vertices:
      ip2d.append(iv)
      p2d.append(plan[iv])
    #endfor iv in hull2d.vertices

    ip2d.append(ip2d[0])
    p2d.append(p2d[0])

    iplanes.append(ip2d)
    planes.append(p2d)

  #endfor i in range(ntria):

  corns = []
  Hull3D = []

  for ipl in range(nplan):
    plan = planes[ipl]
    Hull3D.append(plan)
    iplan = iplanes[ipl]
    for ip in range(len(plan)):
      ipoi = iplan[ip]
      x = plan[ip][0]
      y = plan[ip][1]
      z = plan[ip][2]
      corns.append([ipoi+1,ipl+1,x,y,z])
    #endfor ip in range(len(plan))
  #endfor ipl in range(nplan)

  corns = pd.DataFrame(corns)
  corns.columns = ['ipoi','iplan','x','y','z']

  return corns

#enddef hull3d(points)

def mhull3d(nt='?',varlis='',select='',isame=0,
            facecolor='blue',edgecolor='black', alpha=0.5, iplot=1):

  global Hull3D

  if type(nt) == str and nt == '?':
    print("\nUsage: vertices = mhull3d(nt='?',varlis='',select='', isame=0, facecolor='blue',edgecolor='black', alpha=0.5, iplot=1)")
    return
  #endif type(nt) == str and nt == '?'

  idn = GetIndexN(nt)

  if idn == -1:
    if type(nt) != str:
      print("*** Error in nhull3d: Nt must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #if type(nt) != str

  vx, vy, vz = ncopv(nt,varlis,select)
  points = np.array([vx,vy,vz]).T

  vertices =  hull3d(points)

  if iplot:
    plothull3d(isame,facecolor=facecolor,edgecolor=edgecolor,alpha=alpha)

  return vertices

#enddef mhull3d(nt='?',varlis='',select='', plopt='',iplot=1)

def nhull3d(nt='?',varlis='',select='', plopt='',iplot=1, iretval=0,color='!',
           mcolor='!'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: hull = nhull3d(nt,varlis,select,iplot=1)")
    return
  #endif type(nt) == str and nt == '?'

  if type(nt) == Tdf:
    pass
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      print("*** Error in nhull3d: Unknown Ntuple ***")
      return -1
    #endif
    nt = Ntup[ind]
  elif type(nt) == int and nt > 0:
    nt = Ntup[idn]
  else:
    print("*** Error in nhull3d: Unknown Ntuple ***")
    return -2
  #endif nt >= 0:

  if len(select):
    N = nt.query(select)
    nt = N
  #endif len(select)

  if not len(nt):
    print("*** Error in nhull3d: No data, check ntuple and selection ***")
    return -1
  #endif

  if color != '!':
    setlinecolor(color)
  #endif

  if mcolor != '!':
    setlinecolor(mcolor)
  #endif

  varl = nlistcolon(varlis)

  sx = eval(nparse(nt,varl[0]))
  sy = eval(nparse(nt,varl[1]))
  sz = eval(nparse(nt,varl[2]))

  ntd = ncre("ntd","ntd","x:y:z",1)
  var = nlistcolon(varlis)
  #ntd.x = nt[var[0]]
  #ntd.y = nt[var[1]]
  #ntd.z = nt[var[2]]
  ntd.x = sx
  ntd.y = sy
  ntd.z = sz
  ntd = ntd.drop_duplicates()
  ntd.index = range(len(ntd))
  #nupdate_header("ntd")
  #nt = ncopv(ntd,varlis)

  vx = ntd.x
  vy = ntd.y
  vz = ntd.z

  vxmin = vx.min()
  vxmax = vx.max()
  vymin = vy.min()
  vymax = vy.max()
  vzmin = vz.min()
  vzmax = vz.max()

  dvx = vxmax - vxmin
  if dvx:
    wx = (vx - vxmin) / dvx
  else:
    wx = vx
  #endif
  dvy = vymax - vymin
  if dvy:
    wy = (vy - vymin) / dvy
  else:
    wy = vy
  #endif
  dvz = vzmax - vzmin
  if dvz:
    wz = (vz - vzmin) / dvz
  else:
    wz = vz
  #endif

  dvx = wx.max() - wx.min()
  dvy = wy.max() - wy.min()
  dvz = wz.max() - wz.min()

  points = np.array([wx,wy,wz]).T

  try:

    hull = ConvexHull(points)

  except:

    ntd2 = ncre("ntd2","ntd2","x:y",1)

    p1 = [wx[0],wy[0],wz[0]]
    p2 = [wx[1],wy[1],wz[1]]
    p3 = [wx[2],wy[2],wz[2]]

    vnor = vcross(listsub(p2,p1),listsub(p3,p1))
    vnor = listnorm(vnor)

    dvx = abs(vnor[0])
    dvy = abs(vnor[1])
    dvz = abs(vnor[2])

    dmx = max(dvx,dvy,dvz)

    if dmx == dvx:
      ntd2.x = wy
      ntd2.y = wz
      modez = 1
    elif dmx == dvy:
      ntd2.x = wx
      ntd2.y = wz
      modez = 2
    else:
      ntd2.x = wx
      ntd2.y = wy
      modez = 3
    #endif

    ntd2 = ntd2.drop_duplicates()
    ntd2.index = range(len(ntd2))
    nupdate_header("ntd2")
    nhull = nhull2d("ntd2","x:y","",0,iretval=1)

    if iplot:
      v1 = []
      v2 = []
      v3 = []
      i1 = int(nhull.ipoi[0])
      v1.append(vx[i1])
      v2.append(vy[i1])
      v3.append(vz[i1])
      for i in range(1,len(nhull),2):
        i1 = int(nhull.ipoi[i])
        v1.append(vx[i1])
        v2.append(vy[i1])
        v3.append(vz[i1])
      #endfor
      vplxyz(v1,v2,v3,'sameline')
    #endif

    if iretval:
      return nhull
    else:
      return
    #endif

  #entry

  if iplot:
    plotoptions(plopt)

  nhull = ncre("Nhull3d","Nhull3d","ipoi:iplan:x:y:z",ioverwrite=1)

  # Find planes

  ntria = len(hull.simplices)
  done = np.zeros([ntria],dtype=int)

  nplan = 0

  iplanes = []
  planes = []

  for i in range(ntria):

    if done[i]: continue

    iplan = []
    plan = []

    s1 = hull.simplices[i]

    i1 = s1[0]
    i2 = s1[1]
    i3 = s1[2]

    x1 = wx[i1]
    x2 = wx[i2]
    x3 = wx[i3]
    y1 = wy[i1]
    y2 = wy[i2]
    y3 = wy[i3]
    z1 = wz[i1]
    z2 = wz[i2]
    z3 = wz[i3]

    iplan.append([i1,i2,i3])
    plan.append([vx[i3],vy[i3],vz[i3]])
    plan.append([vx[i1],vy[i1],vz[i1]])
    plan.append([vx[i2],vy[i2],vz[i2]])
    plan.append([vx[i3],vy[i3],vz[i3]])

    v = [x2-x1,y2-y1,z2-z1]
    w = [x3-x1,y3-y1,z3-z1]
    vn = [ v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0] ]
    vn = vn / (vn[0]*vn[0] + vn[1]*vn[1] + vn[2]*vn[2])**0.5

    done[i] = 1

    j=i+1

    while j < ntria:

      s2 = hull.simplices[j]

      ihit = 0
      for k in [0,1,2]:

        l = s2[k]

        xx = wx[l]
        yy = wy[l]
        zz = wz[l]

        dr = [xx-x1,yy-y1,zz-z1]
        dist = dr[0]*vn[0] + dr[1]*vn[1] + dr[2]*vn[2]

        if abs(dist) < 1.e-9: ihit += 1

      #endfor k in [0,1,2]

      if ihit == 3:

        done[j] = 1

        i1 = s2[0]
        i2 = s2[1]
        i3 = s2[2]

        iplan.append([i1,i2,i3])
        plan.append([vx[i1],vy[i1],vz[i1]])
        plan.append([vx[i2],vy[i2],vz[i2]])
        plan.append([vx[i3],vy[i3],vz[i3]])

      #endif ihit == 3:

      j += 1
    #endwhile j <= ntria:

    nplan += 1

    vnmax = abs(vn).max()
    for iax in [0,1,2]:
      if abs(vn[iax]) == vnmax:
        break
    #endfor iax in [0,1,2]

    plan2d = []
    if iax == 2:
      for ip in range(len(plan)):
        plan2d.append([plan[ip][0],plan[ip][1]])
    elif iax == 1:
      for ip in range(len(plan)):
        plan2d.append([plan[ip][0],plan[ip][2]])
    elif iax == 0:
      for ip in range(len(plan)):
        plan2d.append([plan[ip][1],plan[ip][2]])
    #endif i == 0

    hull2d = ConvexHull(plan2d)

    ip2d = []
    p2d = []
    for iv in hull2d.vertices:
      ip2d.append(iv)
      p2d.append(plan[iv])
    #endfor iv in hull2d.vertices

    ip2d.append(ip2d[0])
    p2d.append(p2d[0])

    iplanes.append(ip2d)
    planes.append(p2d)

  #endfor i in range(ntria):

  fill = []

  Hull3D = []
  for ipl in range(nplan):
    plan = planes[ipl]
    Hull3D.append(plan)
    iplan = iplanes[ipl]
    for ip in range(len(plan)):
      ipoi = iplan[ip]
      x = plan[ip][0]
      y = plan[ip][1]
      z = plan[ip][2]
      fill.append([ipoi+1,ipl+1,x,y,z])
    #endfor ip in range(len(plan))
  #endfor ipl in range(nplan)

  fill = np.array(fill)
  nhull = nfill("Nhull3d",fill)

  if iplot:
    if Isame == 0:
      if Imarker: nplot(nhull,"x:y:z","","","marker")
      for ipl in range(nplan):
        sel = '"iplan==' + str(ipl+1) + '"'
        #if ipl: eval('nplot("Nhull3d","x:y:z",' + sel + ',"","sameline")')
        #else: eval('nplot("Nhull3d","x:y:z",' + sel + ',"","line")')
        if ipl: eval('nplot(nhull,"x:y:z",' + sel + ',"","sameline")')
        else: eval('nplot(nhull,"x:y:z",' + sel + ',"","line")')
      #endfor
    else:
      if Imarker: nplot(nhull,"x:y:z","","","samemarker")
      for ipl in range(nplan):
        sel = '"iplan==' + str(ipl+1) + '"'
        #eval('nplot("Nhull3d","x:y:z",' + sel + ',"","sameline")')
        eval('nplot(nhull,"x:y:z",' + sel + ',"","sameline")')
    #if Isame == 0
    #endfor ipl in range(nplan)
  #endif iplot

  if iretval: return nhull
  else: return

#enddef nhull3d(nt='?')

def nplothull3d(nt="Nhull3d",plopt=''):

  idn = GetIndexN(nt)

  if idn == -1:
    if type(nt) != str:
      print("*** Error in nplothull3d: Nt must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #if type(nt) != str

  nt = N
  nplan = int(nt.iplan.max())

  nplot(nt,"x:y:z","","",plopt)

  for ipl in range(nplan):
    sel = '"iplan==' + str(ipl+1) + '"'
    eval('nplot(nt,"x:y:z",' + sel + ',"","sameline")')
  #endfor ipl in range(nplan)

#enddef nplothull3d(nt="Nhull3d",plopt=''

def nappend(nt='?', nt2=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?' or len(nt2) == 0:
    print("\nUsage: Ntup1 = nappend(Ntup1, Ntup2)")
    return
  #endif

  id1 = GetIndexN(nt)
  id2 = GetIndexN(nt2)

  n1 = nget(id1)
  n2 = nget(id2)

  if id1 == -1:
    if type(nt) != str:
      print("*** Error in nappend: nt must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #if type(nt) != str

  nn = pd.concat([n1,n2])

  Ntup[id1] = nn
  nupdate_header(Ntup[id1],reindex=1)

  return Ntup[id1]
#enddef nappend(nt='?', nt2='')

def nfill(nt='?', data=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?' or len(data) == 0:
    print("\nUsage: Ntup = nfill(Ntup, data)")
    return
  #endif

  idn = GetIndexN(nt,1)

  if idn == -1:
    if type(nt) != str:
      print("*** Error in nfill: nt must be Ntuple or Ntuple name ***")
      return -1
    #endif type(nt) != str
  #if type(nt) != str

  nhead = Nhead[idn]
  nvar = nhead[3]

  varlis = []
  for k in range(nvar):
    l = 4 + k
    var = nhead[l][0]
    varlis.append(var)
  #endfor

  nt = N
  ndum = ncre("ndum","ndum",varlis,ioverwrite=1)

  if type(data) == list:
    dat = []
    dat.append(data)
    data = np.array(dat)
  #endif

  ndat = len(data)
  dat = np.array(data).T

  for k in range(nvar):
    l = 4 + k
    var = nhead[l][0]
    ndum[var] = dat[k]
  #endfor k in range(nvar):

  nn = pd.concat([nt,ndum])

  idn = GetIndexN(nt,1)
  Ntup[idn] = nn
  nupdate_header(Ntup[idn])

  ind = Ind

#  for i in range(nvar):
#    nhead[4+i][1] = nn[varlis[i]].min()
#    nhead[4+i][2] = nn[varlis[i]].max()
#  #endfor i in range(nvar):

#  nhead[4+nvar]=len(nn)

  return nn
#enddef nfill(nt='?', data='')

def npeaksabs(nt='?', varlis='', select='', pkmin=0.5,nsmooth=0,isilent=0,iretval=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: npeaksabs(Ntup,varlis='', select='', pkmin=0.5, nsmooth=0,isilent=0, iretval=0")
    return
  #endif type(nt) == str and nt == '?'

  if type(varlis) == str:
    if not len(varlis):
      print("*** Error in npeaksabs: No columns specified in var ***")
      return
    #endif not len(varlis)
    varlis = nlistcolon(varlis)
  #if type(varlis) == str

  if Ncolon != 1:
    print("*** Error in npeaksabs(...): Variable must contain exactly two terms, i.e. one colon ***")
    return -1
  #endif Ncolon != 1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  svar = ""
  svl = []

  nt = N
  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
  #endfor k in range(Ncolon)

  sv = nparse(nt,varlis[Ncolon])
  svar += nparse(nt,varlis[Ncolon])

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  Nlines = len(N)
  N.columns = ['x','y']
  N.index = range(Nlines)

  ixpeaks, xpeaks,ypeaks,sigma = \
  vpeaksabs(np.array(N.x),np.array(N.y),pkmin=0.5,nsmooth=0,isilent=0)

  if iretval:
    return ixpeaks, xpeaks,ypeaks,sigma
  else:
    return
  #endif iretval

#enddef npeaksabs(nt='?', select='', pkmin=0.5,nsmooth=0,isilent=0,iretval=0)

def hpeaks(h='?', select='', pkmin=0.5,nsmooth=0,isilent=0,iretval=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(h) == str and h == '?':
    print("\nUsage: hpeaks(H1, select='', pkmin=0.5, nsmooth=0,isilent=0, iretval=0")
    return
  #endif type(h) == str and h == '?'

  nhpeaks = hcopn(h, 'nhpeaks','x:y:ey',kweedzero=1)

  if iretval:
    ixpeaks, xpeaks, ypeaks, sigma = npeaks(nhpeaks,"x:y",select,pkmin,nsmooth,isilent,iretval)
    return ixpeaks, xpeaks, ypeaks, sigma
  else:
    npeaks(nhpeaks,"x:y",select,pkmin,nsmooth,isilent,iretval)
    return
  #endif iretval

#enddef hpeaks(nt='?', varlis='', select='', pkmin=0.5,nsmooth=0,isilent=0,ireval=0)

def npeaks(nt='?', varlis='', select='', pkmin=0.5,nsmooth=0,isilent=0,iretval=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: npeaks(Ntup,varlis='', select='', pkmin=0.5, nsmooth=0,isilent=0, iretval=0")
    return

  if type(varlis) == str:
    if not len(varlis):
      print("*** Error in npeaks: No columns specified in var ***")
      return
    varlis = nlistcolon(varlis)

  if Ncolon != 1:
    print("*** Error in npeaks(...): Variable must contain exactly two terms, i.e. one colon ***")
    return -1
  #endif

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  svar = ""
  svl = []

  nt=N
  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
  sv = nparse(nt,varlis[Ncolon])
  svar += nparse(nt,varlis[Ncolon])

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  Nlines = len(N)
  N.columns = ['x','y']
  N.index = range(Nlines)

  vx = np.array(N.x)
  vy = np.array(N.y)
  ixpeaks, xpeaks,ypeaks,sigma = vpeaks(vx,vy,pkmin=0.5,nsmooth=0,isilent=isilent)

  if iretval:
    return ixpeaks, xpeaks,ypeaks,sigma
  else:
    return
  #endif

#enddef npeaks(nt='?', varlis='', select='', pkmin=0.5,nsmooth=0,isilent=0,ireval=0)

def nstat(nt='?',var='',select='', iretval=1, isilent=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: nstat(Ntup,var='', select='', ireval=0, isilent=0")
    return
  #endif type(nt) == str and nt == '?'

  if type(var) == str:
    if not len(var):
      print("*** Error in nstat: No columns specified in var ***")
      return
    #endif not len(var):
    var = nlistcolon(var)
  #endif type(var) == str

  if Ncolon > 1:
    print("*** Error in nstat(...): Variable must contain one or two terms, i.e. no or one colon ***")
    return -1
  #endif Ncolon > 0:

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N

  svarx = nparse(nt,var[0])

  if len(var) > 1:

    svary = nparse(nt,var[1])
    scom = "global N; N = pd.DataFrame([" + svarx + "," + svary + "]).T"
    #print(scom)
    exec(scom)

    N.columns = ['x','y']
    N.index = range(len(N))

    res = vstat(N.x,N.y)

    Nmin = res[0]
    Nmax = res[1]
    Nmean = res[2]
    Nrms = res[3]
    Nxopt = res[4]
    Nyopt = res[5]

    if not isilent:
      print("\nNmin:",Nmin)
      print("Nmax:",Nmax)
      print("Nmean:",Nmean)
      print("Nrms:",Nrms)
      print("Nxopt:",Nxopt)
      print("Nyopt:",Nyopt)
    #endif not isilent

    if iretval:
      return Nmin,Nmax,Nmean,Nrms,Nxopt,Nyopt
    else:
      return
    #endif iretval

  else:

    scom = "global N; N = pd.DataFrame([" + svarx + "]).T"
    #print(scom)
    exec(scom)

    Nlines = len(N)
    scom = "global Nmin; Nmin = (" + svarx + ").min()"; exec(scom)
    scom = "global Nmax; Nmax = (" + svarx + ").max()"; exec(scom)
    scom = "global Nmean; Nmean = (" + svarx + ").mean()"; exec(scom)
    scom = "global Nrms; Nrms = (" + svarx + ").std()"; exec(scom)

    if not isilent:
      print("\nNmin:",Nmin)
      print("Nmax:",Nmax)
      print("Nmean:",Nmean)
      print("Nrms:",Nrms)
    #endif not isilent

    if iretval:
      return Nmin,Nmax,Nmean,Nrms
    else:
      return
    #endif iretval

  #endif len(var) > 1

#enddef nstat(nt='?',var='',select='',iretval,isilent)

def nmax(nt='?',var='',select='',iretval=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: nmax(Ntup,var='',select='', iretval=1")
    return

  if type(var) == str:
    if not len(var):
      print("*** Error in nmax: No columns specified in var ***")
      return
    var = nlistcolon(var)

  if Ncolon > 0:
    print("*** Error in nmax(...): Variable must contain exactly one term, i.e. no colons ***")
    return -1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N

  nt=N
  svar = nparse(nt,var[0])

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  Nlines = len(N)

  scom = "global Nmax; Nmax = " + svar + ".max()"
  #print(scom)
  exec(scom)

  if iretval:
    return Nmax
  else:
    return
#enddef nmax(nt='?',var='',select='',iretval=1)

def nmin(nt='?',var='',select='', iretval=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: nmin(Ntup,var='',select='', iretval=1")
    return

  if type(var) == str:
    if not len(var):
      print("*** Error in nmin: No columns specified in var ***")
      return
    var = nlistcolon(var)

  if Ncolon > 0:
    print("*** Error in nmin(...): Variable must contain exactly one term, i.e. no colons ***")
    return -1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N

  nt=N
  svar = nparse(nt,var[0])

  Nlines = len(N)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  scom = "global Nmin; Nmin = " + svar + ".min()"
  #print(scom)
  exec(scom)

  if iretval:
    return Nmin
  else:
    return
#enddef nmin(nt='?',var='',select='', iretval=1)

def nminmax(nt='?',var='',select='',iretval=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: nminmax(Ntup,var='',select='', iretval=1")
    return

  if type(var) == str:
    if not len(var):
      print("*** Error in nminmax: No columns specified in var ***")
      return
    var = nlistcolon(var)

  if Ncolon > 0:
    print("*** Error in nminmax(...): Variable must contain exactly one term, i.e. no colons ***")
    return -1
  #endif Ncolon > 0

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  svar = nparse(nt,var[0])

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  Nlines = len(N)

  scom = "global Nmin; Nmin = " + svar + ".min()"
  exec(scom)
  scom = "global Nmax; Nmax = " + svar + ".max()"
  exec(scom)

  if iretval:
    return Nmin, Nmax
  else:
    return
#enddef nmax(nt='?',var='',select='',iretval=1)

def rlist(filei="wlist.dat"):

  f = open(filei,"r")
  lis = []
  while True:
    line = f.readline()
    if not line: break
    words = line.split()
    lis.append(float(words[0]))
  f.close()

  return lis

def wlist(lis='?',fileo="wlist.dat"):

  if type(lis) == str: print("\nUsage: wlist(lis,file")

  f = open(fileo,"w")
  for i in range(len(lis)):
    f.write(str(lis[i])+"\n")
  #endfor i in range(len(lis)):
  f.close()

def lilo():
  global LogX, LogY
  LogX = 0
  LogY = 1

def lolo():
  global LogX, LogY
  LogX = 1
  LogY = 1

def loli():
  global LogX, LogY
  LogX = 1
  LogY = 0

def lili():
  global LogX, LogY
  LogX = 0
  LogY = 0

def Line(x1,y1,x2,y2,color=-9, width=2.):
  global Linewidth,Linecolor
  linecolor = Linecolor
  linewidth = Linewidth
  if color != -9: setlinecolor(color)
  setlinewidth(width)
  xpl = [x1,x2]
  ypl = [y1,y2]
  vplxy(xpl,ypl,'linesame')
  setlinecolor(linecolor)
  setlinewidth(linewidth)
#enddef Line

def LineNDC(x1,y1,x2,y2,color=-9, width=2.):
  global Linewidth,Linecolor

  linecolor = Linecolor
  linewidth = Linewidth

  if color != -9: setlinecolor(color)

  setlinewidth(width)

  ax = getax()

  if ax.get_xscale() == 'log': LogX = 1
  else: LogX = 0
  if ax.get_yscale() == 'log': LogY = 1
  else: LogY = 0

  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()

  if LogX:
    xlmin = log(xmin)
    xlmax = log(xmax)
    x1 = xlmin + x1 * (xlmax - xlmin)
    x2 = xlmin + x2 * (xlmax - xlmin)
    xpl = [exp(x1),exp(x2)]
  else:
    xpl = [xmin+x1*(xmax-xmin),xmin+x2*(xmax-xmin)]
  #endif LogX

  if LogY:
    ylmin = log(ymin)
    ylmax = log(ymax)
    y1 = ylmin + y1 * (ylmax - ylmin)
    y2 = ylmin + y2 * (ylmax - ylmin)
    ypl = [exp(y1),exp(y2)]
  else:
    ypl = [ymin+y1*(ymax-ymin),ymin+y2*(ymax-ymin)]
  #endif LogY

  vplxy(xpl,ypl,'linesame')

  setlinecolor(linecolor)
  setlinewidth(linewidth)

#enddef LineNDC

def Box(x1,x2,y1,y2,color=-9, width=2.):
  Line(x1,y1,x2,y1,color=-9, width=2.)
  Line(x1,y2,x2,y2,color=-9, width=2.)
  Line(x1,y1,x1,y2,color=-9, width=2.)
  Line(x1,y1,x1,y2,color=-9, width=2.)
#endif

def BoxNDC(x1,x2,y1,y2,color=-9, width=2.):
  LineNDC(x1,y1,x2,y1,color=-9, width=2.)
  LineNDC(x1,y2,x2,y2,color=-9, width=2.)
  LineNDC(x1,y1,x1,y2,color=-9, width=2.)
  LineNDC(x1,y1,x1,y2,color=-9, width=2.)
#endif

def textWC(x,y,text,fontsize='!', ishow=1, color='black',halign='center',
            valign='center',angle=0.):

  global TextFontSize, Textcolor

  if fontsize == '!': fontsize=TextFontSize

  col = color
  if color == -9: col = Textcolor

  Ax.text(x,y,text,fontsize=fontsize, color=col,
          horizontalalignment=halign,verticalalignment=valign,rotation=angle)

  if ishow: showplot()

#enddef textWC

def NDC_to_WC(x,y):

  # My version referes to plot not to the device

  ax = getax()

  if ax.get_xscale() == 'log': LogX = 1
  else: LogX = 0
  if ax.get_yscale() == 'log': LogY = 1
  else: LogY = 0

  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()

  if LogX:
    xlmin = log(xmin)
    xlmax = log(xmax)
    x = exp(xlmin + x * (xlmax - xlmin))
  else:
    x = xmin + x * (xmax - xmin)
  #endif LogX

  if LogY:
    ylmin = log(ymin)
    ylmax = log(ymax)
    y = exp(ylmin + y * (ylmax - ylmin))
  else:
    y = ymin + y * (ymax - ymin)
  #endif LogY

  return x,y
#enddef NDC_to_WC(xNDC)

def textNDC(x,y,text,fontsize='!', ishow=1, color='black',halign='center',
            valign='center',angle=0.):

  # My version referes to plot not to the device

  global TextFontSize, Textcolor

  # Compare to textndc()

  if fontsize == '!': fontsize=TextFontSize

  ax = getax()

  if ax.get_xscale() == 'log': LogX = 1
  else: LogX = 0
  if ax.get_yscale() == 'log': LogY = 1
  else: LogY = 0

  col = color
  if color == -9: col = Textcolor

  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()

  if LogX:
    xlmin = log(xmin)
    xlmax = log(xmax)
    x = exp(xlmin + x * (xlmax - xlmin))
  else:
    x = xmin + x * (xmax - xmin)
  #endif LogX

  if LogY:
    ylmin = log(ymin)
    ylmax = log(ymax)
    y = exp(ylmin + y * (ylmax - ylmin))
  else:
    y = ymin + y * (ymax - ymin)
  #endif LogY

  Ax.text(x,y,text,fontsize=fontsize, color=col,transform=Ax.transAxes,
          horizontalalignment=halign,verticalalignment=valign,rotation=angle)

  if ishow: showplot()

#enddef textNDC

def nrenvars(nt,varlis):
# Renames variable names of Ntuple to varlis
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  varlis = nlistcolon(varlis)

  idx= GetIndexN(nt)

  if idx == -1:
    print('*** Error in nrenvars(...): Index of Ntuple not found ***')
    return -1
  #endif type(nt) == Tdf or type(nt) == int:

  N.columns = varlis

  nhead = Nhead[idx]
  nvar = nhead[3]

  if nvar != len(varlis):
    print('*** Error in nrenvars(...): Wrong number of variables ***')
    return -1
  #endif type(nt) == Tdf or type(nt) == int:

  for k in range(nvar):
    l = 4 + k
    nhead[l][0] = varlis[k]
  #endfor k in range(nvar):

  return
#enddef nrenvars():

def nparse(nt,varlis):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global Gdebug
  #Gdebug = 1

  svar=''

  if type(nt) == str or type(nt) == int:

    i= GetIndexN(nt)

    if i == -1:
      print('*** Error in nparse(...): Index of Ntuple not found ***')
      return svar
    #endif

    nhead = Nhead[i]
    nnam = nhead[1]
    nvar = nhead[3]

    variables = []
    for k in range(nvar):
      l = 4 + k
      var = nhead[l][0]
      variables.append(var)
    #endfor k in range(nvar):

  elif type(nt) == Tdf:

    nnam = 'nt'
    variables = []
    for var in nt.columns: variables.append(var)

  else:
    print('*** Error in nparse(...): Unknown type of Ntuple ***')
    return svar

  #endif type(nt) == Tdf or type(nt) == int:

  variables.sort(key=len)
  variables.reverse()

  if Gdebug: print(variables)

  a = ''; itake = 0; atoms = []; iatoms = []

  ic = 0
  c = ''
  c2 = ''

  isnum = 0

  for c in varlis:

    c2 = c
    iwasnum = isnum

    if Gdebug: print("1::c,a,itake,atoms:",c,a,itake,atoms)
    isnum = 0

    if c in ['0','1','2','3','4','5','6','7','8','9','.'] \
    or iwasnum and c.lower() == 'e' \
    or iwasnum and (c2.lower() + c) == 'e-' \
    or iwasnum and (c2.lower() + c) == 'e+':
      isnum = 1
    #if c in ['0','1','2','3','4','5','6','7','8','9','.']

    idelim = 0
    if c in ['(','+','-','*','/',' ',')']:
      idelim = 1
    #endif c in ['(','+','-','*','/',' ',')']:

    if a == '' and idelim == 0 and isnum == 0: itake=1
    if Gdebug: print("2::isnum,idelim,itake:",isnum,idelim,itake)

    if itake and idelim or itake and ic == len(varlis) - 1:
      itake = 0
      if c != '(':
        atoms.append(a)
        iatoms.append(ic)
      a = ''
    if itake:
      a += c
    if Gdebug: print("3::a, atoms, itake:",a,atoms,itake,'\n--------------------------------\n')
    ic += 1
  #for c in varlis

  if Gdebug: print("4::atoms:",atoms)
  ia = -1; svl = 0
  for i2 in iatoms:
    ia += 1
    i1 = i2 - len(atoms[ia])
    sv = varlis[svl:i1] + nnam + '.' + varlis[i1:i2]
    if Gdebug: print(varlis[i1:i2],sv)
    svar += sv
    svl = i2
  sv = varlis[svl:]
  svar += sv
  if Gdebug: print(svar)

  return svar
#enddef nparse():

def make_empty_dataframe(varlis):
  global Tdf
  varlis = nlistcolon(varlis)
  df = pd.DataFrame(columns=varlis)
  Tdf = type(df)
  return df
#def make_empty_dataframe(varlis)

def make_dataframe(varlis='x',x='',y='',z='',s='', t='',bx='',by='',bz=''):

  varlis = nlistcolon(varlis)
  nvar = len(varlis)

  nt = make_empty_dataframe(varlis)

  nt[varlis[0]] = x
  if nvar > 1: nt[varlis[1]] = y
  if nvar > 2: nt[varlis[2]] = z
  if nvar > 3: nt[varlis[3]] = s
  if nvar > 4: nt[varlis[4]] = t
  if nvar > 5: nt[varlis[5]] = bx
  if nvar > 6: nt[varlis[6]] = by
  if nvar > 7: nt[varlis[7]] = bz

  return nt
#def make_dataframe(nt,varlis='x:y:z:bx:by:bz',x,y='',z='',bx='',by='',bz='')

def set_linecolor(lcol='r'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  Linecolor = lcol
#enddef set_linecolor(lcol='r')

def h2fill(idh='?', x=1.e30, y=1.e30, w=1.):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == str and idh == '?':
    print("h2fill(idh=-1, x=1.e30, y=1.e30, w=1.)")
    return
  #endif

  ind = GetIndexH2(idh)

  if ind == -1:
    print("*** Non existing histogram ***")
    return
  #endif

  head2 = H2head[ind]

  idh = head2[0]
  tit = head2[1]
  nx = head2[2]
  xmin = head2[3]
  xmax = head2[4]
  dx = head2[5]
  ny = head2[6]
  ymin = head2[7]
  ymax = head2[8]
  dy = head2[9]
  nent = head2[10]
  zmin = head2[11]
  zmax = head2[12]
  sumz = head2[13]
  meanx = head2[14]
  rmsx = head2[15]
  meany = head2[16]
  rmsy = head2[17]
  nunderx = head2[18]
  noverx = head2[19]
  underx = head2[20]
  overx = head2[21]
  nundery = head2[22]
  novery = head2[23]
  undery = head2[24]
  overy = head2[25]
  nunderxuy = head2[26]
  noverxuy = head2[27]
  underxuy = head2[28]
  overxuy = head2[29]
  nunderxoy = head2[30]
  noverxoy = head2[31 ]
  underxoy = head2[32]
  overxoy = head2[33]

  sumx = head2[34]
  sumx2 = head2[35]
  sumy = head2[36]
  sumy2 = head2[37]

  nent += 1
  head2[10] = nent
  sumz += w
  head2[13] = sumz

  if x < xmin:
    nunderx+=1
    underx+=w
    head2[18] = nunderx
    head2[20] = underx
    if y < ymin:
      nunderxuy+=1
      underxuy+=w
      head2[26] = nunderxuy
      head2[28] = underxuy
    elif y >= ymax:
      nunderxoy+=1
      underxoy+=w
      head2[30] = nunderxoy
      head2[32] = underxoy
    #endif
  elif x >= xmax:
    noverx+=1
    overx+=w
    head2[19] = noverx
    head2[21] = overx
    if y < ymin:
      noverxuy+=1
      overxuy+=w
      head2[27] = noverxuy
      head2[29] = overxuy
    elif y >= ymax:
      noverxoy+=1
      overxoy+=w
      head2[31] = noverxoy
      head2[33] = overxoy
    #endif
  #endif

  if y < ymin:
    nundery+=1
    undery+=w
    head2[22] = nundery
    head2[24] = undery
  elif y >= ymax:
    novery+=1
    overy+=w
    head2[23] = novery
    head2[25] = overy

  else:

    sumx2+=x**2*w #different from h1fill
    meanx = (meanx * sumx + x * w)/(sumx+w)
    sumx+=x*w #different from h1fill
    rmsx = (sumx2/sumx - meanx**2)**(0.5)

    head2[34] = sumx
    head2[35] = sumx2
    head2[14] = meanx
    head2[15] = rmsx

    sumy2+=y**2*w
    meany = (meany * sumy + y * w)/(sumy+w)
    sumy+=w
    rmsy = (sumy2/sumy - meany**2)**(0.5)

    head2[36] = sumy
    head2[37] = sumy2
    head2[16] = meany
    head2[17] = rmsy

    ichax = int((x-xmin)/dx) # numbering starts with zero
    ix = ichax + 1
    ichay = int((y-ymin)/dy) # numbering starts with zero
    iy = ichay + 1

    icha = ichay*nx + ichax

    h2 = H2[ind]

    h2['ix'][icha] = ix
    h2['iy'][icha] = iy
    h2['z'][icha] += w
    h2['z2'][icha] += w*w
    h2['n'][icha] += 1.
    h2['ave'][icha] = h2['z'][icha]/h2['n'][icha]
    h2['ez'][icha] = (h2['z2'][icha]*h2['n'][icha]-h2['y'][icha]**2)**0.5

    head2[11] = min(h2.z)
    head2[12] = max(h2.z)

  #endif x < xmin...

  return

#enddef h2fill

def h1fill(idh=-1, x=1.e30, wei=1.):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == int and idh <0:
    print("h1fill(idh=-1, x=1.e30, wei=1.)")
    return
  #endif

  ind = GetIndexH1(idh)

  if ind == -1:
    print("*** Non existing histogram ***")
    return

  head1 = H1head[ind]

  nx = head1[2]
  xmin = head1[3]
  xmax = head1[4]
  dx = head1[5]

  ny = head1[6]

  contmn = head1[7]
  contmx = head1[8]

  nent = head1[9]

  sumy = head1[10]
  xsqsum = head1[11]
  xmean = head1[12]
  xrms = head1[13]

  nunder = head1[14]
  nover = head1[15]
  under = head1[16]
  over = head1[17]

  sumx = head1[18]
  sumx2 = head1[19]

  if x < xmin:
    nunder+=1
    under+=wei
    head1[14] = nunder
    head1[16] = under
  elif x >= xmax:
    nover+=1
    over+=wei
    head1[15] = nover
    head1[17] = over
  else:

    # xrms**2 = sum((xi*wi)**2)/sum(wi) - (sum(xi*wi)/sum(wi))**2
    # xrms**2 * sum(wi) = sum((xi*wi)**2)/sum(wi) - (sum(xi*wi)/sum(wi))**2
    nent += 1

    xsqsum += x**2 * wei
    if sumy + wei != 0.0:
      xmean = (xmean * sumy + x * wei)/(sumy+wei)
    #endif

    sumy+=wei
    if sumy != 0.0:
      xrms = (xsqsum/sumy - xmean**2)**(0.5)
    #endif

    sumx+=x
    sumx2+=x*x

    head1[9] = nent
    head1[10] = sumy
    head1[11] = xsqsum
    head1[12] = xmean
    head1[13] = xrms
    head1[18] = sumx
    head1[19] = sumx2

    icha = int((x-xmin)/dx) # numbering starts with zero

    h1 = H1[ind]
    h1['y'][icha] += wei
    h1['y2'][icha] += wei*wei
    h1['n'][icha] += 1.
    h1['ave'][icha] = h1['y'][icha]/h1['n'][icha]
    h1['ey'][icha] = (h1['y2'][icha]/h1['n'][icha]-(h1['ave'][icha])**2)**0.5

    head1[7] = min(h1.y)
    head1[8] = max(h1.y)

  #endif x < xmin

  return

#enddef h1fill

def hbook2(idh=-1, tit='Histogram2D',
           nx=10, xmin=0., xmax=1., ny=10, ymin=0., ymax=1.0, overwrite=False):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == int and idh <0:
    print("hbook2(idh=-1, tit='Histogram2D', nx=10, xmin=0., xmax=1., ny=10, ymin=0., ymax=1.0, overwrite=False)")
    return 0

  ind = GetIndexH2(idh)

  if ind != -1 and not overwrite and Khdeleted == 0:
    print("*** Already existing histogram " + tit + ", will not overwrite it due to overwrite flag ***")
    return -1
  #endif ind != -1 and not overwrite and Khdeleted == 0

  if xmin >= xmax:
    print("*** Error in hbook2(...): xmin >= xmax ***")
    return -2
  #endif xmin >= xmax

  if ymin >= ymax:
    print("*** Error in bhook2(...): ymin >= ymax ***")
    return -3
  #endif ymin >= ymax

  if nx <= 0:
    print("*** Error in hbook2: Negative or zero number of channels ***")
    return -4
  elif nx == 1:
    dx = (xmax-xmin)
    x = np.array([(xmin+xmax)/2.])
  else:
    dx = (xmax-xmin)/nx
    x = np.linspace(xmin,xmax,nx+1)
    x = x[:-1] + dx/2.
  #endif nx <= 0:

  ix = np.arange(1,nx+1)

  if ny <= 0:
    print("*** Error in hbook2: Negative or zero number of channels ***")
    return -4
  elif ny == 1:
    dy = (ymax-ymin)
    y = np.array([(ymin+ymax)/2.])
  else:
    dy = (ymax-ymin)/ny
    y = np.linspace(ymin,ymax,ny+1)
    y = y[:-1] + dy/2.
  #endif ny <= 0:

  iy = np.arange(1,ny+1)

  nent = 0;

  sumx = 0.; sumx2=0.; meanx=0.; rmsx = 0.
  nunderx=0; noverx=0; underx=0.0; overx=0.0
  sumy = 0.; sumy2=0.; meany=0.; rmsy = 0.
  nundery=0; novery=0; undery=0.0; overy=0.0

  nunderxuy=0; noverxuy=0; underxuy=0; overxuy=0;
  nunderxoy=0; noverxoy=0; underxoy=0; overxoy=0;

  sumz = 0.; sumz2=0.; meanz=0.; rmsz = 0.; zmin = 0.; zmax = 0.;
  kdelete = 0

  head2 = [idh,tit, # 0 1
           nx,xmin,xmax,dx, #2 3 4 5
           ny,ymin,ymax,dy, # 6 7 8 9
           nent, zmin, zmax, sumz, meanx, rmsx, meany, rmsy, # 10 11 12 13 14 15 16 17
           nunderx, noverx, underx, overx, # 18 19 20 21
           nundery, novery, undery, overy, # 22 23 24 25
           nunderxuy, noverxuy, underxuy, overxuy, # 26 27 28 29
           nunderxoy, noverxoy, underxoy, overxoy, # 30 31 32 33
           sumx, sumx2, sumy, sumy2, kdelete #34 35 36 37 38
           ]

  xmesh, ymesh = np.meshgrid(x,y)
  ixm, iym = np.meshgrid(ix,iy)

  ixf = ixm.flatten()
  iyf = iym.flatten()
  xf = xmesh.flatten()
  yf = ymesh.flatten()
  zf = 0.0 * xf

  z2 = xf * 0.0
  ent = xf * 0.0
  zave = xf * 0.0
  ez = xf * 0.0

  cont2 = pd.DataFrame([ixf,iyf,xf,yf,zf,zave,ez,z2,ent]).T
  cont2.columns=['ix','iy','x','y','z','ave','ez','z2','n']

  if ind < 0:
    H2head.append(head2)
    H2.append(cont2)
    Nh2+=1
  else:
    H2head[ind] = head2
    H2[ind] = cont2
  #endif ind < 0

  return cont2
#enddef hbook2

def h2reset(idh):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  ind = GetIndexH2(idh)

  if ind < 0:
    print("*** h2reset: Non existing histogram ***")
    return
  #endif

  head2 = H2head[ind]

  H2h.z = 0.0
  H2h.ave = 0.0
  H2h.ez = 0.0
  H2h.z2 = 0.0
  H2h.n = 0.0

  idh = head2[0]
  tit = head2[1]
  nx = head2[2]
  xmin = head2[3]
  xmax = head2[4]
  dx = head2[5]
  ny = head2[6]
  ymin = head2[7]
  ymax = head[8]
  dy = head2[9]

  nent = 0;
  sumx = 0.; sumx2=0.; meanx=0.; rmsx = 0.
  nunderx=0; noverx=0; underx=0.0; overx=0.0
  sumy = 0.; sumy2=0.; meany=0.; rmsy = 0.
  nundery=0; novery=0; undery=0.0; overy=0.0

  nunderxuy=0; noverxuy=0; underxuy=0; overxuy=0;
  nunderxoy=0; noverxoy=0; underxoy=0; overxoy=0;

  sumz = 0.; sumz2=0.; meanz=0.; rmsz = 0.; zmin = 0.; zmax = 0.;
  kdelete = 0

  head2 = [idh,tit, # 0 1
           nx,xmin,xmax,dx, #2 3 4 5
           ny,ymin,ymax,dy, # 6 7 8 9
           nent, zmin, zmax, sumz, meanx, rmsx, meany, rmsy, # 10 11 12 13 14 15 16 17
           nunderx, noverx, underx, overx, # 18 19 20 21
           nundery, novery, undery, overy, # 22 23 24 25
           nunderxuy, noverxuy, underxuy, overxuy, # 26 27 28 29
           nunderxoy, noverxoy, underxoy, overxoy, # 30 31 32 33
           sumx, sumx2, sumy, sumy2, kdelete #34 35 36 37 38
           ]

#enddef h2reset

def get_htitle(h):
  idx = GetIndexH1(h)
  if idx == -1:
    idx = GetIndexH2(h)
    if idx == -1:
      print("*** Non-existing histogram ",h, " ***")
      return None
    return H2hh[1]
  return H1hh[1]

def set_htitle(h,title):
  idx = GetIndexH1(h)
  if idx == -1:
    idx = GetIndexH2(h)
    if idx == -1:
      print("*** Non-existing histogram ",h, " ***")
      return
    else:
      H2hh[1] = title
  else:
    H1hh[1] = title

def hbook1(idh=-1, tit='Histogram1D', nx=10, xmin=0., xmax=1., overwrite=False):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == int and idh <0:
    print("hbook1(idh=-1, tit='Histogram1D', nx=10, xmin=0., xmax=1.)")
    return 0
  #endif type(idh) == int and idh <0

  ind = GetIndexH1(idh)

  if ind != -1 and not overwrite and Khdeleted == 0:
    print("*** Already existing histogram " + ", will not overwrite it due to overwrite flag ***")
    return -1
  #endif ind != -1 and not overwrite and Khdeleted == 0

  if xmin >= xmax:
    if xmin == 0:
      xmin = -1
      xmax = 1
    else:
      dx = xmin / 10.
      xmin -= dx
      xmax += dx
    #endif
    #print("*** Error in hbook1(...): xmin >= xmax ***")
    #return -2
  #endif xmin >= xmax

  if nx <= 0:
    print("*** Error in hbook1: Negative or zero number of channels ***")
    return -3
  elif nx == 1:
    dx = (xmax-xmin)
    x = np.array([(xmin+xmax)/2.])
  else:
    dx = (xmax-xmin)/nx
    x = np.linspace(xmin,xmax,nx+1)
    x = x[:-1] + dx/2.
  #endif nx <= 0:

  nent = 0; nunder = 0; nover = 0
  sum = 0.; sumsq=0.; rms = 0.; contmn = 0.; contmx = 0.; xmean = 0.
  over = 0.; under = 0.; sumx = 0.; sumx2 = 0.
  kdeleted = 0

  ny=0
  head1 = [idh,tit,nx,xmin,xmax,dx,ny,contmn,contmx,
           nent,sum,sumsq,xmean,rms,nunder,nover,under,over,sumx,sumx2,kdeleted]

  ix = np.arange(1,nx+1)

  y = x * 0.0
  y2 = x * 0.0
  ent = x * 0.0
  ave = x * 0.0
  ey = x * 0.0

  cont1 = pd.DataFrame([ix,x,y,ave,ey,y2,ent]).T
  cont1.columns=['ix','x','y','ave','ey','y2','n']

  if ind < 0:
    H1head.append(head1)
    H1.append(cont1)
    Nh1+=1
  else:
    H1head[ind] = head1
    H1[ind] = cont1
  #endif idn < 0

  return cont1

#enddef hbook1

def nscan(nt='?',varlis='',select='',isilent=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: nscan(nt,varlis,select)")
    return
  #if type(nt) == str and nt == '?'

  if type(nt) == int and nt < 0:
    pass
  else:
    ind = GetIndexN(nt)
    nupdate_header(Ntup[ind])
    if ind == -1: return -1
  #endif type(nt) == int and nt < 0

  if type(varlis) == str:
    if varlis == '':
      varlis = list(N.columns)
      Ncolon = len(varlis) - 1
    else:
      varlis = nlistcolon(varlis)
    #endif varlis == ''
  #endif type(varlis) == str

  nhead = Nhead[ind]

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  svar = ""

  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
  #endfor k in range(Ncolon)

  sv = nparse(nt,varlis[Ncolon])
  svar += nparse(nt,varlis[Ncolon])
  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  cols = N.columns
  nc = 0

  svar = []

  for k in range(len(cols)):
    try:
      svar.append(varlis[varlis.index(cols[k])])
    except:
      nc += 1
      svar.append('C' + str(nc))
    #endtry
  #endfor k in range(len(cols))

  N.columns = svar
  if not isilent: print(N.to_string())

  Nlines = len(N)
  return Nlines

#enddef nscan

def nfitxy(nt='?',varlis='',select='',fitfun=None, absolute_sigma='default',
           parstart=None, bounds=None, method=None,isilent=0, ninter=101,
           iretval=1,iplot=0,fitcol='!',kweedzero=1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  # Histograms and Ntuples
  global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
  Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
  H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
  Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
  Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
  H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
  FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman

  #global N,Ney,Nfit #,FitFit, FitPar #,Nfitint ,Nfitxy

  if type(nt) == str and nt == '?':
    print("\nUsage: nfitxy(nt='',varlis='',select='',fitfun=None, absolute_sigma='default'," + \
    "parstart=None, bounds=None, method=None,isilent=0, ninter=101," + \
    "iretval=0)")
    print("varlis = x:y or x:y:ey")
    print("If iretval: return par,sigma,chi2ndf,f")
    return
  #endif

  if fitfun == None:
    print("*** Error in nfitxy(...): NO fitfun provided ***")
    return -1,None,None,None
  #endif not fitfun

  if type(varlis) == str: varlis = nlistcolon(varlis)

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1,None,None,None
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2,None,None,None
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select:

  nt=N
  var = []

  svar = ""
  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor

  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  if Ncolon == 1:
    svar += "," + sv + "*0.0"
    scom = "global N; N = pd.DataFrame([" + svar + "]).T"
    exec(scom)
  elif Ncolon == 2:
    scom = "global N; N = pd.DataFrame([" + svar + "]).T"
    exec(scom)
  else:
    print("*** Error in nfitxy(...): Two or three items must be given in varlis ***")
    return -1,None,None,None
  #endif Ncolon == 1

  N.columns = ['x','y','ey']
  Nlines = len(N)
  N.index = range(Nlines)

  if Ncolon == 1:
    ey = ''
  else:
    ey = N.ey
  #endif

  igauss = 0

  if type(fitfun) == str:
    if fitfun == 'expo':
      fitfun = expo
      par,sigma,chi2ndf,f = vfitexp(N.x,N.y,ey,\
      absolute_sigma, parstart, \
      bounds, method,isilent, ninter,kweedzero)
    elif fitfun == 'expo2':
      fitfun = expo2
      par,sigma,chi2ndf,f = vfitexp2(N.x,N.y,ey,\
      absolute_sigma, parstart, \
      bounds, method,isilent, ninter,kweedzero)
    elif fitfun == 'gauss':
      igauss = 1
      fitfun = gauss
      if parstart == None:
        gmin,gmax,gmean,grms, gxopt, gyopt = nstat(N,var='x:y',iretval=1,isilent=1)
        parstart = [gyopt,gmean,grms]
      #endif
      par,sigma,chi2ndf,f = vfitgauss(N.x,N.y,ey,\
      absolute_sigma, parstart, \
      bounds, method,isilent, ninter,kweedzero)
    elif fitfun == 'cosh':
      fitfun = fcosh
      par,sigma,chi2ndf,f = vfitcosh(N.x,N.y,ey,\
      absolute_sigma, parstart, \
      bounds, method,isilent, ninter,kweedzero)
    elif fitfun == 'cos':
      fitfun = fcos
      par,sigma,chi2ndf,f = vfitcos(N.x,N.y,ey,\
      absolute_sigma, parstart, \
      bounds, method,isilent, ninter,kweedzero)
    #endif
  elif type(fitfun) == int:
    nord = fitfun
    par,sigma,chi2ndf,f = vfitpoly(nord,N.x,N.y,ey, cov='default',
                                   isilent=isilent, ninter=ninter,
                                   iretval=iretval,kweedzero=kweedzero)
  else:
    par,sigma,chi2ndf,f = vfit(fitfun,N.x,N.y,ey,\
    absolute_sigma, parstart, \
    bounds, method,isilent, ninter,kweedzero)
  #endif type(fitfun) == str

  x = N.x
  y = N.y
  ey = N.ey

  Nfitxy = ncre("Nfitxy","x:y:ey:f",ioverwrite=1)

  Nfitxy.x = x
  Nfitxy.y = y
  Nfitxy.ey = ey
  Nfitxy.f = f

  nw = Nfitxy.drop_duplicates()

  Nfitxy.x = nw.x
  Nfitxy.y = nw.y
  Nfitxy.ey = nw.ey
  Nfitxy.f = nw.f

  nupdate_header(Nfitxy)

  if iplot:

    if fitcol == '!': fitcol = 'black'

#    Kold = Kstat
#    if Kstat: Kstat = False

    npl(Nfitxy,"x:y")

    if Nfitxy.ey.abs().max() > 0:
      lcol = getlinecolor()
      plt.errorbar(Nfitxy.x,Nfitxy.y,Nfitxy.ey, \
      ls='',marker=Markertype, fillstyle=Fillstyle, mfc=lcol, \
      mec=lcol, ms=Markersize, mew=1, c=lcol)
    #endif

    x = Nfitxy.x.drop_duplicates()
    if len(x) == len(Nfitxy.x):
      vplxy(x,Nfitxy.f[:len(x)],"samespline",color=fitcol)
    else:
      vplxy(x,Nfitxy.f[:len(x)],"sameline",color=fitcol)
    #endif

#    Kstat = Kold

    if Kfit:
      if StatFontSize < 0:
        dpi = Fig.dpi
        nxy = max(Nxzone,Nyzone)
        sfs = 12 - nxy * 2
      else:
        sfs = StatFontSize
      #endif StatFontSize < 0

      sig = FitSig

      if igauss:
        tex = "A = " + '{:.4g}'.format(par[0]) + " +/- " + '{:.4g}'.format(sig[0]) + \
        "\n$\mu$ = " + '{:.4g}'.format(par[1])  + " +/- " + '{:.4g}'.format(sig[1]) + \
        "\n$\sigma$ = " + '{:.4g}'.format(par[2])  + " +/- " + '{:.4g}'.format(sig[2]) + "\n"
      else:
        tex = ""
        ip = 0
        for p in par:
          tex += "P" + str(ip) + " = " + '{:.4g}'.format(p)  + '{:.4g}'.format(sig[ip]) + "\n"
          ip += 1
        #endfor
      #endif

      tex += "$\chi^2/Ndf$" + " = " + '{:.4g}'.format(chi2ndf) + "\n"
      tex += "$\chi^2 prob$" + " = " + '{:.4g}'.format(FitChi2Prob) + "\n"
      text(Xstat,Ystat,tex,halign='left')

    #endif Kstat

  #endif iplot

  if iretval: return par,sigma,chi2ndf,f
  else: return
#enddef nfitxy(nt='?',varlis='',select='',fitfun=None, absolute_sigma='default',

def nfitg(nt='?',varlis='',select='',fitfun=None, absolute_sigma='default',
           parstart=None, bounds=None, method=None,isilent=0, ninter=101,
           iretval=1,iplot=0,fitcol='!'):

  if iretval:
    return nfitxy(nt,varlis,select,'gauss',absolute_sigma,parstart,bounds,method,
         isilent, ninter,iretval,iplot,fitcol)
  else:
    nfitxy(nt,varlis,select,'gauss',absolute_sigma,parstart,bounds,method,
         isilent, ninter,iretval,iplot,fitcol)
  #endif
#enddef

def nfitp1(nt='?',varlis='',select='',fitfun=None, absolute_sigma='default',
           parstart=None, bounds=None, method=None,isilent=0, ninter=101,
           iretval=1,iplot=0,fitcol='!'):
  if iretval:
    par,sigma,chi2ndf,f = nfitxy(nt,varlis,select,1, absolute_sigma, \
           parstart, bounds, method,isilent, ninter, \
           iretval,iplot,fitcol)
    return par,sigma,chi2ndf,f
  else:
    nfitxy(nt,varlis,select,1, absolute_sigma, \
           parstart, bounds, method,isilent, ninter, \
           iretval,iplot,fitcol)
    return
  #endif


def nintern(nt='?',varlis='',select='',xint='!'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and ( nt == '?' or nt == '' ):
    print("nintern(nt='?',varlis='',select='',xint='!')")
    return

  if type(varlis) == str:
    varlis = nlistcolon(varlis)

  if Ncolon != 1:
    print("*** Error in nintern(...): Two items must be given in varlis ***")
    return -1
  #endif Ncolon != 1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  var = []
  svar = ""

  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor k in range(Ncolon)

  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  N.columns = ['x','y']
  Nlines = len(N)
  N.index = range(Nlines)

  yint = vintern(N.x,N.y,xint)

  return

#enddef nintern(nt='',varlis='',select='',xint='!')

def ninter(nt='?',varlis='',select='',xint='!'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and ( nt == '?' or nt == '' ):
    print("ninter(nt='?',varlis='',select='',xint='!')")
    return
  #endif type(nt) == str and ( nt == '?' or nt == '' )

  if type(varlis) == str:
    varlis = nlistcolon(varlis)

  if Ncolon != 1:
    print("*** Error in ninter(...): Two items must be given in varlis ***")
    return -1
  #endif Ncolon != 1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  var = []
  svar = ""

  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor k in range(Ncolon)

  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  Ninter = ncre("Ninter","x:y:yp:ypp:yint",ioverwrite=1)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  N.columns = ['x','y']
  Nlines = len(N)
  N.index = range(Nlines)

  if type(xint) == list:
    xint = np.array(xint)
  elif type(xint)==str and xint=='!':
    xmin = N.x.min()
    xmax = N.x.max()
    xint = vec(xmin,xmax,(max(2,Nlines)-1)*10)
  elif type(xint) == float or type(xint) == int:
    xint = np.array([xint])
  #endif type(xint) == list

  yint = vinter(N.x,N.y,xint)

  Ninter.x = xint
  Ninter.y = yint

  yp = deepcopy(xint)
  ypp = deepcopy(xint)
  yinteg = deepcopy(xint)

  n = len(xint)-1

  yp[0] = (yint[1]-yint[0]) / (xint[1]-xint[0])
  yp[n] = (yint[n]-yint[n-1]) / (xint[n]-xint[n-1])
  yinteg[0] = 0.0

  for i in range(1,n):
    yinteg[i] = yinteg[i-1] + (yint[i]+yint[i-1])/2. * (xint[i]-xint[i-1])
    yp[i] = (yint[i+1]-yint[i-1]) / (xint[i+1]-xint[i-1])
  #endfor

  yinteg[n] = yinteg[n-1] + (yint[n]+yint[n-1])/2. * (xint[n]-xint[n-1])

  ypp[0] = (yp[1]-yp[0]) / (xint[1]-xint[0])
  ypp[n] = (yp[n]-yp[n-1]) / (xint[n]-xint[n-1])

  for i in range(1,n):
    ypp[i] = (yp[i+1]-yp[i-1]) / (xint[i+1]-xint[i-1])
  #endfor

  Ninter.yp = yp
  Ninter.ypp = ypp
  Ninter.yint = yinteg

  nupdate_header(Ninter)

  return

#enddef ninter(nt='',varlis='',select='',xint='!')

def nspline(nt='?',varlis='',select='',xspl='!',periodic=False):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and ( nt == '?' or nt == '' ):
    print("\nnspline(nt='?',varlis='',select='',xspl='!',periodic=False,mode='!')")
    print("\nIf x[len(x)-1] = x[0] vspline_index is used instead of vspline.")
    print("In this case xspl should be an integer, giving the number of spline points returned,")
    print("Otherwize the default of 1001 points is used. Vspline_index is called")
    print("with periodic=True.")
    return
  #endif type(nt) == str and ( nt == '?' or nt == '' )

  if type(varlis) == str: varlis = nlistcolon(varlis)

  if Ncolon != 1:
    print("*** Error in nspline(...): Two items must be given in varlis ***")
    return -1
  #endif Ncolon != 1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  var = []
  svar = ""

  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor k in range(Ncolon)

  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  N.columns = ['x','y']
  Nlines = len(N)
  N.index = range(Nlines)

  dx = abs(N.x[Nlines-1]-N.x[0]) / N.x.std()

  if dx < 1.0e-9:
    nspl = xspl
    if type(xspl) != int:
      print("*** Warning in nspline(...): First and last x are the same,")
      print("so, vspline_index is called, but xspl is not integer, Will use")
      print("default value 1001.")
      nspl = 1001
    yspl = vspline_index(N.x,N.y,nspl, periodic=True)
  else:
    yspl = vspline(N.x,N.y,xspl, periodic=False)
  #endif

  return

#enddef nspline(nt='',varlis='',select='',xspl='!',periodic=False)

def nsolve(nt='?',varlis='',select='',val=0.0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and ( nt == '?' or nt == '' ):
    print("nsolve(nt='?',varlis='',select='',val=0.0")
    return
  #endif type(nt) == str and ( nt == '?' or nt == '' )

  if type(varlis) == str: varlis = nlistcolon(varlis)

  if Ncolon != 1:
    print("*** Error in nsolve(...): Two items must be given in varlis ***")
    return -1
  #endif Ncolon != 1

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  var = []
  svar = ""
  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor k in range(Ncolon)
  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  N.columns = ['x','y']
  Nlines = len(N)
  N.index = range(Nlines)

  if Nlines < 100:
    xspl = N.x
    xspl = vec(N.x[0],N.x[Nlines-1],100)
    yspl = vspline(N.x,N.y,xspl)
    xsolve = vsolvelin(Nspline.x,Nspline.y,val)
  else:
    xsolve = vsolve(N.x,N.y,val)
  #endif Nlines < 100:

  return float(xsolve)

#enddef nsolve(nt='',varlis='',select='',val=0.0)

def ndump(nt='',varlis='',select='',fout='ndump.dat', sep=' ',floatform='%.5e',
          k_header=True, k_index=False,decimal="."):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  varliso = varlis

  if not len(nt):
    print("nlines=ndump(nt='',varlis='',select='',fout='ndump.dat',sep=' ',floatform='%.5e',k_header=True, k_index=False)")
    print("ihead=1: print variables as header")
    print("else: no header")
    print("nt < 0: Ntuple to dump is N")
    print("Examples for float format: 6g")
    return 0
  #endif not len(nt)


  if type(nt) == Tdf:
    nupdate_header(nt)
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(varlis) == 0:
    varl = list(N.columns)
    varlis = varl[0]
    for i in range(len(varl)-1):
      varlis += ':' + varl[i+1]
    #endfor
    varliso = varlis
  #endif

  varlis = nlistcolon(varlis)

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  var = []
  svar = ""

  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor

  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  if k_header:
    Fout = open(fout,'w')
    Fout.write("* " + str(varliso) + NL)
    Fout.close()
    mod = 'a'
  else:
    mod = 'w'
  #endif

  if floatform == '!':
    N.to_csv(fout,header=False,index=False,sep=sep,decimal=decimal,mode=mod)
  else:
    N.to_csv(fout,header=False,index=False,sep=sep,mode=mod,
             float_format=floatform,decimal=decimal)
  #endif

  Nlines = len(N)
  return Nlines

#enddef ndump----------------------------------------------------------------

def nlistcolon(varlis):
  global Ncolon
# converts variable list like 'x:y' to ['x','y']
  if type(varlis) == str:
    vari = varlis.split(':')
    ncolon = len(vari) - 1
    Ncolon = ncolon
    varlis=[]
    for i in range(ncolon+1):
      varlis.append(vari[i])
    #endfor i in range(ncolon+1)
  elif type(varlis) == list:
    Ncolon = len(varlis) - 1
  return varlis
#enddef nlistcolon(varlis)

def nreset(nt='?', varlis=''):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: Ntup = nreset(Ntup, varlis")
    return
  #endif

  idx = GetIndexN(nt)

  if idx == -1:
      print("*** Error in nreset: Non-existing Ntuple ***")
      return -1
  #endif

  nt= N

  if len(varlis) == 0:
    nvar = Nhead[idx][3]
    varlis = []
    for k in range(nvar):
      varlis.append(Nhead[idx][4+k][0])
      Nhead[idx][4+k][1] = None
      Nhead[idx][4+k][2] = None
    Nhead[idx][5+k] = 0
  else:
    nvar = len(varlis)
    for k in range(nvar):
      Nhead[idx][4+k][1] = None
      Nhead[idx][4+k][2] = None
      print(k)
    #endfor
    Nhead[idx][5+k] = 0
  #endif len(varlis) == 0:

  contnt = []
  Ntup[idx] = pd.DataFrame(contnt,columns=varlis)
  nt = Ntup[idx]

  NLast = N
  NiLast = Ind
  N = nt

  Ind = Nntup
  Ni = Ind

  return nt

#enddef nreset(...)

def ndelete(nt='?',isilent=0):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str and nt == '?':
    print("\nUsage: ndelete(Ntup")
    return
  #endif type(nt) == str and nt == '?'

  idx = GetIndexN(nt,isilent)

  if idx == -1:
    if isilent == 0: print("*** Error in ndelete: Non-existing Ntuple ***")
    return
  #endif idx == -1:

  Nntup -= 1
  if Nntup == 0:
    Nhead = []
    Ntup = []
    Nind = []
    return
  #endif

  Nhead[idx] = Nhead[Nntup-1]
  Ntup[idx] = Ntup[Nntup-1]
  Nind[idx] = Nhead[idx][1]

  kdx = Nntup - 1

  Nhead.pop(kdx)
  Ntup.pop(kdx)
  Nind.pop(kdx)

  for i in range(Nntup): Nhead[i][0] = i

#enddef ndelete(...)

def ncre(ntname='', nttit='', varlis='', ioverwrite=0):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  #reakpoint()
  if len(ntname) == 0:
    print("ncre(Name, Title, List of Variables, e.g. 'x:y', ioverwrite=0)")
    return -1
  #endif len(ntname) == 0

  if nttit != '' and varlis == '':
    varlis = nttit
    nttit = ntname
  #if nttit != '' and varlis == ''

  if type(varlis) == str:
    varlis = nlistcolon(varlis)

  for var in varlis:
    if var == 'mode':
      print("*** Error in ncre(...) for Ntuple ", ntnam)
      Quit("*** Variables must not be named 'mode'")
    #endif
  #endfor

  idx = -1

  for i in range(Nntup):
    if Nind[i] == ntname:
      idx = i
      break
  #endfor i in range(Nntup):

  if idx != -1:
    if not ioverwrite:
      print("*** Error in ncre: Already existing Ntuple " + ntname + "  ***")
      return -1
    else:
      ndelete(ntname)
#      print("vor nreset")
#      nt = nreset(ntname,varlis)
#      print("nach nreset")
#      return nt
    #endif not ioverwrite
  #endif idx != -1

  headnt = [Nntup,ntname]
  if len(nttit) == 0: nttit = ntname

  nvar = len(varlis)

  headnt.append(nttit)
  headnt.append(nvar)

  for i in range(nvar):
    if (varlis[i] == "imag"):
      print("*** Error in ncre(...) for Ntuple ", ntname)
      Quit("*** Variable name: 'imag' is reserved ***")
    if (varlis[i] == "int"):
      print("*** Error in ncre(...) for Ntuple ", ntname)
      Quit("*** Variable name: 'int' is reserved ***")
    if (varlis[i] == "nint"):
      print("*** Error in ncre(...) for Ntuple ", ntname)
      Quit("*** Variable name: 'nint' is reserved ***")
    #endif
    var = [varlis[i],None,None]
    headnt.append(var)
  #endfor i in range(nvar)

  nent = 0
  headnt.append(nent)

  Nhead.append(headnt)

  contnt = []
  contnt = pd.DataFrame(contnt,columns=varlis)

  Ntup.append(contnt)
  Nind.append(ntname)

  nt = contnt
  NLast = N
  NiLast = Ind
  N = nt

  Ind = Nntup
  Ni = Ind

  Nntup+=1

  return nt

#enddef ncre(...)

def GetIndexH2(idh='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == str and idh == '?':
    print("\nUsage: index = GetIndexH2(idh), H2h, H2I and Ind are also set, \nreturns -1 if histogram has not been found")
    return
  #endif type(idh) == str:

  if Nh2 <= 0: return -1

  if type(idh) == int:
    idx = idh
  elif type(idh) == Tdf:
    idx = -1
    for i in range(len(H2)):
      if id(H2[i]) == id(idh):
        idx = i
        break
      #endif id(H2[i]) == id(idh):
    #endfor i in range(H2):
  else:
    idx = -1
    for i in range(Nh2):
      if H2head[i][0] == idh:
        idx = i
        break
    #endfor i in range(Nh1):

  #endif type(idh) == int:

  Ind = idx
  H2I = Ind
  H2h = H2[idx]
  H2hh = H2head[idx]
  H2HLast = H2head[idx]
  H2Last = H2h
  Khdeleted = H2hh[38]

  return idx
#def GetIndexH2(idh='?')

def GetIndexN(nt='?', isilent=0):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  reset_status()

  idx = -1

  if not len(Ntup): return idx

  if type(nt) == int:
    idx = nt
  else:

    if type(nt) != Tdf:

      if nt == '?':
        print("\nUsage: index = GetIndex(ntname or nt), N, Ni and Ind are also set, \nreturns -1 if histogram has not been found")
        return idx
      #endif nt == '?':

      for i in range(Nntup):
        if Nind[i] == nt:
          idx = i
          break
      #endfor i in range(Nntup):

    else:

      #if Gdebug == 1: quit()

      for i in range(Nntup):
        if id(Ntup[i]) == id(nt):
          idx = i
          break
        #endif id(Ntup[i]) == id(nt):
      #endfor i in range(Nntup):

    #endif not type(nt) == Tdf:

    if idx == -1:
      if type(nt) != Tdf:
        set_error("*** GetIndexN: Index not found. Ntuple: " + str(nt) + " ***")
      else:
        set_error("*** GetIndexN: Index not found ***")
      #endif
      if not isilent: print(ErrorText)
      return -1
    #endif idx == -1:

  #endif type(nt) == int:

  Ind = idx
  Ni = Ind
  NLast = N
  N = Ntup[Ind]

  return idx

#enddef GetIndexN(nt='?', isilent=0)

def GetIndexNct(idh='?'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == int:
    idx = idh
  else:
    if idh == '?':
      print("\nUsage: index = GetIndex(idh), NctI and Ind are also set, \nreturns -1 if histogram has not been found")
      return
    else:
      idx = -1
      for i in range(Nnctup):
        if Nctup[i][0][0] == idh:
          idx = i
          break
        #endif Nctup[i][0][0] == idh:
      #endfor i in range(Nh1):
    #endif idh == '?':
  #endif type(idh) == int:

  Ind = idx
  NctI = Ind

  return idx

def GetIndexH1(idh='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  #nreakpoint()
  if type(idh) == str:
    if idh == '?':
      print("\nUsage: index = GetIndex(idh), returns -1 if histogram has not been found")
      return
  #endif type(idh) == str:

  Ind = -1
  H1HLast = H1h
  H1I = Ind
  H1h = None

  if Nh1 <= 0: return -1

  if type(idh) == int:
    idx = idh
  elif type(idh) == Tdf:
    idx = -1
    for i in range(len(H1)):
      if id(H1[i]) == id(idh):
        idx = i
        break
      #endif id(H1[i]) == id(idh):
    #endfor i in range(H1):
  else:
    idx = -1
    for i in range(Nh1):
      if H1head[i][0] == idh:
        idx = i
        break
      #endfor i in range(Nh1):
  #if type(nt) == int

  Ind = idx
  H1HLast = H1h
  H1I = Ind
  H1h = H1[Ind]
  H1hh = H1head[Ind]
  Khdeleted = H1hh[20]

  return idx
#enddef GetIndexH1(idh='?')

def GetIndex(idh='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if idh == '?':
    print("\nUsage: index = GetIndex(idh), H1h, H1I and Ind are also set, \nreturns [-1,None] item has not been found")
    return
  #endif

  Ind = GetIndexH1(idh)

  if Ind != -1:
    H1I = Ind
    return [Ind,'h1']
  #endif

  Ind = GetIndexH2(idh)

  if Ind != -1:
    H2I = Ind
    return [Ind,'h2']
  #endif

  Ind = GetIndexN(idh)

  if Ind != -1:
    Ni = Ind
    return [Ind,'nt']
  #endif

  Ind = GetIndexNct(idh)

  if Ind != -1:
    NctI = Ind
    return [Ind,'nct']
  #endif

  return [-1,None]
#enddef GetIndex(idh='?')

def getn(idn):
  if type(idn) == int: return Ntup[idn]
  else: return(Ntup[GetIndexN(idn)])

def h1opt(idh):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  xopt = None
  yopt = None

  if type(idh) == Tdf:

    h1 = idh

  else:

    if type(idh) == int:
      idx = idh
    else:
      idx = GetIndexH1(idh)
    #if type(idh) == int

    H1HLast = H1head[idx]
    H1Last = H1[idx]
    h1h = H1HLast
    h1 = H1h
  #endif type(idh) == Tdf

  ymin = 1.e30; ymax = -1.e30

  #nx = h1h[2]
  nx = len(h1)

  if nx < 3: return xopt, yopt

  ixmin = h1.y.idxmin() + 1
  ixmax = h1.y.idxmax() + 1
  ymin = h1.y.min()
  yxmax = h1.y.max()
#  for ix in range(nx):
#    y = h1.y[ix]
#    if y < ymin: ixmin = ix+1; ymin = y
#    if y > ymax: ixmax = ix+1; ymax = y
#  #endfor ix in range(h1h[2])

  if abs(ymin) <= abs(ymax):
    if ixmax < nx and ixmax > 1:  ix = ixmax - 1; ixm = ix - 1; ixp = ix + 1
    elif ixmax == nx:  ixp = nx - 1; ix = ixp - 1; ixm = ix - 1
    elif ixmax == 1:  ixp = 2; ix = 1; ixm = 0
    x = [h1.x[ixm],h1.x[ix],h1.x[ixp]]
    y = [h1.y[ixm],h1.y[ix],h1.y[ixp]]
  else:
    if ixmin < nx and ixmin > 1: ix = ixmin - 1; ixm = ix - 1; ixp = ix + 1
    elif ixmin == nx: ixp = nx - 1; ix = ixp - 1; ixm = ix - 1
    elif ixmin == 1:  ixp = 2; ix = 1; ixm = 0
    x = [h1.x[ixm],h1.x[ix],h1.x[ixp]]
    y = [h1.y[ixm],h1.y[ix],h1.y[ixp]]
  #endif abs(ymin) <= abs(ymax)

  a, yp, opt, ifail = util_parabel(x,y)

  if ifail:
    xopt = None
    yopt = None
  else:
    xopt = opt[0]
    yopt = opt[1]
  #endif

  return xopt,yopt
#enddef h1opt(idh)

def voptpar(vx,vy):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  xopt = None
  yopt = None

  vx = pd.Series(vx)
  vy = pd.Series(vy)

  nx = len(vx)
  mx = nx - 1

  if len(vy) != nx:
    print("*** Error in voptpar(x,y): Length if vectors x and y differ ***")
    return
  #if len(vy) != nx

  if nx < 3: return xopt, yopt

  ixmin = vy.idxmin() + 1
  ixmax = vy.idxmax() + 1
  ymin = vy.min()
  ymax = vy.max()

  if abs(ymin) <= abs(ymax):
    if ixmax < nx and ixmax > 1:  ix = ixmax - 1; ixm = ix - 1; ixp = ix + 1
    elif ixmax == nx:  ixp = nx - 1; ix = ixp - 1; ixm = ix - 1
    elif ixmax == 1:  ixp = 2; ix = 1; ixm = 0
    x = [vx[ixm],vx[ix],vx[ixp]]
    y = [vy[ixm],vy[ix],vy[ixp]]
  else:
    if ixmin < nx and ixmin > 1: ix = ixmin - 1; ixm = ix - 1; ixp = ix + 1
    elif ixmin == nx: ixp = nx - 1; ix = ixp - 1; ixm = ix - 1
    elif ixmin == 1:  ixp = 2; ix = 1; ixm = 0
    x = [vx[ixm],vx[ix],vx[ixp]]
    y = [vy[ixm],vy[ix],vy[ixp]]
  #endif abs(ymin) <= abs(ymax)

  x = [vx[ixm],vx[ix],vx[ixp]]
  y = [vy[ixm],vy[ix],vy[ixp]]

  a, yp, opt, ifail = util_parabel(x,y)

  if ifail:
    xopt = None
    yopt = None
  else:
    xopt = opt[0]
    yopt = opt[1]
  #if ifail

  VoptX = xopt
  VoptY = yopt

  return xopt,yopt
#enddef voptpar(x,y)

def h1print(idh='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  if type(idh) == str and idh == '?':
    print("\nUsage: h1print(idh), H1Last, H1I, and Ind are also set, \nreturns -1 if histogram has not been found")
    return

  if type(idh) == int:
    idx = idh
  else:
    idx = GetIndexH1(idh)

  if idx == -1: return idx

  H1HLast = H1head[idx]
  H1Last = H1[idx]
  h1h = H1HLast

  print(H1[idx].to_string())

#def h1print(idh='?')

def H1Info(idh='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == str and idh == '?':
    print("\nUsage: H1Info(idh), H1Last, H1I, and Ind are also set, \nreturns -1 if histogram has not been found")
    return
  #endif

  if type(idh) == int:
    idx = idh
  else:
    idx = GetIndexH1(idh)
  #endif

  if idx == -1: return idx

  H1HLast = H1head[idx]
  H1Last = H1[idx]
  h1h = H1HLast

  #for ix in range(len(H1h)):

  print("idh: ",h1h[0])
  print("index: ",idx)
  print("title: ",h1h[1])
  print("number of bins: ",h1h[2])
  print("xmin: ",h1h[3])
  print("xmax: ",h1h[4])
  print("dx: ",h1h[5])
  print("ymin: ",h1h[7])
  print("ymax: ",h1h[8])
  print("number of entries: ",h1h[9])
  print("sum: ",h1h[10])
  print("mean x: ",h1h[12])
  print("rms x: ",h1h[13])
  xopt, yopt = h1opt(idh)
  print("xopt: ",xopt)
  print("yopt: ",yopt)
  print("number of underflows: ",h1h[14])
  print("number of overflows: ",h1h[15])
  print("underflows: ",h1h[16])
  print("overflows: ",h1h[17])

def H2Info(idh='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) == str and idh == '?':
    print("\nUsage: H2Info(idh), H2Last, H2I, and Ind are also set, \nreturns -1 if histogram has not been found")
    return

  if type(idh) == int:
    idx = idh
  else:
    idx = GetIndexH2(idh)

  if idx == -1: return idx

  H2Last = H2[idx]
  h2h = H2head[idx]

  print("idh: ",h2h[0])
  print("title: ",h2h[1])

  print("number of x-bins: ",h2h[2])
  print("xmin: ",h2h[3])
  print("xmax: ",h2h[4])
  print("dx: ",h2h[5])

  print("number of y-bins: ",h2h[6])
  print("ymin: ",h2h[7])
  print("ymax: ",h2h[8])
  print("dy: ",h2h[9])

  print("number of entries: ",h2h[10])

  print("zmin: ",h2h[11])
  print("zmax: ",h2h[12])
  print("sum: ",h2h[13])

  print("meanx: ",h2h[14])
  print("rmsx: ",h2h[15])

  print("meany: ",h2h[16])
  print("rmsy: ",h2h[17])

  print("number of underflows in x: ",h2h[18])
  print("number of overflows in x: ",h2h[19])
  print("underflows in x: ",h2h[20])
  print("overflows in x: ",h2h[21])

  print("number of underflows in y: ",h2h[22])
  print("number of overflows in y: ",h2h[23])
  print("underflows in y: ",h2h[24])
  print("overflows in y: ",h2h[25])

  print("number of underflows in x for y < ymin:  ",h2h[26])
  print("number of overflows in x for y < ymin: ",h2h[27])
  print("underflows in x for y < ymin: ",h2h[28])
  print("overflows in x for y < ymin: ",h2h[29])

  print("number of underflows in x for y > ymax:  ",h2h[30])
  print("number of overflows in x for y > ymax: ",h2h[31])
  print("underflows in x for y > ymax: ",h2h[32])
  print("overflows in x for y > ymax: ",h2h[33])

#enddef H2Info(idh='?')

def H1List():
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  print("\n Index in H1, ID, Title,  Nbins Xmin Xmax, dx, ContMin, ContMax, Status")
  print("-------------------------------------------------------------------")

  for i in range(Nh1):
    h1 = H1head[i]
    print(i," ",h1[0]," ",h1[1]," ",h1[2]," ",h1[3]," ",h1[4]," ",h1[5]," ",
          h1[7]," ",h1[8], " ", h1[20])
  #endfor i in range(Nh1):

def H2List():
  # Histograms and Ntuples
  global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
  Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
  H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
  Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
  Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
  H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
  FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman
  print("\n Index in H2, ID, Title  Nxbins Xmin Xmax, Nybins Ymin, Ymax, ContMin, ContMax, Status")
  print("--------------------------------------------------------------")

  for i in range(Nh2):
    h2 = H2head[i]
    print(i," ",h2[0]," ",h2[1]," ",h2[2]," ",h2[3]," ",h2[4]," ",h2[5],\
    " ",h2[6]," ",h2[7]," ",h2[8]," ",h2[9], " ",h2[38])
  #endfor i in range(Nh2):

def hinfo(idh):

  idx = GetIndexH1(idh)

  if idx != -1:
    H1Info(idh)
    return
  #endif

  idx = GetIndexH2(idh)
  if idx != -1:
    H2Info(idh)
    return
  #endif

  print("*** Non-exixsting histogram ***")
  return

#enddef hinfo(idh)

def nentry(nt='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) != Tdf:
    if nt == '?':
      print("\nUsage: index = nentry(nt), Ni and Ind are also set, \nreturns -1 if histogram has not been found")
      return
  #endif type(nt) == np.DataFrame:

  i= GetIndexN(nt)

  Ni = i
  Ind = i

  if i == -1: return i

  nupdate_header(nt)

  nhead = Nhead[i]

  return nhead[-1]
#enddef ninfo(nt='?')

def ninfo(nt='?'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) != Tdf:
    if nt == '?':
      print("\nUsage: index = ninfo(nt), Ni and Ind are also set, \nreturns -1 if histogram has not been found")
      return
  #endif type(nt) == np.DataFrame:

  i= GetIndexN(nt)

  Ni = i
  Ind = i

  if i == -1: return i

  nupdate_header(nt)

  nhead = Nhead[i]

  print("\nIndex and Id:",nhead[0]," ",nhead[1])
  print("Title: ",nhead[2],"\n")
  print("\n Num  Var         Min        Max        Mean       Std")

  nvar = nhead[3]
  Nlines = nhead[-1]

  for k in range(nvar):
    l = 4 + k
    sh = nhead[l]
    try:
      print('{0:3}'.format(k),"  ",'{0:7}'.format(sh[0]), \
      '{0:10}'.format(pg5(sh[1])), '{0:10}'.format(pg5(sh[2])), \
      '{0:10}'.format(pg5(sh[3])), '{0:10}'.format(pg5(sh[4])))
    except:
      print('{0:3}'.format(k),"  ",'{0:7}'.format(sh[0]),
      "   -          -          -          -")
  #endfor k in range(nvar):

  print("\nnumber of entries",nhead[l+1])
#enddef ninfo(nt='?')

def nlist():

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  print("\n Index in Ntup, ID, Title, Number of variables, Nevent")
  print("---------------------------------------------------------")

  for i in range(Nntup):
    nhead = Nhead[i]
    print(nhead[0]," ",nhead[1]," ",nhead[2]," ",nhead[3]," ",nhead[4+nhead[3]])
  #endfor i in range(Nntup):

def NctList():
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  print("\n Index in Nctup, ID Title, Number of Lines")
  for i in range(Nnctup):
    nhead = Nctup[i][0]
    nent = len(Nctup[i][1])
    print(i," ",nhead[0]," ",nhead[1]," ",nent)
  #endfor i in range(Nnctup):

def ncolumns(fname='ntuple.dat', skiphead=-1, sep=' '):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global VarHead

  if fname == '?':
    print("ncols = ncolumns(nt, file='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0, comment='*', sep=' ')")
    return None
  #endif fname == '?'

  VarHead = ''

  flook = open(fname,'r')

  if skiphead < 0:
    while True:
      cline = flook.readline()
      c1 = cline[0]
      if c1 != '*' and c1 != '%' and c1 != '!':
        if sep == ' ':
          ncols = len(cline.split())
        else:
          ncols = len(cline.split(sep))
        #endif
        break
      elif cline[0:12] == '*Variables::':
        varlis= cline.strip().split('::')
        VarHead = varlis[1]
      #endif
    #endwhile
  else:
    for l in range(skiphead+1):
      line = flook.readline().strip()
      if sep == ' ':
        ncols = len(line.split())
      else:
        ncols = len(line.split(sep))
      #endfor l in range(skip+1)
  #endif

  flook.close()

  return ncols

#enddef ncolumns

def ncolumnsguess(fname='ntuple.dat'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global VarHead

  if fname == '?':
    print("nlines, nheader, ncols = ncolumnsguess(file='ntuple.dat')")
    return None
  #endif fname == '?'

  flook = open(fname,'r')
  clines = flook.readlines()
  flook.close()

  ncols = 0
  nmin = 1000000
  nmax = 0
  nncols = 0
  nhead = 0
  llines = len(clines)

  if not llines: return llines, nhead, ncols

  for cline in clines:
    if cline[0:12] == '*Variables::':
      varlis= cline.strip().split('::')
      VarHead = varlis[1]
    #endif
    if cline[0] == '*' or cline[0] == '%' or cline[0] == '!': continue
    n = len(cline.split())
    if n < nmin:
      nmin = n
    elif n > nmax:
      nmax = n
      # else:
    #endif
    nncols += 1
    ncols += n
    #endif
  #endfor

  ncols = int(ncols/nncols+0.5)

  for cline in clines:
    if cline[0] == '*' or cline[0] == '%' or cline[0] == '!' or \
    len(cline.split()) != ncols:
      nhead += 1
    else:
      break
    #endif
  #endfor

  return llines, nhead, ncols

#enddef ncolumns

def nlook(nt='Nlook',fname='ntuple.dat',header=None, skiphead=-1, skipfoot=0,\
silent=0, comment='*', sep=' ',iguessncols=1, iplot=1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if len(nt) == 0 or type(nt)== str and nt == '?':
    print("nt = nlook(nt, file='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0, comment='*', sep=' ',iguessncols=1, iplot=1)")
    print("Consider mlook(...):")
    print("mlook(file='ntuple.dat',header=None, skiphead=0, skipfoot=0, silent=0, comment='*', sep=' ',iguessncols=1,, iplot=1)")
    print("mlook creates Nlook")
    return None
  #endif len(nt) == 0

  if iguessncols: nlines, skiphead,  ncols = ncolumnsguess(fname)
  else: ncols = ncolumns(fname,skiphead,sep)

  cols = 'c1'
  for i in range(2,ncols+1): cols += ':c' + str(i)

  Nt = ncread(nt,cols,fname,header,skiphead,skipfoot,silent,comment,sep)

  if not silent: ninfo(Nt)

  if iplot: npl(Nt,"c1:c2")

  return Nt
#enddef nlook

def mlook(fname='ntuple.dat',header=None, skiphead=-1, skipfoot=0,\
silent=0, comment='*', sep=' ', iguessncols=1, iplot=1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if len(fname) == 0:
    print("mlook(file='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0, comment='*', sep=' ',iplot=1)")
    print("Consider nlook(...):")
    print("nt = nlook(nt, file='ntuple.dat',header=None, skiphead=0, skipfoot=0, silent=0, comment='*', sep=' ', iguessncols=1, iplot=1)")
  #endif len(nt) == 0

  if iguessncols: nlines, skiphead,  ncols = ncolumnsguess(fname)
  else: ncols = ncolumns(fname,skiphead,sep)

  cols = 'c1'
  for i in range(2,ncols+1): cols += ':c' + str(i)

  Nlook = ncread('Nlook',cols,fname,header,skiphead,skipfoot,silent,comment,sep)

  if not silent: ninfo(Nlook)
  if iplot: npl(Nlook,"c1:c2")

  return
#enddef mlook

def nread(nt='',fname='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0,\
comment='*', sep=' '):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if len(nt) == 0:
    print("nt = nread(nt, file='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0, comment='*', sep=' ')")
    return None
  #endif

  ind = GetIndexN(nt)
  if ind == -1:
    print("*** Error in nread: Ntuple not found ***")
    return ind
  #endif ind == -1:

  if fname[0] == '~':
    try: fname = os.environ['HOME'] + fname[1:]
    except: pass
  #endif

  if not os.path.exists(fname):
    if not silent:
      print("\n*** Error nread(...): No file " + fname)
    #endif
    return nt
  #endif

  if not fsize(fname):
    ninfo(nt)
    Quit("hallo",nt,"holla")
    return nt
  #endif

  if skiphead < 0:
    try:
      Fnt = open(fname,'r')
      ncom = 0
      while True:
        cline = Fnt.readline()
        c1 = cline[0]
        if c1 == '*' or c1 == '%' or c1 == '!':
          if cline[0:12] != "*Variables::": print(NL,cline)
          ncom += 1
        else:
          break
        #endif
      #endwhile
      Fnt.close()
    except:
      if not silent:
        print("\n*** Error nread(...): Failed to read " + fname)
      #endif
    #endtry
    skiphead = ncom
  #endif

  try:
    if sep == ' ' or sep == '':
      nt = pd.read_csv(fname,header=header, delim_whitespace=True, skiprows=skiphead, comment='*', skipfooter=skipfoot)
    else:
      nt = pd.read_csv(fname,header=header, sep=sep, skiprows=skiphead, comment='*', skipfooter=skipfoot)
    #endif
  except:
    if not silent:
      print("\n*** Error nread(...): Failed to read " + fname)
    #endif
    return nt
  #endtry

  Ntup[ind] = nt

  nhead = Nhead[ind]
  nvar=nhead[3]

  varlis = []
  for i in range(nvar):
    varlis.append(nhead[4+i][0])
  #endfor i in range(nvar):

  nt.columns = varlis

  for i in range(nvar):
    try:
      nhead[4+i][1] = nt[varlis[i]].min()
    except:
      nhead[4+i][1] = None
    #endtry
    try:
      nhead[4+i][2] = nt[varlis[i]].max()
    except:
      nhead[4+i][2] = None
    #endtry
  #endfor i in range(nvar):

  nhead[4+nvar]=len(nt)
  Nlines = len(nt)

  return nt
#enddef nread

def ncread(ntname='', varlis='', fname='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0, \
comment='*', sep=' ',iguessncols=1, ioverwrite=1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global VarHead
  #reakpoint()

  if len(ntname) == 0:
    print("ncread(ntname='', varlis='x:y:z:bx:by:bz', file='ntuple.dat',header=None, skiphead=-1, skipfoot=0, silent=0, comment='*', sep=' ',iguessncols=1, ioverwrite=1")
    print("If no variables are provided, the first valid line is taken as header.")
    return 0
  #endif

  if fname[0] == '~':
    try: fname = os.environ['HOME'] + fname[1:]
    except: pass
  #endif

  if not os.path.exists(fname):
     if silent == 0:
       print("\n*** Error in ncread(...): File ",fname," not found ***")
     return -1
  #endif

  if skiphead != -1: iguessncols = 0
  fs = fsize(fname)
  varl = nlistcolon(varlis)
  ncol = Ncolon + 1

  if fs:
    try:
      if iguessncols: nlines, skiphead,  ncol = ncolumnsguess(fname)
      else: ncol = ncolumns(fname,skiphead=skiphead,sep=sep)
    except:
      print("\n*** Error ncread: Could not recognise number of columns ***")
      print("\n*** Ntuple and file ",ntname, fname, "***\n")
      return None
    #endtry
  #endif

  if VarHead != '': varlis = VarHead

  spvar = varlis.split(":")
  nvar = len(spvar)

  if nvar == 0:
    varlis=''
    for i in range(ncol-1):
      varlis += "c" + str(i+1) + ":"
    #endfor
    varlis += "c" + str(ncol)
  elif nvar > ncol:
    varlis=''
    for i in range(ncol-1):
      varlis += spvar[i] + ":"
    #endfor
    varlis +=  spvar[ncol-1]
  elif nvar < ncol:
    varlis=''
    for i in range(nvar):
      varlis += spvar[i] + ":"
    #endfor
    for i in range(nvar,ncol-1):
      varlis += "c" + str(i+1) + ":"
    #endfor
    varlis += "c" + str(ncol)
  #endif

  ntc = ncre(ntname,fname,varlis,ioverwrite)
  if not fs: return ntc

  nt=nread(ntname, fname ,header, skiphead, skipfoot, silent, comment, sep)

  return nt

#enddef ncread(ntname='', varlis=[], file='ntuple.dat',header=None, \

def h1list():
  H1List()
#enddef h2list()
def h2list():
  H2List()
#enddef h2list()

def hlist():
  H1List()
  H2List()
#enddef hlist()

def hindex(idh):
  idim = 0
  idx = GetIndexH1(idh)
  if idx == -1:
    idx = GetIndexH2(idh)
    if idx >= 0:
      idim = 2
  else:
    idim = 1
  return idim,idx
#enddef hindex():

def nproj2(nt='?', xy='', weight=1., select='',
           scalex=1., scaley=1., scalez=1.0, nx=51, ny=51, idh=-1,
           ioverwrite=1):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(nt) == str:
    if nt == "" or nt == "?":
      print("\nUsage: nproj2(nt='', xy='', weigth=1. select='', scalex=1., scaley=1., scalez=1., nx=51, ny=51, idh=-1)")
      print("if histogram does not exist, it is created with nx,ny channels")
      return 0
  #endif

  try:
    typ = type(nt)
  except:
    print("*** Non-existing Ntuple in nproj2 ***")
    return -1
  #endtry

  idn = -1

  if typ != Tdf:
    idn = GetIndexN(nt)
    if idn == -1:
      print("*** Non-existing Ntuple in nproj2: " + str(nt) + " ***")
      return -1
    #endif idn == -1
    nt = N
  #endif typ != Tdf

  if len(nt) == 0:
    idn = GetIndexN(nt)
    print("*** Empty Ntuple in nproj2: " + str(Nind[idn]) + " ***")
    return -1
  #endif len(nt) == 0

  if Kecho:
    if type(idh) == str: sidh = "'" + idh + "'"
    elif type(idh) == int: sidh = str(idh)
    else: sidh = "'" + H2head[idx][1] + "'"

    if type(nt) == str: snt = "'" + nt + "'"
    elif type(nt) == int: snt = str(nt)
    elif idn != -1: snt = "'" + Nhead[idn][1] + "'"
    if idn != -1: sj2 = "nproj2(nt='" + snt + "', xy='"
    else: sj2 = "nproj2(nt,"

    print(sj2 + xy + "', weight='" + str(weight)\
    + "', select='" + select + "', scalex=" + str(scalex) + ", scaley="\
    + str(scaley) + ", scalez=" + str(scalez)\
    + ", nx=" + str(nx) + ", ny=" + str(ny) + ", idh=" + sidh + ")" )
  #endif Kecho

  if select:
    nt = nt.query(select)
    Nsel = nt
    select=''
    if len(nt) == 0:
      print("*** No data survived selection in nproj2  ***")
      return -1
    #endif len(nt) == 0
  #endif select:

  if type(xy) == str: xy = nlistcolon(xy)

  if type(xy) != list or len(xy) == 0:
    print("*** Error in nproj2: No list given for variables ***")
    return -2
  #endif type(xy) != list or len(xy) == 0

  x = "(" + nparse(nt,xy[0]) + ") * " + str(scalex)
  y = "(" + nparse(nt,xy[1]) + ") * " + str(scaley)

  if weight == 1. or weight == '' or weight == '!':
    w = "(" + nparse(nt,xy[0]) + ") * 0.0 + 1.0"
    wt = ''
  else:
    w = nparse(nt,weight)
    wt = weight
  #endif weight == 1. or weight == '' or weight == '!'

  scom = "global N; N = pd.DataFrame([" + x + "," + y + "," + w + "]).T"
  exec(scom)

  N.columns = ['x','y','w']
  Nlines = len(N)
  N.index = range(Nlines)

  x = N.x
  y = N.y
  w = N.w

  z2 = x * 0.0
  ent = x * 0.0
  zave = x * 0.0
  ez = x * 0.0

  if type(idh) != int:
    idx = GetIndexH2(idh)
  else:
    idx = idh
  #endif type(idh) != int

  if idx > -1 and ioverwrite and nx > 0 and ny > 0:
    hname = H2head[idx][0]
    hdelete(hname)
  #endif idx > -1:

  if idx == -1 or Khdeleted:

    xmin=x.min(); xmax=x.max()

    if xmin == xmax:
      dx = 1.
      nx = 1
    else:
      if nx == 1:
        dx = (xmax-xmin)
      else:
        dx = (xmax-xmin)/(nx-1)
    #endif xmin == xmax

    xmin -= dx/2.; xmax += dx/2.

    ymin=y.min(); ymax=y.max()

    if ymin == ymax:
      dy = 1.
      ny = 1
    else:
      if ny == 1:
        dy = (ymax-ymin)
      else:
        dy = (ymax-ymin)/(ny-1)
    #endif ymin == ymax

    ymin -= dy/2.; ymax += dy/2.

    if type(idh) == str: hname = idh
    else: hname = 'HnPlot'

    htit = hname + " :: " + xy[0] + ":" + xy[1] + ":" + wt

    if type(weight) == str and len(weight):
      htit += ":" + weight
    elif weight == 1.0:
      pass
    else:
      htit += ":" + str(weight)
    #endif type(weight) == str and if len(weight):

    if len(select): htit += " [" + select + "]"

    hret = hbook2(hname,htit,nx,xmin,xmax,ny,ymin,ymax)

    idx = GetIndexH2(hname)

    if hname != 'HnPlot' and Kecho:
      print("--- Have created histogram: ", idx, " ",hname," ",Nhead[idn][2])

  #endif idx == -1

  nx = H2hh[2]
  xmin = H2hh[3]
  xmax = H2hh[4]
  dx = H2hh[5]

  ny = H2hh[6]
  ymin = H2hh[7]
  ymax = H2hh[8]
  dy = H2hh[9]

  hn, xedg, yedg  = np.histogram2d(x, y, bins=[nx,ny], range=[[xmin,xmax],[ymin,ymax]])
  hz, xedg, yedg  = np.histogram2d(x, y, bins=[nx,ny], range=[[xmin,xmax],[ymin,ymax]], weights=w)
  hz2, xedg, yedg  = np.histogram2d(x, y, bins=[nx,ny], range=[[xmin,xmax],[ymin,ymax]], weights=w*w)

  xl = np.linspace(xmin+dx/2.,xmax-dx/2.,nx)
  yl = np.linspace(ymin+dy/2.,ymax-dy/2.,ny)

  hz = hz * scalez
  hz2 = hz2 * scalez * scalez

  x,y = np.meshgrid(xl,yl)
  x = x.T.flatten()
  y = y.T.flatten()
  hn = hn.flatten()
  hz = hz.flatten()
  hz2 = hz2.flatten()

  h = pd.DataFrame([x,y,hz,hz2,hn]).T
  h.columns=['x','y','z','z2','n']

  h.z[np.isnan(h.z)] = 0.0
  h.z2[np.isnan(h.z2)] = 0.0
  h.n[np.isnan(h.n)] = 0.0

  h['ave'] = h.z/h.n
  h.ave[np.isnan(h.ave)] = 0.0

  h['ez'] = (h.z2/h.n-h.ave**2)**0.5
  h.ez[np.isnan(h.ez)] = 0.0

  head2 = H2head[idx]

  sumz = h.z.sum()

  head2[10] = h.n.sum()
  head2[11] = h.z.min()
  head2[12] = h.z.max()
  head2[13] = sumz
  xmean = (h.x*h.z).sum()/sumz
  head2[14] = xmean
  head2[15] = max(0.0,(h.x**2*h.z).sum()/sumz-xmean**2)**0.5
  ymean = (h.y*h.z).sum()/sumz
  head2[16] = ymean
  head2[17] = max(0.0,(h.y**2*h.z).sum()/sumz-ymean**2)**0.5

  xmin = head2[3]
  xmax = head2[4]
  ymin = head2[7]
  ymax = head2[8]

  Nux = N.query('x<' + str(xmin))
  nux = len(Nux)
  ux = Nux.w.sum()

  Nuxoy = Nux.query('y>' + str(ymax))
  nuxoy = len(Nux)
  uxoy = Nux.w.sum()

  Nuxuy = Nux.query('y<' + str(ymin))
  nuxuy = len(Nuxuy)
  uxuy = Nuxuy.w.sum()

  Noy = N.query('y>' + str(ymax))
  noy = len(Noy)
  oy = Noy.w.sum()

  Nox = N.query('x>' + str(xmax))
  nox = len(Nox)
  ox = Nox.w.sum()

  Noxov = Nox.query('y>' + str(ymax))
  noxoy = len(Nox)
  oxoy = Nox.w.sum()

  Noxuy = Nox.query('y<' + str(ymin))
  noxuy = len(Noxuy)
  oxuy = Noxuy.w.sum()

  Nuy = N.query('y<' + str(ymin))
  nuy = len(Nuy)
  uy = Nuy.w.sum()

  head2[18] = nux
  head2[19] = nox
  head2[20] = ux
  head2[21] = ox
  head2[22] = nuy
  head2[23] = noy
  head2[24] = uy
  head2[25] = oy
  head2[26] = nuxuy
  head2[27] = noxuy
  head2[28] = uxuy
  head2[29] = oxuy
  head2[30] = nuxoy
  head2[31] = noxoy
  head2[32] = oxuy
  head2[33] = oxoy

  H2h=h
  H2[idx]=h

  return 0
#enddef nproj2()

def nproj1(nt='?', var='', weight=1., select='', scalex=1., scaley = 1,
           nx=101, idh=-1, ioverwrite=1):

  import numpy as np
  import pandas as pd

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  #reakpoint()
  idn = -1

  varl = nlistcolon(var)

  if len(varl) > 1:
    if weight != 1.:
      print("*** Error in nproj1: Weight not 1. and var contains more then one item ***")
      return -1
    #endif
    weight = varl[1]
    var = varl[0]
  #endif

  try:
    typ = type(nt)
  except:
    print("*** Non-existing Ntuple in nproj1 ***")
    return -1
  #endtry

  if typ == str:
    if nt == '' or nt == '?':
      print("nproj1(nt='', var='', weigth=1. select='', scalex=1., scaley=1. nx=100, idh=-1)")
      print("if histogram does not exist, it is created with nx channels")
      print("if nx <= 0: nx = min(100,len(nt))")
      return 0
  #endif type(nt) != Tdf and len(nt) == 0

  if typ != Tdf:
    idn = GetIndexN(nt)
    if idn == -1:
      print("*** Non-existing Ntuple in nproj1: " + str(nt) + "  ***")
      return -1
    #endif idn == -1
    nt = N
  #endif typ != Tdf

  if len(nt) == 0:
    idn = GetIndexN(nt)
    print("*** Empty Ntuple in nproj1: " + str(Nind[idn]) + " ***")
    return -1
  #if len(nt) == 0

  idx = GetIndexH1(idh)

  if idx > -1:
    hname = H1head[idx][0]
    if ioverwrite and nx > 0: h1reset(hname)
  #endif idx > -1:

  if Kecho:

    if type(idh) == str: sidh = "'" + idh + "'"
    elif type(idh) == int: sidh = str(idh)
    else: sidh = "'" + H1head[idx][1] + "'"

    if type(nt) == str: snt = "'" + nt
    elif type(nt) == int: snt = str(nt)
    elif idn != -1: snt = "'" + Nhead[idn][1] + "'"
    if idn != -1: sj1 = "nproj1(nt=" + snt + ", var='"
    else: sj1 = "nproj1(nt, '"

    print(sj1 + var + "', weight='" + str(weight)\
    + "', select='" + select + "', scalex=" + str(scalex) + ", scaley="\
    + str(scaley) + ", nx=" + str(nx) + ", idh=" + sidh + ", ioverwrite=" + str(ioverwrite) + ")" )
  #endif Kecho

  if type(var) == list:
    if len(var) != 1:
      print("*** Error in nproj1: Length of variable list not one ***")
      return -1
    else:
      var = var[0] * scalex
  #endif len(var) != 1

  if len(var) == 0:
      print("*** Error in nproj1: No variable specified ***")
      return -1
  #endif len(var) == 0

  Nsel = nt
  if select:
    nt = nt.query(select)
    Nsel = nt
    if len(nt) == 0:
      print("*** No data survived selection in nproj1  ***")
      return -1
    #endif len(nt) == 0
    select=''
  #endif select:

  nt = Nsel
  x = "(" + nparse(nt,var) + ") * " + str(scalex)

  try: w = str(float(weight))
  except: w = weight

  if w == "1.0" or w == '' or w == '!':
    w = "(" + nparse(nt,var) + ") * 0.0 + 1.0"
  else:
    w = nparse(nt,weight)
  #endif weight == 1. or weight == '' or weight == '!'

  scom = "global N; N = pd.DataFrame([" + x + "," + w + "]).T"

  exec(scom)

  N.columns = ['x','w']
  Nlines = len(N)
  N.index = range(Nlines)

  x = N.x
  w = N.w

  if nx <= 0: nx = min(101,len(N.x.drop_duplicates()))

  #reakpoint()
  if idx == -1 or Khdeleted:

    xmin=x.min(); xmax=x.max()

    if xmin == xmax:
      if xmin == 0.0:
        dx = 1.
        nx = 1
        xmin -= dx/2.; xmax += dx/2.;
      else:
        dx = abs(xmax)*0.05
        nx = 1
        xmin -= dx/2.; xmax += dx/2.;
    else:
      if nx == 1:
        dx = (xmax-xmin)
      else:
        dx = (xmax-xmin)/(nx-1)
      #endif nx == 1
    #endif xmin == xmax

    if type(idh) == str: hname = idh
    else: hname = 'HnPlot'

    htit = hname + " :: " + var

    if type(weight) == str:
      if len(weight) > 0: htit += ":" + weight
    elif weight != 1.0:
      htit += ":" + str(weight)
    #endif

    if len(select): htit += " [" + select + "]"

    hret = hbook1(hname,htit,int(nx+0.5),xmin,xmax,ioverwrite)

    if hname != 'HnPlot' and Kecho:
      print("--- Have created histogram: ", idx, " ",hname," ",htit)

  #endif idx == -1

  idx = GetIndexH1(hname)

  nx = H1hh[2]
  xmin = H1hh[3]
  xmax = H1hh[4]
  dx = H1hh[5]

  hy, xedg = np.histogram(x, bins=nx, range=[xmin,xmax], weights=w)
  hy2, xedg = np.histogram(x, bins=nx, range=[xmin,xmax], weights=w*w)
  hn, xedg = np.histogram(x, bins=nx, range=[xmin,xmax])

  hy = hy * scaley
  hy2 = hy2 * scaley*scaley

  x = np.linspace(xmin+dx/2.,xmax-dx/2.,nx)
  h = pd.DataFrame([x,hy,hy2,hn]).T

  h.columns=['x','y','y2','n']

  h.y[np.isnan(h.y)] = 0.0
  h.y2[np.isnan(h.y2)] = 0.0
  h.n[np.isnan(h.n)] = 0.0

  h['ave'] = h.y/h.n
  h.ave[np.isnan(h.ave)] = 0.0

  if w.min() == w.max() and w.min() == 1.:
    h['ey'] = h.n**0.5
  else:
    h['ey'] = (h.y2/h.n-h.ave**2)**0.5
  #endif

  h.ey[np.isnan(h.ey)] = 0.0

  H1h = h
  H1[idx] = h

  head1 = H1head[idx]
  head1[7] = min(h.y)
  head1[8] = max(h.y)
  head1[9] = h.n.sum()

  sumy = h.y.sum()
  ya = abs(h.y)
  sumya = ya.sum()

  xmean = 0.
  xrms = 0.
  xsqsum = 0.

  if sumya:
    xsqsum = ((h.x)**2*ya).sum()
    xmean = (h.x*ya).sum()/sumya
    xrms = np.sqrt(max(0.0,xsqsum/sumya-xmean**2))
  #endif

  head1[10] = sumy

  head1[11] = xsqsum
  head1[12] = xmean
  head1[13] = xrms

  quer = 'x < ' + str(xmin)
  hund = h.query(quer)
  nunder = hund.n.sum()
  under = hund.y.sum()

  quer = 'x >= ' + str(xmax)
  hov = h.query(quer)
  nover = hov.n.sum()
  over = hov.y.sum()

  head1[14] = nunder
  head1[16] = under
  head1[15] = nover
  head1[17] = over

  head1[18] = h.x.sum()
  head1[19] = ((h.x)**2).sum()

#enddef nproj1()

def setnoempty(no=1):
  global Inoempty
  Inoempty = no
  return
#enddef setnoempty(no=1)

def getnoempty():
  global Inoempty
  return Inoempty
#enddef getnoempty

def hstat1d(idh='?'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  global Debug, Ical

  NxBinMax=0

  if type(idh) == str and idh == '?':
    print("\nUsage: sumn,sumy,xmean,xrms,xopt,yopt = hstat1d(idh='?')")
    return
  #endif

  idx = GetIndexH1(idh)
  if idx == -1:
    print('*** Error in hstat1d: Non-existing histogram ',idh)
    return -1
  elif H1hh[20]:
    print('*** Warning in hplot1d: Histogram has been marked as deleted ',idh)
  #endif idx == -1:

  if Kecho:
    print("hstat1d(" + str(idh) + ")")
  #endif Kecho

  hq = H1[idx]

  x = hq['x']
  y = hq['y']
  ya = abs(y)
  yave = hq['ave']
  ey = hq['ey']
  ny = hq['n']
  ly = len(y) - 1

  sumy = hq.y.sum()
  sumn = hq.n.sum()
  ymaxa = abs(y).max()
  sumya = ya.sum()

  xmean = -9999.
  xrms = -9999.
  xsqsum = 0.
  lstat = 1

  if sumya:
    xsqsum = (x**2*ya).sum()
    xmean = (x*ya).sum()/sumya
    xrms = np.sqrt(max(0.0,xsqsum/sumya-xmean**2))
  #endif sumy

  # 24.9.2021 stdy = ey/ny # RMS of bin
  stdy = deepcopy(ey) # RMS of bin
  stdy[np.isnan(stdy)] = 0.0

  # 24.9.2021 stdyprof = ey/ny/(ny-1)**0.5 # estimated error of mean
  stdyprof = stdy/(ny-1)**0.5 # estimated error of mean
  stdyprof[np.isnan(stdyprof)] = 0.0

  xopt, yopt = h1opt(hq)

  return [sumn,sumy,xmean,xrms,xopt,yopt]

#enddef hstat1d(idh,plopt='')

def vstat(x='?',y=''):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  #nreakpoint()

  if type(x) == str:
    print("\nvstat(x,y) returns [xmin, xmax, xmean, xrms, xopt, yopt]")
    print("If y is missing, y is treated as unity.\n")
    return
  #if type(x) == str

  ione = 0
  if type(y) == str:
    ione = 1
    y = x*0 + 1.
  #endif type(y) == str

  ya = abs(y)
  ny = len(y)
  ly = len(y) - 1

  sumy = y.sum()
  ymaxa = abs(y).max()
  sumya = ya.sum()

  xmean = -9999.
  xrms = -9999.
  xsqsum = 0.

  if sumya:
    xsqsum = (x**2*ya).sum()
    xmean = (x*ya).sum()/sumya
    xrms = np.sqrt(max(0.0,xsqsum/sumya-xmean**2))
  #endif sumy

  xopt = None
  yopt = None

  if not ione: xopt, yopt = voptpar(x,y)

  return [x.min(),x.max(),xmean,xrms,xopt,yopt]

#enddef vstatxy(x,y)

def vmean(x='?',y=''):
  res = vstat(x,y)
  return res[2]
#enddef

def vrms(x='?',y=''):
  res = vstat(x,y)
  return res[3]
#enddef

def hplot1d(idh='?', plopt='2d', Tit='!', xTit='', yTit='', legend='',
            barwidth=-9, block=False):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  global Debug, Ical

  NxBinMax=0

  if type(idh) == str and idh == '?':
    print("\nUsage: hplot1d(idh='?', plopt='2d', Tit='!', xTit='', yTit='', legend='', barwidth=-9, block=False): ")
    return
  #endif type(idh) == str and idh == '?'

  idx = GetIndexH1(idh)

  if idx == -1:
    print('*** Error in hplot1d: Non-existing histogram ',idh)
    return -1
  elif H1hh[20]:
    print('*** Warning in hplot1d: Histogram has been marked as deleted ',idh)
  #endif idx == -1:

  if Kecho:
    s = "hplot1d('" + H1head[idx][0] + "', plopt='" + plopt \
    + "', Tit='" + Tit + "', xTit='" + xTit + "', yTit='" + yTit \
    + "', legend='" + legend + "', barwidth=" + str(barwidth) \
    + ", block=" + str(block) \
    + ")"
    print(s)
  #endif Kecho

  plotopt(plopt)

  if Isame:
    Tit = ''
    xTit = ''
    yTit = ''
  #endif Isame

  getzone()

  if Isame:
    querx = 'x >= ' + str(ZoomXmin) + ' and x <= ' + str(ZoomXmax)
    if Iprof:
      query = 'ave >= ' + str(ZoomYmin) + ' and ave <= ' + str(ZoomYmax)
    else:
      query = 'y >= ' + str(ZoomYmin) + ' and y <= ' + str(ZoomYmax)
    #endif
    squer = querx + " and " + query
    if Inoempty: squer += " and y != 0.0"
    hq = H1[idx].query(squer)
  else:
    if not Inoempty:
      hq = H1[idx]
    else:
      hq = H1[idx].query("y != 0.0")
    #endif not Inoempty
  #endif Isame

  if len(hq) == 0 and len(H1[idx]) != 0:
    print("--- Warning in hplot1d: Nothing to plot, due to zoom clipping ---")
    return
  #endif len(hq) == 0 and len(H1[idx]) != 0

  hq = hq.reset_index()

  x = hq['x']
  y = hq['y']
  ya = abs(y)
  yave = hq['ave']
  ey = hq['ey']
  ny = hq['n']
  ly = len(y) - 1

  sumy = hq.y.sum()
  sumn = hq.n.sum()
  ymaxa = abs(y).max()
  sumya = ya.sum()

  xmean = -9999.
  xrms = -9999.
  xsqsum = 0.
  lstat = 1

  if sumya:
    xsqsum = (x**2*ya).sum()
    xmean = (x*ya).sum()/sumya
    xrms = np.sqrt(max(0.0,xsqsum/sumya-xmean**2))
  #endif sumy

  # 24.9.2021 stdy = ey/ny # RMS of bin
  stdy = deepcopy(ey) # RMS of bin
  stdy[np.isnan(stdy)] = 0.0
  # 24.9.2021 stdyprof = ey/ny/(ny-1)**0.5 # estimated error of mean
  stdyprof = stdy/(ny-1)**0.5 # estimated error of mean
  stdyprof[np.isnan(stdyprof)] = 0.0
  stdyprof[np.isinf(stdyprof)] = 0.0

  nx = len(x)
  xmin = x.min()
  xmax = x.max()
  dx = (xmax-xmin) / max(nx,1)

  if re.search('2d',plopt.lower()): plopt = Mode2d

  Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist)
  plt.rcParams['axes.labelsize'] = Atitfontsize

  #kfig = getzone()

  if re.search('2d',plopt): plopt = Mode2d

  if plopt == '' or plopt == 'same' or plopt == 'S':
    if plopt == 'same' or plopt == 'S':
      if ey.max() == 0:
        plopt = 'h'
      else:
        plopt = 'e'
      #endif ey.max() == 0
    else:
      if ey.max() == 0:
        plopt = 'hsame'
      else:
        plopt = 'errsame'
      #endif ey.max() == 0
    #endif
  #endif plopt == '':

  if plopt == '!':
    if nx <= 101: plopt = 'h'
    else: plopt = 'line'
  #endif plopt == '!'

  plotopt(plopt)

  if barwidth == -9: barwidth = H1head[idx][5] * 1.0
  if type(barwidth) == str and barwidth == '!': barwidth = Histbarwidth

  if nx > 1 and dx/max(abs(xmax),abs(xmin)) < 1.e-5: NxBinMax = 5

  iplot = 0

  if plopt == 'e' or Ierr:

    lincol = Markercolor

    xpl = []
    ypl = []
    epl = []

    for i in range(len(x)):
      if ny[i] == 0: continue
      xpl.append(x[i])
      if Iprof:
        ypl.append(yave[i])
      else:
        ypl.append(y[i])
      #endif
      epl.append(ey[i])
    #endfor

    plt.errorbar(xpl, ypl, epl, ls='',marker=Markertype,fillstyle=Fillstyle, mfc=Markercolor, mec=Markercolor, ms=Markersize, mew=1, c=lincol)

    iplot=1

    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexyz(x,yave,stdy,fout)
      print("\nData written to ",fout)
    #endif

  elif Iprof:

    lincol = Markercolor

    xpl = []
    ypl = []
    epl = []

    for i in range(len(x)):
      if ny[i] == 0: continue
      xpl.append(x[i])
      ypl.append(yave[i])
      epl.append(stdyprof[i])
    #endfor
    plt.errorbar(xpl,ypl,epl, ls='',marker=Markertype,fillstyle=Fillstyle, mfc=Markercolor, mec=Markercolor, ms=Markersize, mew=1, c=lincol)

    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexyz(x,yave,stdyprof,fout)
      print("\nData written to ",fout)
    #endif

    iplot=1

  elif plopt == 'h' or Ihist:
    plt.bar(x,y,width=barwidth,color=Histcolor,edgecolor=Histedgecolor, fill=Ifill1d)
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexy(x,y,fout)
      print("\nData written to ",fout)
    #endif
    iplot=1
  elif plopt == 'he' or Ihist and Ierr:
    plt.bar(x,y,yerr=ey,width=barwidth,color=Histcolor,edgecolor=Histedgecolor)
    iplot=1
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexyz(x,y,ey,fout)
      print("\nData written to ",fout)
    #endif
  elif Iline and Imarker:
    plt.plot(x,y,c=Linecolor,ls=Linestyle,lw=Linewidth,marker=Markertype,fillstyle=Fillstyle, mfc=Markercolor, mec=Markercolor, ms=Markersize, mew=1)
    iplot=1
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexy(x,y,fout)
      print("\nData written to ",fout)
    #endif
  elif Ispline:
    xmin = x.min(); xmax = x.max()
    xspl = vec(xmin,xmax,1001)
    yspl = vspline(x,y,xspl)
    plt.plot(xspl,yspl,c=Linecolor,ls=Linestyle,lw=Linewidth)
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexy(splx,yspl,fout)
      print("\nData written to ",fout)
    #endif
    if Imarker: plt.plot(x,y,ls='',marker=Markertype,fillstyle=Fillstyle, mfc=Markercolor, mec=Markercolor, ms=Markersize, mew=1)
    iplot=1
  elif Imarker:
    plt.plot(x,y,ls='',marker=Markertype,fillstyle=Fillstyle, mfc=Markercolor, mec=Markercolor, ms=Markersize, mew=1)
    iplot=1
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexy(x,y,fout)
      print("\nData written to ",fout)
    #endif
  elif Iline:
    plt.plot(x,y,c=Linecolor,ls=Linestyle,lw=Linewidth)
    iplot=1
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexy(x,y,fout)
      print("\nData written to ",fout)
    #endif
  elif Imarker:
    plt.plot(x,y,ls='',marker=Markertype,fillstyle=Fillstyle, mfc=Markercolor, mec=Markercolor, ms=Markersize, mew=1)
    iplot=1
    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      vwritexy(x,y,fout)
      print("\nData written to ",fout)
    #endif Kdump
  #endif plopt == 'e' or Ierr

  if iplot == 0:
    plt.plot(x,y)
  #endif plopt

  if Tit == "!": Tit = H1head[idx][1]

  if Kstat and lstat:

    if StatFontSize < 0:
      dpi = Fig.dpi
      nxy = max(Nxzone,Nyzone)
      sfs = 12 - nxy * 2
    else:
      sfs = StatFontSize
    #endif StatFontSize < 0

    if not Iprof:
      xopt, yopt = h1opt(hq)
    else:
      res = vstat(hq.x,hq.ave)
      xopt = res[4]
      yopt = res[5]
    #endif Iprof

    tex = "N: " + str(int(sumn)) + \
    "\nSum: " + '{:.4g}'.format(sumy) + \
    "\nMean: " + '{:.4g}'.format(xmean) + \
    "\nRMS: " + '{:.4g}'.format(xrms)

    if yopt:
      tex += "\nxOpt: " + '{:.4g}'.format(xopt) + \
      "\nyOpt: " + '{:.4g}'.format(yopt)
    #endif

    text(Xstat,Ystat,tex,halign='left')
    # Latex makes trouble with exponents
    #    latex(Xstat,Ystat,tex,color='black',fontsize=sfs,halign='left',valign='top', \
    #    bbox='none',bbstyle='round',fc='white',ec='black',pad=0.2)

  #endif Kstat

  if len(legend): Legend.append(legend)

  txyz(Tit,xTit,yTit)
  showplot()

  if Kpdf:
    Npdf += 1
    fout = WaveFilePrefix + str(Npdf) + ".pdf"
    pplot(fout,'A4','landscape')
    #print("\nPlot written to ",fout)
  #endif Kdump

  Kplots[Kzone-1] = 1

#enddef hplot1d(idh,plopt='')

def hplot(idh, plopt='!', Tit='!', xTit='', yTit='', zTit = '', legend='', block=False):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  idx = GetIndexH1(idh)

  if idx >= 0:
    hplot1d(idh,plopt,Tit,xTit,yTit,legend,block)
  else:
    idx = GetIndexH2(idh)
    if idx == -1:
      print('*** Error in hplot: Non-existing histogram ',idh)
      return -1
    #endif idx == -1:
    hplot2d(idh,plopt,block)
  #endif idx >= 0:

#enddef hplot(idh,plopt='h')

def hplave(idh, plopt='!', Tit='!', xTit='', yTit='', zTit = '', legend='', block=False):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  idx = GetIndexH1(idh)

  if idx >= 0:
    h1 = hget(idh)
    h1y = deepcopy(h1.y)
    h1.y = h1.ave
    hplot1d(idh,plopt,Tit,xTit,yTit,legend,block)
    h1.y = h1y
  else:
    idx = GetIndexH2(idh)
    if idx == -1:
      print('*** Error in hplot: Non-existing histogram ',idh)
      return -1
    #endif idx == -1:
    h2 = hget(idh)
    h2z = deepcopy(h2.z)
    h2.z = h2.ave
    hplot2d(idh,plopt,block)
    h2.z = h2z
  #endif idx >= 0:

#enddef hplave(idh,plopt='h')

def geo_proz(sgeo):

  global ScreenWidth, ScreenHeight

  swid = int(ScreenWidth)
  shei = int(ScreenHeight)

  wh,x,y = sgeo.split('+')
  w,h = wh.split('x')

  if w[-1] == '%': w = int(swid*float(w[:-1])/100.)
  if h[-1] == '%': h = int(shei*float(h[:-1])/100.)
  if x[-1] == '%': x = int(swid*float(x[:-1])/100.)
  if y[-1] == '%': y = int(shei*float(y[:-1])/100.)

  return str(w) + 'x' + str(h) + '+' + str(x) + '+' + str(y)

#enddef

def read_window_geometry(fname='ntupplot.cfg'):

  global Figgeom, Figgeoms, Figgeom2, FiggeomR, FiggeomL, Fwidth, Fheight, Fxoff, Fyoff, XtermGeo

  global WavesMode

  if WavesMode == 'WAVES' or WavesMode == 'WPLOT': fname = 'waveplot.cfg'
  elif WavesMode == 'UNDUMAG': fname = 'undugui.cfg'
  elif WavesMode == 'BRILL': fname = 'brill.cfg'
  else:
    if fexist("waveplot.cfg"): fname = 'waveplot.cfg'
    else: fname = 'ntupplot.cfg'
  #endif

  #print("WavesMode:",WavesMode)

  if not Figgeom:

    Figgeom =  "50%x60%+5%+5%"
    Figgeom2 = "40%x40%+5%+5%"
    FiggeomR = "40%x40%+55%+5%"
    FiggeomL = "40%x40%+5%+5%"
    XtermGeo = "30%x20%+60%+50%"

    Figgeoms.append(Figgeom)
    Figgeoms.append(Figgeom2)
    Figgeoms.append(FiggeomR)
    Figgeoms.append(FiggeomL)

    if os.path.exists(fname):
      Figgeoms = []
      Fgeom = open(fname,"r")
      Figgeoms.append(Fgeom.readline().split('!')[0].strip())
      Figgeoms.append(Fgeom.readline().split('!')[0].strip())
      Figgeoms.append(Fgeom.readline().split('!')[0].strip())
      Figgeoms.append(Fgeom.readline().split('!')[0].strip())
      Fgeom.close()
      Figgeom =  Figgeoms[0].split('!')[0].strip()
      Figgeom2 = Figgeoms[1].split('!')[0].strip()
      FiggeomR = Figgeoms[2].split('!')[0].strip()
      FiggeomL = Figgeoms[3].split('!')[0].strip()
    #endif os.path.exists("'ntupplot.cfg'")

  #endif not Figgeom

#def read_window_geometry()

def window_geometry(geom='', fig=-1, set=True):
  global Figgeom, Fig, Tfig, Figman, ScreenWidth, ScreenHeight, Kplots, Wmaster

  if geom == '' or geom == '!':
    if Figgeom == '':
      read_window_geometry()
    #endif Figgeom != '':
    geom = Figgeom
  #endif geom == ''

  geom = geo_proz(geom)

  if type(fig) == int and fig == -1:
    fig = plt.gcf()
    Tfig = type(fig)
    Fig = fig
    Figman =  plt.get_current_fig_manager()
  #endif type(fig) == int and fig == -1

  if set:
    fig.canvas.manager.window.wm_geometry(geom)
  else:
    geom = fig.canvas.manager.window.wm_geometry()
  #endif set

  return geom
#enddef window_geometry(geom='700x600+1200+100', fig=-1)

def get_window_geometry(fig=-1):
  return window_geometry('', fig, False)

def get_window_width(fig=-1):
  geo = get_window_geometry(fig)
  s = geo.split('x')
  return s[0]

def get_window_height(fig=-1):
  geo = get_window_geometry(fig)
  s = geo.split('x')
  h = s[1].split('+')
  return h[0]

def window_get_title(fig=-1):
  global Fig, Tfig, Figman

  if type(fig) == int and fig == -1:
    fig = plt.gcf()
    Tfig = type(fig)
    Fig = fig
  #endif

  Figman = plt.get_current_fig_manager()
  return Figman.get_window_title()

#enddef window_get_title()

def window_set_title(Title='',fig=-1):

  global Fig, Tfig, Figman

  if type(fig) == int and fig == -1:
    fig = plt.gcf()
    Tfig = type(fig)
    Fig = fig
  #if type(fig) == int and fig == -1

  #fig.canvas.set_window_title(Title)
  Figman = plt.get_current_fig_manager()
  Figman.set_window_title(Title)

#enddef window_set_title()

def gui_key_press(ev):
  if ev.key in ['q', 'Q']: Quit()
#enddef

def window(title='', geom="!", block=False, projection = '2d',
           getconsole=True, visible=False):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  global Tfig, Tax2d, Tax3d, IsameGlobal, ScreenWidth, ScreenHeight, Tdate, \
  Figman, Wmaster

  #nreakpoint()

  if Nwins <= 0:
    read_window_geometry()
    Nfigs = len(plt.get_fignums())
    #if Nfigs > Nwins: Nwins = Nfigs
    Nwins = Nfigs
  #endif

  Nfigs+=1
  Nwins+=1
  Kplots[Nwins-1] = 0

  Kfig = Nwins

  tit = title

  if title == '?':
    print("title='?', geom='700x600+1200+100', block=False")
  elif title == '':
    tit = 'Win_' + str(Nwins)
  #endif title == '?':

  Fig = plt.figure(tit)
  #Fig = plt.gcf()
  Tfig = type(Fig)
  Figman =  plt.get_current_fig_manager()

  w,h = getplotsize()
#  print(w,h,ScaleSizeX,ScaleSizeY)
  setplotsize(Fig,w*ScaleSizeX,h*ScaleSizeY)
#  siz = Fig.get_size_inches()
#  Fig.set_size_inches(siz[1]*ScaleSizeX, siz[0]*ScaleSizeY)

  MPLmain = Fig
  MPLmaster = MPLmain.canvas.toolbar.master
  Wmaster = MPLmaster

  plt.connect('key_press_event', gui_key_press)

  ScreenWidth = MPLmaster.winfo_screenwidth()
  ScreenHeight = MPLmaster.winfo_screenheight()

  if type(Klegend) == int or type(Klegend) == bool:
    iledg = Klegend
    Klegend = IntVar(MPLmaster)
    Klegend.set(iledg)
  #if type(Klegend) == int or type(Klegend) == bool

  if type(Icmap) == int:
    icmap = Icmap
    Icmap = IntVar(MPLmaster)
    Icmap.set(icmap)
  #endif type(Icmap) == int

  if type(IsameGlobal) == int:
    isame = IsameGlobal
    IsameGlobal = IntVar(MPLmaster)
    IsameGlobal.set(isame)
  #endif type(IsameGlobal) == int

  geo = window_geometry(geom,Fig)

  if projection.lower() != '3d':
    Ax = Fig.add_subplot(111,visible=visible,label='111')
    #Ax = Fig.add_subplot(111,visible=False,label='111')
    Tax2d = type(Ax)
  else:
    Ax = Fig.add_subplot(111,projection='3d',visible=visible,label='111')
    #Ax = Fig.add_subplot(111,projection='3d',visible=False,label='111')
    Tax3d = type(Ax)
  #endif projection.lower() != '3d'

  #Figs.append(Fig)

  Axes = Fig.axes
  Kzone = 1
  Nyzone = 1
  Nxzone = 1

  if len(Figs) < Nwins:
    Figs.append(Fig)
  else:
    Figs[Nwins-1] = Fig
  #endif len(Figs) < Nwins

  z = []
  z.append(Fig)
  z.append(Nxzone)
  z.append(Nyzone)
  z.append(Kzone)
  z.append(1) # date on figure
  Zones.append(z)

  date_on_figure()
  if Krun != None: run_on_figure()

  shpl()
  return

#enddef window(title='?', xoff=0., yoff=0., wid = 700., hig=500.)

def win2(title='Win_2', geom="!", block=False, projection = '2d', getconsole=True, visible=True):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  if geom == '!': geom = Figgeom2
  window(title=title, geom=geom, block=block, projection =projection, getconsole=getconsole, visible=visible)
  shpl()
#enddef win2(title='', geom="!", block=False, projection = '2d')

def winr(title='Win_r', geom="!", block=False, projection = '2d', getconsole=True, visible=False):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  if Nwins < 1: read_window_geometry(fname='ntupplot.cfg')
  if geom == '!':
    if FiggeomR: geom = FiggeomR
    else: FiggeomR = Figgeom
  #endif
  window(title=title, geom=geom, block=block, projection =projection, getconsole=getconsole, visible=visible)
  shpl()
#enddef winr(title='', geom="!", block=False, projection = '2d')

def winl(title='Win_l', geom="!", block=False, projection = '2d', getconsole=True, visible=False):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  if Nwins < 1: read_window_geometry(fname='ntupplot.cfg')
  if geom == '!': geom = FiggeomL
  window(title=title, geom=geom, block=block, projection =projection, getconsole=getconsole, visible=visible)
  shpl()
#enddef winl(title='', geom="!", block=False, projection = '2d')

def showplot(visible=True):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global NxBinMax, Ical, Aspect, Figman

  Nfigs = len(plt.get_fignums())
  if not Nfigs: return

  if Nwins <= 0 and Nfigs:
    plt.grid(Kgrid)
    if Igetconsole: get_console()
    plt.show(block=False)
    return
  #endif

  Isame = 0
  is3d = 0

  if hasattr(Ax,'zaxis'): is3d = 1

  if Itight:
    plt.tight_layout(pad=Tightpad, w_pad=Xtightpad, h_pad=Ytightpad)
  else:
    plt.subplots_adjust(left=LeftMargin, right=RightMargin, bottom=BottomMargin, top=TopMargin, wspace=Xspace, hspace=Yspace)
  #endif Itight

  if Icallfromoverview:
    if LogX: Ax.set_xscale('log')
    else: Ax.set_xscale('linear')
    if LogY: Ax.set_yscale('log')
    else: Ax.set_yscale('linear')
    #Ax.set_visible(visible)
    #plt.show(block=False)
    xmin,xmax = Ax.get_xlim()
    if (xmax+xmin) > 0 and (xmax-xmin)/(xmax+xmin) < 1.e-3: NxBinMax = 5
    if NxBinMax:
      plt.locator_params(axis='x', nbins=NxBinMax)
    return
  #endif Icallfromoverview

  if is3d:
    if NXtick3d >= 0:
      _set_number_of_ticks_3d()
  else:
    if NXtick >= 0:
      _set_number_of_ticks()
  #endif hasattr(Ax,'zaxis')

  if LogX: Ax.set_xscale('log')
  else: Ax.set_xscale('linear')
  if LogY: Ax.set_yscale('log')
  else: Ax.set_yscale('linear')

  Fig.set_visible(visible)

  Ax.set_visible(visible)
  Ax.set_aspect(Aspect)

  if type(Klegend) == int or type(Klegend) == bool:
    if Klegend and len(Legend): plt.legend(Legend)
  else:
    if Klegend.get() and len(Legend): plt.legend(Legend)
    plt.show(block=False)
  #endif type(Klegend) == int or type(Klegend) == bool

  nax = 0
  for ax in Axes:

    if hasattr(ax,'get_label'): nax += 1

    if hasattr(ax,'xaxis'):
      try: ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3),
                               useOffset=True, useMathText=True)
      except: pass
    #endif hasattr(ax,'xaxis')

    if hasattr(ax,'yaxis'):
      try: ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useMathText=True)
      except: pass
    #endif hasattr(ax,'yaxis')

    if hasattr(ax,'zaxis'):
      ax.ticklabel_format(style='sci',scilimits=(-LexpPow,LexpPow))
    #endif hasattr(ax,'zaxis')

  #endfor ax in Axes

  if Itight:
    plt.tight_layout(pad=Tightpad, w_pad=Xtightpad, h_pad=Ytightpad)
  else:
    plt.subplots_adjust(left=LeftMargin, right=RightMargin, bottom=BottomMargin, top=TopMargin, wspace=Xspace, hspace=Yspace)
  #endif

  if not is3d:
    xmin,xmax = Ax.get_xlim()
    if (xmax+xmin) > 0 and (xmax-xmin)/(xmax+xmin) < 1.e-3: NxBinMax = 5
    if NxBinMax: plt.locator_params(axis='x', nbins=NxBinMax)
  #endif not is3d

  plt.grid(Kgrid)

  Figman =  plt.get_current_fig_manager()
  if Igetconsole: get_console()
#  getplotsize()

#enddef showplot()

def optconsole(con=True): global Igetconsole; Igetconsole = con

def hplot2d(idh, plopt='3d', block=False, scalex=1., scaley=1., scalez=1.,
            cmap='', surfcolor='', tit='', xtit='', ytit='', ztit=''):

  import numpy as np

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idh) != int:
    idx = GetIndexH2(idh)
  else:
    if idh < 0: idh = H2I
    idx = idh
  #endif

  if idx == -1:
    print('*** Error in hplot2d: Non-existing histogram ',idh)
    return -1
  elif H2hh[38]:
    print('*** Warning in hplot2d: Histogram has been marked as deleted ',idh)
  #endif idx == -1:

  if Kecho:
    s = "hplot2d('" + H2head[idx][0] + "', plopt='" + plopt + "', block=" + str(block) + \
    ", scalex=" + str(scalex) + ", scaley=" + str(scaley) + ", scalez=" + str(scalez) \
    + ", cmap='" + cmap + "', surfcolor='" +  surfcolor + ", tit='" + \
    tit + "', xtit='" + xtit + "', ytit='" + ytit + "', ztit='" + ztit +"')"
    print(s)
  #endif

  nx = H2head[idx][2]
  ny = H2head[idx][6]

  x = H2[idx]['x'] * scalex
  y = H2[idx]['y'] * scaley

  z = H2[idx]['z'] * scalez
  ez = H2[idx]['ez'] * scalez

  xmin = x.min(); xmax = x.max()
  if xmin == xmax:
    xmin-=abs(xmin*0.1)
    xmax+=abs(xmax*0.1)
  #endif

  ymin = y.min(); ymax = y.max()
  if ymin == ymax:
    ymin-=abs(ymin*0.1)
    ymax+=abs(ymax*0.1)
  #endif

  if nx == 1:
    dx = H2head[idx][5] * np.ones_like(x) * scalex
  else:
    dx = H2head[idx][5] * np.ones_like(x) * scalex
  #endif

  if ny == 1:
    dy = H2head[idx][9] * np.ones_like(y) * scaley
  else:
    dy = H2head[idx][9]  * np.ones_like(y) * scaley
  #endif

  dz = z
  set_plot_params_3d()

  if not plopt: plopt = Mode3d
  #print("*** plopt:",plopt)
  plotoptions(plopt)

  if Ihist:

    if surfcolor == '': surfcolor = Surfcolor

    getzone(projection='3d')

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)
    plt.rcParams['axes.labelsize'] = Atitfontsize3d

    if not Ierr:

      x = x - dx / 2.
      y = y - dy / 2.
      z = z - dz

      Ax.bar3d(x,y,z,dx,dy,dz, color=surfcolor)

    else:

      dx = 0.1 * dx
      dy = 0.1 * dy
      dz = ez

      x = x - dx / 2.
      y = y - dy / 2.
      z = z - dz

      Ax.bar3d(x,y,z,dx,dy,dz, color=surfcolor)

    #endif

    txyz(tit,xtit,ytit,ztit)

  elif Isurf:

    getzone(projection='3d')

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)
    plt.rcParams['axes.labelsize'] = Atitfontsize3d

    xsh = x.values.reshape(nx,ny)
    ysh = y.values.reshape(nx,ny)
    zsh = z.values.reshape(nx,ny)

    if not cmap and not surfcolor:
      surfcolor = None
      if Icmap.get():
        if cmap == '' or cmap == '!': cmap=Cmap
        surfcolor = None
      else:
        cmap = None
        if surfcolor == '': surfcolor = Surfcolor
      #endif
    #endif

    if surfcolor: cmap = None

    Ax.plot_surface(xsh,ysh,zsh,rstride=1, cstride=1, shade=True, cmap=cmap,color=surfcolor)

    txyz(tit,xtit,ytit,ztit)

  elif Icont3d:

    getzone(projection='3d')

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)
    plt.rcParams['axes.labelsize'] = Atitfontsize3d

    xsh = x.values.reshape(nx,ny)
    ysh = y.values.reshape(nx,ny)
    zsh = z.values.reshape(nx,ny)

    if cmap == '' or cmap == '!': cmap=Cmap

    Ax.contourf(xsh,ysh,zsh,100,cmap=cmap)

    txyz(tit,xtit,ytit,ztit)

  elif Itrisurf:

    getzone(projection='3d')

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)
    plt.rcParams['axes.labelsize'] = Atitfontsize3d

    if not cmap and not surfcolor:
      if Icmap.get():
        if cmap == '' or cmap == '!': cmap=Cmap
        surfcolor = None
      else:
        cmap = None
        if surfcolor == '': surfcolor = Surfcolor
      #endif
    #endif

    Ax.plot_trisurf(x,y,z,shade=True,cmap=cmap,color=surfcolor)
    #Ax.plot_surface(xsh,ysh,zsh,rstride=1, cstride=1, cmap=cm.viridis)

    txyz(tit,xtit,ytit,ztit)

  elif Iprof:

    Quit("Ende in hplot2d")

  elif Iinter == 0 :

    getzone()

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist)
    plt.rcParams['axes.labelsize'] = Atitfontsize3d

    if cmap == '' or cmap == '!': cmap=Cmap

    zsh = z.values.reshape(nx,ny)
    zsh = zsh.T

    txyz(tit,xtit,ytit)

    colmap = Ax.imshow(zsh
                       ,cmap=cmap, aspect="auto"
                       ,extent=[xmin,xmax,ymin,ymax] # plot size
             )

    if Colorbarpad != '!':
      fcm = Fig.colorbar(colmap, pad=Colorbarpad)
    else:
      fcm = Fig.colorbar(colmap)
    #endif

    fcm.ax.tick_params(labelsize=Axislabelsize)

    Axes.append(fcm)

  else:

    getzone()

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist)
    plt.rcParams['axes.labelsize'] = Atitfontsize

    methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
               'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
               'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

    zsh = z.values.reshape(nx,ny)
    zsh = zsh.T

    if cmap == '' or cmap == '!': cmap=Cmap

#    cmap = plt.cm.copper
#    ls = LightSource(315, 45)
#    zsh = ls.shade(zsh,plt.cm.copper)

    txyz(tit,xtit,ytit)

    colmap = Ax.imshow(zsh,
                       interpolation='spline36'
                       ,cmap=cmap, aspect="auto"
                       ,extent=[xmin,xmax,ymin,ymax] # plot size
             )

    if Colorbarpad != '!':
      fcm = Fig.colorbar(colmap, pad=Colorbarpad)
    else:
      fcm = Fig.colorbar(colmap)
    #endif

    fcm.set_label(label=ztit, labelpad=Axistitledist, size=Axislabelsize)
    fcm.ax.tick_params(labelsize=Axislabelsize)
    fcm.ax.get_yaxis().get_offset_text().set(size=Axislabelsize,
                                             position=(1+0.5*Nxzone,0.))

    Axes.append(fcm)

  #endif Ihist:

  plt.show(block=block)

  if Kdump:
    Ndump += 1
    fout = WaveFilePrefix + str(Ndump) + ".dat"
    vwritexyz(x,y,z,fout)
    print("\nData written to ",fout)
  #endif Kdump

  if Kpdf:
    Npdf += 1
    fout = WaveFilePrefix + str(Npdf) + ".pdf"
    pplot(fout)
    #print("\nPlot written to ",fout)
  #endif Kdump

  Kplots[Kzone-1] = 1

#enddef hplot2d(idh,plopt='')

def nextcolor():
  global Kecho
  if Kecho:
    print("nextcolor()")
  nextlinecolor()
  nextmarkercolor()
  nexthistcolor()
#enddef nextcolor()

def nexthistcolor():
  global Colors
  col = gethistcolor()
  icol = getcolorindex(col)
  if icol > len(Colors) - 3: icol = 1
  else: icol += 1
  sethistcolor(Colors[icol])
  if Kecho:
    print("nexthistcolor:",Colors[icol])
#enddef nexthistcolor()

def nextmarkercolor():
  global Colors
  col = getmarkercolor()
  icol = getcolorindex(col)
  if icol > len(Colors) - 3: icol = 1
  else: icol += 1
  setmarkercolor(Colors[icol])
  if Kecho:
    print("nextmarkercolor:",Colors[icol])
#enddef nextmarkercolor()

def nextlinecolor():
  global Colors
  col = getlinecolor()
  icol = getcolorindex(col)
  if icol > len(Colors) - 3: icol = 1
  else: icol += 1
  setlinecolor(Colors[icol])
  if Kecho:
    print("nextlinecolor:",Colors[icol])
#enddef nextlinecolor()

def samezone(isame=1):
  global Isame, Kecho
  if Kecho:
    print("samezone(isame=" + str(isame) +")")
  IsameGlobal.set(isame)
  if not isame: Isame = 0
#enddef samezone(isame=1)

def nextzone(projection='2d', visible=True):

  global Nxzone, Nyzone, Kzone, Isame, Kecho, Kplots, Kzone

  if Kecho:
    print("nextzone(projection=" + projection + ", visible=" + str(visible) +")")

  try:
    IsameGlobal.set(0)
  except:
    return
  #endtry

  Isame = 0
  Kplots[Kzone-1] = 0
  reset_zoom()

  nx = Nxzone
  ny = Nyzone
  kzone = Kzone

  if kzone < nx * ny:
    kzone += 1
    zone(nx,ny,kzone,'s',projection=projection,visible=visible)
  else:
    kzone = 1
    zone(nx,ny,kzone,'',projection=projection,visible=visible)
  #endif kzone == nx * ny

#def nextzone(projection='2d', visible=True)

def zone(nx=1, ny=1, kzone=1, isame='', projection='2d', visible=True):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug


  global Tax2d, Tax3d, Debug

  if Nwins <= 0:

    Nfigs = len(plt.get_fignums())

    if Nfigs > Nwins:

      Nwins = Nfigs

      Fig = plt.gcf()

      if not Axes:
        Kzone = 1
        Nyzone = 1
        Nxzone = 1
      #endif

      Ax = plt.gca()
      Axes = Fig.axes

      z = []
      z.append(Fig)
      z.append(Nxzone)
      z.append(Nyzone)
      z.append(Kzone)
      z.append(1) # date on figure
      Zones.append(z)

#      return

    else:

      window(projection=projection)
      reset_zoom()

    #endif Nfigs > Nwins

  #endif Nwins <= 0

  if kzone > nx * ny:
    print("*** Error in zone: kzone > nxzone * nyzone ***")
    return
  #endif kzone > nx * ny

  if Nwins <= 0:
    Nfigs = len(plt.get_fignums())
    if Nfigs > Nwins: Nwins = Nfigs
    window(projection=projection)
    reset_zoom()
  #endif Nwins <= 0

  nxyk = ny*100+nx*10+kzone
  label = str(nx) + "_" + str(ny) + "_" + str(kzone)

  Legend = []

  if type(IsameGlobal) == int:
    ksame = IsameGlobal
  else:
    ksame = IsameGlobal.get()
  #endif type(IsameGlobal) == int

  if not isame and Isame == 0 and ksame == 0:

    Fig.clear()
    Tdate = None
    reset_zoom()
    date_on_figure()

    #if Waveplot: run_on_figure(iforce=1)

    for z in Zones:
      if z[0] == Fig:
        z[4] = 0
        break
    #endfor z in Zones:

    if projection.lower() == '2d':
      Ax = Fig.add_subplot(ny,nx,kzone,visible=visible,label=label)
      Tax2d = type(Ax)
    elif projection.lower() == '3d':
      Ax = Fig.add_subplot(ny,nx,kzone,projection='3d',visible=visible,label=label)
      Tax3d = type(Ax)
      set_plot_params_3d()
    else:
      Ax = Fig.add_subplot(ny,nx,kzone,projection=projection,visible=visible,label=label)
    #endif

  else:

    if nx != Nxzone or ny != Nyzone or kzone != Kzone:
      if projection.lower() == '2d':
        Ax = Fig.add_subplot(ny,nx,kzone,visible=visible,label=label)
        Tax2d = type(Ax)
      elif projection.lower() == '3d':
        Ax = Fig.add_subplot(ny,nx,kzone,projection='3d',visible=visible,label=label)
        Tax3d = type(Ax)
        set_plot_params_3d()
      else:
        Ax = Fig.add_subplot(ny,nx,kzone,projection=projection,visible=visible,label=label)
      #if projection.lower() == '2d'
    #endif nx == Nxzone and ny == Nyzone and kzone == Kzone
  #endif isame == '' and Isame == 0 and IsameGlobal.get() == 0

  Nxzone = nx
  Nyzone = ny
  Kzone = kzone

  ifound = 0
  for z in Zones:
    if z[0] == Fig:
      z[1] = Nxzone
      z[2] = Nyzone
      z[3] = Kzone
      ifound = 1
      break
  #endfor z in Zones:

  if ifound == 0:
    z = []
    z.append(Fig)
    z.append(Nxzone)
    z.append(Nyzone)
    z.append(Kzone)
    Zones.append(z)
  #endif ifound == 0

  Ax.set_visible(visible)

  #Debug += 1
  if hasattr(Ax,'zaxis'):
    set_plot_params_3d()
  else:
    set_plot_params()
  #endif hasattr(Ax,'zaxis')

#enddef zone(nx,ny)

def window_close(win=-1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  if Nwins == 1:
    Nwins = 0
    Nfigs = 0
    Figs = []
    Legend = []
    Zones = []
    Nxzone = 0
    Nyzone = 0
    Kzone = 0
    Kplot = []
    Fig = None
    Ax = None
    plt.close(plt.gcf())
    return
  #endif

  if win == -1: win = plt.gcf()

  figs = []
  Legend = []

  for w in Figs:
    if w != win: figs.append(w)
  #endfor w in Figs

  plt.close(win)

  Figs = figs
  Nwins = len(Figs)
  Nfigs = Nwins

  Fig = None
  Ax = None

  if Nwins:
    Fig = plt.gcf()
    Ax = plt.gca()
  #endif

#enddef window_close(win=None):

def window_clear(win=-1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  global Tfig, Tdate, Figman

  if win == -1:
    Fig = plt.gcf()
    TFig = type(Fig)
  else:
    Fig = win
  #endif win == -1

  Fig.clear()
  Tdate = None

  zone(1,1)
  Figman =  plt.get_current_fig_manager()

#enddef window_clear(win=-1)

def figure_clear(kfig=-1):
  global Fig, Tfig,Tdate,Figman
  if kfig == -1:
    Fig = plt.gcf()
    Fig.clear()
    TFig = type(Fig)
  else:
    kfig.clf()
  #endif kfig == -1
  Tdate = None
  Figman =  plt.get_current_fig_manager()
#enddef figure_clear(kfig=-1)

def figure_close(kfig=-1):
  if kfig == -1:
    plt.close()
  else:
    plt.close(kfig)
  #endif kfig == -1
#enddef figure_close(kfig=-1)

def set_title(title='Title',tfs=-9.,titx=-9.,tity=-9):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if tfs == -9.:
    tfs = TitFontSize
    if Nxzone*Nyzone > 1:
      if type(TextFontSize) == str: tfs = 'small'
      else: tfs = TextFontSize*0.75
  if titx == -9.: titx = Xtitle
  if tity == -9.: tity = Ytitle

  Ptit = title
  Ax.set_title(title,fontsize=tfs,x=titx,y=tity)

  showplot()
#def set_title(title='Title')

def set_x_title(xtit='xTit',pos=0.5):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  Xtit = xtit
  Ax.set_xlabel(xtit,x=pos)
  shpl()
#def set_x_title(title='xTit')
def get_x_title(): return Xtit

def set_y_title(ytit='yTit', pos=0.5):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  Ytit = ytit
  Ax.set_ylabel(ytit,y=pos)
  shpl()
#def set_y_title(title='yTit')
def get_y_title(): return Ytit

def set_z_title(ztit='zTit',pos=0.5):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  if hasattr(Ax,'zaxis'):
    Ztit = ztit
    Ax.set_zlabel(ztit,z=pos)
    shpl()
#def set_z_title(title='zTit')
def get_z_title(): return Ztit

def set_titles(gtit='',pltit='Title',xtit='xTit', ytit='yTit', ztit=''):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if gtit != '':
    Fig.suptitle(gtit, fontsize=GtitFontSize)

  txyz(pltit,xtit,ytit,ztit)

  plt.show(block=False)
#enddef set_titles(pltit='Title',xtit='xTit', ytit='yTit', ztit='')

def set_global_title(gtit='', fontsize='!'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global Figman

  if fontsize == '!': fontsize=GtitFontSize

  Fig = plt.gcf()
  Tfig = type(Fig)
  Figman =  plt.get_current_fig_manager()

  Gtit = gtit

  if gtit != '':
    Fig.suptitle(gtit, fontsize=fontsize)

  plt.show(block=False)
#def set_global_title(gtit='Run_and_Code', fontsize=16)

def txyz(pltit='Title',xtit='', ytit='', ztit='', tfs=-9., xyzfs=-9,
        titx=-9, tity=-9, titlesize=-9):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  global Ktitles
  if not Ktitles: return

  if not Ax: return

  if tfs == -9.:
    tfs = TitFontSize
    if Nxzone*Nyzone > 1:
      if type(TextFontSize) == str: tfs = 'small'
      else: tfs = TextFontSize*0.75
    #endif Nxzone*Nyzone > 1
  #endif tfs == -9.

  if tity == -9: tity = Ytitle
  if titx == -9: titx = Xtitle

  if xyzfs == -9:
    xyzfs = Atitfontsize
    if Nxzone*Nyzone > 1: xyzfs *= 0.75
  #endif xyzfs == -9

  if titlesize == -9:
    Titlesize = Atitfontsize
    if Nxzone*Nyzone > 1: Titlesize *= 0.75
  #endif titlesize == -9

  Xtit = xtit
  Ytit = ytit
  Ztit = ztit
  Ptit = pltit

  #plt.rc('font',size=xyzfs) # for the exponent

  #if pltit != '': Ax.set_title(pltit,fontsize=tfs,position=(titx,tity),pad=TitPad)

  #if pltit != '': Ax.set_title(pltit,x=titx,y=tity)
  #if xtit != '': Ax.set_xlabel(xtit)
  #if ytit != '': Ax.set_ylabel(ytit)

  if hasattr(Ax,'zaxis') and ztit != '':
    plt.rcParams['axes.labelpad'] = Axistitledist3d
    plt.rcParams['axes.labelsize'] = Atitfontsize3d
  else:
    plt.rcParams['axes.labelpad'] = Axistitledist
    plt.rcParams['axes.labelsize'] = Atitfontsize
  #endif hasattr(Ax,'zaxis') and ztit != ''

  if pltit != '': Ax.set_title(pltit,fontsize=tfs,x=titx,y=tity)

  if xtit != '': Ax.set_xlabel(xtit,fontsize=xyzfs)
  if ytit != '': Ax.set_ylabel(ytit,fontsize=xyzfs)

  txexp = Ax.xaxis.get_offset_text()
  txexp.set_size(Axislabelsize)
  tyexp = Ax.yaxis.get_offset_text() # for the exponent
  tyexp.set_size(Axislabelsize)

  fcm = Axes[len(Axes)-1]

  if hasattr(Ax,'zaxis') and ztit != '':

    tzexp = Ax.zaxis.get_offset_text()

    # Bug in matplotlib:
    #tzexp.set_size(AxisLabelSize)
    #tzexp.set_position(xy=(LexpX,LexpY))
    #tzexp.set_rotation(LexpRot)

    # Workaround: See showplot()

    tzexp.set(visible=False)
    zticks = Ax.get_zticks()
    zamax = max(abs(zticks[:-1]))

    zexp = np.log10(zamax)

    if zexp >= 0: izexp = int(zexp)
    else: izexp = int(zexp) - 1

    if abs(izexp) >= 2:
      ztit += '       1.e' + str(izexp)
      #ztit += '                      '
    #endif abs(izexp) >= 2

    Ax.set_zlabel(ztit,fontsize=xyzfs,linespacing=0.5) #,fontname='Helvetica',family='serif')
    #Ax.set_zlabel(ztit,fontsize=xyzfs) #,fontname='Helvetica',family='serif')

  elif type(fcm) == Tcmap:

    fcm.set_label(label=ztit, labelpad=Axistitledist, size=Axislabelsize)
    fcm.ax.tick_params(labelsize=Axislabelsize)
    fcm.ax.get_yaxis().get_offset_text().set(size=Axislabelsize,
                                             position=(1+0.5*Nxzone,0.))
  #endif hasattr(Ax,'zaxis') and ztit != ''

  plt.show(block=False)

#enddef set_txyz(pltit='Title',xtit='xTit', ytit='yTit', ztit='')

def null3d(xmin=-10., xmax=10., ymin=-10., ymax=10., zmin=-10., zmax=10.):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  global Tfig, Tax3d

  if Nwins <= 0: window()

  getzone('3d')

  Ax.set_xlim3d(xmin, xmax)
  Ax.set_ylim3d(ymin, ymax)
  Ax.set_zlim3d(zmin, zmax)

  Tax3d = type(Ax)

  showplot()
#enddef null3d(xmin=-10., xmax=10., ymin=-10., ymax=10., zmin=-10., zmax=10.)

def null(xmin=-10., xmax=10., ymin=-10., ymax=10.):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  global Tfig,Tax2d
  #nreakpoint()

  if Nwins <= 0: window()
  else: getzone()

  plt.axis([xmin,xmax,ymin,ymax])
  Kplots[Kzone-1] = 1
  Tax2d = type(Ax)
  showplot()
#enddef null(xmin=-10., xmax=10., ymin=-10., ymax=10.)

def nnull(xmin=-10., xmax=10., ymin=-10., ymax=10.):
  nextzone()
  null(xmin, xmax, ymin, ymax)
#enddef

def figtext(x,y,text,fontsize='!', ishow=1):
  global Tfig, TextFontSize
  # text for normalized axes coordinates
  if fontsize == '!': fontsize=TextFontSize
  Fig = plt.gcf()
  Tfig = type(Fig)
  Fig.text(x,y,text,fontsize=fontsize)
  if ishow: showplot()
#enddef figtext(x,y,text,fontsize=12, ishow=1)

def textndc(x,y,text,fontsize='!', ishow=1, color='black',halign='center',
            valign='center',angle=0.):

  global TextFontSize, Textcolor

  # text for normalized figure coordinates
  # Compare to textNDC()

  if fontsize == '!': fontsize=TextFontSize

  ax = getax()

  col = color
  if color == -9: col = Textcolor

  Ax.text(x,y,text,fontsize=fontsize, color=col,
          transform=Fig.transFigure,
          horizontalalignment=halign,verticalalignment=valign,rotation=angle)

  if ishow: showplot()

#enddef textndc(x,y,text,fontsize=12, ishow=1)

def text3d(x,y,z,text,fontsize='!', ishow=1, fontname='', fontfamily='',color='black',
         halign='center', valign='center',angle=0.):
  global TextFontSize, Textcolor, Tax2d, Tax3d

  # text for normalized figure coordinates

  col = color
  if color == -9: col = Textcolor

  if fontsize == '!':
    if Nxzone*Nyzone == 1:
      fontsize=TextFontSize
    else:
      if type(TextFontSize) == str: fontsize = 'small'
      else: fontsize=TextFontSize*0.75
  #endif fontsize == '!'

  ax = getax()

  if type(Ax) == Tax2d:
    text(x,y,text,fontsize, ishow, fontname, fontfamily,color,halign, valign,angle)
    return
  #endif

  if fontname and fontfamily:
    Ax.text(x,y,z,text,transform=Ax.transAxes,fontsize=fontsize,
            fontname=fontname, fontfamily=fontfamily, color=col,
            horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  elif fontname:
    Ax.text(x,y,z,text,transform=Ax.transAxes,fontsize=fontsize,
            fontname=fontname,color=col,
            horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  elif fontfamily:
    Ax.text(x,y,z,text,transform=Ax.transAxes,fontsize=fontsize,
            fontfamily=fontfamily,color=col,
            horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  else:
      Ax.text(x,y,z,text,transform=Ax.transAxes,fontsize=fontsize,color=col,\
      horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  #endif
#endif fontname and fontfamily

  if ishow: showplot()

#enddef text(x,y,text,fontsize=12, ishow=1)

def text(x,y,text,fontsize='!', ishow=1, fontname='', fontfamily='',color='black',
         halign='center', valign='center',angle=0.):
  global TextFontSize, Textcolor, Tax2d, Tax3d

  # text for normalized figure coordinates

  if type(Ax) == Tax3d:
    print('\n --- for 3d use text3d(...) ---')
    return
  #endif

  col = color
  if color == -9: col = Textcolor

  if fontsize == '!':
    if Nxzone*Nyzone == 1:
      fontsize=TextFontSize
    else:
      if type(TextFontSize) == str: fontsize = 'small'
      else: fontsize=TextFontSize*0.75
  #endif fontsize == '!'

  ax = getax()

  if fontname and fontfamily:
    Ax.text(x,y,text,transform=Ax.transAxes,fontsize=fontsize,
            fontname=fontname, fontfamily=fontfamily, color=col,
            horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  elif fontname:
    Ax.text(x,y,text,transform=Ax.transAxes,fontsize=fontsize,
            fontname=fontname,color=col,
            horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  elif fontfamily:
    Ax.text(x,y,text,transform=Ax.transAxes,fontsize=fontsize,
            fontfamily=fontfamily,color=col,
            horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  else:
    Ax.text(x,y,text,transform=Ax.transAxes,fontsize=fontsize,color=col,\
    horizontalalignment=halign,verticalalignment=valign,rotation=angle)
  #endif
#endif fontname and fontfamily

  if ishow: showplot()

#enddef text(x,y,text,fontsize=12, ishow=1)

def get_run_on_figure():
  global Krun,Kruns
  return Krun
#endif

def set_run_on_figure(krun=True):
  global Krun,Kruns
  Krun = krun
#endif

def run_on_figure(x=0.03,y=0.95,fontsize='!',ishow=1, iforce=0):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  global Krun,Kruns

  if ROFx: x = ROFx
  if ROFy: y = ROFy

  if Waveplot == 1: #WAVE

    if Krun == None or Krun == False: return

    global Icallfromoverview
    if Icallfromoverview: return

    try:
      if krunini: pass
    except:
      krunini = -1
      krunold = False
    #endtry

    if fontsize == '!': fontsize=DateFontSize

    srun = "Run " + str(Wrun)
    ftext = plt.gcf()

    if not Trun or iforce:
      Trun = ftext.text(x,y,srun,fontsize=fontsize)
    elif Trun:
      Trun.set_x(x)
      Trun.set_y(y)
      Trun.set_text(srun)
      Trun.set_fontsize(fontsize)
    #endif

    if Krun:
      Trun.set_visible(True)
      if TrunOv: TrunOv.set_visible(True)
    else:
      Trun.set_visible(False)
      if TrunOv: TrunOv.set_visible(False)
    #endif Krun

    if (krunini == -1 or Krun != krunold) and ishow == 1: plt.show(block=False)

    krunold = Krun

  elif Waveplot == 2: #UNDUMAG

    if Krun == None or Krun == False: return

    try:
      if krunini: pass
    except:
      krunini = -1
      krunold = False
    #endtry

    if fontsize == '!': fontsize=DateFontSize

    #if len(Mcomment):
    srun = "Run " + str(Mrun) + ": " + str(Mcomment) + "\n" + str(Mdate)
    #else:
    #  srun = "Run " + str(Mrun) + ": " + "\n" + str(Mdate)

    ftext = plt.gcf()

    if not Trun or iforce:
      #print("ROF",x)
      #x = len(srun)
      Trun = ftext.text(x,y,srun,fontsize=fontsize)
    #endif

    if Krun:
      Trun.set_visible(True)
      if TrunOv: TrunOv.set_visible(True)
    else:
      Trun.set_visible(False)
      if TrunOv: TrunOv.set_visible(False)
    #endif Krun

    if (krunini == -1 or Krun != krunold) and ishow == 1: plt.show(block=False)

    krunold = Krun

  #endif Waveplot == 1

#enddef run_on_figure():

def date_on_figure(x=0.04,y=0.02,fontsize='!',ishow=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  global Kdate

  if Kdate == None or Kdate == False: return

  try:
    if Tdate:
      Tdate.remove()
      Tdate = None
  except:
    Tdate = None
  #endtry

  try:
    if kdateini: pass
  except:
    kdateini = -1
    kdateold = False
  #endtry

  if fontsize == '!': fontsize=DateFontSize

  localtime = time.asctime( time.localtime(time.time()) )
  if Author: localtime = Author + ', ' + localtime

  if Icallfromoverview:

    if Krun:
      if TrunOv: TrunOv.set_visible(True)
    else:
      if TrunOv: TrunOv.set_visible(False)
    if Kdate:
      if TdateOv: TdateOv.set_visible(True)
    else:
      if TdateOv: TdateOv.set_visible(False)
    #endif

  else:

    if not Tdate and Fig != None: Tdate = Fig.text(x,y,localtime,fontsize=fontsize)

    if Krun:
      if Trun: Trun.set_visible(True)
      if TrunOv: TrunOv.set_visible(True)
    else:
      if Trun: Trun.set_visible(False)
      if TrunOv: TrunOv.set_visible(False)
    #endif Krun

    if Kdate:
      if Tdate: Tdate.set_visible(True)
      if TdateOv: TdateOv.set_visible(True)
    else:
      if Tdate: Tdate.set_visible(False)
      if TdateOv: TdateOv.set_visible(False)
    #endif Kdate

  if (kdateini == -1 or Kdate != kdateold) and ishow == 1:
    plt.show(block=False)

  kdateold = Kdate
  #run_on_figure()

#enddef date_on_figure():

def optnrun(krun=False):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  Krun = krun
  run_on_figure()
#enddef optnrun()

def optrun(krun=True):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  Krun = krun
  run_on_figure()
#enddef optrun()

def optndate(kdate=False):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  Kdate = kdate
  date_on_figure()
#enddef optndate()

def optdate(kdate=True):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  Kdate = kdate
  date_on_figure()
#enddef optdate()

def getkdate():
  global Kdate
  return Kdate
#enddef

def getkrun():
  global Krun,Kruns
  return Krun
#enddef

def set_author(author=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  if author != '': Author = author
#enddef set_author(author='')

def hcopy1d(idh,idnew,tit='',scalex=1.,scaley=1., reset=0, overwrite=True):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idnew) != str: idh = 'h' + str(idnew)

  idx = GetIndexH1(idh)

  if idh == -1:
    print('*** Error in hcopy1d: Non-existing source histogram ',idh,' ***')
    return -1
  #endif idh == -1

  if tit == '': tit = H1hh[1]

  head1 = deepcopy(H1hh)

  nx = int(H1hh[2])
  xmin = H1hh[3] * scalex
  xmax = H1hh[4] * scalex

  ny = int(H1hh[6])
  ymin = H1hh[7] * scaley
  ymax = H1hh[8] * scaley

  hret = hbook1(idnew,tit,nx,xmin,xmax,overwrite=overwrite)

  idxnew = GetIndexH1(idnew)
  if reset: return idxnew

  head1[0] = idnew
  head1[1] = tit

  head1[3] *= scalex
  head1[4] *= scalex
  head1[5] *= scalex
  head1[11] *= scalex**2
  head1[12] *= scalex
  head1[13] *= scalex
  head1[18] *= scalex
  head1[19] *= scalex**2

  head1[7] *= scaley
  head1[8] *= scaley
  head1[10] *= scaley

  H1[idxnew] = deepcopy(H1[idx])
  H1[idxnew].y *= scaley
  H1[idxnew].ey *= scaley
  H1[idxnew].y2 *= scaley**2

  H1head[idxnew] = head1

  return idxnew
#def hcopy1d(idh,idnew,tit='',scalex=1.,scaley=1., reset=0)

def hcopy2d(idh,idnew,tit='',scalex=1.,scaley=1., scalez=1., reset=0, overwrite=True):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  idx = GetIndexH2(idh)

  if idh == -1:
    print('*** Error in hcopy2: Non-existing source histogram ',idh,' ***')
    return -1
  #endif idh == -1

  if tit == '': tit = H2hh[1]

  head2 = deepcopy(H2hh)

  nx = int(H2hh[2])
  xmin = H2hh[3] * scalex
  xmax = H2hh[4] * scalex

  ny = int(H2hh[6])
  ymin = H2hh[7] * scaley
  ymax = H2hh[8] * scaley

  hret = hbook2(idnew,tit,nx,xmin,xmax,ny,ymin,ymax,overwrite=overwrite)

  idxnew = GetIndexH2(idnew)
  if reset: return idxnew

  H2[idxnew] = deepcopy(H2[idx])

  H2[idxnew].z = H2[idx].z * scalez
  H2[idxnew].z2 = H2[idx].z2 * scalez
  H2[idxnew].ave = H2[idx].ave * scalez
  H2[idxnew].ez = H2[idx].ez * scalez

  H2head[idxnew][14] = head2[14] * scalex
  H2head[idxnew][15] = head2[15] * scalex
  H2head[idxnew][16] = head2[16] * scaley
  H2head[idxnew][17] = head2[17] * scaley
  H2head[idxnew][20] = head2[20] * scalez
  H2head[idxnew][21] = head2[21] * scalez
  H2head[idxnew][24] = head2[24] * scalez
  H2head[idxnew][25] = head2[25] * scalez
  H2head[idxnew][28] = head2[28] * scalez
  H2head[idxnew][29] = head2[29] * scalez
  H2head[idxnew][32] = head2[32] * scalez
  H2head[idxnew][33] = head2[33] * scalez
  H2head[idxnew][34] = head2[34] * scalex
  H2head[idxnew][35] = head2[35] * scalex * scalex
  H2head[idxnew][36] = head2[36] * scaley
  H2head[idxnew][37] = head2[37] * scaley * scaley

#def hcopy2d(idh,idnew,tit='',scalex=1.,scaley=1., scalez=1.)

def nplc(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcr(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcb(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcc(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcm(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplgr(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcs(nt='?',varlis='',select='',weights='',plopt='samespline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcrs(nt='?',varlis='',select='',weights='',plopt='samespline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcbs(nt='?',varlis='',select='',weights='',plopt='samespline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplccs(nt='?',varlis='',select='',weights='',plopt='samespline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcms(nt='?',varlis='',select='',weights='',plopt='samespline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplcgs(nt='?',varlis='',select='',weights='',plopt='samespline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nnplc(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nnplcr(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nnplcb(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nnplcc(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nnplcm(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def nplgr(nt='?',varlis='',select='',weights='',plopt='spline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef

def npll(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef npll(

def nnpll(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef npll(

def nplm(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nplm(

def nnplm(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nnplm(

def npllr(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef npllr(

def nnpllr(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nnpllr(

def nplmr(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nplmr(

def nnplmr(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nnplmr(

def npllb(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nplot(nt,varlis,select,weights,plopt, legend, \
  scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllr(

def nplmb(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nplmb(

def nnpllb(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
  scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nnpllr(

def nnplmb(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nnplmb(

def npllc(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef npllc(

def nplmc(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nplmc(

def npllc(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef npllc(

def nplmy(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='yellow'):
  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nplmy(

def nnplly(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='yellow'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nnplly(

def nnplmc(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):
  nnplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)
#endef nnplmc(

def npllm(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllm(

def nplmm(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmm(

def npllg(nt='?',varlis='',select='',weights='',plopt='line', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllg(

def nplmg(nt='?',varlis='',select='',weights='',plopt='marker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmg(

def nplls(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplls(

def nplms(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplms(

def npllrs(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllrs(

def nplmrs(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='red'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmrs(

def npllbs(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllrs(

def nplmbs(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='blue'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmbs(

def npllgs(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllgs(

def npllms(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllms(

def nplmms(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmms(

def npllcs(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllcs(

def nplmcs(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='cyan'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmcs(

def npllls(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllls

def nplmls(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='magenta'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmls

def npllys(nt='?',varlis='',select='',weights='',plopt='sameline', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='yellow'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef npllys(

def nplmys(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='yellow'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmys(

def nplmgs(nt='?',varlis='',select='',weights='',plopt='samemarker', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='green'):

  nplot(nt,varlis,select,weights,plopt, legend, \
        scalex, scaley, scalez, scalet, cmap, hist, color)

#endef nplmgs(

def nplot(nt='?',varlis='',select='',weights='',plopt='', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='!',
          color='default',isort=0,nx=101,ny=101):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  #reakpoint()

  NxBinMax = 0
  nto = nt

  kstato = Kstat

  global FillColor

  #LastPlot = ['nplot',varlis,select,weights,plopt,legend,scalex,scaley,scalez,scalet,cmap]
  if type(nt) == str and nt == '?':
    print("\nnplot(nt='', varlis='', select='', weights='', plopt='', legend='', scalex=1., scaley=1., scalez=1., scalet=1., hist='Hnplot', color='default', isort=0)")
    print("\nFor weights != '', the markers are colored accoring to the colormap and a color bar is drawn.")
    print("For a single variable and plopt == '' or '!' the data are projected into histgram first.")
    print("For two variables and plopt == 'prof' the data are projected into a profile histogram first.")
    return 0
  #endif type(nt) == str and len(nt) == 0

  if type(weights) != str:
    print("*** Error in nplot: Argument 'weights' must be a string ***")
    return -1
  #endif type(weights) != str

  idn = GetIndexN(nt)

  if idn == -1:
    print("*** Non-existing Ntuple in nplot ***")
    return -1
  #endif idn == -1

  nt = N
  ntname = Nhead[idn][1]
  nttitle = Nhead[idn][2]
  nvar = Nhead[idn][3]

  if color == 'default':
    mcol = Markercolor
    lcol = Linecolor
  else:
    mcol = color
    lcol = color
  #endif color == 'default'

  if hist == '!':
    if hexist("HnPlot"):
      hdelete("HnPlot")
    #endif
    hist = "HnPlot"
  #endif

  idx = GetIndexH1(hist)

  if Kecho:
    if type(hist) == str: shis = "'" + hist + "'"
    elif type(hist) == int: shis = str(hist)
    else: shis = H1hh[1]
    if idn != -1:
      s = "nplot('" + Nhead[idn][1] + "', '" + varlis + "', select='" + select + "', weights='" + weights + "', plopt='" + plopt + "', legend='" + \
      legend + "', scalex=" + str(scalex) + ", scaley=" + str(scaley) + ", scalez=" + str(scalez) + ", scalet=" + str(scalet) + ", cmap='" + cmap + "'," + shis + ")"
    else:
      s = "nplot(nt, '" + varlis + "', select='" + select + "', weights='" + weights + "', plopt='" + plopt + "', legend='" + \
      legend + "', scalex=" + str(scalex) + ", scaley=" + str(scaley) + ", scalez=" + str(scalez) + ", scalet=" + str(scalet) + ", cmap='" + cmap + "'," + shis + ")"
    print(s)
  #endif Kecho

  if Nwins <= 0: window()

  plotoptions(plopt)

  varliso =varlis

  if type(varlis) == str:
    varlis = nlistcolon(varlis)
    if varlis[0] == '':
      varlis = nt.columns
      vdum = []
      vdum.append(varlis[0])
      if len(varlis) > 1:
        vdum.append(varlis[1])
      varlis = vdum
    #endif varlis[0] == ''
  #endif type(varlis) == str

  if select == '!': select = ''

  if select:
    nt = nt.query(select)
    Nsel = nt
  #endif select:

  if len(varlis) == 1:

    getzone()
    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist)

    sx = "(" +nparse(nt,varlis[0]) + ") * " + str(scalex)

    iplot = 0

    if Ihist or plopt == '' or plopt == '!': Ihist = 1

    if plopt == 'same' or plopt == 'S':
      Ihist = 1
      plopt = 'errhistsame'
    #endif

    if Ihist or Iprof or Ierr or Imarker:

      idx = GetIndexH1(hist)
      nproj1(nt,varlis[0],weights,select,scalex=scalex,nx=nx,idh=hist)
      if idx > -1: hist = H1[idx]
      ocol=getmarkercolor()
      hcolo=gethistcolor()
      setmarkercolor(mcol)
      if color != 'default': sethistcolor(color)

      hplot1d(hist,plopt,legend=legend)
      setmarkercolor(ocol)
      sethistcolor(hcolo)

      return
    #endif Ihist:

    if not iplot:
        sopt = ", c='" + mcol + "',marker='" + Markertype + "',ms=" + str(Markersize) + ",ls=''"
        if FillColor != 'none': plt.fill(eval(sx),color=FillColor)
        eval('plt.plot(' + sx + sopt + ')')
        #plt.plot(x,c=mcol,marker=Markertype,fillstyle=Fillstyle, ms=Markersize,ls='')
    #endif not plopt:

    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      print("\nData written to ",fout)
    #endif

  elif len(varlis) == 2:

    getzone()

    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist)

    if len(weights) == 0:

      sx = "(" + nparse(nt,varlis[0]) + ") * " + str(scalex)
      sy = "(" + nparse(nt,varlis[1]) + ") * " + str(scaley)

      iplot = 0

      try:

        imark = 0
        isplin = 0

        if Ispline == 1:
          nspline(nt,varliso,select,1001)
          imark = Imarker
          if Isame: vplxy(Nspline.x,Nspline.y,'sameline',color=lcol)
          else: vplxy(Nspline.x,Nspline.y,'line',color=lcol)
          iplot = 1
          Kstat = 0
          if imark:
            iplot = 0
            Imarker = imark
          #if Imarker: iplot = 0

        elif Iline == 1 or Iclosed:

          sopt = ", c='" + lcol + "',ls='" + Linestyle + "',lw=" + str(Linewidth)

          if isort:

            scom = "global VsortX, VsortY; VsortX, VsortY = vsortxy(" + sx + "," + sy + ")"
            exec(scom)
            if Iclosed:
              sx = list(VsortX)
              sy = list(VsortY)
              sx.append(sx[0])
              sy.append(sy[0])
#              VsortX = sx
#              VsortY = sy
            #endif

          else: #isort

            scom = "global VxyzX, VxyzY; VxyzX, VxyzY = vxy(" + sx + "," + sy +  ")"
            exec(scom)
            sx = list(VxyzX)
            sy = list(VxyzY)
            if Iclosed:
              sx.append(sx[0])
              sy.append(sy[0])
            #endif #Iclosde

          #endif isort

          if FillColor != 'none':
            sval = ' plt.fill(sx,sy,color=FillColor)'
          else:
            sval = 'plt.plot(sx,sy' + sopt + ')'
          #endif

          exec(sval)

          iplot = 1
          if Imarker: iplot = 0

        elif Iprof or Ierr:
          idx = GetIndexH1(hist)
          nproj1(nt,varlis[0],varlis[1],select,scalex=scalex,nx=nx,idh=hist)
          if idx > -1: hist = H1[idx]
          ocol=getmarkercolor()
          setmarkercolor(mcol)
          hplot1d(hist,plopt,legend=legend)
          setmarkercolor(ocol)
          return
        elif Ihist or Isurf or Itrisurf or Iinter or Iboxes :
          nproj2(nto,varliso,weights,select,scalex=scalex,scaley=scaley,nx=nx,ny=ny,idh=hist)
          hplot2d(hist,plopt)
          iplot = 1
        #endif Iline == 1

        if not iplot:
          if FillColor != 'none': plt.fill(eval(sx),eval(sy),color=FillColor)
          sopt = ", ls='',  c='" + mcol + "',marker='" + Markertype + "',ms=" + str(Markersize) + ", fillstyle='" + str(Fillstyle) + "'"
          exec('plt.plot(' + sx + ',' + sy + sopt + ')')
        #endif not iplot

        Ispline = isplin

      except:
        print("\n*** Error in nplot(...) while parsing of list of variables, due to wrong variable or due to floats with exponents. \nIn this case try scaling factor instead ***\n")
        return
      #endtry

      if Kdump:
        Ndump += 1
        fout = WaveFilePrefix + str(Ndump) + ".dat"
        #eval("vwritexy(" + sx + "," + sy + ",'" + fout + "')")
        eval("vwritexy(sx,sy,'" + fout + "')")
        print("\nData written to ",fout)
      #endif

      if Kstat:

        try:

          xmin,xmax,xmean,xrms,xopt,yopt = nstat(nt,varlis[0]+":"+varlis[1],"",
                                                 iretval=1,isilent=1)

          tex = "Mean: " + '{:.4g}'.format(xmean) + \
          "\nRMS: " + '{:.4g}'.format(xrms)

          if yopt != None:
            tex += \
            "\nxOpt: " + '{:.4g}'.format(xopt) + \
            "\nyOpt: " + '{:.4g}'.format(yopt)
          else:
            scom = "yopt = nt." + varlis[1] + ".max()"
            exec(scom)
          #endif

          text(Xstat,Ystat,tex,halign='left')
        except:
          pass
        #endtry

      #endif

    else:

        w = nt[weights] * scalez

        sx = "(" + nparse(nt,varlis[0]) + ") * " + str(scalex)
        sy = "(" + nparse(nt,varlis[1]) + ") * " + str(scaley)
        sw = "(" + nparse(nt,weights) + ") * " + str(scalez)

        if Ihist or Isurf or Itrisurf or Iinter or Iboxes:

          nproj2(nto,varliso,select,weights,scalex=scalex,scaley=scaley,nx=nx,ny=ny,idh=hist)
          hplot2d(hist,plopt)
          iplot = 1

        else:

          s = np.ones_like(w) * Markersize**2

          if cmap == '' or cmap == '!': cmap=Cmap

          sopt = ",s=s, c=" + sw + ", marker='" + Markertype + "',cmap='" + cmap + "'"
          scom = 'plt.scatter(' + sx + ',' + sy + sopt + ')'
          img = eval(scom)

          if Colorbarpad != '!':
            fcm = Fig.colorbar(img, pad=Colorbarpad)
          else:
            fcm = Fig.colorbar(img)
          #endif Colorbarpad != '!'

          fcm.ax.tick_params(labelsize=Axislabelsize)

          Axes.append(fcm)

        #endif Ihist or Isurf or Itrisurf or Iinter or Iboxes

        if Kdump:
          Ndump += 1
          fout = WaveFilePrefix + str(Ndump) + ".dat"
          eval("vwritexyz(" + sx + "," + sy + "," + sw + ",'" + fout + "')")
          print("\nData written to ",fout)
        #endif

    #endif len(weights) = 0

  elif len(varlis) == 3:

    getzone(projection='3d')
    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)

    sx = "(" + nparse(nt,varlis[0]) + ") * " + str(scalex)
    sy = "(" + nparse(nt,varlis[1]) + ") * " + str(scaley)
    sz = "(" + nparse(nt,varlis[2]) + ") * " + str(scalez)

    if Iline == 1 or Iclosed == 1:

      sopt = ", c='" + lcol + "',ls='" + Linestyle + "',lw=" + str(Linewidth)

      if Iclosed:
        scom = "global VxyzX, VxyzY, VxyzZ; VxyzX, VxyzY, VxyzZ = vxyz(" + sx + "," + sy + "," + sz + ")"
        exec(scom)
        sx = list(VxyzX)
        sy = list(VxyzY)
        sz = list(VxyzZ)
        sx.append(sx[0])
        sy.append(sy[0])
        sz.append(sz[0])
        scom='plt.plot(sx,sy,sz' + sopt + ')'
      else:
        scom='plt.plot(' + sx + ',' + sy + ',' + sz + sopt + ')'
      #endif Iclosed

      eval(scom)

    elif Itrisurf:

      getzone(projection='3d')

      Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)
      plt.rcParams['axes.labelsize'] = Atitfontsize3d

      if cmap != '':
        if Icmap.get():
          if cmap == '' or cmap == '!': cmap=Cmap
          surfcolor = None
          scom = 'Ax.plot_trisurf(' + sx + ',' + sy + ',' + sz + ',shade=True, cmap="' +  str(Cmap)  + '")'
        else:
          cmap = None
          if surfcolor == '': surfcolor = Surfcolor
          scom = 'Ax.plot_trisurf(' + sx + ',' + sy + ',' + sz + ',shade=True, color="' +  str(surfcolor)  + '")'
        #endif
      else:
        if color == '' or color == '!' or color == 'default':
          surfcolor = Surfcolor
        else:
          surfcolor = color
        #endif
        scom = 'Ax.plot_trisurf(' + sx + ',' + sy + ',' + sz + ',shade=True, color="' +  str(surfcolor)  + '")'
      #endif

      #Ax.plot_trisurf(x,y,z,shade=True,cmap=cmap,color=surfcolor)
      eval(scom)

    else:
      sopt = ", ls='',  c='" + mcol + "',marker='" + Markertype + "',ms=" + str(Markersize)
      eval('plt.plot(' + sx + ',' + sy + ',' + sz + sopt + ')')
    #endif Iline == 1

    if Kdump:
      Ndump += 1
      fout = WaveFilePrefix + str(Ndump) + ".dat"
      eval("vwritexyz(" + sx + "," + sy + "," + sz + ",'" + fout + "')")
      print("\nData written to ",fout)
    #endif Kdump

  elif len(varlis) == 4:

    getzone(projection='3d')
    Ax.tick_params(labelsize=Axislabelsize, pad=Axislabeldist3d)

    sx = "(" + nparse(nt,varlis[0]) + ") * " + str(scalex)
    sy = "(" + nparse(nt,varlis[1]) + ") * " + str(scaley)
    sz = "(" + nparse(nt,varlis[2]) + ") * " + str(scalez)

    st = "(" + nparse(nt,varlis[3]) + ") * " + str(scalet)

    if cmap == '' or cmap == '!': cmap=Cmap

    #s = np.ones_like(nparse(nt.varlis[3])) * Markersize**2
    s = np.arange(1,len(nt)+1)
    s = np.ones_like(s) * Markersize
    sopt = ",s=s ,c=" + st + ",cmap='" + cmap + "',marker='" + Markertype + "'"
    scom = 'Ax.scatter(' + sx + ',' + sy + ',' + sz + sopt + ')'
    img = eval('Ax.scatter(' + sx + ',' + sy + ',' + sz + sopt + ')')

    if Colorbarpad != '!':
      fcm = Fig.colorbar(img, pad=Colorbarpad)
    else:
      fcm = Fig.colorbar(img)
    #endif

    ttit = varlis[3]
    if scalet != 1: scalet += ' x ' + str(scalet)
    ttit = ""

    fcm.set_label(label=ttit, labelpad=Axistitledist, size=Axislabelsize)
    fcm.ax.tick_params(labelsize=Axislabelsize)

    Axes.append(fcm)

  else:
    print("*** Error in nplot: Don't know how to plot ",varlis," ***")
    return -9
  #if len(varlis) == 0

  if len(legend): Legend.append(legend)
#  print("Klegend in nplot:",Klegend.get())

  Kplots[Kzone-1] = 1
  showplot()

  if Kpdf:
    Npdf += 1
    fout = WaveFilePrefix + str(Npdf) + ".pdf"
    pplot(fout)
    #print("\nPlot written to ",fout)
  #endif Kdump

  Kstat = kstato

#enddef nplot(...) nt idn

def nprof(nt='?',varlis='',select='',weights='',plopt='sprof', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default',isort=0):
          nplot(nt,varlis,select,weights,plopt, legend,scalex, scaley, scalez, scalet,
                cmap, hist,color,isort)
#enddef
def nprofs(nt='?',varlis='',select='',weights='',plopt='samesprof', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default',isort=0):
          nplot(nt,varlis,select,weights,plopt, legend,scalex, scaley, scalez, scalet,
                cmap, hist,color,isort)
#enddef

def nnplot(nt='?',varlis='',select='',weights='',plopt='', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default',isort=0):
  nextzone()
  nplot(nt,varlis,select,weights,plopt, legend,scalex, scaley, scalez, scalet,
        cmap, hist,color,isort)
#enddef nnplot(

def nplt(nt='?',varlis='',select='',weights='',plopt='', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default',isort=0):
  nplot(nt,varlis,select,weights,plopt, legend,scalex, scaley, scalez, scalet,
        cmap, hist,color,isort)

  if type(varlis) == str:
    varlis = nlistcolon(varlis)
    if varlis[0] == '':
      varlis = nt.columns
      vdum = []
      vdum.append(varlis[0])
      if len(varlis) > 1:
        vdum.append(varlis[1])
      varlis = vdum
    #endif varlis[0] == ''
  #endif type(varlis) == str

  tit = ""
  for v in varlis: tit += ":" + v
  tit = tit[1:]

  if select: tit += " (" + select + ")"
  if weights: tit += " * (" + weights + ")"
  txyz(tit)
#enddef nplt

def nnplt(nt='?',varlis='',select='',weights='',plopt='', legend='',
          scalex=1., scaley=1., scalez=1., scalet=1., cmap='', hist='HnPlot',
         color='default',isort=0):
  nextzone()
  nplt(nt,varlis,select,weights,plopt, legend,scalex, scaley, scalez, scalet,
        cmap, hist,color,isort)
#enddef nnplt

def vprint(v):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  Nv.x = v
  print(Nv)
  print("\nmin, max:",Nv.x.min(),Nv.x.max())
  print("sum, mean, sigma:",Nv.x.sum(),Nv.x.mean(),Nv.x.std())
#enddef vprint(xmin,xmax,nx)

def vprintxy(x,y):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  Nxy.x = x
  Nxy.y = y
  print(Nxy)
#enddef vprintxy(xmin,xmax,nx)

def vec(xmin,xmax,nx):
  import numpy as np
  if type(nx) == int:
    return np.linspace(xmin,xmax,nx)
  else:
    print("*** Expected arguments: xmin,xmax,nx")
    return None
    eps = (xmax-xmin) / 1.e10
    return np.arange(xmin,xmax+eps,nx)
  #endif
#enddef vec(xmin,xmax,nx)

def vplxy(x='!',y='!',plopt='',label='',color='!',fillcolor='none'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ


  if type(x) == str and x == '!':
    print("vplxy(x='!',y='!',plopt='',label='',color='!',fillcolor='none')")
    return
  #endif

  if Nwins <=0: window()

  if type(x) == str and x == '!': x = Vxint
  if type(y) == str and y == '!': y = Vyint

  x = np.array(x)
  y = np.array(y)

  plotopt(plopt)
  getzone()

  if Iclosed:
    xl = list(x)
    yl = list(y)
    xl.append(xl[0])
    yl.append(yl[0])
    Iline=1
  #endif

  if color == '!': color = Linecolor

  iplot = 0
  if not Imarker and not Ispline and not Ihist: Imarker=1

  if fillcolor == 'none': fillcolor = getfillcolor()
  else: setfillcolor(fillcolor)

  if Iline == 1:
    # For more sophisticated filling, check plt.fill(...) and polygon patches
    if color == '!': color = Linecolor
    if fillcolor != 'none': plt.fill(x,y,color=fillcolor)
    if Iclosed:
      plt.plot(xl,yl,c=color,ls=Linestyle,lw=Linewidth,label=label)
    else:
      plt.plot(x,y,c=color,ls=Linestyle,lw=Linewidth,label=label)
    #endif
    if Ispline: iplot=0
    else: iplot=1

  elif Imarker == 1:
    # For more sophisticated filling, check plt.fill(...) and polygon patches
    if color == '!': color = Markercolor
    if fillcolor != 'none': plt.fill(x,y,color=fillcolor)
    plt.plot(x,y,c=color,marker=Markertype,fillstyle=Fillstyle, \
    ms=Markersize,ls='')
    if Ispline: iplot=0
    else: iplot=1

  elif plopt == 'h' or Ihist:
    if color == '!': color = Histcolor
    plt.bar(x,y,width=barwidth,color=color,edgecolor=Histedgecolor, fill=Ifill1d)
    if Ispline: iplot=0
    else: iplot=1
  #endif

  if Ispline == 1:

    if type(x) == 'list':  x = np.array(x)
    if type(y) == 'list':  y = np.array(y)
    xmin, xmax = x.min(), x.max()

    nx = len(x) - 1

    if nx < 1000:
      if abs(x[nx]-x[0]) < x.std() *1.0e-9: #RMS
        xspl, yspl = vspline_index(x,y,1001,True)
      else:
        xspl = vec(xmin,xmax,1001)
        yspl = vspline(x,y,xspl)
      #endif
    else:
      xspl = x
      yspl = y
    #endif len(x) < 1001
    # For more sophisticated filling, check plt.fill(...) and polygon patches

    if x.std() > 0.0:
      if color == '!': color = Linecolor
      if fillcolor != 'none': plt.fill(xspl,yspl,color=fillcolor)
      plt.plot(xspl,yspl,c=color,ls=Linestyle,lw=Linewidth,label=label)
      iplot=1
    #endif
  #endif Iline == 1

  if iplot == 0:
    # For more sophisticated filling, check plt.fill(...) and polygon patches
    if fillcolor != 'none': plt.fill(x,y,color=fillcolor)
    plt.plot(x,y,plopt,label=label)
  #endif iplot == 0:

  if label:
    #Klegend.set(1)
    if len(label): Legend.append(label)
  #endif label

  if Kstat:

    if StatFontSize < 0:
      dpi = Fig.dpi
      nxy = max(Nxzone,Nyzone)
      sfs = 12 - nxy * 2
    else:
      sfs = StatFontSize
    #endif StatFontSize < 0

    if Ispline:
      xmin, xmax, xmean, xrms, xopt, yopt = vstat(xspl,yspl)
    else:
      xmin, xmax, xmean, xrms, xopt, yopt = vstat(x,y)
    #endif

    tex = \
    "N, Sum: " + str(int(len(x))) + ", " + '{:.4g}'.format(y.sum()) + \
    "\nMean: " + '{:.4g}'.format(xmean) + \
    "\nMean: " + '{:.4g}'.format(xmean) + \
    "\nRMS: " + '{:.4g}'.format(xrms)

    if xopt != None and yopt != None:
      tex += "\nxOpt: " + '{:.4g}'.format(xopt) + \
      "\nyOpt: " + '{:.4g}'.format(yopt)
    #endif

    text(Xstat,Ystat,tex,halign='left')
    # Latex makes trouble with exponents
    #    latex(Xstat,Ystat,tex,color='black',fontsize=sfs,halign='left',valign='top', \
    #    bbox='none',bbstyle='round',fc='white',ec='black',pad=0.2)

  #endif Kstat

  Kplots[Kzone-1] = 1

  showplot()
#enddef vplxy(x,y,plopt='line',label=''):

def vpll(x='!',y='!',plopt='line',label=''):
  vplxy(x,y,plopt,label,tit,xtit,ytit)
def vplm(x='!',y='!',plopt='marker',label=''):
  vplxy(x,y,plopt,label,tit,xtit,ytit)
def vplls(x='!',y='!',plopt='sameline',label=''):
  vplxy(x,y,plopt,label,tit,xtit,ytit)
def vplms(x='!',y='!',plopt='samemarker',label=''):
  vplxy(x,y,plopt,label,tit,xtit,ytit)

def vplxyey(x,y,ey='',plopt='o',label='',
            marker='o', mfc='', mec='', ms='', mew=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if Nwins <=0: window()

  if len(ey) == 0: ey = y * 0.0
  plotopt(plopt)
  getzone()

  if mfc == '': mfc = Markercolor
  if mec == '': mec = Markercolor # marker edge color
  if ms == '': ms = Markersize

  lincol = mec
  plt.errorbar(x, y, ey, ls='',c=lincol, marker=marker, mfc=mfc, mec=mec, ms=ms, mew=mew)

  if len(label): Legend.append(label)
  showplot()

#enddef vplxyey(x,y,ey,plopt='',label='',...

def vplxyerr(x,y,ey='',ex='',plopt='o',label='',
            marker='o', mfc='', mec='', ms='', mew=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if Nwins <=0: window()

  if len(ey) == 0: ey = y * 0.0
  if len(ex) == 0: ex = x * 0.0

  plotopt(plopt)
  getzone()

  if mfc == '': mfc = Markercolor
  if mec == '': mec = Markercolor # marker edge color
  if ms == '': ms = Markersize

  lincol = mec
  plt.errorbar(x, y, ey, ex, ls='',c=lincol, marker=marker, mfc=mfc, mec=mec, ms=ms, mew=mew)

  if len(label): Legend.append(label)
  showplot()

#enddef vplxyerr(x,y,ey,plopt='',label='',tit='',xtit='',ytit='',...

def vinter(x,y,xint='!'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  n = len(x)

  if type(xint) == list:
    xint = np.array(xint)
  elif type(xint)==str and xint=='!':
    xmin = x.min()
    xmax = x.max()
    xint = vec(xmin,xmax,(max(2,n)-1)*10)
  elif type(xint) == float or type(xint) == int:
    xint = np.array([xint])
  #endif type(xint) == list

  if xint.min() < x.min() or xint.max() > x.max():
    print("*** Error in vinter: x out of range:")
    print("*** x.min(), x.max():",x.min(), x.max())
    print("*** xint.min(), xint.max():",xint.min(), xint.max())
    return -1
  #endif xint.min() < x.min() or xint.max() > x.max():

  nint = len(xint)

  yint = np.zeros_like(xint)

  reset_status()

  if n < 2:
    ErrorText = "*** Error in vinter: Less than two points given ***"
    Istatus = -1
    return yint
  #endif n < 2

  if xint.max() > x.max() or xint.min() < x.min():
    WarningText="*** Warning in vinter: Interpolation interval is out  of range ***"
    Istatus = 1
  #endif xint.max() > x.max() or xint.min() < x.min()

  f = interp1d(x,y)
  yint = f(xint)   # use interpolation function returned by `interp1d`

  Vxint = xint
  Vyint = yint

  return yint
#enddef vinter(x,y,xint)

def vintern(x,y,xint='!'):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ


  n = len(x)

  if type(xint) == list:
    xint = np.array(xint)
  elif type(xint)==str and xint=='!':
    xint=np.linspace(x.min(),x.max(),n*10)
  elif type(xint) == float or type(xint) == int:
    xint = np.array([xint])
  #endif type(xint) == list

  if xint.min() < x.min() or xint.max() > x.max():
    print("*** Error in vintern: x out of range:")
    print("*** x.min(), x.max():",x.min(), x.max())
    print("*** xint.min(), xint.max():",xint.min(), xint.max())
    return -1
  #endif xint.min() < x.min() or xint.max() > x.max():

  nint = len(xint)

  yint = np.zeros_like(xint)
  yintp = np.zeros_like(xint)
  yintpp = np.zeros_like(xint)
  yintint = np.zeros_like(xint)

  Ninter = pd.DataFrame(columns=['x','y','yp','ypp','yint'])

  reset_status()

  if n < 3:
    ErrorText = "*** Error in vintern: Less than three points given ***"
    Istatus = -1
    return yint
  #endif n < 2

  if xint.max() > x.max() or xint.min() < x.min():
    WarningText="*** Warning in vintern: Interpolation interval is out  of range ***"
    Istatus = 1
  #endif xint.max() > x.max() or xint.min() < x.min()

  m = n - 1
  kold = 0

  for k in range(nint):
    xx = xint[k]
    for i in range(m):
      if xx >= x[i] and xx <= x[i+1]: break
    #for i in range(n):
    kold = i
    yint[k] = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i]) * (xx-x[i])
  #endfor i in range(nint):

  kold = 0
  for k in range(0,nint):
    xx = xint[k]
    for i in range(kold,m):
      if xx >= x[i] and xx <= x[i+1]: break
    #for i in range(n):
    kold = i
    if i == 0:
      yintp[k] = (y[1]-y[0])/(x[1]-x[0])
    elif i == nint-1:
      yintp[nint-1] = (y[nint-1]-y[nint-2])/(x[nint-1]-x[nint-2])
    else:
      yintp[k] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
  #endfor i in range(nint):

  kold = 0
  for k in range(kold,nint):
    xx = xint[k]
    for i in range(m):
      if xx >= x[i] and xx <= x[i+1]: break
    #for i in range(n):
    kold = i
    if i == 0:
      yintpp[k] = \
      ((y[2]-y[1])/(x[2]-x[1]) - (y[1]-y[0])/(x[1]-x[0])) / \
      ((x[2]+x[1])/2. - (x[1]+x[0])/2.)
    elif i == nint-2:
      yintpp[k] = \
      ((y[nint-1]-y[nint-2])/(x[nint-1]-x[nint-2]) - (y[nint-2]-y[nint-3])/ \
      (x[nint-2]-x[nint-3])) / \
      ((x[nint-1]+x[nint-2])/2. - (x[nint-2]+x[nint-3])/2.)
    else:
     yintpp[k] = \
      ((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/ (x[i]-x[i-1])) / \
      ((x[i+1]+x[i])/2. - (x[i]+x[i-1])/2.)
   #endfor i in range(nint):

  for k in range(nint-1):
    dx = xint[k+1] - xint[k]
    yintint[k+1] = yintint[k]+dx*(yint[k]+dx*(yintp[k]/2.+yintpp[k]*dx/6.))
  #endfor k in range(nint)

  Ninter.x = xint
  Ninter.y = yint
  Ninter.yp = yintp
  Ninter.ypp = yintpp
  Ninter.yint = yintint

  return yint
#enddef vintern(x,y,xint)

def vspline_index(x,y,nspl=1001, periodic=False, ypp1=0.0, yppn=0.0):

  import numpy as np

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ


  if type(nspl) != int:
    print("*** Error in vspline_index: Argument nspl must be integer!")
    return x*0
  #endif

  nspl = max(3,nspl)

  nx = len(x)
  ind = np.linspace(0.,float(nx-1),nx)

  yspl = vspline(ind,y,nspl, periodic, ypp1, yppn)
  yp = deepcopy(Nspline.yp)
  ypp = deepcopy(Nspline.ypp)

  xspl = vspline(ind,x,nspl, periodic, ypp1, yppn)

  xp = Nspline.yp+(1.0e-15*Nspline.yp.std())
  xpp = deepcopy(Nspline.ypp)

  ssl = deepcopy(VsplX)

  VsplX = deepcopy(xspl)
  VsplY = deepcopy(yspl)
  Vspl1 = yp / xp
  Vspl2 = np.zeros_like(yp)
  VsplI = np.zeros_like(yp)


  Vspl2 = (ypp*xp-yp*xpp) / xp

  VsplI[0] = 0.0
  for i in range(len(xspl)-1):
    VsplI[i+1]=VsplI[i] + (VsplX[i+1]-VsplX[i])*0.5 *(VsplY[i]+VsplY[i+1]) \
    -(VsplX[i+1]-VsplX[i])**3/24.*(Vspl2[i]+Vspl2[i+1])
  #endfor i in range(len(xspl)-1)

  Nspline = ncre("Nspline","x:y:yp:ypp:yint",ioverwrite=1)

  Nspline.x = VsplX
  Nspline.y = VsplY
  Nspline.yp = Vspl1
  Nspline.ypp = Vspl2
  Nspline.yint = VsplI

  Nspline.index = range(len(Nspline))
  nupdate_header("Nspline")

  return xspl,yspl
#enddef vspline_index(...

def vspline(x,y,xspl='!', periodic=False, ypp1=0.0, yppn=0.0):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ


  if SplineMode.lower() != 'new': return vspline_old(x,y,xspl,periodic)

  import numpy as np

  #nreakpoint()
  n = len(x)

  if n < 2:
    ErrorText = "*** Error in vspline: Less than two points given ***"
    Istatus = -1
    return y
  #endif n < 2

  if type(x) == list: x = np.array(x)
  if type(y) == list: y = np.array(y)

  if x.std() == 0.0:
    ErrorText = "*** Error in vspline: Zero RMS of x values ***"
    Istatus = -1
    return y
  #endif

  if type(xspl) == list:
    xspl = np.array(xspl)
  elif type(xspl) == int:
    nspl = xspl
    if nspl < 10: nspl = 10
    xspl=np.linspace(x.min(),x.max(),nspl)
  elif type(xspl) == str and xspl == '!':
    xspl=np.linspace(x.min(),x.max(),n*10)
  #endif type(xspl) == list

  try:
    nspl = len(xspl)
  except:
    xpd = pd.DataFrame([xspl])
    xpd.columns = ['x']
    nspl = len(xpd)
    xspl = xpd.x
  #endtry

  xs = x
  ys = y
  x,y = vsortxy(xs,ys)

  xspl = vsortx(xspl)

  VsplX = xspl
  VsplY = np.zeros_like(xspl)
  Vspl1 = np.zeros_like(xspl)
  Vspl2 = np.zeros_like(xspl)
  VsplI = np.zeros_like(xspl)

  reset_status()

  if n < 2:
    ErrorText = "*** Error in vspline: Less than two points given ***"
    Istatus = -1
    return VsplY
  #endif n < 2

  if VsplX.max() > x.max() or VsplX.min() < x.min():
    WarningText="*** Warning in vspline: Interpolation interval is out of range ***"
    Istatus = 1
  #endif VsplX.max() > x.max() or VsplX.min() < x.min()

  if periodic:
    kfail, y2p = util_spline_coef_periodic(x,y)
  else:
    kfail, y2p = util_spline_coef(x,y,ypp1,yppn)
  #endif periodic

  Istatus += kfail

  yy,yyp,yypp = util_spline_inter(x,y,y2p,xspl[0],-1)
  for i in range(len(xspl)):
    VsplY[i],Vspl1[i],Vspl2[i] = util_spline_inter(x,y,y2p,xspl[i],0)
    #print(i,xspl[i],VsplY[i])
  #endfor i in range(len(xspl)):

  VsplI[0] = 0.0
  for i in range(len(xspl)-1):
    VsplI[i+1]=VsplI[i] + (VsplX[i+1]-VsplX[i])*0.5 *(VsplY[i]+VsplY[i+1]) \
    -(VsplX[i+1]-VsplX[i])**3/24.*(Vspl2[i]+Vspl2[i+1])
  #endfor i in range(len(xspl)-1)

  Nspline = ncre("Nspline","x:y:yp:ypp:yint",ioverwrite=1)

  Nspline.x = VsplX
  Nspline.y = VsplY
  Nspline.yp = Vspl1
  Nspline.ypp = Vspl2
  Nspline.yint = VsplI

  Nspline.index = range(len(Nspline))
  nupdate_header("Nspline")

  VsplCoef = y2p

  return VsplY
#enddef vspline(x,y,xspl)

def vspline_old(x,y,xspl='!', periodic=False):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ


  import numpy as np

  n = len(x)

  if type(xspl) == list:
    xspl = np.array(xspl)
  elif type(xspl) == int:
    if xspl < 10: xspl = 10
    xspl=np.arange(x.min(),x.max(),xspl)
  elif type(xspl) == str and xspl == '!':
    xspl=np.arange(x.min(),x.max(),xspl*10)
  #endif type(xspl) == list

  try:
    nspl = len(xspl)
  except:
    xpd = pd.DataFrame([xspl])
    xpd.columns = ['x']
    nspl = len(xpd)
    xspl = xpd.x
  #endtry

  xs = x
  ys = y
  x,y = vsortxy(xs,ys)
  xspl = vsortx(xspl)

  VsplX = xspl
  VsplY = np.zeros_like(xspl)
  Vspl1 = np.zeros_like(xspl)
  VsplI = np.zeros_like(xspl)

  reset_status()

  if n < 2:
    ErrorText = "*** Error in vspline_old: Less than two points given ***"
    Istatus = -1
    return VsplY, istat
  #endif n < 2

  if VsplX.max() > x.max() or VsplX.min() < x.min():
    WarningText="*** Warning in vspline_old: Interpolation interval is out of range ***"
    Istatus = 1
  #endif VsplX.max() > x.max() or VsplX.min() < x.min()

  if n > 3:

    VsplCoef = interpolate.splrep(x,y, s=0, per=periodic) # from scipy

    VsplY = interpolate.splev(xspl, VsplCoef, der=0)
    Vspl1 = interpolate.splev(xspl, VsplCoef, der=1)
    Vspl2 = interpolate.splev(xspl, VsplCoef, der=2)

    if nspl > 1:
      const = interpolate.splint(0, xspl[0], VsplCoef)
      for i in range(nspl):
        VsplI[i] = interpolate.splint(0, xspl[i], VsplCoef) - const
      #endfor i in range(nspl)
    #endif nspl > 1

  else:

    if n > 2:
      f2 = interp1d(x, y, kind='quadratic')
    elif n > 1:
      f2 = interp1d(x, y)
    #endif n > 2

    x4 = vec(x[0],x[n-1],4)
    y4 = f2(x4)

    VsplCoef = interpolate.splrep(x4,y4, s=0) # from scipy

    VsplY = interpolate.splev(xspl, VsplCoef, der=0)
    Vspl1 = interpolate.splev(xspl, VsplCoef, der=1)
    Vspl2 = interpolate.splev(xspl, VsplCoef, der=2)

    if nspl > 1:
      const = interpolate.splint(0, xspl[0], VsplCoef)
      for i in range(nspl):
        VsplI[i] = interpolate.splint(0, xspl[i], VsplCoef) - const
      #endfor i in range(nspl)
    #endif nspl > 1

  #endif n > 3

  #Nspline = pd.DataFrame(columns=['x','y','yp','ypp','yint'])
  Nspline = ncre("Nspline","x:y:yp:ypp:yint",ioverwrite=1)

  Nspline.x = VsplX
  Nspline.y = list(VsplY)
  Nspline.yp = Vspl1
  Nspline.ypp = Vspl2
  Nspline.yint = VsplI

  Nspline.index = range(len(Nspline))
  nupdate_header("Nspline")

  return VsplY
#enddef vspline_old(x,y,xspl)

def vcopn(nt,varlis='x',x='',y='',z='',s='', t='',bx='',by='',bz=''):

  idn = GetIndexN(nt,1)
  varlis = nlistcolon(varlis)

  if idn == -1:
    if type(nt) != str:
      print("*** Error in vcopn: nt must be Ntuple or Ntuple name ***")
      return -1
    else:
      nt = ncre(nt,nt,varlis)
    #endif type(nt) != str
  #if type(nt) != str

  nvar = len(varlis)

  nt[varlis[0]] = x
  if nvar > 1: nt[varlis[1]] = y
  if nvar > 2: nt[varlis[2]] = z
  if nvar > 3: nt[varlis[3]] = s
  if nvar > 4: nt[varlis[4]] = t
  if nvar > 5: nt[varlis[5]] = bx
  if nvar > 6: nt[varlis[6]] = by
  if nvar > 7: nt[varlis[7]] = bz

  ind = Ind

  nhead = Nhead[ind]
  nvar=nhead[3]

  varlis = []
  for i in range(nvar):
    varlis.append(nhead[4+i][0])
  #endfor i in range(nvar):

  nt.columns = varlis

  for i in range(nvar):
    nhead[4+i][1] = nt[varlis[i]].min()
    nhead[4+i][2] = nt[varlis[i]].max()
  #endfor i in range(nvar):

  nhead[4+nvar]=len(nt)

  return nt
#def vcopn(nt,varlis='x:y:z:bx:by:bz',x,y='',z='',bx='',by='',bz='')

def nupdate_header(nt,reindex=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  idn = -1
  idn = GetIndexN(nt,0)

  if idn == -1: return -1
  nheado = Nhead[idn]
  nt = Ntup[idn]

  if reindex:
    #nt = nt.reset_index()
    nt.index = range(len(nt))
    Ntup[idn] = nt
  #endif reindex:

  headnt = nheado[0:3]
  varlis = nt.columns
  nvar = len(varlis)
  headnt.append(nvar)

  for i in range(nvar):
    try:
      v = nt[varlis[i]];
      vmin=v.min(); vmax=v.max()
      vmean=v.mean(); vstd=v.std()
      if str(vstd) == 'nan': vstd = 0.0
      var = [varlis[i],vmin,vmax,vmean,vstd]
    except:
      var = [varlis[i],None,None,None,None]
    #endtry
    headnt.append(var)
  #endfor i in range(nvar)

  nent = len(nt)
  headnt.append(nent)

  Nhead[idn] = headnt

  return
#enddef nupdate_header(nt)

def ndelrow(nt,irow1,irow2,reindex=1):
  idn = GetIndexN(nt)
  if idn == -1: return -1

  nt = Ntup[idn]

  nt = nt.drop(nt.index[irow1:irow2+1])
  Ntup[idn] = nt

  nupdate_header(Ntup[idn],reindex)

  return nt

#def ndelrow(nt,irow1,irow2,reindex=1)

def naddcol(nt,newnam='',newdef=''):

  if not newnam: return
  if not newdef: return

  ind = GetIndexN(nt)
  if ind == -1: return -1

  ntnam = Nhead[ind][1]
  nt = Ntup[ind]

  snew = nparse(nt,newdef)
  print(nparse(nt,newdef))
  com = "nt['" + newnam+ "'] = " + snew
  exec(com)

  nupdate_header(nt)

  print(nt)

#enddef

def ndelcol(nt,varlis=''):

  if not varlis: return

  ind = GetIndexN(nt)
  if ind == -1: return -1

  nt = Ntup[ind]
  varlis = nlistcolon(varlis)

  for i in range(len(varlis)): del nt[varlis[i]]

  nupdate_header(nt)

  Ntup[ind] = nt

  return nt
#enddef ndelcol(nt,varlis='')

def vsolve(x,y,val=0.0,xmin=-1.0e30,xmax=1.0e30):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  reset_status()
  #nreakpoint()

  x = np.array(x)
  y = np.array(y)

  n = len(x)

  if n < 2:
    WarningText = "*** Error in vsolve: Not enough data ***"
    Istatus = 1
    return -9999.
  elif n == 2:
    try:
      a = (x[1]-x[0])/(y[1]-y[0])
      xs = x[0] + a * (val - y[0])
    except:
      print("\n *** Error in vsolve: Problems with vsolve. Are x,y monton?")
      xs = -9999.
    #endtry
    return float(xs)
  elif n == 3:
    try:
      a, yp, opt, ifail = util_parabel(y,x)
      xs = a[0] + a[1]*val + a[2]*val*val
    except:
      print("\n *** Error in vsolve: Problems with vsolve. Are x,y monton?")
      xs = -9999.
    #endtry
    return float(xs)
  #endif

  if x.min() >= xmin and x.max() >= xmax:
    xspl = x
    yspl = y
  else:
    xpl = []; ypl = []
    for i in range(len(x)):
      if x[i] >= xmin and x[i] <= xmax:
        xpl.append(x[i])
        ypl.append(y[i])
      #endif
    #endfor
  #endif

  n = len(xpl)
  ysort,xsort = vsortxy(ypl,xpl)
  #nreakpoint()

  try:
    coef = interpolate.splrep(ysort,xsort, s=0) # from scipy
    if val < y.min() or val > y.max():
      WarningText = "*** Error in vsolve: Value out of range ***"
      Istatus = 1
    #endif val < y.min() or val > y.max()
    xs = interpolate.splev(val, coef, der=0)
  except:
    print("\n *** Error in vsolve: Problems with vspline. Are x,y monton?")
    xs = -9999.
  #endtry

  return float(xs)
#enddef vsolve(x,y,val):

def vsolvelin(x,y,val=0.0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  reset_status()

  f = interp1d(y,x)
  xs = f(val)

  return xs
#enddef vsolve(x,y,val):

def voptspl(x,y):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(x) == list: x = np.array(x)
  if type(y) == list: y = np.array(y)

  n = len(x)

  xspl = x
  yspl = y

  if n < 100:
    xspl = vec(x.min(),x.max(),100)
    yspl = vspline(x,y,xspl)
  #endif n < 100

  yspl = vspline(xspl,yspl,xspl)
  xopt, yopt = voptpar(xspl,yspl)

  return xopt,yopt
#enddef voptspl(x,y):

def set_error(text,istat=-1):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  Istatus = istat; ErrorText=text
#def set_error()

def set_warning(text,istat=1):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  Istatus = istat; WarningText=text
#def set_error()

def reset_status():
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  Istatus = 0; ErrorText=''; WarningText=''
#def reset_status()

def ncopn(nt,ncnam,varlis='',select='',ioverwrite=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  reset_status()

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt = N
  nt.index = range(len(nt))

  v = nt.columns

  nc = nclone(nt,ncnam,ncnam,ioverwrite)

  if varlis == '':
    nc.columns = list(v)
  else:
    Quit(nc)

    nc.columns = nlistcolon(varlis)
    Quit()
  #endif varlis == ''

  for k in range(len(v)):
    nc[v[k]] = nt[v[k]]
  #endfor k in range(ncn)

  nupdate_header(nc)

  return nc

#enddef ncopn(nt,"nc",varlis,select='',ioverwrite=0)

def ncopv(nt,varlis,select=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  reset_status()

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(varlis) == 0: varlis = list(N.columns)
  varlis = nlistcolon(varlis)
  n = len(varlis)

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif len(select)

  nt=N
  var = []
  svar = ""

  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
    var.append(sv)
  #endfor

  sv = nparse(nt,varlis[Ncolon])
  svar += sv
  var.append(sv)

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  exec(scom)

  N.index = range(len(N))
  v = N.columns

  if n == 1: return N[v[0]]
  elif n == 2: return N[v[0]], N[v[1]]
  elif n == 3: return N[v[0]], N[v[1]],N[v[2]]
  elif n == 4: return N[v[0]], N[v[1]],N[v[2]], N[v[3]]
  elif n == 5: return N[v[0]], N[v[1]],N[v[2]], N[v[3]],N[v[4]],
  elif n == 6: return N[v[0]], N[v[1]],N[v[2]], N[v[3]],N[v[4]], \
  N[v[5]]
  elif n == 7: return N[v[0]], N[v[1]],N[v[2]], N[v[3]],N[v[4]], \
  N[v[5]],N[v[6]]
  elif n == 8: return N[v[0]], N[v[1]],N[v[2]], N[v[3]],N[v[4]], \
  N[v[5]],N[v[6]], N[v[7]]

#def ncopv()

def nclone(nt,ncnam,nctit='',ioverwrite=0):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  reset_status()

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if nctit == '': nctit = ncnam
  vlis = list(N.columns)

  varlis = ""
  for i in range(len(vlis)-1) :
    v = vlis[i]
    varlis += v + ':'
  #endfor v in vlis
  varlis += vlis[len(vlis)-1]

  nc = ncre(ncnam,nctit,varlis,ioverwrite)

  return nc
#def nclone()

def nsort(nt,var):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug


  reset_status()

  if type(nt) == str:
    idn = GetIndexN(nt)
    if idn == -1: return -1
    nt = N
  #endif type(nt) == str

  nso =  nt.sort_values(by=var)
  nso.index = range(len(nt))
  return nso

#def nsort()

def vxy(x,y):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()
  nt = make_dataframe('x:y',x,y)
  nt.index = range(len(nt))
  return nt.x, nt.y
#def vxy()

def vxyz(x,y,z):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()
  nt = make_dataframe('x:y:z',x,y,z)
  nt.index = range(len(nt))
  return nt.x, nt.y, nt.z
#def vxyz()

def vsortxy(x,y):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  nt = make_dataframe('x:y',x,y)
  nt = nt.sort_values(by='x')
  nt.index = range(len(nt))

  return nt.x, nt.y

#def vsortxy()

def vflip(v):
  w = deepcopy(v)
  n = len(w)
  for i in range(n): w[i] = v[n-i-1]
  return w
#enddef

def vsortx(x):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  nt = make_dataframe('x',x)
  nt = nt.sort_values(by='x')
  nt.index = range(len(nt))

  return nt.x
#def vsortx()

def vwritexy(x,y,filename='vxy.dat',ifirst=-1,ilast=-1,nlines=-1):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug


  reset_status()

  if ifirst == -1: ifirst = 0

  if ilast == -1:
    if nlines == -1:
      ilast = len(x)
    else:
      ilast = ifirst + nlines
    #endif
  #endif

  nt = make_dataframe('x:y',x[ifirst:ilast],y[ifirst:ilast])
  nt.to_csv(filename,header=False,index=False,sep=' ')
#def vwritexy()

def vreadx(fname='ndump.dat',header=None, skiphead=0, skipfoot=0, comment='*', sep=' '):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  if type(fname) == str:
    if sep == ' ' or sep == '':
      nt = pd.read_csv(fname,header=header, delim_whitespace=True, skiprows=skiphead, comment='*', skipfooter=skipfoot)
    else:
      nt = pd.read_csv(fname,header=header, sep=sep, skiprows=skiphead, comment='*', skipfooter=skipfoot)
  else:
    nt = pd.DataFrame(fname)

  nt.columns = ['x']
  return nt.x
#def vreadxy()

def vreadxy(fname='ndump.dat',header=None, skiphead=0, skipfoot=0, comment='*', sep=' '):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  if sep == ' ' or sep == '':
    nt = pd.read_csv(fname,header=header, delim_whitespace=True, skiprows=skiphead, comment='*', skipfooter=skipfoot)
  else:
    nt = pd.read_csv(fname,header=header, sep=sep, skiprows=skiphead, comment='*', skipfooter=skipfoot)

  nt.columns = ['x','y']
  return nt.x, nt.y
#def vreadxy()

def vwritex(x,filename='vx.dat'):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()
  nt = make_dataframe('x',x)
  nt.to_csv(filename,header=False,index=False,sep=' ')
#def vwritexy()

def vwritexyz(x,y,z,filename='vxyz.dat',ifirst=-1,ilast=-1,nlines=-1):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()
  if ifirst == -1: ifirst = 0

  if ilast == -1:
    if nlines == -1:
      ilast = len(x)
    else:
      ilast = ifirst + nlines
    #endif
  #endif

  nt = make_dataframe('x:y:z',x[ifirst:ilast],y[ifirst:ilast],z[ifirst:ilast])
  nt.to_csv(filename,header=False,index=False,sep=' ')
#def vwritexyz()

def vwritexyzt(x,y,z,t,filename='vxyzt.dat',ifirst=-1,ilast=-1,nlines=-1):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  if ilast == -1:
    if nlines == -1:
      ilast = len(x)
    else:
      ilast = ifirst + nlines
    #endif
  #endif

  nt = make_dataframe('x:y:z:t',x[ifirst:ilast],y[ifirst:ilast],z[ifirst:ilast],t[ifirst:ilast])
  nt.to_csv(filename,header=False,index=False,sep=' ')
#def vwritexyz()

def vsmooth(vx,nsmooth=0):

# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  ndim = len(vx)
  x = vx
#  x.index = range(ndim)

  smooth = x * 0.0
  #print(ndim,len(x),smooth[0],smooth[ndim-1])

  nsmooth=min(nsmooth,ndim)
  nsmooth=int(nsmooth/2)
  nsmooth1=int(2*nsmooth+1)

  sum=0.0

  for i in range(nsmooth1):
    sum=sum+x[i]
    smooth[i] = x[i]
  #endfor i in range(nsmooth1)

  i = ndim - nsmooth - 1
  while i <= ndim - 1:
    smooth[i]=x[i]
    i+=1
  #endwhile i <= ndim - 1

  n1 = - 1
  n2=nsmooth1 - 1
  i=nsmooth
  while n2 < ndim - 1 :
    i=i+1
    n1=n1+1
    n2=n2+1
    sum=sum-x[n1]+x[n2]
    smooth[i]=sum/nsmooth1
  #endwhile n2 < ndim - 1

  return smooth

#enddef vsmooth(x,nsmooth=0)

def vpeaks(x,y,pkmin=0.5,nsmooth=0,isilent=0):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  global Ical

  fmxtot=-1.0e30
  ndim = len(x)

  nsmooth=min(nsmooth,ndim)
  nsmooth=int(nsmooth/2)
  nsmooth1=int(2*nsmooth+1)

  smooth = y * 0.0
  sum=0.0

  for i in range(nsmooth1):
    sum=sum+y[i]
    smooth[i] = y[i]
  #endfor i in range(nsmooth1)

  i = ndim - nsmooth - 1
  while i <= ndim - 1:
    smooth[i]=y[i]
    i+=1
  #endwhile i <= ndim - 1

  n1 = - 1
  n2=nsmooth1 - 1
  i=nsmooth
  while n2 < ndim - 1 :
    i=i+1
    n1=n1+1
    n2=n2+1
    sum=sum-y[n1]+y[n2]
    smooth[i]=sum/nsmooth1
    if smooth[i] > fmxtot: fmxtot=smooth[i]

  if not isilent:
    print('\nPeaks (Number, center, height, width, distance to previous peak):')
    print("(The width sigma is calculate from the peak parabola:  f(xpeak +/- sigma) = max/2.)")
  #endif not isilent

  npeaks=0
  thresh=fmxtot*pkmin

  ixpeaks =  []
  xpeaks =  []
  ypeaks =  []
  sigma =  []

  i = 1
  while i < ndim-1:

# calculate s=a0+a1*(x-x0)+a2*(x-x0)**2
# change system: (x0,s0)->(0,0), i.e.
# calculate s=a1*dx+a2*dx**2
#  ds/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

    x0=x[i]
    dxm=x[i-1]-x0
    dxp=x[i+1]-x0

    s0=smooth[i]
    sm=smooth[i-1]-s0
    sp=smooth[i+1]-s0

    if s0 >= thresh and sm < 0.0 and sp <= 0.0:

      npeaks=npeaks+1

      dxm2=dxm*dxm
      dxp2=dxp*dxp

      det=dxm*dxp2-dxp*dxm2
      if det:
        a1=(sm*dxp2-sp*dxm2)/det
        a2=(sp*dxm-sm*dxp)/det
      else:
        a2=0.0

      if a2 < 0.0:
        dxmax=-a1/(2.0*a2)
        dymax=(a1+a2*dxmax)*dxmax
        ixpeaks.append(i)
        xpeaks.append(x0+dxmax)
        ypeaks.append(s0+dymax)
        sigma.append((-abs(s0)/(2.0*a2))**0.5)
      else:
        ixpeaks.append(i)
        xpeaks.append(x[i])
        ypeaks.append(smooth[i])
        sigma.append(-9999.0)
      #endif a2 < 0.0

      n1 = npeaks -1
      if not isilent:
        if npeaks == 1:
          print(npeaks,"%10.5g"%xpeaks[n1], "%10.5g"%ypeaks[n1], "%10.5g"%sigma[n1])
        else:
          print(npeaks,"%10.5g"%xpeaks[n1], "%10.5g"%ypeaks[n1], "%10.5g"%sigma[n1],"%10.5g"%(xpeaks[n1]-xpeaks[n1-1]))
      #endif not isilent

    #endif s0 >= thresh and sm < 0.0 and sp <= 0.0

    i+=1
  #end while i < ndim-1

  return ixpeaks,xpeaks,ypeaks,sigma
#enddef vpeaks()

def vpeaksabs(x,y,pkmin=0.5,nsmooth=0,isilent=0):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  fmxtot=-1.0e30
  ndim = len(x)

  nsmooth=min(nsmooth,ndim)
  nsmooth=int(nsmooth/2)
  nsmooth1=int(2*nsmooth+1)

  y = abs(y)
  smooth = y

  sum=0.0
  for i in range(nsmooth1):
    sum=sum+y[i]
    smooth[i]=y[i]
  #endfor i in range(nsmooth1)

  i = ndim - nsmooth - 1
  while i <= ndim - 1:
    smooth[i]=y[i]
    i+=1

  n1 = - 1
  n2=nsmooth1 - 1
  i=nsmooth
  while n2 < ndim - 1 :
    i=i+1
    n1=n1+1
    n2=n2+1
    sum=sum-y[n1]+y[n2]
    smooth[i]=sum/nsmooth1
    if smooth[i] > fmxtot: fmxtot=smooth[i]

  if not isilent:
    print('\nPeaks (center, height, width, distance to previous peak):')
    print("Width sigma:  f(xpeak +/ -sigma) = max/2.")

  npeaks=0
  thresh=fmxtot*pkmin

  ixpeaks = []
  xpeaks =  []
  ypeaks =  []
  sigma =  []

  i = 1
  while i < ndim-1:

# calculate s=a0+a1*(x-x0)+a2*(x-x0)**2
# change system: (x0,s0)->(0,0), i.e.
# calculate s=a1*dx+a2*dx**2
#  ds/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

    x0=x[i]
    dxm=x[i-1]-x0
    dxp=x[i+1]-x0

    s0=smooth[i]
    sm=smooth[i-1]-s0
    sp=smooth[i+1]-s0

    if s0 >= thresh and sm < 0.0 and sp <= 0.0:

      npeaks=npeaks+1

      dxm2=dxm*dxm
      dxp2=dxp*dxp

      det=dxm*dxp2-dxp*dxm2
      if det:
        a1=(sm*dxp2-sp*dxm2)/det
        a2=(sp*dxm-sm*dxp)/det
      else:
        a2=0.0

      if a2 < 0.0:
        dxmax=-a1/(2.0*a2)
        dymax=(a1+a2*dxmax)*dxmax
        ixpeaks.append(i)
        xpeaks.append(x0+dxmax)
        ypeaks.append(s0+dymax)
        sigma.append((-abs(s0)/(2.0*a2))**0.5)
      else:
        ixpeaks.append(i)
        xpeaks.append(x[i])
        ypeaks.append(smooth[i])
        sigma.append(-9999.0)
      #endif a2 < 0.0

      n1 = npeaks -1
      if not isilent:
        if npeaks == 1:
          print(npeaks,"%10.5g"%xpeaks[n1], "%10.5g"%ypeaks[n1], "%10.5g"%sigma[n1])
        else:
          print(npeaks,"%10.5g"%xpeaks[n1], "%10.5g"%ypeaks[n1], "%10.5g"%sigma[n1],"%10.5g"%(xpeaks[n1]-xpeaks[n1-1]))

    #endif s0 >= thresh and sm < 0.0 and sp <= 0.0

    i+=1
  #end while i < ndim-1

  return ixpeaks,xpeaks,ypeaks,sigma
#def vpeaksabs()

def hfwhm(h='?',nsmooth=0,isilent=0):

  if type(h) == str and h == '?':
    print("\nUsage: fwhm = hfwhm(histo, nsmooth=0, isilent=0")

  vxy = h1copv(h)
  fwhm = vfwhm(vxy.x,vxy.y,nsmooth,isilent)

  return fwhm
#enddef hfwhm(h='?',nsmooth=0,isilent=0):

def hfwhmabs(h='?',nsmooth=0,isilent=0):

  if type(h) == str and h == '?':
    print("\nUsage: fwhm = hfwhmabs(histo, nsmooth=0, isilent=0")

  vxy = h1copv(h)
  fwhm = vfwhm(vxy.x,abs(vxy.y),nsmooth,isilent)

  return fwhm
#enddef hfwhmabs(h='?',nsmooth=0,isilent=0):

def nfwhmabs(nt='?',varlis='', select='', nsmooth=0,isilent=0):

  if type(nt) == str and nt == '?':
    print("\nUsage: fwhm = nfwhmabs(Ntup, varlis='', select='', nsmooth=0, isilent=0")
    return
  #endif type(nt) == str and nt == '?'

  vx,vy = ncopv(nt,varlis,select)
  fwhm = vfwhmabs(vx,vy,nsmooth,isilent)

  return fwhm
#enddef nfwhmabs(h='?',nsmooth=0,isilent=0):

def nfwhm(nt='?',varlis='', select='', nsmooth=0,isilent=0):

  global N

  if type(nt) == str and nt == '?':
    print("\nUsage: fwhm = nfwhm(Ntup, varlis='', select='', nsmooth=0, isilent=0")
    return
  #endif type(nt) == str and nt == '?'

  if type(varlis) == str:
    if not len(varlis):
      print("*** Error in nfwhm: No columns specified in var ***")
      return
    varlis = nlistcolon(varlis)

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N
  #endif

  svar = ""
  svl = []

  nt = N
  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
  svar += nparse(nt,varlis[Ncolon])

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  Nlines = len(N)
  N.columns = ['x','y']
  N.index = range(Nlines)

  fwhm = vfwhm(N.x,N.y,nsmooth,isilent)

  return fwhm
#enddef nfwhm(h='?',nsmooth=0,isilent=0):

def nfwhmtot(nt='?',varlis='', select='', nsmooth=0,isilent=0):

  global N

  if type(nt) == str and nt == '?':
    print("\nUsage: fwhm = nfwhmtot(Ntup, varlis='', select='', nsmooth=0, isilent=0")
    return
  #endif type(nt) == str and nt == '?'

  if type(varlis) == str:
    if not len(varlis):
      print("*** Error in nfwhmtot: No columns specified in var ***")
      return
    varlis = nlistcolon(varlis)

  if type(nt) == Tdf:
    N = nt
  elif type(nt) == str:
    ind = GetIndexN(nt)
    if ind == -1:
      return -1
  elif type(nt) == int and nt < 0:
    pass
  else:
    return -2
  #endif nt >= 0:

  if len(select):
    N = N.query(select)
    Nsel = N

  svar = ""
  svl = []

  nt = N
  for k in range(Ncolon):
    sv = nparse(nt,varlis[k])
    svar += sv + ","
  svar += nparse(nt,varlis[Ncolon])

  scom = "global N; N = pd.DataFrame([" + svar + "]).T"
  #print(scom)
  exec(scom)

  Nlines = len(N)
  N.columns = ['x','y']
  N.index = range(Nlines)

  fwhm = vfwhmtot(N.x,N.y,nsmooth,isilent)

  return fwhm
#enddef nfwhmtot(h='?',nsmooth=0,isilent=0):

def vfwhmabs(x='?',y='',nsmooth=0,isilent=0):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug


  reset_status()

  if type(x) == str and x == '?':
    print("\nUsage: fwhm = vfwhmabs(x,y,nsmooth=0, isilent=0")

  return vfwhm(x,abs(y),nsmooth,isilent)

#enddef vfwhmabs(x='?',y='',nsmooth=0,isilent=0)

def vfwhm(x='?',y='',nsmooth=0,isilent=0):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug

  reset_status()

  if type(x) == str and x == '?':
    print("\nUsage: fwhm = vfwhm(x,y,nsmooth=0, isilent=1")

  ixpeaks,xpeaks,ypeaks,sigma = vpeaks(x,y,0.5,nsmooth,1)

  fwhm = []
  npeaks = len(ixpeaks)
  npoi = len(x)

  for i in range(npeaks):

    if sigma[i] == -9999.: continue

    im = -9999
    ip = -9999

    k = ixpeaks[i]
    ymax = ypeaks[i]
    half = ymax / 2.

    for l in range(0,k+1):
      if y[l] >= half and l > 0:
        im = l - 1
        break
    #endif y[l] >= half

    for l in range(k,npoi):
      if y[l] <= half and l < npoi:
        ip = l
        break
    #endif y[l] >= half

    xm = -9999.
    xp = -9999.
    w = -9999.

    if im != -9999 and y[im] <= half:
      x1 = x[im]
      x2 = x[im+1]
      y1 = y[im]
      y2 = y[im+1]
      xm = x1 + (x2-x1)/(y2-y1) * (half-y1 )
    #endif im != -9999 and ip != -9999

    if ip != -9999 and y[ip]  < half:
      x1 = x[ip-1]
      x2 = x[ip]
      y1 = y[ip-1]
      y2 = y[ip]
      xp = x1 + (x2-x1)/(y2-y1) * (half-y1 )
    #endif im != -9999 and ip != -9999

    if xm != -9999. and xp != -9999: w = xp - xm

    fwhm.append([i,w,xm,xpeaks[i],xp])

  #endfor i in range(len(ixpeaks))

  if not isilent:

    print("\n  N      FHWH         Xm       Xpeak       Xp       Wm         Wp")
    print("----------------------------------------------------------------------")

    for k in range(npeaks):
      if fwhm[k][0] == -9999: continue
      wm = -9999.
      wp = -9999.
      if fwhm[k][3] != -9999. and fwhm[k][2] != -9999.:
        wm = fwhm[k][3]-fwhm[k][2]
      if fwhm[k][3] != -9999. and fwhm[k][4] != -9999.:
        wp = fwhm[k][4]-fwhm[k][3]
      print(" ",k+1,\
      '     {:.4g}'.format(fwhm[k][1]),\
      '     {:.4g}'.format(fwhm[k][2]),\
      '     {:.4g}'.format(fwhm[k][3]),\
      '     {:.4g}'.format(fwhm[k][4]),\
      '     {:.4g}'.format(wm),\
      '     {:.4g}'.format(wp)\
      )
    #endfor k in range(npeaks)
  #endif not isilent

  return fwhm

#enddef vfwhm(x,y,nsmooth=0,isilent=0)

def vfwhmtot(x='?',y='',nsmooth=0,isilent=0):

  """ Referes to only one total maximum, not to peaks like vfwhm(...) """

# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug


  reset_status()

  if type(x) == str and x == '?':
    print("\nUsage: fwhm = vfwhm(x,y,nsmooth=0, isilent=1")

  ymax = y.max()
  k = -1
  for i in range(len(y)):
    if y[i] == ymax:
      k = i
      break
    #endif y[i] == ymax
  #endfor i in range(len(y))

  npoi = len(x)

  im = -9999
  ip = -9999

  half = ymax / 2.

  for l in range(0,k+1):
    if y[l] >= half and l > 0:
      im = l - 1
      break
    #endif y[l] >= half
  #endfor l in range(0,k+1)

  for l in range(k,npoi):
    if y[l] <= half and l < npoi:
      ip = l
      break
    #endif y[l] >= half
  #endfor l in range(k,npoi)

  xm = -9999.
  xp = -9999.
  x1m = -9999.
  x1p = -9999.
  x2m = -9999.
  x2p = -9999.
  y1m = -9999.
  y1p = -9999.
  y2m = -9999.
  y2p = -9999.
  w = -9999.

  if im != -9999 and y[im] <= half:
    x1m = x[im]
    x2m = x[im+1]
    y1m = y[im]
    y2m = y[im+1]
    xm = x1m + (x2m-x1m)/(y2m-y1m) * (half-y1m)
  #endif im != -9999 and ip != -9999

  if ip != -9999 and y[ip]  < half:
    x1p = x[ip-1]
    x2p = x[ip]
    y1p = y[ip-1]
    y2p = y[ip]
    xp = x1p + (x2p-x1p)/(y2p-y1p) * (half-y1p)
  #endif im != -9999 and ip != -9999

  if xm != -9999. and xp != -9999: w = xp - xm

         # 0 1   2   3   4   5   6   7   8  9    10  11
  fwhm = ([w,xm,x[k],xp,x1m,y1m,x2m,y2m,x1p,y1p,x2p,y2p])

  if not isilent:
    print("\n       FHWH        Xm       Xpeak        Xp        Wm       Wp")
    print("--------------------------------------------------------------------")

    wm = -9999.
    wp = -9999.

    if fwhm[3] != -9999. and fwhm[2] != -9999.:
      wm = fwhm[3]-fwhm[2]
    if fwhm[3] != -9999. and fwhm[4] != -9999.:
      wp = fwhm[4]-fwhm[3]

    print(" "\
    '     {:.4g}'.format(fwhm[0]),\
    '     {:.4g}'.format(fwhm[1]),\
    '     {:.4g}'.format(fwhm[2]),\
    '     {:.4g}'.format(fwhm[3]),\
    '     {:.4g}'.format(wm),\
    '     {:.4g}'.format(wp)\
    )
  #endif not isilent

  return fwhm

#enddef vfwhmtot(x,y,nsmooth=0,isilent=0)

def vfwhmtotabs(x='?',y='',nsmooth=0,isilent=0):
    return vfwhmtot(x,abs(y),nsmooth,isilent)
#enddef vfwhmtotabs(x='?',y='',nsmooth=0,isilent=0)

def getzone(projection=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  Fig = plt.gcf()
  Axes = Fig.get_axes()

  if not len(Axes):
    Ax = plt.gca()
    Axes = Fig.get_axes()
    Nwins = 1
    Nxzone = 1
    Nyzone = 1
    Kzone = 1
    return
  #endif not len(Fig.get_axes())

  if Nwins <= 0:
    window(projection=projection)
    zone(projection=projection)
    return
  #endif Nwins <= 0

  Axes = Fig.axes
  Ax = plt.gca()

  ifound = 0

  for z in Zones:
    if z[0] == Fig:
      Nxzone = z[1]
      Nyzone = z[2]
      Kzone = z[3]
      ifound = 1
      break
  #endfor z in Zones:

  if not ifound: print("*** Error in getzone(...): No zones found for this window ***")

  isame = Isame

  if type(IsameGlobal) == int:
    if IsameGlobal == 1: isame += 1
  else:
    if IsameGlobal.get() == 1: isame += 1
  #endif type(IsameGlobal) == int

  if isame: return

  try:
    if Tdate:
      Tdate.remove()
      Tdate = None
  except:
    Tdate = None
  #endtry

  i = 0; iax = -1
  for ax in Axes:
    if ax == Ax:
      iax=i
      break
    #endif ax == Ax
    i += 1
  #endfor ax in Axes

  if iax >= 0 and iax < len(Axes)-1 and  hasattr(Axes[iax+1],'get_label'):
    Fig.delaxes(Axes[iax+1])

  Fig.delaxes(Ax)
  Legend = []

  label = str(Nxzone) + "_" + str(Nyzone) + "_" + str(Kzone)

  set_plot_params()

  if len(projection):
    if projection.lower() == '3d': set_plot_params_3d()
    Ax = Fig.add_subplot(Nyzone,Nxzone,Kzone,projection=projection,
                        label = label)
  else:
    Ax = Fig.add_subplot(Nyzone,Nxzone,Kzone,label=label)
  #endif len(projection)

  Axes = Fig.axes
  if Kdate: date_on_figure()

  return
#enddef getzone()

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  Console = console

  if platform.system() == 'Linux':
    sys.stdout.write("\x1b]2;" + console + "\x07")
  elif platform.system() == 'Windows':
    #ctypes.windll.kernel32.SetConsoleTitleW(console)
    os.system("title "+console)
  #endif

#enddef set_console_title()

def get_console(console=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if console == '': console = Console

  if platform.system() == 'Linux':

    if FirstConsole:
      com = "(sleep 3; wmctrl -a " + console +str(") &")
      FirstConsole = 0
    else:
      com = "(sleep 3; wmctrl -a " + console +str(") &")
    #endif

    stat = os.system(com)
    if stat: print('... Could not raise console, wmctrl is not installed ...')
  #endif

#enddef get_console()

def getax(visible=True):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  Fig = plt.gcf()
  ax = plt.gca()
  ax.set_visible(visible)
  Ax = ax
  Tax = type(Ax)
  return ax
#def getax():

def vplbxy(x,y,u,v,scale=-9999.0,plopt='',tit='',xtit='',ytit='',ztit='',label='',
           color='default'):

#There are many options more..., not yet implemented here

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if Nwins <=0: window()

  plotopt(plopt)
  if not Isame: getzone()

  if color == 'default':
    lcol = Linecolor
  else:
    lcol = color
  #endif color == 'default'

  if scale == -9999:
    Ax.quiver(x, y, u, v, color=lcol)
  else:
    Ax.quiver(x, y, u, v, scale=1./scale, color=lcol)
  #endif

  txyz(tit,xtit,ytit)
  showplot()

#enddef vplbxy

def nplbxy(nt='?',varlis='x:y:bx:by',select='',scale=-9999.0,plopt='',
            tit='',xtit='',ytit='',label='',color='default'):
  x,y,bx,by = ncopv(nt,varlis,select)
  vplbxy(x,y,bx,by,scale,plopt,tit,xtit,ytit,label,color)
#enddef nplbxyz

def nplbxyz(nt='?',varlis='x:y:z:bx:by:bz',select='',scale=1.,plopt='',
            tit='',xtit='',ytit='',ztit='',label='',color='default'):
  x,y,z,bx,by,bz = ncopv(nt,varlis,select)
  vplbxyz(x,y,z,bx,by,bz,scale,plopt,tit,xtit,ytit,ztit,label,color)
#enddef nplbxyz

def vplbxyz(x,y,z,u,v,w,scale,plopt='',tit='',xtit='',ytit='',ztit='',label='',
           color='default'):

#There are many options more..., not yet implemented here

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if Nwins <=0: window()

  plotopt(plopt)
  if not Isame: getzone(projection='3d')

  if color == 'default':
    lcol = Linecolor
  else:
    lcol = color
  #endif color == 'default'

  Ax.quiver(x, y, z, u*scale, v*scale, w*scale, color=lcol)

  txyz(tit,xtit,ytit,ztit)
  showplot()

#enddef vplbxyz

def vplxyz(x,y,z,plopt='',tit='',xtit='',ytit='',ztit='',label='',
           color='default',cmap='', surfcolor=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  #if label: Klegend.set(1)

  if Nwins <=0: window()

  plotopt(plopt)
  if not Isame: getzone(projection='3d')

  iplot = 0

  if Iscatter:
    img = Ax.scatter(x,y,z,marker=Markertype,c=z,label=label)
    if Colorbarpad != '!':
      cbar = Fig.colorbar(img, pad=Colorbarpad)
    else:
      cbar = Fig.colorbar(img)
    #endif
    cbar.set_label(label=ztit, labelpad=Axistitledist, size=Axislabelsize)
    cbar.ax.tick_params(labelsize=Axislabelsize)
    Axes.append(cbar)
    txyz(tit,xtit,ytit)
    cbar.set_label(ztit, rotation=90)
    showplot()

    return
  #endif Iscatter

  if Iclosed:

    if type(x) == 'list':
      x.append(x[0])
      y.append(y[0])
      z.append(z[0])
    else:
      l = len(x)
      x[l] = x[0]
      y[l] = y[0]
      z[l] = z[0]
    #endif

    Iline = 1
  #endif

  if Iline:

    if color == 'default': color = Linecolor
    plt.plot(x,y,z,c=color,ls=Linestyle,lw=Linewidth,label=label)
    txyz(tit,xtit,ytit,ztit)
    return
    iplot=1

  if Itrisurf:

    if cmap == '' or cmap == '!': cmap=Cmap

    if surfcolor == '':
      cmap = None
      surfcolor = Surfcolor

    if cmap == None:
      Ax.plot_trisurf(x,y,z,shade=True,color=surfcolor,label=label)
    else:
      Ax.plot_trisurf(x,y,z,shade=True,cmap=cmap,label=label)
    #endif

    iplot=1

  if Iscat3d:
    Ax.scatter(x,y,z,marker=Markertype,c=Markercolor,label=label)
    iplot = 1
  #endif

  if not iplot:
    Ax.scatter(x,y,z,marker=Markertype,c=Markercolor,label=label)
  #endif iplot == 0:

  txyz(tit,xtit,ytit,ztit)
  showplot()
#enddef vplxyz

def vplxyzt(x,y,z,t,plopt='',tit='',xtit='',ytit='',ztit='', label='',
           color='default',cmap=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if Nwins <=0: window()
  #if label: Klegend.set(1)

  plotopt(plopt)
  getzone(projection='3d')

  iplot = 0

  if cmap == '' or cmap == '!': cmap=Cmap

  if Iscat3d:

    Ax.scatter(x,y,z,cmap=cmap,c=t,marker=Markertype,fillstyle=Fillstyle,label=label)
    iplot = 1

  if not iplot:
    Ax.scatter(x,y,z,cmap=cmap,c=t,marker=Markertype,fillstyle=Fillstyle,label=label)
  #endif iplot == 0:

  txyz(tit,xtit,ytit,ztit)

  showplot()
#enddef vplxyzt

def textbox(text,x=0.05, y=0.95, tcolor=None, bgcolor='white', alpha=0.9,
           rotation='horizontal'):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

  props = dict(facecolor=bgcolor, alpha=alpha)
  if type(Ax) == Tax3d:
    Ax.text2D(x,y,text,transform=Ax.transAxes,
            verticalalignment='top', bbox=props, rotation=rotation)
  else:
    Ax.text(x,y,text,transform=Ax.transAxes,
            verticalalignment='top', bbox=props, rotation=rotation)
#enddef textbox()

def shpl():
  plt.show(block=False)
#def shpl()

def seed(iseed=0, plopt=''):

  from pickle import dump
  import numpy as np

  if iseed >= 0: np.random.seed(iseed)

  if plopt == 'w':
    with open('m_hbook.seed', 'wb') as f:
      dump(np.random.get_state(), f)
  elif plopt == 'r':
    with open('m_hbook.seed', 'rb') as f:
      np.random.set_state(load(f))
  #endif plopt='w'

#def seed(iseed)

#----------------------------------------------------------------------

def vfitpoly(nord,x,y, ey='', cov='default', isilent=0, ninter=101, iretval=1,
             kweedzero=1):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  # weights are applied to y. For Gaussian uncertenties use 1/sigma not
  # 1/sigma**2

  if len(ey):
    iey = 1
    weight = 1.0/ey
    if cov == 'default': cov = 'unscaled'
  else:
    iey = 0
    weight = None
    if cov == 'default': cov = 'scaled'
  #endif len(ey)

  npar = nord + 1
  ndat = len(x)

  if (ndat == npar):
    # Fuer ndat = npar funktioniert der Fit nicht
    for i in range(ndat):
      x[ndat+i] = x[i]
      y[ndat+i] = y[i]
      if len(ey):
        ey[ndat+i] = ey[i]
    #endfor
    if len(ey):
      weight = 1.0/ey
  #endif

  par,covmat = np.polyfit(x, y, nord, w=weight, cov=cov)
  sigma = [covmat[i,i]**0.5 for i in range(npar)]
  f = np.polyval(par,x[:ndat])

  if not isilent:
    print("\n")
    for i in range(nord+1):
      k = nord - i
      pari = "P" + str(k) + ':'
      print(pari,par[i]," +/- ",sigma[i])
    #endfor
  #endif not isilent

  chi2 = 0.0
  chi2ndf = 0.0
  FitChi2Prob = 0.0

  ndf = ndat-npar-1

  if not iey: ey = 1.

  if ndat > npar:
    ch2 = ((y-f)/ey)**2
    chi2=ch2.sum()
    chi2ndf = chi2/ndf
    FitChi2Prob = chi2ndf_prob(chi2ndf,ndf)
  #endif ndat > npar:

  print("\nChi2, NDF, Chi2/NDF, Chi2Prob:",chi2,ndf,chi2ndf,FitChi2Prob)

  if cov != 'unscaled':
    print("\nSince cov != 'unscaled', sigmas are scaled such that chi2/ndf = 1\n")

  FitPar = list(par)
  FitSig = sigma

  Fitp=open("vfit.par","w")
  for i in range(nord+1): Fitp.write(str(par[i]) + "  " + str(sigma[i]) + "\n")
  Fitp.close()

  xmin = x.min()
  xmax = x.max()
  dx = (xmax-xmin)/(ninter-1)
  xx = xmin-dx

  Ffit = open("vfit.dat","w")
  for i in range(ninter):
    ff = np.polyval(par,xx)
    Ffit.write(str(xx) + " " + str(ff) + " \n")
  #endfor
  Ffit.close()

  if not isilent:
    print("\nFit parameters and sigmas are stored in FitPar and FitSig and also\n")
    print("written to vfit.par")
    print("The interpolated fit curve is written to vfit.dat")
  #endif

  if iretval: return par,sigma,chi2ndf,f
  else: return

#enddef vfitpoly(x,y,ey)

def vfitp1(x,y, ey='', cov='default', isilent=0, ninter=101, iretval=1,
           kweedzero=1):
  par,sigma,chi2ndf,f = vfitpoly(1,x,y, ey, cov, isilent, ninter, iretval,kweedzero)
  if iretval: return par,sigma,chi2ndf,f
#enddef

def gauss(x, A=0.3989422804014, mu=0.0, sig=1.0):
  import numpy as np
  return A * np.exp(-(((x-mu)/sig)**2)/2.)
#def gauss(x, mu=0.0, sig=1.0)

def fcos(x, A=1., x0=0.0, L=6.283185307179586):
  import numpy as np
  k = 2.*pi/L
  return A * cos(k*(x-x0))

def fcosh(x, A=1., x0=0.0, L=6.283185307179586):
  import numpy as np
  k = 2.*pi/L
  return A * cosh(k*(x-x0))

def b_to_K(bv='?',lam=None,bh=0.0):
  global \
  clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
  emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
  grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
  radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
  ecdipev,ecdipkev

  #Bh, Bv in Tesla, lam in mm
  if type(bv) == str:
    print("\nUsage: b_to_K(Bv/T, lamba/mm [, Bh=0.0]")
    return

  if type(bv) == list:
    K = []
    for i in range(len(bv)):
      bvv = bv[i]
      if type(bh) == list: bhh = bh[i]
      else: bhh = bh
      K.append(echarge1*(bhh**2+bvv**2)**0.5*lam/1000.0/(2.*pi1*emasskg1*clight1))
  else:
    K=echarge1*(bh**2+bv**2)**0.5*lam/1000.0/(2.*pi1*emasskg1*clight1)
  return K

def K_to_b(K='?',lam=None):
  global \
  clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
  emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
  grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
  radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
  ecdipev,ecdipkev

  #B in Tesla, lam in mm
  if type(K) == str:
    print("\nUsage: K_to_b(K, lamba/mm B,h=0.0)")
    return
  if type(K) == list:
    b = []
    for k in K:
      b.append(k/(echarge1*lam/1000.0/(2.*pi1*emasskg1*clight1)))
  else:
    b=K/(echarge1*lam/1000.0/(2.*pi1*emasskg1*clight1))
  return b

def K_to_harm(K='?',lam=None,ebeam=None):
  global \
  clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
  emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
  grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
  radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
  ecdipev,ecdipkev

  if type(K) == str:
    print("\nUsage: K_to_harm(K, lamba/mm, Ebeam/GeV")
    return
  gam = ebeam/emassg1
  if type(K) == list:
    harm = []
    for k in K:
      wlen1=(1.+k**2/2.)/2./gam**2*lam/1000.0e0*1.e9
      harm.append(wtoe1/wlen1)
  else:
    wlen1=(1.+K**2/2.)/2./gam**2*lam/1000.0e0*1.e9
    harm = wtoe1/wlen1
  return harm

def b_to_harm(b='?',lam=None,ebeam=None):
  global \
  clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
  emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
  grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
  radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
  ecdipev,ecdipkev

  if type(b) == str:
    print("\nUsage: b_to_harm(B/T, lamba/mm, Ebeam/GeV")
    return
  K = b_to_K(b,lam)
  return K_to_harm(K,lam,ebeam)

def harm_to_K(ebeam='?',lam=None,nharm=None,harm=None):
  global \
  clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
  emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
  grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
  radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
  ecdipev,ecdipkev


  if type(ebeam) == str:
    print("\nUsage: harm_to_K(Ebeam/GeV, lamba/mm, Nharm, Harmonic /eV")
    return

  h = harm/nharm
  rlamb1=wtoe1/h*1.0e-6
  gamma=ebeam/emassg1
  K=2.0*(rlamb1/lam*2.0*gamma**2-1.0)

  if K > 0.0 :
    K=K**0.5
  else:
    K=0.0

  return K, K_to_b(K,lam)

def expo(x, a=1.0, lam=1.0):
  import numpy as np
  return a * np.exp(lam*x)
#enddef expo(x, a=1.0, lam=1.0)

def expo2(x, a=1.0, lam=1.0, lam2=0.0):
  import numpy as np
  return a * np.exp(x*(lam+lam2*x))
#enddef expo2(x, a=1.0, lam=1.0, lam2=0.0)

def vgauss(x, mu=0.0, sig=1.0):
  import numpy as np
  return 3.989422804014e-1/sig * np.exp(-(((x-mu)/sig)**2)/2.)
#def vgauss(x, mu=0.0, sig=1.0)

def getfitpar():
  global FitPar
  return FitPar
#enddef getfitpar()

def getfitndf():
  global FitNdf
  return FitNdf
#enddef getfitndf

def getfitchi2ndf():
  global FitChi2ndf
  return FitChi2ndf
#enddef getfitchi2ndf

def getfitchi2prob():
  global FitChi2Prob
  return FitChi2Prob
#enddef getfitchi2ndf

def getfitsig():
  global FitSig
  return FitSig
#enddef getfitsig()

def hfit(idh, fitfun, select='',absolute_sigma='default', parstart=None,
         bounds=None, method=None,isilent=0, ninter=101, iretval=1,
         iplot=1,fitcol='!',kweedzero=1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  # Histograms and Ntuples
  global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
  Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
  H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
  Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
  Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
  H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
  FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman

#  global Kold, Kstat, Mode2d, Nfitxy, Kplots, Kzone,

  nt = hcopn(idh,'nHfit',kweedzero=1)
  if H1I < 0: return None,None,None,None

  par,sigma,chi2ndf,f = nfitxy(nt,'x:y:ey',select,fitfun, absolute_sigma, \
  parstart,bounds,method,isilent, ninter,1,0,fitcol,kweedzero)

  igauss = 0
  if type(fitfun) == str and fitfun == 'gauss': igauss = 1

  if iplot:

    if fitcol == '!': fitcol = 'black'

    Kold = Kstat
    if Kstat: Kstat = False

    if Nfitxy.ey.abs().max() > 0: hplot1d(idh,'e')
    else: hplot1d(idh,'hist')

#      plt.errorbar(Nfitxy.x,Nfitxy.f,Nfitxy.ey, \
#    ls='',marker=Markertype, fillstyle=Fillstyle, mfc=fitcol, \
#    mec=fitcol, ms=Markersize, mew=1, c=fitcol)

    vplxy(Nfitxy.x,Nfitxy.f,"samespline",color=fitcol)
    Kstat = Kold

    if Kstat:
      if StatFontSize < 0:
        dpi = Fig.dpi
        nxy = max(Nxzone,Nyzone)
        sfs = 12 - nxy * 2
      else:
        sfs = StatFontSize
      #endif StatFontSize < 0

      sig = FitSig

      if igauss:
        tex = "A = " + '{:.4g}'.format(par[0]) + " +/- " + '{:.4g}'.format(sig[0]) + \
        "\n$\mu$ = " + '{:.4g}'.format(par[1])  + " +/- "  + '{:.4g}'.format(sig[1]) + \
        "\n$\sigma$ = " + '{:.4g}'.format(par[2])  + " +/- " + '{:.4g}'.format(sig[2]) + "\n"
      else:
        tex = ""
        ip = 0
        for p in par:
          tex += "P" + str(ip) + " = " + '{:.4g}'.format(p)  + '{:.4g}'.format(sig[ip]) + "\n"
          ip += 1
        #endfor
      #endif

      tex += "$\chi^2/Ndf$" + " = " + '{:.4g}'.format(chi2ndf) + "\n"
      tex += "$\chi^2 prob$" + " = " + '{:.4g}'.format(FitChi2Prob) + "\n"
      text(Xfit,Yfit,tex,halign='left')

    #endif Kstat
  #endif iplot

  if iretval: return par,sigma,chi2ndf,f
  else: return

#enddef

def hfitg(idh, select='',absolute_sigma='default', parstart=None, bounds=None,
          method=None,isilent=0, ninter=101, iretval=1, iplot=1, fitcol='!',
          kweedzero=1):

  par,sigma,chi2ndf,f = hfit(idh,'gauss', select,absolute_sigma, parstart, \
  bounds,method,isilent,ninter,1,iplot,fitcol,kweedzero)

  if iretval: return par,sigma,chi2ndf,f
  else: return

#enddef

def vfit(fitfun, x, y, ey = '', absolute_sigma='default', parstart=None,
         bounds=None, method=None,isilent=0, ninter=101, iretval=1,
         kweedzero=1):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  # Histograms and Ntuples
  global H1h, H1hh, H2h, H2hh, H1, H2, H1head, H2head, H1HLast, Nhead, Ntup, \
  Nctup, Nh1, Nh2, Nntup, Nnctup, Hdir, Ndir, Kdir, Cdir, Fdir, \
  H1Last, H2Last, NLast, H1h, H2h, N, Nct, Ind, IndLast, \
  Nmin, Nmax, Nmean, Nrms, Nxopt, Nyopt, Nlook, \
  Tdf, Tfig, Tax, Tax3d, Tax2d , H1ind, H2ind, Ncind, \
  H1ILast, NiLast, H1I, H2I, H2ILast, Ni, NctI, Nind, Nsel, Nlines, Ncolon, \
  FitPar, FitFit, FitSig, FitChi2ndf, FitNdf, FitChi2Prob,Figman

  from scipy.optimize import curve_fit

  if not bounds: bounds = (-np.inf,np.inf)

  iey = 1

  if not len(ey) or ey.abs().max() == 0.0:
    iey=0
    ey = y
    ef = None
    if absolute_sigma == 'default':
      absolute_sigma = False
  else:
    if absolute_sigma == 'default':
      absolute_sigma = True
  #endif not len(ey)

  if kweedzero:
    xf = []
    yf = []
    if iey: ef = []
    for i in range(len(x)):
      if y[i] !=0 and ey[i] != 0:
        xf.append(x[i])
        yf.append(y[i])
        if iey: ef.append(ey[i])
      #endif
    #endfor
    xf = pdf(xf,'x').x
    yf = pdf(yf,'y').y
    if iey: ef = pdf(ef,'ey').ey
    else: ey = ey * 0.0
  else:
    xf = x
    yf = y
    if iey: ef = ey
  #endif

  par, covmat = curve_fit(fitfun,xf,yf,p0=parstart,sigma=ef, \
  absolute_sigma=absolute_sigma)

  ndat = len(xf)
  npar = len(par)

  sigma = np.sqrt(np.diag(covmat))
  f = fitfun(x,*par) # Example: fitfun(x,p1,p2)

  chi2ndf = 0.0
  chi2 = 0.0
  ndf = ndat - npar - 1

  if not iey: ef = 1.

  if ndat > npar:
    ch2 = ((yf-f)/ef)**2
    ch2[np.isinf(ch2)] = 0.0
    chi2=ch2.sum()
    chi2ndf = chi2/ndf
  #endif ndat > npar:

  FitPar = par
  FitSig = sigma
  FitChi2ndf = chi2ndf
  FitNdf = ndf
  FitFit = f
  FitChi2Prob = chi2ndf_prob(chi2ndf,ndf)

  if not isilent:
    print("\n")
    for i in range(npar):
      pari = "P" + str(i) + ':'
      print(pari,par[i]," +/- ",sigma[i])
    #endif
    print("Chi2, NDF, Chi2/NDF, Chi2Prob: ",chi2, ndf, chi2ndf,FitChi2Prob,NL)
  #endif not isilent

  if not absolute_sigma and not isilent:
      print("Since absolute_simga=True, sigmas are scaled such that chi2/ndf = 1")
      print("However, the return value referes to the unscaled ey!\n")
  #endif not absolute_sigma and not isilent

  if not isilent:
    print("Return values: par, sigma, chi2ndf, fit (if iretval is set)")
    print("Also stored in FitPar, FitSig, FitChi2ndf, FitNdf, FitChi2Prob, FitFit\n")
  #endif not isilent

  if ninter < 1: ninter = 2

  Fres = open("fit.res","w")
  for i in range(ndat):
    xx = xf[i]
    ff = fitfun(xx,*par)
    Fres.write(str(xx) + " " + str(y[i]) + " " +str(ff)+"\n")
  #endfor
  Fres.close()

  xmin = x.min()
  xmax = x.max()
  dx = (xmax-xmin)/(ninter-1)
  xx = xmin-dx

  Ffit = open("fit.int","w")
  for i in range(ninter):
    xx += dx
    ff = fitfun(xx,*par)
    Ffit.write(str(xx) + " " + str(ff)+"\n")
  #endfor
  Ffit.close()

  Fpar = open("fit.par","w")
  for i in range(npar):
    Fpar.write(str(par[i]) + " " + str(sigma[i]) +"\n")
  #endfor
  Fpar.close()

  if not isilent:
    print("\nFit parameters written to fit.par.")
    print("Interpolated fit curve written to fit.int.")
    print("fit data points written to fit.res.")
    print("\nNtuple Nfitxy contains x:y:ey:f data.")
    print("To get Parameters etc., use getfitpar(), getfitsig(), and getchi2ndf(), getchi2prob()")
  #endif not isilent

  Nfitxy = ncre("Nfitxy","x:y:ey:f",ioverwrite=1)

  Nfitxy.x = x
  Nfitxy.y = y
  Nfitxy.ey = ey
  Nfitxy.f = f

  nupdate_header(Nfitxy)

  if iretval: return par,sigma,chi2ndf, f
  else: return

#enddef vfit(fitfun,x,y, ey = None, parstart=None, bounds=None, method=None,isilent=0)

def vfitexp(x,y, ey = '', absolute_sigma='default', parstart=None,
            bounds=None, method=None,isilent=0,ninter=101, iretval=1,
            kweedzero=1):
  # See also vfit(...)
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  from scipy.optimize import curve_fit

  par, sigma, chi2ndf, f = vfit(expo,x,y,ey,'default',parstart, bounds,
                                   method,isilent,ninter,kweedzero)

  if not isilent: print("\nvfitexp: fit = par[0] * exp(par[1]*x)\n")
  return  par, sigma, chi2ndf, f

#enddef vfitexp(x,y, ey = '', absolute_sigma='default', parstart=None, bounds=None, method=None,isilent=0)

def vfitexp2(x,y, ey = '', absolute_sigma='default', parstart=None,
             bounds=None, method=None,isilent=0,ninter=101, iretval=1,
             kweedzero=1):
  # See also vfit(...)
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  from scipy.optimize import curve_fit

  if parstart == None: parstart = [1.,1.,0]

  par, sigma, chi2ndf, f = vfit(expo,x,y,ey,absolute_sigma,parstart[:2],bounds,
                                method,1,ninter,iretval,kweedzero)
  parstart[0:2] = par[0:2]

  par, sigma, chi2ndf, f = vfit(expo2,x,y,ey,absolute_sigma,parstart,bounds,
                                method,1,ninter,iretval,kweedzero)

  if not isilent: print("\nvfitexp2: fit = par[0] * exp(x*(par[1]+par[2]*x)\n")
  return  par, sigma, chi2ndf, f

#enddef vfitexp2(x,y, ey = '', absolute_sigma='default', parstart=None, bounds=None, method=None,isilent=0)

def vfitgauss(x,y, ey = '', absolute_sigma='default',
              parstart=None, bounds=None, method=None,isilent=0,ninter=101,
              kweedzero=1):
  # See also vfit(...)
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  from scipy.optimize import curve_fit

  par, sigma, chi2ndf, f = vfit(gauss,x,y,ey,'default',parstart, bounds,
                                   method,isilent,ninter,kweedzero)
  return  par, sigma, chi2ndf, f

#enddef vfitgauss(x,y, ey = '', absolute_sigma='default', parstart=None, bounds=None, method=None,isilent=0)

def vfitcosh(x,y, ey = '', absolute_sigma='default',
             parstart=None, bounds=None, method=None,isilent=0,ninter=101,
             kweedzero=1):
  # See also vfit(...)
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  from scipy.optimize import curve_fit

  par, sigma, chi2ndf, f = vfit(fcosh,x,y,ey,'default',parstart, bounds,
                                   method,isilent,ninter,kweedzero)
  return  par, sigma, chi2ndf, f

#enddef vfitcosh(x,y, ey = '', absolute_sigma='default', parstart=None, bounds=None, method=None,isilent=0)

def vfitcos(x,y, ey = '', absolute_sigma='default',
            parstart=None, bounds=None, method=None,isilent=0,ninter=101,
            kweedzero=1):
  # See also vfit(...)
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  from scipy.optimize import curve_fit

  par, sigma, chi2ndf, f = vfit(fcos,x,y,ey,'default',parstart, bounds,
                                   method,isilent,ninter,kweedzero)
  return  par, sigma, chi2ndf, f

#enddef vfitcos(x,y, ey = '', absolute_sigma='default', parstart=None, bounds=None, method=None,isilent=0)

def hget(idh=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz

  if type(idh) == int:
    return H1[idh]
  else:
    idx = GetIndexH2(idh)
    if idx != -1:
      return H2[idx]
    else:
      idx = GetIndexH1(idh)
      if idx != -1:
        return H1[idx]
      else:
        ErrorText = "*** Error in hget: Non-existing histogram " "idh"
        Istatus = -1
    return idx
#enddef hget(idh='h')

def nget(idn=''):

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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux

#+PATCH,//WAVES/PYTHON
#+KEEP,vecglobind,T=PYTHON.

  global VsortX, VsortY, VoptX, VoptY, VsplX, VsplY, Vspl1, Vspl2, VsplI, \
  VsplCoef, Nspline,Ninter, Nfitxy, Nfitint, Vxint, Vyint, SplineMode, \
  VxyzX,VxyzY,VxyzZ

#+KEEP,nxyzglobind,T=PYTHON.
#*CMZ :          29/09/2019  11.11.01  by  Michael Scheer
  global N1, N2, N3, N4, N5, N6, N7,N8,N9,Nv, Nx, Nxy, Nxyz


  if type(idn) == int:
    nupdate_header(Ntup[idn])
    return Ntup[idn]
  else:
    idx = GetIndexN(idn)
    if idx != -1:
      nupdate_header(Ntup[idx])
      return Ntup[idx]
    else:
        ErrorText = "*** Error in nget: Non-existing Ntuple " "idn"
        Istatus = -1
    return idx
#enddef nget(idn='')

def hadd(h1,h2,hres,tit=''):

  idx1 = GetIndexH1(h1)
  tit1 = H1hh[1]
  nx1 = H1hh[2]

  if idx1 == -1:
    print('*** Error in hadd: Non-existing histogram ',h1)
    return -1

  idx2 = GetIndexH1(h2)
  tit2 = H1hh[1]
  nx2 = H1hh[2]

  if tit == '': tit = tit1 + " x " + tit2

  if idx2 == -1:
    print('*** Error in hadd: Non-existing histogram ',h2)
    return -1

  if nx1 != nx2:
    print('*** Error in hadd: Number of channels nx1, nx2 do not match', nx1, nx2)
    return -1

  e1 = H1[idx1].ey
  e2 = H1[idx2].ey

  x = H1[idx1].x
  y = H1[idx1].y + H1[idx2].y
  ey = (e1**2 + e2**2)**0.5

  nt = make_dataframe('x:y',x,y)

  nt["x2"]=nt.x * nt.x
  nt["xy"]=nt.x * nt.y
  nt["x2y"]=nt.x * nt.x * nt.y

  sumx2 = nt.x2.sum()
  sumxy = nt.xy.sum()
  sumx2y = nt.x2y.sum()
  sumy = nt.y.sum()
  sumx = nt.x.sum()
  sumx2 = nt.x2.sum()

  xmean = 0.0
  x2mean = 0.0
  if sumy:
    xmean = sumxy / sumy
    x2mean = sumx2y / sumy
  #endif
  xrms = np.sqrt(max(0.0,x2mean-xmean**2))

  hret = hbook1(hres,tit,nx2,H1hh[3],H1hh[4])
  idxres = GetIndexH1(hres)

  H1head[idxres][7] = nt.y.min()
  H1head[idxres][8] = nt.y.max()
  H1head[idxres][9] = H1head[idx1][9] #entries
  H1head[idxres][10] = sumy
  H1head[idxres][11] = sumx2y
  H1head[idxres][12] = xmean
  H1head[idxres][13] = xrms
  H1head[idxres][18] = sumx
  H1head[idxres][19] = sumx2

  H1[idxres].n = H1[idx1].n
  H1[idxres].y = y
  H1[idxres].ave = y
  H1[idxres].ey = ey
  H1[idxres].y2 = ( ey**2 + y**2) * H1[idx1].n

  return idxres
#def hadd(h1,h2,hres)

def hsub(h1,h2,hres,tit=''):

  idx1 = GetIndexH1(h1)
  tit1 = H1hh[1]
  nx1 = H1hh[2]

  if idx1 == -1:
    print('*** Error in hadd: Non-existing histogram ',h1)
    return -1

  idx2 = GetIndexH1(h2)
  tit2 = H1hh[1]
  nx2 = H1hh[2]

  if tit == '': tit = tit1 + " x " + tit2

  if idx2 == -1:
    print('*** Error in hadd: Non-existing histogram ',h2)
    return -1

  if nx1 != nx2:
    print('*** Error in hadd: Number of channels nx1, nx2 do not match', nx1, nx2)
    return -1

  e1 = H1[idx1].ey
  e2 = H1[idx2].ey

  x = H1[idx1].x
  y = H1[idx1].y - H1[idx2].y
  ey = (e1**2 + e2**2)**0.5

  nt = make_dataframe('x:y',x,y)

  nt["x2"]=nt.x * nt.x
  nt["xy"]=nt.x * nt.y
  nt["x2y"]=nt.x * nt.x * nt.y

  sumx2 = nt.x2.sum()
  sumxy = nt.xy.sum()
  sumx2y = nt.x2y.sum()
  sumy = nt.y.sum()
  sumx = nt.x.sum()
  sumx2 = nt.x2.sum()

  xmean = 0.0
  x2mean = 0.0
  if sumy:
    xmean = sumxy / sumy
    x2mean = sumx2y / sumy
  xrms = np.sqrt(max(0.0,x2mean-xmean**2))

  hret = hbook1(hres,tit,nx2,H1hh[3],H1hh[4])
  idxres = GetIndexH1(hres)

  H1head[idxres][7] = nt.y.min()
  H1head[idxres][8] = nt.y.max()
  H1head[idxres][9] = H1head[idx1][9] #entries
  H1head[idxres][10] = sumy
  H1head[idxres][11] = sumx2y
  H1head[idxres][12] = xmean
  H1head[idxres][13] = xrms
  H1head[idxres][18] = sumx
  H1head[idxres][19] = sumx2

  H1[idxres].n = H1[idx1].n
  H1[idxres].y = y
  H1[idxres].ave = y
  H1[idxres].ey = ey
  H1[idxres].y2 = ( ey**2 + y**2) * H1[idx1].n

  return idxres
#def hsub(h1,h2,hres)

def hdiv(h1,h2,hres,tit=''):

  idx1 = GetIndexH1(h1)
  tit1 = H1hh[1]
  nx1 = H1hh[2]

  if idx1 == -1:
    print('*** Error in hdiv: Non-existing histogram ',h1)
    return -1

  idx2 = GetIndexH1(h2)
  tit2 = H1hh[1]
  nx2 = H1hh[2]

  if tit == '': tit = tit1 + " x " + tit2

  if idx2 == -1:
    print('*** Error in hdiv: Non-existing histogram ',h2)
    return -1

  if nx1 != nx2:
    print('*** Error in hdiv: Number of channels nx1, nx2 do not match', nx1, nx2)
    return -1

  if H1[idx1].y != 0.0:
    e1r = abs(H1[idx1].ey / H1[idx1].y)
  else:
    e1r = 0.0

  if H1[idx2].y != 0.0:
    e2r = abs(H1[idx2].ey / H1[idx2].y)
  else:
    e2r = 0.0

  x = H1[idx1].x
  if H1[idx2].y == 0.0 and H1[idx1].y == 0.0:
    y = 0.0
  elif H1[idx2].y == 0.0:
    y = inf
  else:
    y = H1[idx1].y / H1[idx2].y

  ey = abs(y * (e1r**2 + e2r**2)**0.5)

  nt = make_dataframe('x:y',x,y)

  nt["x2"]=nt.x * nt.x
  nt["xy"]=nt.x * nt.y
  nt["x2y"]=nt.x * nt.x * nt.y

  sumx2 = nt.x2.sum()
  sumxy = nt.xy.sum()
  sumx2y = nt.x2y.sum()
  sumy = nt.y.sum()
  sumx = nt.x.sum()
  sumx2 = nt.x2.sum()

  xmean = 0.0
  x2mean = 0.0
  if sumy:
    xmean = sumxy / sumy
    x2mean = sumx2y / sumy
  xrms = np.sqrt(max(0.0,x2mean-xmean**2))

  hret = hbook1(hres,tit,nx2,H1hh[3],H1hh[4])
  idxres = GetIndexH1(hres)

  H1head[idxres][7] = nt.y.min()
  H1head[idxres][8] = nt.y.max()
  H1head[idxres][9] = H1head[idx1][9] #entries
  H1head[idxres][10] = sumy
  H1head[idxres][11] = sumx2y
  H1head[idxres][12] = xmean
  H1head[idxres][13] = xrms
  H1head[idxres][18] = sumx
  H1head[idxres][19] = sumx2

  H1[idxres].n = H1[idx1].n
  H1[idxres].y = y
  H1[idxres].ave = y
  H1[idxres].ey = ey
  H1[idxres].y2 = ( ey**2 + y**2) * H1[idx1].n

  return idxres
#def hdiv(h1,h2,hres)

def hmul(h1,h2,hres,tit=''):

  idx1 = GetIndexH1(h1)
  tit1 = H1hh[1]
  nx1 = H1hh[2]

  if idx1 == -1:
    print('*** Error in hmul: Non-existing histogram ',h1)
    return -1

  idx2 = GetIndexH1(h2)
  tit2 = H1hh[1]
  nx2 = H1hh[2]

  if tit == '': tit = tit1 + " x " + tit2

  if idx2 == -1:
    print('*** Error in hmul: Non-existing histogram ',h2)
    return -1

  if nx1 != nx2:
    print('*** Error in hmul: Number of channels nx1, nx2 do not match', nx1, nx2)
    return -1
  #endif

  e1r = abs(H1[idx1].ey / H1[idx1].y)
  e2r = abs(H1[idx2].ey / H1[idx2].y)

  x = H1[idx1].x
  y = H1[idx1].y * H1[idx2].y
  ey = abs(y * (e1r**2 + e2r**2)**0.5)

  nt = make_dataframe('x:y',x,y)

  nt["x2"]=nt.x * nt.x
  nt["xy"]=nt.x * nt.y
  nt["x2y"]=nt.x * nt.x * nt.y

  sumx2 = nt.x2.sum()
  sumxy = nt.xy.sum()
  sumx2y = nt.x2y.sum()
  sumy = nt.y.sum()
  sumx = nt.x.sum()
  sumx2 = nt.x2.sum()

  xmean = 0.0
  x2mean = 0.0
  if sumy:
    xmean = sumxy / sumy
    x2mean = sumx2y / sumy
  xrms = np.sqrt(max(0.0,x2mean-xmean**2))

  hret = hbook1(hres,tit,nx2,H1hh[3],H1hh[4])
  idxres = GetIndexH1(hres)

  H1head[idxres][7] = nt.y.min()
  H1head[idxres][8] = nt.y.max()
  H1head[idxres][9] = H1head[idx1][9] #entries
  H1head[idxres][10] = sumy
  H1head[idxres][11] = sumx2y
  H1head[idxres][12] = xmean
  H1head[idxres][13] = xrms
  H1head[idxres][18] = sumx
  H1head[idxres][19] = sumx2

  H1[idxres].n = H1[idx1].n
  H1[idxres].y = y
  H1[idxres].ave = y
  H1[idxres].ey = ey
  H1[idxres].y2 = ( ey**2 + y**2) * H1[idx1].n

  return idxres
#def hmul(h1,h2,hres)

def df_to_ntup(df,nt,tit='!'):

  if type(df) != Tdf:
    print("*** Error in df_to_ntupdf must be of type DataFrame ***")

  idn = GetIndexN(nt,1)

  if idn != -1:
    print("*** Error in df_to_ntupdf: Ntuple dose already exist ***")
    return -1

  if type(nt) != str: nt = "N_" + str(nt)
  if tit == '!': tit = nt


  cols = df.columns
  nvar = len(cols)

  varlis = cols[0]
  for v in cols[1:]: varlis += ':'+ v

  nt = ncre(nt,tit,varlis)

  if type(nt) == int and nt == -1:
    print("*** Error in df_to_ntupdf: Bad return from ncre ***")
    return -1

  idn = GetIndexN(nt)

  Ntup[idn] = pd.concat([nt,df])
  nt = Ntup[idn]

  nhead = Nhead[idn]
  for i in range(nvar):
    nhead[4+i][1] = nt[cols[i]].min()
    nhead[4+i][2] = nt[cols[i]].max()

  nhead[4+nvar] = len(df)

  return idn
#def dftontup(df,nt)

def h1copv(h, varlis='x:y', sel=''):
# +PATCH,//WAVES/PYTHON
# +KEEP,statusglobind,T=PYTHON.
  global Istatus, WarningText, ErrorText, Gdebug


  reset_status()

  if type(h) == str:
    idx = GetIndexH1(h)
    h = H1h
  #endif type(nt) == str:

  vxy = make_empty_dataframe(varlis)
  varxy = varlis.split(':')

  vxy[varxy[0]] = h.x
  vxy[varxy[1]] = h.y

  if len(sel):
    return  vxy.query(sel)
  else: return vxy

def vlocate(command='new'):

  global Vlocate, CanVlocate, Fig, \
  VlocX,VlocY, VlocDistX,VlocDistY, VlocDist, \
  VlocXO,VlocYO, VlocDistXO,VlocDistYO, VlocDistO

  print("--- Get values with left mouse click, delete with right click ---")
  print("--- Terminate with double click ---")
  print("--- Values are in list Vlocate  ---\n")

  print("--- n or new clears list ---")
  print("--- d or deletes last item ---")

  if command == 'new' or command == 'n':
    Vlocate = []
    VlocX = 0
    VlocX = 0
    VlocDistX = 0
    VlocDistY = 0
    VlocDist = 0
  elif command == 'delete' or command == 'd':
    print(Vlocate.pop(-1)," deleted")
  #endif command

  CanVlocate = Fig.canvas.mpl_connect('button_press_event',_vlocate)

#enddef vlocate()

def setlinewidth(lw=1.5):
  global Linewidth
  mpl.rcParams['lines.linewidth'] = lw
  Linewidth = lw

def getlinewidth(): return mpl.rcParams.get('lines.linewidth')

def settextcolor(tc='black'):
  global Textcolor
  Textcolor = tc

def gettextcolor(): return Textcolor

def setlinecolor(lc='red'):
  global Linecolor
  mpl.rcParams['lines.color'] = lc
  Linecolor = lc

def getlinecolor(): return mpl.rcParams.get('lines.color')

def sethistcolor(hc='blue'):
  global HistEdgeColor
  HistEdgeColor=hc
  _sethistcolor()

def gethistcolor(): return HistEdgeColor

def shcolblack(): sethistcolor('black')
def shcolblue(): sethistcolor('blue')
def shcolcyan(): sethistcolor('cyan')
def shcolgray(): sethistcolor('gray')
def shcolgreen(): sethistcolor('green')
def shcolmagenta(): sethistcolor('magenta')
def shcolred(): sethistcolor('red')
def shcolwhite(): sethistcolor('white')
def shcolyellow(): sethistcolor('yellow')

def slcolblack(): setlinecolor('black')
def slcolblue(): setlinecolor('blue')
def slcolcyan(): setlinecolor('cyan')
def slcolgray(): setlinecolor('gray')
def slcolgreen(): setlinecolor('green')
def slcolmagenta(): setlinecolor('magenta')
def slcolred(): setlinecolor('red')
def slcolwhite(): setlinecolor('white')
def slcolyellow(): setlinecolor('yellow')

def setmarkercolor(mc='red'):
  global Markercolor
  Markercolor = mc
#enddef setmarkercolor(mc='red')

def getmarkercolor():
  global Markercolor
  return Markercolor
#def getmarkercolor()

def smcolblack(): setmarkercolor('black')
def smcolblue(): setmarkercolor('blue')
def smcolcyan(): setmarkercolor('cyan')
def smcolgray(): setmarkercolor('gray')
def smcolgreen(): setmarkercolor('green')
def smcolmagenta(): setmarkercolor('magenta')
def smcolred(): setmarkercolor('red')
def smcolwhite(): setmarkercolor('white')
def smcolyellow(): setmarkercolor('yellow')

smcol = setmarkercolor
smblack = smcolblack
smblue = smcolblue
smcyan = smcolcyan
smgray = smcolgray
smgreen = smcolgreen
smmagenta = smcolmagenta
smred = smcolred
smwhite = smcolwhite
smbyellow = smcolyellow

shcol = sethistcolor
shblack = shcolblack
shblue = shcolblue
shcyan = shcolcyan
shgray = shcolgray
shgreen = shcolgreen
shmagenta = shcolmagenta
shred = shcolred
shwhite = shcolwhite
shbyellow = shcolyellow

slcol = setlinecolor
slblack = slcolblack
slblue = slcolblue
slcyan = slcolcyan
slgray = slcolgray
slgreen = slcolgreen
slmagenta = slcolmagenta
slred = slcolred
slwhite = slcolwhite
slbyellow = slcolyellow

def setlinestyle(ls='solid'):
  global Linestyle
  mpl.rcParams['lines.linestyle'] = ls
  Linestyle = ls
#enddef setlinestyle(ls='solid')

def setmarkertype(ms='o'):
  global Markertype
  #mpl.rcParams['lines.marker'] = ms
  Markertype = ms
#def setmarkertype(ms='o')

def setfillstyle(fs='none'):
  global Fillstyle
  Fillstyle = fs
#def setfillstyle(fs='none')

#def dot(): setmarkertype('m00')
def dot(): setmarkersize(0.1)

def bull():
  setmarkertype('o')
  setfillstyle('full')
#enddef bull():
def circle():
  setmarkertype('o')
  setfillstyle('none')
#enddef bull():
def square():
  setmarkertype('s')
#enddef square()
def star():
  setmarkertype('*')
#enddef star():

def nextmarkertype():
  global Markertype
  marker = ['o','s','*','^','v','8','p','x','+']
  n = len(marker)
  for i in range(n):
    if Markertype == marker[i]:
      k = i + 1
      break
    #if Markertype == marker[i]
  #for i in len(marker)
  if k > n - 1: k=0
  Markertype = marker[k]
#enddef nextmarkertype()

def setmarkersize(ms=6.):
  global Markersize
  mpl.rcParams['lines.markersize'] = ms
  Markersize = ms
#enddef setmarkersize(ms=6.)

def getmarkersize(): return mpl.rcParams.get('lines.markersize')

def getaxislabelsize(): return Axislabelsize

def getaxistitlesize3d(): return Atitfontsize3d

def getaxistitledist(): return mpl.rcParams['axes.labelpad']

def getaxistitledist3d(): return mpl.rcParams['axes.labelpad']

def getaxislabeldist(): return Axislabeldist

def getaxislabeldist3d(): return AxisLabelDist3d

def setaxislabelsize(size):
  global AxisLabelSize
  AxisLabelSize = size

def setaxistitlesize3d(size):
  global AtitFontSize3d
  AtitFontSize3d = size

def setaxistitlesize(size):
  global AtitFontSize
  AtitFontSize = size

def setaxistitledist(dist):
  global AxisTitleDist
  AxisTitleDist = dist

def setaxistitledist3d(dist):
  global AxisTitleDist3d
  AxisTitleDist3d = dist

def setaxislabeldist(dist):
  global AxisLabelDist
  AxisLabelDist = dist

def setaxislabeldist3d(dist):
  global AxisLabelDist3d
  AxisLabelDist3d = dist

def showparams(): print(mpl.rc_params())

def optgrid(g=True):
    global Kgrid
    Kgrid = g
    plt.rcParams['axes.grid'] = g
    if not plt.get_fignums(): return
    plt.grid(g)
    showplot()
#enddef optgrid(g=True)

def grid(g=True):
    global Kgrid
    Kgrid = g
    plt.grid(g)
    showplot()
#enddef grid(g=True)

def nogrid(g=False):
    global Kgrid, Ax
    Kgrid = g
    if Ax:
      plt.grid(g)
      showplot()
    #endif
#enddef nogrid(g=False)

def getlegend():
  global Klegend
  return Klegend

def optlegend(ilege=1,pos='upper right'):
  global Klegend
  if type(Klegend) == int: Klegend=ilege
  else: Klegend.set(ilege)
  setlegendposition(pos=pos)

def setlegendposition(pos='upper right'):
  mpl.rcParams['legend.loc'] = pos

def gettitlepad():
  return mpl.rcParams['axes.titlepad']

def settitlepad(pad=6.0):
  mpl.rcParams['axes.titlepad'] = pad

def getlinestyleindex(lsty):
  global Linestyles
  idx = -1
  for i in range(len(Linestyles)):
    if Linestyles[i] == lsty:
      idx=i
      break
  #endfor i in range(len(Linestyles))
  return idx

def getmarkertypeindex(mty):
  global Markertypes
  idx = -1
  for i in range(len(Markertypes)):
    if Markertypes[i] == mty:
      idx=i
      break
  #endfor i in range(len(Markertypes))
  return idx

def getcolorindex(col):
  global Colors
  idx = -1
  for i in range(len(Colors)):
    if Colors[i] == col:
      idx=i
      break
  #endfor i in range(len(Colors))
  return idx

def getsurfcolorindex(col):
  global Surfcolors
  idx = -1
  for i in range(len(Surfcolors)):
    if Surfcolors[i] == col:
      idx=i
      break
  #endfor i in range(len(Surfcolors))
  return idx

def getcmapindex(cma):
  global Cmaps
  idx = -1
  for i in range(len(Cmaps)):
    if Cmaps[i] == cma:
      idx=i
      break
  #endfor i in range(len(Cmaps))
  return idx

def getmode3dindex(cmo):
  global Mode3ds
  idx = -1
  for i in range(len(Mode3ds)):
    if Mode3ds[i] == cmo:
      idx=i
      break
  #endfor i in range(len(mode3ds))
  return idx

def getwin():
  #tit = plt.gcf().canvas.get_window_title()
  tit = window_get_title()
  if tit == 'Figure 1':
    plt.close(plt.gcf())
    window()
  #endif
  #return plt.gcf().canvas.get_window_title()
  return window_get_title()

def set_win_visible(vis=True):
  plt.gcf().set_visible(vis)
  plt.show(block=False)
#enddef

def setwin(wintit):

  global Fig, Ax, Figman, Nwins, Nfigs


  ifig = -1
  fnums = plt.get_fignums()
  Nwins = len(fnums)
  Nfigs = Nwins

  if type(wintit) == str:

    tlis = plt.get_figlabels()

    for i in range(Nwins):
      if tlis[i] == wintit:
        ifig = fnums[i]
        break
      #endif
    #endfor

  elif type(wintit) == int:

    ifig = wintit

  #endif

  if ifig == -1:
    print("*** Window " + wintit + " not found  ***")
    return
  #endif

  Fig = plt.figure(ifig)

  if len(Fig.axes) == 0:
    #print("no axes")
    Ax = Fig.add_subplot(111,visible=False,label='111')
  else:
    #print("axes!")
    Ax = Fig.axes[0]
  #endif

  plt.sca(Ax)
  Figman = plt.get_current_fig_manager()

#  for win in Figs:
#    #if win.canvas.get_window_title() == wintit:
#    if window_get_title() == wintit:
#      for ax in win.axes:
#        plt.sca(ax)
#      break
#    #endif
#  #for win in Figs

#enddef setwin(wintit)

def getgeo():

    fig = plt.gcf()

    geo=fig.canvas.manager.window.wm_geometry()

    git = geo.split('+')
    wh = git[0].split('x')

    x = int(git[1]); y = int(git[2])
    wid = int(wh[0]); h = int(wh[1]);

    return wid,h,x,y
#def getgeo()

def nexist(nt):
  idx = GetIndexN(nt,1)
  if idx < 0: return 0
  else: return 1
#enddef nexist(nt)

def printg3(x): print('{:.3g}'.format(x))
def printg4(x): print('{:.4g}'.format(x))
def printg5(x): print('{:.5g}'.format(x))

def g8(x): return '{:.8g}'.format(x)
def g7(x): return '{:.7g}'.format(x)
def g6(x): return '{:.6g}'.format(x)
def g5(x): return '{:.5g}'.format(x)
def g4(x): return '{:.4g}'.format(x)
def g3(x): return '{:.3g}'.format(x)

def set_y_title(y='!'):
  global YTitle, Ytitle, Nxzone, Nyzone
  if y == '!':
    y = YTitle
  if Nyzone > 1:
    y += (Nyzone-1)*0.05
  Ytitle = y
def get_y_title():
  return Ytitle

def get_y_stat(): return Ystat

def set_number_of_ticks(n=-9):
  global NXtick
  NXtick = n
#enddef

def set_number_of_ticks_3d(n=-9):
  global NXtick3d
  NXtick3d = n
#enddef

def opttitles(k=1): global Ktitles; Ktitles=k

def optecho(k=1): global Kecho; Kecho = k
def getecho(): global Kecho; return Kecho

def getcolormap():
  global CMap
  return CMap
#enddef

def setcolormap(cmap='jet'):
  global CMap, Cmap
  if cmap == '!': CMap = cmap
  Cmap = CMap
#enddef

def setcolorbarpad(pad='!'):
  global ColorbarPad, Colorbarpad
  ColorbarPad = pad
  _setcolorbarpad()
#enddef

def getcolorbarpad(): global ColorbarPad; return ColorbarPad

def optfit(k=1):
  global Kfit, Kecho
  if Kecho: print("optfit(k=" + str(k) + ")")
  Kfit=k
#enddef

def optnfit(k=1):
  global Kfit, Kecho
  if Kecho: print("opntfit(k=" + str(k) + ")")
  if k: optfit(False)
  else: optfit(True)
#enddef

def getfit(): global Kfit; return Kfit

def optstat(k=1):
  global Kstat, Kecho
  if Kecho: print("optstat(k=" + str(k) + ")")
  Kstat=k
#enddef

def optnstat(k=1):
  global Kstat, Kecho
  if Kecho: print("opntstat(k=" + str(k) + ")")
  if k: optstat(False)
  else: optstat(True)
#enddef

def getstat(): global Kstat; return Kstat

def settopmargin(mar=0.86):
  global TopMargin
  TopMargin=mar
  plt.subplots_adjust(top=TopMargin)
  shpl()
#enddef

def gettopmargin(): global TopMargin; return TopMargin

def setbottommargin(mar=0.92):
  global BottomMargin
  BottomMargin=mar
  plt.subplots_adjust(bottom=BottomMargin)
  shpl()
#enddef

def getbottommargin(): global BottomMargin; return BottomMargin

def setrightmargin(mar=0.12):
  global RightMargin
  RightMargin=mar
  plt.subplots_adjust(right=RightMargin)
  shpl()
#enddef

def getrightmargin(): global RightMargin; return RightMargin

def setleftmargin(mar=0.12):
  global LeftMargin
  LeftMargin=mar
  plt.subplots_adjust(left=LeftMargin)
  shpl()
#enddef

def getleftmargin(): global LeftMargin; return LeftMargin

def setklegend(kl=True): global Klegend; Klegend=kl
def getklegend(): global Klegend; return Klegend
def setktitles(kt=True): global Ktitles; Ktitles=kt
def getktitles(): global Ktitles; return Ktitles

def set_plot_params():
  global Debug
  _setcolormap()
  _sethistcolor()
  _setaxislabelsize()
  _setaxistitlesize()
  _setaxistitledist()
  _setaxislabeldist()
  set_x_fit()
  set_x_stat()
  set_y_fit()
  set_y_stat()
  set_y_title()
  _setmode2d()
  _setmode3d()
  _setcolorbarpad()
#enddef set_plot_params()

def set_plot_params_3d():
  set_plot_params()
  _setaxistitlesize3d()
  _setaxistitledist3d()
  _setaxislabeldist3d()
  _setcolorbarpad()
#enddef set_plot_params()

def mred():
  setmarkercolor('red')
def mblue():
  setmarkercolor('blue')
def myellow():
  setmarkercolor('yellow')
def mcyan():
  setmarkercolor('cyan')
def mgreen():
  setmarkercolor('green')
def mmagenta():
  setmarkercolor('magenta')
def mgray():
  setmarkercolor('gray')

def lred():
  setlinecolor('red')
def lblue():
  setlinecolor('blue')
def lyellow():
  setlinecolor('yellow')
def lcyan():
  setlinecolor('cyan')
def lgreen():
  setlinecolor('green')
def lmagenta():
  setlinecolor('magenta')
def lgray():
  setlinecolor('gray')

def tred():
  settextcolor('red')
def tblue():
  settextcolor('blue')
def tyellow():
  settextcolor('yellow')
def tcyan():
  settextcolor('cyan')
def tgreen():
  settextcolor('green')
def tmagenta():
  settextcolor('magenta')
def tgray():
  settextcolor('gray')

def lmred():
  setlinecolor('red')
  setmarkercolor('red')
def lmblue():
  setlinecolor('blue')
  setmarkercolor('blue')
def lmyellow():
  setlinecolor('yellow')
  setmarkercolor('yellow')
def lmcyan():
  setlinecolor('cyan')
  setmarkercolor('cyan')
def lmgreen():
  setlinecolor('green')
  setmarkercolor('green')
def lmmagenta():
  setlinecolor('magenta')
  setmarkercolor('magenta')
def lmgray():
  setlinecolor('gray')
  setmarkercolor('gray')

def lhmred():
  setlinecolor('red')
  sethistcolor('red')
  setmarkercolor('red')
def lhmblue():
  setlinecolor('blue')
  sethistcolor('blue')
  setmarkercolor('blue')
def lhmyellow():
  setlinecolor('yellow')
  sethistcolor('yellow')
  setmarkercolor('yellow')
def lhmcyan():
  setlinecolor('cyan')
  sethistcolor('cyan')
  setmarkercolor('cyan')
def lhmgreen():
  setlinecolor('green')
  sethistcolor('green')
  setmarkercolor('green')
def lhmmagenta():
  setlinecolor('magenta')
  sethistcolor('magenta')
  setmarkercolor('magenta')
def lhmgray():
  setlinecolor('gray')
  sethistcolor('gray')
  setmarkercolor('gray')

def setfillcolor(fillcolor='none'):
  global FillColor
  FillColor = fillcolor
#enddef setfillcolor()

def getfillcolor():
  global FillColor
  return FillColor
#enddef getfillcolor()

#----------------------------------------------------------------------
set_plot_params_3d()
set_plot_params()
setlinecolor()
setmarkercolor()
sethistcolor()
#----------------------------------------------------------------------
#begin of aliases in m_hbook

# aliases

h1info = H1Info
h2info = H2Info

hfill2 = h2fill
hfill1 = h1fill

h1index = GetIndexH1
h2index = GetIndexH2
hindex1 = h1index
hindex2 = h2index

h1list = H1List
hlist1 = H1List
h2list = H2List
hlist2 = H2List

h1book = hbook1
h2book = hbook2

hplot1 = hplot1d
h1plot = hplot1d
hplot2 = hplot2d
h2plot = hplot2d
hpl = hplot

wc = window_close
wcl = window_clear
cc = _clearCanvas

pp = pplot
plotopt = plotoptions
ncolonlist = nlistcolon

vpl = vplxy

nupdh = nupdate_header
rstat = reset_status

medf = make_empty_dataframe
mdf = make_dataframe
vton = make_dataframe
vw = vwritex
vwxy = vwritexy
vwxyz = vwritexyz
vwxyt = vwritexyzt

setecho = optecho

vrxy = vreadxy

tbox = textbox

getzones = getzone
win = window

print3g = printg3
print4g = printg4
print5g = printg5

getconsole = get_console
cd = mhb_cd
mkdir = mhb_mkdir

nlc = nextlinecolor
nmc = nextmarkercolor
nhc = nexthistcolor

pp = pplot
Q = Quit
vplotxy = vplxy
setmarkerstyle = setmarkertype

nexists = nexist
npl = nplot
nnpl = nnplot
vexpand = vspline
vgrid = vspline
nexpand = nspline
ngrid = nspline
atan = np.arctan
polygon = polygon2d
nconcat = nappend

p2g = pg2
p3g = pg3
p4g = pg4
p5g = pg5
p6g = pg6
p7g = pg7
p8g = pg8
p9g = pg9
p10g = pg10

lmblau = mblue
lmrot = mred
lmhellblau = lmcyan
lmlila = lmmagenta

mblau = mblue
mrot = mred
mhellblau = mcyan
mlila = mmagenta

tblau = tblue
trot = tred
thellblau = tcyan
tlila = tmagenta

lblau = mblue
lrot = mred
lhellblau = lcyan
llila = lmagenta


lhmhellblau = lhmcyan
lhmlila = lhmmagenta

vfitxypoly = vfitpoly
circ = circle
vgrafit = vfitp1
ngrafit = nfitp1
shell = os.system
she = os.system

setdump = optdump
ndelet = ndelete
pmark = vplm
fexists = fexist

vfitg = vfitgauss

cdu = mhb_cd_up
cdd = mhb_cd_down
cup = mhb_cd_up
cdo = mhb_cd_down
wcd = mhb_cd
wpwd = mhb_pwd

nentries = nentry
nplp = nprof
setxstat = set_x_stat
setystat = set_y_stat
getxstat = get_x_stat
getystat = get_y_stat
#end of aliases in m_hbook

#end of m_hbook
def plotoptions_unklar(plopt=''):
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
  Textcolor, WaveFilePrefix, \
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
  ZoomXmin,ZoomXmax, ZoomYmin, ZoomYmax,\
  Tdate, TdateOv, Trun, TrunOv, Icallfromoverview,\
  LogX,LogY, LogZ, NxBinMax, Khdeleted, Waveplot, \
  Mrun, Mcomment, Mdate, ROFx, Rofy, Hull3D, Kgrid, KxAxis,KyAxis,KzAxis,Kbox, \
  FillColor,WisLinux


  Iplotopt = 0

  if type(plopt) != str: return

  if re.search('2d',plopt.lower()): plopt = re.sub('2[dD]',Mode2d,plopt)
  elif re.search('3d',plopt.lower()): plopt = re.sub('3[dD]',Mode3d,plopt)

  Ihist = 0
  Iprof = 0
  Ispline = 0
  Iboxes = 0
  Isame = 0
  Iscatter = 0
  Iscat3d = 0
  Icont3d = 0
  Itrisurf = 0
  Iclosed = 0
  Ierr = 0
  Isurf = 0
  Iline = 0
  Iinter = 0
  Ifill1d = 0
  Imarker = 0
  Inoempty = 0

  if not Fig: window()

  if re.search('marker',plopt) or re.search('M',plopt) or re.search('P',plopt):
    Iplotopt = 1;  Imarker = 1
  if re.search('spline',plopt) or re.search('C',plopt):
    Iplotopt = 1;  Ispline = 1; plopt = re.sub("spline","",plopt)
  if re.search('same',plopt) or re.search('S',plopt) or IsameGlobal.get() == 1:
    Iplotopt = 1;  Isame = 1
  if re.search('scatter',plopt):
    Iplotopt = 1;  Iscatter = 1

  if re.search('boxes',plopt): Iplotopt = 1;  Iboxes = 1

  if not Iboxes and re.search('scat3d',plopt): Iplotopt = 1;  Iscat3d = 1
  if not Iboxes and re.search('cont3d',plopt): Iplotopt = 1;  Icont3d = 1
  if not Iboxes and re.search('trisurf',plopt): Iplotopt = 1;  Itrisurf = 1
  if not Iboxes and re.search('inter',plopt): Iplotopt = 1;  Iinter = 1
  if not Iboxes and re.search('hist',plopt) or re.search('H',plopt): Iplotopt = 1;  Ihist = 1

  if re.search('err',plopt) or re.search('E',plopt): Iplotopt = 1;  Ierr = 1
  if Itrisurf == 0 and re.search('surf',plopt): Iplotopt = 1;  Isurf = 1
  if re.search('line',plopt) or re.search('L',plopt): Iplotopt = 1;  Iline = 1
  if re.search('closed',plopt): Iplotopt = 1;  Iclosed = 1
  if re.search('fill1d',plopt) or re.search('F',plopt): Iplotopt = 1;  Ifill1d = 1
  if re.search('prof',plopt) or re.search('P',plopt): Iplotopt = 1;  Iprof = 1
  if re.search('noempty',plopt) or re.search('N',plopt): Iplotopt = 1;  Inoempty = 1

  if not Linestyle or Linestyle == 'none' or Linewidth <= 0.: Iline = 0
  if not Markertype or Markertype == 'none' or Markersize <= 0.: Imarker = 0

#enddef plotoptions(plopt=''):

def waveray(nt="nray",fray='wave_ray.dat',iplot=1):

  nray = ncread(nt,"z:y:s0:p1:p2:p3",fray,skiphead=1)

  Fray = open(fray,'r')
  headray = Fray.readline().split()
  Fray.close()

  nx=int(headray[0])
  ny=int(headray[1])
  krun=int(headray[2])
  ener=float(headray[3])
  code=headray[4]
  xmin = nray.z.min()
  xmax = nray.z.max()
  dx = (xmax-xmin)/(nx-1)
  ymin = nray.y.min()
  ymax = nray.y.max()
  dy = (ymax-ymin)/(ny-1)

  htit = fray + ", run " + str(krun) + ", E=" + '{:.5g}'.format(ener) + " eV"
  htitnob = fray + "_" + str(krun) + "_" + '{:.5g}'.format(ener) + "_eV"
  hray = hbook2("hray",htit,
                nx,xmin-dx/2.,xmax+dx/2.,
                ny,ymin-dy/2.,ymax+dy/2.,
                1)
  nproj2(nray,"z:y","s0",idh="hray")

  if iplot:
    hplot2("hray")
    txyz(htit,"Phi [mrad]","Theta [mrad]","S0")
    pp(htitnob + ".pdf")
  #endif

#enddef

def write(*args):
  File = args[0]
  line = ""
  for a in args[1:]: line += " " + str(a)
  File.write(line+"\n")
#enddef

def hexist(idh):
  i = GetIndexH1(idh)
  if i == -1: i = GetIndexH2(idh)
  if i == -1: return 0
  else: return 1
#enddef

def optitight(itight=0):
  global Itight
  if itight: Itight = 1
  else: Itight = 0
#endif

def setxspace(xs=0.4):
  # XWIN
  global Xspace
  Xspace = xs
#enddef

def setyspace(ys=0.5):
  # YWIN
  global Yspace
  Yspace = ys
#enddef

def setleftmargin(m=0.12):
  #PAD
  global LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace
  LeftMargin = m
#enddef
def settopmargin(m=0.86):
  #PAD
  global LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace
  TopMargin = m
#enddef
def setbottommargin(m=0.13):
  #PAD
  global LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace
  BottomMargin = m
#enddef
def setrightmargin(m=0.92):
  #PAD
  global LeftMargin,RightMargin,TopMargin,BottomMargin, Xspace, Yspace
  RightMargin = m
#enddef

# chi2.ppf: Returns the chi2-value that should not be exceeded for a good fit
#           and a given convidence level
# chi2.cdf: Returns the integral of the chi2-distribution up to the given chi2
# chi2.sf:  Returns the probability that the hypothesis is true


#Example: Ndf = 77
#chi2.ppf(0.95,77) = 98.48438345934042 the chi2 > 98.48 rejects the hypothesis
#chi2.sf(98.48438345934042,77) =  0.05
#chi2.cdf(98.48438345934042,77) = 0.95

def chi2_crit(crit,ndf):
  #find Chi-Square critical value
  return chi2.ppf(crit, df=ndf)
#enddef
def chi2_prob(ch2,ndf):
  return chi2.sf(ch2,ndf)
#enddef
def chi2ndf_prob(ch2n,ndf):
  return chi2.sf(ch2n*ndf,ndf)
#enddef
def pdf(vals,cnam='x'):
  # see also make_dataframe etc....
  return pd.DataFrame(vals,columns=[cnam])
#enddef

def setsurfcolor(col):
  global Surfcolor
  Surfcolor = col
#enddef setsurfcolor(col)

def vcross(x,y):
  z = deepcopy(x)
  z[0] = x[1]*y[2] - x[2]*y[1]
  z[1] = x[2]*y[0] - x[0]*y[2]
  z[2] = x[0]*y[1] - x[1]*y[0]
  return z
#enddef vcross(x,y)

def setisame(isa=1):
  global Isame,IsameGlobal
  Isame = isa
  IsameGlobal.set(isa)
#enddef setisame(isa=1)

def vzero(n): return vcre(n,0.,0.)

def listsub(l1,l2):
  if len(l1) != len(l2):
    print("*** Error in listsub: Different length of lists")
    return None
  #endif
  l21 = []
  for i in range(len(l1)):
    l21.append(l2[i]-l1[i])
  #endfor
  return l21
#enddef

def listnorm(l):
  ln = []
  s2 = 0
  for i in range(len(l)):
    s2 += l[i]**2
  #endfor
  s = s2**0.5
  for i in range(len(l)):
    ln.append(l[i]/s)
  #endfor
  return ln
#enddef
def pnl(): print("\n")

def ivmin(v): return pd.Series(v).idxmin()
def vmin(v): return v.min(v)
def ivmax(v): return pd.Series(v).idxmax()
def vmax(v): return v.max(v)

def vminmax(x='?',y=''):

  if type(x) == str and x == '?':
    print("vminmax(x,[y]) returns idxmin, x[idxmin], y[idxmin], idxmax, x[idxmax], y[idxmax]\n or only values for x if y is missing. ")
  #endif

  if type(y) == str and y == '':
    im = ivmin(x)
    ip = ivmax(x)
    return im, x[im], ip, x[ip]
  else:
    im = ivmin(y)
    ip = ivmax(y)
    return im, x[im], y[im],ip, x[ip], y[ip]
  #endif

#enddef
# End of sequence m_hbook
########################################################

