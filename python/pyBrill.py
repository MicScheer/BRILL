
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
#+DECK,pyBrill,T=PYTHON.

import os,sys,platform,shutil,time,re

import tkinter as tk
from tkinter import *

import numpy as np
from scipy import special

import matplotlib as mpl
import matplotlib.pyplot as plt

import m_hbook as m
from m_hbook import *


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

def fqnke(n,K,Kyx):

    #NOTE: x and y reversed compared to brill_ellip.kumac

    # n is integer order of Fn(K)
    nm = int((n-1)/2)
    np = int((n+1)/2)

    K2 = K*K
    Kx2 = K2 / (1.0+Kyx**2)
    Ky2 = K2 - Kx2
    Kx = Kx2**0.5
    Ky = 0.0
    if Ky2.min() > 0.0: Ky = Ky2**0.5

    K221 = 1. + K2/2.

    x = n*(Kx2-Ky2) / (4.*K221)

    Jm = special.jv(nm,x)
    Jp = special.jv(np,x)

    Ax = Kx * (Jp-Jm)
    Ay = Ky * (Jp+Jm)

    fnke = (n/K221)**2 * (Ax*Ax+Ay*Ay)

    qnke = K221 * fnke / n

    return fnke,qnke

#enddef fqnk(n,K)

def fqnk(n,K):

    # n is integer order of Fn(K)
    nm = int((n-1)/2)
    np = int((n+1)/2)

    K2 = K*K
    K221 = 1. + K2/2.

    x = n*K2 / (4.*K221)

    Jm = special.jv(nm,x)
    Jp = special.jv(np,x)

    fnk = n * K / K221 * (Jm - Jp)
    fnk = fnk * fnk

    qnk = K221 * fnk / n

    return fnk,qnk

#enddef fqnk(n,K)

def calc_brill(l='?', nKvals=101, Kmin=0.5, Kmax=3., n=100, ebeam=1.722, curr=0.1,
          emitx=4.4, emity=0.066, betx=14., bety=3.4,
          sige=0.001, mode=1):

    global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList

    global Calculated,Nmax

    if type(l) == str and l == '?' :
        print("calc_brill(l='?', nKvals=101,Kmin=0.5, Kmax=3., n=100, ebeam=1.722, curr=0.1,")
        print("emitx=4.4, emity=0.066, betx=14., bety=3.4,sige=0.001, mode=1)")
        return

    dK = (Kmax-Kmin)/(nKvals-1)
    Kvals = np.arange(Kmin,Kmax+dK,dK)
    b0 = Kvals/(echarge1 * L /1000./(2.*pi1*emasskg1*clight1))

    if type(KyxList) != int:
        print("\n***********************************************************")
        print("Attention: For elliptical undulators K is shift-dependend!")
        print("Thus K must actually set for each harmonic...\n")
        print("***********************************************************\n")
        if Nmax > len(KyxList)*2-1:
            print("*** Uptonow, the max. harmonic is limited to len(KyxList)*2-1 ***\n")
            Nmax = len(KyxList)*2-1
        #endif
        time.sleep(3)
    #endif type(KyxList) != int

    Fkn = []
    Qn = []
    F = []
    FB = []
    FC = []
    FD = []
    B = []
    Lam = []
    Harm = []
    Sigr = []
    Sigrp = []

    if type(KyxList) != int:
        if Nmax > len(KyxList)*2-1:
            print("*** Uptonow, the max. harmonic is limited to len(KyxList)*2-1 ***\n")
            Nmax = len(KyxList)*2-1
        #endif
    #endif

    for k in range(Nmax+1):

        i = k
        if i == 0: i = 1

        if type(KyxList) == int:
            fk, q  = fqnk(i,Kvals)
        else:
            fk, q  = fqnke(i,Kvals,KyxList[int(i/2)])
        #endif type(KyxList) == int:

        Fkn.append(fk)
        Qn.append(q)

        F.append(0.0)
        FB.append(0.0)
        FC.append(0.0)
        FD.append(0.0)
        B.append(0.0)
        Lam.append(0.0)
        Harm.append(0.0)
        Sigr.append(0.0)
        Sigrp.append(0.0)
    #endfor k in range(10)

    #for i in [1,3,5,7,9,11]

    sigx=np.sqrt(emitx*1.e-9*betx)*1000. #mm
    sigxp=np.sqrt(emitx*1.e-9/betx)*1000. #mrad
    sigy=np.sqrt(emity*1.e-9*bety)*1000.
    sigyp=np.sqrt(emity*1.e-9/bety)*1000.

    if mode == 1:
        print("\n  Mode=1 (Kim):")
        print("\n  sigr := sqrt(lambda*length)/4/pi")
        print("  sigrp := sqrt(lambda/length)")
        print("  sigr*sigrp := lambda/(4*pi)\n\n")
    elif mode == 2:
        print("\n  Mode=2 (Walker, recommended):")
        print("\n  sigr := sqrt(2*lambda*length)/2/pi")
        print("  sigrp := sqrt(lambda/2/length)")
        print("  sigr*sigrp := lambda/(2*pi)\n\n")
    elif mode == -1:
        print("\n  Mode=-1 Erik Wallen, MAX-IV:")
        print("\n  sigw[i]=0.36/[i]/[N] | width of harmonic")
        print("  mulam[i]=sqrt(1.+(2.*[sigE]*[i]*[N]/0.36)**2)")
        print("  sigr[i]=sqrt(lam1/[i]*[N]*[l]/(8.*pi**2))") # like Kim
        print("  sigrp[i]=sqrt(mulam[i]*lam1/[i]/(2.*[N]*[l]))") # rad like Walker
        print("  sigr*sigrp := lambda/(4*pi)\n\n") # like Kim
    elif mode == 3:
        print("\n  Mode=3:")
        print("\n  sigr := sqrt(sqrt(2)*lambda*length)/2/pi")
        print(" 	sigrp := sqrt(lambda/2/sqrt(2)/length)")
        print(" 	sigr*sigrp := lambda/sqrt(2)/(2*pi)\n\n")
    else:
        print("\n *** Mode has no meaning, set to 1 (Kim) *** \n")
        mode = 1
    #endif

    time.sleep(1)

    if not Calculated_Brill:
        print('\n For the flux calculation the enery-spread is not taken into account, since')
        print(' it has no effect for the brillant flux and a smaller effect for the max. flux')
        print(' compared to the effect on the flux-density.')
        print(" The flux is always calculated according to Kim's formular (xray")
        print(' data booklet eqn. 17). This overestimates the flux by a factor')
        print(' of two, but since the max. is about a this factor higher than')
        print(' the on-resonant flux, this seems to be alright.\n\n')
        Calculated = True

    time.sleep(1)

    pi = pi1

    for i in range(1,Nmax+1,2):

        lami = 13.056*(1.+Kvals**2/2.)*L/10./ebeam**2/1.e10*1.e3/i # mm
        Lam[i] = lami

        harmi = i*0.950/(1.+Kvals**2/2.)/(l/10.)*ebeam**2 # keV
        Harm[i] = harmi

        rsigphi = \
        (1./(i*n*2.*np.sqrt(2.0))) / \
        np.sqrt((1.0/(i*n*2.*np.sqrt(2.0)))**2+(2.0*sige)**2)

        if mode == 1:
            sigrpi = np.sqrt(lami/(n*l)) # rad
            sigri = lami/4./pi/sigrpi # mm
        elif mode == 2:
            sigrpi = np.sqrt(lami/2./(n*l)) # rad
            sigri = lami/2./pi/sigrpi # mm
        elif mode == -1:
            #MAX IV formulas (Erik Wallen?)
            sigwi = 0.36/i/n # width of harmonic
            mulami = np.sqrt(1.+(2.*sigE*i*n/0.36)**2)
            rsigphi = 0.
            sigri = np.sqrt(lam1/i*n*l/(8.*pi**2)) # like Kim
            sigrpi = np.sqrt(mulami*lam1/i/(2.*n*l)) # rad # like Walker

        elif mode == 3:
            sigrpi = np.sqrt(lami/(n*l*2.*np.sqrt(2.))) #rad
            sigri=np.sqrt(np.sqrt(2.)*lami*(n*l))/2./pi #mm
        #endif mode == ...

        Sigr[i] = sigri
        Sigrp[i] = sigrpi

        sigrpi=sigrpi*1000 # mrad

        sigrxi = np.sqrt(sigx**2+sigri**2)
        sigrpxi = np.sqrt(sigxp**2+sigrpi**2)
        sigryi = np.sqrt(sigy**2+sigri**2)
        sigrpyi = np.sqrt(sigyp**2+sigrpi**2)

        rsigxyi = sigrpi**2/(sigrpxi*sigrpyi)

        F[i] = 1.431e14 * n * Qn[i] * curr # flux, data booklet Kim
        FD[i] = 1.744e14 * n**2 * ebeam**2 * curr * Fkn[i] # spectral brightness, data booklet Kim
        FB[i] = FD[i]*2. * pi * sigrpi**2 # brilliant flux
        FD[i] = FD[i]*rsigphi*rsigxyi # d.h. keine Emittanzwirkung auf FB!

        if mode == 1:
            B[i] = F[i] * rsigphi / (2.*pi)**2 / (sigrxi*sigrpxi*sigryi*sigrpyi)
        elif mode == 2:
            B[i] = FB[i] * rsigphi / (2.*pi)**2 /(sigrxi*sigrpxi*sigryi*sigrpyi)
        elif mode == -1:
            # MAX IV formulas (Erik Wallen?)
            rsigxyi = sigrpi**2 / (sigrpxi*sigrpyi)
            F[i] = 1.431e14 * n * Qn[i] * curr/2. # flux, data booklet Kim divived by 2 (MAX IV formula)
            FB[i] = 1.431e14 * n * Qn[i] * curr/2. # flux, data booklet Kim divived by 2 (MAX IV formula)
            B[i] = FB[i] / (2.*pi)**2 / (sigrxi*sigrpxi*sigryi*sigrpyi)
        elif mode == 3:
            B[i] = FB[i] * rsigphi/(2.*pi)**2 / (sigrxi*sigrpxi*sigryi*sigrpyi)
        #endif mode

        FC[i]= B[i] * (lami/2.)**2 * 1.0e6 # mrad**2 -> rad**2!!

    #endfor i in range(1,Nmax+1,2)

#enddef calc_brill(l='?', n=100, ebeam=1.722, curr=0.1,

# Wrting
def write_brill(fout="brilliance.dat"):

    global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList

    global Nmax

    Fout = open(fout,'w')
    Fout.write("* n, Harm, Flux-dens., Flux, Brilliance, Fbrill, Fcoh, SigR, SigRp\n")

    ll = len(Harm[1]) - 1

    for i in range(1,Nmax+1,2):
        for l in range(len(Harm[i])):
            k = ll - l
            sw = str(i)
            sw += " " + '{:.4g}'.format(Harm[i][k])
            sw += " " + '{:.4g}'.format(FD[i][k])
            sw += " " + '{:.4g}'.format(F[i][k])
            sw += " " + '{:.4g}'.format(B[i][k])
            sw += " " + '{:.4g}'.format(FB[i][k])
            sw += " " + '{:.4g}'.format(FC[i][k])
            sw += " " + '{:.4g}'.format(Sigr[i][k])
            sw += " " + '{:.4g}'.format(Sigrp[i][k])
            sw += "\n"
            Fout.write(sw)
        #endfor k in range(len(F)))
    #endfor i in range(1,Nmax+1,2)

    Fout.close()
#enddef write_brill()
global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,ScreenWidth, ScreenHeight,Vsetup_Plot

global Esel
global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


global Unamelist,Useed
global LastPlot; LastPlot = ['','']

global MShWelcome
MShWelcome = False

def get_mshwelcome():
  global MShWelcome
  return MShWelcome
#enddef

def mshwelcome(program='pyBrill',year='2023'):

  global Kdate, Kbox, KxAxis, KyAxis, MShWelcome

  KdateO, KboxO, KxAxisO, KyAxisO = Kdate, Kbox, KxAxis, KyAxis

  MShWelcome = True

  optdate(False)

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

def _spec_key_press(ev):

  global LastPlot
  global Esel,IEsel,S_Esel,S_IEsel

  if LastPlot[0] == 'FdPin':

    Nepho = int(Dsetup['Nepho'][1])
    EphMin = float(Dsetup['EphMin'][1])
    EphMax = float(Dsetup['EphMax'][1])

    if Nepho > 1:
      dE = (EphMax-EphMin)/(Nepho-1)
    else:
      dE = 0
    #endif

    if ev.key == 'up' or ev.key == 'pageup' or ev.key == 'right':
      IEsel += 1
      if IEsel > Nepho: IEsel = 1
    elif ev.key == 'down' or ev.key == 'pagedown' or ev.key == 'left':
      IEsel -= 1
      if IEsel < 1: IEsel = Nepho
    #endif

    Esel = EphMin + (IEsel-1)*dE
    S_IEsel.set(IEsel)
    S_Esel.set(Esel)

    _pFdPin(LastPlot[1])

  #endif LastPlot

#enddef _spec_key_press(ev)

def _set_uname():

  global Unamelist,Useed,Dsetup

  Unamelist = [ \
  'Mthreads','Ebeam','Curr','Step','Nelec','Noranone','Icohere','Ihbunch', \
  'Bunchlen','BunchCharge','Modebunch','PinX','PinY','PinZ','PinW','PinH', \
  'NpinY','NpinZ','Modepin','ModeSphere','Perlen','Shift','Nper','Nharm', \
  'Harm','Beffv','Beffh','Nepho','EphMin','EphMax','Espread','BetaH', \
  'BetaV','EmitH','EmitV','Disph','Dispph','Dispv','Disppv','Modeph','Pherror', \
  'IFieldProp','PinXprop','PinWprop','PinHprop','NpinYprop','NpinZprop']

  Useed = [376577121, 52147852, -1273034815, -1963249100, 1195262240, \
  -1718716574, -224354675, 432587481, 1692325775, 1934175653, \
  -107106772, 648589804, -1919014861, -1763988460, -1039845022, \
  1414926465, -1214705659, 560082688, 527470902, -1078636718, \
  272932485, -356992740, -2013991490, -588501795, -1010120436, \
  -1558306344, -1116776222, 794926823, -1157173406, 63711032, \
  -1870802148, -674825931, -690546468, 1671737514, -224394481, \
  2026233226, -1141752469, 2061158685, -1225625467, -147464566, \
  1692325775, 1934175653, -107106772, 648589804, -1919014861, \
  -1763988460, -1039845022, 1414926465, -1214705659, 560082688, \
  527470902, -1078636718, 272932485, -356992740, -2013991490, \
  -588501795, -1010120436, -1558306344, -1116776222, 794926823, \
  -1157173406, 63711032, -1870802148, -674825931]

  Dsetup['Mthreads'] = ['Number of threads (<0: Use all cores)',-1]
  Dsetup['Icohere'] = ['',0]
  Dsetup['BunchCharge'] = ['Not yet',0.0]
  Dsetup['Bunchlen'] = ['Not yet',0.0]
  Dsetup['BunchCharge'] = ['Not yet',0.0]
  Dsetup['Ihbunch'] = ['Each Ihbunch_th bunch is recorded',1]
  Dsetup['Modebunch'] = ['Not yet',0]
  Dsetup['Modeph'] = ['Mode for phase-error',0]
  Dsetup['Noranone'] = ['No random change for first e-',1]
  Dsetup['Noranone'] = ['No random change for first e-',0]

#enddef _set_uname()

global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
global MBrill, Omenu, NMBrill, NOmenu, Myfont, Toolbar, \
Calculated_Brill, Fig, Ax, Grid, Calculated_Spec

global MSetup,MBrill,MSpec
global Kellip

BeamPar = ['Ebeam','Curr','EmitH','EmitV','BetaH','BetaV','SigE', \
'Disph','Dispph','Dispv','Disppv']
UnduPar = ['Perlen','Nper','Beffv','Beffh','Nharm','Harmonic','Shift']
BrillPar = ['nKvals','Kmin','Kmax','Nmin','Nmax','Mode']
SpecPar = ['Nelec','Modepin','ModeSphere','Nepho','EphMin','EphMax','PinX', \
'PinY','PinZ','PinW','PinH','NpinZ','NpinY','Step','Pherror', \
'IFieldProp','PinXprop','PinWprop','PinHprop','NpinYprop','NpinZprop','Ifixseed']
PlotPar = ['Mode3d','Markersize','Linewidth','Linecolor']

global LastSetUp_Esel
LastSetUp_Esel = 0

def _wesel(ev):

  global Wesel,CanWesel

  global Esel,IEsel,S_Esel,S_IEsel

  Nepho = int(Dsetup['Nepho'][1])
  EphMin = float(Dsetup['EphMin'][1])
  EphMax = float(Dsetup['EphMax'][1])

  Esel = ev.xdata

  if Nepho > 1:
    dE = (EphMax-EphMin)/(Nepho-1)
  else:
    dE = 0
  #endif

  if Esel <= 0.0:
    IEsel = 1
    Esel = EphMin
  elif Esel > EphMax:
    IEsel = Nepho
    Esel = EphMax
  #endif

  if Nepho > 1:
    IEsel = int((Esel-EphMin)/dE)+1
    if IEsel <=0:
      IEsel = 1
    elif IEsel > Nepho:
      IEsel = Nepho
    #endif
  else:
    IEsel = 1
  #endif

  Esel = EphMin + (IEsel-1)*dE

  S_IEsel.set(IEsel)
  S_Esel.set(Esel)

  Wesel.canvas.mpl_disconnect(CanWesel)
  window_close()

#enddef wesel()

def _sel_Esel():
  global Dsetup,Vsetup,Vsetup_Beam,Vsetup_Spec,Vsetup_Brill,Vsetup_Undu, \
  Vsetup_Plot, IEsel,Esel

  global Curr,EmitX,EmitV,BetaH,BetaV,SigE,Disph,Dispph,Dispv,Disppv,L,N,Beffv,Beffh,Nharm, \
  Harm,Shift,nKvals,Kmin,Kmax,Nmin,Nmax,Mode,Nelec,modepin,modesphere,Nepho, \
  EphMin,EphMax,PinX,PinY,PinZ,PinW,PinH,NpinZ,NpinY,Step,Pherror, \
  Ifixseed,Kellip

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global CanWesel, Wesel

  SetUp_Esel.destroy()

  window('Select Photon Energy','')
  if Calculated_Spec == False or nexist("nflx") == 0: _calc_spec()

  global Dsetup,Nepho,EphMax,EphMin,WmainMaster

  Nepho = int(Dsetup['Nepho'][1])
  EphMin = float(Dsetup['EphMin'][1])
  EphMax = float(Dsetup['EphMax'][1])

  Wesel = plt.gcf()
  CanWesel = Wesel.canvas.mpl_connect('button_press_event',_wesel)
  #Fig = plt.gcf()
  #Wmaster = Fig.canvas.toolbar.master

  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()

  if xm <= ScreenWidth/2:
    window_geometry(str(int(wm*0.75)) + 'x' + str(int(hm*0.75)) + '+' + str(int(xm+wm*1.1)) + '+' + str(int(ym-hm*0.1)))
  else:
    window_geometry(str(int(wm*0.75)) + 'x' + str(int(hm*0.75)) + '+' + str(int(xm-wm*1.1)) + '+' + str(int(ym-hm*0.1)))
  #endif

  optnstat()

  if Modepin != 0:
    zone(1,1)
  else:
    zone(1,2)
  #endif

  mso = float(Vsetup_Plot[1][1][1])
  lwo = float(Vsetup_Plot[2][1][1])
  lco = Vsetup_Plot[3][1][1]
  mso = getmarkersize()

  setmarkersize(5.)
  setlinewidth(2.)
  setlinecolor('b')

  if Modepin == 0:
    selzy = "abs(z-" + str(PinZ) + ") < 1.0e-10 and abs(y-" + str(PinY) + ") < 1.e-10"
    npl(nfld,"egam:s0",selzy,plopt='line')
    xtit="photon energy [eV]"
    ytit = 'N$_{\gamma}$' + '/mm$^2$/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
    titp = "\nS0 (x={:.3g}m, y={:.3g}mm, z={:.3g}mm)". \
    format(PinX/1000.,PinY,PinZ)
    txyz(titp,'',ytit)
    if Nepho < 100:
      npl(nfld,"egam:s0",selzy,plopt='samemarker')
    nextzone()
  #endif

  npl(nflx,"egam:s0",plopt='line')

  if Nepho < 100:
    npl(nflx,"egam:s0",plopt='markersame')

  xtit="photon energy [eV]"

  if NpinZ > 1 and NpinY > 1:
    ytit = 'N$_{\gamma}$' + '/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
  else:
    ytit = 'N$_{\gamma}$' + '/mm/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
  #endif NpinZ, NpinY

  titp = "\nS0 (w={:.3g}mm, h={:.3g}mm, x={:.3g}m, y={:.3g}mm, z={:.3g}mm)". \
  format(PinW,PinH,PinX/1000.,PinY,PinZ)

  txyz(titp,xtit,ytit)

  Vsetup_Plot[1][1][1] = mso
  Vsetup_Plot[2][1][1] = lwo
  Vsetup_Plot[3][1][1] = lco
  setmarkersize(mso)

#enddef _sel_Esel()

def _SetUpIn_Esel(event,kvar):
  global LastSetUp_Esel
  LastSetUp_Esel = [event,kvar]
#enddef _SetUpInSpec(event,kvar)

def _SetUpOut_Esel(event,kvar):

  global SetUp_Esel,LastSetUp_Esel, IEsel,Esel, S_Esel,S_IEsel
  global Dsetup,Nepho,EphMax,EphMin

  ev = LastSetUp_Esel[0].widget
  val = ev.get()
  if len(val.split('.')) > 1:
    val = float(val)
  else:
    val = int(val)
  #endif

  Nepho = int(Dsetup['Nepho'][1])
  EphMin = float(Dsetup['EphMin'][1])
  EphMax = float(Dsetup['EphMax'][1])

  if Nepho > 1:
    dE = (EphMax-EphMin)/(Nepho-1)
  else:
    dE = 0
  #endif

  if kvar == 1:
    IEsel = int(S_IEsel.get())
    if IEsel <=0:
      IEsel = 1
    elif IEsel > Nepho:
      IEsel = Nepho
    #endif
  elif kvar == 2:
    Esel = float(S_Esel.get())
    if Esel <= 0.0:
      IEsel = 1
      Esel = EphMin
    elif Esel > EphMax:
      IEsel = Nepho
      Esel = EphMax
    #endif
    if Nepho > 1:
      IEsel = int((Esel-EphMin)/dE)+1
      if IEsel <=0:
        IEsel = 1
      elif IEsel > Nepho:
        IEsel = Nepho
      #endif
    else:
      IEsel = 1
    #endif
  #endif
  Esel = EphMin + (IEsel-1)*dE
  S_Esel.set(Esel)
  S_IEsel.set(IEsel)
#enddef _SetUpOutSpec(event,kvar)

def _closeSetUp_Esel():
  global SetUp_Esel,LastSetUp_Esel,IEsel,Esel,S_Esel,S_IEsel
  LastSetUp_Esel = 0
  IEsel = int(S_IEsel.get())
  Esel = float(S_Esel.get())
  SetUp_Esel.destroy()
#def _closeSetUp_Esel

def _setup_esel():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  global SetUp_Esel,IEsel,Esel

  SetUp_Esel = Toplevel()
  SetUp_Esel.title('Select Photon Energy')
  SetUp_Esel.attributes('-topmost',1)
  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()
  SetUp_Esel.geometry('+' + str(int(xm+wm/2)) + '+' + str(int(ym+hm*0.7)))

  f = Frame(SetUp_Esel)
  flab = Label(f,text="Index of E_photon")
  fent =  Entry(f,text=S_IEsel)
  flab.pack(side=LEFT)
  fent.pack(side=RIGHT)
  fent.bind('<FocusIn>',lambda event,kvar=1:_SetUpIn_Esel(event,kvar))
  fent.bind('<FocusOut>',lambda event,kvar=1:_SetUpOut_Esel(event,kvar))
  fent.bind('<Return>',lambda event,kvar=1:_SetUpOut_Esel(event,kvar))
  f.pack(fill='x')

  f = Frame(SetUp_Esel)
  flab = Label(f,text="E_photon")
  fent =  Entry(f,text=S_Esel)
  flab.pack(side=LEFT)
  fent.pack(side=RIGHT)
  fent.bind('<FocusIn>',lambda event,kvar=2:_SetUpIn_Esel(event,kvar))
  fent.bind('<FocusOut>',lambda event,kvar=2:_SetUpOut_Esel(event,kvar))
  fent.bind('<Return>',lambda event,kvar=2:_SetUpOut_Esel(event,kvar))
  f.pack(fill='x')

  bSel = Button(SetUp_Esel,text='Select from spectrum',command=_sel_Esel)
  bSel.pack()
  bClose = Button(SetUp_Esel,text='Close',command=_closeSetUp_Esel)
  bClose.pack()

#enddef _setup_esel()

def _ini_Esel():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global Esel,IEsel,S_Esel,S_IEsel

  Nepho = int(Dsetup['Nepho'][1])
  EphMin = float(Dsetup['EphMin'][1])
  EphMax = float(Dsetup['EphMax'][1])

  if nexist("nbun"):
    s0max = nflx.s0.max()
    EphMaxS0 = nflx.query("s0=="+str(s0max)).egam.max()
    Esel = EphMaxS0
    IEsel = int(nflx.query("abs(egam-"+str(EphMaxS0)+")<1.e-10").iegam.max())
  else:
    IEsel = int((Nepho+1)/2)
    if IEsel > 1:
      Esel = EphMin + (EphMax-EphMax)/(IEsel-1)
    else:
      Esel = EphMin
    #endif
  #endif

  S_IEsel.set(IEsel)
  S_Esel.set(Esel)

#enddef _ini_Esel()

def _pFdProp(key='s0'):

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global Esel,IEsel,S_Esel,S_IEsel
  global LastPlot; LastPlot = ['FdProp',key]
  global nfdp

  if Calculated_Spec == False or nexist("nbun") == 0 \
  or nexist("nfdp") == 0: _calc_spec()

  s0max = nfdp.s0.max()
  if np.isnan(s0max) == True: return

  #getzone()
  optnstat()
  _set_plot_spec()

  keyu = key.upper()
  keyl = key.lower()

  if Esel <= 0 : _ini_Esel()
  elif Esel < EphMin :
    Esel = EphMin
    IEsel = 1
  elif Esel > nfld.egam.max() :
    Esel = EphMax
    IEsel = Nepho
  #endif

  S_IEsel.set(IEsel)
  S_Esel.set(Esel)

  #selgam = "abs(egam-" + str(Esel) + ")<1.0e-10"
  selgam = "iegam==" + str(IEsel)

  ymin = -PinHprop/2.
  ymax = PinHprop/2.
  zmin = -PinWprop/2.
  zmax = PinWprop/2.

  set_plot_params_3d()

  if keyu == 'S0' or keyu == 'S1' or keyu == 'S2' or keyu == 'S3' or keyu == 'P':

    if keyu == 'S0': htit = 'Distribution of S$_0$'
    elif keyu == 'S1': htit = 'Distribution of S$_1$'
    elif keyu == 'S2': htit = 'Distribution of S$_2$'
    elif keyu == 'S3': htit = 'Distribution of S$_3$'
    elif keyu == 'P': htit = 'Distribution of Power'

    plopt = Vsetup_Plot[0][1][1]

    if plopt == 'surf' or plopt == 'boxes' or plopt == 'inter':
      hnam = 'HpinProp_' + keyu
      hbook2(hnam,htit,NpinZprop,zmin,zmax,NpinYprop,ymin,ymax,overwrite=1)
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nproj2(nbun,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      else:
        nproj2(nfdp,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      #endif
      hplave(hnam,plopt)
    elif plopt == 'scat3d':
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nplot(nbun,"z:y:"+keyl+":"+keyl,selgam)
      else:
        nplot(nfdp,"z:y:"+keyl+":"+keyl,selgam)
      #endif
    else:
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nplot(nbun,"z:y",selgam,keyl)
      else:
        nplot(nfdp,"z:y",selgam,keyl)
      #endif
    #endif

    tunit = 'N$_{\gamma}$' + '/mm$^2$/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"

    cbp = getcolorbarpad()
    xuni = 0.97 + cbp
    yuni = 0.5
    auni = 90.

    ax = plt.gca()
    if type(ax) == Tax2d:
      zunit = ''
    else:
      zunit = tunit
    #endif
    if keyu == 'S0':
      txyz("Dens. of S$_0$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'S1':
      txyz("Dens. of S$_1$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'S2':
      txyz("Dens. of S$_2$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'S3':
      txyz("Dens. of S$_3$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'P':
      xuni = 1.0 + cbp
      yuni = 1.05
      auni = 0.0
      zunit = '[W/mm$^2$]'
      txyz("Dens. of Power for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    #endif

    if type(ax) == Tax2d:
      text(xuni,yuni,tunit,halign='left',angle=auni)
    #endif

  elif keyu == 'EYI' or keyu == 'EYR' or keyu == 'EZR' or keyu == 'EZI':

    if keyu == 'EYR': htit = 'Distribution of Ey_real'
    elif keyu == 'EYI': htit = 'Distribution of Ey_imag'
    elif keyu == 'EZR': htit = 'Distribution of Ez_real'
    elif keyu == 'EZI': htit = 'Distribution of Ez_imag'

    plopt = Vsetup_Plot[0][1][1]

    if plopt == 'surf' or plopt == 'boxes' or plopt == 'inter':
      hnam = 'HpinProp_' + keyu
      hbook2(hnam,htit,NpinZprop,zmin,zmax,NpinYprop,ymin,ymax,overwrite=1)
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nproj2(nbun,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      else:
        nproj2(nfdp,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      #endif
      hplave(hnam,plopt)
    elif plopt == 'scat3d':
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nplot(nbun,"z:y:"+keyl+":"+keyl,selgam)
      else:
        nplot(nfdp,"z:y:"+keyl+":"+keyl,selgam)
      #endif
    else:
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nplot(nbun,"z:y",selgam,keyl)
      else:
        nplot(nfdp,"z:y",selgam,keyl)
      #endif
    #endif

    tunit = 'Sqrt(N$_{\gamma}$' + '/mm$^2$/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA)"

    cbp = getcolorbarpad()
    xuni = 0.97 + cbp
    yuni = 0.5
    auni = 90.

    ax = plt.gca()
    if type(ax) == Tax2d:
      zunit = ''
    else:
      zunit = tunit
    #endif
    if keyu == 'EYR':
      txyz("Dens. of Ey_real for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'EYI':
      txyz("Dens. of Ey_imag for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'EZR':
      txyz("Dens. of Ez_real for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'EZI':
      txyz("Dens. of Ez_imag for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    #endif

    if type(ax) == Tax2d:
      text(xuni,yuni,tunit,halign='left',angle=auni)
    #endif

  elif keyu == 'P0' or keyu == 'P1' or keyu == 'P2' or keyu == 'P3':
    Quit("Baustelle Pprop")

    if keyu == 'P0': htit = 'Distribution of P$_0$'
    elif keyu == 'P2': htit = 'Distribution of P$_2$'
    elif keyu == 'P3': htit = 'Distribution of P$_3$'

    plopt = Vsetup_Plot[0][1][1]

    if plopt == 'surf' or plopt == 'boxes' or plopt == 'inter':
      for p in ['0','1','2','3']:
        hnam = 'Hpin_S' + p
        kl = 's' + p
        hbook2(hnam,htit,NpinZ,zmin,zmax,NpinY,ymin,ymax,overwrite=1)
        if Modepin != 0:
          print("\n***Modepin != 0 not yet for propagated fields available ***")
          return
          nproj2(nbun,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
        else:
          nproj2(nfld,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
        #endif
      #endfor
      hnam = 'Hpin_' + keyu
      if keyu == 'P1':
        htit = 'Distribution of P$_1$'
        hdiv('Hpin_S1','Hpin_S0',hnam,htit)
      #endif
      hplave(hnam,plopt)
    elif plopt == 'scat3d':
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nplot(nbun,"z:y:"+keyl+":"+keyl,selgam)
      else:
        nplot(nfld,"z:y:"+keyl+":"+keyl,selgam)
      #endif
    else:
      if Modepin != 0:
        print("\n***Modepin != 0 not yet for propagated fields available ***")
        return
        nplot(nbun,"z:y",selgam,keyl)
      else:
        nplot(nfld,"z:y",selgam,keyl)
      #endif
    #endif

    if keyu == 'P0': \
    txyz("P$_0$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")
    elif keyu == 'P1': \
    txyz("P$_1$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")
    elif keyu == 'P2': \
    txyz("P$_2$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")
    elif keyu == 'P3': \
    txyz("P$_3$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")

  #endif key

#enddef _pFdProp()

def _pFdPin(key='s0'):

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global Esel,IEsel,S_Esel,S_IEsel
  global LastPlot; LastPlot = ['FdPin',key]

  if Calculated_Spec == False or nexist("nbun") == 0 \
  or nexist("nfld") == 0: _calc_spec()

  s0max = nflx.s0.max()
  if np.isnan(s0max) == True: return

  #getzone()
  optnstat()
  _set_plot_spec()

  keyu = key.upper()
  keyl = key.lower()

  if Esel <= 0 : _ini_Esel()
  elif Esel < EphMin :
    Esel = EphMin
    IEsel = 1
  elif Esel > nfld.egam.max() :
    Esel = EphMax
    IEsel = Nepho
  #endif

  S_IEsel.set(IEsel)
  S_Esel.set(Esel)

  #selgam = "abs(egam-" + str(Esel) + ")<1.0e-10"
  selgam = "iegam==" + str(IEsel)

  ymin = PinY - PinH/2.
  ymax = PinY + PinH/2.
  zmin = PinZ - PinW/2.
  zmax = PinZ + PinW/2.

  set_plot_params_3d()

  if keyu == 'S0' or keyu == 'S1' or keyu == 'S2' or keyu == 'S3' or keyu == 'P':

    if keyu == 'S0': htit = 'Distribution of S$_0$'
    elif keyu == 'S1': htit = 'Distribution of S$_1$'
    elif keyu == 'S2': htit = 'Distribution of S$_2$'
    elif keyu == 'S3': htit = 'Distribution of S$_3$'
    elif keyu == 'P': htit = 'Distribution of Power'

    plopt = Vsetup_Plot[0][1][1]

    if plopt == 'surf' or plopt == 'boxes' or plopt == 'inter':
      hnam = 'Hpin_' + keyu
      hbook2(hnam,htit,NpinZ,zmin,zmax,NpinY,ymin,ymax,overwrite=1)
      if Modepin != 0:
        nproj2(nbun,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      else:
        nproj2(nfld,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      #endif
      hplave(hnam,plopt)
    elif plopt == 'scat3d':
      if Modepin != 0:
        nplot(nbun,"z:y:"+keyl+":"+keyl,selgam)
      else:
        nplot(nfld,"z:y:"+keyl+":"+keyl,selgam)
      #endif
    else:
      if Modepin != 0:
        nplot(nbun,"z:y",selgam,keyl)
      else:
        nplot(nfld,"z:y",selgam,keyl)
      #endif
    #endif

    tunit = 'N$_{\gamma}$' + '/mm$^2$/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"

    cbp = getcolorbarpad()
    xuni = 0.97 + cbp
    yuni = 0.5
    auni = 90.

    ax = plt.gca()
    if type(ax) == Tax2d:
      zunit = ''
    else:
      zunit = tunit
    #endif
    if keyu == 'S0':
      txyz("Dens. of S$_0$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'S1':
      txyz("Dens. of S$_1$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'S2':
      txyz("Dens. of S$_2$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'S3':
      txyz("Dens. of S$_3$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'P':
      xuni = 1.0 + cbp
      yuni = 1.05
      auni = 0.0
      zunit = '[W/mm$^2$]'
      txyz("Dens. of Power for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    #endif

    if type(ax) == Tax2d:
      text(xuni,yuni,tunit,halign='left',angle=auni)
    #endif

  elif keyu == 'EYI' or keyu == 'EYR' or keyu == 'EZI' or keyu == 'EZR':

    if keyu == 'EYR': htit = 'Distribution of Ey_real'
    elif keyu == 'EYI': htit = 'Distribution of Ey_imag'
    elif keyu == 'EZR': htit = 'Distribution of Ez_real'
    elif keyu == 'EZI': htit = 'Distribution of Ez_imag$'

    plopt = Vsetup_Plot[0][1][1]

    if plopt == 'surf' or plopt == 'boxes' or plopt == 'inter':
      hnam = 'Hpin_' + keyu
      hbook2(hnam,htit,NpinZ,zmin,zmax,NpinY,ymin,ymax,overwrite=1)
      if Modepin != 0:
        nproj2(nbun,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      else:
        nproj2(nfld,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
      #endif
      hplave(hnam,plopt)
    elif plopt == 'scat3d':
      if Modepin != 0:
        nplot(nbun,"z:y:"+keyl+":"+keyl,selgam)
      else:
        nplot(nfld,"z:y:"+keyl+":"+keyl,selgam)
      #endif
    else:
      if Modepin != 0:
        nplot(nbun,"z:y",selgam,keyl)
      else:
        nplot(nfld,"z:y",selgam,keyl)
      #endif
    #endif

    tunit = 'Sqrt(N$_{\gamma}$' + '/mm$^2$/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA)"

    cbp = getcolorbarpad()
    xuni = 0.97 + cbp
    yuni = 0.5
    auni = 90.

    ax = plt.gca()
    if type(ax) == Tax2d:
      zunit = ''
    else:
      zunit = tunit
    #endif
    if keyu == 'EYR':
      txyz("Dens. of Ey_real for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'EYI':
      txyz("Dens. of Ey_imag for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'EZR':
      txyz("Dens. of Ez_real for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    elif keyu == 'EZI':
      txyz("Dens. of Ez_imag for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]",zunit)
    #endif

    if type(ax) == Tax2d:
      text(xuni,yuni,tunit,halign='left',angle=auni)
    #endif

  elif keyu == 'P0' or keyu == 'P1' or keyu == 'P2' or keyu == 'P3':
    Quit("Baustelle Pspec")

    if keyu == 'P0': htit = 'Distribution of P$_0$'
    elif keyu == 'P2': htit = 'Distribution of P$_2$'
    elif keyu == 'P3': htit = 'Distribution of P$_3$'

    plopt = Vsetup_Plot[0][1][1]

    if plopt == 'surf' or plopt == 'boxes' or plopt == 'inter':
      for p in ['0','1','2','3']:
        hnam = 'Hpin_S' + p
        kl = 's' + p
        hbook2(hnam,htit,NpinZ,zmin,zmax,NpinY,ymin,ymax,overwrite=1)
        if Modepin != 0:
          nproj2(nbun,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
        else:
          nproj2(nfld,"z:y",keyl,selgam,idh=hnam,ioverwrite=0)
        #endif
      #endfor
      hnam = 'Hpin_' + keyu
      if keyu == 'P1':
        htit = 'Distribution of P$_1$'
        hdiv('Hpin_S1','Hpin_S0',hnam,htit)
      #endif
      hplave(hnam,plopt)
    elif plopt == 'scat3d':
      if Modepin != 0:
        nplot(nbun,"z:y:"+keyl+":"+keyl,selgam)
      else:
        nplot(nfld,"z:y:"+keyl+":"+keyl,selgam)
      #endif
    else:
      if Modepin != 0:
        nplot(nbun,"z:y",selgam,keyl)
      else:
        nplot(nfld,"z:y",selgam,keyl)
      #endif
    #endif

    if keyu == 'P0': \
    txyz("P$_0$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")
    elif keyu == 'P1': \
    txyz("P$_1$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")
    elif keyu == 'P2': \
    txyz("P$_2$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")
    elif keyu == 'P3': \
    txyz("P$_3$ for E$_{\gamma}$  = " + pg5(Esel) + " eV","z [mm]","y [mm]")

  #endif key

#enddef _pFdPin()

def _pElec(key='zizpi'):

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global MSpec
  global LastPlot; LastPlot = ['Elec',key]

  MSpec.unpost()

  if Calculated_Spec == False or nexist("nbun") == 0 \
  or nexist("nfld") == 0: _calc_spec()

  s0max = nflx.s0.max()
  if np.isnan(s0max) == True: return

  #getzone()
  optnstat()
  _set_plot_spec()

  keyu = key.upper()
  keyl = key.lower()

  set_plot_params_3d()

  #zone(1,1)

  if keyl == 'eel':
    htit = 'Energy Distribution'
    npl(nbun,"eel","iegam==1")
    txyz(htit,"E[GeV]")
  elif keyl == 'zizpi':
    nplot(nbun,"rzi1:zpi1","iegam==1")
    txyz('Horizontal Phase-space at Entrance',"z[mm]","z'[mrad]")
  elif keyl == 'ziyi':
    nplot(nbun,"rzi1:ryi1","iegam==1")
    txyz('Phase-space at Entrance',"z[mm]","y[mm]")
  elif keyl == 'zpiypi':
    nplot(nbun,"zpi1:ypi1","iegam==1")
    txyz('Phase-space at Entrance',"z'[mrad]","y'[mrad]")
  elif keyl == 'yiypi':
    nplot(nbun,"ryi1:ypi1","iegam==1")
    txyz('Vertical Phase-space at Entrance',"y[mm]","y'[mrad]")
  #endif

#enddef _pElec()

def _pFdSpec(key='s0'):

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global LastPlot; LastPlot = ['FdSpec',key]

  if Modepin != 0: return

  if Calculated_Spec == False or nexist("nflx") == 0: _calc_spec()

  keyu = key.upper()
  keyl = key.lower()

  s0max = nfld.s0.max()
  if np.isnan(s0max): return

  optnstat()

  setmarkersize(float(Vsetup_Plot[1][1][1]))
  setlinewidth(float(Vsetup_Plot[2][1][1]))
  setlinecolor(Vsetup_Plot[3][1][1])

  selzy = "abs(z-" + str(PinZ) + ") < 1.0e-10 and abs(y-" + str(PinY) + ") < 1.e-10"

  if keyl[0] == 's':
    npl(nfld,"egam:"+keyl,selzy,plopt='line')
    xtit="photon energy [eV]"
    ytit = 'N$_{\gamma}$' + '/mm$^2$/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
    #endif NpinZ, NpinY
    titp = "\n" + keyu + " (x={:.3g}m, y={:.3g}mm, z={:.3g}mm)". \
    format(PinX/1000.,PinY,PinZ)
  elif keyl[0] == 'p':
    Quit("Baustelle PFdSpec")
  #endif key

  txyz(titp,xtit,ytit)

#enddef _pFdSpec()

def _pFluxSpec(key='s0'):

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global LastPlot; LastPlot = ['FluxSpec',key]

  if Calculated_Spec == False or nexist("nflx") == 0: _calc_spec()

  kplot = 0
  keyu = key.upper()
  keyl = key.lower()

  s0max = nflx.s0.max()
  if not np.isnan(s0max): kplot = 1

  optnstat()
  #getzone()
  #zone(1,1)

  setmarkersize(float(Vsetup_Plot[1][1][1]))
  setlinewidth(float(Vsetup_Plot[2][1][1]))
  setlinecolor(Vsetup_Plot[3][1][1])

  if keyl[0] == 's':
    npl(nflx,"egam:"+keyl,plopt='line')
    xtit="photon energy [eV]"
    if NpinZ > 1 and NpinY > 1:
      ytit = 'N$_{\gamma}$' + '/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
    else:
      ytit = 'N$_{\gamma}$' + '/mm/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
    #endif NpinZ, NpinY
    titp = "\n" + keyu + " (w={:.3g}mm, h={:.3g}mm, x={:.3g}m, y={:.3g}mm, z={:.3g}mm)". \
    format(PinW,PinH,PinX/1000.,PinY,PinZ)
  elif keyl[0] == 'p':
    Quit("Baustelle P")
    npl(nflx,"egam:"+keyl,plopt='line')
    xtit="photon energy [eV]"
    if NpinZ > 1 and NpinY > 1:
      ytit = 'N$_{\gamma}$' + '/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
    else:
      ytit = 'N$_{\gamma}$' + '/mm/s/0.1' + '%BW/' + str(int(Curr*1000)) + "mA"
    #endif NpinZ, NpinY
    titp = "\n" + keyu + " (w={:.3g}mm, h={:.3g}mm, x={:.3g}m, y={:.3g}mm, z={:.3g}mm)". \
    format(PinW,PinH,PinX/1000.,PinY,PinZ)
  #endif key

  txyz(titp,xtit,ytit)

#enddef _pFluxSpec()

def _write_urad_phase_nam():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global Unamelist,Useed

  try:
    shutil.copyfile("urad_phase.nam","urad_phase.nam.bck")
  except: pass

  Fnam = open("urad_phase.nam",'w')

  Fnam.write(" $uradphasen\n\n")
  for var in Unamelist:
    if var == 'Harm':
      Fnam.write("  " 'harm='+ str(Dsetup['Harmonic'][1]) + '          !' + Dsetup['Harmonic'][0] + '\n')
    else:
      Fnam.write("  " + var+'='+ str(Dsetup[var][1]) + '          !' + Dsetup[var][0] + '\n')
  #endfor
  Fnam.write("\n $end\n\n")

  Fnam.write(" $seedn\n\n")
  Fnam.write('  Ifixseed=' + str(Dsetup['Ifixseed'][1]) + \
  '         !' + Dsetup['Ifixseed'][0] + '\n')
  i=0
  for var in Useed:
    i += 1
    Fnam.write('  irnseed(' + str(i) + ')=' + str(Useed[i-1]) + '\n')
  #endfor
  Fnam.write("\n $end\n\n")

  Fnam.close()

#enddef _write_urad_phase_nam()

def _write_urad_phase_nam_alt():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  shutil.copyfile("urad_phase.nam","urad_phase.nam.bck")

  Fnam = open("urad_phase.nam",'r')
  fin = Fnam.readlines()
  Fnam.close()

  Fnam = open("urad_phase.nam",'w')
  n = 0
  for line in fin:
    if line[0] == '$' or len(line.strip()) == 0:
      Fnam.write(line)
      continue
    #endif
    scom = line.split('!')
    if len(scom) > 1: com = " !" + scom[1].strip()
    else: com = ''
    sline = scom[0].split('=')
    if len(sline) > 1:
      svar = sline[0].strip().lower()
    else:
      Fnam.write(line)
      continue
    #endif
    ifound = 0

    i = 0
    for key in list(Dsetup):
      i += 1
      if key.lower() == svar:
        sval = str(Dsetup[key][1])
        Fnam.write("  " + key + "=" + sval + com + '\n')
        ifound = 1
      #endif
    #endfor
    if ifound == 0: Fnam.write(line)
  #endfor
  Fnam.close()
#enddef _write_urad_phase_nam()

def _mouse_click(ev):
  global Xmouse,Ymouse

  print("ev:",ev)
  #print("ev.xy:",ev.xy)
  #print(ev.xy[0])
  return
  Xmouse = ev.xdata
  Ymouse = ev.ydata
  if ev.button is MouseButton.LEFT:
    pass
  #endif
#enddef _mouse_click(event)

def _dvsetup():
  global Dsetup,Vsetup,Vsetup_Beam,Vsetup_Spec,Vsetup_Brill,Vsetup_Undu, \
  Vsetup_Plot, IEsel,Esel

  global Curr,EmitX,EmitV,BetaH,BetaV,SigE,Disph,Dispph,Dispv,Disppv,L,N,Beffv,Beffh,Nharm, \
  Harm,Shift,nKvals,Kmin,Kmax,Nmin,Nmax,Mode,Nelec,modepin,modesphere,Nepho, \
  EphMin,EphMax,PinX,PinY,PinZ,PinW,PinH,NpinZ,NpinY,Step,Pherror, \
  Ifixseed,Kellip

  Vsetup_Beam = []
  for key in BeamPar:
    Vsetup_Beam.append([key,Dsetup[key]])
  #endfor

  Vsetup_Undu = []
  for key in UnduPar:
    Vsetup_Undu.append([key,Dsetup[key]])
  #endfor

  Vsetup_Brill = []
  for key in BrillPar:
    Vsetup_Brill.append([key,Dsetup[key]])
  #endfor

  Vsetup_Spec = []
  for key in SpecPar:
    Vsetup_Spec.append([key,Dsetup[key]])
  #endfor

  Vsetup_Plot = []
  for key in PlotPar:
    Vsetup_Plot.append([key,Dsetup[key]])
  #endfor

  Vsetup = []

  for vs in Vsetup_Beam:
    Vsetup.append(vs)
  #endif
  for vs in Vsetup_Undu:
    Vsetup.append(vs)
  #endif
  for vs in Vsetup_Brill:
    Vsetup.append(vs)
  #endif
  for vs in Vsetup_Spec:
    Vsetup.append(vs)
  #endif
  for vs in Vsetup_Plot:
    Vsetup.append(vs)
  #endif

  for key in list(Dsetup):
    sval = Dsetup[key][1]
    if type(sval) == str:
      scom = "global " + key + "; " + key + " = '" + sval + "'"
    else:
      scom = "global " + key + "; " + key + " = " + str(sval)
    #endif
    exec(scom)
  #endfor

#enddef _dvsetup()

def debug(s=''):
  if s: print("debug:",s)
#enddef debug()

def _get_spec():
  global Calculated_Spec, NcalcSpec
  global nsto,nflx,nfld,nbun,nfdp

  #nsto = ncread("nsto","x:y:z:iegam:egam:s0:s1:s2:s3","urad_phase.sto")
  nflx = ncread("nflx","iegam:egam:s0:s1:s2:s3","urad_phase.flx")
  if fexist("urad_phase.fdp"):
    nfdp = ncread("nfdp","x:y:z:iegam:egam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fdp")
  nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
  nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt","urad_phase.bun")

  Calculated_Spec = True

  Frun = open("pybrill_spec.run",'r')
  line = Frun.readline().strip().split()
  Frun.close()
  Run_pyBrill = int(line[0])

  nlist()
  print('\n Spectra read from Run',Run_pyBrill,'\n')

  NcalcSpec += 1

#enddef _get_spec()

def _calc_spec():
  global Calculated_Spec, NcalcSpec
  global nsto,nflx,nfld,nbun,nfdp

  _UpdateVars()

  cwd = os.getcwd()
  _write_urad_phase_nam()
  t0 = time.time()

  localtime = time.asctime( time.localtime(time.time()) )

  if platform.system() == 'Windows':
    print('\n',"Starting spectrum calculation with urad_phase_win32.exe")
    print(localtime)
    os.system(cwd + '\\..\\bin\\urad_phase_win32.exe')
  else:
    print('\n',"Starting spectrum calculation with urad_phase.exe")
    print(localtime)
    os.system(cwd + "/../bin/urad_phase.exe")
  #endif platform.system() == 'Windows'

  t1 = time.time()
  localtime = time.asctime( time.localtime(time.time()) )
  print(localtime)
  print(" Calculation time [sec]:",int(t1-t0+0.5))
  print(' Loading N-tuples\n')

  #nsto = ncread("nsto","x:y:z:iegam:egam:s0:s1:s2:s3","urad_phase.sto")
  nflx = ncread("nflx","iegam:egam:s0:s1:s2:s3","urad_phase.flx")
  nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
  nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt","urad_phase.bun")
  if fexist("urad_phase.fdp"):
    nfdp = ncread("nfdp","x:y:z:iegam:egam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fdp")
  nlist()
  print("\nFinished")

  Calculated_Spec = True

  Frun = open("pybrill_spec.run",'w')
  Frun.write(str(Run_pyBrill) + " " + str(time.time()) + '\n')
  Frun.close()

  NcalcSpec += 1
#enddef _calc_spec()

def _set_plot_spec():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global SetUp_Plot, Vsetup_Plot, LastSetUp_Plot, Dsetup

  setmarkersize(float(Vsetup_Plot[1][1][1]))
  setlinewidth(float(Vsetup_Plot[2][1][1]))
  setlinecolor(Vsetup_Plot[3][1][1])

#enddef _set_plot_spec()

def _SetUpIn_Plot(event,kvar):
  global SetUp_Plot, Vsetup_Plot, LastSetUp_Plot
  LastSetUp_Plot = [event,kvar]
#enddef _SetUpInPlot(event,kvar)

def _SetUpOut_Plot(event,kvar):
  global SetUp_Plot, Vsetup_Plot, LastSetUp_Plot,Dsetup
  ev = LastSetUp_Plot[0].widget
  val = ev.get()
  try:
    if len(val.split('.')) > 1:
      val = float(val)
    else:
      val = int(val)
    #endif
  except:
    pass
  #entry
  Vsetup_Plot[kvar][1][1] = val
  vs = Vsetup_Plot[kvar]
  Dsetup[vs[0]] = vs[1]
  _set_plot_spec()
#enddef _SetUpOutPlot(event,kvar)

def _closeSetUp_Plot():

  global SetUp_Plot, Vsetup_Plot, LastSetUp_Plot, Dsetup

  if LastSetUp_Plot:
    ev = LastSetUp_Plot[0].widget
    kvar = LastSetUp_Plot[1]
    val = ev.get()
    try:
      if len(val.split('.')) > 1:
        val = float(val)
      else:
        val = int(val)
      #endif
    except: pass
    Vsetup_Plot[kvar][1][1] = val
    vs = Vsetup_Plot[kvar]
    Dsetup[vs[0]] = vs[1]
  else:
    #_set_plot_spec()
    for vs in Vsetup_Plot:
      Dsetup[vs[0]] = vs[1]
    #endfor
  #endif LastSetup

  _set_plot_spec()

  SetUp_Plot.destroy()

#def _closeSetUp_Plot(win)

def _vsetup_plot_ini():
  global SetUp_Plot, Vsetup_Plot, LastSetUp_Plot, Dsetup

  LastSetUp_Plot = 0
  SetUp_Plot = 0

  Dsetup['Mode3d'] = ["Mode3d [boxes, surf, scat2d, scat3d, inter]",'scat2d']
  Dsetup['Markersize']= ["Markersize",5.]
  Dsetup['Linewidth']= ["Linewidth",2.]
  Dsetup['Linecolor'] = ["Linecolor",'b']

  Vsetup_Plot = []
  for key in PlotPar:
    Vsetup_Plot.append([key,Dsetup[key]])
  #endfor

  _set_plot_spec()
  #print("ini:",Vsetup_Plot[0][1],Dsetup['Mode3d'])
#enddef _vsetup_plot_ini()

def _setup_plot():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  global SetUp_Plot, Vsetup_Plot

  SetUp_Plot = Toplevel()
  SetUp_Plot.title('Plotting Options for Spectra')
  SetUp_Plot.attributes('-topmost',1)
  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()
  SetUp_Plot.geometry('+' + str(int(xm+wm/2)) + '+' + str(int(ym+hm*0.7)))

  if not len(Vsetup_Plot): _vsetup_plot_ini()

  for i in range(len(Vsetup_Plot)):
    f = Frame(SetUp_Plot)
    flab = Label(f,text=Vsetup_Plot[i][1][0])
    fent =  Entry(f)
    fent.insert(1,Vsetup_Plot[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Plot(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Plot(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Plot(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Plot,text='Close',command=_closeSetUp_Plot)
  bClose.pack()

#enddef _setup_plot()

def _SetUpIn_Spec(event,kvar):
  global LastSetUp_Spec
  LastSetUp_Spec = [event,kvar]
#enddef _SetUpInSpec(event,kvar)

def _SetUpOut_Spec(event,kvar):
  global SetUp_Plot, Vsetup_Plot, LastSetUp_Plot
  global Vsetup_Spec, LastSetUp_Spec,Dsetup
  ev = LastSetUp_Spec[0].widget
  val = ev.get()
  if len(val.split('.')) > 1:
    val = float(val)
  else:
    val = int(val)
  #endif
  Vsetup_Spec[kvar][1][1] = val
  vs = Vsetup_Spec[kvar]
  Dsetup[vs[0]] = vs[1]
  _dvsetup()
#enddef _SetUpOutSpec(event,kvar)

def _closeSetUp_Spec():
  global SetUp_Spec,Vsetup_Spec, LastSetUp_Spec,Dsetup

  global Calculated_Spec, Calculated_Brill
  dso = deepcopy(Dsetup)

  if LastSetUp_Spec:
    ev = LastSetUp_Spec[0].widget
    kvar = LastSetUp_Spec[1]
    val = ev.get()
    if len(val.split('.')) > 1:
      val = float(val)
    else:
      val = int(val)
    #endif
    Vsetup_Spec[kvar][1][1] = val
    vs = Vsetup_Spec[kvar]
    Dsetup[vs[0]] = vs[1]
  else:
    for vs in Vsetup_Spec:
      Dsetup[vs[0]] = vs[1]
    #endfor
  #endif LastSetup

  _dvsetup()
  SetUp_Spec.destroy()

  if dso != Dsetup:
    Calculated_Spec = False
    Calculated_Beam = False
  #endif

#enddef _closeSetUp_Spec(win)

def _setup_spec():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  if not len(Vsetup_Spec):
    Vsetup_Spec = []
    for key in SpecPar:
      Vsetup_Spec.append([key,Dsetup[key]])
    #endfor
  #endif

  SetUp_Spec = Toplevel()
  SetUp_Spec.title('Spectra Calculation')
  SetUp_Spec.attributes('-topmost',1)

  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()

  SetUp_Spec.geometry('+' + str(int(xm+wm/3)) + '+' + str(int(ym+hm*0.3)))

  LastSetUp_Spec = 0

  for i in range(len(Vsetup_Spec)):
    f = Frame(SetUp_Spec)
    flab = Label(f,text=Vsetup_Spec[i][1][0])
    fent =  Entry(f)
    fent.insert(1,Vsetup_Spec[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Spec(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Spec(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Spec(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Spec,text='Close',command=_closeSetUp_Spec)
  bClose.pack()

#enddef _setup_spec()

def _SetUpIn_Brill(event,kvar):
  global LastSetUp_Brill
  LastSetUp_Brill = [event,kvar]
#enddef _SetUpIn_Brill(event,kvar)

def _SetUpOut_Brill(ev,kvar):
  global Vsetup_Brill, LastSetUp_Brill,Dsetup
  ev = LastSetUp_Brill[0].widget
  val = ev.get()
  if len(val.split('.')) > 1:
    val = float(val)
  else:
    val = int(val)
  #endif
  Vsetup_Brill[kvar][1][1] = val
  vs = Vsetup_Brill[kvar]
  Dsetup[vs[0]] = vs[1]
  _dvsetup()
#enddef _SetUpOut_Brill(event,kvar)

def _closeSetUp_Brill():
  global SetUp_Brill,Vsetup_Brill, LastSetUp_Brill,Dsetup

  global Calculated_Spec, Calculated_Brill
  dso = deepcopy(Dsetup)

  if LastSetUp_Brill:
    ev = LastSetUp_Brill[0].widget
    kvar = LastSetUp_Brill[1]
    val = ev.get()
    if len(val.split('.')) > 1:
      val = float(val)
    else:
      val = int(val)
    #endif
    Vsetup_Brill[kvar][1][1] = val
    vs = Vsetup_Brill[kvar]
    Dsetup[vs[0]] = vs[1]
  else:
    for vs in Vsetup_Brill:
      Dsetup[vs[0]] = vs[1]
    #endfor
  #endif LastSetup

  _dvsetup()
  SetUp_Brill.destroy()

  if dso != Dsetup:
    #Calculated_Spec = False
    Calculated_Beam = False
  #endif

#def _closeSetUp_Brill(win)

def _setup_brill():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  if not len(Vsetup_Brill):
    Vsetup_Brill = []
    for key in BrillPar:
      Vsetup_Brill.append([key,Dsetup[key]])
    #endfor
  #endif

  SetUp_Brill = Toplevel()
  SetUp_Brill.title('Brilliance Calculation')
  SetUp_Brill.attributes('-topmost',1)
  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()
  SetUp_Brill.geometry('+' + str(int(xm+wm/3.5)) + '+' + str(int(ym+hm*0.6)))

  LastSetUp_Brill = 0

  for i in range(len(Vsetup_Brill)):
    f = Frame(SetUp_Brill)
    flab = Label(f,text=Vsetup_Brill[i][1][0])
    fent =  Entry(f)
    fent.insert(1,Vsetup_Brill[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Brill(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Brill(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Brill(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Brill,text='Close',command=_closeSetUp_Brill)
  bClose.pack()

#enddef _setup_brill()

def _SetUpIn_Undu(event,kvar):
  global LastSetUp_Undu
  LastSetUp_Undu = [event,kvar]
#enddef _SetUpIn_Undu(event,kvar)

def _SetUpOut_Undu(event,kvar):
  global Vsetup_Undu, LastSetUp_Undu,Dsetup
  ev = LastSetUp_Undu[0].widget
  val = ev.get()
  if len(val.split('.')) > 1:
    val = float(val)
  else:
    val = int(val)
  #endif
  Vsetup_Undu[kvar][1][1] = val
  vs = Vsetup_Undu[kvar]
  Dsetup[vs[0]] = vs[1]
  _dvsetup()
#enddef _SetUpOutUndu(event,kvar)

def _closeSetUp_Undu():
  global SetUp_Undu,Vsetup_Undu, LastSetUp_Undu

  global Calculated_Spec, Calculated_Brill
  dso = deepcopy(Dsetup)

  if LastSetUp_Undu:
    ev = LastSetUp_Undu[0].widget
    kvar = LastSetUp_Undu[1]
    val = ev.get()
    if len(val.split('.')) > 1:
      val = float(val)
    else:
      val = int(val)
    #endif
    Vsetup_Undu[kvar][1][1] = val
    vs = Vsetup_Undu[kvar]
    Dsetup[vs[0]] = vs[1]
  else:
    for vs in Vsetup_Undu:
      Dsetup[vs[0]] = vs[1]
    #endfor
  #endif LastSetup

  _dvsetup()
  SetUp_Undu.destroy()

  if dso != Dsetup:
    Calculated_Spec = False
    Calculated_Beam = False
  #endif
#enddef _closeSetUp_Undu(win)

def _setup_undu():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global SetUp_Undu,Vsetup_Undu, LastSetUp_Undu

  if not len(Vsetup_Undu):
    Vsetup_Undu = []
    for key in UnduPar:
      Vsetup_Undu.append(Dsetup[key])
    #endfor
  #endif

  SetUp_Undu = Toplevel()
  SetUp_Undu.title('Undulator Parameters')
  SetUp_Undu.attributes('-topmost',1)

  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()
  SetUp_Undu.geometry('+' + str(int(xm+wm/3.8)) + '+' + str(int(ym+hm*0.6)))

  LastSetUp_Undu = 0

  for i in range(len(Vsetup_Undu)):
    f = Frame(SetUp_Undu)
    flab = Label(f,text=Vsetup_Undu[i][1][0])
    fent =  Entry(f)
    fent.insert(1,Vsetup_Undu[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Undu(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Undu(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Undu(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Undu,text='Close',command=_closeSetUp_Undu)
  bClose.pack()

#enddef _setup_undu()

def _SetUpIn_Beam(event,kvar):
  global LastSetUp_Beam
  LastSetUp_Beam = [event,kvar]
#enddef _SetUpIn_Beam(event,kvar)

def _SetUpOut_Beam(ev,kvar):
  global Vsetup_Beam, LastSetUp_Beam,Dsetup
  ev = LastSetUp_Beam[0].widget
  val = ev.get()
  if len(val.split('.')) > 1:
    val = float(val)
  else:
    val = int(val)
  #endif
  Vsetup_Beam[kvar][1][1] = val
  vs = Vsetup_Beam[kvar]
  Dsetup[vs[0]] = vs[1]
  _dvsetup()
#enddef _SetUpOut_Beam(event,kvar)

def _closeSetUp_Beam():
  global SetUp_Beam,Vsetup_Beam, LastSetUp_Beam,Dsetup

  global Calculated_Spec, Calculated_Brill

  dso = deepcopy(Dsetup)

  if LastSetUp_Beam:
    ev = LastSetUp_Beam[0].widget
    kvar = LastSetUp_Beam[1]
    val = ev.get()
    if len(val.split('.')) > 1:
      val = float(val)
    else:
      val = int(val)
    #endif
    Vsetup_Beam[kvar][1][1] = val
    vs = Vsetup_Beam[kvar]
    Dsetup[vs[0]] = vs[1]
  else:
    for vs in Vsetup_Beam:
      Dsetup[vs[0]] = vs[1]
    #endfor
  #endif LastSetup

  _dvsetup()
  SetUp_Beam.destroy()

  if dso != Dsetup:
    Calculated_Spec = False
    Calculated_Beam = False
  #endif
#enddef _closeSetUp_Beam(win)

def _setup_beam():

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  if not len(Vsetup_Beam):
    for key in BeamPar:
      Vsetup_Beam.append([key,Dsetup[key]])
      #endfor
  #endif

  SetUp_Beam = Toplevel()
  SetUp_Beam.attributes('-topmost',1)
  SetUp_Beam.title('Beam Parameters')

  xm = Wmaster.winfo_x()
  ym = Wmaster.winfo_y()
  wm = Wmaster.winfo_width()
  hm = Wmaster.winfo_height()

  SetUp_Beam.geometry('+' + str(int(xm+wm/4)) + '+' + str(int(ym+hm*0.4)))

  LastSetUp_Beam = 0

  for i in range(len(Vsetup_Beam)):
    f = Frame(SetUp_Beam)
    flab = Label(f,text=Vsetup_Beam[i][1][0])
    fent =  Entry(f)
    fent.insert(1,Vsetup_Beam[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Beam(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Beam(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Beam(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Beam,text='Close',command=_closeSetUp_Beam)
  bClose.pack()
#enddef _setup_beam()

def _writelastrun():

  global Dsetup

  _UpdateVars()

  fl = ".pybrill_last.dat"

  Fl = open(fl,"w")
  for key in list(Dsetup):
    Fl.write(key + ': ' + str(Dsetup[key]) + '\n')
  #endfor
  Fl.close()

#enddef _writelastrun()

def _readlastrun():
  global Dsetup,Vsetup,Vsetup_Beam,Vsetup_Spec,Vsetup_Brill,Vsetup_Undu, \
  Vsetup_Plot, IEsel,Esel

  global Curr,EmitX,EmitV,BetaH,BetaV,SigE,Disph,Dispph,Dispv,Disppv,L,N,Beffv,Beffh,Nharm, \
  Harm,Shift,nKvals,Kmin,Kmax,Nmin,Nmax,Mode,Nelec,modepin,modesphere,Nepho, \
  EphMin,EphMax,PinX,PinY,PinZ,PinW,PinH,NpinZ,NpinY,Step,Pherror, \
  Ifixseed,Kellip

  fl = ".pybrill_last.dat"

  if os.path.exists(fl):

    Fl = open(fl,"r")
    lines = Fl.readlines()
    Fl.close()

    global Dummy
    for line in lines:
      words = line.strip().split(':',1)
      key = words[0]
      scom = "global Dummy; Dummy=" + words[1] + "[0]"
      exec(scom)
      vs = Dummy
      scom = "global Dummy; Dummy=" + words[1] + "[1]"
      exec(scom)
      sval = Dummy
      if type(sval) == str:
        scom = "global " + key + "; " + key + "='" + sval + "'"
      else:
        scom = "global " + key + "; " + key + "=" + str(sval)
      #endif
      exec(scom)
      Dsetup[key] = [vs,sval]
    #endfor

  #endif

#enddef _readlastrun

def pwplot(y,n,ftit):

    global L, nKvals, Kvals, b0, Kmin, Kmax, Kellip, N, Ebeam, Curr, EmitH, EmitV
    global BetaH, BetaV, SigE, Mode, KyzList

    x = Harm[n]
    plt.plot(x,y)

    fn = ftit + "_Harm_" + str(n) + ".dat"

    Fn = open(fn,'w')

    Fn.write("* Period-length [mm]" + str(L) + "\n")
    Fn.write("* Beam energy [GeV], Curr [A] " + str(Ebeam) + " " + str(Curr) + "\n")
    Fn.write("* Hori. and vert. Emittance [nm-rad] " +  str(EmitH) + " " + str(EmitV) + "\n")
    Fn.write("* Hori. and vert. Beta-functions [m] " + str(BetaH) + " " + str(BetaV) + "\n")
    Fn.write("* Rel. beam energy spread " + str(SigE) + "\n")
    Fn.write("* Columns: \n")
    Fn.write("* Energy photons Keff Beff " + str(SigE) + "\n")

    for j in range(nKvals):
        i = nKvals - j - 1
        line = str(j+1) + " " + str(x[i]) + " " + str(y[i])
        line +=  " " + str(Kvals[i]) + " " + str(b0[i])
        Fn.write(line + "\n")
    Fn.close()

    print("Written to " + fn)
#enddef pwplot(x,y,fn='pyBrill.dat')

def grid(alpha=0.5):
    global Fig, Ax, Grid
    if Grid:
        plt.grid(Grid)
        Ax.grid(which='minor',alpha=alpha)
    #endif Grid
#enddef grid()

def zoom(xmin,xmax,ymin,ymax):
    global Fig, Ax
    Ax.axis([xmin,xmax,ymin,ymax])

    grid()
    plt.show(block=False)
#enddef zoom()

def _pcohflux():
  global LastPlot; LastPlot = ['cohflux','']

  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
  global Fig,Ax, Curr

  if Calculated_Brill == False: _calc_brill()

  Fig = plt.gcf()
  Fig.clear()

  Nmin = max(Nmin,1)
  Nmax = max(Nmax,Nmin)

  for i in range(Nmin,Nmax+1,2):
    if i >= Nmin and i <= Nmax:
      FCmax = max(FCmax,max(FC[i]))
      pwplot(FC[i],i,'Coh_Flux')
    #endif
  #endfor

  Fig = plt.gcf()

  Ax = plt.gca()
  Ax.set_xscale('log')
  Ax.set_yscale('log')

  Ax.set_title("Coherent Flux")

  Ax.set_xlabel("photon enery [keV]")
  Ax.set_ylabel("N [1/s/mm$^{2}$/0.1%BW/" + str(Curr) + "A]")

  Ax.set_ylim(FCmax/1.e4,FCmax*2.)

  grid()
  plt.show(block=False)
#enddef _pcohflux()

def _pbrillflux():
  global LastPlot; LastPlot = ['brillflux','']

  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
  global Fig,Ax, Curr

  if Calculated_Brill == False: _calc_brill()

  Fig = plt.gcf()
  Fig.clear()

  Nmin = max(Nmin,1)
  Nmax = max(Nmax,Nmin)

  for i in range(Nmin,Nmax+1,2):
    if i >= Nmin and i <= Nmax:
      pwplot(FB[i],i,'Brilliant_Flux')
      FBmax = max(FBmax,max(FB[i]))
    #endif
  #endfor

  Fig = plt.gcf()

  Ax = plt.gca()
  Ax.set_xscale('log')
  Ax.set_yscale('log')

  Ax.set_title("Brilliant Flux")

  Ax.set_xlabel("photon enery [keV]")
  Ax.set_ylabel("N [1/s/mm$^{2}$/0.1%BW/" + str(Curr) + "A]")

  Ax.set_ylim(FBmax/1.e4,FBmax*2.)

  grid()
  plt.show(block=False)
#enddef _pbrillflux()

def _pflux():
  global LastPlot; LastPlot = ['flux','']

  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
  global Fig,Ax, Curr

  if Calculated_Brill == False: _calc_brill()

  Fig = plt.gcf()
  Fig.clear()

  Nmin = max(Nmin,1)
  Nmax = max(Nmax,Nmin)

  for i in range(Nmin,Nmax+1,2):

    if i >= Nmin and i <= Nmax:
      Fmax = max(Fmax,max(F[i]))
      pwplot(F[i],i,'Flux')
    #endif

    Fig = plt.gcf()

    Ax = plt.gca()
    Ax.set_xscale('log')
    Ax.set_yscale('log')

    Ax.set_title("Flux")

    Ax.set_xlabel("photon enery [keV]")
    Ax.set_ylabel("N [1/s/0.1%BW/" + str(Curr) + "A]")

    Ax.set_ylim(Fmax/1.e4,Fmax*2.)

    grid()
    plt.show(block=False)
#enddef _pflux()

def _pfluxden():
  global LastPlot; LastPlot = ['fluxden','']

  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
  global Fig,Ax, Curr

  if Calculated_Brill == False: _calc_brill()

  Fig = plt.gcf()
  Fig.clear()

  Nmin = max(Nmin,1)
  Nmax = max(Nmax,Nmin)

  for i in range(Nmin,Nmax+1,2):

    if i >= Nmin and i <= Nmax:
      FDmax = max(FDmax,max(FD[i]))
      pwplot(FD[i],i,'Flux-density')
    #endif

    Fig = plt.gcf()

    Ax = plt.gca()
    Ax.set_xscale('log')
    Ax.set_yscale('log')

    Ax.set_title("Flux-density")

    Ax.set_xlabel("photon enery [keV]")
    Ax.set_ylabel("N [1/s/mrad$^{2}$/0.1%BW/" + str(Curr) + "A]")

    Ax.set_ylim(FDmax/1.e4,FDmax*2.)

    grid()
    plt.show(block=False)
#enddef _pfluxden()

def _pbrill():
  global LastPlot; LastPlot = ['brill','']

  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
  global Fig,Ax, Curr, Vsetup

  if Calculated_Brill == False: _calc_brill()

  Fig = plt.gcf()
  Fig.clear()

  Nmin = max(Nmin,1)
  Nmax = max(Nmax,Nmin)

  for i in range(Nmin,Nmax+1,2):
    if i >= Nmin and i <= Nmax:
      Bmax = max(Bmax,max(B[i]))
      pwplot(B[i],i,'Brilliance')
    #endif
  #endfor

  Fig = plt.gcf()

  Ax = plt.gca()
  Ax.set_xscale('log')
  Ax.set_yscale('log')

  Ax.set_title("Brilliance")

  Ax.set_xlabel("photon enery [keV]")
  Ax.set_ylabel("N [1/s/mm$^{2}$/mrad$^{2}$/0.1%BW/" + str(Curr) + "A]")

  Ax.set_ylim(Bmax/1.e4,Bmax*2.)

  grid()
  plt.show(block=False)
#enddef _pbrill()

def _UpdateVars():
  global \
  clight1,cgam1,cq1,alpha1,dnull1,done1,sqrttwopi1,\
  emassg1,emasse1,echarge1,emasskg1,eps01,erad1,\
  grarad1,hbar1,hbarev1,hplanck1,pol1con1,pol2con1,\
  radgra1,rmu01,rmu04pi1,twopi1,pi1,halfpi1,wtoe1,gaussn1,ck934,\
  ecdipev,ecdipkev,g1const,g1max,h2const,h2max

  global Dsetup,Vsetup,Vsetup_Beam,Vsetup_Spec,Vsetup_Brill,Vsetup_Undu, \
  Vsetup_Plot, IEsel,Esel

  global Curr,EmitX,EmitV,BetaH,BetaV,SigE,Disph,Dispph,Dispv,Disppv,L,N,Beffv,Beffh,Nharm, \
  Harm,Shift,nKvals,Kmin,Kmax,Nmin,Nmax,Mode,Nelec,modepin,modesphere,Nepho, \
  EphMin,EphMax,PinX,PinY,PinZ,PinW,PinH,NpinZ,NpinY,Step,Pherror, \
  Ifixseed,Kellip

  global Calculated_Brill, LastSetup, Kellip, Calculated_Spec, NcalcSpec
  global L, nKvals, Kvals, b0, Kmin, Kmax, Kellip, N, Ebeam, Curr, \
  EmitH, EmitV, BetaH, BetaV, SigE, Mode, Espread

  dsetupold = deepcopy(Dsetup)
  _dvsetup()

  N = Dsetup['Nper'][1]
  L = Dsetup['Perlen'][1]
  Dsetup['Espread'] = Dsetup['SigE']

  Calculated_Brill = False
  Calculated_Spec = False

  if dsetupold == Dsetup:
    if NcalcSpec > 0: Calculated_Spec = True
    if NcalcBrill > 0: Calculated_Brill = True
  #endif

  if Calculated_Brill == True: return

  if Beffh != 0: Kellip = True
  else: Kellip == False

  if Kellip: KyxList = [1.0, 0.42, 0.32, 0.27, 0.24, 0.22]
  else: KyxList = 0

  dK = (Kmax-Kmin)/(nKvals-1)
  Kvals = np.arange(Kmin,Kmax+dK,dK)
  b0 = Kvals/(echarge1 * L /1000./(2.*pi1*emasskg1*clight1))

#enddef _UpdateVars()

def _calc_brill():

    global Calculated_Brill, LastSetup, Kellip, NcalcBrill
    global L, nKvals, Kvals, b0, Kmin, Kmax, Kellip, N, Ebeam, Curr, \
    EmitH, EmitV, BetaH, BetaV, SigE, Mode

    _UpdateVars()

    calc_brill(\
    L, nKvals, Kmin, Kmax, N, Ebeam, Curr, \
    EmitH, EmitV, BetaH, BetaV, SigE, Mode \
    )

    Calculated_Brill = True
    NcalcBrill += 1

#enddef calc

def _exit():

  import platform,sys
  global Run_pyBrill, Run_pyBrill_Time

  Frun = open("pybrill.run",'w')
  Frun.write(str(Run_pyBrill) + " " + str(time.time()) + '\n')
  Frun.close()

  _writelastrun()
  stat = os.system("kill " + str(os.getpid()))

  try: _writelastrun()
  except: pass

  if platform.system() == 'Windows':
    stat = os.system("taskkill /F /PID " + str(os.getpid()))
  else:
    stat = os.system("kill " + str(os.getpid()))
    #endif platform.system() == 'Windows'
  #enddef _exit()

#enddef _exit()

def _SetUpIn(event,kvar):
    global LastSetup
    LastSetup = [event,kvar]
#enddef _SetUpIn(event,kvar)

def _SetUpOut(ev,kvar):
    global Vsetup,Vsetup_Beam,Vsetup_Undu,Vsetup_Brill,Vsetup_Spec, \
    Vsetup_Cont, LastSetup

    ev = LastSetup[0].widget
    val = ev.get()
    if len(val.split('.')) > 1:
      val = float(val)
    else:
      val = int(val)
    #endif
    Vsetup[kvar][1] = val

    Quit("_SetUpOut is obsolete!")

    if kvar == 4:
        v = Vsetup[4][1]
        if type(v) == str:
            if v.lower() in ['true','1','yes','y','j','ja']: Kellip = True
            else: Kellip = False
        else: Kellip = bool(Kellip)

#enddef _SetUpOut(event,kvar)

def _closeSetUp():
    global Vsetup,Vsetup_Beam,Vsetup_Undu,Vsetup_Brill,Vsetup_Spec, \
    Vsetup_Cont,LastSetup, Kellip, KyxList

    if LastSetup:
      ev = LastSetup[0].widget
      kvar = LastSetup[1]
      val = ev.get()
      Vsetup[kvar][1] = val
      vs = Vsetup[kvar]
      Dsetup[vs[0]] = vs[1]
    else:
      for vs in Vsetup:
        Dsetup[vs[0]] = vs[1]
      #endfor
    #endif LastSetup

    v = Vsetup[4][1]
    if type(v) == str:
        if v.lower() in ['true','1','yes','y','j','ja']: Kellip = True
        else: Kellip = False
    else: Kellip = bool(Kellip)

    SetUp.destroy()
#def _closeSetUp(win)

def _closeSetup_Menu():
    global Vsetup,Vsetup_Beam,Vsetup_Undu,Vsetup_Brill,Vsetup_Spec, \
    Vsetup_Cont,LastSetup, Kellip, KyxList

    Setup_Menu.destroy()
#def _closeSetUp(win)

def _setup():


  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global L, nKvals, Kvals, b0, Kmin, Kmax, Kellip, N, Ebeam, Curr, \
  EmitH, EmitV, BetaH, BetaV, SigE, Mode, Nmax
  global MBrill, Myfont
  global SetUp, Setup_Menu, Vsetup,Vsetup_Beam,Vsetup_Undu,Vsetup_Brill, \
  Vsetup_Spec,Vsetup_Cont, LastSetup

  Quit("*** _setup() is obsolete!")
  SetUp = Toplevel()
  LastSetup = 0

  for i in range(len(Vsetup)):
    f = Frame(SetUp)
    flab = Label(f,text=Vsetup[i][0])
    fent =  Entry(f)
    fent.insert(1,Vsetup[i][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut(event,kvar))
    f.pack(fill='x')
  #ENDfor i in range(len(V))

  bClose = Button(SetUp,text='Ok',command=_closeSetUp)
  bClose.pack()

  v = Vsetup[4][1]
  if type(v) == str:
    if v.lower() == ['true','1','yes','ja']: Kellip = True
    else: Kellip = False
  else:
    Kellip = bool(Kellip)
  #endif

  if Kellip:
    print("\n***********************************************************")
    print("Attention: For elliptical undulators K is shift-dependend!")
    print("Thus K must actually set for each harmonic...")
    print("***********************************************************\n")
  #endif Kellip

#enddef _setup

def _showMenu(menu,name):

  global BSetup,BBrill,BSpec,ScreenWidth,ScreenHeight,WmainMaster,Modepin, \
  ModeSphere,mPlotSpec

  fs = int(Myfont[1])

  xw = WmainMaster.winfo_x()
  yw = WmainMaster.winfo_y() - 3*fs
  ScreenWidth = WmainMaster.winfo_screenwidth()
  ScreenHeight = WmainMaster.winfo_screenheight()

  if name == 'MSetup':
    x = xw + 19*fs
    y = yw + 37*fs
    MBrill.unpost()
    MSpec.unpost()
    menu.post(x,y)
  elif name == 'MBrill':
    x = xw + 24*fs
    y = yw + 40*fs
    MSetup.unpost()
    MSpec.unpost()
    menu.post(x,y)
  elif name == 'MSpec':
    x = xw + 33*fs
    y = yw + 37*fs
    MSetup.unpost()
    MBrill.unpost()
    menu.post(x,y)
    if Modepin != 0:
      mPlotSpec.entryconfig(1,foreground='gray')
    else:
      mPlotSpec.entryconfig(1,foreground='black')
    #endif
  #endif

#enddef _showMenu()

def _canbutpybrill(ev):
  global MBrill,MSpec,MSetup
  #print("_canbutpybrill",ev)
  MBrill.unpost()
  MSpec.unpost()
  MSetup.unpost()
#enddef _canbutpybrill(ev)

def window_set_title(Title='',fig=-1):
  global Fig, Tfig, Figman
  if type(fig) == int and fig == -1:
    fig = plt.gcf()
    Tfig = type(fig)
    Fig = fig
  #endif type(fig) == int and fig == -1
  Figman = plt.get_current_fig_manager()
  Figman.set_window_title(Title)
#enddef window_set_title()

#######################################################

Perlen = 50.
nKvals = 101
Kmin = 0.5
Kmax = 3.
Kellip = False
Nmin = 1
Nmax = 11
Nper = 100
Ebeam = 1.722
Curr = 0.1
EmitH = 4.4 #nm-rad
EmitV = 0.066
BetaH = 14. #m
BetaV = 3.4
SigE = 0.001
Disph = 0.0 #m
Dispph = 0.0 #rad
Dispv = 0.0 #m
Disppv = 0.0 #rad
Mode = 2 # Walker

Esel = 0.0
IEsel = 0

Beffv = 1.0
Beffh = 0.0
Nharm = 1
Harmonic = 100.0
Shift = 0.0

Dsetup = {}

Dsetup['Ebeam'] =  ["Beam energy [GeV]",Ebeam]
Dsetup['Curr'] =  ["Current [A]",Curr]
Dsetup['EmitH'] =  ["Hor. Emit. [nm-rad]",EmitH]
Dsetup['EmitV'] =  ["Ver. Emit. [nm-rad]",EmitV]
Dsetup['BetaH'] =  ["Hori. Beta function",BetaH]
Dsetup['BetaV'] =  ["Vert. Beta function",BetaV]
Dsetup['SigE'] =  ["Rel. energy spread",SigE]
Dsetup['Disph'] =  ["Hori. Dispersion [mm]",Disph]
Dsetup['Dispph'] =  ["Derivative of Hori. Dispersion [mrad]",Dispph]
Dsetup['Dispv'] =  ["Vert. Dispersion [mm]",Dispv]
Dsetup['Disppv'] =  ["Derivative of Vert. Dispersion [mrad]",Disppv]

Dsetup['Perlen'] = ["Period-length [mm]",Perlen]
Dsetup['Nper'] = ["Number of periods",Nper]
Dsetup['Beffv'] = ["B0 vert. [T]",Beffv]
Dsetup['Beffh'] = ["B0 hori.. [T]",Beffh]
Dsetup['Nharm'] = ["Nharm (>0: overwrites Beff)",Nharm]
Dsetup['Harmonic'] = ["Harmonic [eV]",Harmonic]
Dsetup['Shift'] = ["Shift [mm] (for spectra only)",Shift]

Dsetup['nKvals'] = ["Number of K values",nKvals]
Dsetup['Kmin'] = ["Kmin",Kmin]
Dsetup['Kmax'] = ["Kmax",Kmax]
Dsetup['Nmin'] = ["Lowest harmonic",Nmin]
Dsetup['Nmax'] = ["Highest harmonic",Nmax]
Dsetup['Mode'] = ["Brilliance Mode [-1,1,2,3]",Mode]

Nelec = 1
Modepin = 0
ModeSphere = 0
Nepho = 11
EphMin = 90.
EphMax = 110.
PinW = 1.0
PinH = 1.0
PinX = 10000.0
PinY = 0.0
PinZ = 0.0
NpinZ = 11
NpinY = 11
Step = 0.2
Pherror = 0.0
IFieldProp = 0
PinXprop = 0.0
PinWprop = 0.1
PinHprop = 0.1
NpinZprop = 21
NpinYprop = 21
Ifixseed = 0

Dsetup['Nelec'] = ["Nelec",Nelec]
Dsetup['Modepin'] = ["Monte-Carlo mode [0,1]",Modepin]
Dsetup['Nepho'] = ["Number of Photon Energies",Nepho]
Dsetup['EphMin'] = ["Min. Photon Energy [eV]",EphMin]
Dsetup['EphMax'] = ["Max. Photon Energy [eV]",EphMax]
Dsetup['PinX'] = ["X of PinHole [mm]",PinX]
Dsetup['PinY'] = ["Y of PinHole [mm]",PinY]
Dsetup['PinZ'] = ["Z of PinHole [mm]",PinZ]
Dsetup['PinW'] = ["Width of PinHole [mm]",PinW]
Dsetup['PinH'] = ["Height of PinHole [mm]",PinH]
Dsetup['NpinZ'] = ["Number of hori. points",NpinZ]
Dsetup['NpinY'] = ["Number of vert. points",NpinY]
Dsetup['ModeSphere'] = ["Arrange grid points on sphere [0,1]",ModeSphere]
Dsetup['Step'] = ["Tracking step size [mm]",Step]
Dsetup['Pherror'] = ["Phase errors",Pherror]

Dsetup['IFieldProp'] = ["Propagate radiation field back to origin",IFieldProp]
Dsetup['PinXprop'] = ["Long. Position of Plane [mm]",PinXprop]
Dsetup['PinWprop'] = ["Width of PinHole in Plane [mm]",PinWprop]
Dsetup['PinHprop'] = ["Height of PinHole in Plane [mm]",PinHprop]
Dsetup['NpinZprop'] = ["Number of hori. points in Plane",NpinZprop]
Dsetup['NpinYprop'] = ["Number of vert. points in Plane",NpinYprop]

Dsetup['Ifixseed'] = ["Fix Seeds [0,1]",Ifixseed]

_vsetup_plot_ini()
_dvsetup()

Fkn = []
F = []
FD = []
FC = []
FB = []
Qn = []
B = []

Lam = []
Harm = []
Sigr = []
Sigrp = []

Calculated_Brill = False
Calculated_Spec = False

MBrill = 0
Omenu = 0
NMBrill = 0
NOmenu = 0

Emin = -1.
Emax = -1.

Bmin = -1.
Bmax = -1.
Fmin = -1.
Fmax = -1.
FDmin = -1.
FDmax = -1.
FCmin = -1.
FCmax = -1.
FBmin = -1.
FBmax = -1.

KyxList = 0

Grid = True
LastSetup = 0

global Run_pyBrill, Run_pyBrill_Time, NcalcBrill, NcalcSpec
NcalcBrill = 0
NcalcSpec = 0
try:
  Frun = open("pybrill.run",'r')
  line = Frun.readline().strip().split()
  Frun.close()
  Run_pyBrill = int(line[0])
  Run_pyBrill_Time = float(line[1])
  Run_pyBrill += 1
except:
  Run_pyBrill = 1
  Run_pyBrill_Time = 0
#endtry

_set_uname()
_vsetup_plot_ini()
_readlastrun()
_dvsetup()

global iVerbose
iVerbose = 1

mpl.use('TkAgg')

window()

global WmainMaster
Wmain = plt.gcf()
WmainMaster = Wmain.canvas.toolbar.master
Wmaster = WmainMaster
#print(id(Wmaster))

Toolbar = Wmain.canvas.toolbar
Myfont = ('arial',13)

CanButpyBrill = Wmain.canvas.mpl_connect('button_press_event',_canbutpybrill)
CanKeySpec = plt.connect('key_press_event', _spec_key_press)

window_set_title("pyBrill")
#zone(1,1)

if get_mshwelcome() == False:
    mshwelcome("pyBrill",2023)

global FigMain,AxMain
FigMain = plt.gcf()
AxMain = plt.gca()

MBrill = Menu(Toolbar,tearoff=1,font=Myfont)
mPlot = Menu(MBrill,tearoff=1,font=Myfont)

NMBrill += 1
MBrill.add_command(label='Calculate', command=_calc_brill)
NMBrill += 1
MBrill.add_cascade(label='Plot', menu=mPlot)
##########

mPlot.add_command(label='Brilliance', command=_pbrill)
mPlot.add_command(label='Flux-density', command=_pfluxden)
mPlot.add_command(label='Flux', command=_pflux)
mPlot.add_command(label='Coherent Flux', command=_pcohflux)
mPlot.add_command(label='Brilliant Flux', command=_pbrillflux)

##########

global BSetup, MSetup

MSetup = Menu(Toolbar,tearoff=1,font=Myfont)
BSetup = Button(Toolbar,text='Set-Up',font=Myfont, \
command= lambda menu = MSetup, name ='MSetup' : _showMenu(menu,name))
BSetup.pack(side=LEFT)

MSetup.add_command(label='Beam', command=_setup_beam)
MSetup.add_command(label='Undulator', command=_setup_undu)
MSetup.add_command(label='Brilliance', command=_setup_brill)
MSetup.add_command(label='Spectra', command=_setup_spec)

global BBrill
BBrill = Button(Toolbar,text='Brilliance',font=Myfont, \
command= lambda menu = MBrill, name = 'MBrill': _showMenu(menu,name))
BBrill.pack(side=LEFT)

global BSpec, MSpec, MFlux

MSpec = Menu(Toolbar,tearoff=1,font=Myfont)
mPlotSpec = Menu(MSpec,tearoff=0,font=Myfont)

MSpec.add_command(label='Calculate', command=_calc_spec)
MSpec.add_command(label='Load Previous Run', command=_get_spec)
BSpec = Button(Toolbar,text='Spectra',font=Myfont, \
command= lambda menu = MSpec, name = 'MSpec': _showMenu(menu,name))
BSpec.pack(side=LEFT)

MSpec.add_command(label='Plotting Options', command=_setup_plot)
MSpec.add_cascade(label='Plot', menu=mPlotSpec)

MFlux = Menu(MSpec,tearoff=0,font=Myfont)
MDist = Menu(MSpec,tearoff=1,font=Myfont)
MProp = Menu(MSpec,tearoff=1,font=Myfont)

global MFd
MFd = Menu(MSpec,tearoff=1,font=Myfont)

global MElec
MElec = Menu(MFlux,tearoff=0,font=Myfont)

mPlotSpec.add_cascade(label='Flux', menu=MFlux)
mPlotSpec.add_cascade(label='Central Flux-density', menu=MFd)
mPlotSpec.add_cascade(label='Distributions in Pinhole', menu=MDist)
mPlotSpec.add_cascade(label='Distributions of Propagated Fields', menu=MProp)
mPlotSpec.add_cascade(label='Electrons', menu=MElec)

MElec.add_command(label="Zi-Zi'", command= lambda key='zizpi': _pElec(key))
MElec.add_command(label="Yi-Yi'", command= lambda key='yiypi': _pElec(key))
MElec.add_command(label="Zi-Yi", command= lambda key='ziyi': _pElec(key))
MElec.add_command(label="Zi'-Yi'", command= lambda key='zpiypi': _pElec(key))
MElec.add_command(label="E", command= lambda key='eel': _pElec(key))

global MFdStokes
MFdStokes = Menu(MFd,tearoff=0,font=Myfont)
MFd.add_cascade(label='Stokes', menu=MFdStokes)
MFdStokes.add_command(label='S0', command= lambda key='s0': _pFdSpec(key))
MFdStokes.add_command(label='S1', command= lambda key='s1': _pFdSpec(key))
MFdStokes.add_command(label='S2', command= lambda key='s2': _pFdSpec(key))
MFdStokes.add_command(label='S3', command= lambda key='s3': _pFdSpec(key))

MFluxStokes = Menu(MFlux,tearoff=0,font=Myfont)
MFlux.add_cascade(label='Stokes', menu=MFluxStokes)
MFluxStokes.add_command(label='S0', command= lambda key='s0': _pFluxSpec(key))
MFluxStokes.add_command(label='S1', command= lambda key='s1': _pFluxSpec(key))
MFluxStokes.add_command(label='S2', command= lambda key='s2': _pFluxSpec(key))
MFluxStokes.add_command(label='S3', command= lambda key='s3': _pFluxSpec(key))

#MFluxPola = Menu(MFlux,tearoff=0,font=Myfont)
#MFlux.add_cascade(label='Polarisation', menu=MFluxPola)
#MFluxPola.add_command(label='P0', command= lambda key='p0': _pFluxSpec(key))
#MFluxPola.add_command(label='P1', command= lambda key='p1': _pFluxSpec(key))
#MFluxPola.add_command(label='P2', command= lambda key='p2': _pFluxSpec(key))
#MFluxPola.add_command(label='P3', command= lambda key='p3': _pFluxSpec(key))

MDistStokes = Menu(MDist,tearoff=0,font=Myfont)
MDistFields = Menu(MDist,tearoff=0,font=Myfont)

MDist.add_cascade(label='Stokes', menu=MDistStokes)
MDist.add_cascade(label='Field Amplitudes', menu=MDistFields)
MDist.add_command(label='Power', command= lambda key='p': _pFdPin(key))
MDist.add_command(label='Select E_photon', command=_setup_esel)

MDistStokes.add_command(label='S0', command= lambda key='s0': _pFdPin(key))
MDistStokes.add_command(label='S1', command= lambda key='s1': _pFdPin(key))
MDistStokes.add_command(label='S2', command= lambda key='s2': _pFdPin(key))
MDistStokes.add_command(label='S3', command= lambda key='s3': _pFdPin(key))

MDistFields.add_command(label='Ey_real', command= lambda key='eyr': _pFdPin(key))
MDistFields.add_command(label='Ey_imag', command= lambda key='eyi': _pFdPin(key))
MDistFields.add_command(label='Ez_real', command= lambda key='ezr': _pFdPin(key))
MDistFields.add_command(label='Ez_imag', command= lambda key='ezi': _pFdPin(key))

MPropStokes = Menu(MProp,tearoff=0,font=Myfont)
MPropFields = Menu(MProp,tearoff=0,font=Myfont)

MProp.add_cascade(label='Stokes', menu=MPropStokes)
MProp.add_cascade(label='Field Amplitudes', menu=MPropFields)
MProp.add_command(label='Select E_photon', command=_setup_esel)

MPropStokes.add_command(label='S0', command= lambda key='s0': _pFdProp(key))
MPropStokes.add_command(label='S1', command= lambda key='s1': _pFdProp(key))
MPropStokes.add_command(label='S2', command= lambda key='s2': _pFdProp(key))
MPropStokes.add_command(label='S3', command= lambda key='s3': _pFdProp(key))

MPropFields.add_command(label='Ey_real', command= lambda key='EYR': _pFdProp(key))
MPropFields.add_command(label='Ey_imag', command= lambda key='EYI': _pFdProp(key))
MPropFields.add_command(label='Ez_real', command= lambda key='EZR': _pFdProp(key))
MPropFields.add_command(label='Ez_imag', command= lambda key='EZI': _pFdProp(key))

bExit = Button(Toolbar,text='Exit',font=Myfont,command=_exit)
bExit.pack(side=LEFT)

global S_Esel,S_IEsel
S_IEsel = StringVar()
S_Esel = StringVar()
S_IEsel.set(IEsel)
S_Esel.set(Esel)

#_setup_menu()
#_setup()
#_pbrill()
