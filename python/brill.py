
# +PATCH,//BRILL/PYTHON
# +DECK,brill,T=PYTHON.

global WavesMode, Gdebug
WavesMode = 'BRILL'
Gdebug = 0

import pandas as pd
import os

#+seq,mshutil.
#+seq,plotglobal.
#+seq,hisglobal.
#+seq,m_hbook.

import m_hbook as m
from m_hbook import *

def brill(ntodo=100000):


  ncurve = 0

  emin = 1.0e30
  emax = 0
  brillmin = 1.0e30
  brillmax = 0
  istat = 0

  datei = open("brill.ebeam","r")
  line = datei.readline()
  words = line.split()
  ebeam = float(words[0])
  datei.close()

  datei = open("brill.curr","r")
  line = datei.readline()
  words = line.split()
  curr = float(words[0])
  datei.close()

  datei = open("brill_files.lis","r")
  files = datei.readlines()
  datei.close()

  ncurve = 0
  titold = ""

  NL = []
  TITS = []

  Colors = ['black','red','green','blue','yellow','magenta','cyan','gray','purple','brown','olive','pink','orange']

  colors = []
  nfiles = len(files)
  kcolor = 1
  ncolor = len(Colors)

  for i in range(ncolor): print(i+1,Colors[i])
  print("\n")

  for i in range(min(nfiles,ntodo+3)):
    brflfile = files[i]
    brflfile = brflfile.strip()
    split = brflfile.split('.')
    if split[0] == 'brfl' and split[1] != 'dat':
      ncurve += 1
      datei = open(brflfile,"r")
      line = datei.readline() # date
      line = datei.readline() # run
      line = datei.readline() # color
      words = line.split()
      color = int(words[1])
      if abs(color) > ncolor:
        kcolor += 1
        if kcolor > ncolor: kcolor = 1
        if color < 0: color = -kcolor
        else: color = kcolor
      #endif
#      print(ncolor,kcolor,color)
      colors.append(color)
      line = datei.readline()
      datei.close()
      tit = line.strip()
      if len(tit) > 1: titold = tit
      else: tit = titold
      if not len(tit) > 1:
        print("*** Warning: No title for ",file)
      nid = "b" + str(ncurve)
      nt = ncread(nid,"e:brill:flux:bright",brflfile)
      NL.append(nt)
      emn = nt["e"].min()
      if emn < emin : emin = emn
      emx = nt["e"].max()
      if emx > emax : emax = emx
      brillmn = nt["brill"].min()
      if brillmn < brillmin : brillmin = brillmn
      brillmx = nt["brill"].max()
      if brillmx > brillmax : brillmax = brillmx
      titnob = tit[3:]
      titnob = re.sub(" ","_",titnob)
      TITS.append(tit[3:])
      dst = titnob + "_" + brflfile
      print("\nCopying",brflfile,"to",dst,"\n")
      copyfile(brflfile,dst)


    #endif split[0] == 'brfl' and split[1] != 'dat'

  #endfor files



  emin = emin*0.9
  emax = emax*1.1
  #emin = 50.
  #emax = 4000.

  brillmin = brillmin/2.
  brillmax = brillmax*2.
  if brillmax > 1.0e18: brillmin = 1.e16
  #brillmax = 1.e18

  window("wBrill","800x1000+800+10")
  optnstat()
  null(emin,emax,brillmin,brillmax)
  grid()

  gtit = "Brilliance"
  gtit += ", " + "{:.3g}".format(ebeam) + " GeV"
  gtit += ", " + "{:.3g}".format(curr*1000) + " mA"

  xtit = "E$_{ph}$ [eV]"
  ytit = "N$_{ph}$/s/0.1%BW/mm$^{2}$/mrad$^{2}$"
  txyz(gtit,xtit,ytit)

  lolo()

  setmarkersize(4.)
  setlinewidth(2.)
  stitold = ""
  ntext = 0

  for i in range(min(ncurve,ntodo)):
    nt = NL[i]
    if colors[i] >= 1000:
      setlinestyle('dotted')
      col = Colors[abs(int(colors[i]/1000))-1]
    elif colors[i] < 0:
      setlinestyle('dashed')
      col = Colors[abs(colors[i])-1]
    else:
      setlinestyle('solid')
      col = Colors[abs(colors[i])-1]
    #endif
    setlinecolor(col)
    stit = TITS[i]
    nplot(nt,"e:brill",plopt='linesame')
    if len(stit) > 0 and stit != stitold:
      ntext += 1
      ytext = 0.05 + 0.035 * ntext
      text(0.75,ytext,stit,10,color=col,halign=LEFT)
      LineNDC(0.6,ytext+0.0,0.7,ytext+0.0,col)
    #endif len(stit) > 0 and stit != stitold
    stitold = stit
  #for i in range(ncurve)


  localtime = time.asctime( time.localtime(time.time()) )
  fig = plt.gcf()
  Tdate = fig.text(0.02,0.04,localtime,fontsize=9)

  w,h  = 16.,20.

  print("Figure saved to brill.pdf")
  pplot("brill.pdf",w,h)
  print("Figure saved to brill.png")
  pplot("brill.png",w,h)
  print("Figure saved to brill.jpg")
  pplot("brill.jpg",w,h)
  wans()
