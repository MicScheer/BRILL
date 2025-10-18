
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

# +PATCH,//NTUPPLOT/PYTHON
# +KEEP,ntupplot,T=PYTHON.

# Begin of NtupPlot

def _exit(): Quit()

def ngui_key_press(ev):
  if ev.key in ['q', 'Q']: Quit()
#enddef ngui_key_press(ev)

def startup(sfile='ntupplot_startup.py'):
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  global WaveFilePrefix, WavesMode

  if get_mshwelcome() == False:
    mshwelcome("Ntup-Plot",2021)
  if WavesMode == 'WAVES' or WavesMode == 'WPLOT' or WavesMode == 'WSHOP':
    fcfg = 'waveplot.cfg'
  elif WavesMode == 'UNDUMAG':
    fcfg = 'undugui.cfg'
  else:
    fcfg = 'ntupplot.cfg'
  #endif

  print("\n")
  print("\nHints:\n------")
  print("If a file " + sfile + " exists, it will be executed at start.")
  print("If a file " + fcfg + " exists, it will used to set window parameters\nof the first windows.")
  print("To spline data, plot them with the spline option; \na N-tuple 'Nspline' will be created then.")
  print("To leave, use the 'Exit' button, or enter 'q' in the canvas,\nor enter 'quit()' or 'Ctrl+q' in the terminal.\n")

  if os.path.exists(sfile):
    Fst = open(sfile,'r')
    print('\nEvaluating ' + sfile+ ":\n")
    lines = Fst.readlines()
    l = 0
    for line in lines:
      l += 1
      line = line.strip()
      if line.upper() == 'EOF': break
      if len(line) and line[0] == '#': continue
      elif len(line) > 6 and line[:6] != 'print(':
        print(line)
      #print(str(l)+": "+line)
      exec(line)
    Fst.close()
  #endif not os.path.exists(sfile)

  WaveFilePrefix = 'NtupPlot_'
#enddef startup()

def _showMenu(menu):
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if WavesMode == 'WAVES' or WavesMode == 'WPLOT'  or WavesMode == 'WSHOP':
    _showMenuWave(menu)
    return
  #endif WavesMode

  x,y = NPLmaster.winfo_pointerxy()

  KmenuPosted = Nmenu.winfo_ismapped()
  KplotPosted = Nplot.winfo_ismapped()
  KoptPosted = Omenu.winfo_ismapped()

  if menu == Nmenu:
    if KoptPosted:
      Omenu.unpost()
      KoptPosted = 0
    #endif
    if KplotPosted:
      Nplot.unpost()
      KplotPosted = 0
    #endif KplotPosted

    if KmenuPosted:
      Nmenu.unpost()
      KmenuPosted = 0
    else:
      Nplot.unpost()
      #    Omenu.unpost()
      Nmenu.post(x-50,y-50-NNmenu*2*Fontsize)
      KmenuPosted = 1
    #endif

  elif menu == Nplot:

    if KmenuPosted:
      Nmenu.unpost()
      KmenuPosted = 0
    #endif

    if KoptPosted:
      Omenu.unpost()
      KoptPosted = 0
    #endif KplotPosted

    if KplotPosted:
      Nplot.unpost()
      KplotPosted = 0
    else:
      Nmenu.unpost()
      #    Omenu.unpost()
      Nplot.post(x-50,y-50-NNplot*2*Fontsize)
      KplotPosted = 1
    #endif

  elif menu == Omenu:

    if KmenuPosted:
      Nmenu.unpost()
      KmenuPosted = 0
    #endif
    if KplotPosted:
      Nplot.unpost()
      KplotPosted = 0
    #endif KplotPosted

    if KoptPosted:
      Omenu.unpost()
      KoptPosted = 0
    else:
      Omenu.post(x-50,y-50-NOmenu*2*Fontsize)
      KoptPosted = 1
    #endif
  #endif menu == Nmenu

#enddef _showMenu(menu)

def framelabentry(win,text,var,stvar,font,widlab,wident):
  stvar.set(var)
  f = Frame(win)
  l = Label(f,text=text,font=font, width=widlab)
  l.pack(side=LEFT)
  e = Entry(f,text=stvar,width=wident,justify=CENTER,font=font)
  e.pack(side=LEFT)
  f.pack(fill='x')
#enddef framelabentry()

def _nTopLevel(title='TopLevel',att='-topmost',attn=1):
  tl = Toplevel()
  tl.title(title)
  tl.attributes(att,attn)
  return tl
#enddef _nTopLevel

def _clFillColor():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------

  setfillcolor(S_nFillColor.get())
  WnFillColor.destroy()
#enddef _clRead()

def _cnFillColor():
  global WnFillColor
  WnFillColor.destroy()
#enddef _cnFillColor()

def _nFillColor():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  S_nFillColor.set(getfillcolor())

  WnFillColor = _nTopLevel('Fillcolor')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-300) + '+' + str(y)
  WnFillColor.geometry(sgeo)

  widlab = 18
  wident = 18

  framelabentry(WnFillColor,'Color',S_nFillColor.get(),S_nFillColor,MyFont,widlab,wident)

  fbot = Frame(WnFillColor)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnFillColor,width=widlab-2)
  bCancel.pack(side=LEFT)
  bClose = Button(fbot,text='Ok',command=_clFillColor)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nFillColor()

def _clText():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  DictText['Text'] = S_nText.get()
  DictText['X'] = S_nTextX.get()
  DictText['Y'] = S_nTextY.get()
  DictText['NDC'] = S_nTndc.get().upper()
  DictText['Angle'] = S_nAngle.get()
  DictText['Color'] = S_nTcolor.get().lower()
  DictText['Halign'] = S_nHalign.get().upper()
  DictText['Valign'] = S_nValign.get().upper()
  DictText['Size'] = S_nTsize.get()

  if  DictText['NDC'] == 'Y' or DictText['NDC'] == 'J' or \
  DictText['NDC'] == 'YES' or  \
  DictText['NDC'] == '1' or DictText['NDC'] == 'JA' or  \
  DictText['NDC'] == 'TRUE':
    DictText['NDC'] = 'yes'
  elif  DictText['NDC'] == 'N' or DictText['NDC'] == 'NO' or  \
  DictText['NDC'] == '0' or DictText['NDC'] == 'NEIN' or  \
  DictText['NDC'] == 'FALSE':
    DictText['NDC'] = 'no'
  #endif

  if  DictText['Halign'] == 'L' or DictText['Halign'] == 'LEFT':
    DictText['Halign'] = 'left'
  if  DictText['Halign'] == 'R' or DictText['Halign'] == 'RIGHT':
    DictText['Halign'] = 'right'
  if  DictText['Halign'] == 'C' or DictText['Halign'] == 'CENTER':
    DictText['Halign'] = 'center'

  if  DictText['Valign'] == 'C' or DictText['Valign'] == 'CENTER':
    DictText['Valign'] = 'center'
  if  DictText['Valign'] == 'T' or DictText['Valign'] == 'TOP':
    DictText['Valign'] = 'top'
  if  DictText['Valign'] == 'B' or DictText['Valign'] == 'BOTTOM':
    DictText['Valign'] = 'bottom'

  x = float(DictText['X'])
  y = float(DictText['Y'])
  siz =int(DictText['Size'])
  ang = float(DictText['Angle'])

  if DictText['NDC'] == 'yes':
    text(x,y,DictText['Text'],fontsize=siz,color=DictText['Color'],
         halign=DictText['Halign'], valign=DictText['Valign'],angle=ang)
  else:
    textWC(x,y,DictText['Text'],fontsize=siz,color=DictText['Color'],
         halign=DictText['Halign'], valign=DictText['Valign'],angle=ang)
  #endif DictText['NDC'] = 'yes'

  WnText.destroy()
#enddef _clText()

def _cnText():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------

  DictText = deepcopy(DictTextO)
  WnText.destroy()
#enddef _cnText()

def _nText():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  DictTextO = deepcopy(DictText)

  WnText = _nTopLevel('Text')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-250) + '+' + str(y-300)
  WnText.geometry(sgeo)

  widlab = 12
  wident = 32

  framelabentry(WnText,'Text',S_nText.get(),S_nText,MyFont,widlab,wident)
  framelabentry(WnText,'X',S_nTextX.get(),S_nTextX,MyFont,widlab,wident)
  framelabentry(WnText,'Y',S_nTextY.get(),S_nTextY,MyFont,widlab,wident)
  framelabentry(WnText,'Norm. X,Y',S_nTndc.get(),S_nTndc,MyFont,widlab,wident)
  framelabentry(WnText,'Angle',S_nAngle.get(),S_nAngle,MyFont,widlab,wident)
  framelabentry(WnText,'Hori. align.',S_nHalign.get(),S_nHalign,MyFont,widlab,wident)
  framelabentry(WnText,'Vert. align.',S_nValign.get(),S_nValign,MyFont,widlab,wident)
  framelabentry(WnText,'Size',S_nTsize.get(),S_nTsize,MyFont,widlab,wident)
  framelabentry(WnText,'Color',S_nTcolor.get(),S_nTcolor,MyFont,widlab,wident)

  fbot = Frame(WnText)
  #bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnText,width=widlab-2)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnText)
  #bCancel.pack(side=LEFT)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clText)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nText()

def _clDump():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  sfile = S_nFile.get()
  snam = S_nName.get()

  if nexists(snam) == 0:
    nError(snam + " not existing!")
    return
  #endif

  svar = S_nVars.get()

  S_nLastCom.set('ndump')

  ssel = S_nSelect.get().lower()
  if ssel[0] == 'n': ssel = ''

  sind = 'no'
  try:
    sind = str(S_nDumpInd.get()).lower()
  except: pass
  if sind == 'none' or sind == 'False' or sind == '0' or sind == 'n': sind = 'no'
  elif sind == 'True' or sind == '1' or sind == 'y': sind = 'yes'

  shead = 'no'
  try:
    shead = str(S_nDumpHead.get()).lower()
  except: pass
  if shead == 'none' or shead == 'False' or shead == '0' or shead == 'n': shead = 'no'
  elif shead == 'True' or shead == '1' or shead == 'y': shead = 'yes'

  nFile = S_nFile.get()
  if nFile == '': nFile = 'ntuple.dat'

  ndump(snam,svar,ssel,sfile,shead,sind)

  WnDump.destroy()
#enddef _clDump()

def _cnDump():
  global WnDump
  WnDump.destroy()
#enddef _cnDump()

def _nDump():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  WnDump = _nTopLevel('Dump')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-250) + '+' + str(y-220)
  WnDump.geometry(sgeo)

  widlab = 24
  wident = 32

  nNam = Nhead[-1][1]
  nid = GetIndexN(nNam)

  varlis = list(Ntup[nid].columns)
  slis = nlistcolon(varlis)
  svar =''
  for s in slis:
    svar += ":" + s
  #endfor
  svar = svar[1:]

  ssel = S_nSelect.get()
  if ssel == '': ssel = 'none'

  sind = 'no'
  try:
    sind = str(S_nDumpInd.get()).lower()
  except: pass
  if sind == 'none' or sind == 'False' or sind == '0' or sind == 'n': sind = 'no'
  elif sind == 'True' or sind == '1' or sind == 'y': sind = 'yes'

  shead = 'no'
  try:
    shead = str(S_nDumpHead.get()).lower()
  except: pass
  if shead == 'none' or shead == 'False' or shead == '0' or shead == 'n': shead = 'no'
  elif shead == 'True' or shead == '1' or shead == 'y': shead = 'yes'

  framelabentry(WnDump,'Ntuple',nNam,S_nName,MyFont,widlab,wident)
  nFile = S_nFile.get()
  if nFile == '': nFile = 'ntuple.dat'
  framelabentry(WnDump,'Variables',svar,S_nVars,MyFont,widlab,wident)
  framelabentry(WnDump,'File',nFile,S_nFile,MyFont,widlab,wident)
  framelabentry(WnDump,'Selection',ssel,S_nSelect,MyFont,widlab,wident)
  framelabentry(WnDump,'Header',shead,S_nDumpHead,MyFont,widlab,wident)
  framelabentry(WnDump,'Index',sind,S_nDumpInd,MyFont,widlab,wident)

  fbot = Frame(WnDump)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnDump)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clDump)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nDump()

def nError(errtxt='Error',mode='widget'):

  global NPLmaster, WError

  if mode == 'widget':

    WError = Toplevel()
    WError.title('Error')

    x,y = NPLmaster.winfo_pointerxy()
    sgeo = '+' + str(x) + '+' + str(y)

    WError.geometry(sgeo)
    WError.attributes('-topmost', 1)

    lerr = Label(WError,text=errtxt,font=MyFont)
    lerr.pack(fill=X)

    bClose = Button(WError,text='Ok',command=WError.destroy)
    bClose.pack(fill=X)

    NPLmaster.wait_window(WError)

  else:
    print("\n",errtxt,"\n")
  #endif mode == 'widget'

#enddef nError(errtxt='Error')

def _clRead():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  sfile = S_nFile.get()
  snam = S_nName.get()

  if not os.path.exists(sfile):
    nError(sfile + " not found!")
    return
  #endif not os.path.exists(sfile)

  if nexists(snam) == 0:
    nError(snam + " not existing!")
    return
  #endif

  try:
    head = int(S_nHeader.get())
  except:
    head = None
  #endtry

  snsep = S_nSep.get()
  if snsep == "none": snsep = ''

  nread(snam,S_nFile.get(),head,int(S_nSkipHead.get()), \
  int(S_nSkipFoot.get()),0,S_nComment.get(),snsep)

  S_nLastCom.set('nread')

  print(NL)
  ninfo(snam)

  WnRead.destroy()
#enddef _clRead()

def _cnRead():
  global WnRead
  WnRead.destroy()
#enddef _cnRead()

def _nRead():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  WnRead = _nTopLevel('Read')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-250) + '+' + str(y-220)
  WnRead.geometry(sgeo)

  widlab = 24
  wident = 32

  snsep = S_nSep.get()
  if snsep == '': snsep = 'none'

  nNam = Nhead[-1][1]
  framelabentry(WnRead,'Ntuple',nNam,S_nName,MyFont,widlab,wident)
  nFile = 'ntuple.dat'
  framelabentry(WnRead,'File',nFile,S_nFile,MyFont,widlab,wident)
  skiphead = 0
  framelabentry(WnRead,'N of header lines to skip',skiphead,S_nSkipHead,MyFont,widlab,wident)
  skipfoot = 0
  framelabentry(WnRead,'N of footer lines to skip',skipfoot,S_nSkipFoot,MyFont,widlab,wident)
  scom = '*'
  framelabentry(WnRead,'Comment character',scom,S_nComment,MyFont,widlab,wident)
  sep = ' '
  framelabentry(WnRead,'Column seperator',ssep,S_nSep,MyFont,widlab,wident)

  fbot = Frame(WnRead)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnRead)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clRead)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nRead()

def _clMerge():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if nexists(S_nName12.get()) == 1:
    nError(S_nName.get() + " already existing!")
    return
  #endif

  S_nLastCom.set('nmerge')

  WnMerge.destroy()

#enddef _clMerge()

def _cnMerge():
  global WnMerge
  Merge = deepcopy(MergeO)
  WnMerge.destroy()
#enddef _cnMerge()

def _nMerge():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  MergeO = deepcopy(Merge)

  WnMerge = _nTopLevel('Merge')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-120)
  WnMerge.geometry(sgeo)

  widlab = 14
  wident = 30

  nNam = 'ntup' + str(Nntup)
  framelabentry(WnMerge,'Name 1',nNam,S_nName,MyFont,widlab,wident)
  nVars = 'x:y'
  framelabentry(WnMerge,'Variables',nVars,S_nVars,MyFont,widlab,wident)

  fbot = Frame(WnMerge)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnMerge)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clMerge)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nMerge()

def _clCreate():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if nexists(S_nName.get()) == 1:
    nError(S_nName.get() + " already existing!")
    return
  #endif

  nt = ncre(S_nName.get(),S_nTit.get(),S_nVars.get())
  S_nLastCom.set('ncre')

  WnCreate.destroy()

#enddef _clCreate()

def _cnCreate():
  global WnCreate
  WnCreate.destroy()
#enddef _cnCreate()

def _nCreate():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  WnCreate = _nTopLevel('Create')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-120)
  WnCreate.geometry(sgeo)

  widlab = 14
  wident = 30

  nNam = 'ntup' + str(Nntup)
  framelabentry(WnCreate,'Name',nNam,S_nName,MyFont,widlab,wident)

  nTit = 'ntup' + str(Nntup)
  framelabentry(WnCreate,'Title',nTit,S_nTit,MyFont,widlab,wident)

  nVars = 'x:y'
  framelabentry(WnCreate,'Variables',nVars,S_nVars,MyFont,widlab,wident)

  fbot = Frame(WnCreate)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnCreate)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clCreate)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nCreate()

def _clNull():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  widlab = 10
  wident = 10

  xmin = float(S_nXmin.get())
  xmax = float(S_nXmax.get())

  ymin = float(S_nYmin.get())
  ymax = float(S_nYmax.get())

  zmin = float(S_nZmin.get())
  zmax = float(S_nZmax.get())

  if zmax > zmin:
    null3d(xmin,xmax,ymin,ymax,zmin,zmax)
    S_nLastCom.set('null3d')
  else:
    null(xmin,xmax,ymin,ymax)
    S_nLastCom.set('null')
  #endif

  WnNull.destroy()

#enddef _clNull()

def _cnNull():
  global WnNull
  WnNull.destroy()
#enddef _cnNull()

def _nNull():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  WnNull = _nTopLevel('Frame')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-150)
  WnNull.geometry(sgeo)

  widlab = 14
  wident = 30

  xmin = float(S_nXmin.get())
  framelabentry(WnNull,'Xmin',xmin,S_nXmin,MyFont,widlab,wident)
  xmax = float(S_nXmax.get())
  framelabentry(WnNull,'Xmax',xmax,S_nXmax,MyFont,widlab,wident)

  ymin = float(S_nYmin.get())
  framelabentry(WnNull,'Ymin',ymin,S_nYmin,MyFont,widlab,wident)
  ymax = float(S_nYmax.get())
  framelabentry(WnNull,'Ymax',ymax,S_nYmax,MyFont,widlab,wident)

  zmin = float(S_nZmin.get())
  framelabentry(WnNull,'Zmin',zmin,S_nZmin,MyFont,widlab,wident)
  zmax = float(S_nZmax.get())
  framelabentry(WnNull,'Zmax',zmax,S_nZmax,MyFont,widlab,wident)

  fbot = Frame(WnNull)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnNull)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clNull)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

  Nplot.unpost()

#enddef _nNull()

def _clTitle():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  widlab = 10
  wident = 20

  ttit = S_nTitT.get()
  xtit = S_nTitX.get()
  ytit = S_nTitY.get()
  ztit = S_nTitZ.get()

  if hasattr(Ax,'zaxis'):
    txyz(ttit,xtit,ytit,ztit)
    S_n3d.set('yes')
  else:
    txyz(ttit,xtit,ytit)
    S_n3d.set('no')
  #endif

  WnTitle.destroy()

#enddef _clTitle()

def _cnTitle():
  global WnTitle
  WnTitle.destroy()
#enddef _cnTitle()

def _nTitle():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  WnTitle = _nTopLevel('Axis Titles')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-150)
  WnTitle.geometry(sgeo)

  widlab = 14
  wident = 30

  ttit = S_nTitT.get()
  xtit = S_nTitX.get()
  ytit = S_nTitY.get()
  ztit = S_nTitZ.get()

  wident = max([wident,len(ttit),len(xtit),len(ytit),len(ztit)])

  framelabentry(WnTitle,'Global title',ttit,S_nTitT,MyFont,widlab,wident)
  framelabentry(WnTitle,'X title',ttit,S_nTitX,MyFont,widlab,wident)
  framelabentry(WnTitle,'Y title',ttit,S_nTitY,MyFont,widlab,wident)
  framelabentry(WnTitle,'Z title',ttit,S_nTitZ,MyFont,widlab,wident)

  fbot = Frame(WnTitle)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnTitle)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clTitle)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

  Nplot.unpost()

#enddef _nTitle()

def _clInfo():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if nexists(S_nName.get()) == 0:
    nError(snam + " not existing!")
    return
  #endif

  ninfo(S_nName.get())
  WnInfo.destroy()
#enddef _clInfo()

def _cnInfo():
  global WnInfo
  WnInfo.destroy()
#enddef _cnInfo()

def _nInfo():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  WnInfo = _nTopLevel('Info')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-50)
  WnInfo.geometry(sgeo)

  widlab = 10
  wident = 10

  nNam = Nhead[-1][1]
  framelabentry(WnInfo,'Name',nNam,S_nName,MyFont,widlab,wident)

  fbot = Frame(WnInfo)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnInfo)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clInfo)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nInfo()

def _clDelete():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if nexists(S_nName.get()) == 0:
    nError(snam + " not existing!")
    return
  #endif

  ndelete(S_nName.get())
  WnDelete.destroy()
#enddef _clDelete()

def _cnDelete():
  global WnDelete
  WnDelete.destroy()
#enddef _cnDelete()

def _nDelete():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  WnDelete = _nTopLevel('Delete')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-100)
  WnDelete.geometry(sgeo)

  widlab = 10
  wident = 10

  nNam = Nhead[-1][1]
  framelabentry(WnDelete,'Name',nNam,S_nName,MyFont,widlab,wident)

  fbot = Frame(WnDelete)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnDelete)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clDelete)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

#enddef _nDelete()

def _clStat():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  snam = S_nName.get()
  svars = S_nVars.get()

  ssel = S_nSelect.get()
  if ssel == 'none': ssel = ''

  if nexists(snam) == 0:
    nError(snam + " not existing!")
    return
  #endif

  nstat(snam,svars,ssel)

  WnStat.destroy()
#enddef _clStat()

def _cnStat():
  global WnStat
  WnStat.destroy()
#enddef _cnStat()

def _nStat():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  WnStat = _nTopLevel('Stat')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-320)
  WnStat.geometry(sgeo)

  widlab = 15
  wident = 15

  nNam = Nhead[-1][1]
  nid = GetIndexN(nNam)
  nhead = Nhead[nid]

  nvar = nhead[3]
  svar = nhead[4][0]
  if nvar > 1: svar += ":" + nhead[5][0]

  ssel = S_nSelect.get()
  if ssel == '': ssel = 'none'

  framelabentry(WnStat,'Name',nNam,S_nName,MyFont,widlab,wident)
  framelabentry(WnStat,'Variables',svar,S_nVars,MyFont,widlab,wident)
  framelabentry(WnStat,'Selection',ssel,S_nSelect,MyFont,widlab,wident)

  fbot = Frame(WnStat)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnStat)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clStat)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

  Nmenu.unpost()

#enddef _nStat()

def _clPlot():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------

  global Imarker,Iline

  snam = S_nName.get()
  svars = S_nVars.get()

  ssel = S_nSelect.get()
  if ssel == 'none': ssel = ''
  swei = S_nWeight.get()
  if swei == 'none': swei = ''

  scx = float(S_nScaleX.get())
  scy = float(S_nScaleY.get())
  scz = float(S_nScaleZ.get())
  sct = float(S_nScaleT.get())

  splopt = S_nPlopt.get()
  if splopt == '': splopt = Mode2d

  plotoptions(splopt)

  sisame = S_nIsame.get().lower()
  if sisame == 'yes' or sisame == 'y': isame = 1
  else: isame = 0

  sisort = S_nIsort.get().lower()
  if sisort == 'yes' or sisort == 'y': isort = 1
  else: isort = 0

  if Isame == 0 and isame == 1:
    if splopt == '!':
      splopt = 'same'
    else:
      splopt = 'same' + splopt
    #endif
  #endif

  sleg = S_nLegend.get()

  smarker = S_nMark.get()
  imarker = 0
  if yesno(smarker) == 'yes': imarker = 1

  sline = S_nLine.get()
  iline = 0
  if yesno(sline) == 'yes': iline = 1
  #endif

  if iline == 1:
    if splopt == '!':
      splopt = 'line'
    else:
      splopt = 'line' + splopt
    #endif
  #endif

  setfillcolor(S_nFillColor.get())

  scol = S_nColor.get()

  h = hget(snam)

  if type(h) == int and h == -1:
    if nexists(snam) == 0:
      nError(snam + " not existing!")
      return
    #endif
    nplot(snam,svars,ssel,swei,splopt,sleg,scx,scy,scz,sct,'','HnPlot',scol,isort)
  else:
    if scx == 1.0 and scy == 1.0 and scz == 1.0 and ssel == '' and swei == '':
      hplot(snam,splopt,legend=sleg)
    else:
      snamN = snam + "_N"
      nh = hcopn(snam,snamN,svars)
      nplot(snamN,svars,ssel,swei,splopt,sleg,scx,scy,scz,sct,'','HnPlot',scol,isort)
    #endif
  #endif

  WnPlot.destroy()
#enddef _clPlot()

def _cnPlot():
  global WnPlot
  WnPlot.destroy()
#enddef _cnPlot()

def _nPlot():
#---------------------------------------------------------------------------

  global NPLmain, NPLmaster, MyFont,Myfont, Nmenu, NNmenu, CanBut, CanKey, Toolbar, Fontsize, \
  WError, WnCreate, S_nName, S_nTit, S_nVars, WnList, WnInfo, WnStat, \
  WnRead, S_nFile, S_nHeader, S_nIndex, S_nPlotInd, S_nPlotHead,S_nDumpInd, S_nDumpHead, \
  S_nSkipHead, S_nSkipFoot, S_nComment, S_nSep, \
  WnPlot,WnDump,WnDelete,WnTitle, S_nSelect, S_nWeight, S_nScaleX, S_nScaleY, \
  S_nScaleZ, S_nScaleT,S_nLegend, S_nHisto, S_nLine, S_nMark,S_nColor, S_nPlopt,S_nIsort, S_nIsame, \
  S_nLineColor, S_nLineStyle, S_nMarkerColor, S_nMarkerStyle, S_nIsame, \
  WnNull, S_nXmin,S_nYmin,S_nZmin,S_nXmax,S_nYmax,S_nZmax, S_nLastCom, \
  S_nTitT,S_nTitX,S_nTitY,S_nTitZ,S_n3d,Omenu,NOmenu,Wmaster, \
  KmenuPosted,KplotPosted,KoptPosted, \
  WnText,S_nText,S_nAngle,S_nTcolor,S_nTsize,S_nTndc,S_nTextX,S_nTextY, \
  DictText,DictTextO, S_nHalign,S_nValign, WnFillColor, S_nFillColor, \
  WnMerge,S_nName2,S_nVars2,S_nName1,S_nVars1,S_nName12,S_nVars12,Merge,MergeO, \
  WavesMode

  global Nplot, NNplot

#---------------------------------------------------------------------------


  global FillColor

  print(WavesMode)
  if not len(Nhead):
    nError("  No Ntuple defined so far!  ")
    return
  #endif not len(Nhead)

  WnPlot = _nTopLevel('Plot')

  x,y = NPLmaster.winfo_pointerxy()
  sgeo = '+' + str(x-200) + '+' + str(y-320)
  WnPlot.geometry(sgeo)

  widlab = 15
  wident = 15

  nNam = Nhead[-1][1]
  nid = GetIndexN(nNam)
  nhead = Nhead[nid]

  if hasattr(Ax,'zaxis'):
    S_n3d.set('yes')
  else:
    S_n3d.set('no')
  #endif

  nvar = nhead[3]
  svar = nhead[4][0]
  if nvar > 1:
    if nNam == 'n10':
      svar += ":" + nhead[6][0]
    else:
      svar += ":" + nhead[5][0]
    #endif
  #endif

  if S_n3d.get() == 'yes' and nvar > 2 : svar += ":" + nhead[6][0]

  snsep = S_nSep.get()
  if snsep == "none": snsep = 'blank'

  ssel = S_nSelect.get()
  if ssel == '': ssel = 'none'

  swei = S_nWeight.get()
  if swei == '': swei = 'none'

  scx = 1.
  scy = 1.
  scz = 1.
  sct = 1.

  splopt = S_nPlopt.get()
  plotoptions(splopt)

  if Isame: same = 'yes'
  else: same = 'no'

  if S_nLastCom.get() == 'null' or S_nLastCom.get() == 'null3d': same = 'yes'

  isort = S_nIsort.get().lower()
  if isort == '': isort = 'no'

  scol = S_nColor.get()
  if scol == '': scol = 'default'

  smarker = S_nMark.get()
  if smarker == '': smark = 'no'
  sline = S_nLine.get()
  if sline == '': scol = 'yes'

  framelabentry(WnPlot,'Name',nNam,S_nName,MyFont,widlab,wident)
  framelabentry(WnPlot,'Variables',svar,S_nVars,MyFont,widlab,wident)
  framelabentry(WnPlot,'Selection',ssel,S_nSelect,MyFont,widlab,wident)
  framelabentry(WnPlot,'Weights',swei,S_nWeight,MyFont,widlab,wident)

  framelabentry(WnPlot,'Scaling of 1st var.',scx,S_nScaleX,MyFont,widlab,wident)
  framelabentry(WnPlot,'Scaling of 2sd var.',scy,S_nScaleY,MyFont,widlab,wident)
  framelabentry(WnPlot,'Scaling of 3rd var.',scz,S_nScaleZ,MyFont,widlab,wident)
  framelabentry(WnPlot,'Scaling of 4th var.',sct,S_nScaleT,MyFont,widlab,wident)

  framelabentry(WnPlot,'Line',sline,S_nLine,MyFont,widlab,wident)
  framelabentry(WnPlot,'Marker',smarker,S_nMark,MyFont,widlab,wident)
  framelabentry(WnPlot,'Coler',scol,S_nColor,MyFont,widlab,wident)
  framelabentry(WnPlot,'Fill Color',FillColor,S_nFillColor,MyFont,widlab,wident)

  framelabentry(WnPlot,'Same picture',same,S_nIsame,MyFont,widlab,wident)
  framelabentry(WnPlot,'Sort data',isort,S_nIsort,MyFont,widlab,wident)

  fbot = Frame(WnPlot)
  bCancel = Button(fbot,text='Cancel',font=MyFont,command=_cnPlot)
  bCancel.pack(side=LEFT,expand=TRUE,fill=X)
  bClose = Button(fbot,text='Ok',command=_clPlot)
  bClose.pack(side=LEFT,expand=TRUE,fill=X)
  fbot.pack(expand=TRUE,fill=X)

  Nplot.unpost()

#enddef _nPlot()

# End of NtupPlot
#=============================================================================

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

global Kdebug
# not zero: Debugging with debugbreak(...)
# Kdebug = 2: Run urad_phase_debug.exe under gdb
# Kdebug = 3: Run urad_phase_debug.exe under ddd + gdb
if fexist('.pybrill.debug'):
  f = open('.pybrill.debug')
  Kdebug = int(f.readline().strip())
else: Kdebug = 0
#endif

global ClearCanvas
ClearCanvas = 0
set_ClearCanvas(ClearCanvas)

def _clear_Canvas(kclear=0):
  global Dsetup, ClearCanvas
  if kclear or Dsetup['ClearCanvas'][1]:
    window_clear()
    showplot(False)
    ClearCanvas = 0
    Dsetup['ClearCanvas'][1] = ClearCanvas
    set_ClearCanvas(ClearCanvas)
#enddef

def debugbreak(s=''):
  global Kdebug
  if not Kdebug: return
  if s: print("Debugbreak:",s)
  breakpoint()
#enddef debugbreak()


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
Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, Ianalytic, \
Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,ScreenWidth, ScreenHeight,Vsetup_Plot

global Esel
global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


global EbeamMin, EbeamMax, dEbeam, nEfold

global Unamelist,Useed
global LastPlot; LastPlot = ['','']

global MShWelcome
MShWelcome = False

def get_mshwelcome():
  global MShWelcome
  return MShWelcome
#enddef

#reakpoint()
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
  global Esel,IEsel,S_Esel,S_IEsel,dE

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
  'Mthreads','Ebeam','Curr','Step','Nelec','Noranone','IfixPhase','Icohere','Ihbunch', \
  'Bunchlen','BunchCharge','Modebunch','PinX','PinY','PinZ','PinW','PinH', \
  'NpinY','NpinZ','Modepin','ModeSphere','Perlen','Shift','Nper','Nharm', \
  'Harm','Beffv','Beffh','Ianalytic','Nepho','EphMin','EphMax','Espread','BetaH', \
  'BetaV','EmitH','EmitV','Disph','Dispph','Dispv','Disppv','Modeph','Pherror', \
  'IFieldProp','PinXprop','PinWprop','PinHprop','NpinYprop','NpinZprop', \
  'IfixPhase','PhGshift','IWigner','NyTheWig','TheYWig','NzTheWig','TheZWig', \
  'nEfold','NoSplineEfold']

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

  Dsetup['IfixPhase'] = ['Random phase for e-',0]
  Dsetup['PhGshift'] = ['Global phase shift. 9999.: Hori. amplitude is zero for (PinX,PinY,PinZ)) ',9999.]
  Dsetup['GlobPhaseProp'] = ['Global phase shift. of prop. fields. 9999.: Hori. amplitude is zero for (PinX,PinY,PinZ)) ',0.0]

  Dsetup['IWigner'] = ['Calculate Wigner Distributions [-4,-3,-2,-1,0,1]',0]
  Dsetup['NyTheWig'] = ['Number of vert. angle steps',51]
  Dsetup['TheYWig'] = ['Vert. angle range',0.05]
  Dsetup['NzTheWig'] = ['Number of hori. angle steps',51]
  Dsetup['TheZWig'] = ['Hori. angle range',0.05]
  Dsetup['nEfold'] = ['Number of E-spread steps',0]
  Dsetup['NoSplineEfold'] = ['Suppress splines for E-spread folding',0]

#enddef _set_uname()

global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
global Nmin, Nmax,Emin,Emax,Bmin,Bmax,FDmin,FDmax,Fmin,Fmax,FBmin,FBmax,FCmin,FCmax
global MBrill, Omenu, NMBrill, NOmenu, Myfont, Mylabe_font_size, Toolbar, \
Calculated_Brill, Fig, Ax, Grid, Calculated_Spec

global MSetup,MBrill,MSpec
global Kellip

debugbreak('Main 1')

BeamPar = ['Ebeam','Curr','EmitH','EmitV','BetaH','BetaV','SigE', \
'Disph','Dispph','Dispv','Disppv']

UnduPar = ['Perlen','Nper','Beffv','Beffh','Ianalytic','Nharm','Harmonic','Shift']

BrillPar = ['nKvals','Kmin','Kmax','Nmin','Nmax','Mode']

SpecPar = ['Nelec','Modepin','Noranone','ModeSphere','Nepho','EphMin','EphMax','PinX', \
'PinY','PinZ','PinW','PinH','NpinZ','NpinY','Step','Pherror', \
'IFieldProp','PinXprop','PinWprop','PinHprop','NpinYprop','NpinZprop', \
'IfixPhase','PhGshift','IWigner','NyTheWig','TheYWig','NzTheWig','TheZWig','nEfold','NoSplineEfold', \
'Ifixseed']

PlotPar = ['Mode3d','Markersize','Linewidth','Linecolor','NxZones','NyZones','ClearCanvas']

global LastSetUp_Esel
LastSetUp_Esel = 0

def webeam(ev):

  global Canbeam, Wbeam

  global Ebeam,IEbeam,S_Ebeam,S_IEbeamsel,nEfold,dEbeam

  nEfold = int(Dsetup['nEfold'][1])
  EbeamMin = float(Dsetup['EbeamMin'][1])
  EbeamMax = float(Dsetup['EbeamMax'][1])

  Ebeam = ev.xdata

  if nEfold > 1:
    dEbeam = (EBeamMax-EBeamMin)/(nEfold-1)
  else:
    dEbeam = 0
  #endif

  if Ebeamsel <= 0.0:
    IEbeamsel = 1
    Ebeamsel = EphMin
  elif Ebeamsel > EbeamMax:
    IEbeamsel = nEfold
    Ebeamsel = EbeamMax
  #endif

  if nEfold > 1:
    IEbeamsel = int((Ebeamsel-EbeamMin)/dEbeam)+1
    if IEbeamsel <=0:
      IEbeamsel = 1
    elif IEbeamsel > nEfold:
      IEbeamsel = nEfold
    #endif
  else:
    IEbeamsel = 1
  #endif

  Ebeamsel = EbeamMin + (IEbeamsel-1)*dEbeam

  S_IEbeamsel.set(IEbeamsel)
  S_Ebeamsel.set(Esbeamel)

  Webeam.canvas.mpl_disconnect(CanWebeam)
  window_close()

#enddef webeam()

def _wesel(ev):

  global Wesel,CanWesel,Wbeam

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
  Ifixseed,Kellip, Ianalytic

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global CanWesel, Wesel,nEfold

  SetUp_Esel.destroy()

  window('Select Photon Energy','')
  if Calculated_Spec == False or nexist("nflx") == 0: _calc_spec()

  global Dsetup,Nepho,EphMax,EphMin,WmainMaster

  Nepho = int(Dsetup['Nepho'][1])
  nEfold = int(Dsetup['nEfold'][1])
  iefold = int(nEfold/2) + 1

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
    #endif
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  flab = Label(f,text="Index of E_photon",font=('arial',MyLabel_font_size))
  fent =  Entry(f,text=S_IEsel)
  flab.pack(side=LEFT)
  fent.pack(side=RIGHT)
  fent.bind('<FocusIn>',lambda event,kvar=1:_SetUpIn_Esel(event,kvar))
  fent.bind('<FocusOut>',lambda event,kvar=1:_SetUpOut_Esel(event,kvar))
  fent.bind('<Return>',lambda event,kvar=1:_SetUpOut_Esel(event,kvar))
  f.pack(fill='x')

  f = Frame(SetUp_Esel)
  flab = Label(f,text="E_photon",font=('arial',MyLabel_font_size))
  fent =  Entry(f,text=S_Esel,font=('arial',MyLabel_font_size))
  flab.pack(side=LEFT)
  fent.pack(side=RIGHT)
  fent.bind('<FocusIn>',lambda event,kvar=2:_SetUpIn_Esel(event,kvar))
  fent.bind('<FocusOut>',lambda event,kvar=2:_SetUpOut_Esel(event,kvar))
  fent.bind('<Return>',lambda event,kvar=2:_SetUpOut_Esel(event,kvar))
  f.pack(fill='x')

  bSel = Button(SetUp_Esel,text='Select from spectrum',font=MyLabel_font_size,command=_sel_Esel)
  bSel.pack()
  bClose = Button(SetUp_Esel,text='Close',font=MyLabel_font_size,command=_closeSetUp_Esel)
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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

def _pWigner(key='WzzZ',select=''):


  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  global LastPlot; LastPlot = ['Wigner',key]
  global Esel,IEsel,S_Esel,S_IEsel,dE

#  if Modepin != 0: return
  debugbreak('_pWigner')

  if Calculated_Spec == False or nexist("nwig") == 0: _calc_spec()

  keyu = key.upper()
  keyl = key.lower()

  if keyu == 'WHHH': keyu = 'WZZZ'
  elif keyu == 'WHHV': keyu = 'WZZY'
  elif keyu == 'WVVV': keyu = 'WYYY'
  elif keyu == 'WVVH': keyu = 'WYYZ'
  elif keyu == 'WHVH': keyu = 'WZYZ'
  elif keyu == 'WVHH': keyu = 'WYZZ'
  elif keyu == 'WHVV': keyu = 'WZYY'
  elif keyu == 'WVHV': keyu = 'WYZY'

  selgam = "iegam==" + str(IEsel)

  set_plot_params_3d()

  kstat = getstat()
  nxzones = Dsetup['NxZones'][1]
  nyzones = Dsetup['NyZones'][1]

  optnstat()

  nz = int(Dsetup['NpinZprop'][1])
  ny = int(Dsetup['NpinYprop'][1])
  pinw = float(Dsetup['PinWprop'][1])
  pinh = float(Dsetup['PinHprop'][1])
  ymin = -pinh/2.0
  ymax =  pinh/2.0
  zmin = -pinw/2.0
  zmax =  pinw/2.0
  nty = int(Dsetup['NyTheWig'][1])
  ntz = int(Dsetup['NzTheWig'][1])
  tz = float(Dsetup['TheZWig'][1])
  ty = float(Dsetup['TheYWig'][1])

  if ny > 1: dy = pinh/(ny-1)
  else: dy = pinh / 2.
  if nz > 1: dz = pinw/(nz-1)
  else: dz = pinw / 2.

  if nty > 1: dty = ty/(nty-1)
  else: dty = ty / 2.
  if ntz > 1: dtz = tz/(ntz-1)
  else: dtz = tz / 2.

  tymin = -ty/2.0
  tymax =  ty/2.0
  tzmin = -tz/2.0
  tzmax =  tz/2.0

  a = ' and '
  sizcut = "iz==" + str(int(nz/2)+1)
  siycut = "iy==" + str(int(ny/2)+1)
  sitzcut = "itz==" + str(int(ntz/2)+1)
  sitycut = "ity==" + str(int(nty/2)+1)

  debugbreak('_pWigner')
  nwig = nget("nwig")

  plopt = Vsetup_Plot[0][1][1]
  lwo = float(Vsetup_Plot[2][1][1])
  lco = Vsetup_Plot[3][1][1]
  mso = getmarkersize()

#  disttit3d = getaxistitledist3d()
#  distlab3d = getaxislabeldist3d()
  colorbarpad = getcolorbarpad()

#  if select.strip() == '':
#    if Esel < EphMin or Esel > EphMax or IEsel <= 0:
#      if Nepho > 1:
#        dE = (EphMax-EphMin)/(Nepho-1)
#      else:
#        dE = 0
#      #endif
#      IEsel = int((nwig.iegam.max()-nwig.iegam.min())/2+1)
#      Esel = EphMin + (IEsel-1)*dE
#      S_IEsel.set(IEsel)
#      S_Esel.set(Esel)
#      select = 'iegam == ' + str(IEsel)
#    #endif select != ''
#  #endif

  if keyu == 'WZZZ':
    hnam = 'HWIGzzZ'
    htit = 'Wzz in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==1"
  elif keyu == 'WZZY':
    hnam = 'HWIGzzY'
    htit = 'Wzz in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==1"
  elif keyu == 'WYYZ':
    hnam = 'HWIGyyZ'
    htit = 'Wyy in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==2"
  elif keyu == 'WYYY':
    hnam = 'HWIGyyY'
    htit = 'Wyy in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==2"
  elif keyu == 'WZYZ':
    hnam = 'HWIGzyZ'
    htit = 'Wzy in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==3"
  elif keyu == 'WZYY':
    hnam = 'HWIGzyY'
    htit = 'Wzy in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==3"
  elif keyu == 'WYZZ':
    hnam = 'HWIGyzZ'
    htit = 'Wyz in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==4"
  elif keyu == 'WZYY':
    hnam = 'HWIGzyY'
    htit = 'Wzy in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==4"
  #endif keyu

  wtit = TeX_gamma + '/s/0.1' + ' %BW/mm$^{2}$/mrad$^{2}$/' + str(int(Curr*1000.+0.5)) + "mA"

  optstat(kstat)

  sel += ' and ' + selgam
  htit += htit + '  (' + str(Esel) + ' eV)'

  if keyu[3] == 'Z':

    h = hbook2(hnam,htit,
                 nz,zmin-dz/2.,zmax+dz/2.,
                 ntz,tzmin-dtz/2.,tzmax+dtz/2.,
                 overwrite=True)

    istat = nproj2(nwig,'z:tz','wig*1.0e-12',sel,1.0,1.0,1.0,nz,ntz,hnam)

    xtit = 'z [mm]'
    ytit = 'Theta_Z [mrad]'

    if h.y.min() < h.y.max():
      hplot2d(hnam,plopt,tit=htit,xtit=xtit,ytit=ytit,ztit='')
    else:
      npl(nwig,'z:tz',sel,'wig*1.0e-12')
      txyz(htit,xtit,ytit,' ')
    #endif


  else:

    h = hbook2(hnam,htit,
               ny,ymin-dy/2.,ymax+dy/2.,
               nty,tymin-dty/2.,tymax+dty/2.,
               overwrite=True)

    istat = nproj2(nwig,'y:ty','wig*1.0e-12',sel,1.0,1.0,1.0,ny,nty,hnam)

    xtit = 'y [mm]'
    ytit = 'Theta_Y [mrad]'

    if h.y.min() < h.y.max():
      iyty = 1
      hplot2d(h,plopt,tit=htit,xtit=xtit,ytit=ytit,ztit='')
    else:
      npl(nwig,'y:ty',sel,'wig*1.0e-12')
      txyz(htit,xtit,ytit,' ')
    #endif

  #endif keyu

  Dsetup['NxZones'][1] = nxzones
  Dsetup['NyZones'][1] = nyzones

#  setaxistitledist3d(disttit3d)
#  setaxislabeldist3d(distlab3d)
  setcolorbarpad(colorbarpad)

  global ClearCanvas
  ClearCanvas = 1
  Dsetup['ClearCanvas'][1] = ClearCanvas
  set_ClearCanvas(ClearCanvas)

#enddef _pWigner()

def _pWignerE(key='WzzZ',select=''):


  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  global LastPlot; LastPlot = ['Wigner',key]
  global Esel,IEsel,S_Esel,S_IEsel,dE,nEfold

  debugbreak('_pWignerE')
  iwigner = Dsetup['IWigner'][1]

  if Calculated_Spec == False or nexist("nwge") == 0: _calc_spec()

  keyu = key.upper()
  keyl = key.lower()

  if keyu == 'WHHH': keyu = 'WZZZ'
  elif keyu == 'WHHV': keyu = 'WZZY'
  elif keyu == 'WVVV': keyu = 'WYYY'
  elif keyu == 'WVVH': keyu = 'WYYZ'
  elif keyu == 'WHVH': keyu = 'WZYZ'
  elif keyu == 'WVHH': keyu = 'WYZZ'
  elif keyu == 'WHVV': keyu = 'WZYY'
  elif keyu == 'WVHV': keyu = 'WYZY'

  kstat = getstat()
  nxzones = Dsetup['NxZones'][1]
  nyzones = Dsetup['NyZones'][1]

  optnstat()

  nz = int(Dsetup['NpinZprop'][1])
  ny = int(Dsetup['NpinYprop'][1])
  pinw = float(Dsetup['PinWprop'][1])
  pinh = float(Dsetup['PinHprop'][1])
  ymin = -pinh/2.0
  ymax =  pinh/2.0
  zmin = -pinw/2.0
  zmax =  pinw/2.0
  nty = int(Dsetup['NyTheWig'][1])
  ntz = int(Dsetup['NzTheWig'][1])
  tz = float(Dsetup['TheZWig'][1])
  ty = float(Dsetup['TheYWig'][1])

  if ny > 1: dy = pinh/(ny-1)
  else: dy = pinh / 2.
  if nz > 1: dz = pinw/(nz-1)
  else: dz = pinw / 2.

  if nty > 1: dty = ty/(nty-1)
  else: dty = ty / 2.
  if ntz > 1: dtz = tz/(ntz-1)
  else: dtz = tz / 2.

  tymin = -ty/2.0
  tymax =  ty/2.0
  tzmin = -tz/2.0
  tzmax =  tz/2.0

  a = ' and '
  sizcut = "iz==" + str(int(nz/2)+1)
  siycut = "iy==" + str(int(ny/2)+1)
  sitzcut = "itz==" + str(int(ntz/2)+1)
  sitycut = "ity==" + str(int(nty/2)+1)

  debugbreak('_pWignerE')
  nwge = nget("nwge")

#  print("\nninfo('nwge'):\n")
#  ninfo('nwge')

  plopt = Vsetup_Plot[0][1][1]
  lwo = float(Vsetup_Plot[2][1][1])
  lco = Vsetup_Plot[3][1][1]
  mso = getmarkersize()

#  disttit3d = getaxistitledist3d()
#  distlab3d = getaxislabeldist3d()
  colorbarpad = getcolorbarpad()

  #reakpoint()
  nEfold = Dsetup['nEfold'][1]

  if nEfold < 2:
    if select.strip() == '':
      if Esel < EphMin or Esel > EphMax or IEsel <= 0:
        if Nepho > 1:
          dE = (EphMax-EphMin)/(Nepho-1)
        else:
          dE = 0
        #endif
        IEsel = int((nwge.iegam.max()-nwge.iegam.min())/2+1)
        Esel = EphMin + (IEsel-1)*dE
        S_IEsel.set(IEsel)
        S_Esel.set(Esel)
        select = 'iegam == ' + str(IEsel)
    #endif select != ''
  #endif

  if keyu == 'WZZZ':
    hnam = 'HWIGzzZE'
    htit = 'Wzz (E-folded) in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==1"
  elif keyu == 'WZZY':
    hnam = 'HWIGzzYE'
    htit = 'Wzz (E-folded) in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==1"
  elif keyu == 'WYYZ':
    hnam = 'HWIGyyZE'
    htit = 'Wyy (E-folded) in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==2"
  elif keyu == 'WYYY':
    hnam = 'HWIGyyYE'
    htit = 'Wyy (E-folded) in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==2"
  elif keyu == 'WZYZ':
    hnam = 'HWIGzyZE'
    htit = 'Wzy (E-folded) in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==3"
  elif keyu == 'WZYY':
    hnam = 'HWIGzyYE'
    htit = 'Wzy (E-folded) in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==3"
  elif keyu == 'WYZZ':
    hnam = 'HWIGyzZE'
    htit = 'Wyz (E-folded) in Z-Theta_Z Plane'
    if select.strip() == '': sel = siycut + a + sitycut
    else: sel = select + a +siycut + a + sitycut
    sel = sel + a + "kpola==4"
  elif keyu == 'WZYY':
    hnam = 'HWIGzyYE'
    htit = 'Wzy (E-folded) in Y-Theta_Y Plane'
    if select.strip() == '': sel = sizcut + a + sitzcut
    else: sel = select + a +sizcut + a + sitzcut
    sel = sel + a + "kpola==4"
  #endif keyu

  wtit = TeX_gamma + '/s/0.1' + ' %BW/mm$^{2}$/mrad$^{2}$/' + str(int(Curr*1000.+0.5)) + "mA"

  optstat(kstat)

  zone(2,1)

  if keyu[3] == 'Z':

    h = hbook2(hnam,htit,
                 nz,zmin-dz/2.,zmax+dz/2.,
                 ntz,tzmin-dtz/2.,tzmax+dtz/2.,
                 overwrite=True)

    istat = nproj2(nwge,'z:tz','wig*1.0e-12',sel,1.0,1.0,1.0,nz,ntz,hnam)

    xtit = 'z [mm]'
    ytit = 'Theta_Z [mrad]'

    if h.y.min() < h.y.max():
      hplot2d(hnam,plopt,tit=htit,xtit=xtit,ytit=ytit,ztit='')
    else:
#      setaxistitledist3d(disttit3d*2)
      setcolorbarpad(colorbarpad*2)
      npl(nwge,'z:tz:wig*1.0e-12:wig*1.0e-12',sel)
      txyz(htit,'  '+xtit,'\n\n'+ytit,' ')
    #endif

    zone(2,2,2,'same')
    npl(nwge,'z:wig*1.0e-12',sel + a + sitzcut)
    tit = 'Theta_Z = 0'
    xtit = 'z [mm]'
    txyz(tit,xtit,wtit)

    zone(2,2,4,'same')
    npl(nwge,'tz:wig*1.0e-12',sel + a + sizcut)
    tit = 'Z = 0'
    xtit = 'Theta_Z [mrad]'
    txyz(tit,xtit,wtit)

  else:
    h = hbook2(hnam,htit,
               ny,ymin-dy/2.,ymax+dy/2.,
               nty,tymin-dty/2.,tymax+dty/2.,
               overwrite=True)

    istat = nproj2(nwge,'y:ty','wig*1.0e-12',sel,1.0,1.0,1.0,ny,nty,hnam)

    xtit = 'y [mm]'
    ytit = 'Theta_Y [mrad]'

    if h.y.min() < h.y.max():
      iyty = 1
      hplot2d(h,plopt,tit=htit,xtit=xtit,ytit=ytit,ztit='')
    else:
#      setaxistitledist3d(disttit3d*2)
#      setaxislabeldist3d(1)
      setcolorbarpad(colorbarpad*2)
      npl(nwge,'y:ty:wig*1.0e-12:wig*1.0e-12',sel)
      txyz(htit,'  '+xtit,'\n\n'+ytit,' ')
    #endif

    zone(2,2,2,'same')
    npl(nwge,'y:wig*1.0e-12',sel + a + sitycut)
    tit = 'Theta_Y = 0'
    xtit = 'y [mm]'
    txyz(htit,xtit,wtit)

    zone(2,2,4,'same')
    npl(nwge,'ty:wig*1.0e-12',sel + a + siycut)
    tit = 'Y = 0'
    xtit = 'Theta_Y [mrad]'
    txyz(tit,xtit,wtit)

  #endif keyu

  Dsetup['NxZones'][1] = nxzones
  Dsetup['NyZones'][1] = nyzones

#  setaxistitledist3d(disttit3d)
#  setaxislabeldist3d(distlab3d)
  setcolorbarpad(colorbarpad)

  global ClearCanvas
  ClearCanvas = 1
  Dsetup['ClearCanvas'][1] = ClearCanvas
  set_ClearCanvas(ClearCanvas)

#enddef _pWignerE()

def _pFdSpec(key='s0'):

  global Dsetup
  global Vsetup_Beam, LastSetUp_Beam,  SetUp_Beam, \
  Vsetup_Undu, LastSetUp_Undu, SetUp_Undu, \
  Vsetup_Brill, LastSetUp_Brill, SetUp_Brill, Vsetup_Plot, \
  Vsetup_Spec, LastSetUp_Spec, SetUp_Spec, ScreenWidth, ScreenHeight

  global Mthreads,Step,Nelec,Noranone,Icohere,Ihbunch,Bunchlen, \
  Bunchcharge,Modebunch,PinX,PinY,PinZ,PinW,PinH,NpinY,NpinZ,modepin,modesphere, \
  Shift,Nper,Nharm,Harm,Beffv,Beffh,Nepho,EphMin,EphMax, \
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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

  nEfold = int(Dsetup['nEfold'][1])

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global Unamelist,Useed

  debugbreak('_write_urad_nam')

  try:
    shutil.copyfile("urad_phase.nam","urad_phase.nam.bck")
  except: pass

  Fnam = open("urad_phase.nam",'w')

  Fnam.write(" $uradphasen\n\n")

  for var in Unamelist:
    if var == 'Harm':
      Fnam.write("  " 'harm='+ str(Dsetup['Harmonic'][1]) + '          !' + Dsetup['Harmonic'][0] + '\n')
    else:
      val = Dsetup[var][1]
      if (\
      var.strip().upper()[0] == 'I' or \
      var.strip().upper()[0] == 'J' or \
      var.strip().upper()[0] == 'K' or \
      var.strip().upper()[0] == 'L' or \
      var.strip().upper()[0] == 'I' or \
      var.strip().upper()[0] == 'N') \
      :
        val = int(val)
        Dsetup[var][1] = val
      #endif
      Fnam.write("  " + var+'='+ str(val) + '          !' + Dsetup[var][0] + '\n')
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  Ifixseed,Kellip, Ianalytic
  global Ebeam, EbeamMin, EbeamMax, Espread

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

def __get_spec():

  global Calculated_Spec, NcalcSpec, SpecPar
  global nsto,nflx,nfld,nbun,nfdp,nwig,nwge
  global IWigner,nEfold,Esel,IEsel

#  debugbreak('_get_spec')

  if fexist("urad_phase_espread.flx"):
    nflx = ncread("nflx","iegam:iebeam:egam:ebeam:s0:s1:s2:s3","urad_phase_espread.flx")
  elif fexist("urad_phase.flx"):
      nflx = ncread("nflx","iegam:iebeam:egam:ebeam:s0:s1:s2:s3:g","urad_phase.flx")
  #endif

  if fexist("urad_phase_espread.fdp"):
    nfdp = ncread("nfdp","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase_espread.fdp")
  elif fexist("urad_phase.fdp"):
    nfdp = ncread("nfdp","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase.fdp")
  #endif

  if fexist("urad_phase_espread.fld"):
    nfld = ncread("nfld","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase_espread.fld")
  elif fexist("urad_phase.fld"):
    nfld = ncread("nfld","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase.fld")
  #endif

  if fexist("urad_phase.bun"):
    nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:","urad_phase.bun")

  iwig = 0

  fwige = 'urad_phase_espread.wig'
  fwig = 'urad_phase.wig'
  if fexist(fwige):
    nwig = ncread("nwig","kpola:iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:eyr:eyi:ezr:ezi:wig",fwige)
  elif (fexist(fwig)):
    nwig = ncread("nwig","kpola:iz:iy:itz:ity:x:y:z:ty:tz:iegam:iebeam:egam:ebeam:ezr:ezi:eyr:eyi:wig:g",fwig)
  #endif

  fpin = open("urad_phase.pin",'r')

  pin = fpin.readline().strip().split()
  npinz = int(pin[0])
  npiny = int(pin[1])
  pinw = float(pin[2])
  pinh = float(pin[3])

  pincen = fpin.readline().strip().split()
  pinx = pincen[0]
  piny = pincen[1]
  pinz = pincen[2]

  words = fpin.readline().strip().split()
  modepin = int(words[0])
  ifold = int(words[1])
  ifixphase = int(words[2])
  ifieldprop = int(words[3])
  nelec = int(words[4])
  ihbunch=int(words[5])
  Ianalytic = float(words[6])

  words = fpin.readline().strip().split()
  betah = float(words[0])
  emith = float(words[1])
  betav = float(words[2])
  emitv = float(words[3])
  espread = float(words[4])

  words = fpin.readline().strip().split()

  npinzprop = int(words[0])
  npinyprop = int(words[1])
  pinxprop = float(words[2])
  pinwprop = float(words[3])
  pinhprop = float(words[4])

  dzprop = pinwprop/(npinzprop-1)
  dyprop = pinhprop/(npinyprop-1)

  words = fpin.readline().strip().split()
  nzwig = int(words[0])
  nywig = int(words[1])
  tzwig = float(words[2])
  tywig = float(words[3])

  words = fpin.readline().strip().split()
  IWigner = int(words[0])
  nosplineefold = int(words[1])
  nEfold = int(words[2])

  fpin.close()

  Dsetup['IWigner'][1] = IWigner
  Dsetup['nEfold'][1] = nEfold

  _reset_mwigner()

  Calculated_Spec = True

  Frun = open("pybrill_spec.run",'r')
  line = Frun.readline().strip().split()
  Frun.close()
  Run_pyBrill = int(line[0])

  if IEsel <= 0:
    IEsel = int(nfld.iegam.max()/2) + 1
    Esel = int(nfld.iegam.max()/2) + 1
    S_IEsel.set(IEsel)
    EphMin = nfld.egam.min()
    EphMax = nfld.egam.max()
    Esel = (EphMax+EphMin)/2.0
    S_Esel.set(Esel)
  #endif

  nlist()
  print('\n Spectra read from Run',Run_pyBrill,'\n')

  NcalcSpec += 1

#enddef __get_spec()

def _get_spec():
  try:
    __get_spec()
  except:
    print('\n*** Failed to load previous run *** \n')
  #enddef
#enddef __get_spec()

def _calc_spec():

  global Calculated_Spec, NcalcSpec,SpecPar,Dsetup
  global nsto,nflx,nfld,nbun,nfdp,nwig,nwge

  debugbreak('calc_spec()')
  _UpdateVars()

  cwd = os.getcwd()
  _write_urad_phase_nam()
  t0 = time.time()

  localtime = time.asctime( time.localtime(time.time()) )

  debugbreak('_calc_spec')

  furph = open('.urad_phase.cal','w')
  furph.write('pyBrill\n')
  furph.close()

  if platform.system() == 'Windows':
    print('\n',"Starting spectrum calculation with urad_phase_win32.exe")
    print(localtime)
    #os.system(cwd + '\\..\\bin\\urad_phase_win32.exe')
    os.system('%BRILL%\\bin\\urad_phase_win32.exe')
  else:
    if Kdebug == 3:
      print('\n',"Starting spectrum calculation with urad_phase_debug.exe")
      print(localtime)
      os.system("ddd -gdb -geometry +20+10 --command ./startup.ddd $BRILL/bin/urad_phase_debug.exe")
    elif Kdebug == 2:
      print('\n',"Starting spectrum calculation with urad_phase_debug.exe")
      print(localtime)
      os.system("gdb --command ./startup.ddd $BRILL/bin/urad_phase_debug.exe")
    else:
      print('\n',"Starting spectrum calculation with urad_phase.exe")
      print(localtime)
      os.system("$BRILL/bin/urad_phase.exe")
    #endif
  #endif platform.system() == 'Windows'

  t1 = time.time()
  localtime = time.asctime( time.localtime(time.time()) )
  print(localtime)
  print(" Calculation time [sec]:",int(t1-t0+0.5))
  print(' Loading N-tuples\n')

  #nsto = ncread("nsto","x:y:z:iegam:egam:s0:s1:s2:s3","urad_phase.sto")

  ifieldprop = Dsetup['IFieldProp'][1]
  iwigner = Dsetup['IWigner'][1]
  nEfold = Dsetup['nEfold'][1]


  if fexist("urad_phase_espread.flx"):
    nflx = ncread("nflx","iegam:iebeam:egam:ebeam:s0:s1:s2:s3","urad_phase_espread.flx")
  elif fexist("urad_phase.flx"):
      nflx = ncread("nflx","iegam:iebeam:egam:ebeam:s0:s1:s2:s3","urad_phase.flx")
  #endif

  if fexist("urad_phase_espread.fdp"):
    nfdp = ncread("nfdp","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase_espread.fdp")
  elif fexist("urad_phase.fdp"):
    nfdp = ncread("nfdp","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase.fdp")
  #endif

  if fexist("urad_phase_espread.fld"):
    nfld = ncread("nfld","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase_espread.fld")
  elif fexist("urad_phase.fld"):
    nfld = ncread("nfld","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase.fld")
  #endif

  if fexist("urad_phase.bun"):
    nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt:g","urad_phase.bun")
  #endif

  fwige = 'urad_phase_espread.wig'
  fwig = 'urad_phase.wig'
  if fexist(fwige):
    nwig = ncread("nwig","kpola:iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:eyr:eyi:ezr:ezi:wig",fwige)
  elif (fexist(fwig)):
    nwig = ncread("nwig","kpola:iz:iy:itz:ity:x:y:z:ty:tz:iegam:iebeam:egam:ebeam:ezr:ezi:eyr:eyi:wig:g",fwig)
  #endif

  nlist()

  print("\nFinished")

  Calculated_Spec = True

  Frun = open("pybrill_spec.run",'w')
  Frun.write(str(Run_pyBrill) + " " + str(time.time()) + '\n')
  Frun.close()

  fpin = open("urad_phase.pin",'r')

  pin = fpin.readline().strip().split()
  npinz = int(pin[0])
  npiny = int(pin[1])
  pinw = float(pin[2])
  pinh = float(pin[3])

  pincen = fpin.readline().strip().split()
  pinx = pincen[0]
  piny = pincen[1]
  pinz = pincen[2]

  words = fpin.readline().strip().split()
  modepin = int(words[0])
  ifold = int(words[1])
  ifixphase = int(words[2])
  ifieldprop = int(words[3])
  nelec = int(words[4])
  ihbunch=int(words[5])
  Ianalytic = float(words[6])

  words = fpin.readline().strip().split()
  betah = float(words[0])
  emith = float(words[1])
  betav = float(words[2])
  emitv = float(words[3])
  espread = float(words[4])

  words = fpin.readline().strip().split()

  npinzprop = int(words[0])
  npinyprop = int(words[1])
  pinxprop = float(words[2])
  pinwprop = float(words[3])
  pinhprop = float(words[4])

  dzprop = pinwprop/(npinzprop-1)
  dyprop = pinhprop/(npinyprop-1)

  words = fpin.readline().strip().split()
  nzwig = int(words[0])
  nywig = int(words[1])
  tzwig = float(words[2])
  tywig = float(words[3])

  words = fpin.readline().strip().split()
  IWigner = int(words[0])
  nosplineefold = int(words[1])
  nEfold = int(words[2])

  fpin.close()

  Dsetup['IWigner'][1] = IWigner
  Dsetup['nEfold'][1] = nEfold

  _reset_mwigner()

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
  Dsetup['NxZones'] = ["NxZones",1]
  Dsetup['NyZones'] = ["NyZones",1]
  global ClearCanvas
  Dsetup['ClearCanvas'] = ["Clear Canvas",ClearCanvas]
  set_ClearCanvas(ClearCanvas)

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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
    flab = Label(f,text=Vsetup_Plot[i][1][0],font=('arial',MyLabel_font_size))
    fent =  Entry(f,font=('arial',MyLabel_font_size))
    fent.insert(1,Vsetup_Plot[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Plot(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Plot(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Plot(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Plot,text='Close',font=MyLabel_font_size,command=_closeSetUp_Plot)
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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

  SetUp_Spec.geometry('+' + str(int(xm+wm/3)) + '+' + str(int(ym+hm*0.05)))

  LastSetUp_Spec = 0

  for i in range(len(Vsetup_Spec)):
    f = Frame(SetUp_Spec)
    flab = Label(f,text=Vsetup_Spec[i][1][0],font=('arial',MyLabel_font_size))
    fent =  Entry(f,font=('arial',MyLabel_font_size))
    fent.insert(1,Vsetup_Spec[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Spec(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Spec(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Spec(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Spec,text='Close',font=MyLabel_font_size,command=_closeSetUp_Spec)
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
    flab = Label(f,text=Vsetup_Brill[i][1][0],font=('arial',MyLabel_font_size))
    fent =  Entry(f,font=('arial',MyLabel_font_size))
    fent.insert(1,Vsetup_Brill[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Brill(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Brill(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Brill(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Brill,text='Close',font=MyLabel_font_size,command=_closeSetUp_Brill)
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar


  global SetUp_Undu,Vsetup_Undu, LastSetUp_Undu

  debugbreak('_setup_undu')
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
    flab = Label(f,text=Vsetup_Undu[i][1][0],font=('arial',MyLabel_font_size))
    fent =  Entry(f,font=('arial',MyLabel_font_size))
    fent.insert(1,Vsetup_Undu[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Undu(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Undu(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Undu(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Undu,text='Close',font=MyLabel_font_size,command=_closeSetUp_Undu)
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

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
    flab = Label(f,text=Vsetup_Beam[i][1][0],font=('arial',MyLabel_font_size))
    fent =  Entry(f,font=('arial',MyLabel_font_size))
    fent.insert(1,Vsetup_Beam[i][1][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn_Beam(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut_Beam(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut_Beam(event,kvar))
    f.pack(fill='x')
  #endfor

  bClose = Button(SetUp_Beam,text='Close',font=MyLabel_font_size,command=_closeSetUp_Beam)
  bClose.pack()
#enddef _setup_beam()

def _writelastrun():

  global Dsetup

  debugbreak('_writelastrun')
  _UpdateVars()

  fl = ".pybrill_last.dat"
  flb = ".pybrill_last.dat.bck"

  if fexist(fl):
    copyfile(fl,flb)
  #endif

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
  Ifixseed,Kellip, Ianalytic

#  debugbreak('_readlastrun')
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
    Fn.write("* Index Energy Photons Keff Beff " + str(SigE) + "\n")

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
  Ifixseed,Kellip, Ianalytic

  global Calculated_Brill, LastSetup, Kellip, Calculated_Spec, NcalcSpec
  global L, nKvals, Kvals, b0, Kmin, Kmax, Kellip, N, Ebeam, Curr, \
  EmitH, EmitV, BetaH, BetaV, SigE, Mode, Espread, EbeamMin,EbeamMax, dEbeam

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

    debugbreak('calc_brill()')
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
  Disph,Dispph,Dispv,Disppv,Pherror,Ifixseed,Ianalytic

  global nsto,nfld,nflx,nbun,Esel
  global BeamPar,UnduPar,SpecPar,BrillPar,PlotPar



  global Fkn, F,FD, FC, FB, Qn, B, Harm, Lam, Sigr, Sigrp, KyxList
  global L, nKvals, Kvals, b0, Kmin, Kmax, Kellip, N, Ebeam, Curr, \
  EmitH, EmitV, BetaH, BetaV, SigE, Mode, Nmax
  global MBrill, Myfont, MyLabel_font_size
  global SetUp, Setup_Menu, Vsetup,Vsetup_Beam,Vsetup_Undu,Vsetup_Brill, \
  Vsetup_Spec,Vsetup_Cont, LastSetup

  Quit("*** _setup() is obsolete!")
  SetUp = Toplevel()
  LastSetup = 0

  for i in range(len(Vsetup)):
    f = Frame(SetUp)
    flab = Label(f,text=Vsetup[i][0],font=('arial',MyLabel_font_size))
    fent =  Entry(f,font=('arial',MyLabel_font_size))
    fent.insert(1,Vsetup[i][1])
    flab.pack(side=LEFT)
    fent.pack(side=RIGHT)
    fent.bind('<FocusIn>',lambda event,kvar=i:_SetUpIn(event,kvar))
    fent.bind('<FocusOut>',lambda event,kvar=i:_SetUpOut(event,kvar))
    fent.bind('<Return>',lambda event,kvar=i:_SetUpOut(event,kvar))
    f.pack(fill='x')
  #ENDfor i in range(len(V))

  bClose = Button(SetUp,text='Ok',font=Myfont,command=_closeSetUp)
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

debugbreak('Main 2')

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
Espread = SigE
Disph = 0.0 #m
Dispph = 0.0 #rad
Dispv = 0.0 #m
Disppv = 0.0 #rad
Mode = 2 # Walker
Noranone = 1

Esel = 0.0
IEsel = 0

Beffv = 1.0
Beffh = 0.0
Ianalytic = 0
Nharm = 1
Harmonic = 100.0
Shift = 0.0

#reakpoint()
Dsetup = {}

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
Dsetup['Beffh'] = ["B0 hori. [T]",Beffh]
Dsetup['Ianalytic'] = ["Use analytic formular for rad. field",Ianalytic]
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

IfixPhase = 0
PhGshift = 9999.
GlobPhaseProp = 9999.

IWigner = 0
NyTheWig = 31
TheYWig = 0.05
NzTheWig = 31
TheZWig= 0.05
nEfold = 0
NoSplineEfold = 0

EbeamMin = Ebeam*(1.0-3.0*Espread)
EbeamMax = Ebeam*(1.0+3.0*Espread)

if nEfold > 1:
  dEbeam = (EbeamMax-EbeamMin)/(nEfold-1)
else:
  dEbeam = 0
#endif

Dsetup['Nelec'] = ["Nelec",Nelec]
Dsetup['Modepin'] = ["Monte-Carlo mode [0,1]",Modepin]
Dsetup['Noranone'] = ["No randomization of first e- [0,1]",Noranone]
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

Dsetup['IfixPhase'] = ['Random phase for e-',IfixPhase]
Dsetup['PhGshift'] = ['Global phase shift. 9999.: Hori. amplitude is zero for (PinX,PinY,PinZ)) ',PhGshift]
Dsetup['GlobPhaseProp'] = ['Global phase shift. of prop. fields. 9999.: Hori. amplitude is zero for (PinX,PinY,PinZ)) ',GlobPhaseProp]
Dsetup['IWigner'] = ['Calculate Wigner Distributions [-4,-3,-2,-1,0,1]',IWigner]
Dsetup['NyTheWig'] = ['Number of vert. angle steps',NyTheWig]
Dsetup['TheYWig'] = ['Vert. angle range',TheYWig]
Dsetup['NzTheWig'] = ['Number of hori. angle steps',NzTheWig]
Dsetup['TheZWig'] = ['Hori. angle range',51,TheZWig]
Dsetup['Ebeam'] =  ["Beam energy [GeV]",Ebeam]
Dsetup['EbeamMin'] =  ["Min. beam energy [GeV]",EbeamMin]
Dsetup['EbeamMax'] =  ["Max. beam energy [GeV]",EbeamMax]
Dsetup['dEbeam'] =  ["dEbeam",dEbeam]
Dsetup['nEfold'] = ['Number of E-spread steps',0,nEfold]
Dsetup['NoSplineEfold'] = ['Suppress splines for E-spread folding',NoSplineEfold]

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

debugbreak('Main 3')

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
Myfont = ('arial',15)
#MyLabel_font_size = StringVar(value=15)
MyLabel_font_size = 12

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

global MWigner, NWignerMentries
NWignerMentries=0

MWigner = Menu(MSpec,tearoff=1,font=Myfont)

global MFd
MFd = Menu(MSpec,tearoff=1,font=Myfont)

global MElec
MElec = Menu(MFlux,tearoff=0,font=Myfont)

mPlotSpec.add_cascade(label='Flux', menu=MFlux)
mPlotSpec.add_cascade(label='Central Flux-density', menu=MFd)
mPlotSpec.add_cascade(label='Distributions in Pinhole', menu=MDist)
mPlotSpec.add_cascade(label='Distributions of Propagated Fields', menu=MProp)
mPlotSpec.add_cascade(label='Wigner Distributions', menu=MWigner)
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
MPropFields.add_command(label='Ez_imag', command= lambda key='EZI': _pFdProp(key))

IWigner = Dsetup['IWigner'][1]
nEfold = Dsetup['nEfold'][1]

def _reset_mwigner():

  global MWigner, NWignerMentries,nEfold

  #reakpoint()

  try: MWigner.delete('WzzZ')
  except: pass
  try: MWigner.delete('WzzY')
  except: pass

  try: MWigner.delete('WyyZ')
  except: pass
  try: MWigner.delete('WyyY')
  except: pass

  try: MWigner.delete('WzyZ')
  except: pass
  try: MWigner.delete('WzyY')
  except: pass

  try: MWigner.delete('WyzZ')
  except: pass
  try: MWigner.delete('WyzY')
  except: pass

  try: MWigner.delete('Select E_photon')
  except: pass

  NWignerMentries = 0

  if IWigner == -1 or IWigner > 0:
    MWigner.add_command(label='WzzZ', command= lambda key='WzzZ': _pWigner(key))
    MWigner.add_command(label='WzzY', command= lambda key='WzzY': _pWigner(key))
    NWignerMentries += 2
  #endif

  if IWigner == -2 or IWigner > 0:
    MWigner.add_command(label='WyyZ', command= lambda key='WyyZ': _pWigner(key))
    MWigner.add_command(label='WyyY', command= lambda key='WyyY': _pWigner(key))
    NWignerMentries += 2
  #endif

  if IWigner == -3 or IWigner > 0:
    MWigner.add_command(label='WzyZ', command= lambda key='WzyZ': _pWigner(key))
    MWigner.add_command(label='WzyY', command= lambda key='WzyY': _pWigner(key))
    NWignerMentries += 2
  #endif

  if IWigner == -4 or IWigner > 0:
    MWigner.add_command(label='WyzZ', command= lambda key='WyzZ': _pWigner(key))
    MWigner.add_command(label='WyzY', command= lambda key='WyzY': _pWigner(key))
    NWignerMentries += 2
  #endif

  MWigner.add_command(label='Select E_photon', command=_setup_esel)
  NWignerMentries += 1

#enddef

debugbreak('Main 4')

_reset_mwigner()

bExit = Button(Toolbar,text='Exit',font=Myfont,command=_exit)
bExit.pack(side=LEFT)

global S_Esel,S_IEsel
S_IEsel = StringVar()
S_Esel = StringVar()
S_IEsel.set(IEsel)
S_Esel.set(Esel)

def Load_Previous_Run():
  _get_spec()
#enddef

# Execute ntupplot_startup.py
startup()
