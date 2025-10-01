print("--- BRILL ipylogon.py ---\n")

# +PATCH,//BRILL/PYTHON
# +DECK,ipylogon,T=PYTHON.
#plt.style.use('seaborn-dark')

set_console_title("brillPython")

import sys
args=sys.argv
arg1 = args[1]
nargs = len(args)

import m_hbook as m
from m_hbook import *

ntuples = 1
global Nhead,Ntup,Nind,Nntup

if arg1 == "last":
  try:
    Farg = open("ipylogon.arg","r")
    argl = Farg.readlines()
    Farg.close()
    if len(argl):
      sys.argv = [sys.argv[0]]
      for arg in argl: sys.argv.append(arg.strip())
    #endif
  except: pass
#endif

args=sys.argv
arg1 = args[1]
nargs = len(args)

if not arg1 == "last" and nargs > 0:
  Farg = open("ipylogon.arg","w")
  for i in range(1,nargs): Farg.write(args[i] + "\n")
  Farg.close()
#endif

ntuples = 1
histos = 1

def brill():
    import brill as b
    m.Kecho = 0
    m.Kpdf = 0
    m.Kdump = 0
    optnstat()
    b.brill(50)
    optconsole(); #getconsole()
    #Quit()
    return
#enddef brill()

def flux():
    import flux as f
    m.Kecho = 0
    m.Kpdf = 0
    m.Kdump = 0
    m.Kstat = 0
    optnstat()
    f.flux()
    optconsole(); #getconsole()
    #Quit()
    return
#enddef flux()

if nargs > 1:
  
  print("--- Argument:",arg1," ---\n")
  
  if arg1 == "default":
    
    brill()
    #reakpoint()
    npaul=ncread("npaul","e:b13:b0636","b2_bend1.3T_b3_bend0.64T.dat")
    ninfo(npaul)

    null(0.5e-3,2.0e4,1.0e8,5.0e17)
    xtit = "E$_{ph}$ [eV]"
    ytit = "N$_{ph}$/s/0.1%BW/mm$^{2}$/mrad$^{2}$"
    txyz("Brilliances of Dipoles",xtit,ytit)
    
    npllrs("b1","e:brill",legend='BRILL, BESSY II, 1.3T')
    npllrs("b2","e:brill",legend=' ')
    
    npllgs("b3","e:brill",legend='BRILL, BESSY III, 0.636T')
    npllgs("b4","e:brill",legend=' ')
    
    setlinestyle('solid')
    npllbs(npaul,"e:b13",legend='Paul, BESSY II, 1.3T')
    npllcs(npaul,"e:b0636",legend='Paul, BESSY III, 0.636T')
    
    legend()
    
  if arg1 == "default1":
    
    #    default()
    fin = open('pyDipole.in','r')
    flines = fin.readlines()
    sdist = flines[7].split()[-1]
    sdist2 = flines[7].split()[-1] + '**2'
    fin.close()
    ndip = ncread("ndip","y:e:fd:pd","pyDipole.out")
    ninfo(ndip)
    zone(1,2)
    emin = ndip.e.min()
    emax = ndip.e.max()
    npll(ndip,"y/"+sdist+":fd*"+sdist2,"e=="+str(emin))
    print("\nRMS",g3(emin*1000),'ev:')
    nrms(ndip,"y/"+sdist+":fd*"+sdist2,"e=="+str(emin))
    txyz(g3(emin*1000) + " eV","vert. angle [mrad]","Flux-density")
    nex()
    npll(ndip,"y/"+sdist+":fd*"+sdist2,"e=="+str(emax))
    txyz(g3(emax*1000) + " eV","vert. angle [mrad]","Flux-density")
    gtit("Dipole BESSY II")
    print(g3(ndip.fd.max()))
#    y,fd = vsymxy(nemin.y,nemin.fd)
#    vplxy(y,fd)
#    vfwhm(y,fd)
    
  elif arg1 == "brill": 
    brill()
    ndip = ncread("ndip","y:e:fd:pd","pyDipole.out")
    ndip0 = ncread("ndip0","e:fd:f:h2:g1:powd:vsiz","pyDipole_0.out")
    nlist()
    
  elif arg1 == "flux": 
    optnstat()
    flux()
    ndip = ncread("ndip","y:e:fd:pd","pyDipole.out")
    ndip0 = ncread("ndip0","e:fd:f:h2:g1:powd:vsiz","pyDipole_0.out")
    nlist()
    
  elif arg1 == "dipole_vsiz":
    
    nvs = ncread("nvs","e:vs","dipole_vert_beam_size.dat")
    lilo()
    npll(nvs,"e:vs")
    txyz("Vertical Beam Size","E [keV]","Flux / Fluxdensity [mrad]")
    
    
  elif arg1 == "dipole":
    
    ndip = ncread("ndip","y:e:fd:pd","pyDipole.out")
    ndip0 = ncread("ndip0","e:fd:f:h2:g1:powd:vsiz","pyDipole_0.out")
    
    #ninfo("ndip")
    ninfo("ndip0")
    
    lolo()
    optnstat()

    zone(1,2)
    
    npll(ndip0,"e*1000:fd","fd>0")
    txyz("","E [eV]","Flux-density")
    print("Max. Flux-density",g5(ndip0.fd.max()),'\n')

    nex()
    npll(ndip0,"e*1000:f","f>0")
    txyz("","E [eV]","Flux")
    print("Max. Flux",g5(ndip0.f.max()),'\n')
    
    gtit("BESSY II Dipole (pyDipole.py)")
    pp("BESSY_II_Dipole_pyDipole.pdf")
    
    dipole(1.7,1.307,0.3)
    
  elif arg1 == "pyBrill" or arg1 == 'U' or arg1 == 'UE' or arg1 == 'W' or arg1 == 'WE':
    
    if (fexist('urad_phase.fld')):      
      nfld = ncread("nfld","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase.fld")
    #endif
    
    if (fexist('urad_phase_espread.fld')):      
      nflde = ncread("nflde","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase_espread.fld")
    #endif
    
    if fexist('urad_phase.fdp'):
      nfdp = ncread("nfdp","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase.fdp")
    
    if fexist('urad_phase_espread.fdp'):
      nfdpe = ncread("nfdpe","x:y:z:iegam:iebeam:egam:ebeam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz:g","urad_phase_espread.fdp")
      
    fwig = '/scheer/urad_phase.wig'
    if (fexist(fwig)):
      nwig = ncread("nwig","kpola:iz:iy:itz:ity:x:y:z:ty:tz:iegam:iebeam:egam:ebeam::ezr:ezi:eyr:eyi:wig:fdzy:fdtzty:g",fwig)
  
    if (fexist('/scheer/urad_phase_espread.wig')):
      nwge = ncread("nwge","kpola:iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:wig:eyr:eyi:ezr:ezi","/scheer/urad_phase_espread.wig")

    nlist()
    #reakpoint()
    
    sely0 = 'y==0'
    selz0 = 'z==0'
    
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
    ianalytic = float(words[6])
    
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
    iwigner = int(words[0])
    nosplineefold = int(words[1])
    nefold = int(words[2])
    
    fpin.close()
    
else:
    pass
#endif nargs > 1

if ntuples:
    for n in range(len(Nhead)):
        snam = Nhead[n][1]
        exec(snam + ' = nget("' + snam + '")')
    #endfor
#endif

if histos:
  for hh in H1head:
    snam = hh[0]
    exec(snam + ' = hget("' + snam + '")')
  #endfor
  for hh in H2head:
    snam = hh[0]
    exec(snam + ' = hget("' + snam + '")')
  #endfor
#endif
