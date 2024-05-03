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

def default():
    ndip = ncread("ndip","y:e:fd:pd","pyDipole.out")
    n100=ndip.query("e==100")
    y,fd = vsymxy(n100.y,n100.fd)
    vplxy(y,fd)
    pass
#enddef default()

def brill():
    import brill as b
    m.Kecho = 0
    m.Kpdf = 0
    m.Kdump = 0
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
    f.flux()
    optconsole(); #getconsole()
    #Quit()
    return
#enddef flux()

if nargs > 1:
  
  print("--- Argument:",arg1," ---\n")
  
  if arg1 == "default":
    
    #    default()
    ndip = ncread("ndip","y:e:fd:pd","pyDipole.out")
    n100=ndip.query("e==100")
    y,fd = vsymxy(n100.y,n100.fd)
    vplxy(y,fd)
    vfwhm(y,fd)
    
  elif arg1 == "brill": brill()
    
  elif arg1 == "flux": flux()
    
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
    npll(ndip0,"e:f","f>0")
    txyz("","E [keV]","Flux-density")
    
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
