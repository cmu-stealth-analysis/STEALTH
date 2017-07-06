import numpy as np
import ROOT
from scipy.integrate import quad

 Keep time
sw = ROOT.TStopwatch()
sw.Start()

def f_up(p1):
  return lambda st: 1./((st/13000.)**p1)

def f_nm(p1):
  return lambda st: 1./((st/13000.)**(p1*np.log(st)))

def f_dn(p1):
  p1 = np.float128(p1)/13000.
  return lambda st: 1./np.exp(p1*st)

def I_(st_min,st_max,p0,f_):
  I,err = quad(f_,st_min,st_max)
  return I*p0,err

p0_fits = []
p1_fits = []
p1err_fits = []
with open('fit_choice_params.dat','r') as f:
  for ijt in f:
    params = ijt.strip('\n').split(' ')
    p0_fits.append(np.float32(params[0]))
    p1_fits.append(np.float32(params[1]))
    p1err_fits.append(np.float32(params[2]))

print p0_fits
print p1_fits
print p1err_fits
#p0_fits = [1.29174E-04, 5.90186E-07, 2.60714E+03]
p0_up = p0_fits[0] 
p0_nm = p0_fits[1] 
p0_dn = p0_fits[2]

#p1_fits = [5.09093E+00, 1.03868E+00, 5.01548E+01] # up, nm, dn
p1_up = p1_fits[0]
p1_nm = p1_fits[1]
p1_dn = p1_fits[2]
#p1err_nm = 1.06890E-01
p1err_nm = p1err_fits[1]

# CHOICE OF FIT UNCERTAINTY
# WARNING: Use DoubleEG data sample not GJet!

#st_mins = [1200,1300,1400]
st_mins = [1500] # DataB
#st_mins = [1050] # DataA
for st_min in st_mins:
  I_up,_ = I_(st_min,np.inf,p0_up,f_up(p1_up))
  I_nm,_ = I_(st_min,np.inf,p0_nm,f_nm(p1_nm))
  I_dn,_ = I_(st_min,np.inf,p0_dn,f_dn(p1_dn))
  print "Choice of fit | st_min=%d: %f %+.2f/%+.2f %%"%(st_min,I_nm,100.*(I_up-I_nm)/I_nm,100*(I_dn-I_nm)/I_nm)
print "-------------------------------------------------------"

# BKG SHAPE UNCERAINTY
# WARNING: Use DoubleEG data sample not GJet!

def r_xx2nm(st_norm,p1_xx,p1_nm):
  #I_std_xx,_ = I_(st_norm,st_norm+100,1.,f_nm(p1_xx))
  #I_std_nm,_ = I_(st_norm,st_norm+100,1.,f_nm(p1_nm))
  I_std_xx,_ = I_(st_norm,st_norm+150,1.,f_nm(p1_xx))
  I_std_nm,_ = I_(st_norm,st_norm+150,1.,f_nm(p1_nm))
  return I_std_nm/I_std_xx
def r_xx2nm_y(st_norm,p1_xx,p1_nm):
  y_xx = f_nm(p1_xx)(st_norm+100)
  y_nm = f_nm(p1_nm)(st_norm+100)
  return y_nm/y_xx

#st_norms = [1000,1100,1200]
st_norms = [1300] #DataB
#st_norms = [900] #DataA
p1_stds = [p1_nm+p1err_nm, p1_nm, p1_nm-p1err_nm]
for st_min in st_mins:
  I_nm,_ = I_(st_min,np.inf,1.,f_nm(p1_nm))
  for st_norm in st_norms:
    for p1 in p1_stds:
      I_xx,_ = I_(st_min,np.inf,1.,f_nm(p1))
      I_xx *= r_xx2nm(st_norm,p1,p1_nm)
      #I_xx *= r_xx2nm_y(st_norm,p1,p1_nm)
      print r_xx2nm(st_norm,p1,p1_nm)
      print "Bkg shape | st_norm=%d | st_min=%d: %f (%+.2f %%)"%(st_norm,st_min,I_xx,100.*(I_xx-I_nm)/I_nm)
    print "-------------------------------------------------------"

# ST SCALING UNCERTAINTY
# WARNING: Use GJet sample not DoubleEG!

def r_jt2nm(st_norm,p1_jt,p1_nm):
  I_jt,_ = I_(st_norm,st_norm+200,1.,f_nm(p1_jt))
  I_nm,_ = I_(st_norm,st_norm+200,1.,f_nm(p1_nm))
  return I_nm/I_jt
def r_jt2nm_y(st_norm,p1_jt,p1_nm):
  y_jt = f_nm(p1_jt)(st_norm+100)
  y_nm = f_nm(p1_nm)(st_norm+100)
  return y_nm/y_jt

p1_nJt = []
with open('st_scaling_params.dat','r') as f:
  for ijt in f:
    params = ijt.strip('\n').split(' ')
    p1_nJt.append(np.float32(params[1]))

print p1_nJt
iref = 0 #2jet
#iref = 1 #3jet
#p1_nJt = [1.00870E+00, 1.20344E+00, 7.93254E-01, 9.01548E-01] #2,3,4,5+
for st_min in st_mins:
  for st_norm in st_norms:
    for ijt, p1 in enumerate(p1_nJt):
      I_nm,_ = I_(st_min,np.inf,1.,f_nm(p1_nJt[iref]))
      I_jt,_ = I_(st_min,np.inf,1.,f_nm(p1))
      I_jt *= r_jt2nm(st_norm,p1,p1_nJt[iref])
      print "ST | njt=%d | st_norm=%d | st_min=%d: %f (%+.2f %%)"%(ijt+2,st_norm,st_min,I_jt,100.*(I_jt-I_nm)/I_nm)
    print "-------------------------------------------------------"


sw.Stop()
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
