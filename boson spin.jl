using ITensors


let
  N = 7
  mu=1
  U=30
  t=1
  dB=2
  alpha=0.5
  delta=0.85
  beta=0.02
  f=1
  sites = [isodd(n) ? siteind("Boson", n; dim=2, conserve_qns=true) : siteind("Qubit", n; conserve_qns=true) for n in 1:N]
  state = ["1","Up","1", "Dn","1", "Up","0"]
  psi0 = random_mps(sites,state)
  #psi0 = randomMPS(sites; linkdims = 20)
#--------------------------set sweeps----------
maxdim = [200, 400]
mindim = [200]
nsweeps=50
niter=[4]
cutoff = [1E-12]
noise=[1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,1E-11,0,0,0,0,0,0,0]
etol=1E-6
#-----------------set H------------------------------
terms = OpSum()
  for i=1:2:(N-2)
  terms += -alpha, "Adag",i,"Z",i+1,"A",i+2
  terms += -alpha, "A",i,"Z",i+1,"Adag",i+2
  for i=1:2:(N-2) 
    terms += -t, "Adag",i,"A",i+2
    terms += -t, "A",i,"Adag",i+2
  end
  for i=1:2:N
    terms += (U/2), "N",i,"N",i
    terms += - (U/2), "N",i
  end
  for i=1:2:N
    terms += -(mu), "N",i
  end
  for i=1:2:(N-2) 
    terms += (delta/2), "Z",i+1
  end
  for i=1:2:(N-2) 
    terms += beta, "X",i+1
  end
H = MPO(terms, sites)
#----------------------------strat calculate-----------------------
energy0,psi = dmrg(H,psi0,ishermitian=false;nsweeps,maxdim,noise,cutoff)
print(energy0)
end
end

