import ROOT
import numpy as np
import os
import sys


file_name = sys.argv[1]
outfile = sys.argv[2]
reco = sys.argv[3]

if reco == "JV":
    option = 0
if reco == "OEC":
    option = 1

tree_name = "CdEvents"
tree_wp_name = "WpEvents"

prompt_min = 1500
prompt_max = 20000
delayed_min = 3700
delayed_max = 6000

t_min = 5
t_max = 2000 #us

FV = 17200 
z_cut = 11000
rho_cut = 3000
distance = 1500
veto_window = 2000 #(us)

coincidenze = []

#Per tree coinc
energyP_reco = []
energyD_reco = []
xP_reco = []
yP_reco = []
zP_reco = []
xD_reco = []
yD_reco = []
zD_reco = []
TT5_P = []
TT5_D = []
tipo_P = []
tipo_D = []
timeP = []
timeD = []

#Per tree veto
energyP_veto = []
energyD_veto = []
xP_veto = []
yP_veto = []
zP_veto = []
xD_veto = []
yD_veto = []
zD_veto = []
timeP_veto = []
timeD_veto = []
time_veto = []
muD_veto = []
FVD_veto = []
muP_veto = []
FVP_veto = []
multy_veto = []


#per leggere dati
energy = []
true_energy = []
time = []
sec = [] #per salvare fSec
nano = [] #per salvare fNanoSec
x = []
y = []
z = []
r = []
TtT5 = []
mu_tag = []
deltaT_muon = []
triggerType = []
tipo = []

energy_WP = []
sec_WP = []
nano_WP = []


print("Opening:" + file_name)
tfile = ROOT.TFile(file_name, "READ")
tree = tfile.Get(tree_name)
tree_wp = tfile.Get(tree_wp_name)
lt = tfile.Get("lt")
duration = lt.GetBinContent(1)

it = iter(tree)
for event in it:
    true_energy.append(event.OECRecoEnergy)
    sec.append(event.TimeStamp.GetSec())    #(float(event.fSec)) #secondi
    nano.append(event.TimeStamp.GetNanoSec()) #(float(event.fNanoSec)) #nanosecondi
    energy.append(float(event.NPE))
    TtT5.append(event.tailToTot500)
    tipo.append(event.OECTag)
    #mu_tag.append(event.muonTag)
    #deltaT_muon.append(event.deltaT_muon)
    triggerType.append([str(s) for s in event.TriggerType])
    if option == 0:
        x.append(event.JRecoX)
        y.append(event.JRecoY)
        z.append(event.JRecoZ)
    if option == 1:
        x.append(event.OECRecox)
        y.append(event.OECRecoy)
        z.append(event.OECRecoz)
        
sec = np.array(sec, dtype=np.float64)
nano = np.array(nano, dtype=np.float64)
energy = np.array(energy, dtype=np.float64)


it_wp = iter(tree_wp)
for event in it_wp:
    energy_WP.append(float(event.NPE)) 
    sec_WP.append(event.TimeStamp.GetSec())#(float(event.fSec))
    nano_WP.append(event.TimeStamp.GetNanoSec())  #(float(event.fNanoSec))

energy_WP = np.array(energy_WP, dtype=np.float64)
sec_WP = np.array(sec_WP, dtype=np.float64)
nano_WP = np.array(nano_WP, dtype=np.float64)
    

# muon veto ###########################################################################################################################################   
    
mu_tag = np.zeros(len(energy), dtype=np.int32) 
muWP_tag = np.zeros(len(energy_WP), dtype=np.int32)
deltaT_muon = [] #np.full(len(energy), -1, dtype=np.float64)
mu_time = []   
count_mu = 0
count_CD = 0
count_WP = 0
count_CD_WP = 0

for i in range(len(energy)):
    if energy[i] > 30000:
        wp_apo = energy_WP[abs((sec[i]-sec_WP)*1e9+(nano[i]-nano_WP)) < 500]
        
        if np.any(wp_apo > 400):
            mu_tag[i] = 1
            count_CD_WP += 1
            tag_apo = 1
            #print("muone CD+WP trovato")
        else:
            mu_tag[i] = 2
            count_CD += 1
            tag_apo = 2
            #print("muone CD trovato")
            
        mu_time.append((sec[i]*1e9+nano[i], sec[i], nano[i], tag_apo))
        count_mu += 1
            
for i in range(len(energy_WP)):
    if energy_WP[i] > 400:
        cd_apo = energy[abs((sec_WP[i]-sec)*1e9+(nano_WP[i]-nano)) < 500]
        
        if np.all(cd_apo <= 30000):
            muWP_tag[i] = 3
            #print("muone WP trovato")
            mu_time.append((sec_WP[i]*1e9+nano_WP[i], sec_WP[i], nano_WP[i], 3))
            count_mu += 1
            count_WP += 1
        
        
mu_time = sorted(mu_time, key=lambda x: x[0]) #serve a ordinare mu_time rispetto al primo elemento della tripletta 
mu_times = np.array([t[0] for t in mu_time])
mu_taggato = np.array([a for t,s,n,a in mu_time])

#for i in range(20):
#    print(mu_times[i])
#    print(mu_taggato[i])

print("ora mi metto a calcolare i deltaT_muon")
for i in range(len(energy)):
    idx = np.searchsorted(mu_times, sec[i]*1e9+nano[i]) #questa funzione mi da l'indice dell tempo del muone immediatamente dopo al mio evento
    if idx == 0:
    	deltaT_muon.append(-1)
    else:
    	deltaT_muon.append((sec[i]-mu_time[idx-1][1]) + (nano[i]-mu_time[idx-1][2])*1e-9) 
    	
print("ora mi metto a calcolare il dead_time")
last_end = 0
dead_time = 0    	
for t in mu_times:
    start = t*1e-3
    end = t*1e-3 + veto_window
    # conto solo la parte che non si sovrappone con il precedente
    dead_time += max(0, end - max(start, last_end))
    last_end = max(last_end, end)  # aggiorna fine ultimo intervallo    	
    	
 
#######################################################################################################################################################   	
     


tags = np.zeros(len(sec))

for i in range (1,len(sec)):
    apo = 0
    
    cut_mask_delayed = mu_tag[i] == 0 and deltaT_muon[i]*1e6>veto_window and (x[i]**2+y[i]**2+z[i]**2)**0.5<FV and ((x[i]**2+y[i]**2)**0.5>rho_cut or abs(z[i])<z_cut) and ((sec[i]-sec[0])*1e6 + (nano[i]-nano[0])*1e-3)>2000 and (triggerType[i][0] == "Multiplicity" or (len(triggerType[i])>1 and triggerType[i][1] == "Multiplicity"))
    
    if delayed_min < energy[i] < delayed_max and cut_mask_delayed:
        #print("candidato delayed")
        #print(energy[i])
        #print(sec[i])
        deltat = (sec[i]-sec[i-1])*1e6 + (nano[i]-nano[i-1])*1e-3 #us
        j = 1
        while (deltat <= t_min and j<=i):
            j += 1
            deltat = (sec[i]-sec[i-j])*1e6 + (nano[i]-nano[i-j])*1e-3  
        while (t_min < deltat < t_max and j <= i):
            ds = np.sqrt((x[i]-x[i-j])**2+(y[i]-y[i-j])**2+(z[i]-z[i-j])**2)

            cut_mask_prompt = mu_tag[i-j] == 0 and deltaT_muon[i-j]*1e6>veto_window and (x[i-j]**2+y[i-j]**2+z[i-j]**2)**0.5<FV and ((x[i-j]**2+y[i-j]**2)**0.5>rho_cut or abs(z[i-j])<z_cut) and ((sec[i-j]-sec[0])*1e6 + (nano[i-j]-nano[0])*1e-3)>2000 and (triggerType[i-j][0] == "Multiplicity" or (len(triggerType[i-j])>1 and triggerType[i-j][1] == "Multiplicity"))

            if prompt_min < energy[i-j] < prompt_max and tags[i-j] == 0 and ds < distance and cut_mask_prompt:
                #print("candidato prompt")
                
                # multiplicity cut -----------------------------------------------------------------------------------------------------
                k = 1
                while(apo == 0 and k<=i-j and ((sec[i-j]-sec[i-j-k])*1e6 + (nano[i-j]-nano[i-j-k])*1e-3) <2000):
                    if prompt_min < energy[i-j-k] < prompt_max:
                        apo = 1
                    else:
                        k += 1
                k = 1
                while(apo == 0 and i-j+k<len(sec) and (((sec[i-j+k]-sec[i-j])*1e6 + (nano[i-j+k]-nano[i-j])*1e-3)<2000 or ((sec[i-j+k]-sec[i])*1e6 + (nano[i-j+k]-nano[i])*1e-3)<2000)):
                    if prompt_min < energy[i-j+k] < prompt_max and i-j+k != i:
                        apo = 1
                    else:
                        k += 1
                #while(apo == 0 and i+k<len(sec) and ((sec[i+k]-sec[i])*1e6 + (nano[i+k]-nano[i])*1e-3) <2000):
                # i----------------------------------------------------------------------------------------------------------------------

                if apo == 1:
                    #Print("Ho fallito sul multiplicity!")
                    energyD_veto.append(true_energy[i])
                    energyP_veto.append(true_energy[i-j])
                    xD_veto.append(x[i])
                    yD_veto.append(y[i])
                    zD_veto.append(z[i])
                    xP_veto.append(x[i-j])
                    yP_veto.append(y[i-j])
                    zP_veto.append(z[i-j])
                    time_veto.append(deltat)
                    timeP_veto.append((nano[i-j]*1e-9) + sec[i-j])
                    timeD_veto.append((nano[i]*1e-9) + sec[i])
                    muD_veto.append(0)
                    FVD_veto.append(0)
                    muP_veto.append(0)
                    FVP_veto.append(0)
                    multy_veto.append(1)
                    break
                else:
                    print("coincidenza trovata!")
                    #print("Entry", i)
                    #print(triggerType[i])
                    #print(triggerType[i][0])
                    #if (len(triggerType[i])>1):
                    #    print(triggerType[i][1])
                    energyD_reco.append(true_energy[i])
                    energyP_reco.append(true_energy[i-j])
                    xD_reco.append(x[i])
                    yD_reco.append(y[i])
                    zD_reco.append(z[i])
                    xP_reco.append(x[i-j])
                    yP_reco.append(y[i-j])
                    zP_reco.append(z[i-j])
                    TT5_D.append(TtT5[i])
                    TT5_P.append(TtT5[i-j])
                    tipo_D.append(tipo[i])
                    tipo_P.append(tipo[i-j])
                    coincidenze.append((energy[i-j], energy[i]))
                    time.append(deltat)
                    timeP.append((nano[i-j]*1e-9) + sec[i-j])
                    timeD.append((nano[i]*1e-9) + sec[i])
                    #print((nano[i]**1e-9) + sec[i])
                    print(sec[i])
                    print(nano[i])
                    #print(sec[i]*1e9+(nano[i]))
                    tags[i] = 2
                    tags[i-j] = 1
                    break
            
            else:
                if (prompt_min < energy[i-j] < prompt_max and (mu_tag[i-j] != 0 or deltaT_muon[i-j]*1e6<veto_window)):
                    energyD_veto.append(true_energy[i])
                    energyP_veto.append(true_energy[i-j])
                    xD_veto.append(x[i])
                    yD_veto.append(y[i])
                    zD_veto.append(z[i])
                    xP_veto.append(x[i-j])
                    yP_veto.append(y[i-j])
                    zP_veto.append(z[i-j])
                    time_veto.append(deltat)
                    timeP_veto.append((nano[i-j]*1e-9) + sec[i-j])
                    timeD_veto.append((nano[i]*1e-9) + sec[i])
                    muD_veto.append(0)
                    FVD_veto.append(0)
                    muP_veto.append(1)
                    FVP_veto.append(0)
                    multy_veto.append(0)
                if prompt_min < energy[i-j] < prompt_max and mu_tag[i-j] == 0 and deltaT_muon[i-j]*1e6>veto_window and ((x[i-j]**2+y[i-j]**2+z[i-j]**2)**0.5>FV or ((x[i-j]**2+y[i-j]**2)**0.5<rho_cut and abs(z[i])>z_cut)):
                    energyD_veto.append(true_energy[i])
                    energyP_veto.append(true_energy[i-j])
                    xD_veto.append(x[i])
                    yD_veto.append(y[i])
                    zD_veto.append(z[i])
                    xP_veto.append(x[i-j])
                    yP_veto.append(y[i-j])
                    zP_veto.append(z[i-j])
                    time_veto.append(deltat)
                    timeP_veto.append((nano[i-j]*1e-9) + sec[i-j])
                    timeD_veto.append((nano[i]*1e-9) + sec[i])
                    muD_veto.append(0)
                    FVD_veto.append(0)
                    muP_veto.append(0)
                    FVP_veto.append(1)
                    multy_veto.append(0)
  
                j += 1
                deltat = (sec[i]-sec[i-j])*1e6 + (nano[i]-nano[i-j])*1e-3
        
    else: 
        if (delayed_min < energy[i] < delayed_max and (mu_tag[i] != 0 or deltaT_muon[i]*1e6<veto_window)):
            energyD_veto.append(true_energy[i])
            energyP_veto.append(0)
            xD_veto.append(x[i])
            yD_veto.append(y[i])
            zD_veto.append(z[i])
            xP_veto.append(0)
            yP_veto.append(0)
            zP_veto.append(0)
            time_veto.append(0)
            timeP_veto.append(0)    #((nano[i-j]**1e-9) + sec[i-j])
            timeD_veto.append((nano[i]*1e-9) + sec[i])
            muD_veto.append(1)
            FVD_veto.append(0)
            muP_veto.append(0)
            FVP_veto.append(0)
            multy_veto.append(0)
        if delayed_min < energy[i] < delayed_max and mu_tag[i] == 0 and deltaT_muon[i]*1e6>veto_window and ((x[i]**2+y[i]**2+z[i]**2)**0.5>FV or ((x[i]**2+y[i]**2)**0.5<rho_cut and abs(z[i])>z_cut)):
            energyD_veto.append(true_energy[i])
            energyP_veto.append(0)
            xD_veto.append(x[i])
            yD_veto.append(y[i])
            zD_veto.append(z[i])
            xP_veto.append(0)
            yP_veto.append(0)
            zP_veto.append(0)
            time_veto.append(0)
            timeP_veto.append(0)    #((nano[i-j]**1e-9) + sec[i-j])
            timeD_veto.append(sec[i]*1e9+(nano[i]))   #((nano[i]**1e-9) + sec[i])
            muD_veto.append(0)
            FVD_veto.append(1)
            muP_veto.append(0)
            FVP_veto.append(0)
            multy_veto.append(0) 
    
    
        
prompt = [p for p, d in coincidenze]  
delayed = [d for p, d in coincidenze]

tfile_out = ROOT.TFile(outfile, "RECREATE")
print("Writing in: ",outfile)
tree_out = ROOT.TTree("Coinc", "Coinc")
tree_time = ROOT.TTree("Time", "Time")
tree_veto = ROOT.TTree("Veto", "Veto")
tree_muon = ROOT.TTree("Muon", "Muon")

e_d = np.zeros(1, dtype=np.float32)  
e_p = np.zeros(1, dtype=np.float32)
OECEnergy_d = np.zeros(1, dtype=np.float32)
OECEnergy_p = np.zeros(1, dtype=np.float32)
xP = np.zeros(1, dtype=np.float32)
yP = np.zeros(1, dtype=np.float32)
zP = np.zeros(1, dtype=np.float32)
xD = np.zeros(1, dtype=np.float32)
yD = np.zeros(1, dtype=np.float32)
zD = np.zeros(1, dtype=np.float32)
TtT5_d = np.zeros(1, dtype=np.float32)
TtT5_p = np.zeros(1, dtype=np.float32)
OECTag_d = ROOT.std.vector('string')()
OECTag_p = ROOT.std.vector('string')()
dt = np.zeros(1, dtype=np.float32)
time_p = np.zeros(1, dtype=np.double)
time_d = np.zeros(1, dtype=np.double)
Total_time = np.zeros(1, dtype=np.float32)


Ep_v = np.zeros(1, dtype=np.float32)
Ed_v = np.zeros(1, dtype=np.float32)
xP_v = np.zeros(1, dtype=np.float32)
yP_v = np.zeros(1, dtype=np.float32)
zP_v = np.zeros(1, dtype=np.float32)
xD_v = np.zeros(1, dtype=np.float32)
yD_v = np.zeros(1, dtype=np.float32)
zD_v = np.zeros(1, dtype=np.float32)
timeP_v = np.zeros(1, dtype=np.double)
timeD_v = np.zeros(1, dtype=np.double)
dt_v = np.zeros(1, dtype=np.float32)
muD_v = np.zeros(1, dtype=np.int32)
FVD_v = np.zeros(1, dtype=np.int32)
muP_v = np.zeros(1, dtype=np.int32)
FVP_v = np.zeros(1, dtype=np.int32)
multy_v = np.zeros(1, dtype=np.int32)

life_t = np.zeros(1, dtype=np.float32)
dead_t = np.zeros(1, dtype=np.float32)
n_mu = np.zeros(1, dtype=np.int32)
CD_mu = np.zeros(1, dtype=np.int32)
WP_mu = np.zeros(1, dtype=np.int32)
CD_WP_mu = np.zeros(1, dtype=np.int32)

tag_mu = np.zeros(1,dtype=np.int32)
tempo_mu = np.zeros(1,dtype=np.double)

tree_out.Branch("E_d", e_d, "E_d/F")
tree_out.Branch("E_p", e_p, "E_p/F")
tree_out.Branch("OECEnergy_d", OECEnergy_d, "OECEnergy_d/F")
tree_out.Branch("OECEnergy_p", OECEnergy_p, "OECEnergy_p/F")
tree_out.Branch("xP", xP, "xP/F")
tree_out.Branch("yP", yP, "yP/F")
tree_out.Branch("zP", zP, "zP/F")
tree_out.Branch("xD", xD, "xD/F")
tree_out.Branch("yD", yD, "yD/F")
tree_out.Branch("zD", zD, "zD/F")
tree_out.Branch("TtT5_d", TtT5_d, "TtT5_d/F")
tree_out.Branch("TtT5_p", TtT5_p, "TtT5_p/F")
tree_out.Branch("OECTag_d", OECTag_d)
tree_out.Branch("OECTag_p", OECTag_p)
tree_out.Branch("dt", dt, "dt/F")
tree_out.Branch("time_p", time_p, "time_p/D")
tree_out.Branch("time_d", time_d, "time_d/D")
tree_out.Branch("Total_time", Total_time, "Total_time/F")

tree_veto.Branch("Ep_v", Ep_v, "Ep_v/F")
tree_veto.Branch("Ed_v", Ed_v, "Ed_v/F")
tree_veto.Branch("dt_v", dt_v, "dt_v/F")
tree_veto.Branch("xP_v", xP_v, "xP_v/F")
tree_veto.Branch("yP_v", yP_v, "yP_v/F")
tree_veto.Branch("zP_v", zP_v, "zP_v/F")
tree_veto.Branch("xD_v", xD_v, "xD_v/F")
tree_veto.Branch("yD_v", yD_v, "yD_v/F")
tree_veto.Branch("zD_v", zD_v, "zD_v/F")
tree_veto.Branch("timeP_v", timeP_v, "timeP_v/D")
tree_veto.Branch("timeD_v", timeD_v, "timeD_v/D")
tree_veto.Branch("FVD_v", FVD_v, "FVD_v/I")
tree_veto.Branch("muD_v", muD_v, "muD_v/I")
tree_veto.Branch("muP_v", muP_v, "muP_v/I")
tree_veto.Branch("FVP_v", FVP_v, "FVP_v/I")
tree_veto.Branch("multy_v", multy_v, "multy_v/I")

tree_time.Branch("life_t", life_t, "life_t/F")
tree_time.Branch("dead_t", dead_t, "dead_t/F")
tree_time.Branch("n_mu", n_mu, "n_mu/I")
tree_time.Branch("CD_mu", CD_mu, "CD_mu/I")
tree_time.Branch("CD_WP_mu", CD_WP_mu, "CD_WP_mu/I")
tree_time.Branch("WP_mu", WP_mu, "WP_mu/I")

tree_muon.Branch("tag_mu", tag_mu , "tag_mu/I")
tree_muon.Branch("tempo_mu", tempo_mu , "tempo_mu/D")

for i in range(len(time)):

    OECTag_d.clear()
    OECTag_p.clear()

    for j in range(len(tipo_D[i])):
        OECTag_d.push_back(tipo_D[i][j])

    for j in range(len(tipo_P[i])):
        OECTag_p.push_back(tipo_P[i][j])
    
    Total_time[0] = ((sec[len(sec)-1]-sec[0]) + (nano[len(sec)-1]-nano[0])*1e-9)/len(time)

    e_d[0] = delayed[i]
    e_p[0] = prompt[i]
    OECEnergy_d[0] = energyD_reco[i]
    OECEnergy_p[0] = energyP_reco[i]
    xP[0] = xP_reco[i]
    yP[0] = yP_reco[i]
    zP[0] = zP_reco[i]
    xD[0] = xD_reco[i]
    yD[0] = yD_reco[i]
    zD[0] = zD_reco[i]
    TtT5_d[0] = TT5_D[i]
    TtT5_p[0] = TT5_P[i]
    time_p[0] = timeP[i]
    print(time_p[0])
    time_d[0] = timeD[i]
    print(time_d[0])
    dt[0] = time[i]

    tree_out.Fill()

'''
for i in range(len(time_veto)):
    
    Ep_v[0] = energyP_veto[i]
    Ed_v[0] = energyD_veto[i]
    xP_v[0] = xP_veto[i]
    yP_v[0] = yP_veto[i]
    zP_v[0] = zP_veto[i]
    xD_v[0] = xD_veto[i]
    yD_v[0] = yD_veto[i]
    zD_v[0] = zD_veto[i]
    timeP_v[0] = timeP_veto[i]
    timeD_v[0] = timeD_veto[i]
    dt_v[0] = time_veto[i]
    muD_v[0] = muD_veto[i]
    FVD_v[0] = FVD_veto[i]
    muP_v[0] = muP_veto[i]
    FVP_v[0] = FVP_veto[i] 
    multy_v[0] = multy_veto[i]

    tree_veto.Fill()
'''

life_t[0] = duration
dead_t[0] = dead_time*1e-6
n_mu[0] = count_mu
CD_mu[0] = count_CD
WP_mu[0] = count_WP
CD_WP_mu[0] = count_CD_WP
tree_time.Fill()

for i in range(len(mu_times)):
    tag_mu[0] = mu_taggato[i]
    tempo_mu[0] = mu_times[i]
    tree_muon.Fill()


tfile_out.Write()
tfile_out.Close()
       
        
        
        
        
        
        
        
        
        
        
        
        

