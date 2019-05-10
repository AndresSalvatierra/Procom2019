import numpy as np
import matplotlib.pyplot as plt
from commpy.filters import rrcosfilter
from commpy.filters import rcosfilter


# ####
# Ruido
# ####
def noise(SNR_a, N,E):    
    SNR_db = SNR_a
    SNR_veces = float(10**(SNR_db/10))
    No = 1.0/SNR_veces #Potencia de ruido
    sigma = np.sqrt(No*E)
    media = 0
    # noise_vect = np.random.normal(media,sigma,N)
    noise_vect=(sigma*(np.random.randn(N)))+media
    return noise_vect


T  = 1.0/100e6  # Periodo de baudio
T_ch = 1.0/100e9  # Periodo de baudio
prbs_len=511
os    = 4           # Over sampling
os_ch = 4.

## Parametros de la respuesta en frecuencia
Nfreqs = 256          # Cantidad de frecuencias

## Parametros del filtro de caida cosenoidal
beta   = 0.99        # Roll-Off
Nbauds = 6.0          # Cantidad de baudios del filtro
Nbauds_ch = 7.0        # Cantidad de baudios del filtro channel

## Parametros funcionales
Ts = T/os              # Frecuencia de muestreo
Ts_ch = T_ch/os_ch              # Frecuencia de muestreo



# h_rrc = rrcosfilter(int(Nbauds_ch*os),beta,T_ch,1./Ts_ch)[1]  #Root Raised Cosine
# h_rc = rcosfilter(int(Nbauds_ch*os),beta,T,1./Ts)[1]    #Raised Cosine

# titulo = "Respuesta al impulso "
# plt.title(titulo)
# plt.plot(h_rrc,'b-',linewidth=2.0,label=r'$\beta=0.5$')
# plt.plot(h_rc,'r-',linewidth=2.0,label=r'$\beta=0.5$')
# plt.legend()
# plt.grid(True)
# plt.xlim(0,len(h_rrc)-1)
# plt.xlabel('Muestras')
# plt.ylabel('Magnitud')
# plt.show()


vector = noise(2,4000000,3)
for i in range (len(vector)):
    vector[i] = np.abs(vector[i])
    
print np.mean(vector) 