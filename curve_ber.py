# main.py
import   numpy as np
from     scipy.special import erfc
import   matplotlib.pyplot as plt

# ####
# Grafica la curva de BER
# ####
def graf_BER(snr,ber):    
   SNR_db=snr
   theoryBER = np.zeros(len(SNR_db),float) #theoritical error
   for i in range(len(SNR_db)):
      #theoryBER[i] = erfc(np.sqrt(10**(SNR_db[i]/10)))
      theoryBER[i] = erfc(np.sqrt(0.5*(10.0**(SNR_db[i]/10)))) - (1/4)*(erfc(np.sqrt(0.5*(10.0**(SNR_db[i]/10)))))**2
      theoryBER[i] = theoryBER[i]/np.sqrt(4)
   plt.semilogy(snr,ber, 'x')  #Simulation
   plt.semilogy(SNR_db, theoryBER, '--')
   plt.ylim(10e-08)
   plt.xlim(0,13)
   plt.ylabel('log(BER)')
   plt.xlabel('[Db]')
   plt.title('BER Curves')
   plt.legend(['Simulation', 'Theory'], loc='upper right')
   plt.grid()
   #plt.show()

 
SNR_db = np.array(np.arange(0, 18, 2),float)
snr=[16.,14.,12.,10.,8.,6.,4.,2.]

#VALORES DE BER CON RUIDO
#ber_fase_0 = [0.0, 6e-06, 0.000176, 0.002962, 0.006498, 0.024546, 0.059262, 0.10712] #Con 8,6,4,2 ponderados
ber_fase_0 = [0.0, 2e-06, 9.4e-05, 0.001268, 0.00761, 0.026376, 0.060972, 0.10709] #Con 16,12,14,10,8,6,4,2 ponderados
# ber_fase_0=[4e-06, 0.0001955, 0.001852, 0.008551, 0.024828, 0.054447, 0.0972895, 0.148002]
# ber_fase_1=[0.250302,0.250575,0.2508235,0.2510045,0.2530255,0.2614175,0.278696,0.301809]
# ber_fase_2=[ 0.0, 0.0, 3.25e-05, 0.000751, 0.00596, 0.0229845, 0.0564575, 0.103797]
# ber_fase_3=[5e-06,0.000185,0.0018345,0.008479,0.0249475,0.054707,0.0976795,0.147922]


#VALORES DE BER CON CANAL Y RUIDO
# ber_fase_0=[4.5e-06, 0.000159, 0.0015375, 0.0075685, 0.0215345, 0.0443135, 0.076653, 0.1182115]
# ber_fase_1=[0.0,0.0,1.4e-05,0.0003445,0.0028515,0.012574,0.0356785,0.074395]
# ber_fase_2=[3.5e-06,0.0001545,0.00157,0.0075245,0.021423,0.0444345,0.076656,0.117781]

# for i in range (len(ber_fase_0)):
#    ber_fase_0[i] = ber_fase_0[i]*np.sqrt(2)
#    ber_fase_1[i] = ber_fase_1[i]*np.sqrt(2)
#    ber_fase_2[i] = ber_fase_2[i]*np.sqrt(2)
#    ber_fase_3[i] = ber_fase_3[i]*np.sqrt(2)


# #VALORES CON ECUALIZADOR
# equ=[3.5e-06, 2.5e-06, 3.3e-05, 0.000958, 0.421837, 0.4997705,0.501389,0.499458]

# theoryBER = np.zeros(len(SNR_db),float) #theoritical error
# theoryBER_vieja=np.zeros(len(SNR_db),float) #theoritical error

# for i in range(len(SNR_db)):
#     theoryBER[i] = erfc(np.sqrt(0.5*(10.0**(SNR_db[i]/10)))) - (1/4)*(erfc(np.sqrt(0.5*(10.0**(SNR_db[i]/10)))))**2
   
#    # theoryBER_vieja[i] = erfc(np.sqrt(10.0**(SNR_db[i]/10)))

# plt.semilogy(snr,equ, 'x')  #Simulation
# plt.semilogy(SNR_db, theoryBER, '--')
# #plt.semilogy(SNR_db,theoryBER_vieja,'--')
# plt.ylim(10e-08)
# plt.xlim(0,13)
# plt.ylabel('log(BER)')
# plt.xlabel('[Db]')
# plt.title('BER Curves')
# plt.legend(['ber_fase_0', 'Theory'], loc='upper right')
# plt.grid()
# plt.show()

##Grafica curvas Ber-teorica y simulada
plt.figure()
graf_BER(snr,ber_fase_0)
# plt.figure()
# graf_BER(snr,ber_fase_1)
# plt.figure()
# graf_BER(snr,ber_fase_2)
# plt.figure()
# graf_BER(snr,ber_fase_3)


plt.show(block=False)
raw_input('Press Enter to Continue')
plt.close()

