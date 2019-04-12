import numpy as np
from tool.r_cosine import *
from tool.DSPtools import *
from tool._fixedInt import *
from scipy.special import erfc
import matplotlib.pyplot as plt

# ####
# Ruido
# ####
def noise(SNR_a, N,E):    
    SNR_db = SNR_a
    SNR_veces = float(10**(SNR_db/10))
    No = 1.0/SNR_veces #Potencia de ruido
    sigma = np.sqrt(No*E)
    media = 0
    #noise_vect = np.random.normal(media,sigma,N)
    noise_vect=(sigma*(np.random.randn(N)))+media
    return noise_vect

# ####
# Grafica la curva de BER
# ####
def graf_BER(snr,ber):    
    SNR_db=snr
    theoryBER = np.zeros(len(SNR_db),float) #theoritical error
    for i in range(len(SNR_db)):
        theoryBER[i] = erfc(np.sqrt(0.5*(10.0**(float(SNR_db[i]/10))))) - (1/4)*(erfc(np.sqrt(0.5*(10.0**(SNR_db[i]/10)))))**2
    plt.semilogy(snr,ber, 'x')  #Simulation
    plt.semilogy(SNR_db, theoryBER, '--')
    plt.ylabel('log(BER)')
    plt.xlabel('[Db]')
    plt.title('BER Curves')
    plt.legend(['Simulation', 'Theory'], loc='upper right')
    plt.grid()


# ####
# ShiftRegister
# ####
def ShiftReg(vector):   
    '''Funcion que realiza el corrimiento del vector en una posicion.
    '''
    shift_vector = vector
    for i in range(0,len(shift_vector)-1):
        shift_vector[len(shift_vector)-(i+1)] = shift_vector[len(shift_vector)-(i+2)]
    
    return shift_vector

# ####
# Grafica de respuesta al impulso
# ####
def graf_respImpulso(rc,tipo):               
    titulo = "Respuesta al impulso " + tipo
    plt.figure()
    plt.title(titulo)
    plt.stem(rc,'k-',linewidth=2.0,label=r'$\beta=0.5$')
    plt.legend()
    plt.grid(True)
    plt.xlim(0,len(rc)-1)
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')

# ####
# Grafica de respuesta en frecuencia del filtro
# ####
def graf_respFrecuencia(H,A,F,tipo):      
    titulo = "Respuesta en frecuencia " + tipo
    plt.figure()
    plt.title('Respuesta en frecuencia filtro transmisor')
    plt.semilogx(F, 20*np.log10(H),'r', linewidth=2.0, label=r'$\beta=0.5$')
    plt.axvline(x=(1./Ts)/2.,color='k',linewidth=2.0)
    plt.axvline(x=(1./T)/2.,color='k',linewidth=2.0)
    plt.axhline(y=20*np.log10(0.5),color='k',linewidth=2.0)
    plt.legend(loc=3)
    plt.grid(True)
    plt.xlim(F[1],F[len(F)-1])
    plt.xlabel('Frequencia [Hz]')
    plt.ylabel('Magnitud [dB]')

# ####
# Simbolos transmitidos
# ####
def graf_simbTransmitidos(simI,simQ):   
    plt.figure()
    plt.title('Simbolos transmitidos')
    plt.subplot(2,1,1)
    plt.plot(simI,'o')
    plt.xlim(0,20)
    plt.grid(True)
    plt.subplot(2,1,2)
    plt.plot(simQ,'o')
    plt.xlim(0,20)
    plt.grid(True)

# ####
# Convolucion
# ####
def graf_Convolucion(senial,tipo):      
    titulo = "Senial " + tipo
    plt.figure()
    plt.title(titulo)
    plt.plot(senial,'r-',linewidth=2.0)
    plt.xlim(1000,1250)
    plt.grid(True)
    plt.legend()
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')


#############################
# Datos
#############################
T=1.0/100e6 #Periodo de baudio
T_ch=1.0/100e9
os=4 #OverSampling
beta=0.99 #Roll-Off
Nbauds=6.0 #Cantidad de baudios del filtro Tx y Rx
Nsymb=1000000


#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#

def main():
    #Filtro Root Raised Cosine
    (t,rc)=rrcosine(beta, T, os, Nbauds)
    #Normalizacion
    rc=rc[:24]
    rc=rc/np.linalg.norm(rc) #lialg.norm ----> Raiz de la suma de los cuadrados de los coeficientes
    # (t,canal)=rrcosine(beta,T_ch,os,7.0)
    # #canal=canal[:28]
    # #canal=canal/np.linalg.norm(canal)
    # canal=np.zeros(28)
    # canal[(len(canal)/2)+1]=1
    # tx_canal=np.convolve(canal,rc)

    energia_fir_equiv=0
    for w in range(len(rc)):
        energia_fir_equiv=energia_fir_equiv+(rc[w]**2)
    print "ENERGIA DEL FILTRO",energia_fir_equiv

    ## Generacion de simbolos a transmitir
    symbolsI = 2*(np.random.uniform(-1,1,Nsymb)>0.0)-1
    symbolsQ = 2*(np.random.uniform(-1,1,Nsymb)>0.0)-1
    
    ### Over sampling
    os_simbolsI = np.zeros(int(os*Nsymb)); os_simbolsI[1:len(os_simbolsI):int(os)]=symbolsI
    os_simbolsQ = np.zeros(int(os*Nsymb)); os_simbolsQ[1:len(os_simbolsQ):int(os)]=symbolsQ

    ### Snr
    snr=[100.0,100.0,16.0,14.0,12.0,10.0,8.0,6.0,4.0,2.0]
    
    ber_fase0=[]
    ber_fase1=[]
    ber_fase2=[]
    ber_fase3=[]


    delta=0.0115


    #Cantidad de coeficientes FIR del ecualizador
    Ntap=31

    ##
    recib_berI=[]
    recib_berQ=[]

    #Memoria del filtro FIR Ecualizador
    mem_fir_adap_I= np.zeros(Ntap) 
    mem_fir_adap_Q= np.zeros(Ntap) 

    #Coeficientes del filtro FIR del ecualizador adaptativo
    coef_fir_adap_I=np.zeros(Ntap) 
    coef_fir_adap_Q=np.zeros(Ntap) 

    #Pos 15 = 1 inicializacion
    coef_fir_adap_I[(Ntap-1)/2]=1 
    coef_fir_adap_Q[(Ntap-1)/2]=1 
    contador=0
    for j in range (0,len(snr)):
        del recib_berI[:] 
        del recib_berQ[:]
        print "CONTADOR", contador
        contador=0
        ### Convolucion transmisor, senial transmitida
        tx_I = np.convolve(rc,os_simbolsI,'same')
        tx_Q = np.convolve(rc,os_simbolsQ,'same')
        
        # canal_I=np.convolve(tx_I,canal*np.sqrt(np.sqrt(os)/np.sqrt(2)),'same')
        # canal_Q=np.convolve(tx_Q,canal*np.sqrt(np.sqrt(os)/np.sqrt(2)),'same')
        
        # ### NOISE
        noise_vector_I=noise(snr[j],len(os_simbolsI),energia_fir_equiv) #genero senial de ruido
        noise_vector_Q=noise(snr[j],len(os_simbolsQ),energia_fir_equiv)
        
        tx_I_noise=tx_I+noise_vector_I
        tx_Q_noise=tx_Q+noise_vector_Q
        
        ##Convolucion receptor, senial recibida
        rx_I=np.convolve(rc*np.sqrt(np.sqrt(os)/np.sqrt(2)),tx_I_noise,'same')
        rx_Q=np.convolve(rc*np.sqrt(np.sqrt(os)/np.sqrt(2)),tx_Q_noise,'same')
        
        print "VAR", np.var(rx_I)
        print "VAR", np.var(rx_Q)
        #raw_input('Press Enter to Continue')
        for i in range (0,len(rx_I)):

            if((i%2!=0)): #nos quedamos con muestras 1 y 3
                mem_fir_adap_I=ShiftReg(mem_fir_adap_I)
                mem_fir_adap_Q=ShiftReg(mem_fir_adap_Q)

                mem_fir_adap_I[0]=rx_I[i]*0.3
                mem_fir_adap_Q[0]=rx_Q[i]*0.3

            #print "RX_I", rx_I[i]
            
            if(((i+1)%os)==0 and i!=0):
                y_I=np.sum(coef_fir_adap_I*mem_fir_adap_I)           #FIR
                y_Q=np.sum(coef_fir_adap_Q*mem_fir_adap_Q)
                
                # print "Y_I", y_I
                # raw_input('Press Enter to Continue')
            #Slicer
                ak_I=(2*(y_I>=0)-1)
                ak_Q=(2*(y_Q>=0)-1)
 #-----------------------------------------------Algoritmo de adaptacion
                error_adap_I=ak_I-y_I
                error_adap_Q=ak_Q-y_Q
                
                for b in range(0,Ntap):
                    coef_fir_adap_I[b]=coef_fir_adap_I[b] + delta*error_adap_I*mem_fir_adap_I[b]
                    coef_fir_adap_Q[b]=coef_fir_adap_Q[b] + delta*error_adap_Q*mem_fir_adap_Q[b]
                
                # print "COEFICIENTES",coef_fir_adap_I
                # raw_input('Press Enter to Continue')
                
                recib_berI.append(ak_I)
                recib_berQ.append(ak_Q)
                contador=contador+1
        # error_fase=0


        
        # for t in range(0,len(recib_berI)):
        #     if(recib_berI[t]!=symbolsI[t]):
        #         error_fase=error_fase+1
        #     if(recib_berQ[t]!=symbolsQ[t]):
        #         error_fase=error_fase+1

        # ber_fase.append(error_fase/(2*float(Nsymb)))
        
        # print "Ber", ber_fase
        # #Downsample
        of0_I=rx_I[0:len(rx_I):int(os)]
        of1_I=rx_I[1:len(rx_I):int(os)]
        of2_I=rx_I[2:len(rx_I):int(os)]
        of3_I=rx_I[3:len(rx_I):int(os)]

        of0_Q=rx_Q[0:len(rx_Q):int(os)]
        of1_Q=rx_Q[1:len(rx_Q):int(os)]
        of2_Q=rx_Q[2:len(rx_Q):int(os)]
        of3_Q=rx_Q[3:len(rx_Q):int(os)]
      
        #Slicer

        for i in range (0, len(of0_I)):
            of0_I[i] =(2*(of0_I[i]>=0)-1)
            of0_Q[i] =(2*(of0_Q[i]>=0)-1)  
            of1_I[i] =(2*(of1_I[i]>=0)-1)  
            of1_Q[i] =(2*(of1_Q[i]>=0)-1)  
            of2_I[i] =(2*(of2_I[i]>=0)-1)  
            of2_Q[i] =(2*(of2_Q[i]>=0)-1)  
            of3_I[i] =(2*(of3_I[i]>=0)-1)  
            of3_Q[i] =(2*(of3_Q[i]>=0)-1)  

        #BER
        error_fase0=0
        error_fase1=0
        error_fase2=0
        error_fase3=0

        for k in range (0,len(of0_I)):
            if(of0_I[k]!=symbolsI[k]):
                error_fase0=error_fase0+1
            if(of0_Q[k]!=symbolsQ[k]):
                error_fase0=error_fase0+1
            
            if(of1_I[k]!=symbolsI[k]):
                error_fase1=error_fase1+1
            if(of1_Q[k]!=symbolsQ[k]):
                error_fase1=error_fase1+1
            
            if(of2_I[k]!=symbolsI[k]):
                error_fase2=error_fase2+1
            if(of2_Q[k]!=symbolsQ[k]):
                error_fase2=error_fase2+1
            
            if(of3_I[k]!=symbolsI[k]):
                error_fase3=error_fase3+1
            if(of3_Q[k]!=symbolsQ[k]):
                error_fase3=error_fase3+1
        
        print "Error fase 0",  error_fase0
        print "Error fase 1",  error_fase1
        print "Error fase 2",  error_fase2
        print "Error fase 3",  error_fase3
        
        ber_fase0.append(error_fase0/(2*float(Nsymb)))
        ber_fase1.append(error_fase1/(2*float(Nsymb)))
        ber_fase2.append(error_fase2/(2*float(Nsymb)))
        ber_fase3.append(error_fase3/(2*float(Nsymb)))

        print "BER0", ber_fase0
        print "BER1", ber_fase1
        print "BER2", ber_fase2
        print "BER3", ber_fase3


    # plt.figure()
    # graf_BER(snr[2:],ber_fase0[2:])
    # plt.figure()
    # graf_BER(snr[2:],ber_fase1[2:])
    # plt.figure()
    # graf_BER(snr[2:],ber_fase2[2:])
    # plt.figure()
    # graf_BER(snr[2:],ber_fase3[2:])

    # #Correlacion
    # plt.figure()
    # fase0_I=np.correlate(symbolsI,of0_I,'same')
    # fase0_Q=np.correlate(symbolsQ,of0_Q,'same')
    # plt.title("FASE 0")
    # plt.plot(fase0_I)
    # plt.plot(fase0_Q)
    
    # plt.figure()
    # fase1_I=np.correlate(symbolsI,of1_I,'same')
    # fase1_Q=np.correlate(symbolsQ,of1_Q,'same')
    # plt.title("FASE 1")
    # plt.plot(fase1_I)
    # plt.plot(fase1_Q)
    

    # plt.figure()
    # fase2_I=np.correlate(symbolsI,of2_I,'same')
    # fase2_Q=np.correlate(symbolsQ,of2_Q,'same')
    # plt.title("FASE 2")
    # plt.plot(fase2_I)
    # plt.plot(fase2_Q)
    

    # plt.figure()
    # fase3_I=np.correlate(symbolsI,of3_I,'same')
    # fase3_Q=np.correlate(symbolsQ,of3_Q,'same')
    # plt.title("FASE 3")
    # plt.plot(fase3_I)
    # plt.plot(fase3_Q)
    
    plt.show(block=False)
    raw_input('Press Enter to Continue')
    plt.close()
main()