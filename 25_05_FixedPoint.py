import numpy as np
import matplotlib.pyplot as plt
from tool.r_cosine import *
from tool.DSPtools import *
from tool._fixedInt import *
from scipy.special import erfc
from commpy.filters import rrcosfilter
from commpy.filters import rcosfilter

#----------------------------FUNCIONES---------------------------#

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

# ####
# Prbs
# ####
def PRBS9(seed):        
    a = int(seed)
    newbit = (((a >> 9-1) ^ (a >> 5-1)) & 1)
    a = ((a << 1) | newbit) & 0x1ff
    return a

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
# Vector de cuantizacion (punto fijo)
# ####
def VectorCuantizacion(NB,NBF,v,mode):          
        N = arrayFixedInt(NB,NBF,v,'S',mode,'saturate')
        v_cuant = []
        for i in range (len(v)):
                v_cuant.append(N[i].fValue)
        
        return v_cuant

# ####
# Grafica de respuesta al impulso
# ####
def graf_respImpulso(rc,tipo):               
    titulo = "Respuesta al impulso " + tipo
    plt.figure()
    plt.title(titulo)
    #plt.stem(rc,'k-',linewidth=2.0)
    plt.plot(rc,'k-',linewidth=2.0)
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
    #plt.xlim(1000,1250)
    plt.grid(True)
    plt.legend()
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')

# ####
# Grafica la curva de BER
# ####
def graf_BER(snr,ber):    
    SNR_db=snr
    theoryBER = np.zeros(len(SNR_db),float) #theoritical error
    for i in range(len(SNR_db)):
        theoryBER[i] = erfc(np.sqrt(10**(SNR_db[i]/10)))
    plt.semilogy(snr,ber, 'x')  #Simulation
    plt.semilogy(SNR_db, theoryBER, '--')
    plt.ylabel('log(BER)')
    plt.xlabel('[Db]')
    plt.title('BER Curves')
    plt.legend(['Simulation', 'Theory'], loc='upper right')
    plt.grid()
    plt.show()

#---------------------------------------------Parametros generales------------------------------------------------#
T  = 1.0/100e6  # Periodo de baudio
T_ch = 1.0/100e9  # Periodo de baudio
prbs_len=511
os    = 4           # Over sampling
os_ch = 4

## Parametros de la respuesta en frecuencia
Nfreqs = 256          # Cantidad de frecuencias

## Parametros del filtro de caida cosenoidal
beta   = 0.99        # Roll-Off
Nbauds = 6.0          # Cantidad de baudios del filtro
Nbauds_ch = 7.0        # Cantidad de baudios del filtro channel

## Parametros funcionales
Ts = T/os              # Frecuencia de muestreo
Ts_ch = T_ch/os_ch              # Frecuencia de muestreo
iteraciones=4*511*1030
E=0.0                   # Energia del filtro resultante (convolucion) entre tx y canal

def main():

    ## Inicializacion de simbolos a transmitir
    symbolsI=0x1AA
    symbolsQ=0x1FE

    ## Filtro Tx y Rx
    #Filtro Root Raised Cosine
    rc=rrcosfilter(int(Nbauds*os),beta,T,1./Ts)[1]
    #Normalizacion
    rc=rc[:24]
    #rc=rc/np.linalg.norm(rc) #lialg.norm ----> Raiz de la suma de los cuadrados de los coeficientes
    rc=rc/np.sqrt(float(os))
    #print rc
    rc_fp= VectorCuantizacion(9,8,rc,'trunc')
    
    #print "\n", rc_fp

    
    # ### FIR para Canal
    # #Filtro Root Raised Cosine 
    rch=rcosfilter(int(Nbauds_ch*os),beta,T_ch,1./Ts_ch)[1]
    # # #Normalizacion fir del canal
    #rch=rch[:28]
    rch=rch/np.linalg.norm(rch)
    #rch=rch/np.linalg.norm(rch) #lialg.norm ----> Raiz de la suma de los cuadrados de los coeficientes
    rch=rch/np.sqrt(float(os))
    # rch= np.zeros(28)
    # rch[13]=1
    # rch=rch/np.sqrt(float(1))
    rch_fp= VectorCuantizacion(9,8,rch,'trunc')

    canal=np.convolve(rc,rch)
    canal_fp=np.convolve(rc_fp,rch_fp)  ## filtro equivalente, resultado de la convolucion de filtro Tx y FIR canal
    canal_fp=VectorCuantizacion(9,8,canal,'trunc')
    
    ###---------------------------------------------------
    ####Calculos de energia de filtros
    #Tx
    energia_tx=0
    for p in range(len(rc)):
        energia_tx=energia_tx+(rc_fp[p]**2)
    print "Energia Tx = ", energia_tx

    # ##FIR canal
    # energia_rch=0
    # for q in range(len(rch)):
    #     energia_rch=energia_rch+(rch[q]**2)
    # print "Energia rch = ",energia_rch

   
    ##Filtro equivalente de Tx y FIR , canal equivalente
    energia_fir_equiv=0
    for w in range(len(canal)):
        energia_fir_equiv=energia_fir_equiv+(canal[w]**2)

    print "Energia canal equivalente: ", energia_fir_equiv

    ###---------------------------------------------------
    #Grafica de respuesta al impulso del filtro
    # graf_respImpulso(rc_fp,"Filtro Transmisor Tx")
    # graf_respImpulso(rch_fp,"Filtro Canal")
    # graf_respImpulso(canal,"Rc Rch")
    
    plt.figure()
    plt.title("Filtro RRC")
    plt.plot(rc,'k-',linewidth=2.0)
    plt.plot(rc_fp,'b-', linewidth =2.0)
    plt.legend()
    plt.grid(True)
    plt.xlim(0,len(rc)-1)
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')


    plt.figure()
    plt.title("Convolucion TxRx")
    plt.plot(canal,'k-',linewidth=2.0)
    plt.plot(canal_fp,'b-', linewidth =2.0)
    plt.legend()
    plt.grid(True)
    plt.xlim(0,len(canal)-1)
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')

    # plt.figure()
    # plt.title("Canal")
    # plt.plot(canal,'k-',linewidth=2.0)
    # plt.plot(canal_fp,'b-', linewidth =2.0)
    # plt.legend()
    # plt.grid(True)
    # plt.xlim(0,len(canal)-1)
    # plt.xlabel('Muestras')
    # plt.ylabel('Magnitud')

    print "Convolucion RC PuntoFlotante \n", canal , "\n \n"
    print "Convolucion RC PunfoFijo \n" , canal_fp, "\n \n"

    #Respuesta en frecuencia
    [H0,A0,F0]=resp_freq(rc, Ts, Nfreqs)
    [H1,A1,F1]=resp_freq(rc_fp, Ts, Nfreqs)

    #Grafica de respuesta en frecuencia del filtro
    #graf_respFrecuencia(H0,A0,F0,"Filtro Transmisor")

    plt.figure()
    plt.title('Respuesta en frecuencia convolucion filtro transmisor y receptor')
    plt.semilogx(F0, 20*np.log10(H0),'r', linewidth=2.0, label=r'$\beta=0.5$')
    plt.semilogx(F1, 20*np.log10(H1),'b', linewidth=2.0, label=r'$\beta=0.5$')  
    plt.axvline(x=(1./Ts)/2.,color='k',linewidth=2.0)
    plt.axvline(x=(1./T)/2.,color='k',linewidth=2.0)
    plt.axhline(y=20*np.log10(0.5),color='k',linewidth=2.0)
    plt.legend(loc=3)
    plt.grid(True)
    plt.xlim(F0[1],F0[len(F0)-1])
    plt.xlabel('Frequencia [Hz]')
    plt.ylabel('Magnitud [dB]')


    #Grafica de respuesta al impulso del filtro
    # graf_respImpulso(rcoseno,"convolucion filtro transmisor y receptor")

    #Respuesta en frecuencia
    # [H1,A1,F1]=resp_freq(rcoseno, Ts, Nfreqs)

    #Grafica de respuesta en frecuencia del filtro
    # graf_respFrecuencia(H1,A1,F1,"convolucion filtro transmisor y receptor")

    #Shift register del transmisor y receptor
    shift_tx_I=np.zeros(24) 
    shift_tx_Q=np.zeros(24)
    shift_rx_I=np.zeros(24)
    shift_rx_Q=np.zeros(24)

    #Bit transmitidos para comparar en la BER
    trans_berI=np.zeros(1024)
    trans_berQ=np.zeros(1024)

    #Bit recibios para comprar en la BER
    recib_berI=np.zeros(511)
    recib_berQ=np.zeros(511)

    #Salida del filtro transmisor
    out_tx_I=0
    out_tx_Q=0
    out_tx_I_fp= DeFixedInt(9,7,'S','trunc') #Punto Fijo
    out_tx_Q_fp= DeFixedInt(9,7,'S','trunc')

    #Salida de del canal mas el ruido
    out_ch_I=0
    out_ch_Q=0

    out_ch_I_fp= DeFixedInt(9,8) #Punto Fijo
    out_ch_Q_fp= DeFixedInt(9,8) 

    #Ruido
    noise_vector_I_fp = DeFixedInt(9,7,'S','trunc')
    noise_vector_Q_fp = DeFixedInt(9,7,'S','trunc')

    #Canal + Ruido
    sym_t_I_fp=DeFixedInt(10,7,'S','trunc')
    sym_t_Q_fp=DeFixedInt(10,7,'S','trunc')

    #Salida del filtro receptor
    rx_I = 0
    rx_Q = 0 
    rx_I_fp = DeFixedInt(12,9,'S','trunc') #Punto Fijo
    rx_Q_fp = DeFixedInt(12,9,'S','trunc')

    #Cantidad de coeficientes FIR del ecualizador
    Ntap=31

    #Memoria del filtro FIR Ecualizador
    mem_fir_adap_I= np.zeros(Ntap) 
    mem_fir_adap_Q= np.zeros(Ntap)

    #Coeficientes del filtro FIR del ecualizador adaptativo
    coef_fir_adap_I= np.zeros(Ntap) 
    coef_fir_adap_Q= np.zeros(Ntap)
    #Pos 15 = 1 inicializacion
    coef_fir_adap_I[(Ntap-1)/2]=1 
    coef_fir_adap_Q[(Ntap-1)/2]=1

    coef_fir_adap_I_fp=DeFixedInt(20,17,'S','trunc')
    coef_fir_adap_Q_fp=DeFixedInt(20,17,'S','trunc')

    #Error realimentacion Equalizador
    error_adap_I = 0
    error_adap_Q = 0
    error_adap_I_fp = DeFixedInt(12,9,'S','trunc')
    error_adap_Q_fp = DeFixedInt(12,9,'S','trunc')
    
    #Valor obtenido del Slicer (detector)
    ak_I=0
    ak_Q=0
    #Salida del FIR del Equalizador
    y_I=0
    y_Q=0
    y_I_fp=DeFixedInt(13,10,'S','trunc')
    y_Q_fp=DeFixedInt(13,10,'S','trunc')
    
    promedio_error = []
    promedio_error2 = []

    offset=1
    
    #Variables utilizadas en la sincronizacion
    error_actual=0
    error_minimo=99999
    pos_trans=0
    pos_aux=0
    cont_ber=0

    #Variables utilizadas post sincronizacion
    error_final=0
    cant_muestras=0
    
    #Habilitacion Prbs, Ber
    value=0

    delta=0.0115
    # delta=0.005

    diezmado=2

    end_sync=0 #Bandera indica fin etapa de sincronizacion (1)

    graficar=0 #Cuando esta en 1 almacena valores para luego graficar

    # #Memoria del canal
    mem_canal_I=np.zeros(len(rch))
    mem_canal_Q=np.zeros(len(rch))


    #Usados para graficar
    senial_transmitir_I=[]
    senial_transmitir_Q=[]

    senial_recibida_I=[]
    senial_recibida_Q=[]

    media_oRx=[]

    senial_recibida_diezmada_I=[]
    senial_recibida_diezmada_Q=[]

    o_ecu=[]

    coeficientes=[]

    error_tiempo=[] 

    input_equ=[]
    var_tx=[]
    var_rx=[]
    var_ch=[]

    #vect_pond=[0.6407,0.6369,0.6356,0.6331,0.6310,0.6275,0.6238,0.6210,0.6097] #ponderacion segun medias
    vect_pond=[0.93,0.65,0.62,0.53,0.47,0.42,0.38,0.32,0.28] #ponderacion a ojo segun entrada ecu
    # vect_pond=[0.38,0.28,0.24]
    #CURVA BER SIMULADA
    snr=[100.,16.,14.,12.,10.,8.,6.,4.,2.]
    # snr=[6.,4.,2.]
    ber=[]
    # snr_iteraciones=[40000,40000,40000]
    snr_iteraciones=[2097153,400000+1000000,100000+1000000,100000+1000000,100000+1000000,100000+1000000,100000+1000000,100000+1000000,100000+1000000]
    #snr_iteraciones=[100000,100000,100000,100000,100000,100000,100000,100000,100000]
    for t in range(len(snr_iteraciones)):
        error_final=0
        noise_vector_I=noise(snr[t],snr_iteraciones[t],energia_fir_equiv/1.2)#/(os**2)) #genero senial de ruido
        noise_vector_Q=noise(snr[t],snr_iteraciones[t],energia_fir_equiv/1.2)#/(os**2))
  
        print "snr_iter: ",snr_iteraciones[t]
        print "t: ",t
        
        print "media: ", np.mean(media_oRx)
        del media_oRx[:]

        for i in range(0,snr_iteraciones[t]):
            value=i%os
            if(i==0):
                sym_t_I=0
                sym_t_Q=0
            else:
                # print "noise: ", noise_vector_I[i]
                
                noise_vector_I_fp.value=noise_vector_I[i]
                noise_vector_Q_fp.value=noise_vector_Q[i]
                sym_t_I=out_ch_I_fp.fValue +noise_vector_I_fp.fValue
                sym_t_Q=out_ch_Q_fp.fValue +noise_vector_Q_fp.fValue

                sym_t_I_fp.value=sym_t_I
                sym_t_Q_fp.value=sym_t_Q
                # print "noise fp: ", noise_vector_I_fp.fValue
                # raw_input('Press Enter to Continue')

            shift_tx_I = ShiftReg(shift_tx_I)
            shift_tx_Q = ShiftReg(shift_tx_Q)

            if(value==0):   #es multiplo de 4?
                symbolsI=PRBS9(symbolsI)    #Generacion de simbolos
                symbolsQ=PRBS9(symbolsQ)
                if((symbolsI>>8)& 0x001==0):
                    shift_tx_I[0]=1
                else:
                    shift_tx_I[0]=-1

                if((symbolsQ>>8)& 0x001==0):
                    shift_tx_Q[0]=1
                else:
                    shift_tx_Q[0]=-1

                trans_berI = ShiftReg(trans_berI)
                trans_berQ = ShiftReg(trans_berQ)

                trans_berI[0]=shift_tx_I[0]  #Agrego nuevo valor transmitido a arreglo bits transmitidos
                trans_berQ[0]=shift_tx_Q[0]

            else:                           #Oversampling
                shift_tx_I[0]=0
                shift_tx_Q[0]=0


            ## Convolucion transmisor, senial transmitida
            out_tx_I=np.sum(rc_fp*shift_tx_I)
            out_tx_Q=np.sum(rc_fp*shift_tx_Q)
            
            # if(i%4):
            #     var_tx.append(out_tx_I)

            out_tx_I_fp.value=out_tx_I
            out_tx_Q_fp.value=out_tx_Q



            # ###
            # Canal
            # ###
            mem_canal_I = ShiftReg(mem_canal_I) #corrimiento de memoria del canal
            mem_canal_Q = ShiftReg(mem_canal_Q)
            
            mem_canal_I[0]=out_tx_I_fp.fValue #agregamos nuevo valor de salida de Tx a la memoria
            mem_canal_Q[0]=out_tx_Q_fp.fValue

            out_ch_I=np.sum(mem_canal_I*(rch_fp)) #salida Tx pasa por canal
            out_ch_Q=np.sum(mem_canal_Q*(rch_fp))

            out_ch_I_fp.value = out_ch_I
            out_ch_Q_fp.value = out_ch_Q

            # if(i%4):
            #     var_ch.append(out_ch_I)
            
            shift_rx_I = ShiftReg(shift_rx_I)
            shift_rx_Q = ShiftReg(shift_rx_Q)

            shift_rx_I[0]=sym_t_I_fp.fValue            
            shift_rx_Q[0]=sym_t_Q_fp.fValue


            ##Convolucion receptor, senial recibida
            rx_I=np.sum(rc_fp*shift_rx_I)*vect_pond[t]
            rx_Q=np.sum(rc_fp*shift_rx_Q)*vect_pond[t]
            
            
            rx_I_fp.value=rx_I
            rx_Q_fp.value=rx_Q        

            # if(i%4 and t==0):
            #     input_equ.append(abs(rx_I))
            # rx_I_fp = rx_I *vect_pond[t]
            # rx_Q_fp = rx_Q *vect_pond[t]

            # rx_I=rx_I*pond
            # rx_Q=rx_Q*pond

            # if(max_rx[t] < abs(rx_I)):
            #     max_rx[t] = abs(rx_I)

            # print "Salida Rx: ",rx_I
            senial_recibida_I.append(rx_I_fp.fValue)
            media_oRx.append(abs(rx_I_fp.fValue))

            # print "Salida sin ruido:" , rx_I
            # print "Salida con ruido:" , sym_t_I

            
            ####
            ###FIR ECUALIZADOR
            #####
            if((i%2!=0) and (i!=0)): #nos quedamos con muestras 1 y 3

                mem_fir_adap_I=ShiftReg(mem_fir_adap_I)
                mem_fir_adap_Q=ShiftReg(mem_fir_adap_Q)

                mem_fir_adap_I[0]=rx_I_fp.fValue
                mem_fir_adap_Q[0]=rx_Q_fp.fValue

            if((value==0) and i!=0):

                
                #for r in range(0,Ntap-1):
                y_I=np.sum(coef_fir_adap_I*mem_fir_adap_I)           #FIR
                y_Q=np.sum(coef_fir_adap_Q*mem_fir_adap_Q)           #FIR

                y_I_fp.value=y_I
                y_Q_fp.value=y_Q

                o_ecu.append(y_I_fp.fValue)

                #Slicer
                ak_I=(2*(y_I_fp.fValue>=0)-1)
                ak_Q=(2*(y_Q_fp.fValue>=0)-1)

                #------------------------------------------Algoritmo de adaptacion
                error_adap_I=ak_I-y_I_fp.fValue
                error_adap_Q=ak_Q-y_Q_fp.fValue

                # promedio_error.append(abs(error_adap_I))

                error_adap_I_fp.value=error_adap_I
                error_adap_Q_fp.value=error_adap_Q

                # promedio_error2.append(abs(error_adap_I_fp.fValue))
                                
                for b in range(0,Ntap):
                    coef_fir_adap_I[b]=coef_fir_adap_I[b] + delta*error_adap_I_fp.fValue*mem_fir_adap_I[b]
                    coef_fir_adap_Q[b]=coef_fir_adap_Q[b] + delta*error_adap_Q_fp.fValue*mem_fir_adap_Q[b]
                
                # print("COEF I SIN FIX")
                # print ("sin fix",coef_fir_adap_I)
                # print("NOSE",coef_fir_adap_I_fp)
                # for b in range(0,Ntap):
                #     coef_fir_adap_I_fp.value=coef_fir_adap_I[b]
                #     coef_fir_adap_I[b]=coef_fir_adap_I_fp.fValue

                # print("FIX",coef_fir_adap_I)

                for b in range(0,Ntap):
                    coef_fir_adap_I_fp.value=coef_fir_adap_I[b]
                    coef_fir_adap_I[b]=coef_fir_adap_I_fp.fValue
                    coef_fir_adap_Q_fp.value=coef_fir_adap_Q[b]
                    coef_fir_adap_Q[b]=coef_fir_adap_Q_fp.fValue
                 
                if(i%8==0): 
                    coeficientes.append(coef_fir_adap_I.copy())    
                # # raw_input('Press Enter to Continue')
                # print("COEF I CON FIX")
                # for b in range(0,Ntap):
                #     print (coef_fir_adap_I[b])
            #Ber
            if(i>100000):
                if((value==0) and i!=0):
                    recib_berI = ShiftReg(recib_berI)
                    recib_berQ = ShiftReg(recib_berQ)

                    recib_berI[0]=ak_I #Bits recibidos
                    recib_berQ[0]=ak_Q #Bits recibidos

                    if(end_sync==0): ##Etapa de sincronizacion
                        if(cont_ber<511):
                            cont_ber=cont_ber+1
                        else: #cont_ber==511
                            cont_ber=0
                            if(error_actual<error_minimo):
                                print "Error minimo actual: ",error_actual
                                print "Pos error minimo actual: ",pos_trans
                                error_minimo=error_actual
                                error_actual=0
                                pos_aux=pos_trans
                            else:
                                error_actual=0

                            if(pos_trans!=1023):
                                pos_trans=pos_trans+1
                            else:   #finalizo etapa de sincronizacion
                                pos_trans=pos_aux
                                end_sync=1
                                # print "end_sync",end_sync
                                #graficar=1
                        if(recib_berI[0]!=trans_berI[pos_trans]):
                            error_actual=error_actual+1
                        if(recib_berQ[0]!=trans_berQ[pos_trans]):
                            error_actual=error_actual+1

                    else: #end_sync==1 pos sincronizacion
                        cant_muestras=cant_muestras+1
                        # print "end_sync"
                        # if(end_sync==1):
                            # print "Bit recibido en BER: ", ak_I
                            # raw_input('Press Enter to Continue')
                            # plt.close()

                        if(recib_berI[0]!=trans_berI[pos_trans]):
                            error_final=error_final+1
                        if(recib_berQ[0]!=trans_berQ[pos_trans]):
                            error_final=error_final+1
                        
                        if(i==(snr_iteraciones[t]-4)):
                            print "snr",snr[t]
                            print "error_final", error_final
                            ber.append(float(error_final)/((((snr_iteraciones[t]-100000)/os)*2)))
                            print ber

    
    # print "VAR_TX", np.var(var_tx)
    # print "VAR_RX", np.var(var_rx)
    # print "VAR_CH", np.var(var_ch)
    
        # print "SNR = ", snr[t], " Maximo = ", max_rx[t]
        # print "SNR = %d - Maximo = %d" %(snr[t],max_rx[t])
    # suma=0
    # for w in range (len(input_equ)):
    #     suma=suma + input_equ[w]
    print "BER", ber
    # plt.figure()
    # plt.title("Entrada Ecualizador")
    # plt.plot(input_equ,".")
    # print "Cantidad de muestras para la media ", len(input_equ) 
    # print "MEDIA:", suma/len(input_equ)
    # print "Vector de Maximos", max_rx
    #print "BER: ", ber

    #print error_minimo
    #print pos_aux
    
    ##Grafico de simbolos transmitidos
    #graf_simbTransmitidos(shift_tx_I,shift_tx_Q)
    
    ### Convolucion transmisor, senial transmitida
    #graf_Convolucion(senial_transmitir_I, "transmitida")
    
    ##Grafico diagrama ojo transmisor
    #eyediagram(senial_transmitir_I[100:len(senial_transmitir_I)-100],os,0,Nbauds)
    #eyediagram(senial_transmitir_Q[100:len(senial_transmitir_Q)-100],os,0,Nbauds)
    
    ##Convolucion receptor, senial recibida
    #graf_Convolucion(var_rx,"recibida")

    ##Grafica curvas Ber-teorica y simulada
    #graf_BER(snr,ber)

    print "media: ", np.mean(media_oRx)
    del media_oRx[:]

    plt.figure()
    plt.title("Entrada al ecualizador")
    plt.plot(senial_recibida_I)

    plt.figure()
    plt.title("Salida Ecualizador")
    plt.plot(o_ecu)

    ##Grafico diagrama ojo receptor
    #eyediagram(senial_recibida_I[100:len(senial_recibida_I)-100],os,0,Nbauds)
    #eyediagram(senial_recibida_Q[100:len(senial_recibida_Q)-100],os,0,Nbauds)

    ##Constelaciones
    # plt.figure()
    # plt.title("Constelaciones")
    # plt.plot(senial_recibida_diezmada_I[(offset):len(senial_recibida_diezmada_I)-((offset)):int(os)],
    #      senial_recibida_diezmada_Q[(offset):len(senial_recibida_diezmada_Q)-((offset)):int(os)],
    #          'r.',linewidth=2.0)
    # plt.axis('equal')
    # plt.grid(True)
    # plt.xlabel('Real')
    # plt.ylabel('Imag')
    
    # plt.figure()
    # plt.title("Error en la adaptacion")
    # plt.plot(error_tiempo)
    
    plt.figure()
    plt.title("Coeficientes en el tiempo")
    plt.plot(coeficientes)
    
    # print len(coeficientes)
    plt.figure()
    plt.title("Ultimos coeficientes del filtro FIR")
    plt.stem(coeficientes[len(coeficientes)-1])
    
    #equa=np.convolve(rcoseno,np.transpose(coeficientes[len(coeficientes)-1]))
    
    # ## Grafica de respuesta al impulso del filtro
    # graf_respImpulso(equa,"a la salida del equalizador")
    
    # ## Respuesta en frecuencia
    # [H2,A2,F2]=resp_freq(equa, Ts, Nfreqs)
    # ## Grafica de respuesta en frecuencia del filtro
    # graf_respFrecuencia(H2,A2,F2,"a la salida del ecualizador")

    plt.show(block=False)
    raw_input('Press Enter to Continue')
    plt.close()

main()