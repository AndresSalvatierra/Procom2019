import numpy as np
from tool.r_cosine import *
from tool.DSPtools import *
from tool._fixedInt import *


## Parametros generales
T     = 1.0/100e6     # Periodo de baudio
prbs_len=511
os    = 4           # Over sampling
## Parametros de la respuesta en frecuencia
Nfreqs = 256          # Cantidad de frecuencias
## Parametros del filtro de caida cosenoidal
beta   = 0.5          # Roll-Off
Nbauds = 6.0          # Cantidad de baudios del filtro
## Parametros funcionales
Ts = T/os              # Frecuencia de muestreo
iteraciones=4*512*1025
def PRBS9(seed):
    a = int(seed)
    newbit = (((a >> 9-1) ^ (a >> 5-1)) & 1)
    a = ((a << 1) | newbit) & 0x1ff
    return a


def main():

    ## Inicializacion de simbolos a transmitir
    symbolsI=0x1AA
    symbolsQ=0x1FE

    #Filtro Root Raised Cosine
    (t,rc)=rrcosine(beta, T, os, Nbauds)

    #Normalizacion
    rc=rc/np.linalg.norm(rc) #lialg.norm ----> Raiz de la suma de los cuadrados de los coeficientes
    rc=rc[:24]

    # rc_pot = sum([i**2 for i in rc])
    # error_pot = sum([(i-j)**2 for i,j in zip(rc0_t, rc)])
    # snr = -10*np.log10(error_pot)

    #print snr

    ## Grafica de respuesta al impulso del filtro
    plt.figure()
    plt.title('Respuesta al impulso')
    plt.plot(rc,'k-',linewidth=2.0,label=r'$\beta=0.5$')
    plt.legend()
    plt.grid(True)
    plt.xlim(0,len(rc)-1)
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')

    ## Respuesta en frecuencia
    [H0,A0,F0]=resp_freq(rc, Ts, Nfreqs)

    ## Grafica de respuesta en frecuencia del filtro
    plt.figure()
    plt.title('Respuesta en frecuencia')
    plt.semilogx(F0, 20*np.log10(H0),'r', linewidth=2.0, label=r'$\beta=0.5$')
    plt.axvline(x=(1./Ts)/2.,color='k',linewidth=2.0)
    plt.axvline(x=(1./T)/2.,color='k',linewidth=2.0)
    plt.axhline(y=20*np.log10(0.5),color='k',linewidth=2.0)
    plt.legend(loc=3)
    plt.grid(True)
    plt.xlim(F0[1],F0[len(F0)-1])
    plt.xlabel('Frequencia [Hz]')
    plt.ylabel('Magnitud [dB]')

    shift_tx_I=np.zeros(24)
    shift_tx_Q=np.zeros(24)
    shift_rx_I=np.zeros(24)
    shift_rx_Q=np.zeros(24)

    trans_berI=np.zeros(1024)
    trans_berQ=np.zeros(1024)

    recib_berI=np.zeros(511)
    recib_berQ=np.zeros(511)

    b_I=0  #Salida de convolucion transmisor
    b_Q=0  #Salida de convolucion transmisor
    rx_pfi_I=0 #Salida de convolucion receptor
    rx_pfi_Q=0   #Salida de convolucion receptor

    tx_I=0
    tx_Q=0

    offset=0

    phase_I=np.zeros(4)
    phase_Q=np.zeros(4)

    error_actual=0
    error_minimo=99999
    pos_trans=0
    pos_aux=0
    cont_ber=0

    graficar=0              #Cuando esta en 1 almacena valores para luego graficar

    #Usados para graficar
    senial_transmitir_I=[]
    senial_transmitir_Q=[]
    senial_recibida_I=[]
    senial_recibida_Q=[]

    Ntap=15 #cantidad de coeficientes
    mem_filtro_I= np.zeros(Ntap-1) #Memoria del filtro FIR
    mem_filtro_Q= np.zeros(Ntap-1) #Memoria del filtro FIR
    coef_I=np.zeros(Ntap-1) #Coeficientes del filtro FIR del equalizador adaptativo
    coef_Q=np.zeros(Ntap-1) #Coeficientes del filtro FIR del equalizador adaptativo
    coef_I[(Ntap-1)/2]=1 #Pos 7 = 1 inicializacion
    coef_Q[(Ntap-1)/2]=1 #Pos 7 = 1 inicializacion
    diesmado=2

    for i in range(0,iteraciones):

        if(i==0):
            sym_t_I=0
            sym_t_Q=0

        else:
            sym_t_I=b_I
            sym_t_Q=b_Q

        if(graficar==1):
            senial_transmitir_I.append(sym_t_I) #Usado para graficar
            senial_transmitir_Q.append(sym_t_Q) #Usado para graficar

        for j in range (0,len(shift_tx_I)-1):
            shift_tx_I[len(shift_tx_I)-(j+1)]=shift_tx_I[len(shift_tx_I)-(j+2)]
            shift_tx_Q[len(shift_tx_Q)-(j+1)]=shift_tx_Q[len(shift_tx_Q)-(j+2)]

        if(i%os==0):
            symbolsI=PRBS9(symbolsI)    #Generacion de simbolos
            symbolsQ=PRBS9(symbolsQ)
            if((symbolsI>>8)& 0x001==0):
                shift_tx_I[0]=1
                # if(graficar==1):
                #     filep = open('prbsi.txt','a')
                #     filep.write("0\n")
                #     filep.close()
            else:
                shift_tx_I[0]=-1
                # if(graficar==1):
                #     filep = open('prbsi.txt','a')
                #     filep.write("1\n")
                #     filep.close()

            if((symbolsQ>>8)& 0x001==0):
                shift_tx_Q[0]=1
                # if(graficar==1):
                #     filep = open('prbsq.txt','a')
                #     filep.write("0\n")
                #     filep.close()
            else:
                shift_tx_Q[0]=-1
                # if(graficar==1):
                #     filep = open('prbsq.txt','a')
                #     filep.write("1\n")
                #     filep.close()

            for m in range (0,len(trans_berI)-1):
                trans_berI[len(trans_berI)-(m+1)]=trans_berI[len(trans_berI)-(m+2)]
                trans_berQ[len(trans_berQ)-(m+1)]=trans_berQ[len(trans_berQ)-(m+2)]

            trans_berI[0]=shift_tx_I[0]  #Arreglo de bits transmitidos
            trans_berQ[0]=shift_tx_Q[0]

        else:                           #Oversampling
            shift_tx_I[0]=0
            shift_tx_Q[0]=0

        ## Convolucion transmisor, senial transmitida
        tx_I=np.sum(shift_tx_I*rc)
        tx_Q=np.sum(shift_tx_Q*rc)

        b_I=tx_I
        # if(graficar==1):
        #     filep = open('txi.txt','a')
        #     filep.write(str(b_I))
        #     filep.write("\n")
        #     filep.close()
        b_Q=tx_Q



        for k in range (0,len(shift_rx_I)-1):
            shift_rx_I[len(shift_rx_I)-(k+1)]=shift_rx_I[len(shift_rx_I)-(k+2)]
            shift_rx_Q[len(shift_rx_Q)-(k+1)]=shift_rx_Q[len(shift_rx_Q)-(k+2)]

        shift_rx_I[0]=sym_t_I
        shift_rx_Q[0]=sym_t_Q

        ##Convolucion receptor, senial recibida
        rx_I=np.sum(shift_rx_I*rc)
        rx_Q=np.sum(shift_rx_Q*rc)

        rx_pfi_I=rx_I
        #if(graficar==1):
            # filep = open('rxi.txt','a')
            # filep.write(str(rx_pfi_I))
            # filep.write("\n")
            # filep.close()
        rx_pfi_Q=rx_Q
        #if(graficar==1):
            # filep = open('rxq.txt','a')
            # filep.write(str(rx_pfi_Q))
            # filep.write("\n")
            # filep.close()

        if((i%diesmado)==0):

            mem_filtro_I[0]=rx_pfi_I
            mem_filtro_Q[0]=rx_pfi_Q
            y_I=0
            y_Q=0
            for r in range(0,Ntap-1):
                y_I=y_I+coef_I[r]*mem_filtro_I[r]           #FIR
                y_Q=y_Q+coef_Q[r]*mem_filtro_Q[r]           #FIR

            for s in range(len(mem_filtro_I)-1,0,-1):
                mem_filtro_I[s]=mem_filtro_I[s-1]   #Actualizacion de la memoria FIR
                mem_filtro_Q[s]=mem_filtro_Q[s-1]   #Actualizacion de la memoria FIR


                if(graficar==1):
                    senial_recibida_I.append(y_I)  #Usada para graficar
                    senial_recibida_Q.append(y_Q)  #Usada para graficar

            for t in range(0,len(phase_I)-1):
                phase_I[t]=phase_I[t+1]
                phase_Q[t]=phase_Q[t+1]

            #Slicer
            if(y_I>=0):
                y_I=1
            else:
                y_I=-1

            if(y_Q>=0):
                y_Q=1
            else:
                y_Q=-1

            phase_I[3]=y_I
            phase_Q[3]=y_Q

        ##Ber
        if((i%os)==0):
            for m in range (0,len(recib_berI)-1):
                recib_berI[len(recib_berI)-(m+1)]=recib_berI[len(recib_berI)-(m+2)]
                recib_berQ[len(recib_berQ)-(m+1)]=recib_berQ[len(recib_berQ)-(m+2)]

            recib_berI[0]=phase_I[offset]  #Arreglo de bits recibidos
            recib_berQ[0]=phase_Q[offset]  #Arreglo de bits recibidos

            if(cont_ber<511):
                cont_ber=cont_ber+1
            else:
                print error_actual
                print pos_trans
                if(error_actual<error_minimo):
                    error_minimo=error_actual
                    error_actual=0
                    pos_aux=pos_trans

                else:
                    error_actual=0
                cont_ber=0

                if(pos_trans!=1023 and graficar==0):
                    pos_trans=pos_trans+1
                else:
                    pos_trans=pos_aux
                    graficar=1
            if(recib_berI[0]!=trans_berI[pos_trans]):
                error_actual=error_actual+1
            if(recib_berQ[0]!=trans_berQ[pos_trans]):
                error_actual=error_actual+1

    ##Simbolos transmitidos
    print error_minimo
    print pos_aux
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(shift_tx_I,'o')
    plt.xlim(0,20)
    plt.grid(True)
    plt.subplot(2,1,2)
    plt.plot(shift_tx_Q,'o')
    plt.xlim(0,20)
    plt.grid(True)

    ### Convolucion transmisor, senial transmitida
    plt.figure()
    plt.title("Senial transmitida")
    plt.plot(senial_transmitir_I,'r-',linewidth=2.0)
    plt.xlim(1000,1250)
    plt.grid(True)
    plt.legend()
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')

    ##Grafico diagrama ojo transmisor
    eyediagram(senial_transmitir_I[100:len(senial_transmitir_I)-100],os,0,Nbauds)
    eyediagram(senial_transmitir_Q[100:len(senial_transmitir_Q)-100],os,0,Nbauds)

    ##Convolucion receptor, senial recibida
    plt.figure()
    plt.title("Senial recibida")
    plt.plot(senial_recibida_I,'r-',linewidth=2.0)
    plt.xlim(1000,1250)
    plt.grid(True)
    plt.legend()
    plt.xlabel('Muestras')
    plt.ylabel('Magnitud')

    ##Grafico diagrama ojo receptor
    eyediagram(senial_recibida_I[100:len(senial_recibida_I)-100],os,0,Nbauds)
    eyediagram(senial_recibida_Q[100:len(senial_recibida_Q)-100],os,0,Nbauds)

    ##Constelaciones
    plt.figure()
    plt.plot(senial_recibida_I[(offset):len(senial_recibida_I)-((offset)):int(os)],
         senial_recibida_Q[(offset):len(senial_recibida_Q)-((offset)):int(os)],
             'r.',linewidth=2.0)
    plt.axis('equal')
    plt.grid(True)
    plt.xlabel('Real')
    plt.ylabel('Imag')


    plt.show(block=False)
    raw_input('Press Enter to Continue')
    plt.close()

main()
