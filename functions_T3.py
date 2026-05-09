import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numba import njit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def response_spectrum_NEC_2024(Z,soil_type,r_region):
    
    if Z>=0.1 and Z<0.2:
        seismic_zone="I"
    
    if Z>=0.2 and Z<0.3:
        seismic_zone="II"
    
    if Z>=0.3 and Z<0.4:
        seismic_zone="III"

    if Z>=0.4 and Z<0.5:
        seismic_zone="IV"

    if Z>=0.5 and Z<0.65:
        seismic_zone="V"


    if r_region == "Provincias de la Costa":
        r=1.2
    else:
        r=1
    
    Fa=pd.DataFrame({"A":{"I":0.90,"II":0.90,"III":0.90,"IV":0.90,"V":0.90},
                    "B":{"I":1.00,"II":1.00,"III":1.00,"IV":1.00,"V":1.00},
                    "C":{"I":1.40,"II":1.30,"III":1.23,"IV":1.19,"V":1.13},
                    "D":{"I":1.60,"II":1.40,"III":1.25,"IV":1.14,"V":1.00},
                    "E":{"I":1.80,"II":1.4,"III":1.1,"IV":0.90,"V":0.62}})   
   
    Fd=pd.DataFrame({"A":{"I":0.90,"II":0.90,"III":0.90,"IV":0.90,"V":0.90},
                    "B":{"I":1.00,"II":1.00,"III":1.00,"IV":1.00,"V":1.00},
                    "C":{"I":1.36,"II":1.28,"III":1.15,"IV":1.08,"V":1.00},
                    "D":{"I":1.62,"II":1.45,"III":1.28,"IV":1.15,"V":1.00},
                    "E":{"I":2.10,"II":1.75,"III":1.65,"IV":1.52,"V":1.36}})   
   
    Fs=pd.DataFrame({"A":{"I":0.75,"II":0.75,"III":0.75,"IV":0.75,"V":0.75},
                    "B":{"I":0.75,"II":0.75,"III":0.75,"IV":0.75,"V":0.75},
                    "C":{"I":0.85,"II":0.94,"III":1.06,"IV":1.17,"V":1.28},
                    "D":{"I":1.02,"II":1.06,"III":1.19,"IV":1.32,"V":1.44},
                    "E":{"I":1.50,"II":1.6,"III":1.80,"IV":1.94,"V":2.09}})
    
    Fas=Fa[soil_type][seismic_zone]
    Fds=Fd[soil_type][seismic_zone]
    Fss=Fs[soil_type][seismic_zone]
    
    To=0.1*Fss*Fds/Fas
    Tc=0.45*Fss*Fds/Fas
    Tl=2.4*Fds

    T=np.arange(0, 10, 0.005)
    n=len(T)
    Sa=np.zeros(n)
    
    for i in range(n):
        Ti=T[i]
        if Ti<To:
            Sa[i]=Z*Fas*(1+1.4*(Ti/To))
        if Ti>=To and Ti<Tc:
            Sa[i]=2.4*Z*Fas
        if Ti>=Tc and Ti<Tl:
            Sa[i]=2.4*Z*Fas*(Tc/Ti)**r
        if Ti>=Tl:
            Sa[i]=2.4*Z*Fas*((Tc/Tl)**r)*(Tl/Ti)**2
    
    return Sa,T
    


@njit
def CalculoDinamicoBNewmarkNuevo(time, a, Ti, dT, Tf, zi):

    m = 1.0
    g = 9.81
    a = a * g  # convertir a m/s²

    dt = time[1] - time[0]
    ndatos = len(a)
    nT = int((Tf - Ti)/dT)
    Sa = np.zeros(nT)
    Sv = np.zeros(nT)
    Sd = np.zeros(nT)
    Tperiod = np.zeros(nT)

    for j in range(nT):
        T = Ti + j*dT
        w = 2*np.pi / T
        k = w**2 * m
        c = 2*zi*w*m

        # reiniciar en cada periodo
        xn = np.zeros(ndatos)
        xvn = np.zeros(ndatos)
        xan = np.zeros(ndatos)
        for i in range(ndatos-1):
            if i == 0:
                xan[i] = (-a[i] - c*xvn[i] - k*xn[i]) / m
            xn[i+1] = xn[i] + dt*xvn[i] + 0.5*dt**2*xan[i]
            xan[i+1] = (-a[i+1]-k*xn[i+1]-c*(xvn[i]+0.5*dt*xan[i]))/(m + 0.5*c*dt)
            xvn[i+1] = xvn[i] + dt*(0.5*xan[i] + 0.5*xan[i+1])
        at = a + xan
        Sd[j] = np.max(np.abs(xn))
        Sv[j] = np.max(np.abs(xvn))
        Sa[j] = np.max(np.abs(at)) / g   # devolver en g
        Tperiod[j] = T

    return Sa, Sv, Sd, Tperiod

@njit
def SPEC_ROTD100_50_BNewmarkNuevo(time, SGNS, SGEW, Ti, dT, Tf, zi, dTheta):

    dim0 = len(time)
    dim1 = int((Tf - Ti)/dT)
    
    thetas = np.arange(0.0, 180.0, dTheta)
    
    thetas = np.append(thetas, 90.0)
    
    dim2 = len(thetas)

    SGCOMB = np.zeros((dim0, dim2))
    SaComb = np.zeros((dim1, dim2))
    SvComb = np.zeros((dim1, dim2))
    SdComb = np.zeros((dim1, dim2))

    for j, theta in enumerate(thetas):
        theta_r = theta * np.pi / 180.0
        
        SGCOMB[:, j] = SGNS*np.cos(theta_r) + SGEW*np.sin(theta_r)

        Sa,Sv,Sd,Tperiod = CalculoDinamicoBNewmarkNuevo(time,SGCOMB[:, j],Ti, dT, Tf, zi)

        SaComb[:, j] = Sa
        SvComb[:, j] = Sv
        SdComb[:, j] = Sd

    RODT100 = np.zeros(dim1)
    RODT50  = np.zeros(dim1)

    for i in range(dim1):
        max_val = SaComb[i, 0]
        suma = 0.0
        for j in range(dim2):
            val = SaComb[i, j]
            if val > max_val:
                max_val = val
            suma += val
        RODT100[i] = max_val
        RODT50[i]  = suma / dim2

    return SaComb, SvComb, SdComb, RODT100, RODT50,Tperiod

def minimum_square_error(Spec_obje,Spec_record):
    n=len(Spec_obje)
    sum1=0
    sum2=0
    for i in range(n):
        sum1+=Spec_obje[i]*Spec_record[i]
        sum2+=Spec_record[i]**2
    k=sum1/sum2

    return k

def escalamiento_amplitud(sismo):
    if not sismo:
        st.error("⚠️ No existen Eventos ingresados")                                            
        return
    
    col_config, col_plot = st.columns([1,3])
    with col_config:
        Z = st.number_input(f"Ingrese Acelereación de Roca Z (Tr 475 años)", min_value=0.05, max_value=0.80, value=0.30, step=0.01, format="%.2f")
        r_region=st.selectbox("Selecciona la región:",["Provincias de la Costa","Provincias de la Sierra y Oriente"])
        soil_type=st.selectbox("Selecciona el Tipo de Suelo:",["A","B","C","D","E"])
        T_lower= st.number_input(f"Ingrese Periodo en el limite inferior de acuerdo a norma ASCE:", min_value=0.01, max_value=10.0, value=0.5, step=0.01, format="%.2f")
        T_upper= st.number_input(f"Ingrese Periodo en el limite superior de acuerdo a norma ASCE:", min_value=0.01, max_value=10.0, value=1.0, step=0.01, format="%.2f")

        st.divider()

        sf= st.number_input(f"Ingrese Factor de Escala a Tr 2475y:", min_value=0.10, max_value=2.0, value=1.5, step=0.01, format="%.3f")
        Sa,Ti=response_spectrum_NEC_2024(Z,soil_type,r_region)
        MCE=Sa*sf
        mask=(Ti>=T_lower)&(Ti<=T_upper)
        MCE90=0.9*MCE[mask]
        
        st.divider()

        if "resultados_espectro_ROTD100" not in st.session_state:
            st.session_state["resultados_espectro_ROTD100"] = {}

        if "resultados_espectro_ROTD100_escalado" not in st.session_state:
            st.session_state["resultados_espectro_ROTD100_escalado"] = {}
        
        if "resultados_espectro_ROTD100_escalado2" not in st.session_state:
            st.session_state["resultados_espectro_ROTD100_escalado2"] = {}

        if "calcular_escalamiento" not in st.session_state:
            st.session_state["calcular_escalamiento"] = {}

        if "resultados_espectro_0" not in st.session_state:
            st.session_state["resultados_espectro_0"] = {}

        if "resultados_espectro_90" not in st.session_state:
            st.session_state["resultados_espectro_90"] = {}

        if "resultados_espectro_0_escalado" not in st.session_state:
            st.session_state["resultados_espectro_0_escalado"] = {}

        if "resultados_espectro_90_escalado" not in st.session_state:
            st.session_state["resultados_espectro_90_escalado"] = {}

        if "tiempo_arias" in st.session_state:
            for nombre in list(st.session_state["resultados_espectro_ROTD100"].keys()):
                if nombre not in st.session_state["tiempo_arias"]:
                    st.session_state["resultados_espectro_ROTD100"].pop(nombre, None)
                    st.session_state["resultados_espectro_ROTD100_escalado"].pop(nombre, None)

        if "evento_corregido" in st.session_state:

            if st.button("Calcular escalamiento de registros",key=f"Boton_calcular"):
                st.divider()
                st.session_state["calcular_escalamiento"] = True
                for nombre_evento, datos in st.session_state["evento_corregido"].items():
                    time=st.session_state["evento_corregido"][nombre_evento]["Tiempo (s)"]
                    an= st.session_state["evento_corregido"][nombre_evento]["a_N"]/9.8067
                    ae= st.session_state["evento_corregido"][nombre_evento]["a_E"]/9.8067
                    
                    SaComb, _, _, RODT100, _,Tperiod=SPEC_ROTD100_50_BNewmarkNuevo(time,an,ae,0.005,0.005,T_upper+1,0.05,2.0)                     
                    mask1=(Tperiod>=T_lower)&(Tperiod<=T_upper)
                    
                    k=minimum_square_error( MCE[mask],RODT100[mask1])

                    st.session_state["resultados_espectro_0"][nombre_evento] =  SaComb[:, 0]
                    st.session_state["resultados_espectro_90"][nombre_evento] = SaComb[:, -1]
                    st.session_state["resultados_espectro_ROTD100"][nombre_evento] = (Tperiod, RODT100)
                    st.session_state["resultados_espectro_ROTD100_escalado"][nombre_evento] = (Tperiod, k*RODT100,k)
                    st.session_state["resultados_espectro_ROTD100_escalado2"][nombre_evento] = (Tperiod, k*RODT100,k)

                espectros_ROTD100=[]
                
                for nombre_evento, (Tperiod, Spec,k) in st.session_state["resultados_espectro_ROTD100_escalado"].items():
                    espectros_ROTD100.append(Spec)

                espectros_ROTD100=np.array(espectros_ROTD100)
                average_ROTD100 = np.mean(espectros_ROTD100, axis=0)
                mask1=(Tperiod>=T_lower)&(Tperiod<=T_upper)
                st.session_state["mask1"] = mask1
                
                k2=np.min(average_ROTD100[mask1]/MCE90)
                if k2>=1:
                    k2=1
                else:
                    k2=1/k2
                espectros_ROTD100_2=[]
                for nombre_evento2, (Tperiod2, Spec,k) in st.session_state["resultados_espectro_ROTD100_escalado"].items():
                    st.session_state["resultados_espectro_ROTD100_escalado2"][nombre_evento2] = (Tperiod2,k2*Spec,k*k2)
                    espectros_ROTD100_2.append(k2*Spec)

                for nombre, data in st.session_state["resultados_espectro_0"].items():
                    st.session_state["resultados_espectro_0_escalado"][nombre]=data*k*k2

                for nombre, data in st.session_state["resultados_espectro_90"].items():
                    st.session_state["resultados_espectro_90_escalado"][nombre]=data*k*k2

                espectros_ROTD100_2=np.array(espectros_ROTD100_2)
                average_ROTD100 = np.mean(espectros_ROTD100_2, axis=0)
                st.session_state["average_ROTD100"] = average_ROTD100

            if st.session_state["calcular_escalamiento"]:    
                for nombre_evento, datos in st.session_state["evento_corregido"].items():
                    col1, col2 = st.columns([2.5, 2.5])
                    
                    col1.write(f"📊 {nombre_evento}")
                
                    col2.checkbox("Graficar", key=f"plot_{nombre_evento}")
                
    with col_plot:
       
        st.markdown("# *Seccion de Escalado de registros Mediante [MSE]*")
        st.markdown("<h1 style='color:blue;'>Espectro objetivo (Borrador NEC24), unicamente demostrativo</h1>",unsafe_allow_html=True)

        st.markdown("### *La presente sección tiene como objetivo la verificación del numeral 16.2.3.2 (*Amplitude Scaling*). (ASCE 7-22).*")
        fig, ax = plt.subplots(figsize=(8, 4), dpi=500)
        ax.plot(Ti, MCE, label=f"Espectro objetivo", lw=0.7)
        
        ax_inset = inset_axes(ax, width="30%", height="25%", loc="upper center")
        ax_inset.tick_params(axis='both', labelsize=3, direction='in')
        ax_inset.plot(Ti[mask], MCE[mask],lw=0.7)
        
        ax_inset.set_ylabel("Aceleración (g)", fontsize=4)

        if "evento_corregido" in st.session_state:
            for nombre_evento, (Tperiod, Spec,k) in st.session_state["resultados_espectro_ROTD100_escalado2"].items():
                if st.session_state.get(f"plot_{nombre_evento}", False):
                    ax.plot(Tperiod,Spec,label=f"ROTD100-escalado k={k:.3f}, {nombre_evento}", lw=0.5)
                    ax_inset.plot(Tperiod[ st.session_state["mask1"]],Spec[ st.session_state["mask1"]], lw=0.5)

        ax_inset.plot(Ti[mask], MCE90,lw=0.7, label="MCE90",ls='-.')
      
        if st.session_state["calcular_escalamiento"]:
            ax.plot(Tperiod, st.session_state["average_ROTD100"], color="black", lw=1, label="Promedio ROTD100")
            if "average_ROTD100" in st.session_state:
                ax_inset.plot(Tperiod[st.session_state["mask1"]], st.session_state["average_ROTD100"][st.session_state["mask1"]], color="black", lw=1)

        ax.set_title("Sa")
        ax.set_xlabel("Periodo T(sec)")
        ax.set_ylabel("Aceleración (g)")
        ax.legend(loc='upper right', fontsize=5)
        ax.grid(True,lw=0.5)
        ax.set_ylim(0, 2.5)
        ax.set_xlim(0, T_upper*2)
        st.pyplot(fig)

