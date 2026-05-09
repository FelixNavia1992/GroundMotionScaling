import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import zipfile
import io

def guardar_txt(nombre, tiempo, aceleracion):
    lineas = []
    
    for t, a in zip(tiempo, aceleracion):
        lineas.append(f"{t:.5f} {a:.6e}\n")
    
    contenido = "".join(lineas)
    return contenido.encode("utf-8")

def crear_zip(registros, modo_tiempo):
    buffer = io.BytesIO()
    
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        
        for nombre, data in registros.items():
            K=st.session_state["resultados_espectro_ROTD100_escalado2"][nombre][2]
           
            t = data["Tiempo (s)"]
            accN = data["a_N"]*K
            accE = data["a_E"]*K
            
            if modo_tiempo == "Tiempo recortado (Arias 5-95)%":
                info = st.session_state.tiempo_arias[nombre]["Info"]
                i0 = info["t5_pos"]
                i1 = info["t95_pos"]
                t = t[i0:i1+1]
                accN = accN[i0:i1+1]
                accE = accE[i0:i1+1]
                t = t - t[0]
           
            if st.session_state.get(f"check_Y_{nombre}", False):
                accN, accE = accE*K, accN*K
            else:
                accN, accE = accN*K, accE*K
            
            archivo_txtN= guardar_txt(nombre, t, accN)
            archivo_txtE= guardar_txt(nombre, t, accE)
            z.writestr(f"{nombre}_Y.txt", archivo_txtN)
            z.writestr(f"{nombre}_X.txt", archivo_txtE)
    
    buffer.seek(0)
    return buffer


def direccionalidad(sismo):

    if not sismo:
        st.error("⚠️ No existen Eventos ingresados")                                             
        return
    col_conf, col_plt = st.columns([2.5, 7.5]) 
    
    if st.session_state["resultados_espectro_90_escalado"] is None:
        return

    with col_conf:
        for nombre in st.session_state["resultados_espectro_90_escalado"]:
            if f"check_X_{nombre}" not in st.session_state:
                st.session_state[f"check_X_{nombre}"] = True
            if f"check_Y_{nombre}" not in st.session_state:
                st.session_state[f"check_Y_{nombre}"] = False
            
            def activar_x(nombre=nombre):
                st.session_state[f"check_X_{nombre}"] = True
                st.session_state[f"check_Y_{nombre}"] = False

            def activar_y(nombre=nombre):
                st.session_state[f"check_X_{nombre}"] = False
                st.session_state[f"check_Y_{nombre}"] = True

            col1, col2, col3 = st.columns([2.5, 2.5 ,2.5]) 
            col1.write(f"📊 {nombre}")
            col2.checkbox("Sin Rotar",key=f"check_X_{nombre}",on_change=activar_x)
            col3.checkbox("Rotar",key=f"check_Y_{nombre}",on_change=activar_y) 

        n = len(st.session_state["resultados_espectro_90_escalado"])*1.0
       
        mean_component_spectraXY=None
        mean_component_spectraX=None
        mean_component_spectraY=None

        nombre=0

        for nombre, data in st.session_state["resultados_espectro_90_escalado"].items():
            if mean_component_spectraXY is None:
                mean_component_spectraXY=np.zeros_like(data)
            mean_component_spectraXY+=(st.session_state["resultados_espectro_0_escalado"][nombre]+st.session_state["resultados_espectro_90_escalado"][nombre])

            if mean_component_spectraX is None:
                mean_component_spectraX=np.zeros_like(data)
           
            if mean_component_spectraY is None:
                mean_component_spectraY=np.zeros_like(data)

            if st.session_state[f"check_Y_{nombre}"]:
                mean_component_spectraX+=st.session_state["resultados_espectro_0_escalado"][nombre]
                mean_component_spectraY+=st.session_state["resultados_espectro_90_escalado"][nombre]
            elif st.session_state[f"check_X_{nombre}"]:
                mean_component_spectraX+=st.session_state["resultados_espectro_90_escalado"][nombre]
                mean_component_spectraY+=st.session_state["resultados_espectro_0_escalado"][nombre]
            #count+=1
        if nombre !=0:
            ultimo_nombre=nombre
            if mean_component_spectraXY is not None and n > 0:
                mean_component_spectraXY = mean_component_spectraXY /(2*n)
                mean_component_spectraX = mean_component_spectraX /(n)
                mean_component_spectraY = mean_component_spectraY /(n)

            mean90XY=0.90*mean_component_spectraXY
            mean110XY=1.10*mean_component_spectraXY

    st.write(n)

    with col_plt:

        st.markdown("# *Seccion de selección de dirección de registros*")
        st.markdown("### *La presente sección tiene como objetivo la verificación del numeral 16.2.4 (*Application of Ground Motions to the Structural Model*). (ASCE 7-22)."
        "No aplica para Fuentes cercanas.*")

        if nombre !=0:
            Tperiod, _, _ = st.session_state["resultados_espectro_ROTD100_escalado2"][ultimo_nombre]
            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)
        
            ax.plot(Tperiod, mean_component_spectraXY, label=f"Espectro promedio XY", lw=0.3)
            ax.plot(Tperiod,  mean90XY, label=f"Espectro promedio XY -10%", lw=0.7, ls='--', color='black')
            ax.plot(Tperiod,  mean110XY, label=f"Espectro promedio XY +10%", lw=0.7, ls='--', color='grey')
            ax.plot(Tperiod, mean_component_spectraX, label=f"Espectro promedio X", lw=0.7)
            ax.plot(Tperiod, mean_component_spectraY, label=f"Espectro promedio Y", lw=0.7)

            mask=st.session_state["mask1"]

            ax_inset = inset_axes(ax, width="40%", height="30%", loc="upper center")
            ax_inset.tick_params(axis='both', labelsize=3, direction='in')
            ax_inset.plot(Tperiod[mask], mean_component_spectraXY[mask],lw=0.3)
            ax_inset.plot(Tperiod[mask], mean90XY[mask],lw=0.7,color='black',ls='--')
            ax_inset.plot(Tperiod[mask], mean110XY[mask],lw=0.7,color='grey',ls='--')
            ax_inset.plot(Tperiod[mask], mean_component_spectraX[mask], lw=0.7)
            ax_inset.plot(Tperiod[mask], mean_component_spectraY[mask], lw=0.7)

            ax.set_title("Sa")
            ax.set_xlabel("Periodo T(sec)")
            ax.set_ylabel("Aceleración (g)")
            ax.legend(loc='upper right', fontsize=5)
            ax.grid(True,lw=0.5)
            ax.set_ylim(0, 2.5)
            st.pyplot(fig)


    if f"check_X_{nombre}" not in st.session_state:
       return 
    
    #modo_tiempo = st.radio("⏱️ Selección de duración del registro", ["Tiempo completo", "Tiempo recortado (Arias 5-95)%"], horizontal=False)
    modo_tiempo = st.radio("⏱️ Selección de duración del registro", ["Tiempo completo"], horizontal=False)
    zip_buffer = crear_zip(st.session_state["evento_corregido"], modo_tiempo)

    st.download_button(
    label="📦 Descargar registros",
    data=zip_buffer,
    file_name="registros_escalados.zip",
    mime="application/zip")    

    st.markdown("### *Los registros descargados se encuentran escalados en unidad de (m/s^2), el usuario es responsable del uso de la información obtenida.*")