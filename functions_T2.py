import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def desp_vel_acce_Base(df):                                                                     # Function to calculate base displacements from seismic record
    t = df["Tiempo (s)"]                                                                        # Time array from seismic record
    a = df["Aceleración (g)"] * 9.8067                                                          # Acceleration array converted to m/s²     
    n = len(a)                                                                                  # Number of data points
    v = np.zeros(n); r = np.zeros(n)                                                            # Initialize velocity and displacement arrays
    v[0] = 0; r[0] = 0                                                                          # Initial conditions: zero velocity and displacement
    for i in range(n-1):                                                                        # Loop through each time step to perform numerical integration
        dt=t[i+1]-t[i]                                                                          # Calculate time step                             
        v[i+1]=( a[i+1]-a[i]) * dt/2 + a[i]*dt + v[i]                                           # Update velocity      
        r[i+1]=(a[i+1]-a[i]) * (dt**2) / 6.0 + a[i]*(dt**2)/2 + v[i]*dt + r[i]                  # Update displacement 
    df_resultado = pd.DataFrame({
        "Tiempo (s)": t,
        "Aceleración (m/s²)": a,
        "Velocidad (m/s)": v,
        "Desplazamiento (m)": r})
    return df_resultado                                                                           

def arias_intensity(df):
    t = df["Tiempo (s)"]
    a = df["Aceleración (g)"]*9.8067
    dt=t[1]-t[0]
    n = len(a)
    arias=np.zeros(n)
    for i in range(n-1):
        arias[i+1]=np.pi/(2*9.8067)*(a[i+1]**2+a[i]**2)*dt/2+arias[i]
    df_arias=pd.DataFrame({
        "Tiempo (s)": t,
        "Arias (m/s)": arias
    })
    return df_arias

def correccion_linea_base(t, a, grado):
    coef = np.polyfit(t, a, grado)
    tendencia = np.polyval(coef, t)
    a_corr = a - tendencia
    return a_corr, coef


def tratamiento_datos(sismo):
    
    if "tiempo_arias" not in st.session_state:
        st.session_state.tiempo_arias = {}

    if not sismo:
        st.error("⚠️ No existen Eventos ingresados")                                             # Warning
        return

    col_config, col_plot = st.columns([1,3]) 
    
    with col_config:

        st.subheader("⚙️ Corrección de Linea Base")
        nombre_evento=st.selectbox("Selecciona evento a corregir", list(sismo.keys()))

        t = sismo[nombre_evento]['N']['df']["Tiempo (s)"].values
        a_N = sismo[nombre_evento]['N']['df']["Aceleración (g)"].values * 9.8067
        a_E = sismo[nombre_evento]['E']['df']["Aceleración (g)"].values * 9.8067

        corregir_N = st.checkbox("Norte")
        if corregir_N:
            metodo_N = st.radio("Método N", ["Lineal","Parabólica"])
            grado_N = 1 if metodo_N == "Lineal" else 2
            a_corr_N, coef_N = correccion_linea_base(t, a_N, grado_N)                     #aceleracion entra en m/s^2 (sale como vector)
            st.write("Coef Norte:", coef_N)
        else:
           a_corr_N=a_N
        
        corregir_E = st.checkbox("Este")

        if corregir_E:
            metodo_E = st.radio("Método E", ["Lineal","Parabólica"])
            grado_E = 1 if metodo_E == "Lineal" else 2
            a_corr_E, coef_E = correccion_linea_base(t, a_E, grado_E)
            st.write("Coef Este:", coef_E)

        else:
            a_corr_E=a_E
        
        df_soil_dynamic_N=desp_vel_acce_Base(sismo[nombre_evento]['N']['df'])       
        df_soil_dynamic_E=desp_vel_acce_Base(sismo[nombre_evento]['E']['df'])

        #CALCULO DE INTENSIDAD DE ARIAS
        #CALCULO DE 5% AL 95%
        
        df_N = pd.DataFrame({
        "Tiempo (s)": t,
        "Aceleración (g)": a_corr_N / 9.8067
        })

        df_E = pd.DataFrame({
            "Tiempo (s)": t,
            "Aceleración (g)": a_corr_E / 9.8067
        })    
        
        df_soil_dynamic_N_corr=desp_vel_acce_Base(df_N)       
        df_soil_dynamic_E_corr=desp_vel_acce_Base(df_E)

        arias_N=arias_intensity(df_N)
        arias_E=arias_intensity(df_E)

        arias_N_E=arias_N.copy()
        arias_N_E["Arias (m/s)"]= arias_N["Arias (m/s)"]+arias_E["Arias (m/s)"]
        arias_norm=arias_N_E["Arias (m/s)"]/arias_N_E["Arias (m/s)"].iloc[-1]
        
        pos5=arias_norm[arias_norm >= 0.05].index[0]
        pos95=arias_norm[arias_norm >= 0.95].index[0]

        t5=arias_N_E["Tiempo (s)"].iloc[pos5]
        t95=arias_N_E["Tiempo (s)"].iloc[pos95]


        if st.button("✅ Aplicar registro", key="Boton guardar registro"):
            if "evento_corregido" not in st.session_state:
                st.session_state["evento_corregido"] = {}

            st.session_state["evento_corregido"][nombre_evento] = {
                "Tiempo (s)": t,
                "a_N": a_corr_N,
                "a_E": a_corr_E
            }
            st.success("Corrección guardada")

            st.session_state.tiempo_arias[nombre_evento] = {                                         #GUARDAR EN MEMORIA PARA STREAMLIT
            "Info": {"t5_pos": pos5, "t95_pos": pos95}
            }

        if st.session_state.tiempo_arias:
            for nombre, datos in st.session_state.tiempo_arias.items():
                with st.expander(f"📁 {nombre}", expanded=False):

                    st.write("Componentes disponibles:")
                    st.write(f"t5_pos:  {datos['Info']['t5_pos']},  t95_pos:  {datos['Info']['t95_pos']}")

                    # 🔴 botón con key única
                    if st.button("🔄 Borrar", key=f"Boton Borrar{nombre}"):

                        # borrar corrección SOLO de ese evento
                        if "evento_corregido" in st.session_state:
                            st.session_state["evento_corregido"].pop(nombre, None)

                        # borrar info guardada
                        del st.session_state.tiempo_arias[nombre]

                        st.warning(f"{nombre} eliminado")
                        st.rerun()

    with col_plot:

        col1, col2, col3 = st.columns(3)                                                               #CREAMOS TRES COLUMNAS PARA MOSTRAR DATOS
        with col1:
            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)

            ax.plot(t, a_N, label="Original N", lw=0.7)
            if corregir_N:
                ax.plot(t, a_corr_N, label="Corregido N", lw=0.7)
            ax.set_title("Acelerograma")
            ax.set_xlabel("Tiempo (s)")
            ax.set_ylabel("Aceleración (m/s²)")
            ax.legend(loc='upper right')
            ax.grid(True,lw=0.5)
            st.pyplot(fig)

            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)

            ax.plot(t, a_E, label="Original E", lw=0.7)
            if corregir_E:    
                ax.plot(t, a_corr_E, label="Corregido E", lw=0.7)
            ax.set_title("Acelerograma")
            ax.set_xlabel("Tiempo (s)")
            ax.set_ylabel("Aceleración (m/s²)")
            ax.legend(loc='upper right')
            ax.grid(True,lw=0.5)
            st.pyplot(fig)


        with col2:
            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)
            ax.plot(df_soil_dynamic_N["Tiempo (s)"], df_soil_dynamic_N["Velocidad (m/s)"], label='Velocidad Norte',lw=0.8)
            if corregir_N:
                ax.plot(df_soil_dynamic_N_corr["Tiempo (s)"], df_soil_dynamic_N_corr["Velocidad (m/s)"], label='Velocidad Norte corregida',lw=0.8)
            ax.set_title("Velocidad")
            ax.set_xlabel("Tiempo (s)")
            ax.set_ylabel("Velocidad (m/s)")
            ax.legend(loc='upper right')
            ax.grid(True,lw=0.5)
            st.pyplot(fig)

            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)
            ax.plot(df_soil_dynamic_E["Tiempo (s)"], df_soil_dynamic_E["Velocidad (m/s)"], label='Velocidad Este',lw=0.8)
            if corregir_E:
                ax.plot(df_soil_dynamic_E_corr["Tiempo (s)"], df_soil_dynamic_E_corr["Velocidad (m/s)"], label='Velocidad Este corregida',lw=0.8)
            ax.set_title("Velocidad")
            ax.set_xlabel("Tiempo (s)")
            ax.set_ylabel("Velocidad (m/s)")
            ax.legend(loc='upper right')
            ax.grid(True,lw=0.5)
            st.pyplot(fig)
        with col3:
            
            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)
            ax.plot(df_soil_dynamic_N["Tiempo (s)"], df_soil_dynamic_N["Desplazamiento (m)"], label='Desplazamiento Norte',lw=0.8)
            if corregir_N:
                ax.plot(df_soil_dynamic_N_corr["Tiempo (s)"], df_soil_dynamic_N_corr["Desplazamiento (m)"], label='Desplazamiento Norte Corregido',lw=0.8)           
            ax.set_title("Desplazamiento")
            ax.set_xlabel("Tiempo (s)")
            ax.set_ylabel("Desplazamiento (m)")
            ax.legend(loc='upper right')
            ax.grid(True,lw=0.5)
            st.pyplot(fig)

            fig, ax = plt.subplots(figsize=(8, 4), dpi=500)
            ax.plot(df_soil_dynamic_E["Tiempo (s)"], df_soil_dynamic_E["Desplazamiento (m)"], label='Desplazamiento Este',lw=0.8)
            if corregir_E:
                ax.plot(df_soil_dynamic_E_corr["Tiempo (s)"], df_soil_dynamic_E_corr["Desplazamiento (m)"], label='Desplazamiento Este Corregido',lw=0.8) 
            ax.set_title("Desplazamiento")
            ax.set_xlabel("Tiempo (s)")
            ax.set_ylabel("Desplazamiento (m)")
            ax.legend(loc='upper right')
            ax.grid(True,lw=0.5)
            st.pyplot(fig)

            #CORRECCION DE LINEA BASE 

        fig, ax = plt.subplots(figsize=(8, 2), dpi=500)
        ax.plot(arias_N_E["Tiempo (s)"], arias_N_E["Arias (m/s)"], label=f"Arias Bidireccional = {arias_N_E['Arias (m/s)'].iloc[-1]:.2f}",lw=0.8)
        ax.plot(arias_N["Tiempo (s)"], arias_N["Arias (m/s)"], label=f"Arias Norte = {arias_N['Arias (m/s)'].iloc[-1]:.2f}",lw=0.8)
        ax.plot(arias_E["Tiempo (s)"], arias_E["Arias (m/s)"], label=f"Arias Este = {arias_E['Arias (m/s)'].iloc[-1]:.2f}",lw=0.8)
    
        ax.axvline(t5, linestyle='--', linewidth=0.8, label=f"5%= {t5:.2f} seg")
        ax.axvline(t95, linestyle='--', linewidth=0.8, label=f"95%= {t95:.2f} seg")

        t = arias_N_E["Tiempo (s)"]
        arias = arias_N_E["Arias (m/s)"]
        mask = (t >= t5) & (t <= t95)
        ax.fill_between(
            t,
            arias,
            where=mask,
            alpha=0.3,
            label=f"Intensidad (5%-95%) = {(t95-t5):.2f} seg"
    )
        ax.set_title("Intensidad de Arias")
        ax.set_xlabel("Tiempo (s)")
        ax.set_ylabel("Intensidad (m/s)")
        ax.legend(loc='upper left',fontsize=6)
        ax.grid(True,lw=0.5)
        st.pyplot(fig,use_container_width=False) 