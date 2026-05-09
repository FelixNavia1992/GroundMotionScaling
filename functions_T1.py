import streamlit as st
import numpy as np
import pandas as pd

def cargar_registro_individual(file):                      #FUNCION QUE DETECTA EL DT DE CADA REGISTRO CON FORMATO AT2

    nombre_archivo = file.name

    content = file.read().decode("utf-8", errors="ignore")
    lines = content.splitlines()

    dt = None
    for line in lines:
        if "DT=" in line:
            try:
                dt = float(line.split("DT=")[1].split()[0])
                break
            except:
                continue

    if dt is None:
        raise ValueError("No se encontró DT")

    acc_values = []

    for line in lines:
        try:
            nums = [float(x) for x in line.split()]
            acc_values.extend(nums)
        except:
            continue

    acc_values = np.array(acc_values)
    n = len(acc_values)

    time = np.arange(0, n * dt, dt)
    time = time[:n]

    df = pd.DataFrame({
        "Tiempo (s)": time,
        "Aceleración (g)": acc_values
    })

    return df, dt, nombre_archivo,n              #EL DATAFRAME INCLUYE TIEMPO Y ACELERACION


def gestionar_eventos():

    st.header("📦 Gestión de eventos sísmicos")

    # Inicializar almacenamiento
    if "eventos" not in st.session_state:
        st.session_state.eventos = {}

    st.subheader("➕ Crear nuevo evento")

    nombre_evento = st.text_input("Nombre del evento sísmico")                               #CREA EL NOMBRE DEL EVENTO SISIMICO
    col1, col2 = st.columns(2)                                                               #CREAMOS DOS COLUMNAS PARA INGRESO DE REGISTROS
    with col1:
        file_N = st.file_uploader("Registro Norte (N)", type=["AT2"], key="N")
    with col2:
        file_E = st.file_uploader("Registro Este (E)", type=["AT2"], key="E")

    if st.button("💾 Guardar evento"):
        if not nombre_evento:
            st.warning("⚠️ Ingresa un nombre de evento")                                    #CONDICIONES
            return
        
        if file_N is None or file_E is None:
            st.warning("⚠️ Debes cargar ambos registros N y E")
            return
        try:
            df_N, dt_N, _,n_N = cargar_registro_individual(file_N)
            df_E, dt_E, _,n_E = cargar_registro_individual(file_E)
            if n_N!=n_E:
                st.error("⚠️ Deben tener la misma cantidad de datos")
                return
        except:
            st.error("❌ Error al leer archivos")
            return


        st.session_state.eventos[nombre_evento] = {                                         #GUARDAR EN MEMORIA PARA STREAMLIT
            "N": {"df": df_N, "dt": dt_N, "npts":n_N},
            "E": {"df": df_E, "dt": dt_E, "npts":n_E}
        }

        st.success(f"✅ Evento '{nombre_evento}' guardado")

    # =========================
    # MOSTRAR EVENTOS (VENTANA)
    # =========================
    st.subheader("📊 Paquete de eventos")

    if st.session_state.eventos:

        for nombre, datos in st.session_state.eventos.items():

            with st.expander(f"📁 {nombre}", expanded=False):

                st.write("Componentes disponibles:")
                st.write("✔️ Norte (N)")
                st.write("✔️ Este (E)")

                # Mostrar info básica
                st.write(f"dt N: {datos['N']['dt']:.5f}, NPTs: {datos['N']['npts']}")
                st.write(f"dt E: {datos['E']['dt']:.5f}, NPTs: {datos['E']['npts']}")

                # Botón eliminar
                if st.button(f"🗑️ Eliminar {nombre}", key=f"del_{nombre}"):
                    del st.session_state.eventos[nombre]
                    del st.session_state.tiempo_arias[nombre]
                    st.session_state["resultados_espectro_ROTD100"].pop(nombre, None)
                    st.session_state["resultados_espectro_ROTD100_escalado"].pop(nombre, None)
                    st.rerun()

    else:
        st.info("No hay eventos guardados")
    return st.session_state.eventos