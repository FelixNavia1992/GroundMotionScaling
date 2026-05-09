import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, detrend
from scipy.fft import fft, fftfreq
from numba import njit
from functions_T1 import *
from functions_T2 import *
from functions_T3 import *
from functions_T4 import *

st.set_page_config(page_title="REG_SISM",page_icon="⚡",layout="wide",)
st.title("📊 Procesamiento de Registros Sísmicos - ASCE 7-22")
tab1, tab2, tab3, tab4 = st.tabs(["📈 Ingreso de Registros sÍsmico","🧠 Tratamiento de Datos","🖥️ Espectro objetivo ","📐 Direccionalidad"])

# =========================
# APLICACION
# =========================
with tab1:
    sismo = gestionar_eventos()                                                        #INGRESO DE REGISTROS SISMICOS PARA ANALISIS
with tab2:
    sismo_base_line=tratamiento_datos(sismo)
with tab3:   
    escalamiento_amplitud(sismo)
with tab4:
    direccionalidad(sismo)
