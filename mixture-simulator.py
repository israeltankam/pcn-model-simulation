#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import streamlit as st
import hydralit_components as hc
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from model import modelTrajectories, diseaseIncidence, finalDiseaseIncidence, cropYield, distinctCropYield, yieldOptimizer, displayOptimal, displayDiseaseDynamics

# Set page layout to centered and responsive
# st.set_page_config(layout="wide")
st.set_page_config(layout='wide',initial_sidebar_state='collapsed')


# specify the primary menu definition
menu_data = [
    {'icon': "far fa-chart-bar", 'label':"Simulation"},#no tooltip message
    {'icon': "fas fa-seedling", 'label':"Edit varieties"},
    {'icon': "fas fa-tachometer-alt", 'label':"About & Settings"},
]

over_theme = {'txc_inactive': '#FFFFFF', 'menu_background':'#85929E'}
st.markdown("# Cassava Mixture simulator")
main_tab= hc.nav_bar(
    menu_definition=menu_data,
    override_theme=over_theme,
    #home_name='Introduction',
    #login_name='Logout',
    hide_streamlit_markers=False, #will show the st hamburger as well as the navbar now!
    sticky_nav=True, #at the top or not
    sticky_mode='pinned', #jumpy or not-jumpy, but sticky or pinned
)


# Set Plant parameters

# Plant A parameters
st.session_state.setdefault("alpha_A", 0.01)  # Acquisition rate
st.session_state.setdefault("beta_A", 0.02)  # Inoculation rate
st.session_state.setdefault("yield_healthy_A", 40.0) # Average yield when healthy
st.session_state.setdefault("yield_diseased_A", 21.0) # Average yield when diseased
st.session_state.setdefault("gamma_A", 0.05)  # Latency speed

# Plant B parameters
st.session_state.setdefault("alpha_B", 0.02)  # Acquisition rate
st.session_state.setdefault("beta_B", 0.018)  # Inoculation rate
st.session_state.setdefault("yield_healthy_B", 44.0) # Average yield when healthy
st.session_state.setdefault("yield_diseased_B", 21.0) # Average yield when diseased
st.session_state.setdefault("gamma_B", 0.0667)  # Latency speed



# Insect parameters
st.session_state.setdefault("sigma", 0.01) # Dispersal parameter
st.session_state.setdefault("omega", 0.01) # Mortality parameter
st.session_state.setdefault("r", 0.01)  # Recovering parameter

# Insect abundance parameters
st.session_state.setdefault("f_very_low", 0.1) # Insect abundance per plant in very low insect pressure
st.session_state.setdefault("f_low", 1.0) # Insect abundance per plant in low insect pressure
st.session_state.setdefault("f_medium", 5.0) # Insect abundance per plant in medium insect pressure
st.session_state.setdefault("f_high", 40.0) # Insect abundance per plant in high insect pressure
st.session_state.setdefault("f_very_high", 200.0) # Insect abundance per plant in  very high insect pressure

# Disease pressure in plants parameters
st.session_state.setdefault("I_proportion_very_low", 0.01) # Initial infected plant proportion in very low pressure
st.session_state.setdefault("I_proportion_low", 0.1) # Initial infected plant proportion in low pressure
st.session_state.setdefault("I_proportion_medium", 0.2) # Initial infected plant proportion in medium pressure
st.session_state.setdefault("I_proportion_high", 0.5) # Initial infected plant proportion in high pressure
st.session_state.setdefault("I_proportion_very_high", 0.8) # Initial infected plant proportion in very high pressure

# Disease pressure in vectors parameters
st.session_state.setdefault("v_proportion_very_low", 0.01) # Initial viruliferous vector proportion in very low pressure
st.session_state.setdefault("v_proportion_low", 0.1) # Initial viruliferous vector proportion in low pressure
st.session_state.setdefault("v_proportion_medium", 0.2) # Initial viruliferous vector proportion in medium pressure
st.session_state.setdefault("v_proportion_high", 0.5) # Initial viruliferous vector proportion in high pressure
st.session_state.setdefault("v_proportion_very_high", 0.8) # Initial viruliferous vector proportion in very high pressure

# Cassava growing parameters
st.session_state.setdefault("K", 10000) # Field density
st.session_state.setdefault("T", 300) # Season duration

# Selected pressure parameters
st.session_state.setdefault("f", st.session_state.f_very_low) # Insect abundance per plant
st.session_state.setdefault("I_proportion", st.session_state.I_proportion_very_low) # Initial proportion of infected plant
st.session_state.setdefault("v_proportion", st.session_state.v_proportion_very_low) # Initial proportion of viruliferous vectors

step = 0.01
                            

# Set Streamlit app title
# st.title("Cassava mixture")

def main():
    col1, col2, col3 = st.columns([2, 10, 5])
    with col1:
        st.session_state.K = st.slider("Field density (plants/ha):", min_value=8000, max_value=15000, value=st.session_state.K, step=500)
    with col2:
        subcol1, subcol2, subcol3 = st.columns([1, 1, 1])
        with subcol1:
            insect_pressure_option_dic = {'Very low': 0, 'Low': 1, 'Medium': 2, 'High': 3, 'Very high': 4}
            selected_pressure = st.selectbox("Insect abundance (per plant)", options=list(insect_pressure_option_dic.keys()))
            if insect_pressure_option_dic[selected_pressure] == 0:
                st.session_state.f = st.session_state.f_very_low
            elif insect_pressure_option_dic[selected_pressure] == 1:
                st.session_state.f = st.session_state.f_low
            elif insect_pressure_option_dic[selected_pressure] == 2:
                st.session_state.f = st.session_state.f_medium
            elif insect_pressure_option_dic[selected_pressure] == 3:
                st.session_state.f = st.session_state.f_high
            elif insect_pressure_option_dic[selected_pressure] == 4:
                st.session_state.f = st.session_state.f_very_high
        with subcol2:
            plant_disease_pressure_option_dic = {'Very low': 0, 'Low': 1, 'Medium': 2, 'High': 3, 'Very high': 4}
            selected_pressure = st.selectbox("Disease pressure in plants", options=list(plant_disease_pressure_option_dic.keys()))
            if plant_disease_pressure_option_dic[selected_pressure] == 0:
                st.session_state.I_proportion = st.session_state.I_proportion_very_low
            elif plant_disease_pressure_option_dic[selected_pressure] == 1:
                st.session_state.I_proportion = st.session_state.I_proportion_low
            elif plant_disease_pressure_option_dic[selected_pressure] == 2:
                st.session_state.I_proportion = st.session_state.I_proportion_medium
            elif plant_disease_pressure_option_dic[selected_pressure] == 3:
                st.session_state.I_proportion = st.session_state.I_proportion_high
            elif plant_disease_pressure_option_dic[selected_pressure] == 4:
                st.session_state.I_proportion = st.session_state.I_proportion_very_high
        with subcol3:
            vector_disease_pressure_option_dic = {'Very low': 0, 'Low': 1, 'Medium': 2, 'High': 3, 'Very high': 4}
            selected_pressure = st.selectbox("Disease pressure in vectors", options=list(vector_disease_pressure_option_dic.keys()))
            if vector_disease_pressure_option_dic[selected_pressure] == 0:
                st.session_state.v_proportion = st.session_state.v_proportion_very_low
            elif vector_disease_pressure_option_dic[selected_pressure] == 1:
                st.session_state.v_proportion = st.session_state.v_proportion_low
            elif vector_disease_pressure_option_dic[selected_pressure] == 2:
                st.session_state.v_proportion = st.session_state.v_proportion_medium
            elif vector_disease_pressure_option_dic[selected_pressure] == 3:
                st.session_state.v_proportion = st.session_state.v_proportion_high
            elif vector_disease_pressure_option_dic[selected_pressure] == 4:
                st.session_state.v_proportion = st.session_state.v_proportion_very_high
    with col3:
        st.session_state.T = st.slider("Season duration (days):", min_value=150, max_value=365, value=st.session_state.T, step=1)
        
   
    st.markdown("<hr>", unsafe_allow_html=True)
        
    
    col1, col2 = st.columns([8,12])
    
    ##################################################################################################
    # plotting
    theta , _ = yieldOptimizer(st.session_state)
    with col1:
        displayOptimal(theta, st.session_state)
    with col2:
        displayDiseaseDynamics(theta, st.session_state)
    
###################################################################################
if main_tab == "Simulation":
    main()
    
if main_tab == "Edit varieties":
    _, col1, col2, _ = st.columns([1, 5, 5, 1])
    with col1:
        st.markdown("### Variety A")
        st.session_state.alpha_A= st.slider("Acquisition rate A", min_value=0.0, max_value=1.0, value=st.session_state.alpha_A, step=0.01)
        st.session_state.beta_A= st.slider("Inoculation rate A", min_value=0.0, max_value=1.0, value=st.session_state.beta_A, step=0.01)
        st.session_state.gamma_A= 1/st.slider("Latency duration A(days)", min_value=0, max_value=20, value=int(1/st.session_state.gamma_A), step=1)
        st.session_state.yield_healthy_A= st.slider("Yield when healthy A(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_healthy_A, step=0.5)
        st.session_state.yield_diseased_A= st.slider("Yield when infected A(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_diseased_A, step=0.5)
    with col2:
        st.markdown("### Variety B")
        st.session_state.alpha_B= st.slider("Acquisition rate B", min_value=0.0, max_value=1.0, value=st.session_state.alpha_B, step=0.01)
        st.session_state.beta_B= st.slider("Inoculation rate B", min_value=0.0, max_value=1.0, value=st.session_state.beta_B, step=0.01)
        st.session_state.gamma_B= 1/st.slider("Latency duration B(days)", min_value=0, max_value=20, value=int(1/st.session_state.gamma_B), step=1)
        st.session_state.yield_healthy_B= st.slider("Yield when healthy B(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_healthy_B, step=0.5)
        st.session_state.yield_diseased_B= st.slider("Yield when infected B(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_diseased_B, step=0.5)
    
elif main_tab == "About & Settings":
    col1, col2, _ = st.columns([1, 5, 4])
    with col1:
        st.markdown("# About")
    with col2:
        st.markdown("# ")
        st.markdown("# ")
        st.markdown("- The model consider a mixture of two cassava varieties named 'Cultivar A' and 'Cultivar B' ")
        st.markdown("- Each cultivar has an intrinsic virus inoculation rate, virus acquisition rate and virus incubaation duration")
        st.markdown("- We name 'Disease pressure in plants' the initial fraction of plants that are infeted and assume it uniformelly distributed in both varieties")
        st.markdown("- We assume that initially they are no cryptic cassava plants. They are either infected of susceptible.")
        st.markdown("- The vector abundance per plant is constant and defined by the 'Insect pressure'")
        st.markdown("- Among the vectors, those that are initially viruliferous are defined by the 'Disease pressure in vectors' ")
        st.markdown("- We assume an initially equal distribution between those who acquired the virus on A cultivars and B cultivars")
    col1, col2, _ = st.columns([1, 5, 4])
    with col1:
        st.markdown("# Settings")
    with col2:
        st.markdown("# ")
        st.markdown("### Insect parameters ")
        st.session_state.sigma = st.slider("Insect dispersal rate $\sigma$:", min_value=0.0, max_value=1.0, value=st.session_state.sigma, step=0.01)
        st.session_state.omega = st.slider("Insect mortality rate $\omega$:", min_value=0.0, max_value=1.0, value=st.session_state.omega, step=0.01)
        st.session_state.r = st.slider("Insect recovery rate $r$:", min_value=0.0, max_value=1.0, value=st.session_state.r, step=0.01)
        
        st.markdown("### Insect abundance setup")
        st.session_state.f_very_low = st.slider("Insect abundance per plant in very low insect pressure:", min_value=0.0, max_value=5.0, value=st.session_state.f_very_low, step=0.1)
        st.session_state.f_low = st.slider("Insect abundance per plant in low insect pressure:", min_value=0.0, max_value=10.0, value=st.session_state.f_low, step=0.1)
        st.session_state.f_medium = st.slider("Insect abundance per plant in medium insect pressure:", min_value=0.0, max_value=50.0, value=st.session_state.f_medium, step=1.0)
        st.session_state.f_high = st.slider("Insect abundance per plant in high insect pressure:", min_value=0.0, max_value=200.0, value=st.session_state.f_high, step=1.0)
        st.session_state.f_very_high = st.slider("Insect abundance per plant in  very high insect pressure:", min_value=0.0, max_value=500.0, value=st.session_state.f_very_high, step=1.0)
        
        st.markdown("### Disease pressure in plant setup")
        st.session_state.I_proportion_very_low = st.slider("Initial infected plant proportion in very low pressure:", min_value=0.0, max_value=1.0, value=st.session_state.I_proportion_very_low, step=0.01)
        st.session_state.I_proportion_low = st.slider("Initial infected plant proportion in low pressure:", min_value=0.0, max_value=1.0, value=st.session_state.I_proportion_low, step=0.01)
        st.session_state.I_proportion_medium = st.slider("Initial infected plant proportion in medium pressure:", min_value=0.0, max_value=1.0, value=st.session_state.I_proportion_medium, step=0.01)
        st.session_state.I_proportion_high = st.slider("Initial infected plant proportion in high pressure:", min_value=0.0, max_value=1.0, value=st.session_state.I_proportion_high, step=0.01)
        st.session_state.I_proportion_very_high = st.slider("Initial infected plant proportion in very high pressure:", min_value=0.0, max_value=1.0, value=st.session_state.I_proportion_very_high, step=0.01)
        
        st.markdown("### Disease pressure in vector setup")
        st.session_state.v_proportion_very_low = st.slider("Initial viruliferous vector proportion in very low pressure:", min_value=0.0, max_value=1.0, value=st.session_state.v_proportion_very_low, step=0.01)
        st.session_state.v_proportion_low = st.slider("Initial viruliferous vector proportion in low pressure:", min_value=0.0, max_value=1.0, value=st.session_state.v_proportion_low, step=0.01)
        st.session_state.v_proportion_medium = st.slider("Initial viruliferous vector proportion in medium pressure:", min_value=0.0, max_value=1.0, value=st.session_state.v_proportion_medium, step=0.01)
        st.session_state.v_proportion_high = st.slider("Initial viruliferous vector proportion in high pressure:", min_value=0.0, max_value=1.0, value=st.session_state.v_proportion_high, step=0.01)
        st.session_state.v_proportion_very_high = st.slider("Initial viruliferous vector proportion in very high pressure:", min_value=0.0, max_value=1.0, value=st.session_state.v_proportion_very_high, step=0.01)
