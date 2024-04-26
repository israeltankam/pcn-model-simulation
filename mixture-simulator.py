#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import streamlit as st
import hydralit_components as hc
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# Set page layout to centered and responsive
# st.set_page_config(layout="wide")
st.set_page_config(layout='wide',initial_sidebar_state='collapsed')


# specify the primary menu definition
menu_data = [
    {'icon': "far fa-chart-bar", 'label':"Simulation"},#no tooltip message
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
st.session_state.setdefault("gamma_B", 0.05)  # Latency speed



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
        
    
    col1, col2, _ = st.columns([2, 10, 5])
    with col1:
        #Create varieties
        st.markdown("<hr>", unsafe_allow_html=True)
        st.markdown("### Variety A")
        st.session_state.alpha_A= st.slider("Acquisition rate A", min_value=0.0, max_value=1.0, value=st.session_state.alpha_A, step=0.01)
        st.session_state.beta_A= st.slider("Inoculation rate A", min_value=0.0, max_value=1.0, value=st.session_state.beta_A, step=0.01)
        st.session_state.gamma_A= 1/st.slider("Latency duration A(days)", min_value=0, max_value=20, value=int(1/st.session_state.gamma_A), step=1)
        st.session_state.yield_healthy_A= st.slider("Yield when healthy A(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_healthy_A, step=0.5)
        st.session_state.yield_diseased_A= st.slider("Yield when infected A(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_diseased_A, step=0.5)
        st.markdown("<hr>", unsafe_allow_html=True)
        st.markdown("### Variety B")
        st.session_state.alpha_B= st.slider("Acquisition rate B", min_value=0.0, max_value=1.0, value=st.session_state.alpha_B, step=0.01)
        st.session_state.beta_B= st.slider("Inoculation rate B", min_value=0.0, max_value=1.0, value=st.session_state.beta_B, step=0.01)
        st.session_state.gamma_B= 1/st.slider("Latency duration B(days)", min_value=0, max_value=20, value=int(1/st.session_state.gamma_B), step=1)
        st.session_state.yield_healthy_B= st.slider("Yield when healthy B(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_healthy_B, step=0.5)
        st.session_state.yield_diseased_B= st.slider("Yield when infected B(ton/ha)", min_value=0.0, max_value=100.0, value=st.session_state.yield_diseased_B, step=0.5)
    
    ################################################################################################
    def modelTrajectories(theta):
        # Initial values
        lA_0 = 0
        iA_0 = st.session_state.I_proportion * theta
        lB_0 = 0
        iB_0 = st.session_state.I_proportion * (1-theta)
        F = st.session_state.f * st.session_state.K
        VA_0 = st.session_state.v_proportion * F/2
        VB_0 = VA_0
        
        # Define the ODE system
        def ode_system(t, y):
            lA, iA, lB, iB, VA, VB = y
            psi = 1/(st.session_state.sigma + st.session_state.omega + st.session_state.r)
            dlAdt = (psi*st.session_state.sigma/st.session_state.K)*st.session_state.beta_A*(theta - lA - iA)*(VA + VB) - st.session_state.gamma_A*lA
            diAdt = st.session_state.gamma_A * lA
            dlBdt = (psi*st.session_state.sigma/st.session_state.K)*st.session_state.beta_B*(1 - theta - lB - iB)*(VA + VB) - st.session_state.gamma_B*lB
            diBdt = st.session_state.gamma_B * lB
            dVAdt = st.session_state.alpha_A*(iA*F - psi*(st.session_state.sigma*iA*(VA+VB) + (st.session_state.omega+st.session_state.r)*VA)) - (st.session_state.omega + st.session_state.r)*VA
            dVBdt = st.session_state.alpha_B*(iB*F - psi*(st.session_state.sigma*iB*(VA+VB) + (st.session_state.omega+st.session_state.r)*VB)) - (st.session_state.omega + st.session_state.r)*VB
            return [dlAdt, diAdt, dlBdt, diBdt, dVAdt, dVBdt]
        
        # Solve the ODE
        sol = solve_ivp(ode_system, [0, st.session_state.T], [lA_0, iA_0, lB_0, iB_0, VA_0, VB_0], t_eval=np.linspace(0, st.session_state.T, 100))
        
        return sol
        
    #################################################################################################
    def diseaseIncidence(theta):
        sol = modelTrajectories(theta)
        
        # Get the solution
        t_values = sol.t
        disease_incidence  = sol.y[1] + sol.y[3]

        return t_values, disease_incidence
    ##################################################################################################
    def finalDiseaseIncidence(theta):
        t_values, disease_incidence = diseaseIncidence(theta)
        return disease_incidence[-1]
    
    ##################################################################################################
    def cropYield(theta):
        sol = modelTrajectories(theta)
        
        # Get the solution
        t_values = sol.t
        lA_values, iA_values, lB_values, iB_values, VA_values, VB_values = sol.y
        lA = lA_values[-1]
        iA = iA_values[-1]
        sA = theta - (lA + iA)
        lB = lB_values[-1]
        iB = iB_values[-1]
        sB = (1-theta) - (lB + iB)
        
        # Calculate yield per hectare
        y = st.session_state.yield_healthy_A*(sA+lA) + st.session_state.yield_diseased_A*iA + st.session_state.yield_healthy_B*(sB+lB) + st.session_state.yield_diseased_B*iB
        return y
    
    ##################################################################################################
    def distinctCropYield(theta):
        sol = modelTrajectories(theta)
        
        # Get the solution
        t_values = sol.t
        lA_values, iA_values, lB_values, iB_values, VA_values, VB_values = sol.y
        lA = lA_values[-1]
        iA = iA_values[-1]
        sA = theta - (lA + iA)
        lB = lB_values[-1]
        iB = iB_values[-1]
        sB = (1-theta) - (lB + iB)
        
        # Calculate yield per hectare
        yA = st.session_state.yield_healthy_A*(sA+lA) + st.session_state.yield_diseased_A*iA 
        yB = st.session_state.yield_healthy_B*(sB+lB) + st.session_state.yield_diseased_B*iB
        return yA, yB
    
    ##################################################################################################
    def yieldOptimizer():
        # Define bounds for theta
        theta_bounds = (0, 1)

        # Define the negative function function
        def neg_cropYield(theta):
            return -cropYield(theta)

        # Minimize the negative function function
        result = minimize_scalar(neg_cropYield, bounds=theta_bounds, method='bounded')
        return result.x, -result.fun
    
    ##################################################################################################
    def displayOptimal(theta):
        # Calculate the values
        percentageA = theta * 100
        percentageB = (1 - theta) * 100
        yieldA, yieldB = distinctCropYield(theta)

        total_yield = yieldA + yieldB

        # Plotting
        sizes = [theta, (1-theta)]

        fig, ax = plt.subplots(figsize=(9.6, 5), subplot_kw=dict(aspect="equal"))  # Adjusted width to 60% of A4 paper width

        wedges, texts, autotexts = ax.pie(sizes, labels=['', ''], autopct='%1.1f%%', startangle=140, explode=(0.1, 0), shadow=True)

        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        kw = dict(arrowprops=dict(arrowstyle="-"),
                  bbox=bbox_props, zorder=0, va="center")

        tooltips = [r'$\bf{Cultivar \ A}$' + f'\n{percentageA:.2f} %\nyield = {yieldA:.2f}' + ' ton/ha',
                    r'$\bf{Cultivar \ B}$' + f'\n{percentageB:.2f} %\nyield = {yieldB:.2f}' + ' ton/ha']

        for i, p in enumerate(wedges):
            ang = (p.theta2 - p.theta1)/2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = f"angle,angleA=0,angleB={ang}"
            kw["arrowprops"].update({"connectionstyle": connectionstyle})
            ax.annotate(tooltips[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                        horizontalalignment=horizontalalignment, fontsize=10, color='black', **kw)

        # Add boxed message for total yield on the right
        ax.text(1.2, 0, f'Total Yield = {total_yield:.2f}' + ' ton/ha', fontsize=12, color='black', ha='left', va='center', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.5'))

        # Display the plot in Streamlit
        st.pyplot(fig)
        
    ##################################################################################################
    def displayDiseaseDynamics(theta):
        sol = modelTrajectories(theta)

        # Access the solution
        t_values = sol.t
        lA_values, iA_values, lB_values, iB_values, VA_values, VB_values = sol.y

        # Plot the solutions
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # First subplot
        axes[0].plot(t_values, iA_values, label='Infected A cultivar')
        axes[0].plot(t_values, iB_values, label='Infected B cultivar')
        axes[0].set_xlabel('Time')
        axes[0].set_ylabel('Proportions')
        axes[0].set_title('Disease Dynamics Over Time')
        axes[0].legend()
        axes[0].grid(True)

        # Second subplot
        axes[1].plot(t_values, VA_values, label='Acquired on A plants')
        axes[1].plot(t_values, VB_values, label='Acquired on B plants')
        axes[1].set_xlabel('Time')
        axes[1].set_ylabel('Populations')
        axes[1].set_title('Infected Insect Dynamics Over Time')
        axes[1].legend()
        axes[1].grid(True)

        fig.tight_layout()

        # Display the figure in Streamlit
        st.pyplot(fig)
    
    ##################################################################################################
    with col2:
        # plotting
        theta , _ = yieldOptimizer()
        displayOptimal(theta)
        displayDiseaseDynamics(theta)
    
###################################################################################
if main_tab == "Simulation":
    main()

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
