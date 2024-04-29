import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import streamlit as st

def modelTrajectories(theta, session_state):
    # Importing parameters
    I_proportion = session_state.I_proportion 
    f = session_state.f 
    K = session_state.K
    T = session_state.T
    sigma = session_state.sigma
    omega = session_state.omega
    r = session_state.r
    alpha_A = session_state.alpha_A
    alpha_B = session_state.alpha_B
    beta_A = session_state.beta_A
    beta_B = session_state.beta_B
    gamma_A = session_state.gamma_A
    gamma_B = session_state.gamma_B
    # Initial values
    lA_0 = 0
    iA_0 = I_proportion * theta
    lB_0 = 0
    iB_0 = I_proportion * (1-theta)
    F = f * K
    VA_0 = session_state.v_proportion * F/2
    VB_0 = VA_0
    
    # Define the ODE system
    def ode_system(t, y):
        lA, iA, lB, iB, VA, VB = y
        psi = 1/(sigma + omega + r)
        dlAdt = (psi*sigma/K)*beta_A*(theta - lA - iA)*(VA + VB) - gamma_A*lA
        diAdt = gamma_A * lA
        dlBdt = (psi*sigma/K)*beta_B*(1 - theta - lB - iB)*(VA + VB) - gamma_B*lB
        diBdt = gamma_B * lB
        dVAdt = alpha_A*(iA*F - psi*(sigma*iA*(VA+VB) + (omega+r)*VA)) - (omega + r)*VA
        dVBdt = alpha_B*(iB*F - psi*(sigma*iB*(VA+VB) + (omega+r)*VB)) - (omega + r)*VB
        return [dlAdt, diAdt, dlBdt, diBdt, dVAdt, dVBdt]
    
    # Solve the ODE
    sol = solve_ivp(ode_system, [0, T], [lA_0, iA_0, lB_0, iB_0, VA_0, VB_0], t_eval=np.linspace(0, T, 100))
    
    return sol

def diseaseIncidence(theta, session_state):
    sol = modelTrajectories(theta, session_state)
    
    # Get the solution
    t_values = sol.t
    disease_incidence  = sol.y[1] + sol.y[3]

    return t_values, disease_incidence

def finalDiseaseIncidence(theta, session_state):
    t_values, disease_incidence = diseaseIncidence(theta)
    return disease_incidence[-1]

def cropYield(theta, session_state):
    sol = modelTrajectories(theta, session_state)
    
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
    y = session_state.yield_healthy_A*(sA+lA) + session_state.yield_diseased_A*iA + session_state.yield_healthy_B*(sB+lB) + session_state.yield_diseased_B*iB
    return y
def distinctCropYield(theta, session_state):
    sol = modelTrajectories(theta,session_state)
    
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
    yA = session_state.yield_healthy_A*(sA+lA) + session_state.yield_diseased_A*iA 
    yB = session_state.yield_healthy_B*(sB+lB) + session_state.yield_diseased_B*iB
    return yA, yB

def yieldOptimizer(session_state):
    # Define bounds for theta
    theta_bounds = (0, 1)

    # Define the negative function function
    def neg_cropYield(theta):
        return -cropYield(theta, session_state)

    # Minimize the negative function function
    result = minimize_scalar(neg_cropYield, bounds=theta_bounds, method='bounded')
    return result.x, -result.fun

def displayOptimal(theta, session_state):
# Calculate the values
    percentageA = theta * 100
    percentageB = (1 - theta) * 100
    yieldA, yieldB = distinctCropYield(theta, session_state)

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

def displayDiseaseDynamics(theta, session_state):
    sol = modelTrajectories(theta, session_state)

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
