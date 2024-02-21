#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import streamlit as st
import hydralit_components as hc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Set page layout to centered and responsive
# st.set_page_config(layout="wide")
st.set_page_config(layout='wide',initial_sidebar_state='collapsed')


# specify the primary menu definition
menu_data = [
    {'icon': "far fa-copy", 'label':"Model & Parameters"},
    {'icon': "far fa-chart-bar", 'label':"Simulation"},#no tooltip message
    #{'icon': "fas fa-chart-line", 'label':"Genetic Drift"},
    {'icon': "fas fa-tachometer-alt", 'label':"Edit Parameters"},
    #{'icon': "fas fa-download", 'label':"Download report"},
]

over_theme = {'txc_inactive': '#FFFFFF', 'menu_background':'#85929E'}
st.markdown("# Nemo")
main_tab= hc.nav_bar(
    menu_definition=menu_data,
    override_theme=over_theme,
    home_name='Introduction',
    #login_name='Logout',
    hide_streamlit_markers=False, #will show the st hamburger as well as the navbar now!
    sticky_nav=True, #at the top or not
    sticky_mode='pinned', #jumpy or not-jumpy, but sticky or pinned
)


# Define default parameter values
st.session_state.setdefault("init_infest", 100)
st.session_state.setdefault("s", 0.25)
st.session_state.setdefault("m", 0.35)
st.session_state.setdefault("v_freq", 1/(2*(1-st.session_state.m)))
st.session_state.setdefault("h", 0.17)
st.session_state.setdefault("mu", 0.10)
st.session_state.setdefault("c", 0.4)
st.session_state.setdefault("e", 500)
st.session_state.setdefault("detection_threshold", 1)
st.session_state.setdefault("num_years", 15)
st.session_state.setdefault("bc", 0.45)
st.session_state.setdefault("r", 4)
st.session_state.setdefault("num_gen", int(np.ceil(st.session_state.num_years/(st.session_state.r + 1))))
step = 0.01
                            

# Define parameter values for reset
st.session_state.setdefault("reset_init_infest", 100)
st.session_state.setdefault("reset_s", 0.25)
st.session_state.setdefault("reset_m", 0.35)
st.session_state.setdefault("reset_v_freq", 1/(2*(1-st.session_state.reset_m)))
st.session_state.setdefault("reset_h", 0.17)
st.session_state.setdefault("reset_mu", 0.10)
st.session_state.setdefault("reset_c", 0.4)
st.session_state.setdefault("reset_e", 500)
st.session_state.setdefault("reset_detection_threshold", 1)
st.session_state.setdefault("reset_num_years", 15)
st.session_state.setdefault("reset_bc", 0.45)
st.session_state.setdefault("reset_r", 4)

# Set Streamlit app title
#st.title("Nemo")
def dec2(x):
    d = round(x * 100) / 100
    return d
def dec1(number):
    return round(number, 1)   
def minimal4(vec):
    m = np.inf
    i = 0
    while i < len(vec) and vec[i] > 4:
        i += 1
    m = i
    return m
    
def generate_main_plot(N,S,V,R):
    fig, ax = plt.subplots(figsize=(14, 10), dpi=100)
    nb_gen = len(N)
    d = st.session_state.detection_threshold
    ax.plot([0, nb_gen-1], [d, d], 'k--', label='Acceptance threshold')
    ax.plot(np.arange(0, nb_gen), S, '-r', linewidth=3, marker='o', label='Susceptible or blocking resistance')
    ax.plot(np.arange(0, nb_gen), N, '-b', linewidth=3, marker='o', label='Masculinizing Resistance')
    ax.set_xlabel("Time (years)", fontsize=20)
    ax.set_ylabel("Nematode population density (eggs/g of soil)", fontsize=20)
    ax.set_xlim([0, nb_gen-1])
    ax.set_yscale('log')

    # X-axis ticks
    x_ticks_locations = np.arange(0, nb_gen, 1)
    x_ticks_labels = [str((st.session_state.r+1) * i) for i in range(nb_gen)]

    # Y-axis ticks
    y_ticks_locations = [10**(-2), 10**(-1), 10**0, d, 10**1, 10**2, 10**3]
    y_ticks_labels = [str(val) for val in y_ticks_locations]

    # Set X-axis ticks
    ax.set_xticks(x_ticks_locations)
    ax.set_xticklabels(x_ticks_labels)

    # Set Y-axis ticks
    ax.set_yticks(y_ticks_locations)
    ax.set_yticklabels(y_ticks_labels)

    # Legend
    ax.legend(fontsize=15, loc='upper right')

    # Create two columns with widths in the ratio 2:1
    col1, col2 = st.columns([2, 1])
    with col1:
        # Display the plot in Streamlit
        st.pyplot(fig)
    with col2:
        # Basic reproduction number
        st.markdown("Reproduction number $R'$ = "+str(dec2(R)))
        st.markdown("Effective reproduction number $R_0 = (1-m)R' = $ "+str(dec2(R*(1-st.session_state.m))))
        # Upper plot
        with st.expander("Progression of virulence $v$"):
            # Create a new figure and axes
            fig_upper, ax_upper = plt.subplots(figsize=(8, 5), dpi=100)

            # Plot the upper plot data
            ax_upper.plot(np.arange(0, nb_gen), V, linewidth=3)
            ax_upper.set_xlabel("Year", fontsize=30)
            ax_upper.tick_params(axis='both', which='major', labelsize=30)
            x_upper_ticks_locations = np.arange(0, nb_gen, 1)
            x_upper_ticks_labels = [str((st.session_state.r+1) * i) for i in range(nb_gen)]
            ax_upper.set_xticks(x_upper_ticks_locations)
            ax_upper.set_xticklabels(x_upper_ticks_labels)
            # Display the upper plot using Streamlit's pyplot function
            st.pyplot(fig_upper)
        with st.expander("Minimal control required with masculinizing resistance"):
            step = 0.001
            R = st.session_state.e*st.session_state.s*(1-st.session_state.reset_mu)*(1-st.session_state.reset_h)
            b = np.arange(step, 1, step)
            jj = np.ceil(np.log(2. / (R * (1 - b))) / np.log((1-st.session_state.mu) * (1 - st.session_state.h) * (1 - b)))
            # Plotting
            fig_pcb, ax_pcb = plt.subplots(figsize=(14, 10), dpi=100)
            # Plotting grey rectangle
            ax_pcb.fill([b[0], b[-1], b[-1], b[0]], [4, 4, jj[0], jj[0]], color='grey', alpha=0.5, edgecolor='none')
            # Find minimal position
            pos = minimal4(jj) * step
            # Plot jj
            ax_pcb.plot(b, jj, 'b-', linewidth=3)

            # Plot vertical dashed line at minimal position
            ax_pcb.plot([pos, pos], [0, 4], 'r--', linewidth=1)

            # Set plot properties
            ax_pcb.set_xlim([0, 1])
            ax_pcb.set_ylim([0, np.max(jj)])
            ax_pcb.set_ylabel("Minimal number of rotations $r_{min}$", fontsize=30)
            ax_pcb.set_xlabel("Biocontrol efficacy $b$", fontsize=30)
            ax_pcb.tick_params(axis='both', which='major', labelsize=30)
            st.pyplot(fig_pcb)

# Main tab
#with st.sidebar:
#    main_tab = st.radio("Navigation", ["Introduction", "Model & Parameters", "Simulation", "Settings"])

if main_tab == "Introduction":
    col1, col2, col3 = st.columns([1, 8, 1])
    with col2:
        st.markdown("# Introduction")
        st.markdown("The pale cyst nematode, Globodera pallida, is a pest that poses a significant threat to potato crops worldwide. The most effective chemical nematicides are toxic to non-target organisms and are now banned. Alternative control methods are therefore required. Crop rotation and biological control methods have limitations for effectively managing nematodes. The use of genetically resistant cultivars is a promising alternative, but nematode populations evolve, and virulent mutants can break resistance after just a few years. Masculinizing resistances, preventing avirulent nematodes from producing females, are thought to be more durable than blocking resistances, preventing infection. Our demo-genetic model, tracking both nematode population densities and their genetic frequencies, shows that virulence against masculinizing resistance may not fix in the pest population, under realistic conditions. Avirulence may persist despite the uniform use of resistance. This is because avirulent male nematodes may transmit avirulence to their progeny by mating with virulent females. Additionally, because avirulent nematodes do not produce females themselves, they weaken the reproductive rate of the nematode population, leading to a reduction in its density by at least 20\%. This avirulence load can even lead to the collapse of the nematode population in theory. Overall, our model shows that combining masculinizing resistance, rotation, and biocontrol may achieve durable suppression of G. pallida in a reasonable time frame. Our work is supported by an online interactive interface allowing users to test their own control combinations.")    
        st.markdown("# Hightlights")
        st.markdown("- Globodera pallida, commonly known as the Pale Cyst Nematode, poses a significant threat to potato crops on a global scale, making it a serious quarantine pest.")
        st.markdown("- Although utilizing resistant potato cultivars is a widely accepted and sustainable pest control strategy, the continuous evolution of nematode populations towards virulence can undermine the long-term efficacy of resistance-based control methods")
        st.markdown("- Masculinizing resistance aims to prevent avirulent nematodes from generating females, potentially leading to the elimination of avirulent nematodes from the population.")
        st.markdown("- Despite the selection pressure, tracking genotypic frequencies in real conditions reveals that the long-term fixation of the virulence allele does not necessarily occur.")
        st.markdown("- Avirulent nematodes, exclusively male, manage to survive as heterozygotes by mating with virulent females, thereby weakening the nematode's reproductive number.")
        st.markdown("- The biocontrol efficcacy required for long-term nematode suppression is lower when deploying masculinizing resistant plants, but there is a potential challenge of achieving the high biocontrol efficacy needed.")
        st.markdown("- Rotations have been established as an efficient and sustainable method for nematode control, but their implementation necessitates extended periods of cultivating non-host crops or leaving the soil bare.")
        st.markdown("- An effective solution for expediting nematode suppression involves combining resistant cultivars with biocontrol methods and rotations.")
        st.markdown("- The model presented in this study concurrently monitors nematode genetics and dynamics, providing insights into the selection for virulence and determining the required size of rotations and biocontrol efficacy for successful suppression.")

elif main_tab == "Model & Parameters":
    st.markdown("# Model & parameters")
    image1_path = "figs/diagram.png"
    st.image(image1_path)
    checkbox = st.checkbox("See the simple models")
    if checkbox:
        st.markdown("$X_n \longrightarrow AA$ PCNs, $\qquad Y_n \longrightarrow Aa$ PCNs, $\qquad Z_n \longrightarrow aa$ PCNs")
        st.markdown("- When susceptible plants are deployed in every generation, the overall PCN population $N_k$ at generation $k$ is given by the law:")
        st.latex(r'''
        \begin{equation*}
             N_{k+1} = (1-m)R \displaystyle\frac{ N_k}{1 + c N_k}
        \end{equation*}
        ''')
        st.markdown("Where $R$ is the PCN reproduction number and $c>0$ is an interspecific competition parameter. The term $(1-m)R$ is the basic reproduction number of the parasite.")
        st.markdown("If $(1-m)R>1$, the parasite population grows until it reaches carrying capacity")
        st.latex(r'''
            \begin{equation}
            K(R,m,c)=\frac{(1-m)R-1}{c}\,.
            \end{equation}
        ''')
        st.markdown("- When resistance plants are deployed in every generation, the overall PCN population $N_k$ at generation $k$ is given by the law:")
        st.latex(r'''
        \begin{equation*}
            \left\{\begin{aligned}
                N_{k+1} &= N_{k+1} = (1-m)R \displaystyle\frac{ N_k}{1 + c N_k}v_k\\
                v_{k+1} &= \displaystyle\frac{m v_k + \frac{1}{2}(1-v_k)}{m v_k + (1-v_k)}
        \end{aligned}\right.
        \end{equation*}
        ''')
        st.markdown("Where $N_k$ tracks the population dynamics while $v_k$ tracks the frequency of virulent nematodes (aa) and $m$ represents the relative proportion of juveniles that develop into virulent males to those that develop into avirulent males")
    st.markdown("### Basic reproduction number")
    checkbox = st.checkbox("Read the text")
    if checkbox:
        st.markdown("Several factors compose G. pallida's reproduction number $R$, which is the number of secondary infections generated by a single female, in a susceptible host population, at low nematode density, and in the absence of control.")
        st.markdown("The first one is the number of eggs, $e$, produced by a single female during her reproductive lifespan. The others are the viability of eggs inside the cyst, $1-\mu$, the fraction of eggs which survive accidental hatching in-between seasons, $1-h$, and the survival fraction of larvae in the soil, $s$:")
        st.latex(r'''
            \begin{equation}
                R = e(1-\mu)(1-h)s
            \end{equation}
        ''')
        st.markdown("Multiplying the reproduction number, $R$, by the proportion of female, $(1-m)$, yields the \textit{basic} reproduction number $R(1-m)$, that is the average number of daughters generated by a single mother, in a susceptible host population, and at low nematode density.")
        st.markdown("We model rotations as regular potato cultivation breaks of $r$ years, meaning that the potato is grown once every $r+1$ years. We consider that alternative crops to potato (or the absence of crop) have the same effect on nematodes (no trap-crop is used). We will refer to the $r$ as the ``rotation number''.")
        st.markdown("We model the biocontrol efficacy as the percentage of nematodes that do not survive biocontrol application, $b$ and we consider that biocontrol is applied every year, regardless of whether potato is grown or not.")
        st.markdown("Under these control methods, the reproduction number becomes")
        st.latex(r'''
        \begin{equation*} R' = e\big[(1-\mu)(1-h)(1-b)\big]^{r+1}s\,.
        \end{equation*}
        ''')
        
        
        markdown_text = r'''
        
        
        Blocking resistances quickly become obsolete as the resistance gene is fixed by natural selection. Thus, as with susceptible plants, the long-term suppression of PCNs must be ensured by a biocontrol that brings the nematode reprodution number $R'$ below 1/(1-m).

        On the other hand, masculinizing resistances keep a partial resistance indefinitely. This is conditioned by the male allocation rate $m$. When $m < \frac{1}{2}$, which is the case in real setups, it suffices to ensure the long-term suppression of PCNs
        that control efforts bring the reproduction number $R'$ below $2$. The partial resistance is ensured by the survival
        of susceptible phenotype through the pairing of avirulent males with virulent females.
        '''

        st.markdown(markdown_text)
    image2_path = "figs/scenario_diagram.png"
    st.image(image2_path, width=1000)
    st.markdown("### Parameters")
    table_md = r'''
    | Parameter | Description | Value | Range |
    | --- | --- | --- | --- |
    | $s$ | Survival fraction of larvae | $25\%$ | [0,100\%] |
    | $m$ | Male fraction in the progeny | $35\%$ | (0, 35\%] |
    | $e$ | Average number of eggs per cyst | 300 | [200, 500] |
    | $\mu$ | Yearly egg mortality | $10\%$ | [0.01, 20\%] |
    | $h$ | Yearly accidental hatching fraction | 17% | [0, 35\%] |
    | $b$ | Biocontrol efficacy fraction | variable | [0, 99.9\%] |
    | $c$ | Intraspecific competition parameter | 0.4 g | [0.1, 0.9] |
    | $\tau$ | Acceptance threshold | 1 egg/g of soil | [1, 3] eggs/g of soil |
    '''
    st.markdown(table_md)
    st.markdown("The simulation allows to choose the plant breed which is deployed per season. The other parameters, very little variable and estimated on literature data, are in the Settings menu.")
    
elif main_tab == "Simulation":
    st.markdown("# Simulation")
    # Other parameters
    col1, col2, col3 = st.columns([6, 6, 10])
    with col1:
        if st.button("Reset initial values"):
            st.session_state.v_freq = st.session_state.reset_v_freq
            st.session_state.init_infest = st.session_state.reset_init_infest
        st.markdown("### Initial values")
        subcol1, subcol2 = st.columns([1,1])
        with subcol1:
            st.session_state.v_freq = st.slider("Initial frequency of virulence $v_0$ (%):", min_value=0.0, max_value=100.0, value=st.session_state.v_freq*100, step=0.1)/100
        with subcol2:
            st.session_state.init_infest = st.number_input("Initial infestation $N_0$ (eggs/g of soil):", min_value=0, max_value=170, value=st.session_state.init_infest, step=1)
        
    with col2:
        if st.button("Reset set up"):
            st.session_state.num_years = st.session_state.reset_num_years
            st.session_state.detection_threshold = st.session_state.reset_detection_threshold
        st.markdown("### Simulation set up")
        subcol1, subcol2 = st.columns([1,1])
        with subcol1:
            st.session_state.num_years = st.number_input("Enter the number of years of simulation:", min_value=1, max_value=1000, value=st.session_state.num_years, step=1)
        with subcol2:
            st.session_state.detection_threshold = st.number_input(f"Acceptance threshold (eggs/g of soil):", min_value=1, max_value=3, value=st.session_state.detection_threshold, step=1)    
    with col3:
        if st.button("Reset control"):
            st.session_state.bc = st.session_state.reset_bc
            st.session_state.r = st.session_state.reset_r
        st.markdown("### Control set up")
        subcol1, subcol2 = st.columns([1,1])
        with subcol1:
            st.session_state.bc = st.slider("Biocontrol efficacy $b$ (%):", 0.0, 100.0, st.session_state.bc*100, 1.0)/100
        with subcol2:
            st.session_state.r = st.number_input("Rotation number", min_value=0, max_value=14, value=st.session_state.r, step=1)
    ###############################
    st.session_state.num_gen = int(np.ceil(st.session_state.num_years/(st.session_state.r + 1)))
    S = np.zeros(st.session_state.num_gen+1) #Susceptible
    N = np.zeros(st.session_state.num_gen+1)
    V = np.zeros(st.session_state.num_gen+1)
    V[0] = st.session_state.v_freq
    N[0] = st.session_state.init_infest
    S[0] = st.session_state.init_infest
    R = st.session_state.e*st.session_state.s*((1-st.session_state.mu)*(1-st.session_state.h)*(1-st.session_state.bc))**(st.session_state.r+1)
    for k in range(st.session_state.num_gen):
        V[k+1] = (st.session_state.m*V[k] + 0.5*(1-V[k])) / (st.session_state.m*V[k] + (1-V[k]))
        S[k+1] = (1-st.session_state.m)*R*S[k]/(1+st.session_state.c*S[k])
        N[k+1] = (1-st.session_state.m)*R*N[k]*V[k] / (1+st.session_state.c*N[k])
    generate_main_plot(N,S,V,R)
    #generate_min_pcb()
    col1, col2, col3 = st.columns([1, 8, 1])
    #with col2:
    #    link_text = "<a href='https://nemo-simulator.streamlit.app/' style='font-size: 20px;'>For more complex scenario simulation, click here to visit the Nemo Simulator</a>"
    #    st.markdown(link_text, unsafe_allow_html=True)
elif main_tab == "Edit Parameters":
    st.markdown("# Edit Parameters")
    st.markdown("These parameters describe the basic biology of the nematode. They are retrieved from intensive literature review and cautious estimations described in the main paper.")
    st.session_state.s = st.slider("Survival fraction of larvae $s$ (%):", min_value=0.0, max_value=100.0, value=st.session_state.s*100, step=0.1)/100
    st.session_state.m = st.slider("Average male fraction in the progeny $m$ (%):", min_value=0.0, max_value=40.0, value=st.session_state.m*100, step=0.1)/100
    st.session_state.mu = st.slider("Yearly egg mortality fraction $\mu$ (%):", min_value=0, max_value=20, value=int(st.session_state.mu*100), step=1)/100
    st.session_state.h = st.slider("Yearly accidental hatching fraction $h$ (%):", min_value=0, max_value=35, value=int(st.session_state.h*100), step=1)/100
    st.session_state.e = st.slider("Average eggs per cyst $e$:", min_value=200, max_value=500, value=st.session_state.e, step=1)
    st.session_state.c = st.slider("Intraspecific competition parameter $c$:", min_value=0.1, max_value=0.9, value=st.session_state.c, step=0.1)
