import streamlit as st

def main():
    st.title("Dynamic Sliders")

    # Get the number of seasons from the user
    num_seasons = st.number_input("Enter the number of seasons (n)", min_value=1, max_value=10, value=3, step=1)

    # Create a scrolling menu to select the season
    selected_season = st.selectbox("Select Season", range(1, num_seasons + 1))

    # Initialize session state to store slider values
    if 'slider_values' not in st.session_state:
        st.session_state.slider_values = {}

    # Create sliders for each season and show/hide based on the selected season
    for k in range(1, num_seasons + 1):
        season_name = f"Season {k}"

        if k == selected_season:
            # Use the stored value if available, otherwise initialize to 0.0
            progress = st.slider(f"{season_name} Progress", 0.0, 0.999, st.session_state.slider_values.get(k, 0.0), key=f"slider_{k}")
            st.session_state.slider_values[k] = progress  # Store the slider value in session state
        else:
            # If it's not the selected season, show the stored value without the slider
            progress = st.session_state.slider_values.get(k, 0.0)
        st.write(f"{season_name}: {progress:.1%}")

if __name__ == "__main__":
    main()
