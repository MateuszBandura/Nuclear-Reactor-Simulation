import plotly.graph_objects as go
import numpy as np
from dash import Dash, html, dcc, callback, Output, Input, State
from scipy.integrate import solve_ivp

# Parameters (Set A - PWR-like)
Lambda = 1e-5  # s
lambda_i = np.array([0.0124, 0.0305, 0.1110, 0.3010, 1.14, 3.01])  # s^-1
beta_i = np.array([0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273])
beta = np.sum(beta_i)  # ≈ 0.0065
K_p = 3e9  # W (nominal power = 3 GW)
tau_th = 200.0  # s
T_in = 290.0  # °C
T_ref = 300.0  # °C - reference temperature for reactivity feedback
alpha_T = -5e-5  # K^-1

# Calculate H from steady-state: P0 = H*(T0 - T_in)
# At nominal power: 3e9 W = H * (300 - 290) K
# So: H = 3e9 / 10 = 3e8 W/K
P_nominal = K_p * 1.0  # Nominal power at n=1
H = P_nominal / (T_ref - T_in)  # W/K
C_th = tau_th * H  # J/K

# Simulation parameters
t_final = 500.0

# Initial conditions (steady state at nominal power)
n0 = 1.0
C0 = (beta_i / (Lambda * lambda_i)) * n0
T0 = T_ref
y0 = np.concatenate(([n0], C0, [T0]))

def reactor_ode(t, y, rho_ext_func):
    """
    ODE system for reactor dynamics
    y = [n, C1, C2, C3, C4, C5, C6, T]
    """
    n = y[0]
    C = y[1:7]
    T = y[7]
    
    # Get external reactivity at time t
    rho_ext = rho_ext_func(t)
    
    # Total reactivity with temperature feedback
    rho = rho_ext + alpha_T * (T - T_ref)
    
    # Neutron density equation (point kinetics)
    dn = ((rho - beta) / Lambda) * n + np.sum(lambda_i * C)
    
    # Delayed neutron precursor equations
    dC = (beta_i / Lambda) * n - lambda_i * C
    
    # Power
    P = K_p * n
    
    # Core temperature equation
    dT = (P - H * (T - T_in)) / C_th
    
    return np.concatenate(([dn], dC, [dT]))

def step_reactivity(t, rho_value, t_step=1.0):
    """Step reactivity insertion at t_step"""
    return rho_value if t >= t_step else 0.0

app = Dash(__name__)

app.layout = html.Div([
    html.H1('Symulacja Reaktora Jądrowego', style={'textAlign': 'center'}),
    html.Div([
        html.Label('Reaktywność zewnętrzna (Δk/k):'),
        dcc.Slider(
            id='slider-val',
            min=-0.005,
            max=0.005,
            step=0.0001,
            value=0.001,
            marks=None,
            tooltip={"placement": "bottom", "always_visible": True}
        ),
    ], style={'padding': '20px'}),
    html.Button('Uruchom Symulację', id='submit-val', n_clicks=0, 
                style={'margin': '20px', 'padding': '10px 20px', 'fontSize': '16px'}),
    html.Div(id='status-div', style={'textAlign': 'center', 'color': 'blue', 'padding': '10px'}),
    dcc.Graph(id="power-graph"),
    dcc.Graph(id="temperature-graph"),
    dcc.Graph(id="precursors-graph")
])

@callback(
    [Output('power-graph', 'figure'),
     Output('temperature-graph', 'figure'),
     Output('precursors-graph', 'figure'),
     Output('status-div', 'children')],
    Input('submit-val', 'n_clicks'),
    State('slider-val', 'value'),
    prevent_initial_call=True
)
def update_graphs(n_clicks, rho_value):
    """Update all three graphs"""
    
    # Verify steady state
    test_derivs = reactor_ode(0, y0, lambda t: 0.0)
    print(f"\n=== DIAGNOSTICS ===")
    print(f"Initial derivatives (should be ~0 at steady state):")
    print(f"  dn/dt = {test_derivs[0]:.6e}")
    print(f"  dC/dt = {test_derivs[1:7]}")
    print(f"  dT/dt = {test_derivs[7]:.6e}")
    print(f"Initial n = {y0[0]}, C sum = {np.sum(y0[1:7]):.6e}, T = {y0[7]}")
    print(f"beta/Lambda * n0 = {beta/Lambda * n0:.6e}")
    print(f"Sum(lambda_i * C0_i) = {np.sum(lambda_i * C0):.6e}")
    print(f"==================\n")
    
    # Define reactivity function (step at t=1s)
    rho_ext_func = lambda t: step_reactivity(t, rho_value, t_step=1.0)
    
    # Solve ODE
    sol = solve_ivp(
        lambda t, y: reactor_ode(t, y, rho_ext_func),
        [0, t_final],
        y0,
        method='Radau',  # Good for stiff systems
        rtol=1e-6,
        atol=1e-9
    )
    
    if not sol.success:
        error_msg = f"Błąd symulacji: {sol.message}"
        empty_fig = go.Figure()
        return empty_fig, empty_fig, empty_fig, error_msg
    
    t = sol.t
    n = sol.y[0]
    C = sol.y[1:7]
    T = sol.y[7]
    P = K_p * n / 1e9  # Power in GW
    
    # Power graph
    fig_power = go.Figure()
    fig_power.add_trace(go.Scatter(
        x=t, y=P,
        mode='lines',
        line=dict(color='#ff7f0e', width=2),
        name='Moc reaktora'
    ))
    fig_power.add_hline(y=K_p * n0 / 1e9, line_dash="dash", line_color="gray",
                        annotation_text="P₀ (nominalna)")
    fig_power.update_layout(
        title="Moc reaktora",
        xaxis_title="Czas [s]",
        yaxis_title="P(t) [GW]",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#F9F9F9",
        hovermode='x unified',
        showlegend=True
    )
    
    # Temperature graph
    fig_temp = go.Figure()
    fig_temp.add_trace(go.Scatter(
        x=t, y=T,
        mode='lines',
        line=dict(color='#d62728', width=2),
        name='Temperatura rdzenia'
    ))
    fig_temp.add_hline(y=T_ref, line_dash="dash", line_color="gray",
                       annotation_text="T₀")
    fig_temp.add_hline(y=T_in, line_dash="dash", line_color="blue",
                       annotation_text="T_in")
    fig_temp.update_layout(
        title="Temperatura rdzenia",
        xaxis_title="Czas [s]",
        yaxis_title="T(t) [°C]",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#F9F9F9",
        hovermode='x unified',
        showlegend=True
    )
    
    # Delayed neutron precursors graph
    fig_precursors = go.Figure()
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for i in range(6):
        fig_precursors.add_trace(go.Scatter(
            x=t, y=C[i] / C0[i],  # Normalized to initial value
            mode='lines',
            line=dict(color=colors[i], width=2),
            name=f'Grupa {i+1} (λ={lambda_i[i]:.3f} s⁻¹)'
        ))
    fig_precursors.update_layout(
        title="Prekursory opóźnionych neutronów (znormalizowane)",
        xaxis_title="Czas [s]",
        yaxis_title="C_i(t) / C_i(0) [względne]",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#F9F9F9",
        hovermode='x unified',
        showlegend=True,
        legend=dict(x=1.05, y=1)
    )
    
    # Calculate final values
    P_final = P[-1]
    P_initial = P[0]
    T_final = T[-1]
    percent_change = ((P_final / P_initial) - 1) * 100
    
    status = (f"✓ Symulacja zakończona | ρ_ext = {rho_value:.4f} ({rho_value*1e5:.1f} pcm) | "
              f"ΔP = {percent_change:+.2f}% | T_final = {T_final:.2f}°C")
    
    return fig_power, fig_temp, fig_precursors, status

if __name__ == '__main__':
    app.run(debug=True)