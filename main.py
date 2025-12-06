import plotly.graph_objects as go
import numpy as np
from dash import Dash, html, dcc, callback, Output, Input, State
from scipy.integrate import solve_ivp

# Simulation parameters
t_final = 500.0

# Parameter options for dropdowns
REACTOR_TYPES = {
    'PWR (A)': {
        'Lambda': 1e-5,
        'beta_i': np.array([0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]),
        'lambda_i': np.array([0.0124, 0.0305, 0.1110, 0.3010, 1.14, 3.01])
    },
    'Badawczy (B)': {
        'Lambda': 5e-5,
        'beta_i': np.array([0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]),
        'lambda_i': np.array([0.0124, 0.0305, 0.1110, 0.3010, 1.14, 3.01])
    },
    'Szybki (C)': {
        'Lambda': 1e-6,
        'beta_i': np.array([0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]),
        'lambda_i': np.array([0.0124, 0.0305, 0.1110, 0.3010, 1.14, 3.01])
    }
}

POWER_OPTIONS = {
    '1 MW': 1e6,
    '100 MW': 1e8,
    '1 GW': 1e9,
    '3 GW': 3e9
}

TAU_TH_OPTIONS = {
    '50 s': 50.0,
    '200 s': 200.0,
    '500 s': 500.0,
    '1000 s': 1000.0
}

ALPHA_T_OPTIONS = {
    '-1×10⁻⁵ K⁻¹': -1e-5,
    '-5×10⁻⁵ K⁻¹': -5e-5,
    '-1×10⁻⁴ K⁻¹': -1e-4,
    '-2×10⁻⁴ K⁻¹': -2e-4
}

T_IN_OPTIONS = {
    '270°C': 270.0,
    '280°C': 280.0,
    '290°C': 290.0,
    '300°C': 300.0
}

DELTA_T_OPTIONS = {
    '5°C': 5.0,
    '10°C': 10.0,
    '15°C': 15.0,
    '20°C': 20.0
}

def get_original_initial_conditions(reactor_type, K_p, T_in, delta_T):
    """Calculate original initial conditions based on parameters"""
    params = REACTOR_TYPES[reactor_type]
    Lambda = params['Lambda']
    beta_i = params['beta_i']
    lambda_i = params['lambda_i']
    
    T_ref = T_in + delta_T
    n0 = 1.0
    C0 = (beta_i / (Lambda * lambda_i)) * n0
    T0 = T_ref
    
    return np.concatenate(([n0], C0, [T0])), T_ref

def reactor_ode(t, y, rho_ext_func, params):
    """
    ODE system for reactor dynamics
    y = [n, C1, C2, C3, C4, C5, C6, T]
    """
    n = y[0]
    C = y[1:7]
    T = y[7]
    
    # Prevent negative neutron density
    if n < 0:
        n = 0
    
    # Get external reactivity at time t
    rho_ext = rho_ext_func(t)
    
    # Unpack parameters
    Lambda = params['Lambda']
    beta_i = params['beta_i']
    lambda_i = params['lambda_i']
    beta = params['beta']
    K_p = params['K_p']
    alpha_T = params['alpha_T']
    T_ref = params['T_ref']
    T_in = params['T_in']
    H = params['H']
    C_th = params['C_th']
    
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

app = Dash(__name__)

app.layout = html.Div([
    html.H1('Symulacja Reaktora Jądrowego', style={'textAlign': 'center'}),
    
    # Store component to hold current state and last simulation
    dcc.Store(id='current-state', data={
        'y0': None,  # Will be initialized on first run
        'rho_current': 0.0,
        'is_original': True,
        'simulation_count': 0,
        'last_simulation': None,
        'parameters': {
            'reactor_type': 'PWR (A)',
            'power': '3 GW',
            'tau_th': '200 s',
            'alpha_t': '-5×10⁻⁵ K⁻¹',
            't_in': '290°C',
            'delta_t': '10°C'
        }
    }),
    
    # Parameter selection panel
    html.Div([
        html.H3('Parametry Reaktora', style={'textAlign': 'center', 'marginBottom': '20px'}),
        
        html.Div([
            html.Div([
                html.Label('Typ reaktora:', style={'fontWeight': 'bold'}),
                dcc.Dropdown(
                    id='reactor-type',
                    options=[{'label': k, 'value': k} for k in REACTOR_TYPES.keys()],
                    value='PWR (A)',
                    clearable=False,
                    style={'width': '200px'}
                )
            ], style={'display': 'inline-block', 'margin': '10px'}),
            
            html.Div([
                html.Label('Moc nominalna (Kₚ):', style={'fontWeight': 'bold'}),
                dcc.Dropdown(
                    id='power',
                    options=[{'label': k, 'value': k} for k in POWER_OPTIONS.keys()],
                    value='3 GW',
                    clearable=False,
                    style={'width': '150px'}
                )
            ], style={'display': 'inline-block', 'margin': '10px'}),
            
            html.Div([
                html.Label('Stała czasowa (τₜₕ):', style={'fontWeight': 'bold'}),
                dcc.Dropdown(
                    id='tau-th',
                    options=[{'label': k, 'value': k} for k in TAU_TH_OPTIONS.keys()],
                    value='200 s',
                    clearable=False,
                    style={'width': '150px'}
                )
            ], style={'display': 'inline-block', 'margin': '10px'}),
        ], style={'textAlign': 'center'}),
        
        html.Div([
            html.Div([
                html.Label('Współczynnik temp. (αₜ):', style={'fontWeight': 'bold'}),
                dcc.Dropdown(
                    id='alpha-t',
                    options=[{'label': k, 'value': k} for k in ALPHA_T_OPTIONS.keys()],
                    value='-5×10⁻⁵ K⁻¹',
                    clearable=False,
                    style={'width': '180px'}
                )
            ], style={'display': 'inline-block', 'margin': '10px'}),
            
            html.Div([
                html.Label('Temperatura wlotu (Tᵢₙ):', style={'fontWeight': 'bold'}),
                dcc.Dropdown(
                    id='t-in',
                    options=[{'label': k, 'value': k} for k in T_IN_OPTIONS.keys()],
                    value='290°C',
                    clearable=False,
                    style={'width': '150px'}
                )
            ], style={'display': 'inline-block', 'margin': '10px'}),
            
            html.Div([
                html.Label('Przyrost temperatury (ΔT):', style={'fontWeight': 'bold'}),
                dcc.Dropdown(
                    id='delta-t',
                    options=[{'label': k, 'value': k} for k in DELTA_T_OPTIONS.keys()],
                    value='10°C',
                    clearable=False,
                    style={'width': '150px'}
                )
            ], style={'display': 'inline-block', 'margin': '10px'}),
        ], style={'textAlign': 'center'}),
        
    ], style={'backgroundColor': '#f0f0f0', 'padding': '20px', 'marginBottom': '20px', 'borderRadius': '10px'}),
    
    html.Div([
        html.Label('Reaktywność zewnętrzna (Δk/k):'),
        dcc.Slider(
            id='slider-val',
            min=-0.005,
            max=0.005,
            step=0.0001,
            value=0.001,
            marks={i/1000: f'{i} pcm' for i in range(-5, 6)},
            tooltip={"placement": "bottom", "always_visible": True}
        ),
    ], style={'padding': '20px'}),
    
    html.Div([
        html.Button('Uruchom Symulację', id='submit-val', n_clicks=0, 
                    style={'margin': '10px', 'padding': '10px 20px', 'fontSize': '16px',
                           'backgroundColor': '#4CAF50', 'color': 'white', 'border': 'none',
                           'borderRadius': '4px', 'cursor': 'pointer'}),
        html.Button('Reset do Wartości Początkowych', id='reset-val', n_clicks=0,
                    style={'margin': '10px', 'padding': '10px 20px', 'fontSize': '16px',
                           'backgroundColor': '#f44336', 'color': 'white', 'border': 'none',
                           'borderRadius': '4px', 'cursor': 'pointer'}),
    ], style={'textAlign': 'center'}),
    
    html.Div(id='status-div', style={'textAlign': 'center', 'color': 'blue', 'padding': '10px'}),
    html.Div(id='state-info', style={'textAlign': 'center', 'color': 'green', 'padding': '5px', 
                                      'fontSize': '14px'}),
    
    dcc.Graph(id="power-graph"),
    dcc.Graph(id="temperature-graph"),
    dcc.Graph(id="precursors-graph")
])

@callback(
    Output('current-state', 'data'),
    [Input('reactor-type', 'value'),
     Input('power', 'value'),
     Input('tau-th', 'value'),
     Input('alpha-t', 'value'),
     Input('t-in', 'value'),
     Input('delta-t', 'value')],
    State('current-state', 'data'),
    prevent_initial_call=True
)
def reset_on_parameter_change(reactor_type, power, tau_th, alpha_t, t_in, delta_t, state_data):
    """Automatically reset state when any parameter changes"""
    
    # Get current parameters from state
    stored_params = state_data.get('parameters', {})
    
    # Check if any parameter changed
    current_params = {
        'reactor_type': reactor_type,
        'power': power,
        'tau_th': tau_th,
        'alpha_t': alpha_t,
        't_in': t_in,
        'delta_t': delta_t
    }
    
    # If parameters changed, reset to new initial conditions
    if current_params != stored_params:
        K_p = POWER_OPTIONS[power]
        T_in_val = T_IN_OPTIONS[t_in]
        delta_T = DELTA_T_OPTIONS[delta_t]
        
        y0_original, T_ref = get_original_initial_conditions(reactor_type, K_p, T_in_val, delta_T)
        
        return {
            'y0': y0_original.tolist(),
            'rho_current': 0.0,
            'is_original': True,
            'simulation_count': 0,
            'last_simulation': state_data.get('last_simulation'),  # Keep for comparison
            'parameters': current_params
        }
    
    # No change, return unchanged
    return state_data

@callback(
    Output('current-state', 'data', allow_duplicate=True),
    Input('reset-val', 'n_clicks'),
    [State('reactor-type', 'value'),
     State('power', 'value'),
     State('t-in', 'value'),
     State('delta-t', 'value'),
     State('current-state', 'data')],
    prevent_initial_call=True
)
def reset_state(n_clicks, reactor_type, power, t_in, delta_t, state_data):
    """Reset to original initial conditions and clear last simulation"""
    K_p = POWER_OPTIONS[power]
    T_in = T_IN_OPTIONS[t_in]
    delta_T = DELTA_T_OPTIONS[delta_t]
    
    y0_original, T_ref = get_original_initial_conditions(reactor_type, K_p, T_in, delta_T)
    
    return {
        'y0': y0_original.tolist(),
        'rho_current': 0.0,
        'is_original': True,
        'simulation_count': 0,
        'last_simulation': None,  # Clear last simulation on manual reset
        'parameters': state_data.get('parameters', {})
    }

@callback(
    [Output('power-graph', 'figure'),
     Output('temperature-graph', 'figure'),
     Output('precursors-graph', 'figure'),
     Output('status-div', 'children'),
     Output('state-info', 'children'),
     Output('current-state', 'data', allow_duplicate=True)],
    Input('submit-val', 'n_clicks'),
    [State('slider-val', 'value'),
     State('current-state', 'data'),
     State('reactor-type', 'value'),
     State('power', 'value'),
     State('tau-th', 'value'),
     State('alpha-t', 'value'),
     State('t-in', 'value'),
     State('delta-t', 'value')],
    prevent_initial_call=True
)
def update_graphs(n_clicks, rho_value, state_data, reactor_type, power, tau_th, alpha_t, t_in, delta_t):
    """Update all three graphs and save final state"""
    
    # Get parameter values from dropdowns
    reactor_params = REACTOR_TYPES[reactor_type]
    Lambda = reactor_params['Lambda']
    beta_i = reactor_params['beta_i']
    lambda_i = reactor_params['lambda_i']
    beta = np.sum(beta_i)
    
    K_p = POWER_OPTIONS[power]
    tau_th_val = TAU_TH_OPTIONS[tau_th]
    alpha_T = ALPHA_T_OPTIONS[alpha_t]
    T_in = T_IN_OPTIONS[t_in]
    delta_T = DELTA_T_OPTIONS[delta_t]
    T_ref = T_in + delta_T
    
    # Calculate H and C_th from steady state
    P_nominal = K_p * 1.0
    H = P_nominal / delta_T
    C_th = tau_th_val * H
    
    # Package parameters
    params = {
        'Lambda': Lambda,
        'beta_i': beta_i,
        'lambda_i': lambda_i,
        'beta': beta,
        'K_p': K_p,
        'alpha_T': alpha_T,
        'T_ref': T_ref,
        'T_in': T_in,
        'H': H,
        'C_th': C_th
    }
    
    # Get current initial conditions from stored state
    # If first run, initialize with original conditions
    if state_data['y0'] is None:
        y0_original, _ = get_original_initial_conditions(reactor_type, K_p, T_in, delta_T)
        y0_current = y0_original
    else:
        y0_current = np.array(state_data['y0'])
    
    rho_previous = state_data.get('rho_current', 0.0)
    sim_count = state_data.get('simulation_count', 0)
    last_sim = state_data.get('last_simulation', None)
    
    # Define reactivity function
    def rho_ext_func(t):
        if t < 1.0:
            return rho_previous
        else:
            return rho_value
    
    # Solve ODE starting from current state
    try:
        sol = solve_ivp(
            lambda t, y: reactor_ode(t, y, rho_ext_func, params),
            [0, t_final],
            y0_current,
            method='Radau',
            rtol=1e-6,
            atol=1e-9
        )
        
        if not sol.success:
            raise RuntimeError(f"Solver failed: {sol.message}")
            
    except Exception as e:
        error_msg = f"Błąd symulacji: {str(e)}"
        empty_fig = go.Figure()
        empty_fig.update_layout(title=error_msg)
        state_info = "Symulacja nie powiodła się - spróbuj innych parametrów"
        return empty_fig, empty_fig, empty_fig, error_msg, state_info, state_data
    
    t = sol.t
    n = sol.y[0]
    C = sol.y[1:7]
    T = sol.y[7]
    P = K_p * n / 1e9  # Power in GW
    
    # Get FINAL state (end of simulation) for next run
    y0_next = sol.y[:, -1]
    
    # Ensure no negative neutron density
    if y0_next[0] < 0:
        y0_next[0] = 0
    
    # Get original C0 for normalization
    y0_original, _ = get_original_initial_conditions(reactor_type, K_p, T_in, delta_T)
    C0_original = y0_original[1:7]
    
    # Power graph with previous simulation overlay
    fig_power = go.Figure()
    
    # Plot previous simulation if it exists (faded)
    if last_sim is not None:
        fig_power.add_trace(go.Scatter(
            x=last_sim['t'], 
            y=last_sim['P'],
            mode='lines',
            line=dict(color='gray', width=2, dash='dash'),
            opacity=0.5,
            name='Poprzednia symulacja'
        ))
    
    # Plot current simulation (bright)
    fig_power.add_trace(go.Scatter(
        x=t, y=P,
        mode='lines',
        line=dict(color='#ff7f0e', width=2.5),
        name='Bieżąca symulacja'
    ))
    
    fig_power.add_hline(y=K_p * 1.0 / 1e9, line_dash="dash", line_color="lightgray",
                        annotation_text="P₀ (nominalna)")
    fig_power.update_layout(
        title="Moc reaktora (porównanie z poprzednią symulacją)",
        xaxis_title="Czas [s]",
        yaxis_title="P(t) [GW]",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#F9F9F9",
        hovermode='x unified',
        showlegend=True
    )
    
    # Temperature graph with previous simulation overlay
    fig_temp = go.Figure()
    
    # Plot previous simulation if it exists (faded)
    if last_sim is not None:
        fig_temp.add_trace(go.Scatter(
            x=last_sim['t'],
            y=last_sim['T'],
            mode='lines',
            line=dict(color='gray', width=2, dash='dash'),
            opacity=0.5,
            name='Poprzednia symulacja'
        ))
    
    # Plot current simulation (bright)
    fig_temp.add_trace(go.Scatter(
        x=t, y=T,
        mode='lines',
        line=dict(color='#d62728', width=2.5),
        name='Bieżąca symulacja'
    ))
    
    fig_temp.add_hline(y=T_ref, line_dash="dash", line_color="lightgray",
                       annotation_text="T₀")
    fig_temp.add_hline(y=T_in, line_dash="dash", line_color="lightblue",
                       annotation_text="T_in")
    fig_temp.update_layout(
        title="Temperatura rdzenia (porównanie z poprzednią symulacją)",
        xaxis_title="Czas [s]",
        yaxis_title="T(t) [°C]",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#F9F9F9",
        hovermode='x unified',
        showlegend=True
    )
    
    # Delayed neutron precursors graph (current simulation only)
    fig_precursors = go.Figure()
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for i in range(6):
        fig_precursors.add_trace(go.Scatter(
            x=t, y=C[i] / C0_original[i],
            mode='lines',
            line=dict(color=colors[i], width=2),
            name=f'Grupa {i+1} (λ={lambda_i[i]:.3f} s⁻¹)'
        ))
    fig_precursors.update_layout(
        title=f"Prekursory opóźnionych neutronów - Symulacja {sim_count+1}",
        xaxis_title="Czas [s] (względny)",
        yaxis_title="C_i(t) / C_i(0_original) [względne]",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#F9F9F9",
        hovermode='x unified',
        showlegend=True,
        legend=dict(x=1.05, y=1)
    )
    
    # Calculate statistics
    P_initial = P[0]
    P_final = P[-1]
    T_initial = T[0]
    T_final = T[-1]
    n_final = n[-1]
    percent_change = ((P_final / P_initial) - 1) * 100 if P_initial > 0 else 0
    rho_change = rho_value - rho_previous
    
    status = (f"✓ Symulacja {sim_count+1} zakończona | "
              f"ρ_poprzednia = {rho_previous:.4f} → ρ_nowa = {rho_value:.4f} "
              f"(Δρ = {rho_change:+.4f} = {rho_change*1e5:+.1f} pcm) | "
              f"ΔP = {percent_change:+.2f}% | T_final = {T_final:.2f}°C")
    
    # State info
    state_info = (f"Warunki końcowe tej symulacji → Warunki początkowe następnej: "
                  f"n = {n_final:.4f}, T = {T_final:.2f}°C, ρ = {rho_value:.4f}")
    
    # Save current simulation data for next run (to show as "previous")
    current_sim_data = {
        't': t.tolist(),
        'P': P.tolist(),
        'T': T.tolist()
    }
    
    # Update stored state with final values
    new_state = {
        'y0': y0_next.tolist(),
        'rho_current': rho_value,
        'is_original': False,
        'simulation_count': sim_count + 1,
        'last_simulation': current_sim_data,
        'parameters': {
            'reactor_type': reactor_type,
            'power': power,
            'tau_th': tau_th,
            'alpha_t': alpha_t,
            't_in': t_in,
            'delta_t': delta_t
        }
    }
    
    return fig_power, fig_temp, fig_precursors, status, state_info, new_state

if __name__ == '__main__':
    app.run(debug=True)