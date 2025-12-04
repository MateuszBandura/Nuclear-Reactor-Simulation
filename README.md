# Równania do symulacji (dziedzina czasu)

## 1. Neutrony (punktowa kinetyka, 6 grup)

\[
\frac{d n(t)}{d t}
= \frac{\rho(t)-\beta}{\Lambda}\,n(t)
+ \sum_{i=1}^{6}\lambda_i\,C_i(t)
\]

## 2. Prekursory opóźnione \((i=1..6)\)

\[
\frac{d C_i(t)}{d t}
= \frac{\beta_i}{\Lambda}\,n(t)
- \lambda_i\,C_i(t),
\qquad i=1,\dots,6
\]

## 3. Temperatura rdzenia (1-węzeł)

\[
C_{th}\,\frac{d T(t)}{d t}
= P(t) - H\bigl(T(t)-T_{in}\bigr),
\qquad
P(t)=K_p\,n(t)
\]

## 4. Sprzężenie reaktywności

\[
\rho(t) = \rho_{\text{ext}}(t) + \alpha_T\bigl(T(t)-T_0\bigr)
\]

## 5. Warunek odwzorowania mocy

\[
P(t)=K_p\,n(t)
\]

---

# Parametry (symbole i jednostki)

- \(n(t)\) — gęstość neutronów (względna)
- \(C_i(t)\) — koncentracja prekursorów
- \(T(t)\) — temperatura rdzenia
- \(\rho_{\text{ext}}(t)\), \(\rho(t)\) — reaktywności
- \(\alpha_T\) — temperaturowy współczynnik reaktywności
- \(\beta_i\), \(\beta\) — frakcje opóźnionych neutronów
- \(\lambda_i\) — stałe rozpadu
- \(\Lambda\) — czas generacji neutronów
- \(K_p\) — przelicznik gęstości neutronów na moc
- \(C_{th}\), \(H\), \(\tau_{th}\) — parametry cieplne
- \(T_{in}\), \(T_0\)

---

# Przykładowe zestawy parametrów

## Zestaw A — PWR-like

\[
\Lambda = 1\times10^{-5}\,\mathrm{s}
\]

\[
\beta_i = [0.000215,\; 0.001424,\; 0.001274,\; 0.002568,\; 0.000748,\; 0.000273]
\]

\[
\lambda_i = [0.0124,\; 0.0305,\; 0.1110,\; 0.3010,\; 1.14,\; 3.01]\ \mathrm{s^{-1}}
\]

\[
K_p = 3\times10^{9}\,\mathrm{W},
\qquad
\tau_{th}=200\,\mathrm{s},
\qquad
\alpha_T=-5\times10^{-5}\,\mathrm{K^{-1}}
\]

## Zestaw B — reaktor badawczy

\[
\Lambda = 5\times10^{-5}\,\mathrm{s},\quad
K_p = 1\times10^{6}\,\mathrm{W},\quad
\tau_{th}=50\,\mathrm{s},\quad
\alpha_T=-1\times10^{-5}\,\mathrm{K^{-1}}
\]

## Zestaw C — szybki/eksperymentalny

\[
\Lambda = 1\times10^{-6}\,\mathrm{s},\quad
K_p = 1\times10^{8}\,\mathrm{W},\quad
\tau_{th}=1000\,\mathrm{s},\quad
\alpha_T=-2\times10^{-4}\,\mathrm{K^{-1}}
\]

---

# Przykład sygnału sterującego

Skok reaktywności \(100\ \mathrm{pcm}\):

\[
\rho_{\text{ext}}(t) =
\begin{cases}
0, & t<1 \\[4pt]
0.001, & t\ge 1
\end{cases}
\]

---

# Warunki początkowe

\[
n(0)=1
\]

\[
C_i(0)=\frac{\beta_i}{\Lambda\,\lambda_i}
\]

\[
T(0)=T_0
\]

---

# Transmitancja Laplace’a

\[
G_{n\rho}(s)
= \frac{N(s)}{\mathcal{R}(s)}
= \frac{\dfrac{n_0}{\Lambda}}
{s + \dfrac{\beta}{\Lambda}
- \dfrac{1}{\Lambda}\sum_{i=1}^{6}
 \frac{\lambda_i\beta_i}{s+\lambda_i}}
\]

\[
G_{P\rho}(s)=K_p\,G_{n\rho}(s),
\qquad
G_{T\rho}(s)=\frac{K_p/C_{th}}{s+1/\tau_{th}}\,G_{n\rho}(s)
\]
