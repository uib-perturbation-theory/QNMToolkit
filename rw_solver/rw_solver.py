import numpy as np

try:
    from scipy.special import lambertw
    _have_scipy = True
except Exception:
    _have_scipy = False
    import mpmath as mp


# ============================================================
#  TORTOISE COORDINATE
# ============================================================

def r_of_s(s, M=1.0):
    s = np.asarray(s, dtype=float)
    z = np.exp(s/(2.0*M) - 1.0)

    if _have_scipy:
        y = lambertw(z, k=0).real
    else:
        y = np.vectorize(lambda zz: float(mp.lambertw(zz).real))(z)

    return 2.0*M*(1.0 + y)


# ============================================================
#  POTENTIAL
# ============================================================

def SchwarzPot(r, l, M=1.0):
    return (1 - 2*M/r)*(l*(l+1)/r**2 - 6*M/r**3)


# ============================================================
#  WINDOW
# ============================================================

def W_window(r, M=1.0, rmin=2.01, rramp=2.2, rstart=500.0, rcut=1000.0):

    r = np.asarray(r, dtype=float)

    rminM, rrampM = rmin*M, rramp*M
    rstartM, rcutM = rstart*M, rcut*M

    wL = np.zeros_like(r)
    wR = np.zeros_like(r)

    midL = (r > rminM) & (r < rrampM)
    wL[r >= rrampM] = 1.0
    wL[midL] = np.sin(0.5*np.pi*(r[midL]-rminM)/(rrampM-rminM))**2

    wR[r <= rstartM] = 1.0
    midR = (r > rstartM) & (r < rcutM)
    wR[midR] = np.cos(0.5*np.pi*(r[midR]-rstartM)/(rcutM-rstartM))**2

    return wL*wR


# ============================================================
#  DERIVATIVES
# ============================================================

def d1_4th(arr, ds):

    arr = np.asarray(arr)
    d = np.empty_like(arr)

    d[2:-2] = (-arr[4:] + 8*arr[3:-1] - 8*arr[1:-3] + arr[:-4])/(12*ds)

    d[0]  = (-25*arr[0] + 48*arr[1] - 36*arr[2] + 16*arr[3] - 3*arr[4])/(12*ds)
    d[1]  = (-3*arr[0] -10*arr[1] + 18*arr[2] - 6*arr[3] + arr[4])/(12*ds)
    d[-2] = ( 3*arr[-1] +10*arr[-2] -18*arr[-3] + 6*arr[-4] - arr[-5])/(12*ds)
    d[-1] = (25*arr[-1] -48*arr[-2] +36*arr[-3] -16*arr[-4] + 3*arr[-5])/(12*ds)

    return d


def L(u, ds, V_pot):

    h2 = ds**2
    L_ = np.empty_like(u)

    L_[2:-2] = (-u[:-4] + 16*u[1:-3] - 30*u[2:-2]
                + 16*u[3:-1] - u[4:])/(12*h2)

    L_[0]  = (35*u[0] -104*u[1] +114*u[2] -56*u[3] +11*u[4])/(12*h2)
    L_[1]  = (10*u[0] -15*u[1] -4*u[2] +14*u[3] -6*u[4] + u[5])/(12*h2)
    L_[-1] = (35*u[-1] -104*u[-2] +114*u[-3] -56*u[-4] +11*u[-5])/(12*h2)
    L_[-2] = (10*u[-1] -15*u[-2] -4*u[-3] +14*u[-4] -6*u[-5] + u[-6])/(12*h2)

    return L_ - V_pot*u


# ============================================================
#  RK4
# ============================================================

def rk4_wave(u, v, dt, ds, V):

    def F(uu):
        return L(uu, ds, V)

    k1u, k1v = v, F(u)
    k2u, k2v = v + 0.5*dt*k1v, F(u + 0.5*dt*k1u)
    k3u, k3v = v + 0.5*dt*k2v, F(u + 0.5*dt*k2u)
    k4u, k4v = v + dt*k3v,     F(u + dt*k3u)

    u_new = u + dt*(k1u + 2*k2u + 2*k3u + k4u)/6
    v_new = v + dt*(k1v + 2*k2v + 2*k3v + k4v)/6

    return u_new, v_new


# ============================================================
#  EVOLUTION
# ============================================================

def simula_full(N, x_min, x_max, Psi_initial,
                TF=400.0, CFL=0.1, l=2, M=1.0):

    s = np.linspace(x_min, x_max, N)
    ds = (x_max-x_min)/(N-1)

    r_vals = r_of_s(s, M=M)
    V = W_window(r_vals)*SchwarzPot(r_vals, l=l)

    Psi = Psi_initial.copy()
    Pi  = d1_4th(Psi, ds)

    dt = CFL*ds
    steps = int(TF/dt)

    U = []
    T = []

    t = 0.0
    for _ in range(steps):
        U.append(Psi.copy())
        T.append(t)

        Psi, Pi = rk4_wave(Psi, Pi, dt, ds, V)
        t += dt

    return np.array(U), np.array(T), s


# ============================================================
#  GAUSSIAN
# ============================================================

def gaussiana(A, x_min, x_max, N, w0, s0, sigma):
    s = np.linspace(x_min, x_max, N)
    return A*np.exp(-(s-s0)**2/(2*sigma**2))*np.cos(w0*(s-s0))


# ============================================================
#  EXTRACTION
# ============================================================

def extract_observer_4th(U, s, target):

    waveform = []

    for u in U:

        idx = np.searchsorted(s, target)
        idx = max(2, min(idx, len(s)-2))
        ii  = [idx-2, idx-1, idx, idx+1]

        xs = s[ii]
        fs = u[ii]

        num = 0.0
        for j in range(4):
            num_j = 1.0
            den_j = 1.0
            for k in range(4):
                if k != j:
                    num_j *= (target - xs[k])
                    den_j *= (xs[j] - xs[k])
            num += fs[j]*num_j/den_j

        waveform.append(num)

    return np.array(waveform)
