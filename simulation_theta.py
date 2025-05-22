import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Field Parameters
alpha = 7.5 
Ks = 1.06e-4
n = 1.89
m = 1 - 1/n
Theta_s = 0.41
Theta_r = 0.065

# Domain 
L = 60
dx = 2.5
dt = 100
T = 1000*dt
Nt = int(T / dt)  # Number of time steps

x = np.linspace(0, L, int(L/dx) + 1)

# Theta to head
def head(Theta):
    h = np.full_like(Theta, -np.inf)  # initialize with -inf
    inside_bounds = (Theta > Theta_r) & (Theta < Theta_s)
    h[Theta >= Theta_s] = 0
    h[inside_bounds] = - (1/alpha) * (((Theta_s - Theta_r) / (Theta[inside_bounds] - Theta_r))**(1/m) - 1)**(1/n)
    h[Theta <= Theta_r] = -np.inf
    return h

def theta_h(h):
    h = np.array(h)
    Theta = np.full_like(h, np.nan)  # initialize with NaN
    inside_bounds = (h < 0)
    Theta[inside_bounds] = Theta_r + (Theta_s - Theta_r) * (1 + (alpha * np.abs(h[inside_bounds])) ** n) ** (-1 / m)
    Theta[h >= 0] = Theta_s
    return Theta


# initial conditions
#Theta0 = np.random.uniform(Theta_r, Theta_s, len(x))  # initial soil moisture
#Theta0 = Theta_r + (Theta_s - Theta_r) * (0.5 + 0.5 * np.sin(np.pi * x / L))
Theta0 = np.linspace(Theta_r + 0.05, Theta_s - 0.05, len(x))
h0 = head(Theta0)  # initial pressure head

# Saturation
def Se(h):
    h = np.array(h)
    return np.where(h < 0,
                    (1 + (alpha * np.abs(h)) ** n) ** (-m),
                    1.0)

#Capillary Capacity
def C(h):
    h = np.array(h)
    return np.where(h < 0,
                    (Theta_s - Theta_r) * alpha * n * m *
                    ((-alpha * h) ** (n - 1)) *
                    (1 + (-alpha * h) ** n) ** (-(2 - (1 / n))),
                    1.0e-04)

# Hydraulic Conductivity
def K(h):
    h = np.array(h)
    Se_val = Se(h)
    term = (1 - (1 - Se_val ** (1 / m)) ** m) ** 2
    return np.where(h < 0,
                    Ks * Se_val ** 0.5 * term,
                    Ks)

def Picard_iteration(h0,max_iter=10000, tol=1e-4):
    h = h0.copy()
    h_new = np.zeros_like(h)
    for it in range(max_iter):
        if it>0:
            h_old_n = h_old.copy()
        h_old = h.copy()
        if it == 0:
            h_old_n = h.copy()
        A = np.zeros([h.shape[0],h.shape[0]])
        b = np.zeros([h.shape[0]])

        A[-1,-2] = 0
        A[-1,-1] = 1
        b[-1] = 0
        A[0,0] = 1
        A[0,1] = 0
        b[0] = 0     
        C_i = C(h_old)
        for i in range(1,x.shape[0]-1):
            K_plus = (K(h_old[i])+K(h_old[i+1]))/2.0
            K_minus = (K(h_old[i])+K(h_old[i-1]))/2.0
            A[i,i-1] = -K_minus/dx**2
            A[i,i] = K_plus/dx**2 + K_minus/dx**2 +C_i[i]/dt
            A[i,i+1] = -K_plus/dx**2
            b[i] = (-C_i[i]/dt)*(h_old[i]-h_old_n[i]) + (K_plus/dx**2)*(h_old[i+1] - h_old[i]) - (K_minus/dx**2)*(h_old[i] - h_old[i-1]) + (K_plus - K_minus)/dx 
        
        delta_h = np.linalg.solve(A,b)
        h_new = h_old + delta_h
        h = h_new.copy()
        # Check convergence
        if np.linalg.norm(h - h_old, ord=np.inf) < tol:
            print(f"Converged in {it :.2f} iterations.")
            break
        else:
            print(f"Iteration {it :.2f}: max change = {np.linalg.norm(h - h_old, ord=np.inf)}")
        
    return h

h_vec = []
h_vec.append(h0)
theta_history = [theta_h(h0)]
for i in range(Nt):
    X = Picard_iteration(h0)
    h0 = X.copy()
    h_vec.append(X)
    theta_history.append(theta_h(X))


theta_history = np.array(theta_history)

# --- ANIMATION ---
fig, ax = plt.subplots(figsize=(8, 5))
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, L)
ax.set_ylim(Theta_r, Theta_s)
ax.set_xlabel("Position [m]")
ax.set_ylabel("Water Content θ")
ax.set_title(f"1D Richards Equation Solution (θ vs x)")
ax.grid(True)

def init():
    line.set_data([], [])
    return line,

def update(frame):
    line.set_data(x, theta_history[frame])
    ax.set_title(f"Time = {frame*dt:.2f} s")
    return line,

ani = FuncAnimation(fig, update, frames=len(theta_history),
                    init_func=init, blit=True, interval=50)

plt.tight_layout()
plt.show()