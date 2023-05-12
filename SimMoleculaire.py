import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

# Abscisses des minimas
x1 = (-1-np.sqrt(5))/2
x2 = (-1+np.sqrt(5))/2

# Fonction potentiel
def V(x):
    return (3*x**4 + 4*x**3 - 6*x**2)/12

# Dérivée
def Vprime(x):
    return x**3 - x + x**2

def plot_potential():
    # représentation graphique
    x = np.linspace(-2.5, 1.5, 100)
    plt.plot(x, V(x), label='Potentiel')
    plt.title("Fonction Potentiel")
    # plot vertical line for x1 and x2
    plt.axvline(x1, color='r', linestyle='--', label=f'x1 = {round(x1, 3)}')
    plt.axvline(x2, color='r', linestyle='--', label=f'x2 = {round(x2, 3)}')
    plt.xlabel("Position")
    plt.ylabel("Potentiel")
    plt.legend()
    plt.show()

# Calcul de la proba de R={t2 < t1}, P(R) = E(1{t2-t1}) par méthode de Monte-Carlo naïve
def MC_naif(X0, epsilon, h, N):
    E = 0
    global x1
    global x2
    sigma = np.sqrt(2*epsilon*h)

    for _ in range(N):
        Xn = X0
        while Xn >= x1 and Xn <= x2:
            Xn = Xn - Vprime(Xn)*h + sigma*np.random.randn()
        if Xn > x2:
            E += 1

    return E/N

def simulation_trajectoire(X0, h, sigma):
    # simulation d'une trajectoire
    Xn = X0
    traj = [X0]
    while Xn >= x1 and Xn <= x2:
        Xn = Xn - Vprime(Xn)*h + sigma*np.random.randn()
        traj.append(Xn)
    return np.array(traj)

def K(trajectoires):
    # retourne le plus petit maximum des trajectoires
    return np.amin([np.amax(traj) for traj in trajectoires])

def I(trajectoires, K_value):
    # retourne les indices des trajectoires dont le maximum est K_value
    return [i for i, traj in enumerate(trajectoires) if np.amax(traj) == K_value]

def estimateur_AMS(X0, M, h, epsilon, debug=False):
    # variables
    global x2
    trajectoires = []
    I_q = []
    K_q = 0
    sigma = np.sqrt(2*epsilon*h)
    p_estim = 1

    # simulation initiale
    trajectoires = [simulation_trajectoire(X0, h, sigma) for _ in range(M)]
    K_q = K(trajectoires)
    I_q = I(trajectoires, K_q)
    p_estim *= 1 - len(I_q)/M

    q = 0
    while K_q < x2:
        if debug:
            print("Iteration n°{} \t K_q/x2 = {}".format(q, K_q/x2))

        # (i) Branchement
        indices = list(set(range(M)) - set(I_q))
        for m in I_q:
            p = np.random.choice(indices)
            idx = np.where(trajectoires[p] > K_q)[0][0]
            trajectoires[m] = np.concatenate((trajectoires[p][:idx], simulation_trajectoire(trajectoires[p][idx], h, sigma)))

        # (ii) Calcul du niveau
        K_q = K(trajectoires)
        I_q = I(trajectoires, K_q)
        p_estim *= 1 - len(I_q)/M

        q += 1

    return p_estim

def T12_MC_naif(x0, h, epsilon, N):
    # variables
    global x1
    global x2
    sigma = np.sqrt(2*epsilon*h)
    t1 = 0
    t2 = 0
    E = 0.0
    Xn = x0

    # simulation
    i = 0
    for _ in range(N):
        # calcul t1
        while Xn >= x1:
            Xn = Xn - Vprime(Xn)*h + sigma*np.random.randn()
            print(Xn)
            i += 1
        t1 = i
        while Xn <= x2:
            Xn = Xn - Vprime(Xn)*h + sigma*np.random.randn()
            print(Xn)
            i += 1
        t2 = i

        E += (t2 - t1)*h/N

    return E



"""
Partie nécessaire à l'execution du script comme script principal
"""

def plot_estimateurs():
    epsilons = [0.001, 0.01, 0.1, 0.5, 1]
    h = int(1e-3)

    naif = [MC_naif(x1+0.1, epsilon, h, 1000) for epsilon in epsilons]
    print("naif ok")
    AMS = [estimateur_AMS(x1+0.1, 1000, h, epsilon) for epsilon in epsilons]
    print("AMS ok")
    plt.plot(np.log10(epsilons), naif, label='Méthode de Monte-Carlo naïve')
    plt.plot(np.log10(epsilons), AMS, label='Estimateur AMS')
    plt.title("Comparaison des méthodes de Monte-Carlo")
    plt.xlabel("Epsilon (echelle log10)")
    plt.ylabel("Probabilité")
    plt.legend()
    # sauvegarde de l'image
    plt.savefig("comparaison.png")
    plt.show()


if __name__ == '__main__':
    print(estimateur_AMS(x1+0.1, 1000, 0.1, 0.01))