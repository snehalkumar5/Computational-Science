from autograd import grad
import autograd.numpy as np
import autograd


def sep(p1, p2):
    return p2 - p1

def pbc_sep(p1, p2):

    arrow = sep(p1, p2)
    rem = np.mod(arrow, 18)  
    mic_separation_vector = np.mod(rem+18/2, 18)-18/2

    return np.array(mic_separation_vector)

def get_energy(geom):
    faa = 0.238
    check = 3.4
    te = 0
    nexxt = np.array(geom)
    pairs = [(a, b) for idx, a in enumerate(nexxt)
             for b in nexxt[idx + 1:]]
    for p in pairs:
        rij = np.linalg.norm(pbc_sep(p[0], p[1]))
        if rij == 0:
            continue
        te += 4*faa*((check/rij)**12-(check/rij)**6)

    return te


def gradient_descent(factor, iters, w):
    gradient = grad(get_energy)

    weghts = [w]           
    cst = [get_energy(w)]
    for k in range(iters):
        # print("It:", k)
        grad_eval = gradient(w)

        w = w - factor*grad_eval
        print(get_energy(w))
        weghts.append(w)
        cst.append(get_energy(w))
    return weghts, cst