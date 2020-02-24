import numpy as np
a = np.load(file="protonated_clusters/oo_steps_dimer.npy")
print(a)


def get_reduced_mass(m1, m2):
    """
    calculate reduced mass in atomic units given two atoms
    :param m1: mass first atom in amu
    :type m1: float
    :param m2: mass second atom in amu
    :type m2: float
    :return: reduced mass (mu) in atomic units
    :rtype: float
    """
    au = 1822.89
    m1 *= au
    m2 *= au
    mu = (m1 * m2) / (m1 + m2)
    print(mu)
    return mu


"""common masses"""
m_hydrogen = 1.00784
m_oxygen = 15.999

# soln = get_reduced_mass(m_hydrogen, m_oxygen)

