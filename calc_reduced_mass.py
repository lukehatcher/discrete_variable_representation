def get_reduced_mass(m1, m2):
    """
    given two atoms, calculate their reduced mass in atomic units
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


"""commonly used masses"""
m_hydrogen = 1.00784
m_oxygen = 15.999
m_carbon = 12.0107
m_nitrogen = 14.0067

#soln = get_reduced_mass(m_oxygen, m_oxygen)

