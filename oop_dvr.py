import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import argparse


class ParseGaussian:
    """
    parse txt file copy of Summary of Optimized Potential Surface Scans from gaussian

    this class is awful
    """

    def __init__(self, filename, num_oo_steps, num_oh_steps):
        """

        :param filename: name of txt copy of gauss summary
        :type filename: str
        :param num_oo_steps: self explanatory
        :type num_oo_steps: int
        :param num_oh_steps: self explanatory
        :type num_oh_steps: int
        """

        self.filename = filename
        self.num_oo_steps = num_oo_steps
        self.num_oh_steps = num_oh_steps

    def do_parsing(self):
        with open(self.filename, "r") as gauss_data:
            ar = [line.split() for line in gauss_data]
            e_list = []
            r4oh_list = []
            r2oo_list = []
            for i in ar:
                if i[0] == "Eigenvalues":
                    energy = i[2:]
                    z = energy[0]  # only one value
                    energy2 = z.split('-')
                    del energy2[0]  # removes weird first term from bad delimiter
                    # energy2 *= -1
                    e_list.append(energy2)
                elif i[0] == "r4":
                    r4oh_list.append(i[1:])
                elif i[0] == "r2":
                    r2oo_list.append(i[1:])

            """ flatten arrays """
            flat_elist = [x for list in e_list for x in list]
            flat_r4oh_list = [x for list in r4oh_list for x in list]
            flat_r2oo_list = [x for list in r2oo_list for x in list]

            """ sort flat arrays """
            oo_np = np.array([float(x) for x in flat_r2oo_list])
            oh_np = np.array([float(x) for x in flat_r4oh_list])
            eng_np = np.array([float(x) for x in flat_elist])
            eng_np *= -1.0  # fixes removed "delimiter"
            sorted_oo_pos = np.lexsort((oh_np, oo_np))

            new_oo = oo_np[sorted_oo_pos]
            new_oh = oh_np[sorted_oo_pos]
            new_eng = eng_np[sorted_oo_pos]

            """ separate by oo steps and save """
            a = self.num_oh_steps - self.num_oh_steps
            b = self.num_oh_steps
            np.save(arr=new_oh[0:self.num_oh_steps], file="dimer_r4oh3")  # EDIT EACH TRIAL!

            real_oo = np.zeros(self.num_oo_steps)  # real_oo includes one of each OO position
            for i in range(1, self.num_oo_steps + 1):
                np.save(arr=new_eng[a:b], file="dimer_Es3_" + str(i))  # EDIT EACH TRIAL!
                real_oo[i - 1] = new_oo[(i - 1) * self.num_oh_steps]
                a += self.num_oh_steps
                b += self.num_oh_steps
            np.save(arr=real_oo, file="oo_steps_dimer3")  # EDIT EACH TRIAL!

        return None


class DVR:
    def __init__(self, dvrgrid, potentials, mu):
        """
        discrete variable representation
        :param dvrgrid: aka OH steps, UNITS=ATOMIC
        :type dvrgrid: np array
        :param potentials: Energies for a
        :type potentials: np array
        :param mu: reduced mass, calc by hand prior UNITS=ATOMIC, thus *1822.89 for amu --> au
        :type mu: float
        """

        self.dvrgrid = dvrgrid
        self.potentials = potentials
        self.mu = mu
        self.deltax = self.dvrgrid[1] - self.dvrgrid[0]

    def pot_energy(self):
        v_matrix = np.diag(self.potentials)  # whats the units on the loaded self.pot?
        return v_matrix

    def kinetic_energy(self):
        t_matrix = np.zeros((len(self.dvrgrid), len(self.dvrgrid)))
        for i in range(len(self.dvrgrid)):
            for j in range(len(self.dvrgrid)):
                if i == j:
                    t_matrix[i, j] = (np.pi ** 2) / (6 * (self.deltax ** 2) * self.mu)
                else:
                    t_matrix[i, j] = ((-1) ** (i - j)) / (self.mu * ((i - j) ** 2) * (self.deltax**2))
        return t_matrix

    def run_dvr(self):
        v_matrix = self.pot_energy()
        kinetic_energy = self.kinetic_energy()
        evals, wfns_vecs = np.linalg.eigh(v_matrix + kinetic_energy)
        # np.save(file="dimer_dvrwfns_" + str(param.dvr_numb), arr=wfns_vecs)  ##############

        """ now plot """
        # plt.plot((self.dvrgrid * 0.529177), wfns_vecs[:, 0] ** 2)  # plotting in angstrom
        # plt.title("$\\Psi^2$")
        # plt.show()

        return wfns_vecs


class Interpolate1D:
    def __init__(self, x, y, n_points):
        """
        interpolate 1D array to desired size
        :param x: x axis array
        :type x: numpy array
        :param y: y axis array
        :type y: numpy array
        :param n_points: number of points to interpolate to
        :type n_points: int
        """

        self.x = x
        self.y = y
        self.n_points = n_points

    def get_interp(self):
        f = interpolate.interp1d(self.x, self.y, kind="cubic")
        new_xOH = np.linspace(start=np.amin(self.x), stop=np.amax(self.x), num=self.n_points, endpoint=True)
        new_yE = f(new_xOH)
        return new_xOH, new_yE


def use_dimer_files(numb_files, filenameEs, fname_save_interpE, fname_savewfns):
    """
    get dvr wfns at all n of OO distances
    :param numb_files: aka numb of OO distances
    :type numb_files: int
    :param filenameEs: name of energy files EXCLUDE EXTENSION
    :type filenameEs: str
    :param fname_saveE: name wanted for saving interped Es
    :type fname_saveE: str
    :param fname_savewfns: name wanted for saving wfns
    :type fname_savewfns: str
    :return: wfn_data:
    :rtype: numpy array
    """
    for i in range(1, numb_files + 1):
        # parse_ob = ParseGaussian()  #  need to add this

        engy_file = np.load(file=filenameEs + str(i) + ".npy")
        interp_ob = Interpolate1D(grid_arr, engy_file, 2000)  # grid_arr is hard coded
        new_xOH, new_yE = interp_ob.get_interp()
        np.save(file=fname_save_interpE + str(i), arr=new_yE)
        np.save(file="xOH3", arr=new_xOH)  # saved in bohr

        dvr_ob = DVR(new_xOH, new_yE, 1728.3085005881399)
        wfn_data = dvr_ob.run_dvr()
        np.save(file=fname_savewfns + str(i), arr=wfn_data)

    return None


# parse_inst = ParseGaussian("dimer_gaussian_data_run2/gauss_dimer_output2_txt", 18, 16)
# parse_inst = ParseGaussian("dimer_gaussian_data_run3/gauss_dimer_output3_txt", 18, 16)
# parse_inst.do_parsing()
#creates files to be used in call below

""" standard OH grid for passing """
grid_arr = np.load(file="dimer_r4oh3.npy")  # EDIT EACH TRIAL
grid_arr *= 1.88973  # angst -> bohr: gaussian gives angst, dvr needs bohr


run = use_dimer_files(18, "dimer_Es3_", "interpd_dimer_E3_", "dimer_dvrwfns3_")


# if __name__ == '__main__':
