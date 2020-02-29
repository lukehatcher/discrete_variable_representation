import numpy as np
import matplotlib.pyplot as plt


class DVRAnalysis:
    def __init__(self, nwfns,
                 wfn_fname,
                 xOH_fname,
                 oo_fname
                 ):
        """
        calc expectation values <O^hat> and <O^Hat**2> for position operator given wfns
        :param nwfns: number of wave functions
        :type nwfns: int
        :param wfn_fname: filename of wfn data to use EXCLUDE EXTENSION
        :type wfn_fname: str
        :param xOH_fname: filename of OH steps INCLUDE EXT
        :type xOH_fname: str
        :param oo_fname: filename of OO steps
        :type oo_fname : str
        """
        self.nwfns = nwfns
        self.wfn_fname = wfn_fname
        self.xOH_fname = xOH_fname
        self.oo_fname = oo_fname
        self.myexp_val_list = np.zeros(self.nwfns)
        self.myexp_val_list2 = np.zeros(self.nwfns)
        self.oo_vals = np.load(self.oo_fname)
        self.x = np.load(self.xOH_fname) / 1.88973  # bohr -> angst

    def calc_expvals(self):
        for i in range(1, self.nwfns + 1):
            wfns = np.load(self.wfn_fname + str(i) + ".npy")  # "dimer_dvrwfns_"
            psi = wfns[:, 0]
            exp_val = psi @ (self.x * psi)
            exp_val2 = psi @ ((self.x ** 2) * psi)
            self.myexp_val_list[i - 1] = np.copy(exp_val)
            self.myexp_val_list2[i - 1] = np.copy(exp_val2)
        return self.myexp_val_list, self.myexp_val_list2

    def calc_stand_dev(self):
        ct = 0
        for expval, expval2 in zip(self.myexp_val_list, self.myexp_val_list2):
            sigma = np.sqrt(expval2 - (expval ** 2))
            plt.plot(self.oo_vals[ct], sigma, "ro")
            ct += 1
        plt.xlabel("$R_{OO}$ $(\AA)$")
        plt.ylabel("$\sigma_{\Psi{OH}}$ $(\AA)$")
        plt.show()

    def calc_psi_max(self):
        for i in range(1, self.nwfns + 1):
            wfns = np.load(self.wfn_fname + str(i) + ".npy")  # "dimer_dvrwfns_"
            psi = wfns[:, 0]
            psi = np.absolute(psi)  # correct for phase
            psi_max_loc = np.argmax(psi)
            coord_psi_max = self.x[psi_max_loc]
            plt.plot(self.oo_vals[i - 1], coord_psi_max, "bo")
        plt.xlabel("$R_{OO} (\AA)$")
        plt.ylabel("$r_{OH} (\AA)$ of $\Psi^{max}$")
        # plt.savefig(fname="dimer_psimax_vs_oo_fig", dpi=250)
        plt.show()
        return None




exp_ob = DVRAnalysis(18, "dimer_gaussian_data_run3/dimer_dvrwfns3_", "dimer_gaussian_data_run3/xOH3.npy", "dimer_gaussian_data_run3/oo_steps_dimer3.npy")
# exp_ob = DVRAnalysis(18, "dimer_gaussian_data_run2/dimer_dvrwfns2_", "dimer_gaussian_data_run1/xOH2.npy", "dimer_gaussian_data_run2/oo_steps_dimer2.npy")

exp_ob.calc_expvals()
exp_ob.calc_stand_dev()
exp_ob.calc_psi_max()

#  -------------------------------------------------------------------


def wfn_plot(numb_wfns):
    for i in range(1, numb_wfns + 1):
        a = np.load(file="dimer_gaussian_data_run3/dimer_dvrwfns3_" + str(i) + ".npy")
        psi = a[:, 0] ** 2
        x = np.load(file="dimer_gaussian_data_run3/xOH3.npy")
        x /= 1.88973
        plt.ylabel("$P$")
        plt.xlabel("$r_{OH}$ $(\AA)$")
        plt.plot(x, psi)
    plt.show()
    return None

foo = wfn_plot(16)


def potential_plots(n_oos, n_ohs):

    all_engs = np.zeros(n_oos * n_ohs)
    a = n_ohs - n_ohs
    b = n_ohs

    for i in range(1, n_oos + 1):
        eng = np.load(file="dimer_gaussian_data_run3/dimer_Es3_" + str(i) + ".npy")
        # eng = np.load(file="dimer_gaussian_data_run2/dimer_Es2_" + str(i) + ".npy")

        all_engs[a:b] = eng
        a += n_ohs
        b += n_ohs
    min = np.amin(all_engs)

    for i in range(1, n_oos + 1):
        engs = np.load(file="dimer_gaussian_data_run3/dimer_Es3_" + str(i) + ".npy")
        # engs = np.load(file="dimer_gaussian_data_run2/dimer_Es2_" + str(i) + ".npy")

        engs -= min
        engs *= 219474.63
        x = np.load(file="dimer_gaussian_data_run3/dimer_r4oh3.npy")
        # x = np.load(file="dimer_gaussian_data_run2/dimer_r4oh2.npy")

        plt.plot(x, engs)
        plt.xlabel("$r_{OH}$ $(\AA)$")
        plt.ylabel("$cm^{-1}$")
        # plt.ylim([0, 7500])
    plt.show()
    return None


food = potential_plots(18, 16)
print("k")




