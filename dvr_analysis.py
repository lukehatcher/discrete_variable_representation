import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})


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
        :param oo_fname: filename of OO steps  # comes in angst
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

    def expvals_vs_roo(self):
        for expval, oo in zip(self.myexp_val_list, self.oo_vals):
            plt.plot(oo, expval, "ro")
        plt.xlabel("$R_{OO}$ $(\AA)$")
        plt.ylabel("$<r_{OH}>$")
        plt.show()
        return None

    def calc_stand_dev(self):
        ct = 0
        for expval, expval2 in zip(self.myexp_val_list, self.myexp_val_list2):
            sigma = np.sqrt(expval2 - (expval ** 2))
            plt.plot(self.oo_vals[ct], sigma, "ro")
            ct += 1
        plt.xlabel("$R_{OO}$ $(\AA)$")
        plt.ylabel("$\sigma_{\Psi{OH}}$ $(\AA)$")
        # plt.savefig(fname="sigma33", dpi=500, bbox_inches="tight")   ##### delete
        plt.show()
        return None

    def calc_psi_max(self):
        for i in range(1, self.nwfns + 1):
            wfns = np.load(self.wfn_fname + str(i) + ".npy")  # "dimer_dvrwfns_"
            psi = wfns[:, 0]
            psi = np.absolute(psi)  # correct for phase
            psi_max_loc = np.argmax(psi)
            coord_psi_max = self.x[psi_max_loc]
            # plt.rcParams.update({'font.size': 15})
            plt.plot(self.oo_vals[i - 1], coord_psi_max, "bo")
        plt.xlabel("$R_{OO} (\AA)$")
        plt.ylabel("$r_{OH} (\AA)$ of $\Psi^{max}$")
        # plt.savefig(fname="psi_max33", dpi=500, bbox_inches="tight")  ##### delete
        plt.show()
        return None


exp_ob = DVRAnalysis(18, "dimer_gaussian_data_run5/dimer_dvrwfns5_", "dimer_gaussian_data_run5/xOH5.npy", "dimer_gaussian_data_run5/oo_steps_dimer5.npy")
exp_ob.calc_expvals()
# exp_ob.expvals_vs_roo()
exp_ob.calc_stand_dev()
exp_ob.calc_psi_max()

#  -------------------------------------------------------------------
list = [1, 9, 18]
#  -------------------------------------------------------------------


def wfn_plot(numb_wfns):
    # for i in range(1, numb_wfns + 1)
    for i in list:
        a = np.load(file="dimer_gaussian_data_run5/dimer_dvrwfns5_" + str(i) + ".npy")
        psi = a[:, 0] ** 2
        x = np.load(file="dimer_gaussian_data_run5/xOH5.npy")
        x /= 1.88973
        plt.ylabel("$P$")
        plt.xlabel("$r_{OH}$ $(\AA)$")
        plt.plot(x, psi, label="OO dist." + str(i))
    plt.legend()
    # plt.savefig(fname="wfns3_3", dpi=500, bbox_inches="tight")  ##### delete
    plt.show()
    return None

foo = wfn_plot(120)


def potential_plots(n_oos, n_ohs):
    all_engs = np.zeros(n_oos * n_ohs)
    a = n_ohs - n_ohs
    b = n_ohs

    # for i in range(1, n_oos + 1):
    for i in list:
        eng = np.load(file="dimer_gaussian_data_run5/dimer_Es5_" + str(i) + ".npy")  #EDIT EACH TRIAL

        all_engs[a:b] = eng
        a += n_ohs
        b += n_ohs
    min = np.amin(all_engs)

    # for i in range(1, n_oos + 1):
    for i in list:
        engs = np.load(file="dimer_gaussian_data_run5/dimer_Es5_" + str(i) + ".npy")  #EDIT EACH TRIAL
        engs -= min
        engs *= 219474.63
        x = np.load(file="dimer_gaussian_data_run5/dimer_r4oh5.npy")  #EDIT EACH TRIAL

        plt.plot(x, engs, label="OO dist. " + str(i))
        plt.xlabel("$r_{OH}$ $(\AA)$")
        plt.ylabel("$cm^{-1}$")
        # plt.ylim([0, 2000])
        plt.xlim([.7, 1.5])
    plt.legend()
    plt.show()
    return None

food = potential_plots(18, 120)


def minimum_vs_oo(n_points, n):
    """ interp pots to 2k"""
    oo_grid = np.load(file="dimer_gaussian_data_run5/oo_steps_dimer5.npy")
    mins_list = []
    for i in range(1, n + 1):
        pots = np.load(file="dimer_gaussian_data_run5/dimer_Es5_" + str(i) + ".npy")
        a = np.amin(pots)
        mins_list.append(a)
        minE = min(mins_list)

    for i in range(1, n + 1):
        pots_b4 = np.load(file="dimer_gaussian_data_run5/dimer_Es5_" + str(i) + ".npy")
        pots_b4 -= minE
        pots_b4 *= 219474.63

        oh_grid = np.load(file="dimer_gaussian_data_run5/xOH5.npy")
        oh_grid_ext = np.linspace(start=np.amin(oh_grid), stop=np.amax(oh_grid), num=len(pots_b4), endpoint=True)
        f = interpolate.interp1d(oh_grid_ext, pots_b4, kind="cubic")

        new_oh = np.linspace(start=np.amin(oh_grid), stop=np.amax(oh_grid), num=n_points, endpoint=True)
        new_pots = f(new_oh)

        minimum_e_index = np.argmin(new_pots)
        plt.plot(oo_grid[i - 1], new_oh[minimum_e_index], "bo")

    plt.xlabel("$r_{OO}$ $(\AA)$")
    plt.ylabel("$(\AA)$")
    plt.show()

    return new_pots, oh_grid_ext


a = minimum_vs_oo(2000, 18)


def freq_vs_oo(n):
    oo_grid = np.load(file="dimer_gaussian_data_run5/oo_steps_dimer5.npy")
    for i in range(1, n + 1):
        evals = np.load(file="dimer_gaussian_data_run5/dimer_evals5_" + str(i) + ".npy")
        eng = evals[1] - evals[0]
        eng *= 219474.63
        plt.plot(oo_grid[i - 1], eng, "go")
    plt.xlabel("$r_{OO}$ $(\AA)$")
    plt.ylabel("$cm^{-1}$")
    plt.show()
    return

b = freq_vs_oo(18)
