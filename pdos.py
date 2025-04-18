import json
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.build import bulk, surface
from ase.visualize import view
from gpaw import GPAW
from gpaw.occupations import FermiDirac
from ase.units import Hartree
from gpaw.utilities.dos import RestartLCAODOS, fold

def load_config_files():
    with open("input_config.json", "r") as f:
        surface_configs = json.load(f)
    with open("gpaw_p_config.json", "r") as f:
        calc_config = json.load(f)
    return surface_configs, calc_config


def run_pdos_simulation(config, calc_config):
    element = config["element"]
    structure = config["structure"]
    a = config["lattice_constant"]
    h, k, l = config["h"], config["k"], config["l"]
    ort = config["orthorhombic"]

    atoms_bulk = bulk(element, crystalstructure=structure, a=a, orthorhombic=ort)
    layers = 1 if ort else 2
    slab = surface(atoms_bulk, indices=(h, k, l), layers=layers, periodic=True)
    slab *= (2, 2, 1)
    slab.center(vacuum=5, axis=2)
    slab.rotate(180, 'x', center='COU')

    view(slab)

    mode = calc_config["pdos"]["mode"]
    h_real = calc_config["pdos"]["grid_spacing"]
    kpts = tuple(calc_config["pdos"]["kpts"])
    convergence = calc_config["pdos"]["convergence"]
    basis = calc_config["pdos"]["basis"]
    xc = calc_config["pdos"]["xc"]
    symmetry = calc_config["pdos"]["symmetry"]
    fermi_smearing = calc_config["pdos"]["fermi_smearing"]

    calc = GPAW(mode=mode,
                h=h_real,
                kpts={'size': kpts, 'gamma': True},
                xc=xc,
                basis=basis,
                symmetry=symmetry,
                convergence=convergence,
                occupations=FermiDirac(width=fermi_smearing),
                txt=f"{element}_{h}{k}{l}_pdos.txt")

    slab.calc = calc
    slab.get_potential_energy()
    calc.write(f"{element}_{h}{k}{l}_pdos.gpw")

    ef = calc.get_fermi_level()
    dos = RestartLCAODOS(calc)

    energies, weights = dos.get_subspace_pdos(range(51))
    e, w = fold(energies * Hartree, weights, 2000, 0.1)

    e, Ni_s_pdos = dos.get_subspace_pdos([0, 1])
    e, Ni_s_pdos = fold(e * Hartree, Ni_s_pdos, 2000, 0.1)
    e, Ni_p_pdos = dos.get_subspace_pdos([2, 3, 4])
    e, Ni_p_pdos = fold(e * Hartree, Ni_p_pdos, 2000, 0.1)
    e, Ni_d_pdos = dos.get_subspace_pdos([5, 6, 7, 8, 9])
    e, Ni_d_pdos = fold(e * Hartree, Ni_d_pdos, 2000, 0.1)

    e, O_s_pdos = dos.get_subspace_pdos([25])
    e, O_s_pdos = fold(e * Hartree, O_s_pdos, 2000, 0.1)
    e, O_p_pdos = dos.get_subspace_pdos([26, 27, 28])
    e, O_p_pdos = fold(e * Hartree, O_p_pdos, 2000, 0.1)

    w_max = []
    for i in range(len(e)):
        if (-5.5 <= e[i] - ef <= 5.5):
            w_max.append(w[i])

    w_max = np.asarray(w_max)

    # Plot things:
    plt.figure()
    plt.plot(e - ef, w, label='Total', c='k', lw=2, alpha=0.7)
    plt.plot(e - ef, O_s_pdos, label='O-s', c='g', lw=2, alpha=0.7)
    plt.plot(e - ef, O_p_pdos, label='O-p', c='b', lw=2, alpha=0.7)
    plt.plot(e - ef, Ni_s_pdos, label='Ni-s', c='y', lw=2, alpha=0.7)
    plt.plot(e - ef, Ni_p_pdos, label='Ni-p', c='c', lw=2, alpha=0.7)
    plt.plot(e - ef, Ni_d_pdos, label='Ni-d', c='r', lw=2, alpha=0.7)

    plt.axis(ymin=0.0, ymax=np.max(w_max), xmin=-5.5, xmax=5.5)
    plt.xlabel(r'$\epsilon - \epsilon_F \ \rm{(eV)}$')
    plt.ylabel('DOS')
    plt.legend(loc=1)
    plt.savefig(f"{h},{k},{l}_pdos.png")
    #plt.show()


if __name__ == "__main__":
    surface_configs, calc_config = load_config_files()
    for config in surface_configs:
        run_pdos_simulation(config, calc_config)
