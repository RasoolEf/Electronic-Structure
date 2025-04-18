import json
import numpy as np
from ase import Atoms
from ase.build import surface, bulk
from ase.visualize import view
from gpaw import GPAW, PW
from gpaw.occupations import FermiDirac
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid


def load_config_files():
    with open("input_config.json", "r") as f:
        surface_configs = json.load(f)
    with open("gpaw_config.json", "r") as f:
        calc_config = json.load(f)
    return surface_configs, calc_config


def run_simulation(config, calc_config):
    element = config["element"]
    structure = config["structure"]
    a = config["lattice_constant"]
    h, k, l = config["h"], config["k"], config["l"]
    ort = config["orthorhombic"]

    bulk_structure = bulk(name=element, crystalstructure=structure, a=a, orthorhombic=ort)
    slab = surface(bulk_structure, indices=(h, k, l), layers=1 if ort else 2, periodic=True)
    slab *= (2, 2, 1)
    slab.rotate(180, 'x', center='COU')
    slab.center(vacuum=5, axis=2)

    view(slab)

    # Ground state calculation
    mode_type = calc_config["ground_state"]["mode"]["type"]
    ecut = calc_config["ground_state"]["mode"]["ecut"]
    xc = calc_config["ground_state"]["xc"]
    kpts = tuple(calc_config["ground_state"]["kpts"])
    occupation = FermiDirac(calc_config["ground_state"]["fermi_smearing"])

    mode = PW(ecut) if mode_type == "PW" else None  # Add other modes later if needed
    calc = GPAW(
        mode=mode,
        xc=xc,
        kpts=kpts,
        occupations=occupation,
        txt=f'{element}_{h}{k}{l}.txt'
    )

    slab.calc = calc
    energy = slab.get_potential_energy()
    ef = calc.get_fermi_level()
    print(f'Fermi level of {element}({h},{k},{l}) is {ef:.3f} eV')

    energy, pdos = calc.get_orbital_ldos(a=0, angular='d')
    calc.write(f'{element}_{h}{k}{l}.gpw')

    # Band structure calculation
    nbands = calc_config["band_structure"]["nbands"]
    symmetry = calc_config["band_structure"]["symmetry"]
    kpts_path = np.array(config["kpoints"])
    npoints = calc_config["band_structure"]["npoints"]
    convergence = calc_config["band_structure"]["convergence"]

    bs_calc = GPAW(f'{element}_{h}{k}{l}.gpw').fixed_density(
        nbands=nbands,
        symmetry=symmetry,
        kpts={"path": kpts_path, "npoints": npoints},
        convergence=convergence
    )

    bs = bs_calc.band_structure()
    bs.plot(filename=f'{element}_{h}{k}{l}_band_structure.png', show=False, emax=10.0)

    # dband center calculation
    energy -= ef
    I = trapezoid(pdos, energy)
    center = trapezoid(pdos * energy, energy) / I
    width = np.sqrt(trapezoid(pdos * (energy - center) ** 2, energy) / I)

    plt.figure()
    plt.plot(energy, pdos)
    plt.xlabel('Energy (eV)')
    plt.ylabel('d-projected DOS')
    plt.title(f'd-band center = {center:.1f} eV, width = {width:.1f} eV')
    plt.savefig(f'{element}_{h}{k}{l}_dband_center.png')


if __name__ == "__main__":
    surface_configs, calc_config = load_config_files()
    for config in surface_configs:
        run_simulation(config, calc_config)
