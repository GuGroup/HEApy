# ytk
# ytk

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
from matplotlib.ticker import FuncFormatter

def draw_volcano_plot(
    uma_result_path: Path
    ) -> None:
    from ..utils.globals import COLORMAP
    # 1. Common things
    bins = 20
    alpha = 0.7
    plt.rcParams['font.family'] = "Arial"
    
    # 2. Load data
    # UMA 
    energies_UMA = json.load(open(uma_result_path, "r"))

    #NOH->HNOH
    blue = [[-8.0879, -5.2776], [-1.5, -0.3707]]
    red = [[-5.2776, -2.724], [-0.3707, -0.2450]]
    black = [[-2.724, 0.2502], [-0.2450, 0.6523]]
    cyan = [[0.2502, 4.9666], [0.6523, -0.3344]]
    x_dot = [[-9, 5], [0.69, 0.69]]
    y_dot = [[-3.476, -3.476], [-1.5, 2]]
    scatter = [[-3.9770, -4.4978, -4.8479, -5.2032, -5.2591], [-0.2317, -0.3311, -0.7152, -0.3907, -0.3708]]

    # 3. Data processing
    # UMA
    for key in energies_UMA.keys():
        del_idxs = []
        for idx, energy in enumerate(energies_UMA[key]):
            if energy > 1000 or energy < -6000:
                del_idxs.append(idx)
        for del_idx in sorted(del_idxs, reverse=True):
            del energies_UMA[key][del_idx]

    # 4. Drawing Figure
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    # Volcano
    ax.plot(blue[0], blue[1], color="dodgerblue", linewidth=2)
    ax.plot(red[0], red[1], color="red", linewidth=2)
    ax.plot(black[0], black[1], color="black", linewidth=2)
    ax.plot(cyan[0], cyan[1], color="cyan", linewidth=2)
    ax.plot(x_dot[0], x_dot[1], color="dimgray", linewidth=2, linestyle="--")
    ax.plot(y_dot[0], y_dot[1], color="dimgray", linewidth=2, linestyle="--")
    # UMA
    x_grid = np.linspace(-6, 2, 1000)
    for binding_site_atom in energies_UMA.keys():
        data = np.array(energies_UMA[binding_site_atom])

        kde = gaussian_kde(data, bw_method='scott')  # or try bw_method=0.3
        y = kde(x_grid)

        if binding_site_atom == "total":
            y_mul = len(energies_UMA[key]) / 0.7
        else:
            y_mul = len(energies_UMA[key]) / 5

        ax2.plot(
            x_grid, 
            y * y_mul,  # multiply to match histogram height scale
            color=COLORMAP[binding_site_atom],
            label=binding_site_atom,
            linewidth=3,
            alpha=alpha)
        
        # Optional black outline
        ax2.plot(
            x_grid, 
            y * y_mul,
            color='black',
            linewidth=1)
        ''' 
        ax2.legend(loc="upper right", 
                   framealpha=0,
                   prop={"weight": "semibold", "size": 13},
                   bbox_to_anchor=(1, 0.98))
        '''
    symbol_legend = [
        Patch(facecolor='gray', edgecolor='black', label='DFT'),
        Line2D([0], [0], color='black', linewidth=3, label='UMA')
    ]

    #scatter
    ax.scatter(scatter[0], scatter[1], marker="o", color="lavender", edgecolor="black", zorder=10)
    # 5. Settings

    # lim
    ax.set_xlim((-9, 5))
    ax.set_ylim((-1.5, 2.01))
    ax2.set_ylim((0, 800))

    # tick
    xticks = np.arange(-9, 5.1, 2)
    ax.set_xticks(xticks)

    yticks_dft = np.arange(-1.5, 2.1, 0.5)
    ax.set_yticks(yticks_dft)

    yticks_uma = np.arange(0, 801, 200)
    ax2.set_yticks(yticks_uma)

    ax.tick_params(axis="x", which='major', direction='in', length=6, width=2.5, labelsize=17, pad=9)
    ax.tick_params(axis="y", which='major', direction='in', length=6, width=2.5, labelsize=17)
    ax2.tick_params(axis="y", which='major', direction='in', length=6, width=2.5, labelsize=17)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('semibold')

    for label in ax2.get_yticklabels():
        label.set_fontweight("semibold")

    major_xticks = ax.get_xticks()
    if len(major_xticks) >= 2:
        major_step = major_xticks[1] - major_xticks[0]
        ax.xaxis.set_minor_locator(MultipleLocator(major_step / 2))

    major_yticks = ax.get_yticks()
    if len(major_yticks) >= 2:
        major_step = major_yticks[1] - major_yticks[0]
        ax.yaxis.set_minor_locator(MultipleLocator(major_step / 2))

    major_yticks = ax2.get_yticks()
    if len(major_yticks) >= 2:
        major_step = major_yticks[1] - major_yticks[0]
        ax2.yaxis.set_minor_locator(MultipleLocator(major_step / 2))

    ax.tick_params(which='minor', direction='in', length=3.5, width=2.5)
    ax2.tick_params(which='minor', direction='in', length=3.5, width=2.5)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    #ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.1f)}"))

    ax.set_xlabel("Gibbs free energy (eV)", size=25, weight="semibold", labelpad=8)
    ax.set_ylabel("Limiting potential (V)", size=23, weight="semibold", labelpad=10)
    ax2.set_ylabel("Frequency (UMA)", size=25, weight="semibold", labelpad=10)
    for spine in ax.spines.values():
        spine.set_linewidth(2.5)
    plt.tight_layout()
    plt.show()


def draw_gibbs_free_energy_distribution(
    dft_result_path: Path,
    uma_result_path: Path,
    bins: int=20,
    alpha: float=0.7,
    ) -> None:

    from ..utils.globals import COLORMAP
    
    #Common property 
    plt.rcParams['font.family'] = "Arial"
    
    # 1. load data each other
    # DFT
    energies_DFT = json.load(open(dft_result_path, "r"))
    
    # UMA
    energies_UMA = json.load(open(uma_result_path, "r"))
    
    # 2. Data exclusion
    # DFT
    for key in energies_DFT.keys():
        del_idxs = []
        for idx, energy in enumerate(energies_DFT[key]):
            if energy > -1.5:
                del_idxs.append(idx) 
        for del_idx in sorted(del_idxs, reverse=True):
            del energies_DFT[key][del_idx]

    # UMA
    for key in energies_UMA.keys():
        del_idxs = []
        for idx, energy in enumerate(energies_UMA[key]):
            if energy > -1.5 or energy < -6.6: 
                del_idxs.append(idx) 
        for del_idx in sorted(del_idxs, reverse=True):
            del energies_UMA[key][del_idx]
   
    # 3. Draw Figure
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    bins = np.linspace(-6, -2, 50)
    # DFT
    for binding_site_atom in energies_DFT.keys():
        if binding_site_atom == "total":
            ax.hist(energies_DFT[binding_site_atom], 
                    bins=bins,
                    color=COLORMAP[binding_site_atom],
                    density=False,
                    histtype="stepfilled",
                    alpha=alpha,
                    label=binding_site_atom)
            ax.hist(energies_DFT[binding_site_atom],
                    bins=bins,
                    density=False,
                    histtype='step',
                    color='black',
                    linewidth=1)
        else:
            ax.hist(energies_DFT[binding_site_atom],
                    bins=bins,
                    color=COLORMAP[binding_site_atom],
                    density=False,
                    histtype="stepfilled",
                    alpha=alpha,
                    label=binding_site_atom)
            ax.hist(energies_DFT[binding_site_atom],
                    bins=bins,
                    density=False,
                    histtype='step',
                    color='black',
                    linewidth=1)

    # UMA line
    x_grid = np.linspace(-6, -2, 1000)

    for binding_site_atom in energies_UMA.keys():
        data = np.array(energies_UMA[binding_site_atom])

        kde = gaussian_kde(data, bw_method='scott')  # or try bw_method=0.3
        y = kde(x_grid)
        
        if binding_site_atom == "total":
            y_mul = len(energies_UMA[key]) / 0.7
        else:
            y_mul = len(energies_UMA[key]) / 5

        ax2.plot(
            x_grid, 
            y * y_mul,  # multiply to match histogram height scale
            color=COLORMAP[binding_site_atom],
            label=binding_site_atom,
            linewidth=3,
            alpha=alpha)

    # Optional black outline
        ax2.plot(
            x_grid, 
            y * y_mul,
            color='black',
            linewidth=1)
        ax2.legend(
            loc="upper right",
            framealpha=0,
            prop={"weight": "semibold", "size": 13})

    symbol_legend = [
        Patch(facecolor='gray', edgecolor='black', label='DFT'),
        Line2D([0], [0], color='black', linewidth=3, label='UMA')
    ]

    # Add to upper left corner
    ax.legend(handles=symbol_legend,
              loc='upper left',
              bbox_to_anchor=(0.00, 0.98),
              framealpha=0,
              prop={'weight': 'semibold', 'size': 13})
    
    dft_ylim = 80
    uma_ylim = 800
    ax.set_ylim((0, dft_ylim))
    ax2.set_ylim((0, uma_ylim))

    yticks_dft = np.arange(0, dft_ylim + 1 , dft_ylim / 4)
    ax.set_yticks(yticks_dft)

    yticks_uma = np.arange(0, uma_ylim +1 , uma_ylim / 4)
    ax2.set_yticks(yticks_uma)

    ax.tick_params(which='major', direction='in', length=6, width=2.5, labelsize=17)
    ax2.tick_params(axis="y", which='major', direction='in', length=6, width=2.5, labelsize=17)
    #ax.tick_params(axis='x', bottom=False, top=False)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('semibold')

    for label in ax2.get_yticklabels():
        label.set_fontweight("semibold")

    major_xticks = ax.get_xticks()
    if len(major_xticks) >= 2:
        major_step = major_xticks[1] - major_xticks[0]
        ax.xaxis.set_minor_locator(MultipleLocator(major_step / 2))

    major_yticks = ax.get_yticks()
    if len(major_yticks) >= 2:
        major_step = major_yticks[1] - major_yticks[0]
        ax.yaxis.set_minor_locator(MultipleLocator(major_step / 2))

    major_yticks = ax2.get_yticks()
    if len(major_yticks) >= 2:
        major_step = major_yticks[1] - major_yticks[0]
        ax2.yaxis.set_minor_locator(MultipleLocator(major_step / 2))

    ax.tick_params(which='minor', direction='in', length=3.5, width=2.5)
    ax2.tick_params(which='minor', direction='in', length=3.5, width=2.5)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{int(y)}"))

    ax.set_xlabel("Gibbs free energy (eV)", size=25, weight="semibold", labelpad=8)
    ax.set_ylabel("Frequency (DFT)", size=25, weight="semibold", labelpad=10)
    ax2.set_ylabel("Frequency (UMA)", size=25, weight="semibold", labelpad=10)


    for spine in ax.spines.values():
        spine.set_linewidth(2.5)
    plt.tight_layout()
    #plt.savefig("Invitation_dist.png")
    plt.show()


def draw_parity_plot(
    ) -> None:
    from ase.io import read
    
    from ..core.dataset_utils import get_binding_atom
    from ..tests.test_utils import get_gibbs_free_energy
    from ..utils.PATHS import test_paths
    
    target = []
    prediction = []

    # Add data points
    for test_path in test_paths:
        slab_target = read(test_path / "slab" / "OUTCAR", index="-1")
        adslab_target = read(test_path / "adslab" / "OUTCAR", index="-1")
        
        binding_site_atom = get_binding_atom(read(test_path / "slab" / "POSCAR"))

        target_gibbs_free_energy = get_gibbs_free_energy(slab_target.calc.get_potential_energy(),
                                                         adslab_target.calc.get_potential_energy(),
                                                         hea_type="AlFeCoNiCu",
                                                         binding_site_atom=binding_site_atom
                                                         )

        slab_prediction = read(test_path / "slab" / "uma_prediction.traj", index="-1")
        adslab_prediction = read(test_path / "adslab" / "uma_prediction.traj", index="-1")

        prediction_gibbs_free_energy = get_gibbs_free_energy(slab_prediction.calc.get_potential_energy(),
                                                             adslab_prediction.calc.get_potential_energy(),
                                                             hea_type="AlFeCoNiCu",
                                                             binding_site_atom=binding_site_atom
                                                             )

        target.append(target_gibbs_free_energy)
        prediction.append(prediction_gibbs_free_energy)

    target = np.array(target)
    prediction = np.array(prediction)

    training_mae = np.mean(np.abs(target-prediction))
    
    plt.rcParams['font.family'] = "Arial"
    fig, ax = plt.subplots()
    ax.scatter(target, prediction, alpha=0.7, edgecolors=None, c="blue", label='test set')

    # y = x line
    min_val = min(target.min(), prediction.min())
    max_val = max(target.max(), prediction.max())
    buffer = 0.1 * (max_val - min_val)

    ax.plot([min_val-buffer, max_val+buffer], 
            [min_val-buffer, max_val+buffer], 
            linestyle='-', 
            c='black', 
            label='y = x')    
    ax.plot([min_val-buffer, max_val+buffer], 
            [min_val-buffer+0.1, max_val+buffer+0.1], 
            linestyle='--', 
            c='black', 
            label='+-0.1eV')
    ax.plot([min_val-buffer, max_val+buffer], 
            [min_val-buffer-0.1, max_val+buffer-0.1], 
            linestyle='--', 
            c='black')
        
    ax.set_xlim(min_val - buffer, max_val + buffer)
    ax.set_ylim(min_val - buffer, max_val + buffer)

    # Labels and annotations
    plt.xlabel('DFT Gibbs Free Energy (eV)', size=22, weight="semibold", labelpad=13)
    plt.ylabel('UMA Gibbs Free Energy (eV)', size=18, weight="semibold", labelpad=13)
    plt.legend(loc='lower right', prop={"size": 13})

    # R value text
    #plt.text(0.06, 0.98, f'MAE (training set) = {training_mae:.4f} eV\nMAE (test set) = {training_mae:.4f} eV', transform=plt.gca().transAxes, verticalalignment='top', fontsize=14, bbox=dict(facecolor='white', alpha=0))
    plt.text(0.06, 0.98, f'MAE (test set) = {training_mae:.3f} eV', transform=plt.gca().transAxes, verticalalignment='top', fontsize=14, bbox=dict(facecolor='white', alpha=0))

    ax.tick_params(which='major', direction='in', length=6, width=2.5, labelsize=17)
    #ax.tick_params(axis='x', bottom=False, top=False)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('semibold')

    major_xticks = ax.get_xticks()
    if len(major_xticks) >= 2:
        major_step = major_xticks[1] - major_xticks[0]
        ax.xaxis.set_minor_locator(MultipleLocator(major_step / 2))

    major_yticks = ax.get_yticks()
    if len(major_yticks) >= 2:
        major_step = major_yticks[1] - major_yticks[0]
        ax.yaxis.set_minor_locator(MultipleLocator(major_step / 2))

    ax.tick_params(which='minor', direction='in', length=3.5, width=2.5)

    for spine in ax.spines.values():
        spine.set_linewidth(2.5)

    plt.tight_layout()
    plt.show()
