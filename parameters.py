import sympy as sp
import FESTIM
import properties

# IDs for volumes and surfaces (must be the same as in xdmf files)

id_lipb = 6
id_W = 7
id_structure = 8
id_baffle = 9
id_pipe_1 = 10
id_pipe_2 = 11
id_pipe_3 = 12
id_pipe_4 = 13
id_pipe_5 = 14
id_pipe_6 = 15
id_pipe_7 = 16
id_pipe_8 = 17
id_pipe_9 = 18
id_pipe_10 = 19
id_eurofers = [
    id_structure,
    id_baffle,
    id_pipe_1,
    id_pipe_2,
    id_pipe_3,
    id_pipe_4,
    id_pipe_5,
    id_pipe_6,
    id_pipe_7,
    id_pipe_8,
    id_pipe_9,
    id_pipe_10,
]

id_inlet = 20
id_outlet = 21
id_first_wall = 22

id_eurofer_coolant_interface = 39
id_pipe_1_coolant_interface = 40
id_pipe_2_coolant_interface = 41
id_pipe_3_coolant_interface = 42
id_pipe_4_coolant_interface = 43
id_pipe_5_coolant_interface = 44
id_pipe_6_coolant_interface = 45
id_pipe_7_coolant_interface = 46
id_pipe_8_coolant_interface = 47
id_pipe_9_coolant_interface = 48
id_pipe_10_coolant_interface = 49
ids_pipe_coolant_interface = [
    id_pipe_1_coolant_interface,
    id_pipe_2_coolant_interface,
    id_pipe_3_coolant_interface,
    id_pipe_4_coolant_interface,
    id_pipe_5_coolant_interface,
    id_pipe_6_coolant_interface,
    id_pipe_7_coolant_interface,
    id_pipe_8_coolant_interface,
    id_pipe_9_coolant_interface,
    id_pipe_10_coolant_interface,
]

materials_eurofer = []

for id_mat in id_eurofers:
    materials_eurofer.append(
        {
            # # EUROfer
            "D_0": properties.D_0_eurofer,
            "E_D": properties.E_D_eurofer,
            "S_0": properties.S_0_eurofer,
            "E_S": properties.E_S_eurofer,
            "thermal_cond": properties.thermal_cond_eurofer,
            "heat_capacity": properties.Cp_eurofer,
            "rho": properties.rho_eurofer,
            "id": id_mat,
        }
    )

mesh_folder = "meshes"

parameters = {
    "mesh_parameters": {
        "mesh_file": mesh_folder + "/mesh_domains_fine.xdmf",
        "cells_file": mesh_folder + "/mesh_domains_fine.xdmf",
        "facets_file": mesh_folder + "/mesh_boundaries_fine.xdmf",
    },
    "materials": [
        {
            # Tungsten
            "D_0": properties.D_0_W,
            "E_D": properties.E_D_W,
            "S_0": properties.S_0_W,
            "E_S": properties.E_S_W,
            "thermal_cond": properties.thermal_cond_W,
            "heat_capacity": properties.Cp_W,
            "rho": properties.rho_W,
            "id": id_W,
        },
        *materials_eurofer,
        {
            # LiPb
            "D_0": properties.D_0_lipb,
            "E_D": properties.E_D_lipb,
            "S_0": properties.S_0_lipb,
            "E_S": properties.E_S_lipb,
            "thermal_cond": properties.thermal_cond_lipb,
            "heat_capacity": properties.Cp_lipb,
            "rho": properties.rho_lipb,
            "id": id_lipb,
        },
    ],
    "traps": [
        {
            "E_k": properties.E_D_W,
            "k_0": properties.D_0_W / (1.1e-10**2 * 6 * properties.atom_density_W),
            "E_p": 0.87,
            "p_0": 1e13,
            "density": 1.3e-3 * properties.atom_density_W,
            "materials": [id_W],
        },
        {
            "k_0": 4.1e-7 / (1.1e-10**2 * 6 * properties.atom_density_W),
            "E_k": 0.39,
            "p_0": 1e13,
            "E_p": 1.0,
            "density": 4e-4 * properties.atom_density_W,
            "materials": [id_W],
        },
        {
            "E_k": properties.E_D_eurofer,
            "k_0": properties.D_0_eurofer
            / (1.1e-10**2)
            * 0.8165
            / properties.atom_density_eurofer,
            "E_p": properties.trap_energy_eurofer,
            "p_0": 1e13,
            "density": properties.trap_density_eurofer,
            "materials": id_eurofers,
        },
    ],
    "temperature": {
        "type": "solve_stationary",
        "boundary_conditions": [
            {"type": "flux", "value": 0.5e6, "surfaces": id_first_wall},
            {"type": "dc", "value": 598.15, "surfaces": id_inlet},
            {
                "type": "convective_flux",
                "h_coeff": 2.2308e04,
                "T_ext": 585.35,
                "surfaces": id_eurofer_coolant_interface,
            },
            {
                "type": "convective_flux",
                "h_coeff": 1.5865e04,
                "T_ext": 584.9,
                "surfaces": ids_pipe_coolant_interface,
            },
        ],
        "source_term": [
            {
                "value": 23.2e06
                * sp.exp(-71.74 * FESTIM.x)
                * (0.001 * (FESTIM.x < 0.02)),
                "volume": id_W,
            },
            {
                "value": (
                    (9.6209e06 * sp.exp(-12.02 * FESTIM.x)) * (FESTIM.x < 0.15)
                    + (4.7109e06 * sp.exp(-7.773 * FESTIM.x)) * (FESTIM.x >= 0.15)
                )
                * (0.001 * (FESTIM.x < 0.05)),
                "volume": id_eurofers,
            },
            {
                "value": (
                    (3.9108e05 * FESTIM.x ** (-1.213)) * (FESTIM.x < 0.15)
                    + (8.4629e06 * sp.exp(-5.485 * FESTIM.x)) * (FESTIM.x >= 0.15)
                )
                * (
                    0.1 * (FESTIM.x < 0.02)
                    + 1 * (0.02 <= FESTIM.x) * (FESTIM.x <= 0.1)
                    + 1.20 * (0.1 <= FESTIM.x) * (FESTIM.x <= 0.15)
                    + 0.5 * (0.15 <= FESTIM.x) * (FESTIM.x <= 0.25)
                    + 1e-05 * (FESTIM.x > 0.25)
                ),
                "volume": id_lipb,
            },
        ],
    },
    "source_term": [
        {
            "value": 6.022e23
            * 1e6
            * (
                1.044e-11 * sp.exp(-0.2182 * FESTIM.x * 1e2)
                + 6.514e-12 * sp.exp(-0.04106 * FESTIM.x * 1e2)
            ),
            "volumes": id_lipb,
        },
    ],
    "boundary_conditions": [
        {
            "type": "dc",
            "surfaces": [id_inlet],
            "value": 0,
        },
        {
            "type": "recomb",
            "surfaces": [*ids_pipe_coolant_interface, id_eurofer_coolant_interface],
            "Kr_0": properties.Kr_0_eurofer,
            "E_Kr": properties.E_Kr_eurofer,
            "order": 2,
        },
        {
            "type": "dc_imp",
            "implanted_flux": 1e20,
            "implantation_depth": 3e-9,
            "D_0": properties.D_0_W,
            "E_D": properties.E_D_W,
            "surfaces": [id_first_wall],
        },
    ],
    "solving_parameters": {
        "type": "solve_stationary",
        "final_time": 86400 * 15,
        "initial_stepsize": 1,
        "adaptive_stepsize": {
            "stepsize_change_ratio": 1.1,
            "t_stop": 1e618,
            "stepsize_stop_max": 1 / 10,
            "dt_min": 1e-4,
        },
        "newton_solver": {
            "absolute_tolerance": 1e12,
            "relative_tolerance": 1e-8,
            "maximum_iterations": 50,
        },
        "traps_element_type": "DG",
    },
    "exports": {
        "xdmf": {
            "functions": ["solute", "1", "2", "3", "retention"],
            "labels": [
                "solute",
                "trap_W_1",
                "trap_W_2",
                "trap_eurofer_1",
                "retention",
            ],
            "folder": "Results",
            "nb_iterations_between_exports": 1,
        },
        "derived_quantities": {
            "total_volume": [
                {"volumes": [id_W, *id_eurofers, id_lipb], "field": "solute"},
                {"volumes": [id_W], "field": "1"},
                {"volumes": [id_W], "field": "2"},
                {"volumes": [id_eurofers], "field": "3"},
                {"volumes": [id_W, *id_eurofers], "field": "retention"},
            ],
            "surface_flux": [
                {
                    "surfaces": [
                        *ids_pipe_coolant_interface,
                        id_eurofer_coolant_interface,
                        id_first_wall,
                    ],
                    "field": "solute",
                }
            ],
            "total_surface": [{"surfaces": [id_outlet], "field": "solute"}],
            "file": "derived_quantities.csv",
            "folder": "Results",
            "nb_iterations_between_exports": 1,
        },
    },
}
