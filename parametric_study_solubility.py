from parameters import parameters
from solve_H_transport import run_H_transport


if __name__ == "__main__":

    test_2_1 = 4e23
    test_2_2 = 2e21

    test_3_1 = 2.82842712e22

    test_5_1 = 7.52120619e21
    test_5_2 = 1.06365918e23

    test_9_1 = 3.87845489e21
    test_9_2 = 1.45853295e22
    test_9_3 = 5.48496351e22
    test_9_4 = 2.06267708e23

    test_values = [
        test_2_1,
        test_2_2,
        test_3_1,
        test_5_1,
        test_5_2,
        test_9_1,
        test_9_2,
        test_9_3,
        test_9_4,
    ]
    folder = "Results/parametric_studies/varying_lipb_S_0"
    # E_S = 0.01399  # average value
    E_S = 0.133
    for S_0 in test_values:
        parameters["materials"][-1]["S_0"] = S_0
        parameters["materials"][-1]["E_S"] = E_S
        parameters["exports"]["xdmf"]["folder"] = folder + "/S_0={:.1e}".format(S_0)
        parameters["exports"]["derived_quantities"][
            "folder"
        ] = folder + "/S_0={:.1e}".format(S_0)
        print("Current step is S_0 = {:.1e}".format(S_0))
        run_H_transport(parameters, S_0=S_0, E_S=E_S, log_level=40)
