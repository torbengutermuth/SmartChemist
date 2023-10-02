from typing import Type

from smart_chemist.models import PatternMatchingJob, PatternMatchingInputModel


def make_pattern_matching_job(request_data, nof_molecules_allowed: int = 100) -> Type[PatternMatchingJob]:
    """Create a PatternMatchingJob from the request_data.

    :param request_data: dict of request data.
    :return: A new PatternMatchingJob object.
    """
    job = PatternMatchingJob()
    job.save()
    job.input_info = PatternMatchingInputModel(parent_pattern_matching_job=job)

    job.input_info.input_max_nof_molecules_allowed = nof_molecules_allowed
    if "smiles" in request_data:
        job.input_info.input_format = 'smiles_list'
        job.input_info.input_string = request_data["smiles"]
    elif "molecule_file" in request_data:
        file_name = request_data["molecule_file"].name
        if file_name.endswith(".smi"):
            job.input_info.input_format = '.smi'
        elif file_name.endswith(".sdf"):
            job.input_info.input_format = '.sdf'
        else:
            raise RuntimeError("Wrong filetype")

        # write file content to database model
        job.input_info.input_string = request_data["molecule_file"].file.read().decode("ascii")
    else:
        job.input_info.input_format = 'smiles_list'
        job.input_info.input_string = """
        P1CCCCC1 Poster Example,
        S1NCNC1 Poster Example,
        N1C2C(NCCCC2)CC1 Poster Example,
        C(=Cc1ccccc1)c1ccccc1 Poster Example,
        O1C2C(NC1)CNCC2 Poster Example,
        o1c2nnccc2cc1 Poster Example,
        C(=C/C(=O)O)\C(=O)O Poster Example,
        S1NC2N(CN1)CCC2 Poster Example,
        n1ncc2nccccc12 Poster Example,
        CNC(N)=N Poster Example,
        N1C2C(NCC2)CC1 Poster Example,
        O1CC2NCCCC2C1 Poster Example,
        c1cc(cnc1)C(=O)N Poster Example,
        S1CNNCC1 Poster Example,
        O1C=COC=C1 Poster Example,
        S1C2NCCCC2CC1 Example,
        Nc1ncnc2c1ncn2[C@@H]1O[C@H](CSCC[C@H](N)C(=O)O)[C@@H](O)[C@H]1O 3QOX_SAH_A_417,
        S(=O)(=O)(N)c1c(N)cc(c(S(=O)(=O)N)c1)C(F)(F)F 2POW_I7C_A_1000,
        S(=O)(=O)(N)c1cc2c(cc1)CCNC2 1HNN_SKF_A_3001,
        O1c2c(OC)c(OC)cc(c2C=C[C@H]1C3CC3)Cc4c(nc(nc4)N)N 3FYW_XCF_X_300,
        O(c1c(OC)cc2c(c1)C3=NNC(c4ccc(c5ccc(O)cc5)cc4)=C3C2)C 4FSQ_HK3_A_300,
        S1C(=NC=C1)C2=Nc3c(N2c4cc5nc(nc(N)c5cc4)N)cc(OC)c(OC)c3 4LEK_1DN_X_202,
        S(=O)(=O)(N)C=1SC(=NN1)N 2HNC_1SA_A_265,
        S(=O)(=O)(N)c1ccc(cc1)C 6GM9_4J8_A_305,
        S(=O)(=O)(N)c1ccc(cc1)CC 6HQX_4JC_A_305,
        O=C(Nc1cc2c(cc1)cccc2)[C@H]3CN(C=4N5N=CN=C5N=C(C4)C)CCC3 5TZ3_7OM_C_1001,
        O(c1c(O)cc2c(c1)[C@@H]3[C@H]([C@H]4[C@@]([C@@H](O)CC4)(CC3)C)CC2)C 1LHW_ESM_A_301,
        S1(=O)(=O)N=C(O)c2c1cccc2 4RIV_LSA_A_302,
        S(=O)(=O)(N)C=1SC=2S(=O)(=O)[C@H](C[C@H](NCC)C2C1)C 4M2U_ETS_A_302,
        S(=O)(=O)(N)C=1SC(=NN1)N 2HNC_1SA_A_265,
        Clc1c(Cl)cccc1c2nc(nc(NC3CC3)c2)N 5MZG_2GE_A_201,
        S(=O)(=O)(N)c1ccc(cc1)C(=O)O 6RFH_4SO_A_305,
        Clc1c(S(=O)(=O)N)cc(S(=O)(=O)N)c(N)c1 2POV_I7B_A_1000,
        S(=O)(=O)(N)c1ccc(cc1)C(=O)O 6RFH_4SO_A_305,
        O=C1N(Cc2ccc(cc2)CN3C(=N[C@@](C3=O)(CC4CCCCC4)CCC5CCCCC5)N)C[C@@H](N1)CCC 3L5E_BDW_B_455,
        O=C(NC(=O)c1cc2c(cc1)cccc2)N[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO 3ZCT_VMP_A_1840,
        S(=O)(=O)(N)C=1SC(=NC(=O)C)N(N1)C 4K0Z_MZM_A_308,
        Clc1cc2N=C(Oc2cc1)N3CCN(C(=O)c4c(N5N=CC=N5)ccc(c4)C)[C@@H](CC3)C 6TO7_SUV_A_401,
        ClC=1SC(C=2ON=C(C2)CN3c4c(cccc4)C=C3C(=O)NC5CCN(C(C)C)CC5)=CC1 2BOH_IIA_B_1,
        Clc1ccc(cc1)[C@H]2Nc3c(cccc3)C(=O)N2 5AKW_NKI_A_2163,
        S(=O)(=O)(N)c1ccc(cc1)CCN 2NNG_ZYX_A_301,
        S(=O)(=O)(N)c1ccc(cc1)C(=O)O 6RFH_4SO_A_305,
        S(=O)(=O)(c1nc(nc(c1)C)N)c2ccc(cc2)C 3B25_B2K_A_1,
        S1(=O)(=O)N(C[C@H](NCC)C2=C1SC(S(=O)(=O)N)=C2)CCCOC 4M2V_BZ1_A_302,
        Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C 3UU7_2OH_B_1,
        Fc1c(N2N=CC3=C2N(C(=O)C=C3Nc4c(ccc(c4)C(=O)NC5CC5)C)C)ccc(F)c1 3GFE_P37_A_361,
        S(=O)(=O)(N)c1c(N)cc(c(S(=O)(=O)N)c1)C(F)(F)F 2POW_I7C_A_1000,
        S(=O)(=O)(OC[C@@]12OC(O[C@H]2[C@@H]3OC(O[C@@H]3CO1)(C)C)(C)C)N 3HKU_TOR_A_300"""
    job.input_info.save()
    job.save()
    return job
