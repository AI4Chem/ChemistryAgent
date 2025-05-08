import random
import re
from collections import Counter

from chemlib import Element

from utils import generate_random_string, generate_random_filepath
from utils.file_io import read_csv, read_JSON, write_JSON
from pubchempy import ELEMENTS, Compound


def is_valid_cas(cas_number):
    # CAS号的正则表达式
    cas_pattern = re.compile(r'^\d{2,7}-\d{2}-\d$')
    
    # 检查格式是否匹配
    if not cas_pattern.match(cas_number):
        return False

    # 拆分成三个部分
    parts = cas_number.split('-')
    first_part = parts[0]
    second_part = parts[1]
    check_digit = int(parts[2])
    
    # 计算校验位
    digits = list(first_part + second_part)
    digits.reverse()
    
    sum = 0
    for i, digit in enumerate(digits):
        sum += (i + 1) * int(digit)
    
    # 校验位应等于sum对10的余数
    return sum % 10 == check_digit

def get_random_compound():
    Flag = True
    while (Flag):
        compound = Compound.from_cid(random.randint(1,999999))
        for synonym in compound.synonyms:
            if is_valid_cas(synonym):
                compound_CAS = synonym
                Flag = False
                break
        compound_SMILES = compound.canonical_smiles
        compound_CID = compound.cid
        compound_MF = compound.molecular_formula
    return {
        "compound": compound,
        "CAS": compound_CAS,
        "SMILES": compound_SMILES,
        "CID": compound_CID,
        "MF": compound_MF
    }

def adjust_dict_value_sum_to_100(input_dict):
    total = sum(input_dict.values())
    factor = 100 / total
    adjusted_dict = {k: round(v * factor, 1) for k, v in input_dict.items()}
    
    # 计算调整后的总和，以确保总和为100
    adjusted_total = sum(adjusted_dict.values())
    
    # 修正由于舍入引起的误差
    difference = 100 - adjusted_total
    if difference != 0:
        # 找到值最大的键，并调整其值以弥补误差
        max_key = max(adjusted_dict, key=adjusted_dict.get)
        adjusted_dict[max_key] = round(adjusted_dict[max_key] + difference, 1)
    
    return adjusted_dict

def get_tool_param(tool_name: str, tool_param_dict: dict) -> list:
    generated_list = []
    for tool_param_name in tool_param_dict:
        tool_param_type = tool_param_dict[tool_param_name]["type"]
        match tool_name:
            # chemistrytools
            case "get_constants":
                return [random.choice(["avogadro_number","boltzmann_constant","electron_radius","faraday_constant","molar_gas_constant","neutron_mass","plancks_constant","speed_of_light","vacuum_permittivity"])]
            case "calculate_spectrum_similarity":
                original_list = list(range(1, 119))
                mz_top = random.sample(original_list, k=random.randint(10,60))
                intens_top = [round(random.uniform(0.1,100),4) for _ in mz_top]
                mz_bottom = random.sample(original_list, k=random.randint(10,60))
                intens_bottom = [round(random.uniform(0.1,100),4) for _ in mz_bottom]
                return [mz_top, intens_top, mz_bottom, intens_bottom]
            case "get_compound_CID":
                while True:
                    try:
                        compound = Compound.from_cid(random.randint(1,999999))
                        for synonym in compound.synonyms:
                            if not synonym[0].isdigit():
                                return [synonym]
                    except Exception as e:
                        print("get get_compound_CID param Error: ",e)

            # chem_lib
            case "calculate_compound_percentage_composition_by_mass":
                compound = get_random_compound()["MF"]
                all_elements = [ ELEMENTS[element_id] for element_id in ELEMENTS]
                element_list = [ element for element in all_elements if element in compound]
                return [compound, random.choice(element_list)]
            case "convert_compound_stoichiometry_amount":
                if tool_param_name == "unit":
                    generated_list.append(random.choice(["grams","moles","molecules"]))
                    generated_list.append(round(random.uniform(0.1, 100.0),4))
                    return generated_list
            case "get_empirical_formula_by_percent_composition":
                element_count = Counter(get_random_compound()["compound"].elements)
                element_count_dict = dict(element_count)
                element_proposition_dict = { element:element_count_dict[element]*Element(element).AtomicMass for element in element_count_dict}
                element_proposition_dict = adjust_dict_value_sum_to_100(element_proposition_dict)
                return [element_proposition_dict]
            case "analyze_combustion":
                while True:
                    element_count = Counter(get_random_compound()["compound"].elements)
                    element_count_dict = dict(element_count)
                    if "C" in element_count_dict and "H" in element_count_dict:
                        break
                random_num = round(random.uniform(0.1,2),1)
                return [round(element_count_dict["C"]*random_num*44,1), round(element_count_dict["H"]*random_num*18,1)]
            case "reaction_stoichiometry_amounts":
                if tool_param_name == "compound_number":
                    formula = generated_list[0]
                    formula_list = formula.split('-->')
                    compound_list = []
                    compound_list.extend(formula_list[0].split('+'))
                    compound_list.extend(formula_list[1].split('+'))
                    generated_list.append(random.randint(1,len(compound_list)+1))
                    generated_list.append(random.choice(["grams","moles","molecules"]))
                    generated_list.append(round(random.uniform(0.1, 100.0),4))
                    return generated_list
            case "limiting_reagent_of_reaction":
                if tool_param_name == "mode":
                    generated_list.append(random.choice(["grams","moles","molecules"]))
                    formula = generated_list[0]
                    formula_list = formula.split('-->')
                    rectant_list = formula_list[0].split('+')
                    generated_list.append([round(random.uniform(0.1, 100.0),4) for _ in rectant_list])
                    return generated_list
            case "acidity_calculation":
                return [round(random.uniform(0.0,14.0),2)]
            case "make_solution_by_grams_per_liter":
                if tool_param_name == "grams":
                    generated_list.append(round(random.uniform(0.1, 100.0),4))
                    generated_list.append(round(random.uniform(0.1, 100.0),4))
                    return generated_list
            case "perform_electrolysis":
                if tool_param_name == "n":
                    generated_list.append(random.randint(1,10))
                    generated_list.append(round(random.uniform(0.1, 100.0),4))
                    generated_list.append(round(random.uniform(0.1, 100.0),3))
                    return generated_list
            case "electromagnetic_wave_by_wavelength":
                return [round(random.uniform(0.1, 100.0),4)]
            case "electromagnetic_wave_by_frequency":
                return [round(random.uniform(0.1, 100.0),4)]
            case "electromagnetic_wave_by_energy":
                return [round(random.uniform(0.1, 100.0),4)]
            case "electromagnetic_wave_by_rydberg_equation":
                if tool_param_name == 'n1':
                    generated_list.append(random.randint(1,10))
                    generated_list.append(generated_list[-1] + random.randint(1,10))
                    return generated_list
            case "energy_of_hydrogen_orbital":
                return [random.randint(1,10)]

            case _:
                pass
        # 一般参数
        if tool_param_name.startswith("element"):
            generated_list.append(ELEMENTS[random.randint(1, 118)])
        elif tool_param_name.startswith("compound"):
            compound = get_random_compound()
            if tool_param_name.startswith("compound_SMILES"):
                generated_list.append(compound["SMILES"])
            elif tool_param_name.startswith("compound_CAS"):
                generated_list.append(compound["CAS"])
            elif tool_param_name.startswith("compound_CID"):
                generated_list.append(compound["CID"])
            else: # 默认为分子式
                generated_list.append(compound["MF"])
        elif tool_param_name.startswith("reaction_formula"):
            formula_list = read_csv("data/chemistry_reactions.csv")
            generated_list.append(formula_list[random.randint(1,len(formula_list))][0])
        else:
            print("出现未知参数！")
            print("Tool Name: ", tool_name)
            print("Tool Param Name: ", tool_param_name)
            print("Tool Param Type: ", tool_param_type)
    return generated_list

def get_matgen_tool_param(tool_name: str, tool_param_dict: dict) -> list:
    generated_list = []
    for tool_param_name in tool_param_dict:
        tool_param_type = tool_param_dict[tool_param_name]["type"]
        match tool_name:
            # Basic_Functionality
            case "get_atomic_mass_in_unit":
                if tool_param_name == "unit":
                    generated_list.append(random.choice(["kg","g","amu"]))
                    continue
            case "get_composition_properties":
                compound = get_random_compound()["MF"]
                all_elements = [ ELEMENTS[element_id] for element_id in ELEMENTS]
                element_list = [ element for element in all_elements if element in compound]
                return [compound, random.choice(element_list)]
            case "get_species_properties":
                if tool_param_name == "oxidation_state":
                    generated_list.append(random.choice([i for i in range(-7, 8) if i != 0]))  # random.randint(-7, 7)
                    continue
            case "create_cubic_lattice":
                if tool_param_name == "lattice_parameter":
                    generated_list.append(round(random.uniform(3.0, 6.0), 3))
                    continue
            case "create_structure":
                if tool_param_name == "species":
                    num = random.randint(3,7)
                    species = [ELEMENTS[random.randint(1, 118)] for _ in range(num)]
                    coordinates = [[round(random.uniform(-3.0, 3.0), 3) for _ in range(3)] for _ in range(num)]
                    generated_list.append(species)
                    generated_list.append(coordinates)
                    generated_list.append("structure_{}.pkl".format(generate_random_string(10)))
                    return generated_list
            case "modify_structure":
                if tool_param_name == "operation":
                    operation = random.choice(["make_supercell","delete","append","change","shift"])
                    generated_list.append(operation)
                    match operation:
                        case "make_supercell":
                            generated_list.append([random.randint(1, 10) for _ in range(random.randint(3,6))])
                        case "delete":
                            generated_list.append([random.randint(1, 6)])
                        case "append":
                            generated_list.append([ELEMENTS[random.randint(1, 118)], [round(random.uniform(-3.0, 3.0), 3) for _ in range(3)]])
                        case "change":
                            generated_list.append([ELEMENTS[random.randint(1, 118)]])
                        case "shift":
                            generated_list.append([random.randint(1,6), [round(random.uniform(-3.0, 3.0), 3) for _ in range(3)]])
                    return generated_list
                
            # Battery_materials_Analysis
            case "analyze_ion_diffusion_pathways":
                if tool_param_name == "min_slab_size":
                    generated_list.append(round(random.uniform(1.0, 3.0), 3))
                    continue
                elif tool_param_name == "min_vacuum_size":
                    generated_list.append(round(random.uniform(1.0, 3.0), 3))
                    continue

            # Catalysis_Studies
            case "model_catalyst_surface":
                if tool_param_name == "miller_index":
                    generated_list.append([random.randint(1, 10) for _ in range(random.randint(3,6))])
                    continue
                elif tool_param_name == "min_slab_size":
                    generated_list.append(round(random.uniform(1.0, 3.0), 3))
                    continue
                elif tool_param_name == "min_vacuum_size":
                    generated_list.append(round(random.uniform(1.0, 3.0), 3))
                    continue
            case "calculate_adsorption_energy":
                if tool_param_name == "adsorbate":
                    generated_list.append(ELEMENTS[random.randint(1, 118)])
                    continue

            # Compostion_Analysis
            case "mix_elements":
                compound = get_random_compound()
                element_count = Counter(get_random_compound()["compound"].elements)
                element_count_dict = dict(element_count)
                element_sum = sum(element_count_dict.values())
                
                composition_str = compound["MF"]
                element_ratios = {element:element_count_dict[element]/element_sum for element in element_count_dict}

                generated_list.append(composition_str)
                generated_list.append(element_ratios)
                generated_list.append("composition_{}.pkl".format(generate_random_string(10)))
                return generated_list
            
            # Crystal_Structure_Manipulation
            case "create_structure":
                lattice = [[round(random.uniform(-3.0, 3.0), 3) for _ in range(3)] for _ in range(3)]
                num = random.randint(2,7)
                species = [ELEMENTS[random.randint(1, 118)] for _ in range(num)]
                coords = [[round(random.uniform(-3.0, 3.0), 3) for _ in range(3)] for _ in range(num)]
                
                generated_list.append(lattice)
                generated_list.append(species)
                generated_list.append(coords)
                generated_list.append("structure_{}.pkl".format(generate_random_string(10)))
                return generated_list
            case "modify_site_occupancy":
                if tool_param_name == "site_index":
                    site_index = random.randint(1,10)
                    occupancy = round(random.uniform(0.0,1.0),3)

                    generated_list.append(site_index)
                    generated_list.append(occupancy)
                    generated_list.append("structure_{}.pkl".format(generate_random_string(10)))
                    return generated_list
            case "create_vacancy":
                if tool_param_name == "site_index":
                    site_index = random.randint(1,10)
                    generated_list.append(site_index)
                    continue
            case "create_interstitial":
                if tool_param_name == "coords":
                    generated_list.append([round(random.uniform(-3.0, 3.0), 3) for _ in range(3)])
                    continue

            # Defect_Analysis
            case "create_defect":
                if tool_param_name == "defect_type":
                    defect_type = random.choice(['vacancy', 'interstitial', 'substitution'])
                    defect_site  = tuple([round(random.uniform(0.0,1.0),3) for _ in range(random.randint(2,4))])
                    multiplicity = random.randint(1,5)
                    charge = random.randint(1,5)

                    generated_list.append(defect_type)
                    generated_list.append(defect_site)
                    generated_list.append(multiplicity)
                    generated_list.append(charge)
                    generated_list.append("defect_{}.pkl".format(generate_random_string(10)))
                    return generated_list
            case "create_defect_entry":
                if tool_param_name == "energy":
                    generated_list.append(round(random.uniform(1.0,3.0),3))
                    continue

            # Diffusion_Analysis
            case "create_migration_graph":
                if tool_param_name == "migrating_ion":
                    generated_list.append(ELEMENTS[random.randint(1, 118)])
                    continue

            # IO_Operations
            case "get_structure_by_material_id":
                if tool_param_name == "api_key":
                    generated_list.append(generate_random_string(10))
                    continue
                elif tool_param_name == "material_id":
                    generated_list.append(str(random.randint(1,1000)))
                    continue

            # Machine_Learning_Integration
            case "train_model":
                if tool_param_name == "targets":
                    generated_list.append([round(random.uniform(1.0,50.0),3) for _ in range(random.randint(2,5))])
                    continue
                if tool_param_name == "test_size":
                    generated_list.append(round(random.uniform(0.1,0.3),3))
                    continue
                if tool_param_name == "random_state":
                    generated_list.append(random.randint(0,2))
                    continue

            # Reactions_And_Batteries
            case "get_all_entries":
                if tool_param_name == "chemical_system":
                    generated_list.append([ELEMENTS[random.randint(1, 118)] for _ in range(random.randint(2,7))])
                    continue
            
            # Structure_Manipulation_and_Analysis
            case "create_structure_from_spacegroup":
                spacegroups = [
                    "P1", "P-1", "P2", "P21", "C2", "Pm", "Pc", "Cm", "Cc", 
                    "P2/m", "P21/m", "C2/m", "P2/c", "P21/c", "C2/c", 
                    "P222", "P2221", "P21212", "P212121", "C2221", "C222", "F222", "I222", "I212121", 
                    "Pmm2", "Pmc21", "Pcc2", "Pma2", "Pca21", "Pnc2", "Pmn21", "Pba2", "Pna21", "Pnn2", 
                    "Cmm2", "Cmc21", "Ccc2", "Amm2", "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2", "Ima2", 
                    "Pmmm", "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn", "Pbcm", "Pnnm", 
                    "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmce", "Cmmm", "Cccm", "Cmme", "Ccce", 
                    "Fmmm", "Fddd", "Immm", "Ibam", "Ibca", "Imma", 
                    "P4", "P41", "P42", "P43", "I4", "I41", 
                    "P-4", "I-4", 
                    "P4/m", "P42/m", "P4/n", "P42/n", "I4/m", "I41/a", 
                    "P422", "P4212", "P4122", "P41212", "P4222", "P42212", "P4322", "P43212", 
                    "I422", "I4122", 
                    "P4mm", "P4bm", "P42cm", "P42nm", "P4cc", "P4nc", "P42mc", "P42bc", "I4mm", "I4cm", "I41md", "I41cd", 
                    "P-42m", "P-42c", "P-42n", "P-42b", "I-42m", "I-42d", 
                    "P4/mmm", "P4/mcc", "P4/nbm", "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm", "P4/ncc", 
                    "I4/mmm", "I4/mcm", "I41/amd", "I41/acd", 
                    "P3", "P31", "P32", "R3", 
                    "P-3", "R-3", 
                    "P312", "P321", "P3112", "P3121", "P3212", "P3221", "R32", 
                    "P3m1", "P31m", "P3c1", "P31c", "R3m", "R3c", 
                    "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c", 
                    "P6", "P61", "P65", "P62", "P64", "P63", 
                    "P-6", 
                    "P6/m", "P63/m", 
                    "P622", "P6122", "P6522", "P6222", "P6422", "P6322", 
                    "P6mm", "P6cc", "P63cm", "P63mc", 
                    "P-6m2", "P-6c2", "P-62m", "P-62c", 
                    "P6/mmm", "P6/mcc", "P63/mcm", "P63/mmc", 
                    "P23", "F23", "I23", "P213", "I213", 
                    "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3", "Ia-3", 
                    "P432", "P4232", "F432", "F4132", "I432", "P4332", "P4132", "I4132", 
                    "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c", "Fd-3m", "Fd-3c", 
                    "Im-3m", "Ia-3d"
                ]
                spacegroup = random.choice(spacegroups)
                lattice = [[round(random.uniform(-3.0, 3.0), 3) for _ in range(3)] for _ in range(3)]
                num = random.randint(2,7)
                species = [ELEMENTS[random.randint(1, 118)] for _ in range(num)]
                coords = [[round(random.uniform(-3.0, 3.0), 3) for _ in range(3)] for _ in range(num)]
            
                generated_list.append(spacegroup)
                generated_list.append(lattice)
                generated_list.append(species)
                generated_list.append(coords)
                generated_list.append("structure_{}.pkl".format(generate_random_string(10)))
                return generated_list
            case "remove_sites":
                if tool_param_name == "indices":
                    generated_list.append([random.randint(1,100) for _ in range(random.randint(2,5))])
                    continue

            # Surface_Interface_Analysis
            case "generate_surface":
                if tool_param_name == "miller_index":
                    miller_index = [random.randint(1,10) for _ in range(random.randint(2,5))]
                    min_slab_size = round(random.uniform(1.0,5.0),3)
                    min_vacuum_size = round(random.uniform(1.0,5.0),3)
                    generated_list.append(miller_index)
                    generated_list.append(min_slab_size)
                    generated_list.append(min_vacuum_size)
                    return generated_list   
            case "model_interfacial_reactions":
                if tool_param_name == "open_elem":
                    generated_list.append(ELEMENTS[random.randint(1, 118)])
                    continue
                if tool_param_name == "temperature":
                    generated_list.append(round(random.uniform(20.0,500.0),3))
                    continue
            case "perform_adsorption_study":
                if tool_param_name == "adsorbate":
                    generated_list.append(ELEMENTS[random.randint(1, 118)])
                    continue

            case _:
                pass
        # 一般参数
        if "pickle_file_path" in tool_param_name:
            if "band_structure_pickle_file_path" in tool_param_name:
                generated_list.append("band_structure_{}.pkl".format(generate_random_string(10)))
            elif "structure_pickle_file_path" in tool_param_name:
                if "list" in tool_param_name:
                    generated_list.append(["structure_{}.pkl".format(generate_random_string(10)) for _ in range(random.randint(2,5))])
                else:
                    generated_list.append("structure_{}.pkl".format(generate_random_string(10)))
            elif "lattice_pickle_file_path" in tool_param_name:
                generated_list.append("lattice_{}.pkl".format(generate_random_string(10)))
            elif "entries_pickle_file_path" in tool_param_name:
                if "reactant" in tool_param_name:
                    generated_list.append("reactant_entries_{}.pkl".format(generate_random_string(10)))
                elif "product" in tool_param_name:
                    generated_list.append("product_entries_{}.pkl".format(generate_random_string(10)))
                else:
                    generated_list.append("entries_{}.pkl".format(generate_random_string(10)))
            elif "entry_pickle_file_path" in tool_param_name:
                if "defect" in tool_param_name:
                    generated_list.append("defect_entry_{}.pkl".format(generate_random_string(10)))
                else:
                    generated_list.append("entry_{}.pkl".format(generate_random_string(10)))
            elif "slab_pickle_file_path" in tool_param_name:
                generated_list.append("slab_{}.pkl".format(generate_random_string(10)))
            elif "composition_pickle_file_path" in tool_param_name:
                generated_list.append("composition_{}.pkl".format(generate_random_string(10)))
            elif "vacancy_pickle_file_path" in tool_param_name:
                generated_list.append("vacancy_{}.pkl".format(generate_random_string(10)))
            elif "defect_pickle_file_path" in tool_param_name:
                generated_list.append("defect_{}.pkl".format(generate_random_string(10)))
            elif "neb_analysis_pickle_file_path" in tool_param_name:
                generated_list.append("neb_analysis_{}.pkl".format(generate_random_string(10)))
            elif "migration_graph_pickle_file_path" in tool_param_name:
                generated_list.append("migration_graph_{}.pkl".format(generate_random_string(10)))
            elif "dos_pickle_file_path" in tool_param_name:
                generated_list.append("dos_{}.pkl".format(generate_random_string(10)))
            elif "workflow_pickle_file_path" in tool_param_name:
                generated_list.append("workflow_{}.pkl".format(generate_random_string(10)))
            elif "job_pickle_file_path" in tool_param_name:
                generated_list.append("job_{}.pkl".format(generate_random_string(10)))
            elif "cif_pickle_file_path" in tool_param_name:
                generated_list.append("cif_{}.pkl".format(generate_random_string(10)))
            elif "vasp_pickle_file_path" in tool_param_name:
                generated_list.append("vasp_{}.pkl".format(generate_random_string(10)))
            elif "feature_pickle_file_path" in tool_param_name:
                generated_list.append("feature_{}.pkl".format(generate_random_string(10)))
            elif "model_pickle_file_path" in tool_param_name:
                generated_list.append("model_{}.pkl".format(generate_random_string(10)))
            elif "phase_diagram_pickle_file_path" in tool_param_name:
                generated_list.append("phase_diagram_{}.pkl".format(generate_random_string(10)))
            elif "reaction_pickle_file_path" in tool_param_name:
                generated_list.append("reaction_{}.pkl".format(generate_random_string(10)))
            elif "pourbaix_diagram_pickle_file_path" in tool_param_name:
                generated_list.append("pourbaix_diagram_{}.pkl".format(generate_random_string(10)))
            else:
                generated_list.append("{}.pkl".format(generate_random_string(10)))
        elif "file_path" in tool_param_name:
            if "structure_file_path" in tool_param_name:
                generated_list.append("structure_{}.{}".format(generate_random_string(10),random.choice(["cif","vasp","nc","cssr","mson","yaml","xsf","xml","res","pwmat"])))
        elif "str" in tool_param_name:
            if "element_str" in tool_param_name:
                generated_list.append(ELEMENTS[random.randint(1, 118)])
            elif "composition_str" in tool_param_name:
                compound = get_random_compound()
                generated_list.append(compound["MF"])
                # if tool_param_name.startswith("compound_SMILES"):
                #     generated_list.append(compound["SMILES"])
                # elif tool_param_name.startswith("compound_CAS"):
                #     generated_list.append(compound["CAS"])
                # elif tool_param_name.startswith("compound_CID"):
                #     generated_list.append(compound["CID"])
                # else: # 默认为分子式
                #     generated_list.append(compound["MF"])
        elif "dir" in tool_param_name:
            if "neb_dir" in tool_param_name:
                generated_list.append(generate_random_filepath(base_dir="/tmp_neb",depth=random.randint(1,4),filename=""))
            else:
                generated_list.append(generate_random_filepath(depth=random.randint(1,4),filename=""))
        elif "list" in tool_param_name:
            if "number_list" in tool_param_name:
                generated_list.append([ random.randint(-100,100) for _ in range(random.randint(3,10))])
        else:
            print("出现未知参数！")
            print("Tool Name: ", tool_name)
            print("Tool Param Name: ", tool_param_name)
            print("Tool Param Type: ", tool_param_type)
    return generated_list



if __name__ == "__main__":

    # print(get_tool_param("analyze_combustion",{"percent_composition":{"type": "dict"}}))

    # print(get_tool_param("reaction_stoichiometry_amounts",{"reaction_formula":{"type": "str"}, "compound_number":{"type": "int"}}))

    # print(get_tool_param("electromagnetic_wave_by_rydberg_equation",{"element":{"type": "str"}, "n1":{"type": "int"}, "n2":{"type": "int"}}))

    compound = Compound.from_cid(8495)
    element_count = Counter(compound.elements)
    element_count_dict = dict(element_count)
    print(
    element_count_dict
    )