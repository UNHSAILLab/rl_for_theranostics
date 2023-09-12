import numpy as np



##### Units for the simulation
##  Volume  --> L
##  Time    --> min
##  Mass    --> gr
##  Moles   --> nmol


class Patient:

    def __init__(self):

        # Constants related to internalization and physiological clearance rates
        lambda_intern = 0.001  # Rate of internalization (l/min)
        lambda_phys = 7.23 * 1e-5  # Physiological clearance rate (l/min)

        # Constants related to receptor binding and unbinding
        k_on = 0.04 / 0.5  # Binding rate constant (/min/nmol)
        k_off = 0.04  # Unbinding rate constant (1/min)

        # Patient-specific parameters
        H = 0.1  # Hematocrit (fraction of blood that is red blood cells)
        BW = 80  # Body weight of the patient (kg)
        gender = "male"  # Gender of the patient
        BSA = 1.94  # Body Surface Area (m^2)

        # Volume of various compartments in the body
        V_body = BW * 1000  # Total body volume (L), based on the assumption that 1 g = 1 mL
        V_tu = 0.087  # Volume of the tumor (L)
        V_L = 1.811  # Volume of the liver (L)
        V_S = 0.198  # Volume of the spleen (L)
        V_K = 0.193  # Volume of the kidney (L)

        # Receptor densities in various organs (in nmol/L)
        R_tu_density = 15  # Receptor density in the tumor
        R_L_density = 1.4  # Receptor density in the liver
        R_S_density = 8.7  # Receptor density in the spleen
        R_K_density = 6.5  # Receptor density in the kidneys
        R_rest_density = 0.5  # Receptor density in the rest of the body

        # Rates of drug release from non-tumor and tumor tissues (in 1/min)
        lambda_rel_NT = 0.7 * 1e-4  # Rate of drug release from non-tumor tissue
        lambda_rel_TU = 1.1 * 1e-4  # Rate of drug release from tumor tissue

        # Tumor-specific parameters
        tumorType = "NET"  # Type of tumor, either "NET" (Neuroendocrine Tumor) or "MEN" (Multiple Endocrine Neoplasia)
        f_tu = 0.1  # Blood flow rate through the tumor (in L/min/g of tumor)
        k_pr = 4.7 * 1e-4  # Rate constant for protein binding (in 1/min)

        # Kidney function parameter
        GFR = 0.11  # Glomerular Filtration Rate (in L/min), a measure of how well the kidneys are filtering blood

        ## Permeability surface area product
        k_mu = 0.02                                ## L/min/kg | for muscle | --> please see the important note bellow. In a nutshell, /kg is the right unit here

        ## Important Note: In the table of Frenco 2021 paper, the units of k_i are in ml/min/g which is in fact
        ## the permeability surface area product per unit mass of organ. In order to get the permeability of the
        ## organ we need to multipy k_i at the mass of organ (1gr = 1ml). So the calculated PS value will be:
        ## PS = 1000 * k * V  (Note that V is in Litre). So PS will have the units of ml/min. Now to get the standard
        ## units of (L/min) we need to multiply it at 0.001. So PS with standard units (L/min) will be: PS = k*V
        # https://jnm.snmjournals.org/content/jnumed/suppl/2016/04/01/jnumed.115.164699.DC1/164699_Supplemental_Data.pdf
        # Calculate the volume of total body serum (V_p) based on gender
        # Different multipliers are used for males and females to account for physiological differences.
        # These multipliers (2.8 for males and 2.4 for females) are likely based on empirical or statistical models.
        if gender == "male":
            V_p = 2.8 * (1 - H) * BSA  # Calculate volume of total body serum for males
        else:
            V_p = 2.4 * (1 - H) * BSA  # Calculate volume of total body serum for females

        # Calculate the flow of total serum (F) using a multiplier of 1.23
        # This may be based on standard physiological data.
        F = 1.23 * V_p  # Flow of total serum (L/min)

        # Initialize tumor parameters based on the type of tumor (NET or MEN)
        # These values are specific to the type of tumor and are likely based on empirical data.
        if tumorType == "NET":
            v_tu_int = 0.3  # Initial volume of tumor interstitial fluid for NET
            v_tu_v = 0.1  # Initial volume of tumor vascular fluid for NET
            k_tu = 0.2  # Permeability surface area product for NET (L/min/Kg)
        else:
            v_tu_int = 0.23  # Initial volume of tumor interstitial fluid for MEN
            v_tu_v = 0.11  # Initial volume of tumor vascular fluid for MEN
            k_tu = 0.31  # Permeability surface area product for MEN (L/min/Kg)

        # Initialize the rate constants for the release of the drug from the tumor.
        # These are likely empirical constants.
        lambda_rel_tu = 1.5 * 1e-4  # Rate constant for drug release from tumor (1/min)
        lambda_intern_tu = 0.001  # Rate constant for internalization of drug in tumor (1/min)

        # Initialize the Tumor dictionary with various physiological and pharmacological parameters
        self.Tumor = {
            "name": "Tumor",  # Name of the compartment
            "F": f_tu * (1 - H) * V_tu,  # Flow rate of fluid into the tumor (L/min)
            "PS": k_tu * V_tu,  # Permeability surface area product for the tumor (L/min)
            "V_total": V_tu,  # Total volume of the tumor compartment (L)
            "V_v": v_tu_v * (1 - H) * V_tu,  # Volume of the vascular fluid in the tumor (L)
            "V_int": v_tu_int * V_tu,  # Volume of the interstitial fluid in the tumor (L)
            "k_on": k_on,  # Rate constant for drug binding (1/min)
            "k_off": k_off,  # Rate constant for drug unbinding (1/min)
            "lambda_intern": lambda_intern_tu,  # Internal rate constant for drug within tumor (1/min)
            "lambda_rel": lambda_rel_tu,  # Rate constant for drug release from the tumor (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (likely for drug metabolism or clearance) (1/min)
            "R0": R_tu_density * V_tu,  # Initial receptor density in the tumor (nmol)
            "K_on": 0,  # Placeholder for the effective rate constant for drug binding, to be calculated as R * k_on
        }

        # Initialize the Liver dictionary with various physiological and pharmacological parameters
        self.Liver = {
            "name": "Liver",  # Name of the compartment
            "F": 0.065 * F,  # Flow rate of fluid into the liver, 6.5% of the total body fluid flow (L/min)
            "PS": 100 * k_mu * V_L,  # Permeability surface area product for the liver (L/min)
            "V_total": V_L,  # Total volume of the liver compartment (L)
            "V_v": 0.085 * V_L,  # Volume of the vascular fluid in the liver, 8.5% of liver volume (L)
            "V_int": 0.2 * V_L,  # Volume of the interstitial fluid in the liver, 20% of liver volume (L)
            "k_on": k_on,  # Rate constant for drug binding (1/min)
            "k_off": k_off,  # Rate constant for drug unbinding (1/min)
            "lambda_intern": 1.7 * lambda_intern_tu,
            # Internal rate constant for drug within liver, 1.7 times that in the tumor (1/min)
            "lambda_rel": lambda_rel_NT,  # Rate constant for drug release from the liver (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (likely for drug metabolism or clearance) (1/min)
            "R0": R_L_density * V_L,  # Initial receptor density in the liver (nmol)
            "K_on": 0,  # Placeholder for the effective rate constant for drug binding, to be calculated as R * k_on
        }
        # Initialize the Spleen dictionary with various physiological and pharmacological parameters
        self.Spleen = {
            "name": "Spleen",  # Name of the compartment
            "F": 0.03 * F,  # Flow rate of fluid into the spleen, 3% of the total body fluid flow (L/min)
            "PS": 100 * k_mu * V_S,  # Permeability surface area product for the spleen (L/min)
            "V_total": V_S,  # Total volume of the spleen compartment (L)
            "V_v": 0.12 * V_S,  # Volume of the vascular fluid in the spleen, 12% of spleen volume (L)
            "V_int": 0.2 * V_S,  # Volume of the interstitial fluid in the spleen, 20% of spleen volume (L)
            "k_on": k_on,  # Rate constant for drug binding (1/min)
            "k_off": k_off,  # Rate constant for drug unbinding (1/min)
            "lambda_intern": 1.7 * lambda_intern_tu,
            # Internal rate constant for drug within spleen, 1.7 times that in the tumor (1/min)
            "lambda_rel": lambda_rel_NT,  # Rate constant for drug release from the spleen (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (likely for drug metabolism or clearance) (1/min)
            "R0": R_S_density * V_S,  # Initial receptor density in the spleen (nmol)
            "K_on": 0,  # Placeholder for the effective rate constant for drug binding, to be calculated as R * k_on
        }

        # Initialize the Kidney dictionary with various physiological and pharmacological parameters
        self.Kidney = {
            "name": "Kidney",  # Name of the compartment
            "F": 0.19 * F,  # Flow rate of fluid into the kidney, 19% of the total body fluid flow (L/min)
            "V_total": V_K,  # Total volume of the kidney compartment (L)
            "V_v": 0.055 * V_K,  # Volume of the vascular fluid in the kidney, 5.5% of kidney volume (L)
            "V_int": 0.15 * V_K,  # Volume of the interstitial fluid in the kidney, 15% of kidney volume (L)
            "V_intra": -1,  # Placeholder for the intracellular volume in the kidney, will be calculated later
            "k_on": k_on,  # Rate constant for drug binding (1/min)
            "k_off": k_off,  # Rate constant for drug unbinding (1/min)
            "lambda_intern": 1.7 * lambda_intern_tu,
            # Internal rate constant for drug within kidney, 1.7 times that in the tumor (1/min)
            "lambda_rel": lambda_rel_NT,  # Rate constant for drug release from the kidney (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (likely for drug metabolism or clearance) (1/min)
            "R0": R_K_density * V_S,  # Initial receptor density in the kidney (nmol)
            "K_on": 0,  # Placeholder for the effective rate constant for drug binding, to be calculated as R * k_on
            "phi": 1.1,  # Filtration constant
            "GFR": GFR,  # Glomerular Filtration Rate (in L/min)
            "f_exc": 0.98,  # Fraction of fluid that is excreted
            "F_fil": -1,  # Placeholder for the filtration rate, will be calculated as GFR * phi
            "F_R": -1,  # Placeholder for the return flow rate, will be calculated as F_fil * (1-f_exc)
        }

        # Calculating and updating the missing values in the Kidney dictionary
        self.Kidney["F_fil"] = self.Kidney["GFR"] * self.Kidney["phi"]
        self.Kidney["F_R"] = self.Kidney["F_fil"] * (1 - self.Kidney["f_exc"])
        self.Kidney["V_intra"] = (self.Kidney["V_total"] - (self.Kidney["V_int"] + self.Kidney["V_v"])) * 2 / 3

        # Initialize the Prostate (if male) or Uterus (if female) dictionary with relevant physiological and pharmacological parameters
        if gender == "male":  # Prostate Initialization
            V_total_prostateUterus = 0.016 * BW / 71  # Total volume of the prostate in L, scaled with body weight
            V_v_prostateUterus = 0.04 * (
                        1 - H) * V_total_prostateUterus  # Volume of vascular fluid in the prostate in L
            V_int_prostateUterus = 0.25 * V_total_prostateUterus  # Volume of interstitial fluid in the prostate in L
            F_prostateUterus = 0.18 * (1 - H) * V_total_prostateUterus  # Flow rate of fluid into the prostate in L/min
            k_prostateUterus = 0.1  # Permeability-surface area product constant for prostate
            R_density_prostateUterus = R_K_density * 0.26  # Initial receptor density in the prostate (nmol/L)

        if gender == "female":  # Uterus Initialization
            V_total_prostateUterus = 0.08 * BW / 71  # Total volume of the uterus in L, scaled with body weight
            V_v_prostateUterus = 0.07 * (1 - H) * V_total_prostateUterus  # Volume of vascular fluid in the uterus in L
            V_int_prostateUterus = 0.5 * V_total_prostateUterus  # Volume of interstitial fluid in the uterus in L
            F_prostateUterus = 1 * (1 - H) * V_total_prostateUterus  # Flow rate of fluid into the uterus in L/min
            k_prostateUterus = 0.2  # Permeability-surface area product constant for uterus
            R_density_prostateUterus = R_K_density * 0.092  # Initial receptor density in the uterus (nmol/L)

        # Initialize the ProstateUterus dictionary
        self.ProstateUterus = {
            "name": "ProstateUterus",  # Name of the compartment
            "F": F_prostateUterus,  # Flow rate of fluid into the compartment in L/min
            "PS": k_prostateUterus * V_total_prostateUterus,
            # Permeability-surface area product for the compartment in L/min
            "V_total": V_total_prostateUterus,  # Total volume of the compartment in L
            "V_v": V_v_prostateUterus,  # Volume of vascular fluid in the compartment in L
            "V_int": V_int_prostateUterus,  # Volume of interstitial fluid in the compartment in L
            "k_on": k_on,  # Rate constant for drug binding (1/min)
            "k_off": k_off,  # Rate constant for drug unbinding (1/min)
            "lambda_intern": 1.7 * lambda_intern_tu,
            # Internal rate constant for drug within the compartment, 1.7 times that in the tumor (1/min)
            "lambda_rel": lambda_rel_NT,  # Rate constant for drug release from the compartment (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (likely for drug metabolism or clearance) (1/min)
            "R0": R_density_prostateUterus * V_total_prostateUterus,
            # Initial receptor density in the compartment (nmol)
            "K_on": 0,  # Placeholder for the effective rate constant for drug binding, to be calculated as R * k_on
        }

        # Initialize the Adrenal Glands with relevant physiological and pharmacological parameters

        # Total volume of the adrenal glands in L, scaled with body weight
        V_total_adrenal = 0.014 * BW / 71

        # Volume of vascular fluid in the adrenal glands in L
        V_v_adrenal = 0.03 * (1 - H) * V_total_adrenal

        # Volume of interstitial fluid in the adrenal glands in L
        V_int_adrenal = 0.24 * V_total_adrenal

        # Arbitrary factor for flow rate into the adrenal glands
        f_adrenal = 6

        # Flow rate of fluid into the adrenal glands in L/min
        F_adrenal = f_adrenal * (1 - H) * V_total_adrenal

        # Permeability-surface area product constant for adrenal glands
        k_adrenal = k_mu * 100

        # Initial receptor density in the adrenal glands (nmol/L)
        R_density_adrenal = R_K_density * 1.65

        # Initialize the Adrenals dictionary
        self.Adrenals = {
            "name": "Adrenals",  # Name of the compartment
            "F": F_adrenal,  # Flow rate of fluid into the compartment in L/min
            "PS": k_adrenal * V_total_adrenal,  # Permeability-surface area product for the compartment in L/min
            "V_total": V_total_adrenal,  # Total volume of the compartment in L
            "V_v": V_v_adrenal,  # Volume of vascular fluid in the compartment in L
            "V_int": V_int_adrenal,  # Volume of interstitial fluid in the compartment in L
            "k_on": k_on,  # Rate constant for drug binding (1/min)
            "k_off": k_off,  # Rate constant for drug unbinding (1/min)
            "lambda_intern": 1.7 * lambda_intern_tu,
            # Internal rate constant for drug within the compartment, 1.7 times that in the tumor (1/min)
            "lambda_rel": lambda_rel_NT,  # Rate constant for drug release from the compartment (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (likely for drug metabolism or clearance) (1/min)
            "R0": R_density_adrenal * V_total_adrenal,  # Initial receptor density in the compartment (nmol)
            "K_on": 0,  # Placeholder for the effective rate constant for drug binding, to be calculated as R * k_on
        }

        # GI (Gastrointestinal Tract) Initialization
        # ------------------------------------------

        # Calculate the total volume of the GI tract in Liters (L)
        # Volumes of various components of the GI are summed and normalized by body weight.
        V_total_GI = (0.385 + 0.548 + 0.104 + 0.15) * BW / 71  # L

        # Calculate the volume of the vascular compartment in the GI tract in Liters (L)
        V_v_GI = 0.076 * V_p  # L

        # Define the ratio of the interstitial volume to the vascular volume in the GI tract
        alpha_GI = 8.8  # Interstitial to vascular ratio

        # Calculate the volume of the interstitial compartment in the GI tract in Liters (L)
        V_int_GI = alpha_GI * V_v_GI  # L

        # Calculate the flow rate of fluid into the GI tract in Liters per minute (L/min)
        F_GI = 0.16 * F

        # Define the permeability-surface area product constant for the GI tract
        k_GI = k_mu

        # Define the initial receptor density in the GI tract in nanomoles per Liter (nmol/L)
        R_density_GI = R_K_density * 0.16

        # Initialize the GI dictionary to store various parameters
        self.GI = {
            "name": "GI",  # Name of the organ
            "F": F_GI,  # Flow rate (L/min)
            "PS": k_GI * V_total_GI,  # Permeability-surface area product (L/min)
            "V_total": V_total_GI,  # Total volume (L)
            "V_v": V_v_GI,  # Vascular volume (L)
            "V_int": V_int_GI,  # Interstitial volume (L)
            "k_on": k_on,  # Rate of ligand binding (1/min)
            "k_off": k_off,  # Rate of ligand unbinding (1/min)
            "lambda_intern": 1.7 * lambda_intern_tu,  # Rate of internalization (1/min)
            "lambda_rel": lambda_rel_NT,  # Rate of release (1/min)
            "lambda_phys": lambda_phys,  # Physiological rate constant (1/min)
            "R0": R_density_GI * V_total_GI,  # Initial number of receptors
            "K_on": 0,  # Product of receptor density and rate of ligand binding (R*k_on), initially set to zero
        }

        ##
        ####
        ######
        ######## RedMarrow initialization
        V_total_RedMarrow = 1.1 * BW / 71  # L
        V_v_RedMarrow = 0.04 * V_p  # L
        alpha_RedMarrow = 3.7  ## interestitial to vascular ratio
        V_int_RedMarrow = alpha_RedMarrow * V_v_RedMarrow  # L
        F_RedMarrow = 0.03 * F
        k_RedMarrow = k_mu*100
        R_density_RedMarrow = R_K_density * 0.028
        self.RedMarrow = {
            "name": "RedMarrow",
            "F": F_RedMarrow,
            "PS": k_RedMarrow * V_total_RedMarrow,
            "V_total": V_total_RedMarrow,
            "V_v": V_v_RedMarrow,
            "V_int": V_int_RedMarrow,
            "k_on": k_on,
            "k_off": k_off,
            "lambda_intern": 1.7 * lambda_intern_tu,
            "lambda_rel": lambda_rel_NT,
            "lambda_phys": lambda_phys,
            "R0": R_density_RedMarrow * V_total_RedMarrow,
            "K_on": 0,  ## R*k_on
        }

        ##
        ####
        ######
        ######## Muscle initialization
        V_total_Muscle = 30.078 * BW / 71  # L
        V_v_Muscle = 0.14 * V_p # L
        alpha_Muscle = 5.9  ## interestitial to vascular ratio
        V_int_Muscle = alpha_Muscle * V_v_Muscle  # L
        F_Muscle = 0.17 * F
        k_Muscle = k_mu
        R_density_Muscle = R_K_density * 0.0056
        self.Muscle = {
            "name": "Muscle",
            "F": F_Muscle,
            "PS": k_Muscle * V_total_Muscle,
            "V_total": V_total_Muscle,
            "V_v": V_v_Muscle,
            "V_int": V_int_Muscle,
            "k_on": k_on,
            "k_off": k_off,
            "lambda_intern": 1.7 * lambda_intern_tu,
            "lambda_rel": lambda_rel_NT,
            "lambda_phys": lambda_phys,
            "R0": R_density_Muscle * V_total_Muscle,
            "K_on": 0,  ## R*k_on
        }

        ##
        ####
        ######
        ######## Lungs initialization
        V_total_Lungs = 1 * BW / 71  # L
        V_v_Lungs = 0.105 * V_p  # L
        alpha_Lungs = 5.5  ## interestitial to vascular ratio
        V_int_Lungs = V_v_Lungs * alpha_Lungs  # L
        F_Lungs = F
        k_Lungs = 100 * k_mu
        self.Lungs = {
            "name": "Lungs",
            "F": F_Lungs,  # The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "PS": k_Lungs*V_total_Lungs,
            "V_total": V_total_Lungs,
            "V_v": V_v_Lungs,
            "V_int": V_int_Lungs,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Skin initialization
        V_total_Skin = 3.408 * BW / 71  # L
        V_v_Skin = 0.03 * V_p  # L
        alpha_Skin = 8.9  ## interestitial to vascular ratio
        V_int_Skin = alpha_Skin * V_v_Skin  # L
        F_Skin = 0.05 * F
        k_Skin = k_mu
        self.Skin = {
            "name": "Skin",
            "F": F_Skin,  # The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "PS": k_Skin * V_total_Skin,
            "V_total": V_total_Skin,
            "V_v": V_v_Skin,
            "V_int": V_int_Skin,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Adipose initialization
        V_total_Adipose = 13.465 * BW / 71  # L
        V_v_Adipose = 0.05 * V_p  # L
        alpha_Adipose = 15.5  ## interestitial to vascular ratio
        V_int_Adipose = alpha_Adipose * V_v_Adipose  # L
        F_Adipose = 0.05 * F
        k_Adipose = k_mu
        self.Adipose = {
            "name": "Adipose",
            "F": F_Adipose,  # The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "PS": k_Adipose * V_total_Adipose,
            "V_total": V_total_Adipose,
            "V_v": V_v_Adipose,
            "V_int": V_int_Adipose,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Brain initialization
        V_total_Brain = 1.45 * BW / 71  # L
        V_v_Brain = 0.012 * V_p  # L
        alpha_Brain = 1  ## interestitial to vascular ratio
        ## Note that PS is zero for brain. To avoid division by zero I
        ## have intentionally put alpha_Brain to be 1 (so V_int is not zero)
        V_int_Brain = alpha_Brain * V_v_Brain  # L
        F_Brain = 0.04 * F
        k_Brain = 0
        self.Brain = {
            "name": "Brain",
            "F": F_Brain,  # The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "PS": k_Brain * V_total_Brain,
            "V_total": V_total_Brain,
            "V_v": V_v_Brain,
            "V_int": V_int_Brain,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Heart initialization
        V_total_Heart = 0.341 * BW / 71  # L
        V_v_Heart = 0.01 * V_p  # L
        alpha_Heart = 3.7  ## interestitial to vascular ratio
        V_int_Heart = alpha_Heart * V_v_Heart  # L
        F_Heart = 0.04 * F
        k_Heart = k_mu
        self.Heart = {
            "name": "Heart",
            "F": F_Heart,  # The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "PS": k_Heart * V_total_Heart,
            "V_total": V_total_Heart,
            "V_v": V_v_Heart,
            "V_int": V_int_Heart,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Bone initialization
        V_total_Bone = 10.165 * BW / 71 - self.RedMarrow['V_total']  # L
        V_v_Bone = 0.07 * V_p - self.RedMarrow["V_v"]  # L
        alpha_Bone = 9.3  ## interestitial to vascular ratio
        V_int_Bone = alpha_Bone * V_v_Bone  # L
        F_Bone = 0.05 * F
        k_Bone = k_mu
        self.Bone = {
            "name": "Bone",
            "F": F_Bone,  # The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "PS": k_Bone * V_total_Bone,
            "V_total": V_total_Bone,
            "V_v": V_v_Bone,
            "V_int": V_int_Bone,
            "lambda_phys": lambda_phys,
        }


        ##
        ####
        ######
        ######## BloodProtein initialization
        self.BloodProtein = {
            "name": "BloodProtein",
            "k_pr": 5e-4,
            "lambda_phys": lambda_phys
        }

        ##
        ####
        ######
        ######## Vein initialization
        V_total_Vein = (0.18 + 0.045) * V_p
        F_Vein = F
        self.Vein = {
            "name": "Vein",
            "F": F_Vein,  ## The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all
            ## flows
            "V_v": V_total_Vein,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Art initialization
        V_total_Art = (0.06 + 0.045) * V_p  # L
        F_Art = F
        self.Art = {
            "name": "Art",
            "F": F_Art,  ## The real value will be calculated once we have the
            ## the flow of each organ. The Vein and Art flow is the sum of all flows
            "V_v": V_total_Art,
            "lambda_phys": lambda_phys,
        }

        ##
        ####
        ######
        ######## Rest initialization
        self.Rest = {
            "name": "Rest",
            "F": F_Muscle,
            "PS": k_Muscle * V_total_Muscle,
            "V_total": V_v_Muscle,
            "V_v": V_v_Muscle,
            "V_int": V_int_Muscle,
            "k_on": k_on,
            "k_off": k_off,
            "lambda_intern": 1.7 * lambda_intern_tu,
            "lambda_rel": lambda_rel_NT,
            "lambda_phys": lambda_phys,
            "R0": R_density_Muscle * V_total_Muscle,
            "K_on": 0,  ## R*k_on
        }




        self.receptorPositiveList = [self.Tumor, self.Liver, self.Spleen, self.RedMarrow, self.GI,
                                     self.Muscle, self.ProstateUterus, self.Adrenals]
        self.KidneyList = [self.Kidney]
        self.LungsList = [self.Lungs]
        self.receptorNegativeList = [self.Skin, self.Adipose, self.Brain, self.Heart, self.Bone]
        self.BloodProteinList = [self.BloodProtein]
        self.ArtVeinList = [self.Art, self.Vein]

        self.calculateK_on()    ## Will calculate K_on = k_on * R0 for RecPos and Kidney

        self.Organs = {
            "ArtVein": self.ArtVeinList,
            "Lungs": self.LungsList,
            "RecNeg": self.receptorNegativeList,
            "BloodProtein": self.BloodProteinList,
            "RecPos": self.receptorPositiveList,
            "Kidney": self.KidneyList,
        }

        self.addRestOrgan(BW, F, V_p, k_mu, R_rest_density, k_on, k_off, lambda_intern_tu, lambda_rel_NT, lambda_phys)


    def calculateK_on(self):
        # for elem in self.receptorPositiveList:
        #     elem["K_on"] = elem["k_on"] * elem["R0"]
        # for elem in self.KidneyList:
        #     elem["K_on"] = elem["k_on"] * elem["R0"]
        #
        for elem in self.receptorPositiveList:
            elem["K_on"] = 0
        for elem in self.KidneyList:
            elem["K_on"] = 0



    def addRestOrgan(self, BW, F, V_p, k_mu, R_rest_density, k_on, k_off, lambda_intern_tu, lambda_rel_NT, lambda_phys):
        typeList = ["RecNeg", "RecPos", "Kidney"]
        F_allOrgans = 0
        V_total_allOrgans = 0
        V_v_allOrgans = 0
        for type in typeList:
            organsList = self.Organs[type]
            for organ in organsList:
                F_allOrgans += organ["F"]
                V_total_allOrgans += organ["V_total"]
                V_v_allOrgans += organ["V_v"]



        V_total_rest = BW - V_total_allOrgans   ## 1kg = 1lit, ##L
        F_rest = F - F_allOrgans
        V_v_rest = V_p - V_v_allOrgans

        alpha_rest = 3.7
        V_int_rest = alpha_rest * V_v_rest
        k_rest = k_mu
        self.Rest = {
            "name": "Rest",
            "F": F_rest,
            "PS": k_rest * V_total_rest,
            "V_total": V_total_rest,
            "V_v": V_v_rest,
            "V_int": V_int_rest,
            "k_on": k_on,
            "k_off": k_off,
            "lambda_intern": 1.7 * lambda_intern_tu,
            "lambda_rel": lambda_rel_NT,
            "lambda_phys": lambda_phys,
            "R0": R_rest_density * V_total_rest,
            "K_on": 0,  ## R*k_on
        }
        self.receptorPositiveList.append(self.Rest)




    def calculateTotalF(self):
        F = 0.0
        for elem in self.receptorNegativeList:
            F += elem['F']
        for elem in self.receptorPositiveList:
            F += elem['F']
        for elem in self.KidneyList:
            F += elem["F"]

        self.Art['F'] = F
        self.Vein['F'] = F
        self.Lungs["F"] = F
