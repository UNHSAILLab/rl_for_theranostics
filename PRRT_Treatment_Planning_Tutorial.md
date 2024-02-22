
# PRRT Treatment Planning Tutorial

## Overview of Peptide Receptor Radionuclide Therapy (PRRT)

/**
 * Peptide Receptor Radionuclide Therapy (PRRT) is a molecular therapy used for treating neuroendocrine tumors.
 * It involves targeting tumor cells with peptides that bind to somatostatin receptors, which are often overexpressed in these tumors.
 */
Peptide Receptor Radionuclide Therapy (PRRT) is a molecular therapy used for treating neuroendocrine tumors. It involves targeting tumor cells with peptides that bind to somatostatin receptors, which are often overexpressed in these tumors.

## Radiobiological Concepts

### Linear-Quadratic (LQ) Model

The LQ model describes the response of cells to radiation based on two types of damage:

- **Linear Damage (α)**: Proportional to the dose.
- **Quadratic Damage (β)**: Proportional to the square of the dose.

The cell survival fraction (SF) can be calculated as:

```math
SF = e^{-lpha D - eta D^2}
```

### Biologically Effective Dose (BED)

BED accounts for different radiation schedules, normalizing to a biologically equivalent dose.

```math
BED = D (1 + rac{d}{lpha/eta})
```

### Overall Biologically Effective Dose (oBED)

For multiple lesions:

```math
oBED = rac{-1}{lpha} \cdot \ln\left( rac{\sum_{i=1}^{N} m_i \cdot e^{-lpha \cdot BED_i}}{\sum_{i=1}^{N} m_i} 
ight)
```

### Tumor Control Probability (TCP)

Probability of tumor eradication:

```math
TCP = e^{-n_0 \cdot SF}
```

### Overall Tumor Control Probability (oTCP)

For multiple treatment cycles:

```math
oTCP = e^{-n_0 \cdot (oSF_c)^{N_c}}
```

## Two-Phase Optimization Process for PRRT

### First Phase

- Establish multiple peptide amounts.
- Consider cell repair rates for tumors and OARs.
- Determine maximum tolerated BEDs for OARs.
- Use the LQ model to understand tissue response.
- Define the maximum achievable molar activity.
- Use individualized PBPK models for biokinetics.
- Calculate G factors and TIACs for each peptide amount.
- Determine LAACs for OARs and the radiopharmaceutical.
- Form a CLAAC combining OAR and radiopharmaceutical constraints.

### Second Phase

- Calculate tumor oBED for each CLAAC point.
- Identify the optimal plan with the highest tumor oBED.
- Select the optimal peptide molar amount and activity to be administered.

## Goals of Optimization

- Ensure OAR safety.
- Maximize tumor oBED.
- Confirm feasibility based on radiopharmaceutical synthesis.

---

This tutorial provides a conceptual framework for PRRT treatment planning. For actual treatment planning, collaboration with medical physicists, oncologists, and radiopharmacists is essential.

