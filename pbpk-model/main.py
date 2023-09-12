# Importing required modules and classes

# Custom class for data processing
from DataProcessing import DataProcessing
# Custom classes for encoding and solving
from Encoder import Encoder
# Custom classes for patient and therapy models
from Patient import Patient
from StiffSolver import StiffSolver
from Therapy import Therapy

# TODO: Add Albumin compartment to the model and evaluate its effect in the dynamics

# List of injection profiles to be tested
# Uncomment the line below if you want to test multiple injection profiles
# injection_profiles = ["constantInjection60", "constantInjection120", "constantInjection180", ...]
# For this example, only 'bolusInjection' is being tested

#A "bolus injection" is a method of drug administration in which a single, concentrated dose of a medication is given intravenously,
# usually all at once over a short period. This is in contrast to other methods like continuous infusion,
# where the drug is administered over an extended period.
injection_profiles = ["bolusInjection"]

# Loop through each injection profile
for index, profile in enumerate(injection_profiles):
    # Initialize patient and therapy models
    patient_model = Patient()
    therapy_model = Therapy(index)  # 'index' is used as an argument here; adjust as needed

    # Initialize the encoder
    encoder_model = Encoder(patient_model, therapy_model)

    # Initialize the solver
    # Uncomment the solver you want to use
    # solver_model = Solver(encoder_model)
    # solver_model = AdaptiveSolver(encoder_model)
    solver_model = StiffSolver(encoder_model)

    # Solve the model
    solver_model.solve()

    # Process and plot the data
    data_processor = DataProcessing(patient_model, therapy_model, encoder_model, solver_model)
    data_processor.plotter()

    # Save data if needed
    # np.save("data/" + profile, data_processor.results.BigVectList)

# Print a message to indicate script has finished running
print("Finished running the script.........")
