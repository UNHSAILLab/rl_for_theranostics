import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


class Solver:
    def __init__(self, encoder):
        # Copy the system matrix and initial state vector from the encoder object
        self.SystemMat = encoder.SystemMat.copy()
        self.BigVect = encoder.BigVect.copy()

        # Retrieve the organ information and injection profile from the encoder object
        self.organsObj = encoder.organsObj
        self.injectionProfile = self.organsObj.therapy.injectionProfile

        # Initialize simulation configuration and injection settings
        self.setSimConf()
        self.setInjection()

        # Initialize a matrix to store the state of the system at each time step
        # The number of rows equals the length of the state vector
        # The number of columns equals the number of time steps
        self.BigVectList = np.zeros((self.BigVect.shape[0], self.tList.shape[0]))

        # Perform the first injection at t = 0
        self.inject(0)

        # Store the initial state in the first column of BigVectList
        self.BigVectList[:, 0] = self.BigVect.copy()

    def setInjection(self):
        # Initialize total amount of hot and cold injections to zero
        self.totalHot = 0.0
        self.totalCold = 0.0

        # Retrieve vein-related information from the organ dictionary
        Vein_dict = self.organsObj.organsDict["ArtVein"]["Vein"]

        # Get the indices of the 'cold' and 'hot' peptides in the big state vector
        self.Vein_index_cold = Vein_dict["stencil"]["base"] + Vein_dict["bigVectMap"]["P"]
        self.Vein_index_hot = Vein_dict["stencil"]["base"] + Vein_dict["bigVectMap"]["P*"]

        # If the injection type is constant, calculate the amount to inject at each time step
        if self.injectionProfile["type"] == "constant":
            t0 = self.injectionProfile["t0"]
            tf = self.injectionProfile["tf"]
            self.toInjectAtEachTimeStep_Hot = self.injectionProfile["totalAmountHot"] * self.h / (tf - t0)
            self.toInjectAtEachTimeStep_Cold = self.injectionProfile["totalAmountCold"] * self.h / (tf - t0)

        # If the injection type is a single bolus, initialize a flag to ensure only one bolus is injected
        if self.injectionProfile["type"] == "bolus":
            self.singleBolusFlag = 0

        # If the injection type is a train of boluses, initialize variables for multiple injections
        if self.injectionProfile["type"] == "bolusTrain":
            self.multiBolusFlag = 0  # Initialize a flag to count the number of boluses
            self.eachBolusCold = self.injectionProfile["totalAmountCold"] / self.injectionProfile["N"]
            self.eachBolusHot = self.injectionProfile["totalAmountHot"] / self.injectionProfile["N"]

    def setSimConf(self):
        # Initial time for the simulation
        t_0 = 0

        # Final time for the simulation;
        # Note: 75 works best when l=16. Increase l by one if you change t_f to 150.
        t_f = 75

        # Number of time steps, represented as 2 to the power of l
        l = 16

        # Generate a list of time points for the simulation
        self.tList = np.linspace(t_0, t_f, 2 ** (l))

        # Calculate the time step (difference between two consecutive time points)
        self.h = self.tList[1] - self.tList[0]

    def F(self, X, i):
        # Call the getSystemMat method to obtain the system matrix at the given time index 'i'
        # and state vector 'X'
        SystemMat = self.getSystemMat(X, i)

        # Perform matrix-vector multiplication between the system matrix and the state vector
        # to get the rate of change of the system
        return np.matmul(SystemMat, X)

    def getSystemMat(self, X, i):
        # Create a copy of the current system matrix
        SystemMat = self.SystemMat.copy()

        # TODO: Kidney needs to be added as well
        # Iterate over keys in the RecPos (Receptor Positive) organs dictionary (Tumor, Liver, Kidney, etc.)
        for key in self.organsObj.organsDict["RecPos"].keys():
            # Get the organ object using the key
            organ = self.organsObj.organsDict["RecPos"][key]

            # Get the previous state vector for time i-1
            BigVector_pre = self.BigVectList[:, i - 1].copy()

            # Calculate the index of the receptor positive (RP) and RP* (unlabeled) compartments in the BigVector
            RP_index = organ["stencil"]["base"] + organ["bigVectMap"]["RP"]
            RP_unlabeled_index = organ["stencil"]["base"] + organ["bigVectMap"]["RP*"]

            # Calculate K_on (binding rate constant) at time i-1 and i
            K_on_pre = (organ["R0"] - (BigVector_pre[RP_index] + BigVector_pre[RP_unlabeled_index])) * organ["k_on"]
            K_on = (organ["R0"] - (X[RP_index] + X[RP_unlabeled_index])) * organ["k_on"]

            # Update the system matrix based on the change in K_on
            for elem in organ["sysMatMap"]["K_on"]:
                sign = elem[-1]
                pos = np.array(elem[:-1]) + organ["stencil"]["base"]
                SystemMat[pos[0], pos[1]] += sign * (K_on - K_on_pre)
        # Return the updated system matrix
        return SystemMat

    def solve(self):
        # Loop through the time list, starting from the second element
        for j, t in enumerate(self.tList[1:]):
            # Calculate the true index (enumerate starts at 0)
            i = j + 1

            # Inject substances at time t
            self.inject(t)

            # Calculate the function values for Runge-Kutta 4th order method
            f0 = self.F(self.BigVect, i)
            f1 = self.F(self.BigVect + f0 * self.h / 2, i)
            f2 = self.F(self.BigVect + f1 * self.h / 2, i)
            f3 = self.F(self.BigVect + f2 * self.h, i)

            # Update the state vector using Runge-Kutta 4th order method
            self.BigVect = self.BigVect + self.h / 6 * (f0 + 2 * f1 + 2 * f2 + f3)

            # Update the system matrix based on the new state
            self.SystemMat = self.getSystemMat(self.BigVect, i)

            # Store the new state vector into the list of all state vectors
            self.BigVectList[:, i] = self.BigVect.copy()

            # Print the current iteration if it's a multiple of 1000 (for debugging or monitoring)
            if j % 1000 == 0:
                print(j)

        ####  Debugging
        # debugList = [9,11,13,15]
        # peptide = self.BigVectList[0,:]*0
        # for i in debugList:
        #     peptide += self.BigVectList[i,:]
        print("Hello")

    def inject(self, t):
        # Check if the injection type is constant
        if self.injectionProfile["type"] == "constant":
            # Check if the current time is less than the final time for injection
            if t < self.injectionProfile["tf"]:
                # Add a certain amount of cold and hot substances to the vein at each time step
                self.BigVect[self.Vein_index_cold] += self.toInjectAtEachTimeStep_Cold
                self.BigVect[self.Vein_index_hot] += self.toInjectAtEachTimeStep_Hot

                # Keep track of the total amount of hot and cold substances injected
                self.totalHot += self.toInjectAtEachTimeStep_Hot
                self.totalCold += self.toInjectAtEachTimeStep_Cold

        # Check if the injection type is a single bolus
        if self.injectionProfile["type"] == "bolus":
            # Check if the current time is greater than or equal to the initial time for bolus injection
            if t >= self.injectionProfile["t0"] and self.singleBolusFlag == 0:
                # Inject the total amount of cold and hot substances once (bolus)
                self.BigVect[self.Vein_index_cold] += self.injectionProfile["totalAmountCold"]
                self.BigVect[self.Vein_index_hot] += self.injectionProfile["totalAmountHot"]

                # Update the total amounts
                self.totalCold += self.injectionProfile["totalAmountCold"]
                self.totalHot += self.injectionProfile["totalAmountHot"]

                # Set the flag to 1 to avoid another bolus injection
                self.singleBolusFlag = 1

        # Check if the injection type is a train of boluses
        if self.injectionProfile["type"] == "bolusTrain":
            # Check if the number of boluses hasn't reached the specified number 'N'
            if self.multiBolusFlag != self.injectionProfile["N"]:
                # Check if the current time matches the time for the next bolus in the train
                if t >= self.injectionProfile["t"][self.multiBolusFlag]:
                    # Inject a fraction of the total amount (one bolus in the train)
                    self.BigVect[self.Vein_index_cold] += self.eachBolusCold
                    self.BigVect[self.Vein_index_hot] += self.eachBolusHot

                    # Update the total amounts
                    self.totalCold += self.eachBolusCold
                    self.totalHot += self.eachBolusHot

                    # Increment the flag to go to the next bolus
                    self.multiBolusFlag += 1











