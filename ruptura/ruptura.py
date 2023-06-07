import _ruptura
import numpy as np
from typing import Union, Literal
import matplotlib.pyplot as plt

isothermTypes = {
    "Langmuir": 0,
    "Anti-Langmuir": 1,
    "BET": 2,
    "Henry": 3,
    "Freundlich": 4,
    "Sips": 5,
    "Langmuir-Freundlich": 6,
    "Redlich-Peterson": 7,
    "Toth": 8,
    "Unilan": 9,
    "OBrien-Myers": 10,
    "Quadratic": 11,
    "Temkin": 12,
    "BingelWalton": 13,
}

isothermLabels = {
    "Langmuir": ["q_sat", "b"],
    "Anti-Langmuir": ["a", "b"],
    "BET": ["q_sat", "b", "c"],
    "Henry": ["a"],
    "Freundlich": ["a", "nu"],
    "Sips": ["q_sat", "b", "nu"],
    "Langmuir-Freundlich": ["q_sat", "b", "nu"],
    "Redlich-Peterson": ["a", "b", "nu"],
    "Toth": ["q_sat", "b", "nu"],
    "Unilan": ["q_sat", "b", "eta"],
    "OBrien-Myers": ["q_sat", "b", "sigma"],
    "Quadratic": ["q_sat", "b", "c"],
    "Temkin": ["q_sat", "b", "c"],
    "BingelWalton": ["q_sat", "a", "b"]
}



class Components:
    def __init__(self, components: list[dict] = []):
        self.components = []
        self.labels = []
        self.carrierGas = None

        for comp in components:
            self.addComponent(**comp)

    def addComponent(
        self,
        name: str,
        gasPhaseMolFraction: float,
        isotherms: list = [],
        massTransferCoefficient: float = 0.0,
        axialDispersionCoefficient: float = 0.0,
        isCarrierGas: bool = False,
    ):
        # get idx from existing components
        idx = len(self.components)

        for site, isotherm in enumerate(isotherms):
            self.labels += [f"c{idx}_s{site}_t{isothermTypes[isotherm[0]]}_{label}" for label in isothermLabels[isotherm[0]]]

        # add isotherm information
        cpp_isotherms = [
            _ruptura.Isotherm(isothermTypes[isotherm[0]], isotherm[1:], len(isotherm) - 1)
            for isotherm in isotherms
        ]
        comp = _ruptura.Component(
            idx,
            name,
            cpp_isotherms,
            gasPhaseMolFraction,
            massTransferCoefficient,
            axialDispersionCoefficient,
            isCarrierGas,
        )

        if isCarrierGas:
            self.carrierGas = idx

        self.components.append(comp)

    def getLabels(self) -> list[str]:
        return [f"{comp.name} (y_i={comp.gasPhaseMolFraction})" for comp in self.components]



class MixturePrediction:
    def __init__(
        self,
        components: Components,
        displayName: str = "Column",
        temperature: float = 433.0,
        pressureStart: float = -1.0,
        pressureEnd: float = -1.0,
        numberOfPressurePoints: int = 100,
        pressureScale: str = "log",
        predictionMethod: str = "IAST",
        iastMethod: str = "FastIAST",
    ):
        self.shape = (numberOfPressurePoints, len(components.components), 6)
        pressureScale = {"log": 0, "linear": 1}[pressureScale]
        predictionMethod = {"IAST": 0, "SIAST": 1, "EI": 2, "SEI": 3}[predictionMethod]
        iastMethod = {"FastIAST": 0, "NestedLoopBisection": 1}[iastMethod]
        self.components = components

        # on carriergas: ruptura first checks if the numberOfCarrierGases is 0 or more. more is
        # always 1, as carrierGasComponent is size_t. if the carriergas is notpresent, set
        # component to 0 (default), it is not checked.

        self.MixturePrediction = _ruptura.MixturePrediction(
            displayName,
            components.components,
            1 if components.carrierGas is not None else 0,
            components.carrierGas or 0,
            temperature,
            pressureStart,
            pressureEnd,
            numberOfPressurePoints,
            pressureScale,
            predictionMethod,
            iastMethod,
        )
        print(self.MixturePrediction)
        self.data = None

    def compute(self):
        self.data = self.MixturePrediction.compute()
        return self.data
    
    def setComponentsParameters(self, params: np.ndarray):
        self.data = None
        self.MixturePrediction.setComponentsParameters(params)

    def getComponentsParameters(self):
        return self.MixturePrediction.getComponentsParameters()
    
    def plot(self, ax, plot_type: Literal["pure", "mixture", "mixture_molfrac"]):
        select = {"pure": 1, "mixture": 2, "mixture_molfrac": 4}[plot_type]

        if self.data is None:
            raise ValueError("Data not computed yet")

        # set axes
        ax.set_xlabel("Total bulk fluid phase fugacity, f/Pa")
        ylabel = {"pure":"Absolute loading q_i", "mixture": "Absolute loading q_i", "mixture_molfrac": "Adsorbed mol-fraction Y_i"}
        ax.set_ylabel(ylabel[plot_type])
        ax.set_xscale("log")

        # set values for loop
        ncomp = self.data.shape[1]
        labels = self.components.getLabels()
        markers = ["o", "+", "^", "D", "x", "*", "p", "s", "v"]

        # plot all components
        for comp in range(ncomp):
            ax.scatter(self.data[:, comp, 0], self.data[:, comp, select], label=labels[comp], marker=markers[comp], s=8.0)
        ax.legend()


class Breakthrough:
    def __init__(
        self,
        components: Components,
        displayName: str = "Column",
        temperature: float = 433.0,
        numberOfTimeSteps: Union[int, str] = "auto",
        numberOfGridPoints: int = 100,
        printEvery: int = 10000,
        writeEvery: int = 10000,
        totalPressure: float = 1e6,
        columnVoidFraction: float = 0.4,
        pressureGradient: float = 0.0,
        particleDensity: float = 1e3,
        columnEntranceVelocity: float = 0.1,
        columnLength: float = 0.3,
        timeStep: float = 5e-4,
        pulseTime: float = None,
    ):
        # take 1e6 as max number of timesteps if autosteps is used
        autoSteps = numberOfTimeSteps == "auto"
        numberOfTimeSteps = 0 if numberOfTimeSteps == "auto" else int(numberOfTimeSteps)
        pulse = pulseTime is not None
        pulseTime = 0 if pulseTime is None else pulseTime
        carrierGas = components.carrierGas if components.carrierGas else 0
        self.components = components
        mix = MixturePrediction(displayName=displayName, temperature=temperature, components=components)

        self.shape = (
            numberOfTimeSteps // writeEvery,
            numberOfGridPoints,
            len(components.components) + 1,
            6,
        )

        self.Breakthrough = _ruptura.Breakthrough(
            displayName,
            components.components,
            carrierGas,
            numberOfGridPoints,
            printEvery,
            writeEvery,
            temperature,
            totalPressure,
            columnVoidFraction,
            pressureGradient,
            particleDensity,
            columnEntranceVelocity,
            columnLength,
            timeStep,
            numberOfTimeSteps,
            autoSteps,
            pulse,
            pulseTime,
            mix.MixturePrediction,
        )

        self.data = None
        print(self.Breakthrough)

    def compute(self):
        self.data = self.Breakthrough.compute()
        return self.data

    def plot(self, ax, plot_type: Literal["breakthrough", "Dpdt", "Dqdt", "P", "Pnorm", "Pt", "Q", "Qeq", "V"]):
        if self.data is None:
            raise ValueError("Data not computed yet")
        
        # set values for loop
        ncomp = (self.data.shape[-1] - 5) // 6
        labels = self.components.getLabels()
        markers = ["o", "+", "^", "D", "x", "*", "p", "s", "v"]

        if plot_type == "breakthrough":
            # plot all components
            
            x = self.data[:, -1, 0]
            ax.set_xlim(x.min(), x.max())
            ax.set_xticks(np.linspace(x.min(), x.max(), 7))
            ax.set_xlabel("Dimensionless time")
            ax.set_ylabel("Concentration exit gas")

            x2 = self.data[:, -1, 1]
            ax2 = ax.twiny()
            ax2.set_xticks(np.linspace(x2.min(), x2.max(), 7))
            ax2.set_xlabel("Time (min)")

            for comp in range(ncomp):
                ax.scatter(x, self.data[:, -1, 8+comp*6], label=labels[comp], marker=markers[comp%len(markers)], s=8.0 )


        ax.legend()