import _ruptura
import numpy as np
from typing import Union

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


class Components:
    def __init__(self, components: list[dict] = []):
        self.components = []
        self.numberOfCarrierGases = 1
        self.carrierGasComponent = 0

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

        self.components.append(comp)


class MixturePrediction:
    def __init__(
        self,
        displayName: str,
        temperature: float,
        components: Components,
        pressureStart: float,
        pressureEnd: float,
        numberOfPressurePoints: int,
        pressureScale: str,
        predictionMethod: str = "IAST",
        iastMethod: str = "FastIAST",
    ):
        self.shape = (numberOfPressurePoints, len(components.components), 6)
        pressureScale = {"log": 0, "linear": 1}[pressureScale]
        predictionMethod = {"IAST": 0, "SIAST": 1, "EI": 2, "SEI": 3}[predictionMethod]
        iastMethod = {"FastIAST": 0, "NestedLoopBisection": 1}[iastMethod]

        self.MixturePrediction = _ruptura.MixturePrediction(
            displayName,
            components.components,
            components.numberOfCarrierGases,
            components.carrierGasComponent,
            temperature,
            pressureStart,
            pressureEnd,
            numberOfPressurePoints,
            pressureScale,
            predictionMethod,
            iastMethod,
        )
        print(self.MixturePrediction)

    def compute(self):
        self.data = self.MixturePrediction.compute()
        return self.data


class Breakthrough:
    def __init__(
        self,
        displayName: str,
        components: Components,
        mixturePrediction: MixturePrediction,
        temperature: float,
        numberOfTimeSteps: Union[int, str],
        numberOfGridPoints: int = 100,
        printEvery: int = 10000,
        writeEvery: int = 10000,
        totalPressure: float = 1e3,
        columnVoidFraction: float = 0.4,
        pressureGradient: float = 0.0,
        particleDensity: float = 1e3,
        columnEntranceVelocity: float = 0.1,
        columnLength: float = 0.3,
        timeStep: float = 5e-3,
        pulseTime: float = 0,
    ):
        # take 1e6 as max number of timesteps if autosteps is used
        autoSteps = numberOfTimeSteps == "auto"
        numberOfTimeSteps = int(1e6) if numberOfTimeSteps == "auto" else numberOfTimeSteps
        pulse = pulseTime is not None

        self.shape = (
            numberOfTimeSteps // writeEvery,
            numberOfGridPoints,
            len(components.components) + 1,
            6,
        )

        self.Breakthrough = _ruptura.Breakthrough(
            displayName,
            components.components,
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
            mixturePrediction.MixturePrediction,
        )
        print(self.Breakthrough)

    def compute(self):
        self.data = self.Breakthrough.compute()
        return self.data
