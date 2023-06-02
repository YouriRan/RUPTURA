import numpy as np
import matplotlib.pyplot as plt
import ruptura

components = ruptura.Components()
components.addComponent(
    name="Helium",
    gasPhaseMolFraction=0.98,
    isCarrierGas=True
)
components.addComponent(
    name="nC6",
    gasPhaseMolFraction=0.01,
    isotherms=[["Langmuir", 1.1, 2.14e-4]],
    massTransferCoefficient=0.06,
    axialDispersionCoefficient=0.0
)
components.addComponent(
    name="2MP",
    gasPhaseMolFraction=0.01,
    isotherms=[["Langmuir", 1.4, 1.68e-4]],
    massTransferCoefficient=0.06,
    axialDispersionCoefficient=0.0
)

mix = ruptura.MixturePrediction(
    displayName="CoBDP",
    temperature=443.0, 
    pressureStart=-1.0,
    pressureEnd=-1.0,
    numberOfPressurePoints=100,
    pressureScale="log",
    components=components
)

brk = ruptura.Breakthrough(
    displayName="CoBDP",
    components=components,
    mixturePrediction=mix,
    temperature=443.0,
    numberOfTimeSteps="auto",
    numberOfGridPoints=100,
    printEvery=10000,
    writeEvery=10000,
    totalPressure=2.0e6,
    columnVoidFraction=0.4,
    pressureGradient=0.0,
    columnEntranceVelocity=0.5,
    columnLength=0.3, 
    timeStep=0.0005,
)
brk.compute()
print(brk.data.shape)