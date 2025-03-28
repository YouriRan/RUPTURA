import ruptura

components = ruptura.Components()
components.addComponent(
    name="nC6",
    gasPhaseMolFraction=0.2,
    isotherms=[["Langmuir-Freundlich", 1.1, 2.14e-4, 0.6], ["Langmuir-Freundlich", 4.8, 2.39e-5, 1.35]]
)
components.addComponent(
    name="2MP",
    gasPhaseMolFraction=0.2,
    isotherms=[["Langmuir-Freundlich", 1.4, 1.68e-4, 0.6], ["Langmuir-Freundlich", 4.8, 3.3e-5, 1.26]]
)
components.addComponent(
    name="3MP",
    gasPhaseMolFraction=0.2,
    isotherms=[["Langmuir-Freundlich", 1.75, 1.14e-4, 0.62], ["Langmuir-Freundlich", 4.8, 3.22e-5, 1.24]]
)
components.addComponent(
    name="23DMB",
    gasPhaseMolFraction=0.2,
    isotherms=[["Langmuir-Freundlich", 1.5, 8.31e-5, 0.63], ["Langmuir-Freundlich", 4.8, 5.46e-5, 1.13]]
)
components.addComponent(
    name="22DMB",
    gasPhaseMolFraction=0.2,
    isotherms=[["Langmuir-Freundlich", 1.74, 5.65e-5, 0.72], ["Langmuir-Freundlich", 4.8, 4.3e-5, 1.04]]
)

mix = ruptura.MixturePrediction(
    displayName="CoBDP",
    temperature=443.0, 
    pressureStart=1e2,
    pressureEnd=1e6,
    numberOfPressurePoints=100,
    pressureScale="log",
    components=components
)

pred = mix.compute()
print(pred.shape)