"""Original publication credits and the requested RUPTURA 2.0 author roster."""

from pathlib import Path

ABOUT = {
    "name": "Ruptura Lab",
    "version": "2.0",
    "description": "A local workspace for adsorption simulations, breakthrough curves, mixture equilibria, and interactive analysis.",
    "repository": "https://github.com/iRASPA/RUPTURA",
    "license": "MIT",
    "contact": {"name": "David Dubbeldam", "email": "d.dubbeldam@uva.nl"},
    "authors": [
        {"name": "Youri Ran", "affiliation": "University of Amsterdam"},
        {"name": "Shrinjay Sharma", "affiliation": "Delft University of Technology"},
        {"name": "Salvador R. G. Balestra", "affiliation": "Universidad Pablo de Olavide"},
        {"name": "Richard Baur", "affiliation": "Shell Global Solutions International B.V."},
        {"name": "Umang Agarwal", "affiliation": "Shell Global Solutions International B.V."},
        {"name": "Muhammed Hassan", "affiliation": "Pusan National University"},
        {"name": "Erik Zuidema", "affiliation": "Shell Global Solutions International B.V."},
        {"name": "Marcello S. Rigutto", "affiliation": "Shell Global Solutions International B.V."},
        {"name": "Sofia Calero", "affiliation": "Eindhoven University of Technology"},
        {"name": "Thijs J. H. Vlugt", "affiliation": "Delft University of Technology"},
        {"name": "David Dubbeldam", "affiliation": "University of Amsterdam"},
    ],
    "citation": {
        "title": "RUPTURA: Simulation Code for Breakthrough, Ideal Adsorption Solution Theory Computations, and Fitting of Isotherm Models",
        "authors": "Shrinjay Sharma, Salvador R. G. Balestra, Richard Baur, Umang Agarwal, Erik Zuidema, Marcello S. Rigutto, Sofia Calero, Thijs J. H. Vlugt, and David Dubbeldam",
        "publication": "Molecular Simulation 49(9), 893–953 (2023)",
        "doi": "10.1080/08927022.2023.2202757",
        "bibtex": """@article{sharma2023ruptura,
  title = {{RUPTURA}: Simulation Code for Breakthrough, Ideal Adsorption Solution Theory Computations, and Fitting of Isotherm Models},
  author = {Sharma, Shrinjay and Balestra, Salvador R. G. and Baur, Richard and Agarwal, Umang and Zuidema, Erik and Rigutto, Marcello S. and Calero, Sofia and Vlugt, Thijs J. H. and Dubbeldam, David},
  journal = {Molecular Simulation},
  year = {2023},
  volume = {49},
  number = {9},
  pages = {893--953},
  doi = {10.1080/08927022.2023.2202757}
}
""",
    },
    "references": [
        {
            "topic": "Ideal adsorbed solution theory",
            "title": "Thermodynamics of Mixed-Gas Adsorption",
            "authors": "A. L. Myers and J. M. Prausnitz",
            "publication": "AIChE Journal 11, 121–127 (1965)",
            "doi": "10.1002/aic.690110125",
        },
        {
            "topic": "T15 breakthrough validation",
            "title": "Prediction of Non-Isothermal Ternary Gas-Phase Breakthrough Experiments Based on Binary Data",
            "authors": "Dorian Marx, Lisa Joss, Nathalie Casas, Johanna Schell, and Marco Mazzotti",
            "publication": "Adsorption 20, 493–510 (2014)",
            "doi": "10.1007/s10450-013-9593-5",
        },
        {
            "topic": "Monolith validation",
            "title": "On Comparing Packed Beds and Monoliths for CO₂ Capture from Air Through Experiments, Theory, and Modeling",
            "authors": "Valentina Stampi-Bombelli, Alba Storione, Quirin Grossmann, and Marco Mazzotti",
            "publication": "Industrial & Engineering Chemistry Research 63, 11637–11653 (2024)",
            "doi": "10.1021/acs.iecr.4c01392",
        },
    ],
}


ABOUT["licenseText"] = Path(__file__).with_name("LICENSE.txt").read_text()


def example_notes(experiment):
    """Attach references to their experiments without adding metadata to CLI JSON."""
    name = experiment["input"].get("DisplayName", "").lower()
    key = experiment.get("presetId", "")
    references = ABOUT["references"]
    if key == "paper-t15" or "t15" in name:
        selected = [references[0], references[1]]
    elif "monolith" in key or "400 ppm" in name:
        selected = [references[2]]
    elif "bea" in key or "bea" in name:
        selected = [references[0]]
    elif experiment["input"].get("MixturePredictionMethod", "IAST") in ("IAST", "SIAST"):
        selected = [references[0]]
    else:
        return ""
    return "References\n\n" + "\n\n".join(
        f"{item['authors']}. {item['title']}. {item['publication']}. https://doi.org/{item['doi']}"
        for item in selected
    )
