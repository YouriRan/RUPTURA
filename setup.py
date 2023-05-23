from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
    "ruptura",
    sources=glob("src/*.cpp"),
    extra_compile_args=["-std=c++17"]
    )
]

setup(
    name="ruptura",
    version="0.0.1",
    author=[
        "Shrinjay Sharma, Delft University of Technology, The Netherlands",
        "Salvador R.G. Balestra, Pablo de Olavide University, Spain",
        "Richard Baur, Shell Global Solutions International B.V., The Netherlands",
        "Umang Agarwal, Shell Global Solutions International B.V., The Netherlands",
        "Eric Zuidema, Shell Global Solutions International B.V., The Netherlands",
        "Marcello Rigutto, Shell Global Solutions International B.V., The Netherlands",
        "Sofia Calero, Eindhoven University of Technology, The Netherlands",
        "Thijs J.H. Vlugt, Delft University of Technology, The Netherlands",
        "David Dubbeldam, University of Amsterdam, The Netherlands"
    ],
    ext_modules=ext_modules,
    zip_safe=False,
)