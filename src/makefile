CXXFLAGS=-g -O3 -std=c++17 -march=native -ffast-math -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Werror -fomit-frame-pointer -ftree-vectorize -fno-stack-check -funroll-loops $(shell python3 -m pybind11 --includes)

default: ruptura;

random_numbers.o: random_numbers.cpp random_numbers.h
	$(CXX) $(CXXFLAGS) -c random_numbers.cpp

special_functions.o: special_functions.cpp special_functions.h
	$(CXX) $(CXXFLAGS) -c special_functions.cpp

isotherm.o: isotherm.cpp isotherm.h
	$(CXX) $(CXXFLAGS) -c isotherm.cpp

multi_site_isotherm.o: multi_site_isotherm.cpp multi_site_isotherm.h isotherm.h
	$(CXX) $(CXXFLAGS) -c multi_site_isotherm.cpp

component.o: component.cpp component.h isotherm.h
	$(CXX) $(CXXFLAGS) -c component.cpp

mixture_prediction.o: mixture_prediction.cpp mixture_prediction.h
	$(CXX) $(CXXFLAGS) -c mixture_prediction.cpp

PYBIND11_INCLUDE = /Users/youriran/mambaforge/envs/compchem/include/pybind11

ruptura_module.o: ruptura_module.cpp
	$(CXX) $(CXXFLAGS) -I$(PYBIND11_INCLUDE) -c ruptura_module.cpp

ruptura_module.cpp: mixture_prediction.h
	$(CXX) $(CXXFLAGS) -I$(PYBIND11_INCLUDE) -E -P -o ruptura_module.cpp mixture_prediction.h

inputreader.o: inputreader.cpp inputreader.h
	$(CXX) $(CXXFLAGS) -c inputreader.cpp

breakthrough.o: breakthrough.cpp breakthrough.h
	$(CXX) $(CXXFLAGS) -c breakthrough.cpp

fitting.o: fitting.cpp fitting.h
	$(CXX) $(CXXFLAGS) -c fitting.cpp

main.o : main.cpp 
	$(CXX) $(CXXFLAGS) -c main.cpp

ruptura: random_numbers.o special_functions.o isotherm.o multi_site_isotherm.o component.o mixture_prediction.o inputreader.o breakthrough.o fitting.o main.o
	$(CXX) main.o fitting.o breakthrough.o inputreader.o mixture_prediction.o component.o multi_site_isotherm.o isotherm.o special_functions.o random_numbers.o ruptura_module.o -o ruptura

clean:
	rm -f *.pcm *.o *.a ruptura
