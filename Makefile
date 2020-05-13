#OUTPUT := crossmods$(shell python3-config --extension-suffix)
#
#$(OUTPUT): crossmods.cpp
#	c++ -Ofast -Wall -shared -std=c++17 -fPIC `python3 -m pybind11 --includes` -fopenmp crossmods.cpp -o $(OUTPUT)

#crossmods.cpp: crossmods.hpp vddm.hpp bindgen.py
#	python3 bindgen.py crossmods crossmods.hpp > crossmods.cpp

benchmark: benchmark.cpp vddm.hpp tdm.hpp
	c++ -g -Ofast -Wall -std=c++17 -fopenmp $< -o $@

