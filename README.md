# Influence_of_hydrogen_trapping

A repository for the publiction "Influence of hydrogen trapping on WCLL breeding blanket performances", and can be found at: https://doi.org/10.1088/1741-4326/ac28b0   

<br>
All simulation scripts are provided along with associated plotting scripts for results, which are found in the publication


## Run FESTIM in Docker
The FEniCS project provides a [Docker image](https://hub.docker.com/r/fenicsproject/stable/) with FEniCS and its dependencies (python3, UFL, DOLFIN, numpy, sympy...)  already installed. See their ["FEniCS in Docker" manual](https://fenics.readthedocs.io/projects/containers/en/latest/).

Get Docker [here](https://www.docker.com/community-edition).

Pull the Docker image and run the container, sharing a folder between the host and container:

For Windows users:
```python
docker run -ti -v ${PWD}:/home/fenics/shared --name fenics quay.io/fenicsproject/stable:latest
```

Install a FESTIM v0.7
```python
pip install git+https://github.com/RemDelaporteMathurin/FESTIM@v0.7.1
```
To run the FESTIM simulations:
```python
python3 main.py
```
