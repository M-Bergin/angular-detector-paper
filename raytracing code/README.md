# Scattering from Rough Surfaces
## Instructions
Clone [github.com/slambrick/SHeM-Ray-Tracing-Simulation.git](https://github.com/slambrick/SHeM-Ray-Tracing-Simulation):

```git clone https://github.com/slambrick/SHeM-Ray-Tracing-Simulation.git```

Then copy this directory into the root of the cloned source.

## Scattering from 1D Rough Surfaces Simulations
The script is `bentley_specScatter2D`. 

## Scattering from 2D Rough Surfaces Simulations
The script is `bentley_scatterMany.m`. It will save the results from simulations into a directory `bentley_scatterAndSave` as specified in `bentley_scatterAndSave3D.m`. One of the lines in `bentley_scatterMany.m` can be uncommented and another re-commented such that instead of rerunning the simulations, the plots will be produced without rerunning the simulation. Within `bentley_scatterMany.m`, the two scattering distributions that can be called are `specular` and `cosine`. The plotting code can be modified by editing `bentley_plotScatter3D.m` for the specular case or `bentley_plotScatter3D_DIFFUSE.m` for the diffuse case.

## Scattering from gratings
The script is `Gratings_scatter.m`