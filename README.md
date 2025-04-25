# Modeling and Simulation of Selective Laser Sintering with COMSOL Multiphysics

## Introduction

Among the various additive manufacturing (AM) techniques, powder bed fusion methods such as Selective Laser Melting (SLM) and Selective Laser Sintering (SLS) have garnered significant attention due to their ability to fabricate complex three-dimensional structures with high precision. SLM involves the use of a high-powered laser to fully melt metallic powder, creating dense, fully functional parts, whereas SLS typically employs a lower-energy laser to sinter polymer or metal powders without fully melting them. For the purposes of this work, the term SLS will be used to refer to both of these AM techniques.

SLS offers several advantages that make it an attractive choice for advanced manufacturing applications. It provides exceptional design flexibility, allowing for the production of intricate geometries that would be difficult or impossible to achieve using traditional subtractive manufacturing techniques.

## Getting Started

In order to run the simulations, it is important to have installed the software COMSOL Multiphysics 6.1. First, clone the repository into your local filesystem.

```git clone https://github.com/finalbossqc/SLS.git```

```cd SLS/Codes```

### Generating MPH files

There are two main models included in this repository ```Molybdenum2DModel.java``` and ```Polyamide2DModel.java```. The former is for a 2D SLS simulation with Molybdenum and the latter is for a 2D simulation with Nylon 12 (Polyamide 12). 

The Java files have to be compiled using the COMSOL 6.1 compiler. 

If you are on a Linux setup, run the command:

```comsol compile Polyamide2DModel.java```

This will generate ```Polyamide2DModel.class```. Then run the command:

```comsol batch -inputfile Polyamide2DModel.class```

This will generate a .mph file, which can be opened using the COMSOL GUI. 

### Running batch jobs on SLURM with HPC setup

Once you have created the .mph file, you can run a batch job (a job that runs asynchronously) on an HPC system by running the shellscript ```comsolrun.sh```. 

The following statement runs Study 1 of the 2D Polyamide SLS model.

```sbatch comsolrun.sh Polyamide2DModel_Model.mph std1```

You can change the SLURM parameters by editing the shell script ```comsolrun.sh```. By default, the time limit is 8 hours for the job and it requests 10 cpu cores.

For 3D simulations, make sure to use the script ```comsolrunbig.sh```, which requests 20 cpu cores for 12 hours. 

### Opening GUI on SLURM with HPC setup

To open the GUI, first ensure that you have an X11 client enabled. This is usually done by specifying "-X" when ssh-ing into the login node of the HPC cluster. On the UNC Longleaf cluster, this would be

```ssh -X <onyen>@longleaf.unc.edu```

Then import the COMSOL module:

```module load comsol/6.1```

Open the GUI with:

```comsol```

Then navigate to File > Open and open the .mph file generated in the previous steps.

On an HPC cluster, it is usually bad manners to do this on the login node, so with SLURM should do something like this:

```salloc --mem 32G --time 2:00:00 -c 5```

```ssh -X $SLURM_NODELIST```

```module load comsol/6.1```

```comsol```

Then open your model.
