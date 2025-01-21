# Garfield

You will need environment variables defined as shown here.

export DCsim_HOME=/home/jmusser/gitrepo/GarfSim/DCSim

export GARFIELD_HOME=/home/jmusser/gitrepo/Garfield

export GARFIELD_IONDATA=/home/jmusser/gitrepo/Garfield/Data/IonMobility_CO2+_CO2.txt


The executables are designed to be run from scripts which define other local environment variables,  eg, genNormal.csh.

You will need the tcsh shell to run these as is. 

## Working with Docker

You can make an image with root and garfield pre-installed. To do so run:
```bash
docker build -t garfield . 
```
To create the image. Then to run the code in a container:
```bash
docker run --rm -it -v $(pwd):/mount garfield
```

This starts a container with the current directory (`$(pwd)`) mounted to `/mount` within the contain. Once in the container run:
```bash
source /mount/setup_env.sh
```

From there you need to recompile `DCSim`:
```bash
cd $DCsim_HOME ; make clean ; make
```

After that you can compile code in the `Garfield` directory, for example:
```bash
cd /mount/Garfield; 
make MakeGasFile_flight
```

Which can then be ran as:
```bash
./make MakeGasFile_flight
```
