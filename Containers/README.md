## Building containers

In this folder, you can find the .def file for building a container either using Apptainer (http://apptainer.org/) or using Docker (https://www.docker.com/):

## 1. For Apptainer:

You should be able to build the container using the command:

```
apptainer build RHOST_container.sif RHOSTS.def
```

and then run the container using the command:
```
apptainer shell --no-home -e -B ./../:/repo RHOST_container.sif
```

After this you should be able to run all the bash files within the folder ```Example``` (or in apptainer, in ```/repo/Example/```)


## 2. For Docker:

You should be able to build the container using the command:

```
docker build -t rhosts-python . 
```

and then run the container using the command:
```
docker run --rm --mount type=bind,source="$(pwd)"/../,target=/repo -it --entrypoint bash rhosts-python
```

After this you should be able to run all the bash files within the folder ```/repo/Example/```
