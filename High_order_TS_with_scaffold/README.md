# Installation of the packages for the Scaffold computation


The bulk of the code in this folder is exactly the same as the one provided in the folder "High_order_TS". The main difference is that of the computation of the homological scaffold 

After installing the python packages using pip, i.e.:

```
pip install --use-deprecated=legacy-resolver -r requirements.txt
```

You will also need to install and download other softwares:

- Jython 2.7.3 (https://www.jython.org) : After downloading the jython installer (https://repo1.maven.org/maven2/org/python/jython-installer/2.7.3/jython-installer-2.7.3.jar ), it can be installed using  ```java -jar jython-installer-2.7.3.jar```. It is then important to include jython in the system path (e.g. ```export PATH=/home/$USER/jython/bin:$PATH```)


There are other packages/files that are already present in this folder, which are listed below:

- Javaplex (https://github.com/appliedtopology/javaplex). In the folder "javaplex", it is already available the file javaplex.jar which is used for computing the scaffold 

- Holes (an adaptation of the python package developed here: https://github.com/lordgrilo/Holes). 
