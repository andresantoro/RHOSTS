# Installation of the packages

To install the python packages using pip, just run the following command:

```
pip install -r requirements.txt
```

------------------------------

If you have trouble with the previous command, try to upgrade the pip version and then try the following command:
```
pip install --use-deprecated=legacy-resolver -r requirements.txt
```
The usage of "--use-deprecated=legacy-resolver" is actually required only for the package "cechmate", since there are some conflicts in a python dependency (i.e. phat might give you this error: "it has inconsistent version: filename has '1.5.0a0', but metadata has '1.5.0'") 

This is a known issue, but if the issue persists after the previous command, which might be the case in MacOS, then try to manually install the phat software by downloading it from here: https://files.pythonhosted.org/packages/43/82/c14de81dc2953a71a060f72f2bc34c41996307956b162751f2a47e2c78f7/phat-1.5.0a.tar.gz#sha256=51e7fe5e05adf5c7e0895765572fff05b979731234251f13011610d71d4980ab 
and then install it using the following command in MacOS:

```
CFLAGS=-stdlib=libc++ python setup.py install
```

Also, if you get this issue with cechmate: "version `GLIBCXX_3.4.29` not found" and you are using a conda environment, then try to remove this file "/home/$USER/anaconda/lib/libstdc++.so.6". If you are using a specific env "env_name", the file should be located here "/home/$USER/envs/$ENV_NAME/lib/libstdc++.so.6"