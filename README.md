# pyemma-markov-state-models
##start your MSM analysis
```
conda create -n pyemma-env python=3.7.12
```

```
conda activate pyemma-env
```

```
conda install pip
```


python -m pip install pyemma==2.5.7

python -m pip install pandas==0.25.3

python -m pip install notebook

python msm_hlu.py



####################################################
Readmore:
#Make sure your conda environment is conssitent!!!!!
##how to install pyemma2.5.7 in Ubuntu18.04

##attention: Pyemma2.5.7 doesn't exit in conda-forge channel!!!!!!

##Do not install pyemma through commands: 1. conda install --channel "conda-forge" pyemma,
###and then 2. conda install notebook, cause numpy version installed will be the rabbit hole you DONOT want to explore! It is a waste of time!!!

#After creating the env environment with the installation of python3.7.12
python3 -m pip install pyemma==2.5.7

#Use Python 3.7.12 and install pyemma2.5.7, after succesffuly installation, the following packages were also installed:
astunparse-1.6.3 bhmm-0.6.3 cycler-0.11.0 decorator-5.1.1 dill-0.3.7 fonttools-4.38.0 h5py-3.8.0
kiwisolver-1.4.5 matplotlib-3.5.3 mdtraj-1.9.9 msmtools-1.2.6 multiprocess-0.70.15 numpy-1.21.6
packaging-24.0 pathos-0.3.1 pillow-9.5.0 pox-0.3.3 ppft-1.7.6.7 psutil-7.0.0 pyemma-2.5.7 pyparsing-3.1.4
python-dateutil-2.9.0.post0 pyyaml-6.0.1 scipy-1.7.3 six-1.17.0 tqdm-4.67.1 typing-extensions-4.7.1

#Then (yes or no, not mandatory):
python3 -m pip install pandas==pandas0.25.3

#Then install jupyter notebook
python3 -m pip install notebook
Successfully installed Send2Trash-1.8.3 anyio-3.7.1 argon2-cffi-23.1.0 argon2-cffi-bindings-21.2.0
attrs-24.2.0 backcall-0.2.0 beautifulsoup4-4.13.3 bleach-6.0.0 cffi-1.15.1 debugpy-1.7.0 defusedxml-0.7.1
entrypoints-0.4 exceptiongroup-1.2.2 fastjsonschema-2.21.1 idna-3.10 importlib-metadata-6.7.0
importlib-resources-5.12.0 ipykernel-6.16.2 ipython-7.34.0 ipython-genutils-0.2.0 jedi-0.19.2 jinja2-3.1.5
jsonschema-4.17.3 jupyter-client-7.4.9 jupyter-core-4.12.0 jupyter-server-1.24.0 jupyterlab-pygments-0.2.2
markupsafe-2.1.5 matplotlib-inline-0.1.6 mistune-3.0.2 nbclassic-1.2.0 nbclient-0.7.4 nbconvert-7.6.0
nbformat-5.8.0 nest-asyncio-1.6.0 notebook-6.5.7 notebook-shim-0.2.4 pandocfilters-1.5.1 parso-0.8.4
pexpect-4.9.0 pickleshare-0.7.5 pkgutil-resolve-name-1.3.10 prometheus-client-0.17.1 prompt-toolkit-3.0.48
ptyprocess-0.7.0 pycparser-2.21 pygments-2.17.2 pyrsistent-0.19.3 pyzmq-26.2.1 sniffio-1.3.1 soupsieve-2.4.
terminado-0.17.1 tinycss2-1.2.1 tornado-6.2 traitlets-5.9.0 wcwidth-0.2.13 webencodings-0.5.1
websocket-client-1.6.1 zipp-3.15.0
###############################################################

