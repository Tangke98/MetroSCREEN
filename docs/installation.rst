.. highlight:: shell

.. role:: bash(code)
   :language: bash

Installation
------------




System requirements
>>>>>>>>>>>>>>>>>>>

* Linux/Unix
* R >= 4.1.3


We recommend to create an independent conda environment for MetroSCREEN. If users do not have conda, please install Miniconda first:
::
   
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh


Installation
>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for MetroSCREEN.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n MetroSCREEN
   conda activate MetroSCREEN

Step 2 Install MetroSCREEN package from :bash:`conda`.
::::::::::::::::::::::::::::::::::::::::::::::::
::

   conda install tangke::metroscreen