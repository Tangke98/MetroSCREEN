.. highlight:: shell

.. role:: bash(code)
   :language: bash

Installation
------------




System requirements
>>>>>>>>>>>>>>>>>>>

* Linux/Unix
* Python >= 3.10


We recommend to create an independent conda environment for MetroSCREEN. If users do not have conda, please install Miniconda first:
::
   
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh


Install the stable version
>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for MetroSCREEN.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n MetroSCREEN python=3.10
   conda activate MetroSCREEN

Step 2 Install MetroSCREEN package from :bash:`pypi`.
::::::::::::::::::::::::::::::::::::::::::::::::
::

   pip install -U MetroSCREEN


Install the developing version
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for MetroSCREEN.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n MetroSCREEN python=3.10
   conda activate MetroSCREEN

Step 2 Download MetroSCREEN package from github.
:::::::::::::::::::::::::::::::::::::::::::
::

   git clone https://github.com/wanglabtongji/MetroSCREEN.git

Step 3 Install dependencies of MetroSCREEN.
::::::::::::::::::::::::::::::::::::::
::

   cd MetroSCREEN
   pip install -r requirements.txt

Step 4 Install MetroSCREEN.
::::::::::::::::::::::
::
  
   pip install .
