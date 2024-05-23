## Setup python environment for MapMyCells

### Step1: Python environment installation using Miniconda

*Note*: skip the substeps 1 and 2 if you already have a python virtual environment set up. 

#### 1.1. Install Miniconda

In your terminal window, run:
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

*Note*:
* These commands will install Miniconda to your home directory. To install it in a different directory, replace '~' with the path to the folder in all of the commands below. 
* For more help with installing miniconda visit [https://docs.anaconda.com/free/miniconda/miniconda-install.html](https://docs.anaconda.com/free/miniconda/)

#### 1.2. Create a virtual environment

In your terminal window, run:
```
conda create --name <your-env-name> -y
```

More information: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

*Note*:
The conda enviromnet needs to be activated before running cell_type_mapper.

#### 1.3. Install cell_type_mapper library

In your terminal window, run:
```
conda activate <your-env-name>
git clone https://github.com/AllenInstitute/cell_type_mapper.git <path_to_your_directory>
cd <path_to_your_directory>/cell_type_mapper
pip install -r requirements.txt
pip install -e .
```

*The coding steps above are*:
* Activate the conda environment.
* Clone cell_type_mapper repository to the local directory.
* Go to cell_type_mapper folder in the local directory.
* Install the required packages (dependencies).
* Install the cell_type_mapper package itself. Note: this needs to be run from the root directory of this repository, i.e. <path_to_your_directory>/cell_type_mapper (done on line 3).
