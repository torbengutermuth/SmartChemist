# SmartChemist
This repository is the backend of the SmartChemist, which shall enable anyone to sound like a smart Chemist.

Once running, you can supply a string of comma separated SMILES, a SMILES file (.smi) or a SDFile (.sdf) and the names of common functional groups, ringsystems and biologicals will be displayed.
There are three different types of patterns displayed, ringsystems and biologicals (cyclic), functional groups (functional_groups) and patterns that are overshadowed by others (overshadowed).
If the match of a substructure is u subset of another (e.g. hydroxyl group in a carbonic acid), this pattern is overshadowed and shall not be displayed as a default setting.

## Installation

As the SmartChemist is not yet published, this repository remains WIP.  
All SMARTS necessary for setting up the SmartChemist will only be published once the publication is accepted.
In addition, the FrontEnd is still under development and not yet publicly available. 

To setup a testserver and test it using CURL, please execute the following steps:

First, create and activate a [CONDA](https://docs.conda.io/projects/conda/en/latest/index.html) environment using the environment.yml file.

```bash
conda create env create -f environment.yml
conda activate your_conda_environment
```

Now, it is time to start setting up the webserver. This can be done in three steps. 
1. Move into the directory and generate the database necessary for the webserver:
```bash
cd smart_chemist_backend
python manage.py migrate
```

3. Add the SMARTS patterns to the database
```bash
python add_patterns_to_db.py
```

4. Start all necessary Software for task execution using the start.sh script.
```bash
./start.sh
```

5. Start the server
```bash
python manage.py runserver
```

Now, we can test if the server executes properly which takes another three steps

1. Upload the string/file to execute to the webserver
For a SMILES string:
```bash
curl -vvv -d 'smiles=S(=O)(=O)(NC(=O)Nc1nc(OC)cc(OC)n1)Cc2c(cccc2)C(=O)OC' localhost:8000/api/names
```
For a file:
```bash
curl -vvv -X POST -F molecule_file=@'/path/to/your/file/ligands.sdf' localhost:8000/api/names
```
The result should contain a job id like "a7d1ac11-c111-4d4f-a644-b86bd8a5fef8"

2. Query the job-id for execution status
```bash
curl -vvv localhost:8000/jobs/a7d1ac11-c111-4d4f-a644-b86bd8a5fef8/
```
There, you can find an output_info id if it successfully finished

3. Query the output of the successful job
```bash
curl -vvv localhost:8000/jobs_output/213/
```
The result is a bit tedious to look at, especially for higher numbers of molecules, but it should contain an SVG and a matches dictionary at the end.
If you see some SMARTS at the end it is probably fine, please be aware that only few SMARTS are in the database and use a working molecule.

There are some tests which will be extended in the coming weeks and automatic CI/CD will be setup using Github.