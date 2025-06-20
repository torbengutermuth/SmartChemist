# SmartChemist
This repository is the backend of the SmartChemist, which shall enable anyone to sound like a smart Chemist.

Once running, you can supply a string of comma separated SMILES, a SMILES file (.smi) or a SDFile (.sdf) and the names of common functional groups, ringsystems and biologicals will be displayed.
There are three different types of patterns displayed, ringsystems and biologicals (cyclic), functional groups (functional_groups) and patterns that are overshadowed by others (overshadowed).
If the match of a substructure is u subset of another (e.g. hydroxyl group in a carbonic acid), this pattern is overshadowed and shall not be displayed as a default setting.

In scientific literature or written documentation, cite the following publication:
Gutermuth,T. et al., SmartChemist - Simplifying communication about organic chemical structures, Journal of Chemical Information and Modelling, 2025 (in press)

## Installation

To setup a server of your own, please execute the following steps:

First, create and activate a [CONDA](https://docs.conda.io/projects/conda/en/latest/index.html) environment using the environment.yml file.

```bash
conda env create -f environment.yml
conda activate SmartChemist
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

Now, we can test if the server executes properly which takes another three steps.
Alternatively you can skip to the frontend test if you want to include it.
But if there are issues on the frontend, testing this you can assure the backend works as expected.

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
If you see some SMARTS at the end it is probably fine.

There are some tests which will be extended in the coming weeks and automatic CI/CD will be setup using Github.

Now we can include and test the frontend.

```bash
cd smart_chemist_backend
mkdir static
cd static
git clone https://github.com/torbengutermuth/smart_chemist_frontend
```
If your server is still running, you should now be able to access the frontend on the website:

http://127.0.0.1:8000/static/smart_chemist_frontend/index.html

And now everything should work like it does on the live server, maybe a bit slower :D 

We use a postgres database in the live server, which is much faster than the sqlite version for the test server. 
To use an existing postgres database, change the Databases dictionary in the settings.py file in the following way: 

 DATABASES = {
     'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'postgres_db_name',
        'USER': 'postgres_db_user',
        'PASSWORD': 'postgres_db_password',
       'PORT': 'postgres_db_port',
     }
 }
