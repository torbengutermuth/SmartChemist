# SmartChemist
This repository is the backend of the SmartChemist, which shall enable anyone to sound like a smart Chemist.

Once running, you can supply a string of comma separated SMILES, a SMILES file (.smi) or a SDFile (.sdf) and the names of common functional groups, ringsystems and biologicals will be displayed.
There are three different types of patterns displayed, ringsystems and biologicals (cyclic), functional groups (functional_groups) and patterns that are overshadowed by others (overshadowed).
If the match of a substructure is u subset of another (e.g. hydroxyl group in a carbonic acid), this pattern is overshadowed and shall not be displayed as a default setting.

## Installation

First, create a virtual environment
```bash
python -m venv smartchemist
```
Activate your virtual environment, in linux this can be done with the following command:
```bash
source smartchemist/bin/activate
```
Then, install the necessary dependencies in requirements.txt
```bash
pip install -r requirements.txt
```

Now, it is time to start the webserver. This can be done in three steps. 
1. Move into the directory and generate the database necessary for the webserver:
```bash
cd smart_chemist_backend
python manage.py migrate
```

3. Add the SMARTS patterns to the database
```bash
python add_patterns_to_db.py
```

4. Open a separate terminal and start the redis server which is needed by celery which we use for distributed task execution.
```bash
redis-server --port 6378
```

5. Open another separate terminal and start the celery worker (note that this is development mode).
```bash
celery -A smart_chemist_backend worker -O fair --loglevel=INFO
```

6. Start the server
```bash
python manage.py runserver
```

A Frontend is currently being developed but not public yet.

**Have Fun!**

If you find any Problems with the Server or the Patterns, feel free to contact me!