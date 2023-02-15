# SmartChemist
This repository is the backend of the SmartChemist, which shall enable anyone to sound like a smart Chemist.

Once running, you can supply a string of comma separated SMILES, a SMILES file (.smi) or a SDFile (.sdf) and the names of common functional groups and ringsystems will be displayed.
There are three different types of patterns displayed, ringsystems (cyclic), functional groups (functional_groups) and patterns that are overshadowed by others (overshadowed).
If the match of a substructure is u subset of another (e.g. hydroxyl group in a carbonic acid), this pattern is overshadowed and shall not be displayed as a default setting.

Currently, the link between Frontend and Backend only works for SMILES strings.

## Installation

First, create a virtual environment

    python -m venv smartchemist

Then, install the necessary dependencies in requirements.txt

    pip install -r requirements.txt

Activate your virtual environment, in linux this can be done with the following command:

    source smartchemist/bin/activate

After taking care of the python dependencies, we have to include the Frontend to easily view the patterns in your molecules.
The Frontend is included as a submodule inside git and can be activated using the following commands:
    
    git submodule init
    git submodule update

Now, it is time to start the webserver. This can be done in three steps. 
1. Generate the database necessary for the webserver:

        python manage.py migrate
2. Add the SMARTS patterns to the database

        python add_patterns_to_db.py
 
3. Start the server

        python manage.py runserver

You can find the Website [here](http://localhost:8000/static/smart_chemist_frontend/index.html) 

**Have Fun!**

If you find any Problems with the Server or the Patterns, feel free to contact me!