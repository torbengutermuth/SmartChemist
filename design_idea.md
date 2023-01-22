# Smart Chemist Design Idea

The general idea is to sound smart because we aren't.

We accept Molecules as SMILES or SD Files. 
The Server matches all known SMARTS present in the data base to all molecules in the input molecule file. 
The Server sends an SVG file including all matching SMART expressions including notation with their respective atom/bond indexes. 


Frontend allows uploading of molecule file. 
Frontend sends molecule file to backend. 
Backend does magic. 
Backend sends JSON back to frontend. 
JSON is centered around molecules present in file
molecule -> smarts/expression/name -> Atom/Bond Indeces

domain smart_chemist.com

/api/names


json = [{Name:str, svg:str, matches:[{atom_indexes:list, trivial_name:{name:str,smarts:str}}]}]