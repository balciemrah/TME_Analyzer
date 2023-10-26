# TME_Analyzer-3810_test
## Does this run in the Stubbs group?

instructons:
get python version 3.8.10
- you might have to activate your environment for following steps; 
(I needed a specific powershell code to activate `powershell -ExecutionPolicy Bypass -File "./.venv/Scripts/Activate.ps1"`)

run code
`python.exe -m pip install -r requirements.txt`
- note that requirements.txt give the minimal packages, and this should get all the packages in the requirements_extended.txt file

you can run the TME-Analyzer after this;
`python.exe TME_analyzer_test.py`
this should get the interface to run

you should also be able to pack the .py file with Cx_Freeze:
`python setup.py build`

