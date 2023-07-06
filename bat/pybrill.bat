cd %BRILL%\stage

if defined PYTHONPATH (python -i ..\python\pyBrill.py) else (
..\python_home\python -i ..\python\pyBrill.py)
