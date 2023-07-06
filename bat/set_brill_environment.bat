rem @echo off

if not defined BRILL (
  if not exist %cd%\python_home (
    echo '*** not in BRILL home, change dir to your BRILL installation ***'
  ) else (
    set BRILL=%CD%
    set PYTHONHOME=%CD%\python_home;
    set PYTHONPATH=%CD%\python;%CD%\python_home
    bat\extend_path.bat %CD%\python_home\Scripts
    bat\extend_path.bat %CD%\python_home
    bat\extend_path.bat %CD%\bat
    bat\extend_path.bat %CD%\bin
    cd %CD%\stage
  )
)

cd stage
