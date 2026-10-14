@echo off
rem Build the Link Budget Tool executable and the zip to deliver.
rem
rem Double-click this file (or run it from a command prompt) in the project
rem folder. It:
rem   1. creates the .venv the first time and installs the requirements and
rem      PyInstaller into it (later runs only add what is missing);
rem   2. deletes the old build and dist folders and runs PyInstaller;
rem   3. copies the data files next to the executable;
rem   4. zips the folder as dist\Link-Budget-Tool_<version>.zip (version from
rem      version.py).
rem The licence file (license.key) is per customer: add it to the folder or
rem the zip yourself (see README).

setlocal
cd /d "%~dp0"
set APP=Link-Budget-Tool
set PY=.venv\Scripts\python.exe

if not exist %PY% (
    echo === Creating the virtual environment .venv
    python -m venv .venv || goto :error
)

echo === Installing the requirements into .venv
%PY% -m pip install --quiet --upgrade pip || goto :error
%PY% -m pip install --quiet -r requirements.txt pyinstaller || goto :error

for /f "delims=" %%v in ('%PY% -c "import version; print(version.APP_VERSION)"') do set VERSION=%%v
if "%VERSION%"=="" goto :error
echo === Building %APP% %VERSION%

if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist dist (
    echo.
    echo The old dist folder could not be deleted: some of its files are in use.
    echo Close %APP%.exe ^(also in Task Manager^) and any window open on dist, then run build.bat again.
    goto :error
)
%PY% -m PyInstaller --noconfirm gui.spec || goto :error

echo === Copying the data files
for %%f in ("ground_stations.txt" "tle.txt" "parameters.json" "Antenna Pattern Default.txt") do (
    if exist %%f (
        copy /y %%f "dist\%APP%\" >nul || goto :error
    ) else (
        echo WARNING: %%f not found, not copied
    )
)

echo === Creating dist\%APP%_%VERSION%.zip
powershell -NoProfile -Command "Compress-Archive -Path 'dist\%APP%' -DestinationPath 'dist\%APP%_%VERSION%.zip' -Force" || goto :error

echo.
echo Done: dist\%APP%\%APP%.exe and dist\%APP%_%VERSION%.zip
echo Remember to add the customer's license.key next to the executable.
pause
exit /b 0

:error
echo.
echo BUILD FAILED - see the messages above.
pause
exit /b 1
