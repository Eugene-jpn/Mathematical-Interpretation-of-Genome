@echo off
REM DNA Website Local Installation Setup Script
REM For Windows systems

echo 🧬 DNA Tandem Repeat Website - Local Installation
echo ==================================================

REM Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ Python is not installed. Please install Python 3.7 or higher.
    echo    Download from: https://python.org/downloads/
    echo    Make sure to check "Add Python to PATH" during installation.
    pause
    exit /b 1
)

REM Display Python version
echo ✅ Python detected
python --version

REM Create virtual environment
echo 📦 Creating virtual environment...
python -m venv venv

REM Activate virtual environment
echo 🔧 Activating virtual environment...
call venv\Scripts\activate.bat

REM Upgrade pip
echo ⬆️  Upgrading pip...
python -m pip install --upgrade pip

REM Install dependencies
echo 📥 Installing dependencies...
pip install flask flask-cors

if %errorlevel% neq 0 (
    echo ❌ Failed to install dependencies. Please check your internet connection.
    pause
    exit /b 1
)

echo ✅ Dependencies installed successfully!

REM Create a start script
echo 📝 Creating start script...
echo @echo off > start_server.bat
echo echo 🚀 Starting DNA Website Server... >> start_server.bat
echo call venv\Scripts\activate.bat >> start_server.bat
echo cd src >> start_server.bat
echo python main.py >> start_server.bat
echo pause >> start_server.bat

echo.
echo 🎉 Installation Complete!
echo =======================
echo.
echo To start the website:
echo 1. Double-click: start_server.bat
echo 2. Open browser and go to: http://localhost:5000
echo.
echo To start manually:
echo 1. venv\Scripts\activate.bat
echo 2. cd src
echo 3. python main.py
echo.
echo 📚 Documentation:
echo - README.md - Installation guide
echo - TECHNICAL_DOCUMENTATION.md - Technical details
echo - ALIGNMENT_METHOD_EXPLANATION.md - Algorithm details
echo.
echo 🧬 Enjoy your DNA analysis toolkit!
echo.
pause

