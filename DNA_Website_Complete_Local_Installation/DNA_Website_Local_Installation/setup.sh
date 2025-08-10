#!/bin/bash

# DNA Website Local Installation Setup Script
# For Linux and macOS systems

echo "🧬 DNA Tandem Repeat Website - Local Installation"
echo "=================================================="

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.7 or higher."
    echo "   Download from: https://python.org/downloads/"
    exit 1
fi

# Check Python version
python_version=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "✅ Python $python_version detected"

# Create virtual environment
echo "📦 Creating virtual environment..."
python3 -m venv venv

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "⬆️  Upgrading pip..."
pip install --upgrade pip

# Install dependencies
echo "📥 Installing dependencies..."
pip install flask flask-cors

# Check if installation was successful
if [ $? -eq 0 ]; then
    echo "✅ Dependencies installed successfully!"
else
    echo "❌ Failed to install dependencies. Please check your internet connection."
    exit 1
fi

# Create a start script
echo "📝 Creating start script..."
cat > start_server.sh << 'EOF'
#!/bin/bash
echo "🚀 Starting DNA Website Server..."
source venv/bin/activate
cd src
python main.py
EOF

chmod +x start_server.sh

echo ""
echo "🎉 Installation Complete!"
echo "======================="
echo ""
echo "To start the website:"
echo "1. Run: ./start_server.sh"
echo "2. Open browser and go to: http://localhost:5000"
echo ""
echo "To start manually:"
echo "1. source venv/bin/activate"
echo "2. cd src"
echo "3. python main.py"
echo ""
echo "📚 Documentation:"
echo "- README.md - Installation guide"
echo "- TECHNICAL_DOCUMENTATION.md - Technical details"
echo "- ALIGNMENT_METHOD_EXPLANATION.md - Algorithm details"
echo ""
echo "🧬 Enjoy your DNA analysis toolkit!"

