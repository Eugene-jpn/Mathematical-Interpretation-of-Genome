# DNA Tandem Repeat Encoder-Decoder Website

## 🧬 **Complete Local Installation Package**

This package contains everything you need to run the DNA Tandem Repeat analysis website locally on your computer.

## 📋 **Features**

### **Four Main Tabs:**
1. **Encode** - Convert DNA sequences to numerical format using tandem repeat analysis
2. **Decode** - Convert numerical format back to DNA sequences (improved algorithm)
3. **Validate** - Test round-trip encoding/decoding consistency
4. **Mutant** - Compare reference vs sample sequences with mutation analysis

### **Advanced Capabilities:**
- ✅ **Numerical encoding** using formula (A=2, C=3, G=5, T=7)
- ✅ **Tandem repeat detection** and pattern analysis
- ✅ **Two download formats** (complete and fractions-only)
- ✅ **Sequence alignment** with red highlighting for differences
- ✅ **Mutation statistics** and pattern change analysis
- ✅ **Professional UI** with modern design and responsive layout

## 🚀 **Quick Start (5 Minutes)**

### **Prerequisites:**
- Python 3.7 or higher
- pip (Python package installer)

### **Installation Steps:**

1. **Extract this package** to your desired location
2. **Open terminal/command prompt** in the project directory
3. **Run the setup script:**

   **On Windows:**
   ```bash
   setup.bat
   ```

   **On macOS/Linux:**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```

4. **Access the website:**
   - Open your browser
   - Go to: `http://localhost:5000`

## 🔧 **Manual Installation**

If the setup scripts don't work, follow these manual steps:

### **Step 1: Install Dependencies**
```bash
# Create virtual environment (recommended)
python -m venv venv

# Activate virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

# Install required packages
pip install flask flask-cors
```

### **Step 2: Start the Server**
```bash
# Navigate to src directory
cd src

# Run the Flask application
python main.py
```

### **Step 3: Access Website**
- Open browser and go to: `http://localhost:5000`

## 📁 **Project Structure**

```
DNA_Website_Local_Installation/
├── src/                          # Main application code
│   ├── main.py                   # Flask application entry point
│   ├── routes/                   # API endpoints
│   │   └── dna_encoder.py        # DNA encoding/decoding logic
│   └── static/                   # Frontend files
│       └── index.html            # Main website interface
├── venv/                         # Virtual environment (created during setup)
├── requirements.txt              # Python dependencies
├── setup.sh                     # Linux/macOS setup script
├── setup.bat                    # Windows setup script
├── README.md                     # This file
├── TECHNICAL_DOCUMENTATION.md   # Detailed technical docs
└── ALIGNMENT_METHOD_EXPLANATION.md  # Alignment algorithm details
```

## 🧪 **Testing the Website**

### **1. Encode Tab:**
- Input: `AAAACGTTTAACGAACG`
- Expected output: `4 2/9, 3, 5, 3 7/9, 2 2/9, 3, 5, 2 2/9, 3, 5`

### **2. Decode Tab:**
- Input: `4 2/9, 3, 5`
- Expected output: `AAAACG`

### **3. Mutant Tab:**
- Reference: `ATCGATCGATCG`
- Sample: `ATCGATCAATCG`
- Should show mutation at position 8 (G→A) with red highlighting

## 🔧 **Troubleshooting**

### **Common Issues:**

#### **Port 5000 already in use:**
```bash
# Kill process using port 5000
# On Windows:
netstat -ano | findstr :5000
taskkill /PID <PID_NUMBER> /F

# On macOS/Linux:
lsof -ti:5000 | xargs kill -9
```

#### **Python not found:**
- Install Python from: https://python.org/downloads/
- Make sure Python is added to PATH

#### **Permission denied (Linux/macOS):**
```bash
chmod +x setup.sh
sudo ./setup.sh
```

#### **Virtual environment issues:**
```bash
# Remove and recreate virtual environment
rm -rf venv
python -m venv venv
```

## 🌐 **Deployment Options**

### **Local Development:**
- Use the included Flask development server
- Perfect for testing and personal use

### **Production Deployment:**
- Use a production WSGI server like Gunicorn
- Deploy to cloud platforms (Heroku, AWS, etc.)
- Configure proper security settings

## 📚 **Documentation**

- **README.md** - This installation guide
- **TECHNICAL_DOCUMENTATION.md** - Detailed technical specifications
- **ALIGNMENT_METHOD_EXPLANATION.md** - Custom alignment algorithm details

## 🔬 **Algorithm Information**

The sequence alignment method used is a **custom implementation** written specifically for this project:
- **NOT plagiarized** from any external source
- **Simple pairwise alignment** optimized for mutation detection
- **Linear time complexity** O(n) for fast processing
- **Position-based comparison** perfect for DNA mutation analysis

## 📞 **Support**

If you encounter any issues:
1. Check the troubleshooting section above
2. Verify all dependencies are installed correctly
3. Ensure Python 3.7+ is being used
4. Check that port 5000 is available

## 📄 **License**

This project is provided as-is for educational and research purposes.

## 🧬 **Enjoy your DNA analysis toolkit!**

The website provides professional-grade DNA sequence analysis with numerical encoding, mutation detection, and comprehensive visualization tools.

