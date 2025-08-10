# DNA Website - Technical Documentation

## 🏗️ **Architecture Overview**

### **Technology Stack:**
- **Backend:** Flask (Python)
- **Frontend:** HTML5, CSS3, JavaScript (Vanilla)
- **Styling:** Tailwind CSS classes
- **API:** RESTful endpoints with JSON responses

### **Project Structure:**
```
src/
├── main.py                 # Flask application entry point
├── routes/
│   └── dna_encoder.py      # Core DNA processing logic
└── static/
    └── index.html          # Single-page application
```

## 🧬 **Core Algorithms**

### **1. DNA Encoding Algorithm**
```python
# Nucleotide to number mapping
MAP = {'A': '2', 'T': '7', 'C': '3', 'G': '5'}

# Encoding process:
# 1. Convert nucleotides to numbers
# 2. Detect tandem repeats
# 3. Format as fractions (count numerator/denominator)
```

### **2. Tandem Repeat Detection**
- Identifies repeating patterns in DNA sequences
- Calculates repeat counts and frequencies
- Formats as mathematical fractions

### **3. Custom Sequence Alignment**
```python
def simple_sequence_alignment(ref_seq, sample_seq):
    # Pad sequences to same length
    max_len = max(len(ref_seq), len(sample_seq))
    ref_padded = ref_seq.ljust(max_len, '-')
    sample_padded = sample_seq.ljust(max_len, '-')
    
    # Compare position by position
    for i, (ref_base, sample_base) in enumerate(zip(ref_padded, sample_padded)):
        is_different = ref_base != sample_base
        # Record differences and mutations
```

## 🔌 **API Endpoints**

### **Base URL:** `/api/dna`

#### **1. Encode Sequence**
- **Endpoint:** `POST /encode`
- **Input:** `{"sequence": "ATCGATCG"}`
- **Output:** `{"encoded": "2, 7, 3, 5, 2, 7, 3, 5", "statistics": {...}}`

#### **2. Decode Sequence**
- **Endpoint:** `POST /decode`
- **Input:** `{"encoded_sequence": "2, 7, 3, 5"}`
- **Output:** `{"decoded": "ATCG", "success": true}`

#### **3. Validate Round-trip**
- **Endpoint:** `POST /validate`
- **Input:** `{"sequence": "ATCG"}`
- **Output:** `{"original": "ATCG", "encoded": "2, 7, 3, 5", "decoded": "ATCG", "is_valid": true}`

#### **4. Download Fractions**
- **Endpoint:** `POST /download-fractions`
- **Input:** `{"encoded_sequence": "4 2/9, 3, 5, 2 7/9"}`
- **Output:** `{"fractions_only": "4 2/9, 2 7/9"}`

#### **5. Mutant Comparison**
- **Endpoint:** `POST /mutant-compare`
- **Input:** `{"reference_sequence": "ATCG", "sample_sequences": ["ATCG", "ATCC"]}`
- **Output:** Detailed comparison with alignment, statistics, and mutations

## 🎨 **Frontend Architecture**

### **Single Page Application (SPA)**
- **Tab-based navigation** using JavaScript
- **Dynamic content loading** without page refreshes
- **Responsive design** with mobile support

### **Key JavaScript Functions:**
```javascript
// Tab management
function showTab(tabName)

// DNA encoding
async function encodeSequence()

// DNA decoding  
async function decodeSequence()

// Validation
async function validateSequence()

// Mutant comparison
async function compareMutantSequences()

// Sample management
function addSampleInput()
function removeSampleInput()
```

### **CSS Framework:**
- **Tailwind CSS** utility classes
- **Custom gradients** and animations
- **Responsive grid** layouts
- **Color-coded** statistics and results

## 🔧 **Configuration**

### **Flask Configuration:**
```python
app = Flask(__name__)
app.config['DEBUG'] = True  # Development mode
CORS(app)  # Enable cross-origin requests
```

### **Server Settings:**
- **Host:** `0.0.0.0` (all interfaces)
- **Port:** `5000` (default Flask port)
- **Debug:** Enabled for development

## 📊 **Data Processing**

### **Input Validation:**
- **DNA sequences:** Only A, T, C, G nucleotides allowed
- **Encoded sequences:** Validates numerical format and fractions
- **Sample limits:** Maximum 5 sample sequences for comparison

### **Error Handling:**
- **Backend:** Try-catch blocks with JSON error responses
- **Frontend:** User-friendly alert messages
- **Validation:** Input sanitization and format checking

### **Performance Optimizations:**
- **Linear algorithms** for fast processing
- **Minimal dependencies** for quick startup
- **Efficient string operations** for large sequences

## 🧪 **Testing Strategy**

### **Unit Testing Approach:**
1. **Encoding accuracy** - Verify correct nucleotide to number conversion
2. **Decoding consistency** - Ensure round-trip accuracy
3. **Mutation detection** - Validate alignment and difference highlighting
4. **Edge cases** - Empty sequences, invalid characters, long sequences

### **Integration Testing:**
1. **API endpoints** - Test all REST endpoints
2. **Frontend-backend** communication
3. **File download** functionality
4. **Multi-sample** comparison

## 🔒 **Security Considerations**

### **Input Sanitization:**
- **Nucleotide validation** prevents injection attacks
- **Length limits** prevent memory exhaustion
- **Character filtering** removes potentially harmful input

### **CORS Configuration:**
- **Cross-origin** requests enabled for development
- **Production deployment** should restrict origins

## 📈 **Scalability**

### **Current Limitations:**
- **Single-threaded** Flask development server
- **In-memory** processing (no database)
- **Synchronous** request handling

### **Production Recommendations:**
- **WSGI server** (Gunicorn, uWSGI)
- **Database integration** for large datasets
- **Caching** for frequently accessed sequences
- **Load balancing** for high traffic

## 🔍 **Debugging**

### **Backend Debugging:**
- **Flask debug mode** enabled
- **Console logging** for errors
- **JSON error responses** with details

### **Frontend Debugging:**
- **Browser console** for JavaScript errors
- **Network tab** for API request monitoring
- **Element inspector** for UI issues

## 📝 **Code Quality**

### **Python Standards:**
- **PEP 8** style guidelines
- **Descriptive** function and variable names
- **Comprehensive** error handling
- **Modular** code organization

### **JavaScript Standards:**
- **ES6+** modern syntax
- **Async/await** for API calls
- **Clear** function separation
- **Consistent** naming conventions

## 🚀 **Deployment Guide**

### **Development Deployment:**
```bash
python src/main.py
```

### **Production Deployment:**
```bash
# Using Gunicorn
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 src.main:app
```

### **Docker Deployment:**
```dockerfile
FROM python:3.9-slim
COPY . /app
WORKDIR /app
RUN pip install flask flask-cors
EXPOSE 5000
CMD ["python", "src/main.py"]
```

## 📋 **Dependencies**

### **Python Requirements:**
```
Flask==3.1.1
Flask-CORS==6.0.0
```

### **Frontend Dependencies:**
- **No external libraries** required
- **Vanilla JavaScript** for maximum compatibility
- **Tailwind CSS** classes (CDN or local)

## 🔄 **Version Control**

### **Git Workflow:**
- **Main branch** for stable releases
- **Feature branches** for new development
- **Commit messages** following conventional format

### **Release Management:**
- **Semantic versioning** (v1.0.0, v1.1.0, etc.)
- **Changelog** documentation
- **Tagged releases** for major versions

This technical documentation provides comprehensive details for developers working with or extending the DNA analysis website.

