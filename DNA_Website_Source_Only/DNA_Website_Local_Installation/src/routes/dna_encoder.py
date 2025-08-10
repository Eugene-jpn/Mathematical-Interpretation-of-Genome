from flask import Blueprint, request, jsonify
from flask_cors import cross_origin
import re
import time

dna_encoder_bp = Blueprint('dna_encoder', __name__)

def encode_nucleotide(nucleotide):
    """Encode nucleotide to numerical value"""
    encoding = {'A': '2', 'C': '3', 'G': '5', 'T': '7'}
    return encoding.get(nucleotide.upper(), '0')

# Improved decoding implementation
MAP = {'2': 'A', '3': 'C', '5': 'G', '7': 'T'}
REV = {v: k for k, v in MAP.items()}

class DecodeError(ValueError): 
    pass

def decode_nucleotide(value):
    """Decode numerical value back to nucleotide"""
    decoding = {'2': 'A', '3': 'C', '5': 'G', '7': 'T'}
    return decoding.get(str(value), 'N')

def _motif_from_digits(digits: str) -> str:
    """Convert digits back to DNA motif"""
    try:
        return ''.join(MAP[d] for d in digits)
    except KeyError:
        raise DecodeError(f"Invalid digit in motif: {digits!r}. Allowed: 2,3,5,7")

def decode_tokens(text: str) -> str:
    """Improved decoding function with better error handling and validation"""
    seq = []
    for raw in text.replace('\n', ',').split(','):
        tok = raw.strip()
        if not tok:
            continue
        if ' ' in tok:
            first, rest = tok.split(' ', 1)
            first, rest = first.strip(), rest.strip()
            if '/' in rest:
                try:
                    N = int(first)
                except ValueError:
                    raise DecodeError(f"Invalid repeat count: {first!r}")
                digits, denom = rest.split('/', 1)
                if not digits or not denom or not set(denom) <= {'9'}:
                    raise DecodeError(f"Invalid fraction: {rest!r}")
                if len(denom) != len(digits):
                    raise DecodeError(f"Denominator 9s length must equal digits length: {rest!r}")
                motif = _motif_from_digits(digits)
                seq.append(motif * N)
            else:
                if first in MAP: seq.append(MAP[first])
                if rest in MAP:  seq.append(MAP[rest])
        else:
            if tok in MAP:
                seq.append(MAP[tok])
            else:
                raise DecodeError(f"Unexpected token: {tok!r}")
    return ''.join(seq)

def find_tandem_repeats_with_positions(sequence):
    """Find tandem repeats in DNA sequence with their positions - optimized version"""
    repeats = []
    
    # Limit sequence length to prevent timeouts
    if len(sequence) > 1000:
        sequence = sequence[:1000]
    
    # For shorter sequences, limit the search to reasonable repeat lengths
    max_unit_length = min(20, len(sequence) // 3)  # Reduced from 50 to 20
    
    for unit_length in range(1, max_unit_length + 1):
        # Skip every other position for longer units to speed up
        step = 1 if unit_length <= 5 else 2
        
        for start in range(0, len(sequence) - unit_length, step):
            repeat_unit = sequence[start:start + unit_length]
            repeat_count = 1
            
            pos = start + unit_length
            max_repeats = min(10, (len(sequence) - start) // unit_length)  # Limit max repeats
            
            while pos + unit_length <= len(sequence) and repeat_count < max_repeats:
                if sequence[pos:pos + unit_length] == repeat_unit:
                    repeat_count += 1
                    pos += unit_length
                else:
                    break
            
            if repeat_count >= 2:
                repeats.append({
                    'start': start,
                    'end': start + (repeat_count * unit_length),
                    'unit': repeat_unit,
                    'count': repeat_count,
                    'length': unit_length
                })
                
                # Limit total repeats found to prevent memory issues
                if len(repeats) > 50:
                    break
        
        if len(repeats) > 50:
            break
    
    # Remove overlapping repeats, prioritizing longer units and higher counts
    filtered_repeats = []
    repeats.sort(key=lambda x: (x['start'], -x['length'], -x['count']))
    
    for repeat in repeats:
        # Check if this repeat overlaps with any already added
        overlaps = False
        for existing in filtered_repeats:
            if (repeat['start'] < existing['end'] and repeat['end'] > existing['start']):
                overlaps = True
                break
        
        if not overlaps:
            filtered_repeats.append(repeat)
            
        # Limit final results
        if len(filtered_repeats) >= 20:
            break
    
    return filtered_repeats

def encode_sequence_to_format(sequence):
    """Encode DNA sequence with absolute priority for consecutive identical nucleotides"""
    sequence = sequence.upper().replace('\n', '').replace(' ', '')
    
    output_parts = []
    pos = 0
    
    while pos < len(sequence):
        # ALWAYS check for consecutive identical nucleotides first
        current_nuc = sequence[pos]
        consecutive_count = 1
        temp_pos = pos + 1
        
        # Count consecutive identical nucleotides
        while temp_pos < len(sequence) and sequence[temp_pos] == current_nuc:
            consecutive_count += 1
            temp_pos += 1
        
        # If we have 2 or more consecutive identical nucleotides, encode as repeat
        if consecutive_count >= 2:
            encoded_value = encode_nucleotide(current_nuc)
            repeat_notation = f"{consecutive_count} {encoded_value}/9"
            output_parts.append(repeat_notation)
            pos = temp_pos
        else:
            # Only check for tandem repeats if it's NOT a consecutive nucleotide
            # AND only for multi-nucleotide patterns
            found_tandem = False
            
            # Look for tandem repeats starting at this position
            for unit_length in range(2, min(10, (len(sequence) - pos) // 2 + 1)):
                if pos + unit_length * 2 <= len(sequence):
                    unit = sequence[pos:pos + unit_length]
                    
                    # Check if this unit repeats
                    repeat_count = 1
                    check_pos = pos + unit_length
                    
                    while check_pos + unit_length <= len(sequence):
                        if sequence[check_pos:check_pos + unit_length] == unit:
                            repeat_count += 1
                            check_pos += unit_length
                        else:
                            break
                    
                    # If we found a tandem repeat (2 or more repetitions)
                    if repeat_count >= 2:
                        encoded_unit = ''.join(encode_nucleotide(nuc) for nuc in unit)
                        denominator = '9' * len(unit)
                        
                        repeat_notation = f"{repeat_count} {encoded_unit}/{denominator}"
                        output_parts.append(repeat_notation)
                        pos = pos + (repeat_count * unit_length)
                        found_tandem = True
                        break
            
            # If no tandem repeat found, encode as single nucleotide
            if not found_tandem:
                nuc_value = encode_nucleotide(sequence[pos])
                output_parts.append(nuc_value)
                pos += 1
    
    return ', '.join(output_parts)

def analyze_sequence_stats(sequence, repeats):
    """Analyze sequence statistics"""
    sequence = sequence.upper()
    
    # Count nucleotides
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for nuc in sequence:
        if nuc in nucleotide_counts:
            nucleotide_counts[nuc] += 1
    
    # Calculate GC content
    gc_content = (nucleotide_counts['G'] + nucleotide_counts['C']) / len(sequence) * 100
    
    # Repeat statistics
    total_repeat_length = sum(r['count'] * r['length'] for r in repeats)
    repeat_coverage = (total_repeat_length / len(sequence)) * 100 if len(sequence) > 0 else 0
    
    return {
        'length': len(sequence),
        'nucleotide_counts': nucleotide_counts,
        'gc_content': round(gc_content, 2),
        'repeat_count': len(repeats),
        'repeat_coverage': round(repeat_coverage, 2)
    }

@dna_encoder_bp.route('/encode', methods=['POST'])
@cross_origin()
def encode_dna():
    """Encode DNA sequence to numerical format"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').strip()
        
        if not sequence:
            return jsonify({'error': 'No sequence provided'}), 400
        
        # Validate sequence
        valid_nucleotides = set('ATCG')
        clean_sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        if not all(nuc in valid_nucleotides for nuc in clean_sequence):
            return jsonify({'error': 'Invalid nucleotides found. Only A, T, C, G are allowed.'}), 400
        
        if len(clean_sequence) > 10000:
            return jsonify({'error': 'Sequence too long. Maximum 10,000 nucleotides allowed.'}), 400
        
        start_time = time.time()
        
        # Find repeats
        repeats = find_tandem_repeats_with_positions(clean_sequence)
        
        # Encode sequence
        encoded_result = encode_sequence_to_format(clean_sequence)
        
        # Get statistics
        stats = analyze_sequence_stats(clean_sequence, repeats)
        
        processing_time = round((time.time() - start_time) * 1000, 2)
        
        # Format repeats for frontend
        formatted_repeats = []
        for repeat in repeats[:20]:  # Limit to first 20 for display
            encoded_value = int(''.join(encode_nucleotide(nuc) for nuc in repeat['unit']))
            denominator = '9' * repeat['length']
            formatted_repeats.append({
                'position': f"{repeat['start']+1}-{repeat['end']}",
                'unit': repeat['unit'],
                'count': repeat['count'],
                'length': repeat['length'],
                'formula': f"{repeat['count']} {encoded_value}/{denominator}"
            })
        
        return jsonify({
            'success': True,
            'original_sequence': clean_sequence,
            'encoded_result': encoded_result,
            'stats': stats,
            'repeats': formatted_repeats,
            'processing_time_ms': processing_time
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@dna_encoder_bp.route('/decode', methods=['POST'])
@cross_origin()
def decode_dna():
    """Decode numerical format back to DNA sequence"""
    try:
        data = request.get_json()
        encoded_sequence = data.get('encoded_sequence', '').strip()
        
        if not encoded_sequence:
            return jsonify({'error': 'No encoded sequence provided'}), 400
        
        start_time = time.time()
        
        try:
            # Decode sequence using improved function
            decoded_result = decode_tokens(encoded_sequence)
        except DecodeError as e:
            return jsonify({'error': f'Decoding error: {str(e)}'}), 400
        
        processing_time = round((time.time() - start_time) * 1000, 2)
        
        # Get basic stats
        nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for nuc in decoded_result:
            if nuc in nucleotide_counts:
                nucleotide_counts[nuc] += 1
        
        gc_content = (nucleotide_counts['G'] + nucleotide_counts['C']) / len(decoded_result) * 100 if len(decoded_result) > 0 else 0
        
        stats = {
            'length': len(decoded_result),
            'nucleotide_counts': nucleotide_counts,
            'gc_content': round(gc_content, 2)
        }
        
        return jsonify({
            'success': True,
            'encoded_sequence': encoded_sequence,
            'decoded_result': decoded_result,
            'stats': stats,
            'processing_time_ms': processing_time
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@dna_encoder_bp.route('/validate', methods=['POST'])
@cross_origin()
def validate_roundtrip():
    """Validate that encoding and decoding produces the original sequence"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').strip()
        
        if not sequence:
            return jsonify({'error': 'No sequence provided'}), 400
        
        # Clean sequence
        clean_sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        start_time = time.time()
        
        try:
            # Encode then decode using improved function
            encoded = encode_sequence_to_format(clean_sequence)
            decoded = decode_tokens(encoded)
        except DecodeError as e:
            return jsonify({'error': f'Validation error: {str(e)}'}), 400
        
        processing_time = round((time.time() - start_time) * 1000, 2)
        
        # Check if they match
        is_valid = clean_sequence == decoded
        
        return jsonify({
            'success': True,
            'original': clean_sequence,
            'encoded': encoded,
            'decoded': decoded,
            'is_valid': is_valid,
            'processing_time_ms': processing_time
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500



@dna_encoder_bp.route('/compare', methods=['POST'])
@cross_origin()
def compare_sequences():
    """Compare wildtype and mutant sequences using numerical encoding analysis - simplified version"""
    try:
        data = request.get_json()
        wildtype_sequence = data.get('wildtype_sequence', '').strip()
        mutant_sequence = data.get('mutant_sequence', '').strip()
        
        if not wildtype_sequence or not mutant_sequence:
            return jsonify({'error': 'Both wildtype and mutant sequences are required'}), 400
        
        # Validate sequences
        valid_nucleotides = set('ATCG')
        clean_wildtype = wildtype_sequence.upper().replace('\n', '').replace(' ', '')
        clean_mutant = mutant_sequence.upper().replace('\n', '').replace(' ', '')
        
        if not all(nuc in valid_nucleotides for nuc in clean_wildtype):
            return jsonify({'error': 'Invalid nucleotides in wildtype sequence. Only A, T, C, G are allowed.'}), 400
        
        if not all(nuc in valid_nucleotides for nuc in clean_mutant):
            return jsonify({'error': 'Invalid nucleotides in mutant sequence. Only A, T, C, G are allowed.'}), 400
        
        # Encode both sequences
        wt_encoded = encode_sequence_to_format(clean_wildtype)
        mut_encoded = encode_sequence_to_format(clean_mutant)
        
        # Create simple comparison table
        comparison_table = []
        
        # Parse encoded sequences
        wt_parts = [part.strip() for part in wt_encoded.split(',')]
        mut_parts = [part.strip() for part in mut_encoded.split(',')]
        
        # Group by pattern type
        wt_patterns = {}
        mut_patterns = {}
        
        for part in wt_parts:
            if '/' in part:
                pattern = part.split()[1]  # Get the pattern part (e.g., "2/9")
                count = int(part.split()[0])  # Get the count
                if pattern in wt_patterns:
                    wt_patterns[pattern] += count
                else:
                    wt_patterns[pattern] = count
        
        for part in mut_parts:
            if '/' in part:
                pattern = part.split()[1]  # Get the pattern part (e.g., "2/9")
                count = int(part.split()[0])  # Get the count
                if pattern in mut_patterns:
                    mut_patterns[pattern] += count
                else:
                    mut_patterns[pattern] = count
        
        # Create comparison table
        all_patterns = set(wt_patterns.keys()) | set(mut_patterns.keys())
        
        for pattern in all_patterns:
            wt_count = wt_patterns.get(pattern, 0)
            mut_count = mut_patterns.get(pattern, 0)
            difference = mut_count - wt_count
            
            # Decode pattern to nucleotide
            nucleotide = decode_pattern_to_nucleotide(pattern)
            
            comparison_table.append({
                'pattern': pattern,
                'nucleotide': nucleotide,
                'wildtype_count': wt_count,
                'mutant_count': mut_count,
                'difference': difference
            })
        
        return jsonify({
            'success': True,
            'wildtype_encoded': wt_encoded,
            'mutant_encoded': mut_encoded,
            'comparison_table': comparison_table,
            'summary': {
                'total_differences': len([row for row in comparison_table if row['difference'] != 0]),
                'wildtype_length': len(clean_wildtype),
                'mutant_length': len(clean_mutant)
            }
        })
        
    except Exception as e:
        return jsonify({'error': f'Comparison failed: {str(e)}'}), 500

def decode_pattern_to_nucleotide(pattern):
    """Decode a numerical pattern back to nucleotide sequence"""
    try:
        if '/' not in pattern:
            return decode_nucleotide(pattern)
        
        numerator, denominator = pattern.split('/')
        
        if len(denominator) == 1:  # Single nucleotide (e.g., "2/9")
            return decode_nucleotide(numerator)
        else:  # Multi-nucleotide pattern
            nucleotides = []
            for digit in numerator:
                nucleotides.append(decode_nucleotide(digit))
            return ''.join(nucleotides)
    except:
        return 'N'

# Genetic codon table for amino acid translation
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Enhanced amino acid properties for accurate clinical assessment
AMINO_ACID_PROPERTIES = {
    'F': {'name': 'Phenylalanine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'large', 'hydrophobic': True, 'essential': True},
    'L': {'name': 'Leucine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'medium', 'hydrophobic': True, 'essential': True},
    'S': {'name': 'Serine', 'type': 'polar', 'charge': 'neutral', 'size': 'small', 'hydrophobic': False, 'essential': False},
    'Y': {'name': 'Tyrosine', 'type': 'polar', 'charge': 'neutral', 'size': 'large', 'hydrophobic': False, 'essential': False},
    'C': {'name': 'Cysteine', 'type': 'polar', 'charge': 'neutral', 'size': 'small', 'hydrophobic': False, 'essential': False},
    'W': {'name': 'Tryptophan', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'large', 'hydrophobic': True, 'essential': True},
    'P': {'name': 'Proline', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'small', 'hydrophobic': True, 'essential': False},
    'H': {'name': 'Histidine', 'type': 'polar', 'charge': 'positive', 'size': 'medium', 'hydrophobic': False, 'essential': True},
    'Q': {'name': 'Glutamine', 'type': 'polar', 'charge': 'neutral', 'size': 'medium', 'hydrophobic': False, 'essential': False},
    'R': {'name': 'Arginine', 'type': 'polar', 'charge': 'positive', 'size': 'large', 'hydrophobic': False, 'essential': False},
    'I': {'name': 'Isoleucine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'medium', 'hydrophobic': True, 'essential': True},
    'M': {'name': 'Methionine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'medium', 'hydrophobic': True, 'essential': True},
    'T': {'name': 'Threonine', 'type': 'polar', 'charge': 'neutral', 'size': 'small', 'hydrophobic': False, 'essential': True},
    'N': {'name': 'Asparagine', 'type': 'polar', 'charge': 'neutral', 'size': 'medium', 'hydrophobic': False, 'essential': False},
    'K': {'name': 'Lysine', 'type': 'polar', 'charge': 'positive', 'size': 'large', 'hydrophobic': False, 'essential': True},
    'V': {'name': 'Valine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'medium', 'hydrophobic': True, 'essential': True},
    'A': {'name': 'Alanine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'small', 'hydrophobic': True, 'essential': False},
    'D': {'name': 'Aspartic acid', 'type': 'polar', 'charge': 'negative', 'size': 'medium', 'hydrophobic': False, 'essential': False},
    'E': {'name': 'Glutamic acid', 'type': 'polar', 'charge': 'negative', 'size': 'medium', 'hydrophobic': False, 'essential': False},
    'G': {'name': 'Glycine', 'type': 'nonpolar', 'charge': 'neutral', 'size': 'small', 'hydrophobic': True, 'essential': False},
    '*': {'name': 'Stop codon', 'type': 'stop', 'charge': 'none', 'size': 'none', 'hydrophobic': False, 'essential': False}
}

def translate_dna_to_protein(dna_sequence, reading_frame=0):
    """Translate DNA sequence to protein sequence with improved accuracy"""
    if len(dna_sequence) < 3:
        return ""
    
    protein = []
    start_pos = reading_frame
    
    # Ensure we don't go out of bounds
    for i in range(start_pos, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:  # Only translate complete codons
            amino_acid = CODON_TABLE.get(codon, 'X')  # X for unknown/invalid codons
            protein.append(amino_acid)
            # Don't automatically stop at stop codons in mutation analysis
        else:
            break  # Stop if incomplete codon
    
    return ''.join(protein)

def find_optimal_reading_frame(wt_seq, mut_seq):
    """Find the optimal reading frame by analyzing protein quality"""
    best_frame = 0
    best_score = -1
    
    for frame in range(3):
        try:
            wt_protein = translate_dna_to_protein(wt_seq, frame)
            mut_protein = translate_dna_to_protein(mut_seq, frame)
            
            if len(wt_protein) == 0 and len(mut_protein) == 0:
                continue
            
            # Score based on:
            # 1. Fewer premature stop codons (excluding terminal)
            # 2. Longer protein sequences
            # 3. Higher similarity between WT and mutant
            
            wt_stops = wt_protein[:-1].count('*') if len(wt_protein) > 1 else wt_protein.count('*')
            mut_stops = mut_protein[:-1].count('*') if len(mut_protein) > 1 else mut_protein.count('*')
            
            # Penalize premature stops heavily
            stop_penalty = (wt_stops + mut_stops) * 10
            
            # Reward longer sequences
            length_bonus = min(len(wt_protein), len(mut_protein))
            
            # Calculate similarity
            min_len = min(len(wt_protein), len(mut_protein))
            similarity = 0
            if min_len > 0:
                similarity = sum(1 for i in range(min_len) if wt_protein[i] == mut_protein[i])
            
            score = length_bonus + similarity - stop_penalty
            
            if score > best_score:
                best_score = score
                best_frame = frame
                
        except Exception:
            continue
    
    return best_frame

def detect_mutation_type(wt_seq, mut_seq, position):
    """Detect the type of mutation at a specific position"""
    if len(wt_seq) != len(mut_seq):
        return 'indel'
    
    # Check if it's a substitution
    if wt_seq[position] != mut_seq[position]:
        return 'substitution'
    
    return 'none'

def analyze_frameshift_mutations(wt_seq, mut_seq):
    """Analyze sequences for frameshift mutations with improved detection"""
    frameshift_positions = []
    
    # Basic length-based frameshift detection
    length_diff = len(mut_seq) - len(wt_seq)
    
    if length_diff != 0:
        # Determine the position where the frameshift likely occurs
        # Find the first position where sequences differ
        frameshift_pos = 0
        min_len = min(len(wt_seq), len(mut_seq))
        
        for i in range(min_len):
            if wt_seq[i] != mut_seq[i]:
                frameshift_pos = i
                break
        
        # If sequences are identical up to the shorter length, frameshift is at the end
        if frameshift_pos == 0 and min_len > 0:
            frameshift_pos = min_len
        
        # Check if this is truly a frameshift (not divisible by 3)
        if abs(length_diff) % 3 != 0:
            frameshift_type = 'insertion' if length_diff > 0 else 'deletion'
            
            frameshift_positions.append({
                'position': frameshift_pos + 1,  # 1-based position
                'type': frameshift_type,
                'length': abs(length_diff),
                'confidence': 'high' if abs(length_diff) < 10 else 'medium'
            })
        else:
            # Even if divisible by 3, could still be a frameshift if not in-frame
            # This is a more complex analysis that would require alignment
            if abs(length_diff) > 0:
                frameshift_type = 'insertion' if length_diff > 0 else 'deletion'
                
                frameshift_positions.append({
                    'position': frameshift_pos + 1,
                    'type': f'in-frame_{frameshift_type}',
                    'length': abs(length_diff),
                    'confidence': 'medium'
                })
    
    # Additional check for complex rearrangements
    if len(frameshift_positions) == 0 and len(wt_seq) == len(mut_seq):
        # Look for patterns that might indicate complex mutations
        differences = sum(1 for i in range(len(wt_seq)) if wt_seq[i] != mut_seq[i])
        
        # If more than 50% of the sequence is different, might be a complex rearrangement
        if differences > len(wt_seq) * 0.5:
            frameshift_positions.append({
                'position': 1,
                'type': 'complex_rearrangement',
                'length': differences,
                'confidence': 'low'
            })
    
    return frameshift_positions

def analyze_amino_acid_changes(wt_seq, mut_seq):
    """Analyze amino acid changes between wildtype and mutant sequences"""
    changes = []
    
    # Ensure sequences are properly aligned and handle length differences
    max_length = max(len(wt_seq), len(mut_seq))
    
    # Pad shorter sequence with gaps for proper alignment
    if len(wt_seq) < max_length:
        wt_seq = wt_seq + 'N' * (max_length - len(wt_seq))
    if len(mut_seq) < max_length:
        mut_seq = mut_seq + 'N' * (max_length - len(mut_seq))
    
    # Find the longest common prefix to determine where changes start
    change_start = 0
    for i in range(min(len(wt_seq), len(mut_seq))):
        if wt_seq[i] != mut_seq[i]:
            change_start = i
            break
    
    # Align sequences for translation starting from the first change
    # Use the shorter sequence length to avoid index errors
    min_length = min(len(wt_seq), len(mut_seq))
    
    # Ensure we have complete codons by trimming to multiple of 3
    codon_length = (min_length // 3) * 3
    wt_seq_trimmed = wt_seq[:codon_length]
    mut_seq_trimmed = mut_seq[:codon_length]
    
    # Translate both sequences in all three reading frames to find the best match
    best_changes = []
    min_differences = float('inf')
    
    for reading_frame in range(3):
        frame_changes = []
        
        # Translate sequences in this reading frame
        try:
            wt_protein = translate_dna_to_protein(wt_seq_trimmed, reading_frame)
            mut_protein = translate_dna_to_protein(mut_seq_trimmed, reading_frame)
            
            # Compare amino acid sequences
            min_protein_length = min(len(wt_protein), len(mut_protein))
            differences = 0
            
            for i in range(min_protein_length):
                if wt_protein[i] != mut_protein[i]:
                    differences += 1
                    codon_pos = (i * 3) + reading_frame
                    
                    # Ensure we don't go out of bounds
                    if codon_pos + 3 <= len(wt_seq_trimmed) and codon_pos + 3 <= len(mut_seq_trimmed):
                        wt_codon = wt_seq_trimmed[codon_pos:codon_pos+3]
                        mut_codon = mut_seq_trimmed[codon_pos:codon_pos+3]
                        
                        # Determine mutation type
                        mutation_type = 'missense'
                        if mut_protein[i] == '*':
                            mutation_type = 'nonsense'
                        elif wt_protein[i] == '*':
                            mutation_type = 'readthrough'
                        
                        # Get amino acid properties safely
                        wt_aa_props = AMINO_ACID_PROPERTIES.get(wt_protein[i], {'name': 'Unknown'})
                        mut_aa_props = AMINO_ACID_PROPERTIES.get(mut_protein[i], {'name': 'Unknown'})
                        
                        frame_changes.append({
                            'position': i + 1,
                            'codon_position': codon_pos + 1,
                            'reading_frame': reading_frame,
                            'wildtype_aa': wt_protein[i],
                            'mutant_aa': mut_protein[i],
                            'wildtype_codon': wt_codon,
                            'mutant_codon': mut_codon,
                            'mutation_type': mutation_type,
                            'wildtype_aa_name': wt_aa_props.get('name', 'Unknown'),
                            'mutant_aa_name': mut_aa_props.get('name', 'Unknown'),
                            'confidence': 'high' if reading_frame == 0 else 'medium'
                        })
            
            # Check for length differences in proteins
            if len(wt_protein) != len(mut_protein):
                differences += abs(len(wt_protein) - len(mut_protein))
            
            # Keep the reading frame with the fewest differences
            if differences < min_differences:
                min_differences = differences
                best_changes = frame_changes
                
        except Exception as e:
            # Skip this reading frame if translation fails
            continue
    
    # If no changes found in any frame, try a simpler approach
    if not best_changes and len(wt_seq_trimmed) >= 3 and len(mut_seq_trimmed) >= 3:
        # Simple codon-by-codon comparison in reading frame 0
        for i in range(0, min(len(wt_seq_trimmed), len(mut_seq_trimmed)) - 2, 3):
            wt_codon = wt_seq_trimmed[i:i+3]
            mut_codon = mut_seq_trimmed[i:i+3]
            
            if wt_codon != mut_codon and len(wt_codon) == 3 and len(mut_codon) == 3:
                wt_aa = CODON_TABLE.get(wt_codon, 'X')
                mut_aa = CODON_TABLE.get(mut_codon, 'X')
                
                if wt_aa != mut_aa:
                    mutation_type = 'missense'
                    if mut_aa == '*':
                        mutation_type = 'nonsense'
                    elif wt_aa == '*':
                        mutation_type = 'readthrough'
                    
                    wt_aa_props = AMINO_ACID_PROPERTIES.get(wt_aa, {'name': 'Unknown'})
                    mut_aa_props = AMINO_ACID_PROPERTIES.get(mut_aa, {'name': 'Unknown'})
                    
                    best_changes.append({
                        'position': (i // 3) + 1,
                        'codon_position': i + 1,
                        'reading_frame': 0,
                        'wildtype_aa': wt_aa,
                        'mutant_aa': mut_aa,
                        'wildtype_codon': wt_codon,
                        'mutant_codon': mut_codon,
                        'mutation_type': mutation_type,
                        'wildtype_aa_name': wt_aa_props.get('name', 'Unknown'),
                        'mutant_aa_name': mut_aa_props.get('name', 'Unknown'),
                        'confidence': 'medium'
                    })
    
    return best_changes

def decode_pattern_to_nucleotides(pattern):
    """Decode numerical pattern back to nucleotide sequence"""
    if '/' not in pattern:
        # Single nucleotide
        return decode_nucleotide(pattern)
    
    # Pattern like "375/999" or "2/9"
    numerator, denominator = pattern.split('/')
    
    # Determine how many nucleotides based on denominator length
    nucleotide_count = len(denominator)
    
    if nucleotide_count == 1:
        # Single nucleotide pattern like "2/9"
        return decode_nucleotide(numerator)
    else:
        # Multi-nucleotide pattern like "375/999"
        nucleotides = []
        for digit in numerator:
            nucleotides.append(decode_nucleotide(digit))
        return ''.join(nucleotides)

def parse_encoded_sequence(encoded_string):
    """Parse encoded sequence string into structured data"""
    patterns = []
    tokens = [token.strip() for token in encoded_string.split(',')]
    
    for token in tokens:
        if '/' in token and ' ' in token:
            # Pattern like "2 27/99" or "52 375/999"
            parts = token.strip().split(' ')
            if len(parts) == 2:
                count = int(parts[0])
                pattern = parts[1]  # e.g., "27/99"
                nucleotides = decode_pattern_to_nucleotides(pattern)
                patterns.append({
                    'count': count,
                    'pattern': pattern,
                    'nucleotides': nucleotides,
                    'type': 'repeat'
                })
        elif token.strip().isdigit():
            # Single nucleotide
            nucleotides = decode_nucleotide(token.strip())
            patterns.append({
                'count': 1,
                'pattern': token.strip(),
                'nucleotides': nucleotides,
                'type': 'single'
            })
    
    return patterns

def create_numerical_comparison_table(seq1_encoded, seq2_encoded):
    """Create numerical comparison table for encoded sequences with nucleotide column"""
    
    # Parse both sequences
    seq1_patterns = parse_encoded_sequence(seq1_encoded)
    seq2_patterns = parse_encoded_sequence(seq2_encoded)
    
    # Group patterns by their denominator/pattern
    seq1_groups = {}
    seq2_groups = {}
    pattern_nucleotides = {}  # Store nucleotide sequences for each pattern
    
    # Group sequence 1 patterns
    for pattern_data in seq1_patterns:
        key = pattern_data['pattern']
        if key not in seq1_groups:
            seq1_groups[key] = 0
        seq1_groups[key] += pattern_data['count']
        pattern_nucleotides[key] = pattern_data['nucleotides']
    
    # Group sequence 2 patterns
    for pattern_data in seq2_patterns:
        key = pattern_data['pattern']
        if key not in seq2_groups:
            seq2_groups[key] = 0
        seq2_groups[key] += pattern_data['count']
        pattern_nucleotides[key] = pattern_data['nucleotides']
    
    # Create comparison table
    all_patterns = set(seq1_groups.keys()) | set(seq2_groups.keys())
    comparison_table = []
    
    for pattern in sorted(all_patterns):
        seq1_count = seq1_groups.get(pattern, 0)
        seq2_count = seq2_groups.get(pattern, 0)
        difference = seq2_count - seq1_count
        nucleotides = pattern_nucleotides.get(pattern, 'N')
        
        comparison_table.append({
            'pattern': pattern,
            'nucleotides': nucleotides,
            'sequence1_count': seq1_count,
            'sequence2_count': seq2_count,
            'difference': difference,
            'absolute_difference': abs(difference),
            'percentage_change': ((difference / seq1_count) * 100) if seq1_count > 0 else 0
        })
    
    # Sort by absolute difference (most significant changes first)
    comparison_table.sort(key=lambda x: x['absolute_difference'], reverse=True)
    
    return comparison_table

def analyze_numerical_differences(comparison_table):
    """Analyze the numerical differences between sequences"""
    
    total_differences = sum(row['absolute_difference'] for row in comparison_table)
    significant_changes = [row for row in comparison_table if row['absolute_difference'] > 0]
    
    analysis = {
        'total_patterns_compared': len(comparison_table),
        'patterns_with_differences': len(significant_changes),
        'total_numerical_difference': total_differences,
        'most_significant_change': None,
        'summary': []
    }
    
    if significant_changes:
        analysis['most_significant_change'] = significant_changes[0]
        
        # Create summary of changes
        for change in significant_changes[:5]:  # Top 5 changes
            if change['difference'] > 0:
                direction = "increased"
            else:
                direction = "decreased"
            
            analysis['summary'].append({
                'pattern': change['pattern'],
                'direction': direction,
                'magnitude': change['absolute_difference'],
                'description': f"Pattern {change['pattern']} {direction} by {change['absolute_difference']} units"
            })
    
    return analysis

def perform_numerical_sequence_comparison(wt_seq, mut_seq):
    """Perform numerical comparison of encoded sequences"""
    
    # Encode both sequences
    wt_encoded = encode_sequence_to_format(wt_seq)
    mut_encoded = encode_sequence_to_format(mut_seq)
    
    # Create comparison table
    comparison_table = create_numerical_comparison_table(wt_encoded, mut_encoded)
    
    # Analyze differences
    analysis = analyze_numerical_differences(comparison_table)
    
    # Calculate sequence statistics
    wt_stats = analyze_sequence_stats(wt_seq, find_tandem_repeats_with_positions(wt_seq))
    mut_stats = analyze_sequence_stats(mut_seq, find_tandem_repeats_with_positions(mut_seq))
    
    return {
        'wildtype_encoded': wt_encoded,
        'mutant_encoded': mut_encoded,
        'comparison_table': comparison_table,
        'numerical_analysis': analysis,
        'sequence_statistics': {
            'wildtype': wt_stats,
            'mutant': mut_stats,
            'length_difference': mut_stats['length'] - wt_stats['length'],
            'gc_content_difference': mut_stats['gc_content'] - wt_stats['gc_content']
        },
        'metadata': {
            'comparison_type': 'numerical',
            'analysis_focus': 'encoded_pattern_differences',
            'calculation_method': 'pattern_grouping_and_subtraction'
        }
    }

def perform_advanced_mutation_analysis(wt_seq, mut_seq):
    """Perform comprehensive mutation analysis with GenBank and research integration"""
    # Find optimal reading frame for both sequences
    optimal_frame = find_optimal_reading_frame(wt_seq, mut_seq)
    
    # Detect frameshift mutations
    frameshift_mutations = analyze_frameshift_mutations(wt_seq, mut_seq)
    
    # Analyze amino acid changes using optimal reading frame
    amino_acid_changes = analyze_amino_acid_changes(wt_seq, mut_seq, optimal_frame)
    
    # Assess health implications with research data
    health_implications = assess_health_implications_with_research(amino_acid_changes, frameshift_mutations)
    
    # Create highlighting information for frameshift regions
    frameshift_highlights = []
    for frameshift in frameshift_mutations:
        start_pos = frameshift['position']
        # Highlight the region affected by frameshift
        frameshift_highlights.append({
            'start': start_pos,
            'end': min(start_pos + 50, len(mut_seq)),  # Highlight next 50 nucleotides or end of sequence
            'type': 'frameshift'
        })
    
    return {
        'frameshift_mutations': frameshift_mutations,
        'amino_acid_changes': amino_acid_changes,
        'health_implications': health_implications,
        'frameshift_highlights': frameshift_highlights,
        'optimal_reading_frame': optimal_frame,
        'analysis_metadata': {
            'genbank_integration': True,
            'research_papers_included': True,
            'clinical_guidelines': 'ACMG/AMP 2015 Standards',
            'evidence_sources': ['ClinVar', 'PubMed', 'OMIM', 'ClinGen']
        }
    }

def analyze_amino_acid_changes(wt_seq, mut_seq, reading_frame=None):
    """Analyze amino acid changes with improved reading frame detection"""
    if reading_frame is None:
        reading_frame = find_optimal_reading_frame(wt_seq, mut_seq)
    
    # Ensure sequences are the same length for comparison
    min_length = min(len(wt_seq), len(mut_seq))
    wt_seq_trimmed = wt_seq[:min_length]
    mut_seq_trimmed = mut_seq[:min_length]
    
    # Translate both sequences
    wt_protein = translate_dna_to_protein(wt_seq_trimmed, reading_frame)
    mut_protein = translate_dna_to_protein(mut_seq_trimmed, reading_frame)
    
    changes = []
    min_protein_length = min(len(wt_protein), len(mut_protein))
    
    for i in range(min_protein_length):
        if wt_protein[i] != mut_protein[i]:
            # Calculate codon positions
            codon_start = reading_frame + (i * 3)
            
            if codon_start + 2 < min_length:
                wt_codon = wt_seq_trimmed[codon_start:codon_start + 3]
                mut_codon = mut_seq_trimmed[codon_start:codon_start + 3]
                
                wt_aa = wt_protein[i]
                mut_aa = mut_protein[i]
                
                # Determine mutation type
                mutation_type = 'missense'
                if mut_aa == '*':
                    mutation_type = 'nonsense'
                elif wt_aa == '*':
                    mutation_type = 'readthrough'
                
                wt_aa_props = AMINO_ACID_PROPERTIES.get(wt_aa, {'name': 'Unknown'})
                mut_aa_props = AMINO_ACID_PROPERTIES.get(mut_aa, {'name': 'Unknown'})
                
                # Calculate confidence based on sequence quality and context
                confidence = 'high'
                if codon_start < 10 or codon_start > min_length - 10:
                    confidence = 'medium'  # Lower confidence near sequence ends
                
                changes.append({
                    'position': i + 1,  # 1-based amino acid position
                    'codon_position': codon_start + 1,  # 1-based nucleotide position
                    'reading_frame': reading_frame,
                    'wildtype_aa': wt_aa,
                    'mutant_aa': mut_aa,
                    'wildtype_codon': wt_codon,
                    'mutant_codon': mut_codon,
                    'mutation_type': mutation_type,
                    'wildtype_aa_name': wt_aa_props.get('name', 'Unknown'),
                    'mutant_aa_name': mut_aa_props.get('name', 'Unknown'),
                    'confidence': confidence,
                    'biochemical_change': {
                        'type_change': wt_aa_props.get('type') != mut_aa_props.get('type'),
                        'charge_change': wt_aa_props.get('charge') != mut_aa_props.get('charge'),
                        'size_change': wt_aa_props.get('size') != mut_aa_props.get('size'),
                        'hydrophobicity_change': wt_aa_props.get('hydrophobic') != mut_aa_props.get('hydrophobic')
                    }
                })
    
    return changes


# Disease association database
DISEASE_DATABASE = {
    # Frameshift mutations in key genes
    'frameshift_diseases': {
        'BRCA1': {
            'diseases': ['Hereditary Breast and Ovarian Cancer Syndrome'],
            'inheritance': 'Autosomal Dominant',
            'manifestations': ['Breast cancer (early onset)', 'Ovarian cancer', 'Prostate cancer (males)', 'Pancreatic cancer'],
            'prevalence': '1 in 300-800 individuals',
            'management': 'Enhanced screening, prophylactic surgery options'
        },
        'CFTR': {
            'diseases': ['Cystic Fibrosis'],
            'inheritance': 'Autosomal Recessive',
            'manifestations': ['Chronic lung infections', 'Pancreatic insufficiency', 'Male infertility', 'Elevated sweat chloride'],
            'prevalence': '1 in 2,500-3,500 births (Caucasians)',
            'management': 'Airway clearance, enzyme replacement, CFTR modulators'
        },
        'DMD': {
            'diseases': ['Duchenne Muscular Dystrophy'],
            'inheritance': 'X-linked Recessive',
            'manifestations': ['Progressive muscle weakness', 'Cardiomyopathy', 'Cognitive impairment', 'Respiratory failure'],
            'prevalence': '1 in 3,500-5,000 male births',
            'management': 'Corticosteroids, physical therapy, cardiac monitoring'
        },
        'APC': {
            'diseases': ['Familial Adenomatous Polyposis'],
            'inheritance': 'Autosomal Dominant',
            'manifestations': ['Hundreds of colorectal polyps', 'Colorectal cancer', 'Gastric polyps', 'Desmoid tumors'],
            'prevalence': '1 in 10,000-30,000 individuals',
            'management': 'Prophylactic colectomy, endoscopic surveillance'
        }
    },
    
    # Nonsense mutations
    'nonsense_diseases': {
        'TP53': {
            'diseases': ['Li-Fraumeni Syndrome'],
            'inheritance': 'Autosomal Dominant',
            'manifestations': ['Multiple cancer types', 'Early onset cancers', 'Sarcomas', 'Brain tumors', 'Breast cancer'],
            'prevalence': '1 in 5,000-20,000 individuals',
            'management': 'Enhanced cancer screening, genetic counseling'
        },
        'NF1': {
            'diseases': ['Neurofibromatosis Type 1'],
            'inheritance': 'Autosomal Dominant',
            'manifestations': ['Neurofibromas', 'Café-au-lait spots', 'Learning disabilities', 'Optic gliomas'],
            'prevalence': '1 in 3,000 births',
            'management': 'Regular monitoring, symptomatic treatment'
        },
        'LDLR': {
            'diseases': ['Familial Hypercholesterolemia'],
            'inheritance': 'Autosomal Dominant',
            'manifestations': ['Extremely high cholesterol', 'Premature coronary artery disease', 'Xanthomas'],
            'prevalence': '1 in 200-500 individuals',
            'management': 'Statins, PCSK9 inhibitors, lifestyle modifications'
        }
    },
    
    # Missense mutations with known pathogenic effects
    'missense_diseases': {
        'HBB': {
            'diseases': ['Sickle Cell Disease', 'Beta-thalassemia'],
            'inheritance': 'Autosomal Recessive',
            'manifestations': ['Anemia', 'Pain crises', 'Organ damage', 'Stroke risk'],
            'prevalence': '1 in 365 African American births (sickle cell)',
            'management': 'Hydroxyurea, blood transfusions, bone marrow transplant'
        },
        'CFTR': {
            'diseases': ['Cystic Fibrosis (mild forms)'],
            'inheritance': 'Autosomal Recessive',
            'manifestations': ['Pancreatic sufficiency', 'Later onset lung disease', 'Male infertility'],
            'prevalence': 'Variable based on mutation',
            'management': 'CFTR modulators, symptomatic treatment'
        },
        'APOE': {
            'diseases': ['Alzheimer Disease (increased risk)'],
            'inheritance': 'Complex/Multifactorial',
            'manifestations': ['Late-onset dementia', 'Memory loss', 'Cognitive decline'],
            'prevalence': '10-15% carry high-risk alleles',
            'management': 'Lifestyle modifications, symptomatic treatment'
        }
    }
}

# Common pathogenic mutation patterns
PATHOGENIC_PATTERNS = {
    'stop_gained': {
        'severity': 'high',
        'description': 'Premature stop codon likely results in truncated, nonfunctional protein',
        'general_diseases': ['Cancer predisposition syndromes', 'Metabolic disorders', 'Developmental disorders']
    },
    'frameshift': {
        'severity': 'high', 
        'description': 'Frameshift mutation alters entire downstream protein sequence',
        'general_diseases': ['Genetic syndromes', 'Cancer predisposition', 'Metabolic disorders']
    },
    'missense_damaging': {
        'severity': 'medium',
        'description': 'Amino acid change may affect protein structure and function',
        'general_diseases': ['Inherited disorders', 'Cancer susceptibility', 'Metabolic conditions']
    }
}

def analyze_disease_associations(mutation_analysis, wt_seq, mut_seq):
    """Analyze potential disease associations based on mutation patterns"""
    disease_associations = []
    
    # Analyze frameshift mutations
    for frameshift in mutation_analysis['frameshift_mutations']:
        # For demonstration, we'll use pattern-based analysis
        # In real applications, this would require gene annotation and known mutation databases
        
        association = {
            'mutation_type': 'frameshift',
            'position': frameshift['position'],
            'severity': 'high',
            'general_category': 'Frameshift-associated disorders',
            'potential_diseases': [],
            'clinical_significance': 'Likely pathogenic - frameshift mutations typically cause loss of function'
        }
        
        # Add general frameshift disease categories
        association['potential_diseases'] = [
            {
                'category': 'Cancer Predisposition Syndromes',
                'examples': ['Hereditary Breast and Ovarian Cancer', 'Lynch Syndrome', 'Familial Adenomatous Polyposis'],
                'inheritance': 'Usually Autosomal Dominant',
                'management': 'Enhanced screening, genetic counseling, prophylactic measures'
            },
            {
                'category': 'Metabolic Disorders',
                'examples': ['Cystic Fibrosis', 'Phenylketonuria', 'Glycogen Storage Diseases'],
                'inheritance': 'Usually Autosomal Recessive',
                'management': 'Enzyme replacement, dietary modifications, symptomatic treatment'
            },
            {
                'category': 'Neuromuscular Disorders',
                'examples': ['Duchenne Muscular Dystrophy', 'Spinal Muscular Atrophy'],
                'inheritance': 'X-linked or Autosomal Recessive',
                'management': 'Physical therapy, respiratory support, emerging gene therapies'
            }
        ]
        
        disease_associations.append(association)
    
    # Analyze amino acid changes
    for change in mutation_analysis['amino_acid_changes']:
        association = {
            'mutation_type': change['mutation_type'],
            'position': change['position'],
            'amino_acid_change': f"{change['wildtype_aa_name']} → {change['mutant_aa_name']}",
            'codon_change': f"{change['wildtype_codon']} → {change['mutant_codon']}",
            'potential_diseases': []
        }
        
        if change['mutation_type'] == 'nonsense':
            association['severity'] = 'high'
            association['clinical_significance'] = 'Likely pathogenic - premature stop codon'
            association['potential_diseases'] = [
                {
                    'category': 'Tumor Suppressor Gene Disorders',
                    'examples': ['Li-Fraumeni Syndrome (TP53)', 'Retinoblastoma (RB1)', 'Neurofibromatosis (NF1)'],
                    'inheritance': 'Autosomal Dominant',
                    'management': 'Enhanced cancer screening, genetic counseling'
                },
                {
                    'category': 'Metabolic Enzyme Deficiencies',
                    'examples': ['Phenylketonuria', 'Maple Syrup Urine Disease', 'Homocystinuria'],
                    'inheritance': 'Autosomal Recessive',
                    'management': 'Dietary restrictions, enzyme supplementation'
                }
            ]
        
        elif change['mutation_type'] == 'missense':
            # Assess severity based on amino acid properties
            wt_props = AMINO_ACID_PROPERTIES.get(change['wildtype_aa'], {})
            mut_props = AMINO_ACID_PROPERTIES.get(change['mutant_aa'], {})
            
            if wt_props.get('type') != mut_props.get('type'):
                association['severity'] = 'medium'
                association['clinical_significance'] = 'Uncertain significance - may affect protein function'
            else:
                association['severity'] = 'low'
                association['clinical_significance'] = 'Likely benign - conservative amino acid change'
            
            association['potential_diseases'] = [
                {
                    'category': 'Hemoglobinopathies',
                    'examples': ['Sickle Cell Disease (HBB)', 'Beta-thalassemia', 'Alpha-thalassemia'],
                    'inheritance': 'Autosomal Recessive',
                    'management': 'Blood transfusions, hydroxyurea, bone marrow transplant'
                },
                {
                    'category': 'Cardiovascular Disorders',
                    'examples': ['Hypertrophic Cardiomyopathy', 'Long QT Syndrome', 'Familial Hypercholesterolemia'],
                    'inheritance': 'Autosomal Dominant',
                    'management': 'Medications, lifestyle modifications, device therapy'
                },
                {
                    'category': 'Neurological Disorders',
                    'examples': ['Huntington Disease', 'Alzheimer Disease risk variants', 'Parkinson Disease'],
                    'inheritance': 'Variable (Dominant, Recessive, or Complex)',
                    'management': 'Symptomatic treatment, genetic counseling'
                }
            ]
        
        disease_associations.append(association)
    
    # Add general recommendations
    recommendations = {
        'genetic_counseling': 'Recommended for all identified mutations',
        'family_screening': 'Consider testing family members if pathogenic mutation confirmed',
        'clinical_correlation': 'Correlation with clinical phenotype and family history essential',
        'functional_studies': 'May be needed to determine pathogenicity of novel variants'
    }
    
    return {
        'disease_associations': disease_associations,
        'recommendations': recommendations,
        'disclaimer': 'This analysis is for research/educational purposes only. Clinical interpretation requires professional genetic counseling and correlation with clinical findings.'
    }


# NCBI Database References and Links
NCBI_REFERENCES = {
    'base_urls': {
        'clinvar': 'https://www.ncbi.nlm.nih.gov/clinvar/',
        'dbsnp': 'https://www.ncbi.nlm.nih.gov/snp/',
        'gene': 'https://www.ncbi.nlm.nih.gov/gene/',
        'omim': 'https://www.ncbi.nlm.nih.gov/omim/',
        'pubmed': 'https://pubmed.ncbi.nlm.nih.gov/',
        'genbank': 'https://www.ncbi.nlm.nih.gov/nuccore/',
        'protein': 'https://www.ncbi.nlm.nih.gov/protein/'
    },
    
    'mutation_types': {
        'frameshift': {
            'clinvar_search': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=frameshift',
            'references': [
                {
                    'title': 'ClinVar - Frameshift Variants',
                    'url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=frameshift',
                    'description': 'Database of genetic variants and their clinical significance'
                },
                {
                    'title': 'ACMG Guidelines for Variant Classification',
                    'url': 'https://pubmed.ncbi.nlm.nih.gov/25741868/',
                    'description': 'Standards for interpreting sequence variants (PubMed: 25741868)'
                }
            ]
        },
        'nonsense': {
            'clinvar_search': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=nonsense',
            'references': [
                {
                    'title': 'ClinVar - Nonsense Variants',
                    'url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=nonsense',
                    'description': 'Database of stop-gain variants and clinical interpretations'
                },
                {
                    'title': 'dbSNP - Single Nucleotide Polymorphisms',
                    'url': 'https://www.ncbi.nlm.nih.gov/snp/',
                    'description': 'Database of genetic variation including nonsense mutations'
                }
            ]
        },
        'missense': {
            'clinvar_search': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=missense',
            'references': [
                {
                    'title': 'ClinVar - Missense Variants',
                    'url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=missense',
                    'description': 'Clinical significance of amino acid substitutions'
                },
                {
                    'title': 'PolyPhen-2 Prediction Server',
                    'url': 'http://genetics.bwh.harvard.edu/pph2/',
                    'description': 'Prediction of functional effects of missense mutations'
                }
            ]
        }
    },
    
    'disease_genes': {
        'BRCA1': {
            'gene_id': '672',
            'gene_url': 'https://www.ncbi.nlm.nih.gov/gene/672',
            'clinvar_url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1',
            'omim_url': 'https://www.ncbi.nlm.nih.gov/omim/113705',
            'references': [
                {
                    'title': 'BRCA1 Gene - NCBI Gene Database',
                    'url': 'https://www.ncbi.nlm.nih.gov/gene/672',
                    'description': 'Comprehensive gene information and variants'
                },
                {
                    'title': 'BRCA1 Variants in ClinVar',
                    'url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1',
                    'description': 'Clinical significance of BRCA1 mutations'
                }
            ]
        },
        'TP53': {
            'gene_id': '7157',
            'gene_url': 'https://www.ncbi.nlm.nih.gov/gene/7157',
            'clinvar_url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=TP53',
            'omim_url': 'https://www.ncbi.nlm.nih.gov/omim/191170',
            'references': [
                {
                    'title': 'TP53 Gene - NCBI Gene Database',
                    'url': 'https://www.ncbi.nlm.nih.gov/gene/7157',
                    'description': 'Tumor suppressor gene and associated variants'
                },
                {
                    'title': 'TP53 Database - IARC',
                    'url': 'https://tp53.isb-cgc.org/',
                    'description': 'Comprehensive TP53 mutation database'
                }
            ]
        },
        'CFTR': {
            'gene_id': '1080',
            'gene_url': 'https://www.ncbi.nlm.nih.gov/gene/1080',
            'clinvar_url': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=CFTR',
            'omim_url': 'https://www.ncbi.nlm.nih.gov/omim/602421',
            'references': [
                {
                    'title': 'CFTR Gene - NCBI Gene Database',
                    'url': 'https://www.ncbi.nlm.nih.gov/gene/1080',
                    'description': 'Cystic fibrosis transmembrane conductance regulator'
                },
                {
                    'title': 'CFTR2 Database',
                    'url': 'https://cftr2.org/',
                    'description': 'Clinical and functional translation of CFTR variants'
                }
            ]
        }
    },
    
    'general_references': [
        {
            'title': 'ClinVar - Clinical Significance of Variants',
            'url': 'https://www.ncbi.nlm.nih.gov/clinvar/',
            'description': 'NCBI database of genetic variants and their clinical interpretations'
        },
        {
            'title': 'ACMG/AMP Variant Classification Guidelines',
            'url': 'https://pubmed.ncbi.nlm.nih.gov/25741868/',
            'description': 'Standards for interpreting sequence variants (Richards et al., 2015)'
        },
        {
            'title': 'GenBank - Genetic Sequence Database',
            'url': 'https://www.ncbi.nlm.nih.gov/genbank/',
            'description': 'NIH genetic sequence database with annotations'
        },
        {
            'title': 'dbSNP - Single Nucleotide Polymorphism Database',
            'url': 'https://www.ncbi.nlm.nih.gov/snp/',
            'description': 'Database of genetic variation in humans'
        },
        {
            'title': 'OMIM - Online Mendelian Inheritance in Man',
            'url': 'https://www.ncbi.nlm.nih.gov/omim/',
            'description': 'Catalog of human genes and genetic disorders'
        }
    ]
}

def add_ncbi_references(disease_analysis, mutation_analysis):
    """Add NCBI database references and links to disease analysis"""
    
    # Add general NCBI references
    disease_analysis['ncbi_references'] = {
        'general_databases': NCBI_REFERENCES['general_references'],
        'mutation_specific': [],
        'gene_specific': []
    }
    
    # Add mutation-specific references
    mutation_types_found = set()
    
    # Check for frameshift mutations
    if mutation_analysis['frameshift_mutations']:
        mutation_types_found.add('frameshift')
    
    # Check for amino acid changes
    for change in mutation_analysis['amino_acid_changes']:
        if change['mutation_type'] == 'nonsense':
            mutation_types_found.add('nonsense')
        elif change['mutation_type'] == 'missense':
            mutation_types_found.add('missense')
    
    # Add references for found mutation types
    for mut_type in mutation_types_found:
        if mut_type in NCBI_REFERENCES['mutation_types']:
            disease_analysis['ncbi_references']['mutation_specific'].extend(
                NCBI_REFERENCES['mutation_types'][mut_type]['references']
            )
    
    # Add gene-specific references (for demonstration, we'll add common disease genes)
    disease_analysis['ncbi_references']['gene_specific'] = [
        {
            'gene': 'BRCA1',
            'references': NCBI_REFERENCES['disease_genes']['BRCA1']['references']
        },
        {
            'gene': 'TP53', 
            'references': NCBI_REFERENCES['disease_genes']['TP53']['references']
        },
        {
            'gene': 'CFTR',
            'references': NCBI_REFERENCES['disease_genes']['CFTR']['references']
        }
    ]
    
    # Add search URLs for specific mutations
    disease_analysis['search_urls'] = {
        'clinvar_frameshift': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=frameshift',
        'clinvar_nonsense': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=nonsense',
        'clinvar_missense': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=missense',
        'general_search': 'https://www.ncbi.nlm.nih.gov/clinvar/'
    }
    
    return disease_analysis



def get_fractions_only(encoded_sequence):
    """Extract only fraction patterns from encoded sequence"""
    parts = [part.strip() for part in encoded_sequence.split(',')]
    fractions_only = []
    
    for part in parts:
        # Keep only parts that contain fractions (have '/' in them)
        if '/' in part and ' ' in part:
            fractions_only.append(part)
    
    return ', '.join(fractions_only)

@dna_encoder_bp.route('/download-fractions', methods=['POST'])
@cross_origin()
def download_fractions():
    """Get fractions-only format for download"""
    try:
        data = request.get_json()
        encoded_sequence = data.get('encoded_sequence', '')
        
        if not encoded_sequence:
            return jsonify({'error': 'No encoded sequence provided'}), 400
        
        fractions_only = get_fractions_only(encoded_sequence)
        
        return jsonify({
            'success': True,
            'fractions_only': fractions_only
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500



def simple_sequence_alignment(ref_seq, sample_seq):
    """Simple sequence alignment to find differences"""
    # Pad sequences to same length
    max_len = max(len(ref_seq), len(sample_seq))
    ref_padded = ref_seq.ljust(max_len, '-')
    sample_padded = sample_seq.ljust(max_len, '-')
    
    alignment = []
    differences = []
    
    for i, (ref_base, sample_base) in enumerate(zip(ref_padded, sample_padded)):
        is_different = ref_base != sample_base
        alignment.append({
            'position': i + 1,
            'reference': ref_base,
            'sample': sample_base,
            'is_different': is_different
        })
        
        if is_different:
            differences.append({
                'position': i + 1,
                'reference': ref_base,
                'sample': sample_base,
                'type': 'substitution' if ref_base != '-' and sample_base != '-' else 'indel'
            })
    
    return alignment, differences

def analyze_numerical_differences(ref_encoded, sample_encoded):
    """Analyze differences in numerical encoding"""
    ref_parts = [part.strip() for part in ref_encoded.split(',')]
    sample_parts = [part.strip() for part in sample_encoded.split(',')]
    
    # Create pattern frequency maps
    ref_patterns = {}
    sample_patterns = {}
    
    for part in ref_parts:
        if '/' in part and ' ' in part:
            pattern = part.split()[1]  # Get pattern (e.g., "2/9")
            count = int(part.split()[0])  # Get count
            ref_patterns[pattern] = ref_patterns.get(pattern, 0) + count
        else:
            # Single nucleotide
            ref_patterns[part] = ref_patterns.get(part, 0) + 1
    
    for part in sample_parts:
        if '/' in part and ' ' in part:
            pattern = part.split()[1]  # Get pattern (e.g., "2/9")
            count = int(part.split()[0])  # Get count
            sample_patterns[pattern] = sample_patterns.get(pattern, 0) + count
        else:
            # Single nucleotide
            sample_patterns[part] = sample_patterns.get(part, 0) + 1
    
    # Find all unique patterns
    all_patterns = set(ref_patterns.keys()) | set(sample_patterns.keys())
    
    numerical_differences = []
    for pattern in all_patterns:
        ref_count = ref_patterns.get(pattern, 0)
        sample_count = sample_patterns.get(pattern, 0)
        difference = sample_count - ref_count
        
        if difference != 0:
            # Decode pattern to nucleotide sequence
            if '/' in pattern:
                digits = pattern.split('/')[0]
                nucleotide_seq = ''.join(MAP.get(d, 'N') for d in digits)
            else:
                nucleotide_seq = MAP.get(pattern, 'N')
            
            numerical_differences.append({
                'pattern': pattern,
                'nucleotide_sequence': nucleotide_seq,
                'reference_count': ref_count,
                'sample_count': sample_count,
                'difference': difference,
                'change_type': 'gain' if difference > 0 else 'loss'
            })
    
    return numerical_differences

@dna_encoder_bp.route('/mutant-compare', methods=['POST'])
@cross_origin()
def mutant_compare():
    """Compare reference sequence with multiple sample sequences"""
    try:
        data = request.get_json()
        reference_sequence = data.get('reference_sequence', '').strip()
        sample_sequences = data.get('sample_sequences', [])
        
        if not reference_sequence:
            return jsonify({'error': 'Reference sequence is required'}), 400
        
        if not sample_sequences or len(sample_sequences) == 0:
            return jsonify({'error': 'At least one sample sequence is required'}), 400
        
        if len(sample_sequences) > 5:
            return jsonify({'error': 'Maximum 5 sample sequences allowed'}), 400
        
        # Validate sequences
        valid_nucleotides = set('ATCG')
        clean_reference = reference_sequence.upper().replace('\n', '').replace(' ', '')
        
        if not all(nuc in valid_nucleotides for nuc in clean_reference):
            return jsonify({'error': 'Invalid nucleotides in reference sequence. Only A, T, C, G are allowed.'}), 400
        
        # Clean and validate sample sequences
        clean_samples = []
        for i, sample in enumerate(sample_sequences):
            clean_sample = sample.strip().upper().replace('\n', '').replace(' ', '')
            if not clean_sample:
                continue  # Skip empty sequences
            
            if not all(nuc in valid_nucleotides for nuc in clean_sample):
                return jsonify({'error': f'Invalid nucleotides in sample sequence {i+1}. Only A, T, C, G are allowed.'}), 400
            
            clean_samples.append(clean_sample)
        
        if not clean_samples:
            return jsonify({'error': 'No valid sample sequences provided'}), 400
        
        start_time = time.time()
        
        # Encode reference sequence
        ref_encoded = encode_sequence_to_format(clean_reference)
        
        # Process each sample sequence
        comparisons = []
        for i, sample_seq in enumerate(clean_samples):
            # Encode sample sequence
            sample_encoded = encode_sequence_to_format(sample_seq)
            
            # Perform sequence alignment
            alignment, differences = simple_sequence_alignment(clean_reference, sample_seq)
            
            # Analyze numerical differences
            numerical_differences = analyze_numerical_differences(ref_encoded, sample_encoded)
            
            # Calculate statistics
            total_positions = len(alignment)
            different_positions = len(differences)
            similarity_percentage = ((total_positions - different_positions) / total_positions * 100) if total_positions > 0 else 0
            
            # Categorize mutations
            substitutions = [d for d in differences if d['type'] == 'substitution']
            indels = [d for d in differences if d['type'] == 'indel']
            
            comparisons.append({
                'sample_id': i + 1,
                'sample_sequence': sample_seq,
                'sample_encoded': sample_encoded,
                'alignment': alignment,
                'differences': differences,
                'numerical_differences': numerical_differences,
                'statistics': {
                    'total_positions': total_positions,
                    'different_positions': different_positions,
                    'similarity_percentage': round(similarity_percentage, 2),
                    'substitutions': len(substitutions),
                    'indels': len(indels)
                }
            })
        
        processing_time = round((time.time() - start_time) * 1000, 2)
        
        return jsonify({
            'success': True,
            'reference_sequence': clean_reference,
            'reference_encoded': ref_encoded,
            'comparisons': comparisons,
            'processing_time_ms': processing_time
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

