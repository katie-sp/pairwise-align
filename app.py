from flask import Flask, request, jsonify
from flask_cors import CORS
import align

app = Flask(__name__)
CORS(app)

@app.route('/')
def index():
    return 'Web App with Python Flask!'

@app.route('/align', methods=['POST'])
def align_seqs():
    data = request.get_json()  # dictionary
    seq1 = data['seq1']
    seq2 = data['seq2']
    seqtype = data['seqtype']
    seq1name = data['seq1name']
    seq2name = data['seq2name']
    output = align.align_seqs(seq1, seq2, seqtype, seq1name, seq2name)
    if type(output) == str: # encountered an error
        return jsonify({'error':output})
    else:
        return jsonify({'alignment-top':output['alignment-top'],
                        'alignment-symbols':output['alignment-symbols'],
                        'alignment-bottom':output['alignment-bottom'],
                        'top-index':output['top-index'],
                        'bottom-index':output['bottom-index'],
                        'mismatch':output['mismatch'],
                        'gap':output['gap'],
                        'clustal':output['clustal'],
                        'spacer':output['spacer']})
        
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=9000, debug=True)
