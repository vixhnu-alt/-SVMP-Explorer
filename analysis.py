"""Interpretable, literature-assisted SVMP analysis; not a validated diagnostic tool."""
from pathlib import Path
import json, re, math
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parent
DB = ROOT/'database'
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')
VALID_NT = set('ACGTU')


def load_motifs():return json.loads((DB/'motifs.json').read_text())
def load_references():
    p=DB/'reference_records.json'
    return json.loads(p.read_text()) if p.exists() else []

def parse_input(text, kind='auto', frame=0):
    """Only complete CDS is safely translated; frame is 0,1,2 for positive strand."""
    lines=[s.strip() for s in text.splitlines() if s.strip() and not s.startswith('>')]
    seq=''.join(lines).upper().replace(' ','')
    if not seq:raise ValueError('Please enter a FASTA or raw sequence.')
    if kind=='auto':kind='mRNA / CDS' if set(seq)<=VALID_NT else 'Protein'
    if kind=='mRNA / CDS':
        if not set(seq)<=VALID_NT:raise ValueError('Nucleotide sequence contains invalid characters.')
        nt=seq.replace('U','T')[frame:]
        if len(nt)%3:raise ValueError('Coding sequence length after selected frame must be divisible by 3.')
        translated=str(Seq(nt).translate(table=1))
        if '*' in translated[:-1]:raise ValueError('Internal stop codon found. Check reading frame or CDS completeness.')
        if translated.endswith('*'):translated=translated[:-1]
        if not translated:raise ValueError('No translated amino acids.')
        return translated,kind
    if kind!='Protein':raise ValueError('Unsupported input type')
    if seq.endswith('*'):seq=seq[:-1]
    if not set(seq)<=VALID_AA:raise ValueError('Protein must use 20 standard amino acids only (no X or internal *).')
    return seq,kind

def scan_motifs(seq):
    matches=[]
    for category,entry in load_motifs().items():
        for pattern in entry['patterns']:
            for m in re.finditer(pattern if entry['kind']=='regex' else re.escape(pattern),seq):
                matches.append({'category':category,'domain':entry['domain'],'pattern':pattern,
                                'match':m.group(),'start':m.start()+1,'end':m.end(),
                                'references':','.join(entry['references'])})
    # De-duplicate overlapping short and long motif descriptions only if identical start/end/name.
    return sorted({(m['category'],m['start'],m['end'],m['match']):m for m in matches}.values(),
                  key=lambda x:(x['start'],x['end'],x['category']))

def properties(seq):
    p=ProteinAnalysis(seq)
    return {'Length':len(seq),'Molecular weight (Da)':round(p.molecular_weight(),1),
            'Theoretical pI':round(p.isoelectric_point(),2),
            'GRAVY':round(p.gravy(),3),'Instability index':round(p.instability_index(),2),
            'Aromaticity':round(p.aromaticity(),3),'Cysteine %':round(seq.count('C')*100/len(seq),2)}

def align_one(query,target):
    """Local alignment with BLOSUM62 and affine gaps; report coverage and identity."""
    from Bio.Align import substitution_matrices
    aligner=PairwiseAligner(mode='local')
    aligner.substitution_matrix=substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score=-10
    aligner.extend_gap_score=-0.5
    aln=aligner.align(query,target)[0]
    qblocks,tblocks=aln.aligned
    aligned_pairs=sum(b-a for a,b in qblocks)
    identical=0
    for (qa,qb),(ta,tb) in zip(qblocks,tblocks):
        identical+=sum(a==b for a,b in zip(query[qa:qb],target[ta:tb]))
    return {'identity':float(identical/aligned_pairs) if aligned_pairs else 0,
            'query_coverage':float(aligned_pairs/len(query)),
            'reference_coverage':float(aligned_pairs/len(target)),
            'aligned_residues':int(aligned_pairs),
            'query_start':int(qblocks[0][0])+1 if len(qblocks) else None,
            'query_end':int(qblocks[-1][1]) if len(qblocks) else None,
            'reference_start':int(tblocks[0][0])+1 if len(tblocks) else None,
            'reference_end':int(tblocks[-1][1]) if len(tblocks) else None}

def domain_fragments(ref):
    """Only uses UniProt domain annotations; doesn't invent C-domain boundaries."""
    seq=ref['sequence'];out=[]
    for f in ref.get('features',[]):
        if f['type']!='Domain':continue
        name=f['description'].lower()
        domain='M' if ('peptidase m12' in name or 'metalloprote' in name) else ('D_like' if 'disintegrin-like' in name else ('D' if 'disintegrin' in name else ('C' if 'cysteine-rich' in name else None)))
        if domain and f['end']<=len(seq) and f['start']>=1:
            # UniProt often labels both P-II D and P-III D-like as 'Disintegrin'.
            if domain=='D' and ref['class']=='P-III':domain='D_like'
            out.append({'domain':domain,'sequence':seq[f['start']-1:f['end']],
                        'ref_start':f['start'],'ref_end':f['end'],'accession':ref['accession']})
    # In P-III UniProt may annotate M and D but omit the C domain.
    # Never pretend a guessed C region is a curated annotation.
    return out

def compare_domains(seq,refs):
    results=[]
    for ref in refs:
        for frag in domain_fragments(ref):
            if len(frag['sequence'])<35:continue
            a=align_one(seq,frag['sequence'])
            results.append({**a,'domain':frag['domain'],'accession':ref['accession'],
                            'reference_class':ref['class'],'domain_length':len(frag['sequence']),
                            'evidence':'UniProt annotated domain + local alignment'})
    return sorted(results,key=lambda a:(a['reference_coverage']*a['identity'],a['aligned_residues']),reverse=True)

def classify(seq,refs,motifs):
    if not refs:return {'prediction':'Not assessed','confidence':'Unavailable','reason':'Reference database not downloaded. Run python build_database.py first.','domains':[],'top_references':[]}
    scores=[]
    for r in refs:
        try:a=align_one(seq,r['sequence'])
        except Exception:continue
        # Require both query and reference coverage to avoid short local matches being decisive.
        metric=a['identity']*min(a['query_coverage'],a['reference_coverage'])
        scores.append({**a,'accession':r['accession'],'class':r['class'],'species':r.get('species',''), 'score':metric})
    scores.sort(key=lambda a:a['score'],reverse=True)
    domains=compare_domains(seq,refs)
    # Evidence from longer domains, not raw motif hits.
    strong={}
    for d in domains:
        if d['identity']>=.38 and d['reference_coverage']>=.65 and d['aligned_residues']>=45:
            strong.setdefault(d['domain'],d)
    zinc=any(m['category']=='zinc_binding' for m in motifs)
    predicted='Uncertain';reason='Insufficient domain evidence to assign a class.'
    if 'M' in strong:
        if 'D_like' in strong and 'C' in strong:
            predicted='P-III';reason='M, D-like and C domains supported by reference alignments.'
        elif 'D' in strong and 'D_like' not in strong:
            predicted='P-II';reason='M and canonical D domains supported by reference alignments.'
        elif 'D_like' in strong:
            predicted='Uncertain (possible P-III)';reason='M and D-like domains supported, but C domain is not independently confirmed.'
        elif not any(k in strong for k in ('D','D_like','C')):
            # Absence of an annotated D match is not proof of its absence.
            predicted='Uncertain (possible P-I)';reason='M domain found; no independently confirmed accessory domains.'
    # Full-length reference similarity is secondary and cannot establish domain architecture.
    if scores and scores[0]['score']>=.60:
        top=scores[0]
        if predicted.startswith('Uncertain'):
            predicted=f"{top['class']} (reference-assisted, provisional)"
            reason+=f" Strong full-length similarity to {top['accession']} ({top['class']}); not a domain-confirmed classification."
        elif predicted!=top['class']:
            predicted='Uncertain';reason+=' Reference and domain evidence disagree.'
    if not zinc:reason+=' Canonical catalytic zinc motif not detected (may be inactive or truncated).'
    if not any(d['domain']=='C' for d in domains):reason+=' Curated C-domain boundaries were unavailable in downloaded reference annotations.'
    return {'prediction':predicted,'confidence':'Provisional (research prototype)',
            'reason':reason,'domains':domains,'supported_domains':strong,'top_references':scores[:8]}
