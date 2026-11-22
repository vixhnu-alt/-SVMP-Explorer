import json
from pathlib import Path
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


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
        if f['type'] not in ('Domain','Region'):continue
        name=f['description'].lower()
        domain='M' if ('peptidase m12' in name or 'metalloprote' in name) else ('D_like' if 'disintegrin-like' in name else ('D' if 'disintegrin' in name else ('C' if ('cysteine-rich' in name or 'cysteine rich' in name) else None)))
        if f['type']=='Region' and domain!='C':continue
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

def _supported_domain_candidates(domains):
    """Require substantial aligned reference-domain coverage and reliable annotations.

    These are provisional homology calls, not proof of domain boundaries.
    """
    out = {}
    for d in domains:
        if (d['identity'] < 0.42 or d['reference_coverage'] < 0.75
                or d['aligned_residues'] < 55 or d['query_start'] is None):
            continue
        out.setdefault(d['domain'], []).append(d)
    return out


def _ordered_nonoverlapping(hits, names):
    """Find distinct ordered regions, allowing at most 15 residues overlap."""
    def walk(i, chosen):
        if i == len(names):
            return chosen
        for hit in hits.get(names[i], []):
            if chosen and hit['query_start'] < chosen[-1]['query_end'] - 15:
                continue
            found = walk(i + 1, chosen + [hit])
            if found:
                return found
        return None
    ordered = {key: sorted(vals, key=lambda d: (-d['identity']*d['reference_coverage'], d['query_start']))
               for key, vals in hits.items()}
    result = walk(0, [])
    return dict(zip(names, result)) if result else None


def classify(seq,refs,motifs):
    # Domain architecture, not short motifs or cysteine percentage, determines class.

    if not refs:
        return {'prediction':'Not assessed','confidence':'Unavailable',
                'reason':'Reference database not downloaded. Run python build_database.py first.',
                'domains':[],'supported_domains':{},'top_references':[]}
    # Exact accession-independent sequence match to an actual downloaded record.
    # Do not infer class from the candidate label alone: require annotated architecture.
    for ref in refs:
        if seq != ref.get('sequence'):continue
        fragments=domain_fragments(ref)
        kinds={f['domain'] for f in fragments}
        cl=ref['class']
        valid=(cl=='P-II' and {'M','D'}<=kinds and 'C' not in kinds) or (cl=='P-III' and {'M','D_like','C'}<=kinds) or (cl=='P-I' and 'M' in kinds and not ({'D','D_like','C'} & kinds))
        if valid:
            return {'prediction':cl,'confidence':'Exact sequence match + annotated architecture',
                    'reason':f"Exact sequence match to downloaded UniProt entry {ref['accession']} with supporting annotated domain architecture. This is a reference lookup, not de novo domain validation.",
                    'domains':[],'supported_domains':{},'top_references':[{'accession':ref['accession'],'class':cl,'identity':1.0,'score':1.0}]}
    scores=[]
    for r in refs:
        try:
            a=align_one(seq,r['sequence'])
        except Exception:
            continue
        metric=a['identity']*min(a['query_coverage'],a['reference_coverage'])
        scores.append({**a,'accession':r['accession'],'class':r['class'],
                       'species':r.get('species',''),'score':metric})
    scores.sort(key=lambda a:a['score'],reverse=True)
    domains=compare_domains(seq,refs)
    candidates=_supported_domain_candidates(domains)
    p3=_ordered_nonoverlapping(candidates,['M','D_like','C'])
    p2=_ordered_nonoverlapping(candidates,['M','D'])
    m=_ordered_nonoverlapping(candidates,['M'])
    # Prevent false P-III calls when a complete, annotated P-II protein is
    # strongly homologous to the query. A C hit must be separate, sufficiently
    # long and supported by an annotated P-III C-domain reference.
    strong_p2 = [
        s for s in scores if s['class'] == 'P-II'
        and s['identity'] >= 0.85
        and s['query_coverage'] >= 0.85
        and s['reference_coverage'] >= 0.85
        and any(r['accession'] == s['accession']
                and {'M', 'D'} <= {f['domain'] for f in domain_fragments(r)}
                and 'C' not in {f['domain'] for f in domain_fragments(r)}
                for r in refs)
    ]
    strong_p2.sort(key=lambda s:s['score'], reverse=True)
    if strong_p2 and p2:
        prediction='P-II'
        reason=('Strong near-full-length similarity to annotated P-II reference '
                + strong_p2[0]['accession']
                + ' and separately supported M and D domains. '
                'Incidental cysteine-rich or D-like similarities are not sufficient to override this.')
        supported=p2
    elif p3:
        prediction='P-III'
        reason='Three separate ordered regions align to curated M, D-like and C reference domains.'
        supported=p3
    elif p2:
        # If there is a substantial D-like or C match, an apparent P-II may be incomplete.
        if candidates.get('D_like') or candidates.get('C'):
            prediction='Uncertain'
            reason='M+D alignment found, but accessory-domain evidence conflicts; cannot safely exclude P-III.'
            supported={}
        else:
            prediction='P-II'
            reason='Separate ordered M and canonical D domain alignments are supported; no C domain is confirmed.'
            supported=p2
    elif m:
        prediction='Uncertain (possible P-I)'
        reason='M domain is supported, but lack of an accessory-domain match does not establish its absence.'
        supported=m
    else:
        prediction='Uncertain'
        reason='No complete, correctly ordered domain architecture passed the alignment thresholds.'
        supported={}
    if scores:
        top=scores[0]
        reason+=f" Closest full-length reference: {top['accession']} ({top['class']}); similarity is supplementary only."
        if prediction in ('P-II','P-III') and not strong_p2 and top['class'] != prediction and top['score'] >= 0.60:
            prediction='Uncertain'
            reason+=' Strong reference and domain evidence disagree.'
    if not any(x['category']=='zinc_binding' for x in motifs):
        reason+=' Canonical zinc-binding motif not detected.'
    if not candidates.get('C'):
        reason+=' No independently annotated C-domain alignment passed the threshold; cysteine abundance is not used as C-domain evidence.'
    return {'prediction':prediction,'confidence':'Provisional (research prototype)',
            'reason':reason,'domains':domains,'supported_domains':supported,
            'top_references':scores[:8]}


st.set_page_config(page_title='SVMP Explorer',page_icon='🧬',layout='wide')
st.title('🧬 SVMP Explorer')
st.caption('Literature-curated motif screening and reference-assisted P-I / P-II / P-III classification')
st.warning('Research/education prototype. Short motifs do not confirm domains. Predictions are provisional and require independent annotation.')

refs=load_references()
a,b,c=st.columns(3)
a.metric('Downloaded reference proteins',len(refs))
b.metric('Classes represented',len({r['class'] for r in refs}))
c.metric('Motif patterns',sum(len(x['patterns']) for x in load_motifs().values()))
if not refs:
    st.error('Reference FASTA files are not yet downloaded. Open Command Prompt in this folder and run: python build_database.py. Motif analysis is available meanwhile.')

with st.sidebar:
    st.header('Input settings')
    input_type=st.selectbox('Sequence type',['Auto-detect','Protein','mRNA / CDS'])
    frame=st.selectbox('mRNA reading frame (5′→3′)',['1','2','3'])
    st.info('mRNA input must be a coding sequence in the correct reading frame. Noncoding regions and reverse-complement ORFs are not automatically inferred.')
    st.markdown('**Classification rules**')
    st.write('P-I: M only; P-II: M + D; P-III: M + D-like + C. All domain calls require sequence-level support.')

seq_text=st.text_area('Paste a protein or coding mRNA FASTA sequence',height=180,
                     placeholder='>Unknown_SVMP\nM...')
file=st.file_uploader('Or upload a single-sequence FASTA (.fasta, .fa, .faa, .fna, .txt)',type=['fasta','fa','faa','fna','txt'])
if file is not None:seq_text=file.getvalue().decode('utf-8-sig')


def make_map(length,motifs,domain_hits):
    fig,ax=plt.subplots(figsize=(13,3.4))
    ax.hlines(1,1,length,color='gray',lw=5,alpha=.28)
    colors={'M':'#d05d54','D':'#408fc1','D_like':'#408fc1','C':'#4a9d84'}
    # Plot only independently supported domain-alignment regions.
    lane=1.13
    for i,(name,hit) in enumerate(domain_hits.items()):
        s,e=hit['query_start'],hit['query_end']
        if s is None or e is None:continue
        ax.add_patch(Rectangle((s,lane+i*.23),max(1,e-s),.17,color=colors.get(name,'#888'),alpha=.75))
        ax.text((s+e)/2,lane+i*.23+.09,name,ha='center',va='center',color='white',fontsize=8,weight='bold')
    for i,m in enumerate(motifs):
        y=.9 if i%2==0 else .78
        ax.plot([m['start'],m['end']],[y,y],lw=5,color='#e5a83c',solid_capstyle='round')
        ax.text(m['start'],y-.07,f"{m['match']} ({m['start']})",rotation=35,fontsize=7,ha='left',va='top')
    ax.set_xlim(0,length+5);ax.set_ylim(.2,1.95)
    ax.set_yticks([]);ax.set_xlabel('Amino acid position (1-based)')
    ax.set_title('Sequence map: yellow = motif; colored rectangles = supported aligned regions')
    for spine in ['left','right','top']:ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig

if st.button('Analyze sequence',type='primary',use_container_width=True):
    try:
        typ='auto' if input_type=='Auto-detect' else input_type
        seq,actual_type=parse_input(seq_text,typ,int(frame)-1)
        motifs=scan_motifs(seq)
        result=classify(seq,refs,motifs)
        st.session_state['last']={'sequence':seq,'kind':actual_type,'motifs':motifs,'result':result}
    except Exception as ex:st.error(str(ex))

if 'last' in st.session_state:
    data=st.session_state['last'];seq=data['sequence'];motifs=data['motifs'];result=data['result']
    st.subheader('Classification')
    st.metric('Predicted class',result['prediction'])
    st.write('**Evidence assessment:**',result['reason'])
    t1,t2,t3,t4=st.tabs(['Motifs and domain map','Protein properties','Reference comparison','Database and sources'])
    with t1:
        hits=result.get('supported_domains',{})
        st.pyplot(make_map(len(seq),motifs,hits),clear_figure=True)
        if motifs:
            st.dataframe(pd.DataFrame(motifs),use_container_width=True,hide_index=True)
        else:st.info('No catalogued short motifs detected.')
        if hits:
            st.write('**Supported domain alignments**')
            st.dataframe(pd.DataFrame([{'domain':k,**v} for k,v in hits.items()]),use_container_width=True,hide_index=True)
        else:st.info('No domain regions met the prototype alignment criteria.')
        st.caption('Motif matches may overlap. Domain rectangles represent aligned regions, not experimentally verified boundaries in the input.')
    with t2:
        st.dataframe(pd.DataFrame([properties(seq)]).T.rename(columns={0:'Value'}),use_container_width=True)
        st.write('**Translated/analyzed amino acid sequence**')
        st.code('\n'.join(seq[i:i+70] for i in range(0,len(seq),70)))
    with t3:
        if result['top_references']:
            st.dataframe(pd.DataFrame(result['top_references']),use_container_width=True,hide_index=True)
            st.caption('Similarity score = identity × minimum(query coverage, reference coverage). It is not a statistical probability.')
        else:st.info('Run build_database.py on a computer with internet to obtain real reference sequences.')
    with t4:
        st.write('**Motif dictionary**')
        st.json(load_motifs())
        st.write('**References**')
        sources=json.loads((Path(__file__).resolve().parent/'database'/'sources.json').read_text())
        st.dataframe(pd.DataFrame([{'ID':k,**v} for k,v in sources.items()]),use_container_width=True,hide_index=True)
        if refs:st.dataframe(pd.DataFrame([{'Accession':r['accession'],'Class':r['class'],'Species':r['species'],'Length':r['length'],'Name':r['name']} for r in refs]),hide_index=True)
    output={'input_kind':data['kind'],'protein_sequence':seq,'motif_matches':motifs,
            'prediction':result['prediction'],'interpretation':result['reason'],
            'reference_hits':result['top_references']}
    st.download_button('Download JSON report',data=json.dumps(output,indent=2),file_name='svmp_analysis.json',mime='application/json')
