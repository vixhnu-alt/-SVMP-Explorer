import json
from pathlib import Path
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from analysis import parse_input, scan_motifs, properties, load_references, load_motifs, classify

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
