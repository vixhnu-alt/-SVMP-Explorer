import json
from pathlib import Path
import re
import subprocess
import sys
import pandas as pd
import streamlit as st

st.set_page_config(page_title="SVMP Explorer", page_icon="🧬", layout="wide")

BASE_DIR = Path(__file__).resolve().parent
DB_DIR = BASE_DIR / "database"
DB_DIR.mkdir(parents=True, exist_ok=True)
MOTIFS_PATH = DB_DIR / "motifs.json"
REFS_PATH = DB_DIR / "references.json"

DEFAULT_MOTIFS = {
    "P-I": {
        "description": "Metalloproteinase domain only",
        "patterns": ["HE[LIVMF]GH[LIVMF]LG[LIVMF]XH", "CIM[AS]"]
    },
    "P-II": {
        "description": "Metalloproteinase + Disintegrin domain",
        "patterns": ["HE[LIVMF]GH[LIVMF]LG[LIVMF]XH", "RGD", "KGD", "MVD", "ECD"]
    },
    "P-III": {
        "description": "Metalloproteinase + Disintegrin-like + Cysteine-rich domains",
        "patterns": ["HE[LIVMF]GH[LIVMF]LG[LIVMF]XH", "ECD", "D[DE]CD", "C[KR]C"]
    }
}

DEFAULT_REFERENCES = [
    {"accession": "P0C7B0", "name": "Atrolysin-A", "class": "P-III", "species": "Crotalus atrox"},
    {"accession": "P15167", "name": "Atrolysin-C", "class": "P-I", "species": "Crotalus atrox"},
    {"accession": "P0C6E3", "name": "Bilitoxin-1", "class": "P-II", "species": "Agkistrodon bilineatus"},
    {"accession": "Q90495", "name": "Ecarin", "class": "P-III", "species": "Echis carinatus"},
    {"accession": "P83512", "name": "BaP1", "class": "P-I", "species": "Bothrops asper"},
    {"accession": "Q8QG88", "name": "Bothropasin", "class": "P-III", "species": "Bothrops jararaca"}
]


def ensure_database():
    build_script = BASE_DIR / "build_database.py"
    if not MOTIFS_PATH.exists() and build_script.exists():
        subprocess.run([sys.executable, str(build_script)], cwd=str(BASE_DIR), check=False)

    if not MOTIFS_PATH.exists():
        MOTIFS_PATH.write_text(json.dumps(DEFAULT_MOTIFS, indent=2), encoding="utf-8")

    if not REFS_PATH.exists():
        REFS_PATH.write_text(json.dumps(DEFAULT_REFERENCES, indent=2), encoding="utf-8")


ensure_database()


def load_motifs():
    return json.loads(MOTIFS_PATH.read_text(encoding="utf-8"))


def load_references():
    if REFS_PATH.exists():
        return json.loads(REFS_PATH.read_text(encoding="utf-8"))
    return DEFAULT_REFERENCES


def parse_fasta(raw_text):
    sequences = []
    header = "Query_1"
    chunks = []
    for line in raw_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if chunks:
                sequences.append((header, "".join(chunks).upper()))
                chunks = []
            header = line[1:].strip() or f"Query_{len(sequences) + 1}"
        else:
            chunks.append(re.sub(r"[^A-Za-z]", "", line))
    if chunks:
        sequences.append((header, "".join(chunks).upper()))
    return sequences


def screen_sequence(seq, motifs):
    matches = {}
    matched_patterns = []
    for svmp_class, info in motifs.items():
        class_hits = []
        for pat in info.get("patterns", []):
            regex = pat.replace("X", "[A-Z]")
            for m in re.finditer(regex, seq):
                class_hits.append(f"{pat} ({m.start() + 1}-{m.end()})")
                matched_patterns.append((svmp_class, pat, m.group(0), m.start() + 1, m.end()))
        matches[svmp_class] = class_hits

    if matches.get("P-III") and any("ECD" in h or "CD" in h or "C[" in h for h in matches["P-III"]):
        predicted = "P-III"
    elif matches.get("P-II") and any("RGD" in h or "KGD" in h or "MVD" in h for h in matches["P-II"]):
        predicted = "P-II"
    elif any(matches.values()):
        predicted = "P-I"
    else:
        predicted = "Unclassified / Non-SVMP"

    return predicted, matches, matched_patterns


motifs_data = load_motifs()
refs_data = load_references()
classes_count = len({r.get("class") for r in refs_data if r.get("class")}) or len(motifs_data)

st.title("🧬 SVMP Explorer")
st.caption("Literature-curated motif screening and reference-assisted P-I / P-II / P-III classification")

st.warning("Research/education prototype. Short motifs do not confirm domains. Predictions are provisional and require independent annotation.")

c1, c2, c3 = st.columns(3)
c1.metric("Downloaded reference proteins", len(refs_data))
c2.metric("Classes represented", classes_count)
c3.metric("Motif patterns", sum(len(x.get("patterns", [])) for x in motifs_data.values()))

st.divider()

tab1, tab2, tab3 = st.tabs(["Motif Screening & Classification", "Curated Motif Library", "Reference Proteins"])

with tab1:
    uploaded = st.file_uploader("Upload FASTA file (.fasta, .fa, .txt)", type=["fasta", "fa", "txt"])
    sample_fasta = ">Sample_SVMP_Query\nMIILLVCLVFPYQVSSIILESGNVNDYEVVYPRKVTALPKGAVQPKYEDAMQYEFKVNGEPVVLHLEKNKGLFSEDYSETHYSPDGREITTNPPVEDHCYYHGRIQNDADSTASISACNGLKGHFKLRGEMYFIEPLKIPDSEAHAVYKYENVEKEDEAPKMCGVTQNWKSYEPIKKASQLNLTPEQQRYVELFIVVDHGMFMKYNGNSDKIRRIWQMVNIMKEMYLTYIDIALTGVEIWSNKDLINVQPAAPQTLDSFGEWRERDLLNRISHDNAQLLTAIDLDGPTIGLAYVGTMCDPKRSTGIVQDHSAINLLVAVTMAHELGHNLGMDHDTGSCSCGGYSCIMSPVISDEPSKYFSDCSYIQCWDFIMKENPQCILNEPLGTDIVSPPVCGNELLEVGEECDCGTPENCQNECCDAATCKLKSGSQCGHGDCCEQCKFSKSGTECRASMSECDPAEHCTGQSSECPADVFHKNGQPCLDNYGYCYNGNCPIMYHQCYALFGADVYEAEDSCFERNQKGNYYGYCRKENGNKIPCAPEDVKCGRLYCKDNSPGQNNPCKMFYSNDDEHKGMVLPGTKCADGKVCSNGHCVDVATAY"
    seq_input = st.text_area("Or paste protein sequence(s) in FASTA format", value=sample_fasta, height=180)

    if st.button("Run Screening", type="primary"):
        raw_data = uploaded.getvalue().decode("utf-8", errors="ignore") if uploaded else seq_input
        entries = parse_fasta(raw_data)
        if not entries:
            st.error("Please provide at least one valid protein sequence.")
        else:
            rows = []
            detail_hits = []
            for header, seq in entries:
                pred_class, _, hits = screen_sequence(seq, motifs_data)
                rows.append({
                    "Sequence ID": header,
                    "Length (aa)": len(seq),
                    "Predicted Class": pred_class,
                    "Total Motif Hits": len(hits)
                })
                for cls, pat, matched_seq, start, end in hits:
                    detail_hits.append({
                        "Sequence ID": header,
                        "Associated Class": cls,
                        "Motif Pattern": pat,
                        "Matched Subsequence": matched_seq,
                        "Start": start,
                        "End": end
                    })

            st.subheader("Classification Summary")
            st.dataframe(pd.DataFrame(rows), use_container_width=True)

            st.subheader("Matched Motif Positions")
            if detail_hits:
                st.dataframe(pd.DataFrame(detail_hits), use_container_width=True)
            else:
                st.info("No literature-curated SVMP motifs matched in the submitted sequence(s).")

with tab2:
    motif_rows = []
    for cls, info in motifs_data.items():
        for pat in info.get("patterns", []):
            motif_rows.append({
                "SVMP Class": cls,
                "Domain Architecture": info.get("description", ""),
                "Motif Pattern": pat
            })
    st.dataframe(pd.DataFrame(motif_rows), use_container_width=True)

with tab3:
    st.dataframe(pd.DataFrame(refs_data), use_container_width=True)
