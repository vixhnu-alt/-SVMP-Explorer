# SVMP Explorer — literature-assisted Python project

## What it does
- Protein FASTA or coding mRNA FASTA input (reading frame selectable)
- Translation using Biopython
- Literature-linked motif scanning with 1-based residue positions
- Downloads **real** UniProt reference proteins and curated feature annotations on your own internet-connected PC
- Aligns input against full-length references and annotated domain fragments
- Reports provisional P-I/P-II/P-III classification or uncertainty
- Matplotlib domain/motif visualization
- Protein properties and downloadable JSON report

## Windows setup (Command Prompt)
1. Install Python 3.11 or 3.12.
2. Extract this ZIP, then open Command Prompt in the extracted `SVMP_Classifier` directory.
3. `py -m pip install -r requirements.txt`
4. `py build_database.py` **(internet needed)**
5. `py -m streamlit run app.py`
6. Open the local URL shown in Command Prompt (usually http://localhost:8501).

## Scientific limitations (read before presenting)
1. **This is an educational prototype, not a validated SVMP annotation pipeline.**
2. `motifs.json` contains *selected verified literature patterns*, not every known or unknown motif. Short patterns are not unique to SVMPs.
3. The 15 accessions are a **candidate reference panel**. `build_database.py` downloads actual sequences from UniProt; the ZIP deliberately contains no invented sequences. Inspect download_status.json for failures. Verify class labels, fragment status and UniProt annotation for each entry.
4. Some UniProt records lack a separately annotated cysteine-rich domain. The program does **not** invent C-domain coordinates. Accordingly, P-III predictions can remain provisional even when a P-III reference matches strongly.
5. Domain alignments use manually chosen heuristic thresholds; they have not been calibrated or validated on an independent test set. Alignment coverage is not a probability.
6. P-II D and P-III D-like annotations can both be named 'Disintegrin' in UniProt; the code distinguishes them by the reference protein class. This is a proxy and should be replaced with family-specific profiles for publication-quality analysis.
7. A translated CDS is a **precursor**, possibly containing signal peptide and prodomain; it is not automatically the mature venom protein.
8. mRNA translation assumes a complete forward-strand coding sequence in the chosen frame; arbitrary transcripts need ORF detection and possibly reverse-complement analysis.
9. Proteolytic cleavage, noncanonical P-III subclasses, inactive metalloproteinases and non-SVMP ADAM-family homologs require manual review.
10. To make a rigorous classifier, add curated M/D/D-like/C boundaries, HMMER/Pfam/InterPro annotation, homology-aware validation, negative controls, and independent expert review.

## Provenance
- `database/sources.json`: bibliography and source URLs (including the three user-provided reviews).
- `database/motifs.json`: documented motif patterns with source IDs.
- `database/references.json`: 5 proposed accession candidates for each class.
- `database/reference_records.json`: only generated after downloading actual UniProt entries.

## Suggested presentation wording
"We developed a Python-based literature-assisted prototype that translates coding mRNA, scans conserved motifs, compares sequences to curated SVMP reference proteins and annotated domains, and visualizes evidence for provisional P-I, P-II and P-III classification. Predictions are not a substitute for full domain annotation."
