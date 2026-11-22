# SVMP Explorer — integrated corrected project

**Included:** Streamlit app, Python motif and domain analysis, UniProt downloader, accession candidates, original literature CSV datasets, literature sources, tests and deployment instructions.

**Not included:** pre-downloaded UniProt FASTA sequences, because the build environment could not access UniProt. No sequences have been invented. The app will explicitly report 0 references until you use the sidebar download button (or `python build_database.py`).

Run: `pip install -r requirements.txt` then `streamlit run app.py`. See `DEPLOYMENT.md` and `DATA_PROVENANCE.md`.

Predictions are educational/provisional. For actual research, use independent domain annotation (e.g., InterPro/Pfam) and a held-out curated benchmark.
