# Streamlit Community Cloud deployment

1. Extract this ZIP and upload the **contents of `SVMP_Classifier`** to your existing GitHub repository, preserving the `database/` folder and its CSV/JSON files. Do not upload the ZIP as the app entry point.
2. Confirm `app.py`, `analysis.py`, `build_database.py`, `requirements.txt`, and `database/references.json` are in the same repository root.
3. On Streamlit Community Cloud set main file to `app.py` and redeploy/reboot.
4. Open the sidebar and click **Download / refresh UniProt references** once. Wait for completion. Reference counts should be >0, and all three classes should be represented before classification.
5. Test UniProt FASTA P0C6E3 and P34182 (P-II). If either is not P-II, inspect downloaded annotations; do not force the prediction. Test a verified P-III precursor with all three annotated domains.

If UniProt is blocked by your host, run `python build_database.py` on a network-connected PC, then commit generated `database/reference_records.json` and `database/reference_sequences.fasta` to GitHub.

`python -m pytest -q` runs unit/regression tests. Synthetic tests are **not** empirical accuracy validation.
