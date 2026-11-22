# Data provenance and limitations

The `database/literature/` directory preserves the supplied earlier CSV datasets unchanged, including conflicting or provisional claims. These are **literature evidence tables**, not automatically validated domain annotations.

`database/references.json` is a catalogue of accession **candidates**. It is not a sequence FASTA dataset. `build_database.py` downloads the actual UniProt records and their domain features. Only successfully downloaded entries are counted.

The UniProt screenshots supplied by the user annotate Bilitoxin-1 (P0C6E3) as M 7–203, D 211–291, MGD 273–275. This is a positive P-II check. Atrolysin E/D (P34182) is another P-II check. Neither is a P-III reference.

Class assignments are **provisional homology inferences**. Short motifs, cysteine counts, and the presence of a zinc-binding motif cannot independently establish P-I, P-II or P-III. Mature fragments and unusual P-III subtypes may require manual expert annotation.

Primary review: Olaoba et al. 2020, Toxicon X, DOI:10.1016/j.toxcx.2020.100052. Other source details are in `database/sources.json` and literature CSVs.

**Offline limitation:** The authoring environment could not connect to rest.uniprot.org, so no real UniProt FASTA sequence was fabricated or embedded. On Streamlit Cloud, use the in-app 'Download / refresh UniProt references' control; on a computer, run `python build_database.py`. Deployment must have outbound internet for this initial download.
