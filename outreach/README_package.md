# Hot-mitochondria proteome reanalysis: code and results

Contents
- `REPORT.md`: full technical write-up (also as `technical_report.pdf`)
- `scripts/`: Python pipeline, numbered in run order (03_biomart.py was abandoned; BioMart returned HTTP 405)
- `results/`: every intermediate and final number quoted in the report
- `figures/`: main figure (wide and 2x2 versions)

Inputs (all public; not included because of size, and fetched by scripts 01–06):
- Meltome Atlas: Jarzab et al. 2020 Nat Methods, Supplementary Tables S3/S4, plus the FLIP benchmark copy for curves
- MitoCarta 3.0 (Broad Institute), human
- UniProt REST (annotations, sequences)
- Alliance of Genome Resources orthology file
- TemStaPro (Pudžiuvelytė et al. 2024) with ProtT5-XL-U50 embeddings

Main estimates: `16_revised.py`; per-condition CIs and figure: `18_conditions_ci.py`, `14_figures.py`; inner membrane vs matrix by lysate/cells: `19_subcompartments.py`;
paralogs: `17_paralogs_curated.py`; sequence thermometer: `13_temstapro.py`, `15_temstapro_paired.py`.
Analysis carried out with AI assistance (Anthropic's Claude) and independently reviewed by OpenAI Codex.
Contact: Paras Chopra, paras@paraschopra.com
