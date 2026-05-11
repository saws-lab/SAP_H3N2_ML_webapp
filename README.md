---
title: SAP H3N2
sdk: streamlit
sdk_version: 1.27.2
python_version: 3.8.12
app_file: app.py
---

# SAP_H3N2_ML Web Application

Streamlit web application for seasonal antigenic prediction of influenza A virus (IAV) H3N2 using the AdaBoost model from the SAP_H3N2_ML project.

This repository provides an interactive interface to predict normalized hemagglutination titer (NHT)-based antigenic differences between virus–antiserum pairs using HA1 amino-acid sequences and optional metadata.

## Related research and upstream repository

- Paper: [Shah et al., *Seasonal antigenic prediction of influenza A H3N2 using machine learning* (Nature Communications, 2024)](https://doi.org/10.1038/s41467-024-47862-9)
- Upstream research/code repository: [saws-lab/SAP_H3N2_ML](https://github.com/saws-lab/SAP_H3N2_ML)
- Live app (Streamlit Community Cloud): [sap-h3n2-ml.streamlit.app](https://sap-h3n2-ml.streamlit.app/)
- Live app (Hugging Face Space): [huggingface.co/spaces/sawshah/SAP_H3N2](https://huggingface.co/spaces/sawshah/SAP_H3N2)

## What this webapp does

- Loads season-specific pretrained AdaBoost models (`trained_model/`)
- Encodes pairwise HA1 sequence differences using amino-acid mutation matrix `GIAG010101` from AAindex
- Optionally incorporates metadata through one-hot encoding:
	- `virus`
	- `serum`
	- `virusPassCat`
	- `serumPassCat`
- Supports:
	- single pair prediction (interactive form)
	- batch prediction from CSV upload with downloadable output

## Repository structure

- `app.py`: Streamlit app entry point
- `app_utilities.py`: input validation, preprocessing, sequence encoding, file handling
- `aaindex.py`: AAindex parser/helper
- `AAindex/aaindex2`: amino-acid index source file used for mutation matrix lookup
- `trained_model/`: serialized AdaBoost models and metadata encoders for seasons
- `requirements.txt`: Python dependencies

## Supported prediction seasons

The app allows selecting test seasons:

- `2018NH`, `2018SH`
- `2019NH`, `2019SH`
- `2020NH`, `2020SH`
- `2021NH`

For each selected test season, the app loads the model trained up to the previous season.

## Input requirements

### Single prediction mode

Required:

- `virusSeq`: HA1 sequence of virus isolate
- `serumSeq`: HA1 sequence of antiserum/reference isolate

Optional:

- `virusName`, `serumName`
- `virusPassage`, `serumPassage`

### Batch upload mode (CSV)

CSV must include at least:

- `virusSeq`
- `serumSeq`

Optional columns:

- `virusName`
- `serumName`
- `virusPassage`
- `serumPassage`

Each row is treated as one virus–antiserum pair.

### Sequence validation

- HA1 sequence length must be **329** for both virus and serum
- Allowed amino-acid characters: 20 standard amino acids plus `X`
- Invalid rows are reported with their pair index in batch mode

## Local setup

### Option 1: Conda (aligned with upstream project)

1. Create environment (Python 3.8.12):

	 ```bash
	 conda create --name SAP_H3N2_ML_webapp python=3.8.12
	 conda activate SAP_H3N2_ML_webapp
	 ```

2. Install dependencies:

	 ```bash
	 pip install -r requirements.txt
	 pip install streamlit==1.27.2
	 ```

### Option 2: venv

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
# source .venv/bin/activate

pip install -r requirements.txt
pip install streamlit==1.27.2
```

## Run the app

```bash
streamlit run app.py
```

Then open the local URL shown in your terminal (typically `http://localhost:8501`).

## Output

- Single mode: displays predicted NHT for one pair
- Batch mode: displays predictions and enables CSV download with columns including:
	- `virusName`
	- `virusPassage`
	- `serumName`
	- `serumPassage`
	- `predicted NHT`
	- `virusSeq`
	- `serumSeq`

## Notes

- This repo is the deployment-focused webapp companion to the main SAP_H3N2_ML research repo.
- The app is currently available on both Streamlit Community Cloud and Hugging Face Spaces.
- If you use GitHub Actions deployment to Hugging Face Space, configure secret `HF_TOKEN` in GitHub.

## Troubleshooting

For questions/comments about the SAP_H3N2 project:

- [awais.shah@unimelb.edu.au](mailto:awais.shah@unimelb.edu.au)
- [sawshah@connect.ust.hk](mailto:sawshah@connect.ust.hk)

Hugging Face Spaces README front matter reference:
[https://huggingface.co/docs/hub/spaces-config-reference](https://huggingface.co/docs/hub/spaces-config-reference)
