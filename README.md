# Mapping marine aquaculture infrastructure with AlphaEarth Foundations embeddings

Code associated with the study:

**Metcher, S.C., Christofidis, M., Sievers, M. & Kuempel, C.D.  
*Mapping marine aquaculture infrastructure with AlphaEarth Foundations embeddings*.**

This repository contains the Google Earth Engine (GEE) and R scripts used to evaluate AlphaEarth Foundations (AEF) embeddings for mapping marine aquaculture infrastructure under limited-reference conditions.

## Overview

The analysis uses the Google Satellite Embedding V1 Annual dataset:

`GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`

AEF provides 64-dimensional annual embedding vectors at 10 m spatial resolution.

Three 30 × 30 km study areas were analysed:

- Melinka–Repollal, Guaitecas Archipelago, Chile
- Finnøy–Fister, Ryfylke, Norway
- Huon Estuary–D’Entrecasteaux Channel, Tasmania, Australia

Models were developed using 2024 AEF embeddings and applied to the corresponding 2025 embeddings for temporally held-out evaluation.

## Classification approaches

The repository implements three approaches:

1. **Cosine similarity**
   - Regional aquaculture reference embeddings are used to calculate cosine similarity across each study area.
   - Thresholds corresponding to 90%, 95% and 99% retention of 2024 aquaculture reference pixels are evaluated.

2. **Standard Random Forest (SRF)**
   - Uses the 64 AEF embedding dimensions directly as predictor variables.

3. **Class-Relationship Random Forest (CR-RF)**
   - Extends the standard RF with derived variables describing:
     - similarity to aquaculture, marine-background and confounding-marine classes;
     - separation between competing class similarities; and
     - similarity to multiple aquaculture prototypes derived using k-means clustering.

## Reference and validation data

Each regional model was developed using 20 reference sites for each of three classes:

- aquaculture;
- marine background; and
- confounding marine.

Reference labels were developed from 2024 imagery. Independent aquaculture validation footprints were delineated from 2025 imagery and were not used during model development.

Classification and evaluation were restricted to water pixels using the ESA WorldCover 2021 water class.

## Software

The workflow uses:

- [Google Earth Engine](https://earthengine.google.com/) for AEF access, spatial analysis and classification;
- R for prototype selection, feature assessment and supporting analyses.
