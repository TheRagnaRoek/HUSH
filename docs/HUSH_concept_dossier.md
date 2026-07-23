# A Modern, Allele-Aware siRNA/shRNA Design Tool
### Project One-Pager

---

## The Problem

Rational siRNA design has been effectively frozen since 2004. The Ui-Tei, Reynolds, and Amarzguioui rule sets — empirical thermodynamic-asymmetry heuristics — remain the backbone of nearly every modern tool, and the one implementation that meaningfully advanced the field beyond them, **siDirect 2.0**, added seed-duplex Tm-based off-target filtering but has otherwise gone unmaintained: no shRNA output, no chemical-modification guidance, no modern interface, and no path from siRNA candidate to a cloneable construct.

Meanwhile, the tools that *have* been built recently split into two categories that don't overlap with siDirect's role: **efficacy confirmation tools** (DeepSilencer, RNAiSpline, and similar ML models that score a given siRNA rather than generate candidates) and **precomputed genome-scale resources** (Biopredsi's one-siRNA-per-human-gene approach). Neither replaces a flexible, from-sequence design pipeline.

Our project targets a specific, high-value application of this gap: **RNA interference against HLA genes** — relevant to hypoimmunogenic/allogeneic cell engineering and to immune-evasion research more broadly. HLA's extreme allelic polymorphism means a generic design tool cannot simply be pointed at "the HLA gene": conserved-region targeting (for broad/pan-allele knockdown) and allele-specific targeting (for precision knockdown) require sequence-set-aware logic that no current RNAi tool implements.

Notably, this exact problem was recently solved **for CRISPR**: HLA-Knockout (Rutgers, 2026) pulls allele sequences from the IPD-IMGT/HLA database and uses mismatch/PAM-disruption filtering to achieve allele-specific knockout without disturbing non-target alleles. **No RNAi equivalent exists.** Prior HLA-knockdown shRNA work (lentiviral, pan-Class I and allele-specific) demonstrates the biology works — but as bespoke, manually-designed reagents, not a generalizable design tool.

---

## Existing Tool Landscape

| Tool | Category | Strength | Gap (relative to our target) |
|---|---|---|---|
| **siDirect 2.0** | siRNA design | Seed-duplex Tm-based off-target filtering; still the most mechanistically sound classical tool | Unmaintained UI, no shRNA output, no chemistry guidance, no allele-awareness |
| **DSIR** | siRNA design | Purpose-built for Dicer-substrate siRNAs | Narrow scope (Dicer-substrates only), not general-purpose |
| **SSD / siRNA Wizard** | siRNA→shRNA conversion | Template mRNA input directly into mono-/bi-/trivalent shRNA cassettes | Naive templating on top of outdated siRNA cores; no seed-Tm or scaffold-processing optimization |
| **pssRNAit** | Integrated RNAi pipeline (plant) | Chains rule-based + ML + genome-wide off-target scoring into one pipeline | Plant-specific biology and databases; not applicable to mammalian systems |
| **SIREN** (2025) | Integrated dsRNA pipeline | Thermodynamically-informed cumulative off-target modeling, modern & open-source | Built for long dsRNA constructs, not discrete siRNA/shRNA design |
| **DeepSilencer / RNAiSpline** | Efficacy *confirmation* (ML) | State-of-the-art scoring of a *given* candidate | Not generative — require a candidate to already exist |
| **Biopredsi** | Precomputed genome-scale | One scored siRNA per human gene | Fixed, non-customizable; no allele or off-target-set flexibility |
| **HLA-Knockout** (2026) | Allele-specific CRISPR design | Solves the HLA-polymorphism targeting problem — for CRISPR | No RNAi equivalent exists |

---

## Our Approach — How It Fills the Gap

We are building a modernized, open-source successor to siDirect's core logic that closes three gaps simultaneously, rather than one:

1. **Modernized design core** — retain the classical rule filter as a fast pre-screen, paired with a lightweight, interpretable ML scoring layer (logistic regression / gradient boosting trained on public benchmark data such as Huesken) rather than a deep architecture — feasible within our timeline and, per recent comparative work, competitive with deep learning approaches on the same benchmarks.
2. **Native shRNA conversion** — seed-duplex-aware siRNA candidates are carried through to scored shRNA scaffold design (stem-loop vs. miR-30/miR-E-style pri-miRNA context), rather than naive template-stamping.
3. **Allele-awareness** - an MSA-driven module that aligns a user-selected allele panel and scores candidates for either pan-allele conserved-region targeting or single-allele-specific targeting, using seed-duplex thermodynamic differentials. The tool is sequence-agnostic; HLA (via IPD-IMGT/HLA integration) is the proof-of-concept panel that motivated the module, making this the RNAi-appropriate analog of what HLA-Knockout does for CRISPR.

## Conclusion

No existing tool combines a modernized, accessible siRNA design core, native shRNA scaffold conversion, and allele-polymorphism-aware targeting in one open-source pipeline — and no RNAi tool addresses HLA's polymorphism problem at all, despite a directly analogous CRISPR tool already existing and validated wet-lab precedent for HLA shRNA knockdown. This positions the project as a focused, judge-legible contribution: a real, currently-unfilled gap, scoped to what a graduate team can credibly deliver and validate within an iGEM timeline.

---

## References

### Classical siRNA design rules & off-target mechanism
- Ui-Tei, K. et al. (2004). Guidelines for the selection of highly effective siRNA sequences for mammalian and chick RNA interference. *Nucleic Acids Research*, 32(3), 936–948. [Oxford Academic](https://academic.oup.com/nar/article/32/3/936/2904484) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/14769950/)
- Reynolds, A. et al. (2004). Rational siRNA design for RNA interference. *Nature Biotechnology*, 22, 326–330. [Nature](https://www.nature.com/articles/nbt936) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/14758366/)
- Amarzguioui, M. & Prydz, H. (2004). An algorithm for selection of functional siRNA sequences. *Biochemical and Biophysical Research Communications*, 316(4), 1050–1058. DOI: [10.1016/j.bbrc.2004.02.157](https://doi.org/10.1016/j.bbrc.2004.02.157)
- Ui-Tei, K. et al. (2008). Thermodynamic stability and Watson-Crick base pairing in the seed duplex are major determinants of the efficiency of the siRNA-based off-target effect. *Nucleic Acids Research*, 36, 7100–7109. DOI: 10.1093/nar/gkn902 (referenced via [siDirect 2.0, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2791777/))

### Design tools & benchmarks
- Naito, Y. et al. (2004). siDirect: highly effective, target-specific siRNA design software for mammalian RNA interference. *Nucleic Acids Research*, 32, W124–W129.
- Naito, Y., Yoshimura, J., Morishita, S., Ui-Tei, K. (2009). siDirect 2.0: updated software for designing functional siRNA with reduced seed-dependent off-target effect. *BMC Bioinformatics*, 10, 392. [Springer](https://link.springer.com/article/10.1186/1471-2105-10-392) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2791777/) · [Tool](https://sidirect2.rnai.jp/)
- Vert, J.P., Foveau, N., Lajaunie, C., Vandenbrouck, Y. An accurate and interpretable model for siRNA efficacy prediction (DSIR). [Tool: biodev.cea.fr/DSIR](https://biodev.cea.fr/DSIR/)
- Huesken, D. et al. (2005). Design of a genome-wide siRNA library using an artificial neural network (BIOPREDsi). *Nature Biotechnology*, 23, 995–1001. [Nature](https://www.nature.com/articles/nbt1118)
- Silva et al. SSD — a free software for designing multimeric mono-, bi- and trivalent shRNAs. [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197978/)
- InvivoGen siRNA Wizard — siRNA/shRNA target design and hairpin insert tool. [sirnawizard.com](https://sirnawizard.com/)
- pssRNAit: A Web Server for Designing Effective and Specific Plant siRNAs with Genome-Wide Off-Target Assessment. *Plant Physiology*. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7479913/) · [Oxford Academic](https://academic.oup.com/plphys/article/184/1/65/6117836)
- SIREN: Suite for Intelligent RNAi Design and Evaluation of Nucleotide Sequences (2025). *Bioinformatics Advances*. [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.05.26.656188v1.full) · [Oxford Academic](https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag075/8524916)

### ML / modern efficacy prediction
- Monopoli, K.R. et al. (2023). Asymmetric trichotomous partitioning overcomes dataset limitations in building machine learning models for predicting siRNA efficacy. *Molecular Therapy Nucleic Acids*, 33, 93–109.
- Liao, W. & Wang, W. (2025). DeepSilencer: A Novel Deep Learning Model for Predicting siRNA Knockdown Efficiency. [arXiv](https://arxiv.org/html/2503.04200v1) · [ResearchGate](https://www.researchgate.net/publication/389648395_DeepSilencer_A_Novel_Deep_Learning_Model_for_Predicting_siRNA_Knockdown_Efficiency)
- RNAiSpline: A Deep learning model for siRNA efficacy prediction (2026). [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.02.14.705949v1)
- Machine Learning Reveals Intrinsic Determinants of siRNA Efficacy (2025) — logistic regression outperforming DeepSilencer on the Huesken benchmark. [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.08.11.667724v1.full)

### shRNA scaffold engineering
- Fellmann, C. et al. (2013). An Optimized microRNA Backbone for Effective Single-Copy RNAi (miR-E). *Cell Reports*, 5(6), 1704–1713. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2211124713006876) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/24332856/) · [Cell Reports](https://www.cell.com/cell-reports/comments/S2211-1247(13)00687-6)
- Genetically-modified immune cells comprising a microRNA-adapted shRNA (shRNAmiR) — B2M shRNAmiR built on the miR-E backbone. [USPTO patent PDF](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11384335)

### Chemical modification & delivery
- Improving the clinical safety of GalNAc–siRNA conjugates using seed-pairing destabilization. *Nucleic Acids Research*. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9262600/) · [Oxford Academic](https://academic.oup.com/nar/article/50/12/6656/6613918)
- Ribofuranose-based GalNAc-conjugated siRNA enhances liver-targeted delivery and elicits robust RNAi-mediated gene silencing (2025/2026). *Molecular Therapy Nucleic Acids*. [Cell Press](https://www.cell.com/molecular-therapy-family/nucleic-acids/fulltext/S2162-2531(25)00355-5) · [Follow-up](https://www.cell.com/molecular-therapy-family/nucleic-acids/fulltext/S2162-2531(26)00042-9)
- Advanced siRNA Designs Further Improve In Vivo Performance of GalNAc-siRNA Conjugates (ESC chemistry). [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5910670/)

### HLA-specific precedents
- Mattivi, C. et al. (2026). HLA-Knockout: Enabling Allele-Specific Knockout of HLA Class I Genes for Immunogenic Engineering. *Tissue Antigens*. DOI: [10.1111/tan.70548](https://doi.org/10.1111/tan.70548) · [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12820913/) · [Tool: hlaknockout.rutgers.edu](https://hlaknockout.rutgers.edu)
- IPD-IMGT/HLA Database. [ebi.ac.uk/ipd/imgt/hla](https://www.ebi.ac.uk/ipd/imgt/hla/)
- Class-, gene-, and group-specific HLA silencing by lentiviral shRNA delivery. *Journal of Molecular Medicine*. [PubMed](https://pubmed.ncbi.nlm.nih.gov/16520945/) · [Springer](https://link.springer.com/article/10.1007/s00109-005-0024-2)
- Suppression of HLA Expression by Lentivirus-mediated Gene Transfer of siRNA Cassettes and In Vivo Chemoselection to Enhance Hematopoietic Stem Cell Transplantation. [PMC](https://ncbi.nlm.nih.gov/pmc/articles/PMC2938774)
- Selective HLA knockdown and PD-L1 expression prevent allogeneic CAR-NK cell rejection and enhance safety and anti-tumor responses in xenograft mice (2025). *Nature Communications*. [Nature](https://www.nature.com/articles/s41467-025-63863-8)
