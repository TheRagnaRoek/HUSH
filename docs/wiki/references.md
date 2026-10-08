# HUSH - combined reference list

Every reference cited anywhere in the repo, de-duplicated and grouped by what it
is used for. Collected from `docs/HUSH_concept_dossier.md`, `docs/design_notes.md`,
`docs/seed_offtarget_design_note.md`, `docs/offtarget_design_note.md`,
`docs/shrna_design_note.md`, `msa/discrimination.py` and `shrna/loops.py`.

Each entry carries:
- the citation as the repo states it,
- the best identifier available (DOI / PMID / arXiv / URL),
- **Used for:** what HUSH actually takes from it, so you can cite it in the right sentence.

For bulk import, paste the contents of `references_zotero.txt` into Zotero's
*Add Item by Identifier* box (the magic-wand icon) - it accepts many identifiers at
once and will fetch complete, correct metadata. That is more reliable than importing a
hand-written BibTeX file, which is why one is not provided here. The handful of entries
with no DOI/PMID are listed at the end of that file to capture with the Zotero
browser connector instead.

**Verify before citing** - flagged inline with (!): a few entries have a
year/volume mismatch or an uncertain journal name in the repo's own text.

---

## 1. Classical siRNA design rules (the 2004 backbone)

HUSH reimplements all three faithfully in `core.rules`.

- **Ui-Tei, K., Naito, Y., Takahashi, F., Haraguchi, T., Ohki-Hamazaki, H., Juni, A., Ueda, R. & Saigo, K. (2004).** Guidelines for the selection of highly effective siRNA sequences for mammalian and chick RNA interference. *Nucleic Acids Research* 32(3), 936-948. PMID **14769950**
  - **Used for:** the strict 4-rule pass/fail filter. Also the source of the GC-rich tension documented in the seed note (rule (ii) keeps the passenger seed GC-rich, rule (iii) pulls the guide seed Tm down).
- **Reynolds, A., Leake, D., Boese, Q., Scaringe, S., Marshall, W.S. & Khvorova, A. (2004).** Rational siRNA design for RNA interference. *Nature Biotechnology* 22, 326-330. DOI **10.1038/nbt936** - PMID **14758366**
  - **Used for:** the score out of 10, implemented following the original algorithm rather than the common "8 equal binary criteria" simplification. Also the internal-structure criterion (Tm <= 20 C) that `core.hairpin` supplies.
- **Amarzguioui, M. & Prydz, H. (2004).** An algorithm for selection of functional siRNA sequences. *Biochemical and Biophysical Research Communications* 316(4), 1050-1058. DOI **10.1016/j.bbrc.2004.02.157**
  - **Used for:** the additive motif score plus the raw A/U differential ("DT3"). Its Table 3 (46 siRNAs) is HUSH's exact-reproduction validation fixture.

## 2. RNAi mechanism: strand selection and seed-mediated off-target

- **Khvorova, A., Reynolds, A. & Jayasena, S.D. (2003).** Functional siRNAs and miRNAs exhibit strand bias. *Cell* 115(2), 209-216. DOI **10.1016/s0092-8674(03)00801-8** - PMID **14567918**
  - **Used for:** the mechanistic basis of `core.asymmetry` - RISC keeps the strand whose 5' end sits at the duplex's looser end. Khvorova computed internal-stability profiles over 5-nt windows, one reason HUSH's asymmetry window is a parameter rather than a constant.
- **Schwarz, D.S., Hutvagner, G., Du, T., Xu, Z., Aronin, N. & Zamore, P.D. (2003).** Asymmetry in the assembly of the RNAi enzyme complex. *Cell* 115(2), 199-208. DOI **10.1016/s0092-8674(03)00759-1** - PMID **14567917**
  - **Used for:** the same, from the complex-assembly side; Schwarz emphasised the few terminal base pairs.
- **Ui-Tei, K., Naito, Y., Nishi, K., Juni, A. & Saigo, K. (2008).** Thermodynamic stability and Watson-Crick base pairing in the seed duplex are major determinants of the efficiency of the siRNA-based off-target effect. *Nucleic Acids Research* 36(22), 7100-7109. DOI **10.1093/nar/gkn902**
  - **Used for:** the whole premise of `core.seed` - off-target strength tracks seed-duplex stability, so a weakly-pairing seed is largely off-target free even where its matches exist.
- **Salomon, W.E., Jolly, S.M., Moore, M.J., Zamore, P.D. & Serebrov, V. (2015).** Single-molecule imaging reveals that Argonaute reshapes the binding properties of its nucleic acid guides. *Cell* 162(1), 84-95. DOI **10.1016/j.cell.2015.06.029** - PMID **26140592**
  - **Used for:** the critique of the seed filter's physics - the seed duplex does not form in solution, and seed pairing proceeds 50-200x faster than annealing of the equivalent naked RNAs.
- **Wee, L.M., Flores-Jasso, C.F., Salomon, W.E. & Zamore, P.D. (2012).** Argonaute divides its RNA guide into domains with distinct functions and RNA-binding properties. *Cell* 151(5), 1055-1067. DOI **10.1016/j.cell.2012.10.036** - PMID **23178124**
  - **Used for:** the same critique; also underpins the treatment of guide position 1 as an anchor held in Argonaute's MID pocket that does not pair with the target.

## 3. Allele discrimination by single mismatches (the MSA stage's evidence base)

These three are the entire published basis of `msa.discrimination`, used in this order.

- **Huang, H., Qiao, R., Zhao, D., Zhang, T., Li, Y., Yi, F. et al. (2009).** Profiling of mismatch discrimination in RNAi enabled rational design of allele-specific siRNAs. *Nucleic Acids Research* 37(22), 7560-7569. DOI **10.1093/nar/gkp835**
  - **Used for:** Table 1, the lookup keyed by the exact base change at the exact position - the most specific evidence available, and the first tier of the ranking tuple (`H+` / `H-` / `H?` in the `disc` readout).
- **Schwarz, D.S., Ding, H., Kennington, L., Moore, J.T., Schelter, J., Burchard, J. et al. (2006).** Designing siRNA that distinguish between genes that differ by a single nucleotide. *PLoS Genetics* 2(9), e140. DOI **10.1371/journal.pgen.0020140**
  - **Used for:** the general principles that order what Huang endorses - purine:purine mismatches discriminate best, and mismatches 3' of the seed discriminate robustly while seed mismatches are overwhelmed by high siRNA concentration. Also the source of the "consistently discriminating positions" set, including the P5-in-the-seed anomaly HUSH reproduces intact.
- **Ohnishi, Y., Tamura, Y., Yoshida, M., Tokunaga, K. & Hohjoh, H. (2008).** Enhancement of allele discrimination by introduction of nucleotide mismatches into siRNA in allele-specific gene silencing by RNAi. *PLoS ONE* 3(5), e2248. DOI **10.1371/journal.pone.0002248**
  - **Used for:** the rule that sits above both - a second mismatch, even inside the seed, markedly improves discrimination, which is why mismatch count is the top of the lexicographic tuple.

## 4. Nearest-neighbor thermodynamics

- **Xia, T., SantaLucia, J., Burkard, M.E., Kierzek, R., Schroeder, S.J., Jiao, X., Cox, C. & Turner, D.H. (1998).** Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs. *Biochemistry* 37(42), 14719-14735. DOI **10.1021/bi9809425**
  - **Used for:** `core.thermo`'s default parameter set; the asymmetry and hairpin readouts; Reynolds' internal-repeat criterion is pinned to it.
- **Freier, S.M., Kierzek, R., Jaeger, J.A., Sugimoto, N., Caruthers, M.H., Neilson, T. & Turner, D.H. (1986).** Improved free-energy parameters for predictions of RNA duplex stability. *PNAS* 83(24), 9373-9377. DOI **10.1073/pnas.83.24.9373** - PMID **2432595**
  - **Used for:** the second, explicitly-named parameter set - required because siDirect's 21.5 C seed threshold was calibrated on Freier numbers (its -10.8 initiation constant is Freier's).
- **SantaLucia, J. & Hicks, D. (2004).** The thermodynamics of DNA structural motifs. *Annual Review of Biophysics and Biomolecular Structure* 33, 415-440. DOI **10.1146/annurev.biophys.32.110601.141800** - PMID **15139820**
  - **Used for:** the Tm formulation in `core.thermo`.
- **Groebe, D.R. & Uhlenbeck, O.C. (1988).** Characterization of RNA hairpin loop stability. *Nucleic Acids Research* 16(24), 11725-11735. DOI **10.1093/nar/16.24.11725** - PMID **3211748**
  - **Used for:** `core.hairpin` - the concentration-independence of intramolecular melting (identical profiles from 1 to 10 uM, hence `Tm = dH/dS`) and the Table 3 loop destabilisation values, composition-weighted and size-interpolated.

## 5. siRNA duplex geometry

- **Elbashir, S.M., Harborth, J., Lendeckel, W., Yalcin, A., Weber, K. & Tuschl, T. (2001).** Duplexes of 21-nucleotide RNAs mediate RNA interference in cultured mammalian cells. *Nature* 411, 494-498. DOI **10.1038/35078107**
  - **Used for:** the canonical duplex - a 19-nt paired core with symmetric 2-nt 3' overhangs, i.e. 21-nt synthesized strands.
- **Elbashir, S.M., Martinez, J., Patkaniowska, A., Lendeckel, W. & Tuschl, T. (2001).** Functional anatomy of siRNAs for mediating efficient RNAi in Drosophila melanogaster embryo lysate. *EMBO Journal* 20, 6877-6888. DOI **10.1093/emboj/20.23.6877**
  - **Used for:** establishing the 2-nt 3' overhang as most efficient, and that 2'-deoxy (dTdT) is tolerated at the 3' end - the basis of the `dtdt` overhang policy.

## 6. Existing design tools (the prior-art landscape)

- **Naito, Y., Yamada, T., Ui-Tei, K., Morishita, S. & Saigo, K. (2004).** siDirect: highly effective, target-specific siRNA design software for mammalian RNA interference. *Nucleic Acids Research* 32, W124-W129. DOI **10.1093/nar/gkh442** - PMID **15215364** - PMC **PMC441580**
  - **Used for:** the original siDirect, the tool HUSH positions itself as the successor to.
- **Naito, Y., Yoshimura, J., Morishita, S. & Ui-Tei, K. (2009).** siDirect 2.0: updated software for designing functional siRNA with reduced seed-dependent off-target effect. *BMC Bioinformatics* 10, 392. DOI **10.1186/1471-2105-10-392**
  - **Used for:** the single most-cited reference in the repo. The exact seed-Tm formula and all its constants; the three-step selection scheme; the transcriptome off-target step; the "cannot design for 5.3% of RefSeq mRNAs" caveat; the Figure 4 low-Tm/high-abundance anticorrelation; and the published figures HUSH validates against (4488 of 16384 heptamers, -12 to 60 C). Tool at sidirect2.rnai.jp.
- **Vert, J.P., Foveau, N., Lajaunie, C. & Vandenbrouck, Y.** An accurate and interpretable model for siRNA efficacy prediction (**DSIR**). Tool: biodev.cea.fr/DSIR/
  - **Used for:** the prior-art table (purpose-built for Dicer-substrate siRNAs, narrow scope). Also named in the seed note as a tool that computes raw seed frequency but not Tm.
- **Huesken, D., Lange, J., Mickanin, C., Weiler, J., Asselbergs, F., Warner, J. et al. (2005).** Design of a genome-wide siRNA library using an artificial neural network (**BIOPREDsi**). *Nature Biotechnology* 23, 995-1001. DOI **10.1038/nbt1118**
  - **Used for:** two roles - the precomputed genome-scale category in the prior-art table, and the public benchmark dataset the planned ML scoring layer would train on.
- **Silva et al. SSD** - a free software for designing multimeric mono-, bi- and trivalent shRNAs. PMC **PMC7197978**
  - **Used for:** the siRNA-to-shRNA precedent HUSH improves on. Key point taken from it: the paper publishes neither loop sequences nor folding rules; it stamps picks into cassettes without scoring the hairpin.
- **InvivoGen siRNA Wizard** - siRNA/shRNA target design and hairpin insert tool. https://sirnawizard.com/
  - **Used for:** the same naive-templating category in the prior-art table.
- **pssRNAit**: A Web Server for Designing Effective and Specific Plant siRNAs with Genome-Wide Off-Target Assessment. *Plant Physiology* 184(1), 65. PMC **PMC7479913**
  - **Used for:** the integrated-pipeline category - chains rule-based + ML + genome-wide off-target scoring, but plant-specific.
- **SIREN** (Vargas-Mejia, P. & Vega-Arreguin, J., 2026): Suite for Intelligent RNAi Design and Evaluation of Nucleotide Sequences. *Bioinformatics Advances* 6(1), vbag075. DOI **10.1093/bioadv/vbag075** - PMID **41890811** - PMC **PMC13012889**. Preprint DOI **10.1101/2025.05.26.656188**
  - **Used for:** the modern open-source comparator - thermodynamically-informed cumulative off-target modelling, but built for long dsRNA constructs.
  - Note: the repo's own text dates it 2025 (the preprint year); the journal version was published 2026-03-16. Cite as 2026.

## 7. Machine-learning efficacy prediction

Context for the (not yet implemented) ML scoring layer and for the prior-art argument that
modern tools confirm efficacy rather than generate candidates.

- **Monopoli, K.R., Korkin, D. & Khvorova, A. (2023).** Asymmetric trichotomous partitioning overcomes dataset limitations in building machine learning models for predicting siRNA efficacy. *Molecular Therapy Nucleic Acids* 33, 93-109. DOI **10.1016/j.omtn.2023.06.010** - PMID **37456778**
- **Liao, W. & Wang, W. (2025).** DeepSilencer: A Novel Deep Learning Model for Predicting siRNA Knockdown Efficiency. arXiv **2503.04200**
  - **Used for:** the "efficacy confirmation, not generative" category.
- **RNAiSpline**: A Deep learning model for siRNA efficacy prediction (2026). bioRxiv DOI **10.64898/2026.02.14.705949**
  - **Used for:** the same category, as the most recent example.
- **Machine Learning Reveals Intrinsic Determinants of siRNA Efficacy** (2025). bioRxiv DOI **10.1101/2025.08.11.667724**
  - **Used for:** the key justification of pillar 1 - logistic regression outperforming DeepSilencer on the Huesken benchmark, i.e. a lightweight interpretable model is competitive with deep learning.

## 8. shRNA scaffold engineering

- **Brummelkamp, T.R., Bernards, R. & Agami, R. (2002).** A system for stable expression of short interfering RNAs in mammalian cells. *Science* 296, 550-553. DOI **10.1126/science.1068999**
  - **Used for:** the canonical pSUPER stem-loop design HUSH builds on (Pol III promoter, inverted repeat - loop - inverted repeat, poly-T terminator; sense-loop-antisense convention) and the `psuper` loop UUCAAGAGA in the registry.
- **Anthony, R.V. & Cantlon, J.D. (2007).** Ribonucleic acid interference: a new approach to the in vivo study of gene function. *Journal of Animal Science* 85(13 Suppl), E18-19. DOI **10.2527/jas.2006-486**
  - **Used for:** general review context in the shRNA note.
  - Note: PubMed lists online publication on 2006-10-13 but the issue (vol. 85) is 2007; cited here as 2007. The shRNA note was updated to match.
- **Ge, Q., Ilves, H., Dallas, A., Kumar, P., Shorenstein, J., Kazakov, S.A. & Johnston, B.H. (2010).** Minimal-length short hairpin RNAs: the relationship of structure and RNAi activity. *RNA* 16(1), 106-117. DOI **10.1261/rna.1894510** - PMID **19952116** - PMC **PMC2802021**
  - **Used for:** the key stem-length structure-activity study - stems <= 19 bp ("sshRNA") behave distinctly and can be extremely potent, which is what puts HUSH's 19-mer cores on a processing boundary.
  - Note: advance-online 2009-12-01, issue 16(1) January 2010; cited here as 2010. The shRNA note was updated to match.
- **Jensen, S.M.R., Schmitz, A., Pedersen, F.S., Kjems, J. & Bramsen, J.B. (2012).** Functional selection of shRNA loops from randomized retroviral libraries. *PLoS ONE* 7(8), e43095. DOI **10.1371/journal.pone.0043095**
  - **Used for:** HUSH's **default loop** UGUGCUU, chosen because it gives robust knockdown with the least stem-sequence dependence - the safest property for an automated tool that does not know the stem in advance.
- **Wang, X.J., Li, Y., Huang, H., Zhang, X.J., Xie, P.W., Hu, W., Li, D.D. & Wang, S.Q. (2013).** A simple and robust vector-based shRNA expression system used for RNA interference. *PLoS ONE* 8(2), e56110. DOI **10.1371/journal.pone.0056110**
  - **Used for:** the loop-length evidence - 8-nt loops best overall (TTCTAGAA, TTGGCCAA), a 4-nt loop dropping below 50% knockdown, and TTCAAGAGA beating CTCGAG.
  - Note: earlier repo text attributed this paper to "Liu, Y.P. et al."; PubMed lists Wang XJ as first author (PMID 23405258). Corrected throughout.
- **Gao, Z., Berkhout, B. & Herrera-Carrillo, E. (2019).** Boosting AgoshRNA activity by optimized 5'-terminal nucleotide selection. *RNA Biology* 16(7), 890-898. DOI **10.1080/15476286.2019.1599259**
  - **Used for:** the AgoshRNA boundary condition - a very short stem with a small loop is processed by Argonaute-2 rather than Dicer, yielding only the guide strand.
- **Berkhout, B. & Herrera-Carrillo, E. (2021).** Design and evaluation of AgoshRNAs with 3'-terminal HDV ribozymes. *Methods in Molecular Biology* 2167, 225-252. DOI **10.1007/978-1-0716-0716-9_13**
  - **Used for:** the same boundary condition.
- **Fellmann, C., Hoffmann, T., Sridhar, V., Hopfgartner, B., Muhar, M., Roth, M. et al. (2013).** An optimized microRNA backbone for effective single-copy RNAi (**miR-E**). *Cell Reports* 5(6), 1704-1713. PMID **24332856**
  - **Used for:** context only. Explicitly **out of scope** - the miR-30/miR-E pri-miRNA route was descoped on 2026-09-12 as needing research effort and wet-lab infrastructure the team does not have.
- **Genetically-modified immune cells comprising a microRNA-adapted shRNA (shRNAmiR)** - B2M shRNAmiR built on the miR-E backbone. US patent **11384335**. https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11384335
  - **Used for:** evidence that the miR-E route is being commercialised for exactly this application area.
- **Moffat, J., Grueneberg, D.A., Yang, X., Kim, S.Y., Kloepfer, A.M., Hinkle, G. et al. (2006).** A lentiviral RNAi library for human and mouse genes applied to an arrayed viral high-content screen. *Cell* 124(6), 1283-1298. DOI **10.1016/j.cell.2006.01.040** - PMID **16564017**
  - **Used for:** the third registry loop, `trc` CUCGAG (CTCGAG as DNA, embeds an XhoI site), from the pLKO.1 / RNAi Consortium (TRC) library. The loop and flanks (21-bp stem, 6-nt XhoI loop, 5' CCGG, 3' TTTTT) are documented in the TRC lentiviral shRNA manual (Horizon/Dharmacon); Addgene asks users of pLKO.1 to cite this paper. The paper's full text was not read; the abstract does not mention the loop.
- **Stewart, S.A., Dykxhoorn, D.M., Palliser, D., Mizuno, H., Yu, E.Y., An, D.S., Sabatini, D.M., Chen, I.S.Y., Hahn, W.C., Sharp, P.A., Weinberg, R.A. & Novina, C.D. (2003).** Lentivirus-delivered stable gene silencing by RNAi in primary cells. *RNA* 9(4), 493-501. DOI **10.1261/rna.2192803** - PMID **12649500** - PMC **PMC1370415**
  - **Used for:** precedent only - its GFP1 hairpin oligo already uses CTCGAG between the stem arms (retroviral vector; the paper gives no loop sequence for its lentiviral constructs).

## 9. Off-target search infrastructure

- **Langmead, B., Trapnell, C., Pop, M. & Salzberg, S.L. (2009).** Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology* 10, R25. DOI **10.1186/gb-2009-10-3-r25**
  - **Used for:** bowtie1, the engine behind `core.offtarget`'s transcriptome-scale backend - chosen because its `-v` mode (end-to-end, at most N mismatches, all hits) is the open offline equivalent of the GGGenome k-mismatch search siDirect relies on.
- **O'Leary, N.A. et al. (2016).** Reference sequence (RefSeq) database at NCBI. *Nucleic Acids Research* 44(D1), D733-D745. DOI **10.1093/nar/gkv1189**
  - **Used for:** the pinned reference transcriptome, chosen over Ensembl/GENCODE because its NM_/NR_ accessions match the benchmark fixtures, making on-target exclusion by accession and the siDirect column comparison exact.

## 10. HLA-specific precedents

- **Mattivi, C., Wang, S., Ji, L., Xiao, Q. & Cao, J. (2026).** HLA-Knockout: Enabling Allele-Specific Knockout of HLA Class I Genes for Immunogenic Engineering. *HLA* 107(1), e70548. DOI **10.1111/tan.70548** - PMC **PMC12820913** - tool at hlaknockout.rutgers.edu
  - **Used for:** the central framing of the whole project - the same allele-polymorphism targeting problem, solved for CRISPR but with no RNAi equivalent. HUSH's allele stage is explicitly positioned as the RNAi analog.
  - Note: CrossRef confirms the journal is *HLA* (formerly *Tissue Antigens*), published online 2026-01-21. The dossier was updated to match.
- **IPD-IMGT/HLA Database.** https://www.ebi.ac.uk/ipd/imgt/hla/
  - **Used for:** the allele-panel data source; also the reason the MSA stage can assume unambiguous sequences (no IUPAC resolution) and must scale to thousands of members per locus.
  - Note: only the website is cited in the repo. The database has its own citable paper - worth adding for the wiki.
- **Class-, gene-, and group-specific HLA silencing by lentiviral shRNA delivery.** *Journal of Molecular Medicine*. DOI **10.1007/s00109-005-0024-2** - PMID **16520945**
  - **Used for:** the wet-lab precedent that HLA shRNA knockdown works - but as bespoke, manually designed reagents, which is the gap HUSH fills.
- **Suppression of HLA Expression by Lentivirus-mediated Gene Transfer of siRNA Cassettes and In Vivo Chemoselection to Enhance Hematopoietic Stem Cell Transplantation.** PMC **PMC2938774**
  - **Used for:** the same precedent, in a transplantation context.
- **Selective HLA knockdown and PD-L1 expression prevent allogeneic CAR-NK cell rejection and enhance safety and anti-tumor responses in xenograft mice (2025).** *Nature Communications*. DOI **10.1038/s41467-025-63863-8**
  - **Used for:** the current, high-profile application case - HLA knockdown for allogeneic cell therapy.

## 11. Chemical modification and delivery (context, not implemented)

Cited in the dossier as part of the "no chemistry guidance" gap in existing tools. HUSH
does not currently implement any modification guidance, so these support the problem
statement rather than a feature.

- **Improving the clinical safety of GalNAc-siRNA conjugates using seed-pairing destabilization.** *Nucleic Acids Research* 50(12), 6656. PMC **PMC9262600**
  - **Used for:** independent evidence that seed-pairing destabilization is a real clinical safety lever - supports the seed filter's premise from the therapeutic side.
- **Ribofuranose-based GalNAc-conjugated siRNA enhances liver-targeted delivery and elicits robust RNAi-mediated gene silencing (2025/2026).** *Molecular Therapy Nucleic Acids*. Cell Press article S2162-2531(25)00355-5, plus a follow-up S2162-2531(26)00042-9.
- **Advanced siRNA Designs Further Improve In Vivo Performance of GalNAc-siRNA Conjugates** (ESC chemistry). PMC **PMC5910670**

## 12. Cited by name but with no reference in the repo

Add these before the wiki goes live, or remove the mention.

- **GGGenome** - the k-mismatch full-text search siDirect depends on, discussed at length in `offtarget_design_note.md` with no citation. It has its own publication from the same group as siDirect.
- **Gusfield's center-star construction** - `msa.align_star` is built on it and names the factor-2 sum-of-pairs bound, but cites nothing. The primary source is Gusfield's work on multiple sequence alignment with guaranteed error bounds; the textbook *Algorithms on Strings, Trees, and Sequences* is the usual citation.
- **Dharmacon siDESIGN** - named in `seed_offtarget_design_note.md` as a tool computing raw seed frequency. Cite the product page if it stays in the text.
- **bowtie's Artistic License / bowtie2's GPLv3** - the licensing argument is made in `offtarget_design_note.md`; link the license texts if the Compatibility section repeats it.
- **RefSeq `H_sapiens/mRNA_Prot` dataset** - the specific NCBI FTP path the index builder uses. Worth linking directly alongside the O'Leary citation.
