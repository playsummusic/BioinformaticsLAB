# Comprehensive Detailed Review: RNA-Seq Differential Expression and Functional Profiling

## Executive Summary

This document provides an academic-level comprehensive treatment of differential gene expression (DGE) analysis and functional profiling, with specific application to the queuosine tRNA modification study in *Escherichia coli* K-12 MG1655. The analysis integrates RNA-seq data processing, statistical testing with DESeq2, quality control assessment, and Gene Ontology enrichment analysis using clusterProfiler.

---

## SECTION I: BACKGROUND AND METHODOLOGY

### Part 1: Queuosine - Biological Significance

#### 1.1 What is Queuosine?

Queuosine (Q) is a rare hypermodified nucleoside with the following characteristics:

**Chemical Structure and Position:**
- Pseudouridine derivative with unusual ribose-base connectivity
- Located at wobble position (position 34) of tRNA anticodons
- Present in tRNAs with GUN anticodons (primarily histidine tRNAs in bacteria)
- Only synthesized in bacteria and lower eukaryotes (may be a bacterial vitamin precursor)

**Biosynthesis Pathway:**
- Encoded by *tgt* gene (tRNA-guanine transglycosylase)
- Multi-step modification involving:
  - QueA protein (queuosine biosynthesis enzyme A)
  - QueB, QueC, QueD, QueE, QueF, QueG proteins
  - Sequential enzymatic modifications of guanosine in tRNA

**Codon Pairing Rules:**
The wobble position (position 34) determines which codons can be read:

| tRNA Position 34 | Codon Position 3 Pair | mRNA Codons Read |
|-----------------|----------------------|-----------------|
| G (normal) | C (Watson-Crick) | CAC (His) only |
| Q (modified) | C (standard) | CAC (His) |
| Q (modified) | U (wobble) | CAU (His) ALSO |

This dual recognition allows Q to decode both CAC and CAU codons (histidine), improving translation efficiency and reducing errors.

#### 1.2 Role in Translation and Metal Homeostasis

**Translation Function:**
- Improves translation fidelity (reduces frameshift errors)
- Enhances translation efficiency (faster decoding)
- Affects codon-anticodon stability
- Particularly important for histidine codons in metal-binding proteins

**Metal Homeostasis Connection:**
- Many metal-binding proteins require histidines (Zn²⁺, Fe²⁺, Fe³⁺, Cu²⁺, Ni²⁺)
- Q-deficiency → impaired translation of metal-binding domains
- Results in compromised metal homeostasis
- Cells compensate by upregulating metal transporters and regulators

**Pleiotropic Effects of Q Deficiency:**
- Enhanced nickel resistance (paradoxically)
- Impaired iron homeostasis
- Elevated oxidative stress responses
- Perturbed bacterial stress responses
- Growth defects under certain conditions

### Part 2: Experimental Design and Rationale

#### 2.1 Experimental Overview

**Research Question:**
How does queuosine deficiency affect the transcriptome response to nickel stress?

**Experimental Design:**
2 × 2 factorial design:
- **Factor 1 (Genotype):** WT vs. TGT knockout
- **Factor 2 (Treatment):** Control vs. Nickel stress
- **Replicates:** 3 biological replicates per group
- **Total samples:** 4 groups × 3 replicates = 12 samples

**Sample Groups:**
1. WT.none (wild-type control) - baseline
2. WT.Nickel (wild-type + 2 mM Ni²⁺)
3. TGT.none (tgt knockout control)
4. TGT.Nickel (tgt knockout + 2 mM Ni²⁺)

#### 2.2 Biological Hypotheses

**Hypothesis 1 (Nickel Stress Response in WT):**
- Expect upregulation of nickel-responsive genes
- Expect activation of metal homeostasis pathways
- Expect ROS defense mechanism activation
- Expected genes: nikA, nikB, nikC, sodA, sodB, catalase, etc.

**Hypothesis 2 (Q Deficiency Effect):**
- Expect constitutive upregulation of stress responses
- Expect impaired iron metabolism (Q-regulated protein translation)
- Expect elevated basal oxidative stress
- Expected genes: metal transporters, chaperones, antioxidants

**Hypothesis 3 (Interaction Effect):**
- Combined Q deficiency + nickel stress may show:
  - Additive effects (Q-deficiency response + nickel response)
  - Synergistic effects (nickel stress more severe in Q-deficient background)
  - Antagonistic effects (cells already maximally stressed)

#### 2.3 Practical Experimental Considerations

**Sample Selection:**
- 40 million reads per sample (adequate for bacterial genome)
- 75 bp paired-end reads (sufficient for mapping uniqueness)
- 3 biological replicates (standard for DGE studies)
- Not strand-specific (trade-off: cost vs. precision)

**Library Preparation Strategy:**
- rRNA depletion preferred over polyA selection (bacterial RNA lacks polyA)
- Ribo-Zero Magnetic Kit removes ~90% of ribosomal RNA
- Enables detection of non-coding RNAs if relevant
- Allows proper quantification of low-abundance mRNAs

**Quality Control Checkpoints:**
1. **Pre-sequencing:** RNA integrity (RIN > 7), 260/280 > 1.8
2. **Raw reads (FastQC):** Phred quality ≥ 30, no adapter contamination
3. **After mapping (Qualimap):** >90% mapped, <5% multi-mapped
4. **After counting (HTSeq-count):** Consistent counts across samples
5. **Normalized data (PCA):** Replicates cluster together

---

## SECTION II: RNA-SEQ DATA PROCESSING PIPELINE

### Part 3: Read Processing and Quality Metrics

#### 3.1 FastQC Quality Metrics

**What is FastQC?**
Tool analyzing quality of raw FASTQ sequencing files, assessing:

**Per-Base Quality Scores:**
- Phred quality score: Q = -10 × log₁₀(P_error)
- Q ≥ 30: Error probability ≤ 0.1% (99.9% accuracy)
- Q ≥ 20: Error probability ≤ 1% (99% accuracy)
- Quality typically declines toward 3' end of reads

**Interpretation from MultiQC:**
- Trend: Quality scores should remain ≥ 30 across read length
- Problem if: Quality drops sharply (poor library prep) or falls below 20
- Position effect: Common for read quality to decline at 3' end
- Red zones (warning): Quality ≤ 20 indicates problematic bases

**Per-Sequence Quality Distribution:**
- Most reads should cluster at high quality (Q > 30)
- Secondary peaks might indicate:
  - Contamination (mixed species)
  - Biased regions (systematic sequencing errors)
  - Remnants of low-quality sequences that passed filtering

**GC Content Distribution:**
- Bacterial genomes typically have 40-60% GC content
- Should see approximately normal distribution
- Deviation from expected: Suggests contamination or bias
- Two peaks: Possible contamination from different organisms

**Overrepresented Sequences:**
- Adapters (should be removed before sequencing)
- rRNA sequences (even with depletion, some remain)
- Highly repetitive genomic regions
- Contaminating species

**Key Indicators of Quality Issues:**
- Per-base quality < 20: Consider read trimming
- Adapter contamination > 1%: Poor library preparation
- Overrepresented sequences: Possible contamination or repeats
- GC distribution bimodal: Mixture of organisms

#### 3.2 STAR Alignment Algorithm

**What is STAR?**
Splicing-Aware Aligner (STAR) maps RNA-seq reads to reference genome
- Considers intron-exon junction information
- Even though bacteria lack introns, STAR handles edge cases efficiently

**Alignment Algorithm Steps:**

1. **Seed Searching Phase:**
   - Find maximal perfect matches (seeds) of length ≥ 19 bp
   - Scan entire read sequence for matching sequences in reference
   - Store all seed positions

2. **Seed Clustering and Extension:**
   - Cluster nearby seeds into potential alignments
   - Extend seeds ungapped to maximal match length
   - Generates multiple candidate alignments per read

3. **Scoring and Selection:**
   - Gapped Smith-Waterman alignment for gaps/mismatches
   - Score based on: match identity, gap penalties, mapping quality
   - Select best scoring alignment (or multiple if ties)
   - Assign mapping quality (MAPQ) score

4. **Output in SAM/BAM Format:**
   - SAM: Human-readable text (1 Gb per million reads)
   - BAM: Binary compressed (100 Mb per million reads)
   - Contains: Sequence, quality, alignment position, mapping quality

**Mapping Quality Metrics (from MultiQC/Qualimap):**

| Metric | Good Range | Interpretation |
|--------|-----------|-----------------|
| Mapping rate | >90% | Most reads successfully align |
| Uniquely mapped | >85% | Confidence in alignment location |
| Multi-mapped | <5% | Minimal ambiguity |
| Insert size mean | Expected ± 2 SD | Quality indicator for pairs |
| Unmapped reads | <5% | Possible contamination or quality issues |
| Chimeric reads | <1% | Rare fusions or structural issues |

**Interpretation for Queuosine Study:**
Expected excellent mapping due to:
- Well-characterized reference genome (*E. coli* K-12 MG1655)
- No splicing events (bacteria)
- High-quality sequencing
- Paired-end reads improve mapping confidence

#### 3.3 HTSeq-count: Read Quantification

**Counting Problem:**
Assigning reads to genes when ambiguity exists

**Ambiguity Sources:**
1. Read spans multiple features (exons of different genes)
2. Read partially overlaps a feature
3. Read in intergenic region (bacteria mostly lack introns)
4. Read multi-maps to multiple genomic locations

**Counting Modes:**

**Mode 1: UNION (most lenient)**
- Rule: Count if overlaps ANY exon (part or full)
- Advantage: Maximum sensitivity
- Disadvantage: Possible misassignment
- Use case: When maximum counts needed

**Mode 2: INTERSECTION-STRICT (most stringent)**
- Rule: Count ONLY if completely contained within one feature
- Advantage: High confidence assignments
- Disadvantage: Lost counts (sensitivity)
- Use case: When confidence most important

**Mode 3: INTERSECTION-NONEMPTY (default, used here)**
- Rule: Count if overlaps exon; allows small multi-feature overlaps
- Resolution: Requires feature to be longest overlap
- Advantage: Balanced sensitivity/specificity
- Use case: Standard RNA-seq analysis

**Standard Settings for This Study:**
- Mode: INTERSECTION-NONEMPTY
- Stranded: NO (libraries not strand-specific)
- Minimum overlap: 1 bp
- Ambiguous reads: NOT counted
- Multi-mapped reads: NOT counted (separate count category)
- Minimum mapping quality: Not filtered (HTSeq counts all)

**Output Data Structure:**
- Rows: genes (identified by feature ID, name, or locus tag)
- Columns: samples
- Values: Fragment counts (typically 0-100,000 per gene)
- Metadata: Gene annotations (description, protein ID, etc.)

#### 3.4 MultiQC Integration Report

**Purpose:**
Aggregate quality metrics from multiple tools and samples into single interactive report

**Tools Integrated:**
1. **FastQC:** Raw read quality across all samples
2. **STAR:** Mapping statistics per sample
3. **HTSeq-count:** Feature assignment summary
4. **Qualimap:** Post-mapping quality assessment
5. **Other tools:** Adapters, biases, contamination checks

**Key Visualizations:**

**Sample Summary Plots:**
- Heat map showing all samples' metrics
- Quickly identify outliers (unusual mapping, quality, etc.)
- Hierarchical clustering of samples by QC metrics

**Per-Sample Plots:**
- Individual sample quality timeline
- Library size (total reads per sample)
- Mapping rate progression through samples

**Red Flags in MultiQC:**
1. Outlier samples with poor mapping (highlighted)
2. Uneven library sizes (>2-fold difference concerning)
3. Per-sample quality issues (low overall or per-position)
4. Inconsistent GC content between samples
5. Batch effects visible in clustering

**Interpretation for Queuosine Study:**
Expected uniform quality across samples:
- Similar library sizes (40 million ± 10%)
- Consistent mapping rates (>90%)
- Normal GC distribution
- No obvious technical batch effects
- Proper replication structure visible in clustering

---

## SECTION III: DIFFERENTIAL EXPRESSION ANALYSIS

### Part 4: DESeq2 Statistical Framework

#### 4.1 Count Data Distribution and Model Selection

**Why Not Poisson?**

Poisson distribution assumes:
- Variance = Mean
- Independence between counts

RNA-seq reality:
- Variance >> Mean (overdispersion)
- Replicate-to-replicate variation substantial
- Biological variation (genes express differently in replicates)

Example with real data:
- Gene A: Mean = 100 counts, SD = 5 → Variance = 25 (Poisson would expect 100)
- Gene B: Mean = 50 counts, SD = 20 → Variance = 400 (Poisson would expect 50)

**Negative Binomial Solution:**

Negative binomial adds dispersion parameter (α):
- Variance = μ + α·μ²
- α → 0: Approaches Poisson (minimal overdispersion)
- α > 0: Captures biological variation
- α large: Gene shows high replicate-to-replicate variability

**DESeq2 Approach:**
- Estimate α empirically for each gene
- Use empirical Bayes shrinkage toward trend
- More conservative p-values (accounts for biological variation)

#### 4.2 Normalization: Size Factor Estimation

**Problem:**
Sequencing depth (library size) varies between samples

**Example:**
```
Sample 1: 40 million reads
Sample 2: 30 million reads (25% fewer)

Gene A raw counts: Sample 1 = 1000, Sample 2 = 800
Naive conclusion: Gene A upregulated in Sample 1
Correct conclusion: Similar expression (accounting for depth)
```

**DESeq2 Median-of-Ratios Method:**

Algorithm:
1. Calculate geometric mean for each gene (across all samples)
2. For each sample: divide gene count by geometric mean
3. Size factor = median of all ratios per sample

Intuition:
- If sample deeply sequenced: Size factor > 1
- If sample shallowly sequenced: Size factor < 1
- Normalized counts = Raw counts / Size factor

Example calculation:
```
Gene A: counts [1000, 800], geometric mean = 894.4
Gene B: counts [500, 400], geometric mean = 447.2

Sample 1 ratios: [1000/894.4=1.12, 500/447.2=1.12] → median = 1.12
Sample 2 ratios: [800/894.4=0.89, 400/447.2=0.89] → median = 0.89

Size factors: Sample 1 = 1.12, Sample 2 = 0.89
Normalized: Sample 1 gene A = 1000/1.12 = 893
            Sample 2 gene A = 800/0.89 = 899
Conclusion: Similar expression (justified!)
```

**Advantages:**
- Robust to highly-expressed genes (doesn't assume constant normalization)
- Robust to genes with zero counts
- No need for RNA spike-ins (internal normalization)
- Doesn't assume total RNA constant (valid for DGE where some genes change dramatically)

#### 4.3 Dispersion Estimation Process

**Step 1: Per-Gene Maximum Likelihood Estimation**

Estimate α individually:
- For each gene: Calculate likelihood of counts given mean and α
- Find α maximizing likelihood
- Result: Separate α for each gene
- Problem: Few replicates (n=3) → estimates very noisy

**Step 2: Trend Fitting**

Observation: Relationship between mean expression and dispersion
- Highly expressed genes: Lower dispersion (less relative noise)
- Lowly expressed genes: Higher dispersion (more relative noise)

Procedure:
- Plot estimated α vs. mean expression
- Fit smooth curve (loess regression) through data
- Trend describes "expected" dispersion at each expression level

**Step 3: Empirical Bayes Shrinkage**

Shrink per-gene estimates toward trend:
- Gene with unusual dispersion pulled toward expectation
- Amount of shrinkage: proportional to deviation and strength of trend
- Result: More stable estimates, especially for low-count genes

Example:
```
Gene A (mean=10): 
  Estimated α = 5.0 (very high, possibly noise)
  Trend α = 1.0 (expected for this mean)
  Shrunk α = 2.5 (compromise)

Gene B (mean=10):
  Estimated α = 0.95 (very close to trend)
  Trend α = 1.0 (expected)
  Shrunk α ≈ 0.95 (minimal shrinkage)
```

**Impact on Results:**
- Low-count genes: p-values more conservative (larger α → larger variance)
- High-count genes: Minimal shrinkage (estimates already reliable)
- Overall: FDR control improved (fewer false positives)

#### 4.4 GLM Fitting and Hypothesis Testing

**Generalized Linear Model for Count Data:**

DESeq2 uses negative binomial GLM with log link:

log(μᵢⱼ) = β₀ + β₁·genotypej + β₂·treatmentj + ...

Where:
- μᵢⱼ = expected count for gene i in sample j
- β₀ = intercept (baseline expression)
- β₁, β₂ = effect coefficients (log fold changes)

**Likelihood Estimation:**
- For each gene: Find coefficients maximizing likelihood
- Standard errors estimated from second derivative (Fisher information)
- Convergence algorithm (iteratively reweighted least squares)

**Wald Test for Significance:**

For each coefficient:
- Test H₀: β = 0 (no effect)
- Test statistic: z = β̂ / SE(β̂) ~ Normal(0,1)
- Two-tailed p-value: p = 2·P(|Z| > |z_obs|)

**Multiple Testing Correction:**

Challenge: Testing ~4,300 genes
- Probability of false positive in ≥1 test: 1 - (1-0.05)^4300 ≈ 100%
- Need stringent correction

Benjamini-Hochberg FDR Control:
1. Rank p-values from smallest to largest
2. For rank i: Adjusted p-value = p_i · (m / i)
3. Floor the adjusted p-values (enforce monotonicity)
4. Interpretation: p.adjust < 0.05 → FDR ≤ 5%

Effect:
- Less stringent than Bonferroni (9·10⁻⁶ vs 1.16·10⁻⁵)
- Maintains control of false discoveries
- Allows ~5% false positives among reported genes

#### 4.5 Log Fold Change and Effect Size

**Definition:**
LFC = log₂(Expression_B / Expression_A)

**Computation:**
- Fitted values from GLM (expected values)
- Ratio of group B mean to group A mean
- Log₂ transform for symmetry (LFC=+1 equivalent to LFC=-1)

**Interpretation:**
| LFC | Fold Change | Interpretation | Biological |
|-----|-------------|----------------|-----------|
| 0 | 1.0× | No change | Same expression |
| 1 | 2.0× | Upregulated | Twice as high |
| 2 | 4.0× | Strongly up | 4× expression |
| -1 | 0.5× | Downregulated | Half the level |
| -2 | 0.25× | Strongly down | 1/4 expression |

**LFC Shrinkage (Adaptive Shrinkage Estimation):**

Problem with raw LFC:
- Small sample sizes (n=3) → unreliable estimates
- Low-count genes: wide confidence intervals, extreme values
- Example: Gene with counts [0,1,2] vs [5,10,8]
  - Same total (3 counts) but first estimates more unreliably
  - Raw LFC might be extreme (±5) due to noise

Solution: Use empirical Bayes shrinkage
- Prior: Distribution of observed LFC values
- Bayesian update: Shrink unusual estimates
- Effect: More stable ranking, better interpretation
- Mechanism: Genes with uncertain LFC pulled toward zero

Result:
```
Gene with stable high expression:
  Raw LFC = 2.0 → Shrunk ≈ 2.0 (minimal shrinkage)

Gene with few reads:
  Raw LFC = 3.5 (extreme) → Shrunk ≈ 1.2 (substantial shrinkage)
```

**Impact on Results:**
- Better ranking genes by true effect size
- Improved sensitivity (less noise)
- Downstream: Improved biological interpretation

---

## SECTION IV: QUALITY CONTROL AND VALIDATION

### Part 5: Post-Analysis Quality Control

#### 5.1 Principal Component Analysis (PCA)

**Purpose:**
Reduce high-dimensional data (gene expression: thousands of genes) to 2D visualization while preserving sample relationships

**Mathematical Principle:**

PCA finds linear combinations of genes (Principal Components) that maximize variance:

PC1 = w₁₁·Gene₁ + w₁₂·Gene₂ + ... + w₁ₙ·Geneₙ (maximum variance)
PC2 = w₂₁·Gene₁ + w₂₂·Gene₂ + ... + w₂ₙ·Geneₙ (next maximum variance, orthogonal)

Where:
- wᵢⱼ = weights (loadings) for PC_i
- Orthogonality ensures independence between PCs

**Variance Explained:**
- PC1 typically 40-60% variance (largest effect, often genotype or treatment)
- PC2 typically 15-30% variance (secondary effect)
- Remaining variance in PC3+

**Interpretation of Sample Plots:**

**Expected Pattern (Good Quality):**
```
PC2
  ^
  |  TGT.none ·· TGT.Nickel
  |  ··              ··
  |  ·                 ·
  |                    
  +--·----·----·----·---> PC1
  |  ·  WT.none  ·
  |  ··              ··
  |     WT.Nickel ··
  |
```

Replicates cluster → Good reproducibility
Treatment separation → Treatment effects present
Genotype separation → Genotype effects present

**Problem Patterns:**

*Outlier Sample:*
- Single replicate distant from group
- Suggests: contamination, mislabeling, technical failure

*Poor Clustering:*
- Replicates spread across plot
- Suggests: High biological variation or batch effects

*No Separation by Treatment:*
- Nickel-treated samples overlay control
- Suggests: Treatment ineffective or conditions too similar

**Variance Explained Interpretation:**
- High variance on PC1 (>80%): One dominant factor (strong biological effect)
- Balanced variance (40-50% each): Multiple biological factors
- Very low variance explained by PC1 (<20%): Weak or absent biological effects

**Data Transformation (rlog):**
Important: PCA uses transformed data, not raw counts
- Raw count data: Heteroscedastic (variance increases with mean)
- rlog transformation: Regularized log-transformation
  - Counts stabilizes variance
  - Improves visualization
  - Prevents high-count genes from dominating PCA

#### 5.2 Individual Gene Examination

**Purpose:**
Validate specific genes of interest visually

**Gene Selection Strategy:**

1. **Knockout target gene (tgt):**
   - Should be expressed in WT
   - Should be absent/very low in TGT knockout
   - Validates experimental design success

2. **Expected responsive genes:**
   - Metal transporters (nikA, nikB, nikC)
   - Stress response genes (sodA, sodB)
   - Expected to show condition-specific patterns

3. **Housekeeping genes:**
   - Should be relatively constant across conditions
   - Validates that normalization didn't remove real differences

4. **Top DE genes:**
   - Genes with smallest p-value or largest LFC
   - Should show clear patterns across conditions

**Visualization Approach:**

Box plot or dot plot with replicates:
- X-axis: Experimental conditions
- Y-axis: Normalized count
- Points: Individual replicates
- Enables visual assessment of:
  - Effect magnitude
  - Replicate consistency
  - Outliers or unusual patterns

**tgt Gene Expected Pattern:**

| Condition | Expected tgt Expression | Biological Explanation |
|-----------|------------------------|----------------------|
| WT.none | Normal (baseline) | tgt gene functional, normal expression |
| WT.Nickel | Normal or slightly increased | Possible Q-biosynthesis upregulation under stress |
| TGT.none | ~Zero | tgt knockout; no functional transcripts |
| TGT.Nickel | ~Zero | tgt still absent (permanent knockout) |

If tgt shows different pattern → Questions validity of knockout or sample preparation

#### 5.3 Volcano Plot Analysis

**Purpose:**
Simultaneous visualization of effect size AND statistical significance

**Axes:**
- X-axis: log₂ Fold Change (effect size)
  - Left (-): Downregulated genes
  - Right (+): Upregulated genes
  - 0: No change
  
- Y-axis: -log₁₀(p.adjust) (significance)
  - Higher: More significant (p.adjust → 0, -log₁₀ → ∞)
  - Lower: Less significant (p.adjust → 1, -log₁₀ → 0)

**Threshold Lines:**
- Vertical lines: LFC = ±1 (2-fold change cutoff)
- Horizontal line: -log₁₀(0.05) ≈ 1.3 (FDR threshold)

**Expected Shape (Volcano):**
- Narrow, pointed peak: Few large-effect significant genes
- Wide, flat base: Many genes with small effect
- Apex at LFC = 0: Symmetric (unbiased effect)

**Problem Patterns:**

*Flat Distribution:*
- Genes scattered uniformly
- Suggests: Low statistical power or high noise

*Only One Direction:*
- All significant genes upregulated (or downregulated)
- Suggests: Systematic bias or contamination

*Non-symmetric:*
- More upregulated than downregulated (or vice versa)
- May be biological (some processes asymmetric)
- But if extreme: Check for batch effects

**Interpretation for Nickel Stress:**

Expected volcano plot:
- Upregulated cluster (right side): Metal response, oxidative stress genes
- Downregulated cluster (left side): Growth, proliferation genes (inhibited under stress)
- Symmetric: Balanced up/down response
- Top genes: Strong candidates for validation

---

## SECTION V: GENE ONTOLOGY ENRICHMENT

### Part 6: Functional Profiling with GO

#### 6.1 Gene Ontology Structure in Detail

**GO Hierarchical Organization:**

The three GO namespaces (independent hierarchies):

**1. Biological Process (BP) - 45,000+ terms**

Root: "biological_process"

Descendants include:
- "metabolic process" → "energy metabolism" → "glycolysis"
- "response to stimulus" → "response to chemical" → "response to nickel"
- "biological regulation" → "gene expression" → "transcription"

Example for queuosine study: "metal ion homeostasis"
- Path 1: metal ion homeostasis ← ion homeostasis ← homeostatic process ← biological regulation
- Path 2: metal ion homeostasis ← metal cluster assembly ← cofactor biosynthesis
- Result: Gene annotated to "metal ion homeostasis" is implicitly also annotated to all ancestors

**2. Molecular Function (MF) - 13,000+ terms**

Root: "molecular_function"

Examples:
- "binding" → "metal ion binding" → "nickel ion binding"
- "catalytic activity" → "transferase activity" → "kinase activity"
- "transporter activity" → "transmembrane transporter" → "active transmembrane transporter"

**3. Cellular Component (CC) - 3,000+ terms**

Root: "cellular_component"

Examples:
- "membrane" → "plasma membrane" → "outer membrane"
- "organelle" → "mitochondrion" (in eukaryotes)
- "ribosome" → "large ribosomal subunit" → "60S ribosomal subunit"

**Directed Acyclic Graph (DAG) Structure:**

Not a simple tree hierarchy:
- Multiple inheritance: term can have multiple parents
- Example: "iron-sulfur cluster assembly" has parents:
  - "metal cluster assembly"
  - "biosynthetic process"
  - "cellular process"
  All are valid ancestries, not mutually exclusive

This structure allows capturing biological knowledge precisely without forcing artificial hierarchies.

#### 6.2 Over-Representation Analysis Theory

**Basic Concept:**

Question: Are genes with property X more frequent in my gene set than background?

Example: Are metal homeostasis genes overrepresented in upregulated genes?

**Statistical Test: Hypergeometric Distribution**

Setup:
- N = total genes in genome (4,300 for E. coli)
- M = genes annotated to GO term (e.g., 45 genes for "metal ion homeostasis")
- n = genes in query set (upregulated genes, e.g., 150)
- k = overlap between query and term-annotated genes (e.g., 8)

**Null Hypothesis:**
Term genes are randomly distributed among query genes.

Expected k = n × M / N = 150 × 45 / 4300 ≈ 1.6 genes

**Observed k = 8 genes**

Much more than expected!

**Hypergeometric Test:**
P(X ≥ 8) = ?

Formula (probability of observing k or more):
P(X = k) = [C(M,k) × C(N-M, n-k)] / C(N, n)

Where C(n,k) = binomial coefficient = n! / (k!(n-k)!)

**Computation:**
Sum P(X=k) for k = 8 to min(n, M)

**Interpretation:**
- p-value = 0.001: Only 0.1% chance of 8+ overlap if random
- Conclusion: Term genes significantly overrepresented (8 >> 1.6 expected)
- Result: This GO term is "enriched" in upregulated genes

#### 6.3 Multiple Testing in Enrichment Analysis

**Challenge:**
Testing 10,000+ GO terms (every term in ontology)

**Probability of False Positives:**
- Without correction: ~500 false discoveries at α=0.05
- Makes it hard to identify true signals

**Solution: Benjamini-Hochberg FDR Correction**

Same procedure as gene-level testing:
1. Rank p-values (smallest to largest)
2. Adjust: p.adjust_i = p_i × (m / rank_i)
3. Monotonicity: p.adjust must be increasing
4. Result: "adjusted p-value" or "q-value"

**Interpretation:**
- p.adjust < 0.05 → FDR ≤ 5%
- On average, 5% of reported enriched terms are false positives
- Other 95% are true positives

**For enrichment analysis:**
Typical thresholds:
- pvalueCutoff = 0.05 (before correction, for speed)
- qvalueCutoff = 0.05 (after correction, stringent)
- Only terms passing both reported

#### 6.4 Gene Ratio vs. Background Ratio

**GeneRatio Definition:**
= (# query genes annotated to term) / (total # query genes)

Example: 
- Term: "metal ion homeostasis"
- Query set: 150 upregulated genes
- Genes in query + term: 8
- GeneRatio = 8/150 = 0.053 = 5.3%

Interpretation: 5.3% of upregulated genes are metal homeostasis genes

**BgRatio Definition:**
= (# background genes annotated to term) / (total # background genes)

Example:
- Background: All 4,300 genes in genome
- Genes in background + term: 45
- BgRatio = 45/4300 = 0.010 = 1%

Interpretation: 1% of all genes are metal homeostasis genes

**Enrichment Fold Change:**
= GeneRatio / BgRatio = 0.053 / 0.010 = 5.3×

Interpretation: Metal homeostasis genes enriched 5.3-fold in upregulated set

**Visualization (Dotplot X-axis):**
- X-axis: GeneRatio
- Terms with higher GeneRatio more "prominent" in gene set
- Also shows Count (size of dot)
- Larger dots: More genes in term (more robust finding)

---

## SECTION VI: EXPECTED FINDINGS FOR QUEUOSINE STUDY

### Part 7: Biological Predictions and Validation

#### 7.1 Predictions Based on Literature

**From Pollo-Oliveira et al. (2022):**

The original queuosine study found:

**In tgt mutant (Q-deficiency) without nickel:**
- 89 upregulated genes
- 65 downregulated genes
- Strong metal homeostasis dysregulation
- Elevated oxidative stress responses
- Impaired iron metabolism

**Top upregulated genes:**
- nikA, nikB, nikC (nickel transport - 10-30 fold upregulation)
- fecA, fecB (iron transport)
- iscS, iscU, iscA (Fe-S cluster assembly)
- sodA, sodB (superoxide dismutase)
- trx (thioredoxin)
- groEL, groES (chaperones)

**Top downregulated genes:**
- Some FUR-regulated genes (reduced iron availability changes FUR state)
- Certain ribosomal protein genes
- Some translation factors

**Expected GO Enrichments for TGT vs. WT:**

Biological Process (upregulated):
1. "metal ion homeostasis" - Strong
2. "response to oxidative stress" - Strong
3. "iron-sulfur cluster assembly" - Moderate
4. "metal ion transport" - Strong
5. "response to metal ion" - Strong
6. "antioxidant activity" - Moderate

Molecular Function (upregulated):
1. "metal ion binding" - Strong
2. "oxidoreductase activity" - Strong
3. "iron binding" - Moderate
4. "transporter activity" - Strong

Cellular Component (upregulated):
1. "cell membrane" - Strong (transporters locate here)
2. "iron-sulfur cluster" - Moderate
3. "cytoplasm" - Weak (general)

**In WT strain with nickel treatment:**
- ~200-300 upregulated genes
- ~150-250 downregulated genes
- Strong metal ion response
- Oxidative stress response activation
- Growth genes downregulated

**Top upregulated genes:**
- nikA, nikB, nikC (strong induction)
- sodA, sodB (strong induction)
- Various metal transporters
- Heat shock proteins
- Oxidative stress response genes

**Expected GO Enrichments for WT.Nickel vs WT.none:**

Biological Process (upregulated):
1. "response to metal ion stress" - Very strong
2. "response to oxidative stress" - Strong
3. "metal ion transport" - Strong
4. "metal ion homeostasis" - Strong
5. "protein folding" - Moderate (chaperone activation)

#### 7.2 Hypothesis Testing Framework

**Hypothesis 1: Q Regulates Metal Homeostasis**

Prediction:
- TGT.none vs WT.none should show metal homeostasis dysregulation
- nikA, nikB, nikC, fecA, fecB should be upregulated (LFC > 1)
- GO term "metal ion homeostasis" should be significantly enriched
- GO term "metal ion transport" should be significantly enriched

Validation:
- Examine individual gene expression plots
- Check that expected genes are in upregulated DE list
- Verify GO enrichment terms match predictions
- Compare fold changes to literature values

If supported: Q plays crucial role in metal homeostasis regulation

**Hypothesis 2: Q Deficiency Causes Oxidative Stress**

Prediction:
- TGT.none should show elevated oxidative stress genes
- sodA, sodB, catalase, trx should be upregulated
- GO term "oxidative stress response" should enrich
- ROS scavenging pathways should activate

Validation:
- Examine antioxidant genes in DE list
- Check enrichment of ROS-related GO terms
- Look for stress hormone pathway activation

If supported: Q maintains redox balance; deficiency causes ROS elevation

**Hypothesis 3: Metal Stress Response Operates Independently of Q**

Prediction:
- Nickel should induce metal transporters in BOTH WT and TGT
- Same metal response genes should activate regardless of Q status
- Similar GO terms enriched in WT.Nickel and TGT.Nickel comparisons
- Fold changes might differ (TGT baseline already elevated)

Validation:
- Compare upregulated genes between two nickel treatment comparisons
- Examine directionality (should both be upregulated in both)
- Check GO term overlap

If supported: Q-independent metal sensing still functional in Q-deficiency

#### 7.3 Validation Approaches

**1. Cross-Species Comparison**
- Are metal homeostasis genes conserved in other bacteria?
- Do similar genes upregulate in *Salmonella*, *Vibrio*, other enterobacteria?
- Suggests biological importance

**2. Protein-Level Validation**
- Proteomics: Are metal transport proteins elevated?
- Immunoblot: Confirm select proteins upregulated
- Native PAGE: Measure actual metal uptake/homeostasis
- ROS measurement: Direct measurement of superoxide/H₂O₂

**3. Phenotypic Validation**
- Growth curves: TGT slower than WT under nickel stress?
- Tolerance assays: Minimum inhibitory concentration (MIC) of nickel
- ROS sensitivity: Paraquat/H₂O₂ challenge tests
- Metal levels: ICP-MS measurement of cellular metal accumulation

**4. Mechanistic Follow-up**
- ChIP-seq: Which transcription factors regulate metal genes?
- Promoter analysis: Common regulatory motifs?
- 5' degradome-seq: tRNA fragmentation from Q deficiency?
- Ribosome profiling: Translation efficiency on metal protein genes?

**5. Literature Integration**
- PubMed searches for validated hits
- Protein database (UniProt) annotations
- Pathway databases (KEGG, String) functional connections
- Disease relevance (if human orthologs affected in disease)

---

## SECTION VII: COMMON PITFALLS AND TROUBLESHOOTING

### Part 8: Critical Interpretation Issues

#### 8.1 Statistical Fallacies to Avoid

**Fallacy 1: Statistical Significance = Biological Significance**

Problem:
- With large sample sizes and large n, ANY difference becomes "significant"
- p < 0.05 means NOT due to chance, not that it matters biologically

Example:
- Gene upregulated 1.01× fold: LFC = 0.01, p = 0.001
- Statistically significant (p < 0.05)
- Biologically trivial (only 1% change)

Solution:
- Always check LFC magnitude alongside p-value
- Set practical cutoff (LFC > 1 = 2× fold change minimum)
- Combine statistical + effect size criteria

**Fallacy 2: Absence of Evidence = Evidence of Absence**

Problem:
- Non-significant p-value (p > 0.05) doesn't mean no effect
- Could be: True null effect OR insufficient power

Example:
- Gene showed 1.5× fold change but p = 0.08
- Can't conclude "no effect" (power might be low)
- Need confidence interval

Solution:
- Examine confidence intervals, not just p-values
- Consider statistical power (sample size, effect size, variance)
- Report results with uncertainty quantified

**Fallacy 3: Correcting Without Understanding**

Problem:
- FDR correction (BH method) controls false discoveries
- But if apply correction TOO stringently, lose real findings

Example:
- Report only p.adjust < 0.001 instead of p.adjust < 0.05
- More conservative, BUT lose valid discoveries
- FDR 0.05 already controls 95% true positive rate

Solution:
- Understand what correction accomplishes
- Use recommended thresholds (p.adjust < 0.05 standard)
- When uncertain, run sensitivity analysis (try multiple thresholds)

#### 8.2 Data Quality Issues

**Issue: Contamination**

Red Flags:
- Mapping rate < 80%
- GC content bimodal
- Overrepresented sequences from non-target organisms
- Unusually high multi-mapped reads

Investigation:
- BLAST unmapped reads → identify contaminants
- Check sample history (how was RNA extracted?)
- Examine picture of culture plates
- Test for bacteria/phage contamination

Resolution:
- If minor (<5%): Consider analyzing as-is (mapping usually filters some)
- If major (>20%): Consider sample unusable

**Issue: Library Preparation Artifacts**

Red Flags:
- Adapters present in overrepresented sequences
- Unusual insert size distribution
- Reads with low complexity (repeated bases)

Investigation:
- Check trimming parameters (were adapters removed?)
- Verify library construction protocol was followed
- Examine positive controls

Resolution:
- Adapter trimming tools (Trimmomatic, Cutadapt)
- Quality trimming (QC step before mapping)
- Remove problematic samples if severe

**Issue: Batch Effects**

Red Flags:
- PCA: Samples cluster by preparation date, not condition
- Unexpected variance in library size
- Systematic bias in one group of samples

Investigation:
- Examine experimental timeline
- Check if batches correlate with biology or technique
- Statistical test: Estimate batch variance proportion

Resolution:
- Model batch in design formula: design = ~ batch + condition
- Batch correction methods (ComBat, RUVSeq)
- If confounded: May reduce power (can't separate batch from biology)

#### 8.3 Analysis Choices and Sensitivity

**Choice 1: Filtering Threshold**

Parameter: Minimum rowSum for filtering (commonly >10)

Sensitivity:
- rowSum > 5: Retain more genes, fewer reads per gene per sample, less power
- rowSum > 10: Standard choice, balance sensitivity/specificity
- rowSum > 20: More stringent, more confidence but fewer genes

Test impact:
- Run with different thresholds (5, 10, 20)
- Compare number of DE genes
- Check if top genes consistent
- If robust: Findings insensitive to threshold

**Choice 2: LFC Cutoff**

Parameter: Minimum |LFC| to consider "biologically meaningful"

Common choices:
- LFC > 0.5: More inclusive, may capture subtle changes
- LFC > 1.0: Standard, 2-fold change
- LFC > 2.0: Conservative, strong effects only

Test impact:
- Run enrichment with different cutoffs
- Top GO terms should be similar
- If different → consider what threshold captures your biology

**Choice 3: p-value Threshold**

Parameter: Maximum adjusted p-value

Standard: p.adjust < 0.05 (5% FDR)

Alternatives:
- 0.01: More stringent (1% FDR)
- 0.10: More lenient (10% FDR)

Test impact:
- Stricter threshold: fewer genes (higher confidence)
- Looser threshold: more genes (higher sensitivity)
- Check consistency of findings across thresholds

**Choice 4: Gene Set Definition for GO**

Parameter: Which genes to include in ORA

Options:
1. All DE genes (both up and down)
2. Upregulated only (LFC > 1)
3. Downregulated only (LFC < -1)

Impact:
- All DE: Detects any process affected
- Upregulated: Detects activation of processes
- Downregulated: Detects suppression of processes
- Should analyze each separately for complete picture

#### 8.4 Over-Interpretation Risks

**Risk 1: Circular Reasoning**

Fallacy:
"Gene is involved in metal homeostasis because GO enrichment for metal homeostasis"

Reality:
- GO terms assigned by manual annotation (biased toward studied pathways)
- Enrichment detects known annotations, not novel functions
- Gene may have other unknown functions

Proper interpretation:
- GO enrichment supports: "Metal homeostasis gene overrepresented in DE set"
- Implies: "Condition affects metal homeostasis pathway"
- Doesn't prove: "Gene's function is metal homeostasis" (circular)
- Should validate with mechanism studies

**Risk 2: Statistical Overfitting**

Fallacy:
"p < 0.05 proves hypothesis"

Reality:
- Multiple testing: Even random data produces p < 0.05 somewhere
- If test 100 terms: ~5 will be significant by chance
- Multiple hypotheses tested → multiple comparisons problem

Proper practice:
- Pre-register hypotheses (avoids p-hacking)
- Or explicitly correct for multiple comparison (use FDR)
- Validate findings in independent dataset if possible

**Risk 3: Generalization Beyond Study**

Fallacy:
"Results in rich medium (LB) apply to all conditions"

Reality:
- This study: Specific strain (K-12), specific medium (LB), specific temperature
- Results may not generalize to:
  - Natural environmental conditions
  - Other growth media
  - Other strains (genetic background matters)
  - Different stresses

Proper interpretation:
- Claim: "In LB medium at 37°C, Q deficiency upregulates metal homeostasis"
- Not: "Q always regulates metal homeostasis"

---

## CONCLUSION

This comprehensive guide covers:

1. **Biological Context:** Queuosine's role in translation and metal homeostasis
2. **Experimental Design:** 2×2 factorial design rationale and sample preparation
3. **Data Processing:** Quality control, alignment, counting methodologies
4. **Statistical Analysis:** DESeq2 framework from data distribution to hypothesis testing
5. **Enrichment Analysis:** Gene Ontology structure and over-representation statistics
6. **Quality Assessment:** PCA, individual gene plots, volcano plots for validation
7. **Biological Interpretation:** Expected findings, validation approaches, pitfalls

The integration of high-throughput sequencing, rigorous statistical analysis, and biological interpretation enables discovery of gene expression changes underlying cellular responses to genetic and environmental perturbations.

For the queuosine study specifically:

**Key Expected Findings:**
- Q deficiency dysregulates metal homeostasis (constitutive upregulation of metal transport)
- Q-deficient cells show elevated oxidative stress responses
- Nickel stress activates metal response genes in both WT and TGT
- Combined Q-deficiency + nickel stress shows complex interaction effects

**Validation Strategy:**
- Confirm key genes using qRT-PCR
- Measure protein levels (proteomics) for metal transport proteins
- Assess phenotypes (growth, tolerance, ROS levels)
- Investigate molecular mechanisms (promoter analysis, protein interactions)
- Validate GO enrichments through literature cross-reference

**Broader Significance:**
Understanding how tRNA modifications regulate cellular responses has implications for:
- Bacterial physiology and adaptation
- Human diseases with tRNA modification defects
- Therapeutic targets in pathogenic bacteria
- Translation quality control mechanisms
