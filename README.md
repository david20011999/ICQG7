# Exploring Imprinting Phenomena: AlphaSimR's New Functionalities

#### Abstract:
This study simulates a breeding scheme incorporating genomic imprinting as an example of the new AlphaSimR extension to simulate genomic imprinting. The simulation involves creating a founder population, establishing base populations with defined phenotypes, and performing multiple generations of selection and crossbreeding. The genetic values (GV) of individuals are calculated considering additive effects, dominance, and imprinting, and the performance of breeding strategies is compared over several generations. 

#### Introduction:
Genomic imprinting is an epigenetic phenomenon where the expression of a gene depends on its parental origin. This study aims to simulate a breeding scheme that accounts for genomic imprinting and compare its effects on genetic improvement with a scheme that considers only additive genetic effects.

#### Materials and Methods:

**1. Founder Population:**
- **Number of Individuals:** 1000
- **Number of Chromosomes:** 2
- **Segregating Sites:** 1000
- **Number of QTLs (Quantitative Trait Loci):** 300
- **Number of SNPs (Single Nucleotide Polymorphisms):** 700 (1000 - 300)

**2. Genetic Parameters:**
- **Mean and Variance for Additive Genetic Effects (Total Population):** 
  - Mean: (0, 0, 0)
  - Variance: (10, 10, 10)
- **Correlation Matrix for Additive Genetic Effects:**
  \[
  \begin{pmatrix}
  1.00 & 0.80 & 0.40 \\\
  0.80 & 1.00 & 0.40 \\\
  0.40 & 0.40 & 1.00
  \end{pmatrix}
  \]
- **Mean and Variance for Imprinting Effects:**
  - Mean: (0, 0, 0)
  - Variance: (9, 9, 9)
- **Correlation Matrix for Imprinting Effects:**
  \[
  \begin{pmatrix}
  1.00 & 0.80 & 0.90 \\\
  0.80 & 1.00 & 0.90 \\\
  0.90 & 0.90 & 1.00
  \end{pmatrix}
  \]

**3. Simulation Procedure:**

**3.1 Founder Genomes:**
- Founder genomes are generated using the `quickHaplo` function.

**3.2 Simulation Parameters:**
- Parameters are defined using `SimParam$new` and `SP$setSexes("yes_sys")`.

**3.3 Trait and SNP Chip Definition:**
- Traits with additive and imprinting effects are added using `SP$addTraitAI`.
- A SNP chip is added using `SP$addSnpChip`.

**3.4 Base Population:**
- A base population is created from the founder genomes using `newPop`.
- Phenotypes are assigned using `setPheno` with environmental variances (10, 10, 10).

**3.5 Selection and Breeding:**
- Two types of populations are created: paternal and maternal.
- Selection is performed for 50 generations based on phenotypic values.

**3.6 Crossbreeding:**
- Crossbreeding between selected individuals from paternal and maternal lines is performed for 10 generations.

**4. Data Collection and Analysis:**
- Genetic values are calculated considering additive, dominance, and imprinting effects.
- Performance is compared between breeding schemes with and without imprinting.
- Data is aggregated and analyzed using `ggplot2`.

#### Results:
The results show the comparison of genetic values across generations for two breeding strategies:
1. **Additive Only (A)**
2. **Additive + Imprinting (AI)**

**Figure:**
A plot illustrating the mean genetic values over generations with confidence intervals for each strategy.

#### Discussion:
The simulation demonstrates the impact of genomic imprinting on genetic improvement. The results highlight differences in the genetic values achieved under the two breeding strategies, emphasizing the importance of considering epigenetic factors in breeding programs.

#### Conclusion:
Incorporating genomic imprinting in breeding simulations provides a more comprehensive understanding of genetic improvement. This study underscores the necessity of considering both additive and imprinting effects for optimizing breeding schemes.
