# CCS - Cohort Congress System

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%203.6.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.7.3-brightgreen.svg)](https://github.com/yourusername/ccs)

## Overview

**CCS (Cohort Congress System)** is a computational framework for personalized pan-cancer genomic classification. It integrates ensemble machine learning, dimensionality reduction, and clustering algorithms to identify cancer subtypes from multi-omics data with high reproducibility and interpretability.

## Key Features

- **Ensemble Learning**: Combines multiple XGBoost models with gene set-based classifiers (GSClassifier)
- **Multi-omics Integration**: Handles genomic, transcriptomic, and other high-dimensional biological data
- **Robust Subtyping**: Employs consensus clustering with DBSCAN and NbClust for stable subtype identification
- **Parallel Processing**: Supports both ensemble and discrete parallel strategies for computational efficiency
- **Reproducibility**: Built-in seed management and model persistence for reproducible research
- **Visualization**: Nature-quality publication-ready plots for batch effects, dimensionality reduction, and subtype distributions

## Installation

### Prerequisites

- R ≥ 3.6.0
- Required dependencies (automatically installed):
  - `luckyBase` (core dependency)
  - `xgboost`, `caret`, `GSClassifier`
  - Visualization: `ggplot2`, `ComplexHeatmap`, `patchwork`
  - Dimensionality reduction: `umap`, `Rtsne`, `uwot`
  - Deep learning: `torch`, `luz`

### Install from Source

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install CCS from source
devtools::install_github("yourusername/ccs")
```

### Development Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/ccs.git
cd ccs

# Build and check the package
R CMD build .
R CMD check ccs_0.7.3.tar.gz

# Or use devtools for interactive development
Rscript -e "devtools::load_all()"
```

## Quick Start

### Basic Workflow

```r
library(CCS)

# Prepare your data
# data: a list of expression matrices (samples × genes)
# geneSet: a list of gene sets for classification
# geneAnnotation: gene annotation data frame

# Run CCS analysis
result <- ccs(
  data = your_expression_data,
  geneSet = your_gene_sets,
  geneAnnotation = your_gene_annotation,
  method = 'GSClassifier',
  geneid = "ensembl",
  parallel.method = 'ensemble',
  params = list(
    nfold = 5,
    nrounds = 100,
    nthread = 2,
    eta = 0.5,
    max_depth = 10,
    subsample = 0.7
  ),
  seed = 489,
  model.dir = './ccs_output',
  verbose = TRUE,
  numCores = 4
)

# Access results
subtypes <- result@Data$CCS
probabilities <- result@Data$Probability
models <- result@Model
```

### Visualization

```r
# Plot batch effects (Nature-quality style)
plotBatchEffect(
  data = your_data,
  batch = batch_labels,
  method = "umap"
)

# Visualize subtype distributions
# (Additional plotting functions available in the package)
```

## Project Structure

```
ccs/
├── R/                      # Core package code
│   ├── ccs.R              # Main CCS class and function
│   ├── ccs_model.R        # Model training functions
│   ├── ccs_embedding.R    # Dimensionality reduction
│   ├── hyperTuningGS.R    # Hyperparameter tuning
│   ├── importance.R       # Feature importance analysis
│   ├── clustering.R       # Clustering algorithms
│   ├── normalize.R        # Data normalization
│   ├── phenotype.R        # Phenotype processing
│   └── plotBatchEffect.R  # Visualization utilities
├── man/                   # Documentation (auto-generated)
├── test/                  # Test scripts and scenarios
│   ├── 01.test ccs.R
│   ├── 03.test.Classes.R
│   └── test.classifier_performance.R
├── DESCRIPTION           # Package metadata
├── NAMESPACE            # Exported functions
└── README.md           # This file
```

## Documentation

### Building Documentation

```r
# Generate Rd files from roxygen comments
devtools::document()
```

### Accessing Help

```r
# View main function documentation
?ccs

# View class documentation
?`CCS-class`

# Browse package vignettes (if available)
browseVignettes("CCS")
```

## Testing

Run scenario tests to validate functionality:

```bash
# Test core CCS functionality
Rscript test/01.test\ ccs.R

# Test S4 class definitions
Rscript test/03.test.Classes.R

# Test classifier performance
Rscript test/test.classifier_performance.R
```

## Development Guidelines

### Coding Style

- **Indentation**: 2 spaces
- **Assignment**: Use `<-` (not `=`)
- **Naming**:
  - Functions/variables: `snake_case`
  - S4 classes/methods: `CamelCase`
- **Documentation**: Complete roxygen headers (`@description`, `@param`, `@return`, `@examples`)

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Pull Request Checklist

- [ ] Code passes `R CMD check` without errors
- [ ] All relevant test scripts run successfully
- [ ] Documentation updated (roxygen comments + README if needed)
- [ ] CHANGELOG.md updated with changes
- [ ] Screenshots/metrics included for plotting/statistics changes

## Algorithm Overview

### Workflow

1. **Data Preprocessing**: Normalization and quality control
2. **Sub-model Training**: Parallel training of gene set-based classifiers
3. **Probability Aggregation**: Ensemble predictions across sub-models
4. **Dimensionality Reduction**: UMAP/t-SNE/PCA for visualization
5. **Consensus Clustering**: DBSCAN + NbClust for robust subtype identification
6. **Validation**: Cross-validation and performance metrics

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `nfold` | Cross-validation folds | 5 |
| `nrounds` | XGBoost boosting rounds | 100 |
| `eta` | Learning rate | 0.5 |
| `max_depth` | Maximum tree depth | 10 |
| `subsample` | Subsample ratio | 0.7 |
| `numCores` | Parallel cores | 4 |

## Citation

If you use CCS in your research, please cite:

```bibtex
@software{ccs2024,
  author = {Huang, Weibin},
  title = {CCS: Cohort Congress System for Pan-cancer Genomic Classification},
  year = {2024},
  version = {0.7.3},
  url = {https://github.com/yourusername/ccs}
}
```

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Author

**Weibin Huang**
📧 hwb2012@qq.com
🔗 [ORCID: 0000-0003-0159-435X](https://orcid.org/0000-0003-0159-435X)

## Acknowledgments

- Built on top of [xgboost](https://xgboost.readthedocs.io/), [GSClassifier](https://github.com/yourusername/GSClassifier), and [luckyBase](https://github.com/yourusername/luckyBase)
- Inspired by ensemble learning and consensus clustering methodologies in cancer genomics

## References

1. [NbClust Package](https://cran.r-project.org/web/packages/NbClust/index.html) - Determining optimal number of clusters
2. [DBSCAN Tutorial](http://www.sthda.com/english/wiki/wiki.php?id_contents=7940) - Density-based clustering
3. [Seurat Documentation](https://satijalab.org/seurat/articles/get_started_v5_new) - Single-cell analysis workflows

## Support

- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/yourusername/ccs/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/yourusername/ccs/discussions)
- 📖 **Documentation**: [Wiki](https://github.com/yourusername/ccs/wiki)

---

**Status**: Active Development | **Version**: 0.7.3 | **Last Updated**: 2026-03-07
