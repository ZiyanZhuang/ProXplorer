Sure, here's an open-source README for your protein analysis tool, under the MIT License, along with a suggested name:

# ProXplorer: Protein Sequence Analysis Toolbox

## Table of Contents

  * [About ProXplorer](https://www.google.com/search?q=%23about-proxplorer)
  * [Features](https://www.google.com/search?q=%23features)
  * [Getting Started](https://www.google.com/search?q=%23getting-started)
      * [Prerequisites](https://www.google.com/search?q=%23prerequisites)
      * [Installation](https://www.google.com/search?q=%23installation)
      * [Running the Application](https://www.google.com/search?q=%23running-the-application)
  * [Usage](https://www.google.com/search?q=%23usage)
  * [Code Overview](https://www.google.com/search?q=%23code-overview)
  * [Contributing](https://www.google.com/search?q=%23contributing)
  * [License](https://www.google.com/search?q=%23license)
  * [Contact](https://www.google.com/search?q=%23contact)
  * [Acknowledgments](https://www.google.com/search?q=%23acknowledgments)

## About ProXplorer

ProXplorer is an open-source MATLAB-based toolbox designed for comprehensive protein sequence analysis. It provides an intuitive and interactive workflow for researchers and students to perform common bioinformatics tasks, including sequence filtering, similarity searches, motif discovery, and phylogenetic tree construction. The tool aims to simplify complex computational biology analyses through a user-friendly graphical interface.

## Features

  * **FASTA File Loading & Summary**: Load protein sequences from standard FASTA files and get an immediate summary of the dataset, including total sequence count, min, max, mean, and median lengths.
  * **Sequence Length Filtering**: Filter sequences based on specified minimum and maximum length criteria and save the filtered set to a new FASTA file.
  * **Similarity Search**: Find proteins similar to a query sequence within your dataset using Smith-Waterman local alignment. Results are ranked by percent identity and alignment score, with the option to save the top N hits and visualize the similarity.
  * **Motif Discovery**: Search for specific protein motif patterns (e.g., PROSITE-like syntax like `P-X(2)-G`) within your sequences.
  * **Phylogenetic Tree Construction**: Build phylogenetic trees using selected sequences. The workflow includes:
      * Multiple Sequence Alignment (MSA) using a ClustalW-like algorithm (via `multialign`).
      * Pairwise evolutionary distance calculation using the JTT (Jones-Taylor-Thornton) model.
      * Neighbor-Joining (NJ) tree construction.
      * Bootstrap analysis to assess node confidence, with support values displayed on the tree branches.
  * **Interactive Visualizations**: Generate histograms for length distributions, and interactive scatter plots for similarity search results (identity vs. length, colored by alignment score).
  * **Bilingual User Interface**: Most user prompts and messages are provided in both English and Chinese for broader accessibility.

## Getting Started

### Prerequisites

  * MATLAB (R2020a or later recommended)
  * Bioinformatics Toolbox™ (required for `fastaread`, `fastawrite`, `swalign`, `multialign`, `seqpdist`, `seqneighjoin`, `seqbootstrp`, `getbybranches`, `plot(phytree)`)

### Installation

1.  **Clone or Download**:
    Clone this repository to your local machine:

    ```bash
    git clone https://github.com/yourusername/ProXplorer.git
    ```

    Or download the ZIP file and extract it.

2.  **Add to MATLAB Path**:
    Open MATLAB, navigate to the directory where you cloned/extracted ProXplorer.
    Add the folder and its subfolders to your MATLAB path:

    ```matlab
    addpath(genpath(pwd));
    savepath; % Optional: to make the changes permanent
    ```

### Running the Application

Once installed, simply run the `main_protein_analysis` function from the MATLAB command window:

```matlab
main_protein_analysis
```

This will launch the main user interface guiding you through the analysis workflow.

## Usage

The application will guide you through the following steps:

1.  **Load FASTA File**: Select your primary FASTA file containing protein sequences.
2.  **Initial Summary**: View a summary of your loaded dataset and its length distribution.
3.  **Optional Filtering**: Choose whether to filter sequences by length. If yes, provide min/max lengths and an output filename for the filtered sequences.
4.  **Main Analysis Menu**: Select from the following options:
      * **Similarity Search**: Provide a query FASTA file, specify the number of top hits to save, and optionally visualize the results.
      * **Motif Search**: Enter a motif pattern (e.g., `P-X(2)-G`) to search for within your sequences.
      * **Build Phylogenetic Tree**: Select sequences from your dataset and specify bootstrap replicates to construct and visualize a phylogenetic tree.
      * **Exit**: Terminate the program.

## Code Overview

  * `main_protein_analysis.m`: The main entry point of the application, orchestrating the overall workflow and user interactions.
  * `filter_and_save_proteins.m`: Function to filter protein sequences by length and save them to a FASTA file.
  * `find_similar_proteins.m`: Implements Smith-Waterman local alignment to find similar proteins, calculates percent identity, and saves top hits.
  * `find_sequences_by_motif.m`: Searches for user-defined protein motifs using regular expressions.
  * `build_phylogenetic_tree.m`: Handles the entire phylogenetic tree construction process, including MSA, distance calculation, tree building, and bootstrap analysis.
  * `visualize_similarity.m`: Generates an interactive scatter plot to visualize similarity search results.

## Contributing

Contributions are highly welcome\! If you have suggestions for improvements, new features, or bug fixes, please feel free to:

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature/YourFeature` or `bugfix/YourBugFix`).
3.  Make your changes.
4.  Commit your changes (`git commit -m 'Add new feature'`).
5.  Push to the branch (`git push origin feature/YourFeature`).
6.  Open a Pull Request.

Please ensure your code adheres to good MATLAB coding practices and includes appropriate comments.

## License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.

## Contact

For any questions or inquiries, please open an issue on the GitHub repository.

## Acknowledgments

  * The developers of MATLAB and the Bioinformatics Toolbox for providing powerful tools for computational biology.
  * The open-source community for inspiration and collaboration.
