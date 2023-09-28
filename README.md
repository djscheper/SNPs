# SNP Annotation Assignment

This script was developed by Dennis Scheper (373689)

---

### Introduction
The main goal of this script is to introduce a single nucleotide morphism (SNP) on a given position in a gene sequence specified by the user in the form of a FASTA file. The effect of the introduced SNP will be measured in conservation against a multiple sequence alignment (MSA), which includes protein sequences closely aligned with the given gene sequence. The script calculates conservation in percentage and gives a final verdict on whether an SNP is deleterious (<90% conservation) or has no effect (≥90% conservation).

A gene sequence and the accompanying MSA present in this repository can be used as an example for execution (`sequence.FASTA` and `msa.clustal`). These files are about the TP53 gene. The MSA contains ten closely aligned protein sequences from multiple species. [[1]][gene] [[2]][msa] However, users are free to use their own files.

### Deploying
To use this script, please follow the steps stated below.

**Step 1: Acquiring Files**

Either [clone][clone] or [download][download] the source files.

**Step 2: Installing Python**

The script was developed in the language Python (version 3.7.9). Please follow the instructions on how to install Python [here][Python].

**Step 3: Virtual environment**
Setting up a virtual environment is not required, although it is strongly advised. Please execute the following command to set up a virtual environment:

```bash
$ python3 -m venv venv
``

After that, virtual environment activation may differ depending on your operation system. 

Linux:
```{bash}
$ activate venv/Scripts/activate
```

Windows:
```{bash}
your\directories> venv\Scripts\activate
```

**Step 4: Installing Necessary Packages**

In addition, a necessary external package needs to be installed; [BioPython][biopython] (version 1.79) is used to read the FASTA and MSA files. After that, execute the following command on the command line:

```bash
$ pip install -r requirements.txt
```

Now, you are set to use this script.

### Usage
As stated before, the user is free to use their files. Therefore, I included one command for direct execution (copy-paste) and a more general description.

#### General description
```bash
$ python3 snp.py -n {A, C, G, T} -p <POSITION> -m <MSA.CLUSTAL FILE> -s <GENE_SEQUENCE.FASTA FILE>
```

Where `-n` represents the nucleotide—A, C, G or T—that forms the SNP, `-p` the position of the replacement, `-m` the MSA in the form of a Clustal file and `-s` the gene sequence as a FASTA file.

#### Direct execution
```bash
$ python3 snp.py -n A -p 1 -m msa.clustal -s sequence.FASTA
```

### Contact
If any issue or question remains, don't hesitate to get in touch with us at [d.j.scheper@st.hanze.nl](mailto:d.j.scheper@st.hanze.nl)

[msa]: https://www.ncbi.nlm.nih.gov/homologene/41131
[gene]: https://www.ncbi.nlm.nih.gov/nuccore/NM_016399.3
[clone]: https://djscheper@bitbucket.org/djscheper/snp_opdracht_bin3.git
[download]: https://bitbucket.org/djscheper/snp_opdracht_bin3/src/master/
[python]: https://www.python.org/
[biopython]: https://biopython.org/