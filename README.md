# <img src="https://raw.githubusercontent.com/vejnar/FONtools/main/img/logo.svg" alt="FONtools" width="45%" />

FONtools combines:
* The new FON (Feature Object Notation) format to store genomic annotations based on [JSON](https://en.wikipedia.org/wiki/JSON)
* Command-line tools to work with FON files: import (from [Ensembl](http://ensembl.org), [Ensembl Genomes](https://ensemblgenomes.org) or any GFF3), export, merge, mask sequence, filter annotations and more.
* A [Python](https://www.python.org) library to work with FON files

With FONtools:
* Using genomic annotations is easy and program-friendly. No need to write parser to extract specific features from annotations, use it directly from the FON files.
* Using genomic annotations is standardized. FON is based on JSON. Every standard library contains a JSON parser.
* Importing annotations from external source ([Ensembl](http://ensembl.org)) is convenient and included.
* Maintaining your repository of genomic annotations, sequence and indices for read-mapping software is built-in.

## Why the FON format?

FON is a program-friendly and extensible format to store genomic annotations. Since FON is using the JSON format for storing data, and JSON stands for JavaScript Object Notation, we named our format FON for Feature Object Notation.

Genomics annotations are mainly stored in files using the [BED](https://www.genome.ucsc.edu/FAQ/FAQformat.html#format1) or [GFF](https://en.wikipedia.org/wiki/General_feature_format), specifically [GFF3](http://gmod.org/wiki/GFF3), formats.

| Format  | Base format | Parse                   | Hierarchical | Extensible | Coordinates |
| ------- | ----------- | ----------------------- | ------------ | ---------- | ----------- |
| BED     | tab         | Simple                  | No           | Limited    | 0-based     |
| GFF     | tab         | Complex                 | Yes          | Yes        | 1-based     |
| **FON** | JSON        | Existing JSON libraries | Possible     | Yes        | 0-based     |

While BED is simple to parse, it was not designed to store hierarchical annotations, such as exons on a transcript. Instead BED12 "sub-splits" columns using commas instead of tabulations to store exons coordinates of transcripts. Alternatively, GFF allows hierarchical annotations, but is difficult to parse. GFF translates such structures into multiple records linked with a common ID. This approach is generic and describes annotations as a graph, thus requiring more complex code to parse it.

To overcome these limitations, FON format enables simple parsing and hierarchical annotation storage by capitalizing on the strengths of the JSON format:
* Dictionaries and lists are used to store structured annotations (limited to trees compared to graph with GFF), such as list of exons.
* Types are implicit: Quoted names are imported as strings while coordinates are imported as numbers. Describing the type of each column or attribute is no longer necessary. More properties can be added to any annotations without extending the format.

**FON isn't intended to replace GFF to share genomic annotations, but rather to simplify, ease and streamline the use of annotations within programs and pipelines.**

## FON1

FON1 is the first version of the Feature Object Notation format. It stores features in a list, and each feature is a dictionary with a set of defined keys. New keys for each feature can be freely added or removed, none of them are required. Programs using specific key(s) should provide the option to select by their name which key(s) to use (for example the `--key` option of `fon_mask_fasta`).

Example for one zebrafish transcript (one feature):
```json
{
    "fon_version": 1,
    "features": [
        {
            "transcript_stable_id": "ENSDART00000171909",
            "gene_stable_id": "ENSDARG00000099339",
            "gene_name": "pacsin3",
            "protein_stable_id": "ENSDARP00000138886",
            "chrom": "7",
            "strand": "+",
            "transcript_version": "2",
            "gene_version": "4",
            "transcript_biotype": "protein_coding",
            "gene_biotype": "protein_coding",
            "exons": [[54260003, 54260129], [54263505, 54263662]],
            "exons_on_transcript": [[0, 126], [126, 283]],
            "cds_exons": [[54260075, 54260129], [54263505, 54263662]],
            "cds_exons_on_transcript": [[72, 126], [126, 283]],
            "cds_exons_frame": [0, 0],
            "cds_exons_frame_on_transcript": [0, 0],
            "utr5_exons": [[54260003, 54260075]],
            "utr5_exons_on_transcript": [[0, 72]],
            "utr3_exons": [],
            "utr3_exons_on_transcript": [],
            "seq": "TTGGTCTCGCGTCTTGTTCTTCACAGTTTGACGACAGCCGCCATCATTCCGTGCTGCAAGGGCGACCCCAAAATGTCTTCCAACGGTGATCTGCAGGACGTTGGGAGTTGGGACAGCTTCTGGGAGCCTGGAAACTACAAGAGGACGGTTAAGCGCATTGACGACGGCTACAAACTTTGCAACGAGCTGGTCAGCTGCTTCCAGGAGCGGGCCAAGATTGAGAAGGGCTATTCCCAGCAGCTGAGCGACTGGGCTAGGAAATGGAGAGGCATTGTGGAGAAAG",
            "go": {
                "GO:0097320": {
                    "term": "plasma membrane tubulation",
                    "domain": "biological_process",
                    "sources": []
                }
            }
        }
    ]
}
```

Field descriptions:
* `exons`, `cds_exons`, `utr5_exons`, `utr3_exons`: Lists of exons, exons from the coding sequence (CDS), and exons from the 5' and 3' UTRs. Each exon is a list of start and end genomic coordinates. All coordinates are 0-based and relative to the **forward** genomic strand. In the example above, the transcript ENSDART00000171909 starts at position 54260003 on chromosome 7, this will be equal to the first exon start.
* `XX_on_transcript`: Lists of exons, CDS etc coordinates translated to **transcript** coordinates. The first exon start is equal to 0 and the last exon end is equal to the length of the transcript. Translation includes the strand of the transcript: coordinates are forward to the transcript making these coordinates directly usable on the transcript sequence stored in `seq`.
* `cds_exons_frame`: Frame of the first nucleotide. With the coding sequence ATGGCA, the following 4 exons would have frame 0, 1, 2 and 0 respectively:
    ```
      0
      |
      ATG-GCA
    
       1
       |
       TG-GCA
    
        2
        |
        G-GCA
    
          0
          |
          GCA
    ```
* `seq` contains the transcript sequence.
* `go` holds [Gene Ontology (GO)](http://geneontology.org) terms, domains and sources. *Optional*: GO is only imported with the `--go` option from the `import_ensembl` script.

## Future FON version

Future versions might address limitation of the currently available FON version. For example, FON1 doesn't allow features to be stored hierarchically; they are stored in a list. Contributions to add new FON versions are welcome.

## Download

See [tags](/../../tags) page.

## Install

```bash
pip3 install fontools
```

If you don't have root permission, install in your home using `--user` option:
```bash
pip3 install fontools --user
```
Scripts are installed in `$HOME/.local/bin`, which should be added to your shell [PATH](https://linuxhint.com/path_in_bash) to run the scripts. After adding for example `$HOME/.local/bin` to your [PATH](https://linuxhint.com/path_in_bash), try:
```bash
import_ensembl -h
```
If you get an error message like `import_ensembl: command not found`, then your [PATH](https://linuxhint.com/path_in_bash) isn't properly configured.

FONtools depend on [pyfaidx](https://github.com/mdshw5/pyfaidx) for reading FASTA and [pyfnutils](https://gitlab.com/vejnar/pyfnutils) for logging.

## Scripts

| Script                                   | Description                                                 |
| ---------------------------------------- | ----------------------------------------------------------- |
| FON                                                                                                    |
| [import_ensembl](#import-ensembl)        | Import Ensembl sequence and annotations                     |
| [fon_import](#import-annotations-to-fon) | Import annotations to FON (from GFF3 for now)               |
| [fon_transform](#transform-fon)          | Transform FON file                                          |
| FON/GFF3/FASTA/TAB                                                                                     |
| [merge_annot](#fon-and-other-formats)    | Merge FON/GFF3/FASTA files                                  |
| [ensembl2ucsc](#fon-and-other-formats)   | Convert names from Ensembl to UCSC (in FASTA, GFF3 and tab) |
| FASTA                                                                                                  |
| [fasta_format](#fasta-tools)             | Format and/or Sort FASTA file (split sequence)              |
| [fon_mask_fasta](#fasta-tools)           | Mask sequence (FASTA) using FON                             |
| [fasta_seq_length](#fasta-tools)         | Create tab file with sequence(s) length from FASTA file     |

## Import Ensembl

The `import_ensembl` script creates and maintains an Ensembl-based annotation repository including:
* Download Ensembl annotations and sequences
* Parse GFF to create FON annotations
* Map chromosome and contig names to create genome tracks in UCSC
* Create indices for read-mapping software.

Annotations are imported using `fon_import`, then `fon_transform` is used:
* To create FON files restricted to a biotype, for example *protein coding* transcripts,
* To create FON files selecting the longest isoform of each gene,
* To create "metagene" FON files obtained by merging all isoforms of a gene together. Example of how a metagene is obtained from 3 isoforms:

    ![Metagene](img/metagene.webp)

    These "metagenes" can be used for counting HTS reads per gene, where reads mapping to any isoforms will map to the metagene.

The script is compatible with [Ensembl](http://ensembl.org) and [Ensembl Genomes](https://ensemblgenomes.org) (see option `--division`/`-n`).

### Config

The `import_ensembl` script aims to maintain a local Ensembl-based repository. Using it requires to set multiple options. But most of these options will be the same each time `import_ensembl` is used. In most cases, the data will always be stored in the same directories and only options specifying the release number or the species will change and be specified on the command-line. To this end, all `import_ensembl` options can be set for convenience in a JSON config file, in addition to the command-line. This config file can be placed:
* Either one of the two following directories (step 2 below):
    * The directory defined by the environment variable `$HTS_CONFIG_PATH`. `$HTS_CONFIG_PATH` can be defined by the user.
    * The directory defined by the environment variable `$XDG_CONFIG_HOME/hts`. `$XDG_CONFIG_HOME` is defined by your desktop [environment](https://wiki.archlinux.org/title/XDG_Base_Directory).
* Or, you can use the `---path_config` option to set the directory where to find a `fontools.json` config file. This option is not used in this tutorial.

To configure `import_ensembl` script:
1. Create the root/main directory. All downloaded, FON, sequences, etc files will be stored in this directory:
    ```bash
    mkdir /data/sai
    mkdir /data/sai/download
    ```
    In the example above, we are using the `/data/sai` directory (*sai* stands for Sequence Annotations & Indices). This is intended for system-wide installation. **Alternatively**, you can change it to the directory of your choice, for example, if you want to store the data in a `sai` directory in your home:
    ```bash
    mkdir ~/sai
    mkdir ~/sai/download
    ```
2. To configure a config directory, add in your `~/.bashrc`:
    ```bash
    export HTS_CONFIG_PATH="/etc/hts"
    ```
    **Alternatively**, in case you created a `sai` directory in your home:
    ```bash
    mkdir ~/sai/config
    ```
    Then, to set the `HTS_CONFIG_PATH` environment variable, add in your `~/.bashrc`:
    ```bash
    export HTS_CONFIG_PATH="$HOME/sai/config"
    ```
3. Create a `fontools.json` config file in `/etc/hts`:
    ```json
    {
        "fontools_path_main": "/data/sai",
        "fontools_path_download": "/data/sai/download"
    }
    ```
    **Alternatively**, in case you created a `sai` directory in your home (please replace `smith` by your username), create a `fontools.json` config file in `~/sai/config`:
    ```json
    {
        "fontools_path_main": "/home/smith/sai",
        "fontools_path_download": "/home/smith/sai/download"
    }
    ```
4. *Optional*. To automatically keep a log of `import_ensembl` actions, you can create the following directory. This will automatically create a different log file per Ensembl release. To specify the location for the log from the script, use the `--path_log`/`-l` option:
    ```bash
    mkdir /data/sai/log
    ```
    **Alternatively**, in case you created a `sai` directory in your home:
    ```bash
    mkdir ~/sai/log
    ```
5. *Optional*. If mapping to UCSC chromosome and contig names are desired, download mapping from the [ChromosomeMappings](https://github.com/dpryan79/ChromosomeMappings) repository:
    ```bash
    mkdir /data/sai/annots
    cd /data/sai/annots
    ```
    **Alternatively**, in case you created a `sai` directory in your home:
    ```bash
    mkdir ~/sai/annots
    cd ~/sai/annots
    ```
    Then:
    ```bash
    wget https://github.com/dpryan79/ChromosomeMappings/archive/refs/heads/master.tar.gz
    tar xvfz master.tar.gz
    rm -f master.tar.gz
    mv ChromosomeMappings-master ChromosomeMappings
    ```
    Add the `fontools_path_mapping` to `fontools.json` config file:
    ```json
    {
        "fontools_path_main": "/data/sai",
        "fontools_path_download": "/data/sai/download",
        "fontools_path_mapping": "/data/sai/annots/ChromosomeMappings"
    }
    ```
    **Alternatively**, in case you created a `sai` directory in your home (please replace `smith` by your username):
    ```json
    {
        "fontools_path_main": "/home/smith/sai",
        "fontools_path_download": "/home/smith/sai/download",
        "fontools_path_mapping": "/home/smith/sai/annots/ChromosomeMappings"
    }
    ```
    Although using name mapping isn't required, it's recommended. Using name mapping, files containing chromosome lengths with UCSC names will be created. These files are essential to create UCSC genome browser tracks.

### Usage

If you haven't set the environment variable `HTS_CONFIG_PATH` (see above), then:
* Use the `---path_config` option to set the directory to find a `fontools.json` config file or,
* Use the `--fontools_path_main` and `--fontools_path_download` on the command line.

To list available species, use (for Ensembl 104):
```bash
import_ensembl -r 104 -s list
```

To get Ensembl 104 data for 4 species using 10 cores:
```bash
import_ensembl -r 104 -s danio_rerio,saccharomyces_cerevisiae,homo_sapiens,mus_musculus -p 10
```

To select what data are generated, use the `--steps`/`-t` option. Currently, the following steps are available:
* `genome` step download FASTA genome sequences, map chromosome/contig names to UCSC names if requested, sort FASTA files, and create chromosome length file.
* `gene` step download GFF annotations, import them to FON files, map chromosome/contig names to UCSC names if requested, and create FON files with:
    * metagenes (isoform union, see above),
    * longest transcripts.
* `bowtie2` and `star` to create indices for [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2) and [STAR](https://github.com/alexdobin/STAR) respectively.
* `all` of the above steps. This is the default. 

To import terms, domains and sources from [Gene Ontology (GO)](http://geneontology.org), add the `--go` option.

### Example

To import from Ensembl (`-n`) release 104 (`-r`), get FASTA and GFF (`-t genome,gene`) and convert to FON:
```bash
import_ensembl -n ensembl -r 104 -s caenorhabditis_elegans -t genome,gene
```

This command will print detailed log that are recorded in `log/ensembl104.log`:
```bash
2021-05-06 13:53:05,501 - import_ensembl - INFO - Starting (caenorhabditis_elegans,104)
2021-05-06 13:53:06,303 - import_ensembl - INFO - Found assembly WBcel235 (toplevel)
2021-05-06 13:53:06,304 - import_ensembl - INFO - Downloading fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
2021-05-06 13:53:49,983 - import_ensembl - INFO - Downloading gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.104.gff3.gz
2021-05-06 13:54:02,400 - import_ensembl - INFO - Downloading fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz
2021-05-06 13:54:22,319 - import_ensembl - INFO - Downloading fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz
2021-05-06 13:54:25,103 - import_ensembl - INFO - Sorting to /data/sai/seqs/caeele_genome_all_ensembl_wbcel235.fa
2021-05-06 13:54:25,104 - import_ensembl - INFO - Start ['fasta_format', '--sort', '--input', '/data/sai/download/ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz', '--output', '/data/sai/seqs/caeele_genome_all_ensembl_wbcel235.fa']
2021-05-06 13:54:28,613 - import_ensembl - INFO - Start ['cp', '/data/sai/download/ftp.ensembl.org/pub/release-104/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.104.gff3.gz', '/data/sai/annots/caeele_cdna_all_ensembl104.gff3.gz']
2021-05-06 13:54:28,617 - import_ensembl - INFO - Start ['gzip', '-d', '/data/sai/annots/caeele_cdna_all_ensembl104.gff3.gz']
2021-05-06 13:54:28,858 - import_ensembl - INFO - Creating chromosome length file /data/sai/annots/caeele_genome_all_ensembl_wbcel235_chrom_length.tab
2021-05-06 13:54:28,859 - import_ensembl - INFO - Start ['fasta_seq_length', '--input', '/data/sai/seqs/caeele_genome_all_ensembl_wbcel235.fa', '--output', '/data/sai/annots/caeele_genome_all_ensembl_wbcel235_chrom_length.tab']
2021-05-06 13:54:29,322 - import_ensembl - INFO - Creating chromosome length file for UCSC /data/sai/annots/caeele_genome_all_ensembl_wbcel235_ucsc_names_chrom_length.tab
2021-05-06 13:54:29,322 - import_ensembl - INFO - Start ['ensembl2ucsc', '--input', '/data/sai/annots/caeele_genome_all_ensembl_wbcel235_chrom_length.tab', '--output', '/data/sai/annots/caeele_genome_all_ensembl_wbcel235_ucsc_names_chrom_length.tab', '--path_mapping', '/data/sai/annots/ChromosomeMappings/WBcel235_ensembl2UCSC.txt']
2021-05-06 13:54:29,344 - import_ensembl - INFO - Importing annotation
2021-05-06 13:54:29,344 - import_ensembl - INFO - Start ['fon_import', '--annotation', '/data/sai/download/ftp.ensembl.org/pub/release-104/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.104.gff3.gz', '--data_source', 'ensembl', '--fasta', '/data/sai/download/ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz', '--fasta', '/data/sai/download/ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz', '--cdna', '--exclude_no_seq', '--biotype', 'all,protein_coding', '--output', '/data/sai/annots/caeele_cdna_${biotype}_ensembl104.fon${version}.json', '--output_format', 'fon']
2021-05-06 13:54:39,531 - import_ensembl - INFO - Transform FON (union,protein_coding)
2021-05-06 13:54:39,531 - import_ensembl - INFO - Start ['fon_transform', '--fon', '/data/sai/annots/caeele_cdna_protein_coding_ensembl104.fon1.json', '--method', 'union', '--output', '/data/sai/annots/caeele_cdna_union2gene_protein_coding_ensembl104.fon${version}.json']
2021-05-06 13:54:42,214 - import_ensembl - INFO - Transform FON (longest,protein_coding)
2021-05-06 13:54:42,214 - import_ensembl - INFO - Start ['fon_transform', '--fon', '/data/sai/annots/caeele_cdna_protein_coding_ensembl104.fon1.json', '--method', 'longest', '--output', '/data/sai/annots/caeele_cdna_longest_transcript_protein_coding_ensembl104.fon${version}.json']
2021-05-06 13:54:44,178 - import_ensembl - INFO - Transform FON (union,all)
2021-05-06 13:54:44,178 - import_ensembl - INFO - Start ['fon_transform', '--fon', '/data/sai/annots/caeele_cdna_all_ensembl104.fon1.json', '--method', 'union', '--output', '/data/sai/annots/caeele_cdna_union2gene_all_ensembl104.fon${version}.json']
2021-05-06 13:54:48,050 - import_ensembl - INFO - Transform FON (longest,all)
2021-05-06 13:54:48,050 - import_ensembl - INFO - Start ['fon_transform', '--fon', '/data/sai/annots/caeele_cdna_all_ensembl104.fon1.json', '--method', 'longest', '--output', '/data/sai/annots/caeele_cdna_longest_transcript_all_ensembl104.fon${version}.json']
```

The following files will be created:
```
├── annots
│   ├── caeele_cdna_all_ensembl104.fon1.json          <-- All transcripts
│   ├── caeele_cdna_all_ensembl104.gff3               <-- All transcripts (GFF3)
│   ├── caeele_cdna_longest_transcript_all_ensembl104.fon1.json              <-- Longest transcript per all gene
│   ├── caeele_cdna_longest_transcript_protein_coding_ensembl104.fon1.json   <-- Longest transcript per protein-coding gene
│   ├── caeele_cdna_protein_coding_ensembl104.fon1.json                      <-- All transcripts of protein-coding gene
│   ├── caeele_cdna_union2gene_all_ensembl104.fon1.json                      <-- Metagenes of all genes
│   ├── caeele_cdna_union2gene_protein_coding_ensembl104.fon1.json           <-- Metagenes of protein-coding genes
│   ├── caeele_genome_all_ensembl_wbcel235_chrom_length.tab                  <-- Chromosome lengths (TAB)
│   ├── caeele_genome_all_ensembl_wbcel235_ucsc_names_chrom_length.tab       <-- Chromosome lengths with UCSC names (TAB)
│   └── ChromosomeMappings                            <-- Ensembl to/from UCSC name mapping
│       .
│       .
│       └── Zv9_UCSC2ensembl.txt
├── config
│   └── fontools.json
├── download
│   └── ftp.ensembl.org
│       └── pub
│           └── release-104
│               ├── fasta
│               │   └── caenorhabditis_elegans
│               │       ├── cdna
│               │       │   └── Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz
│               │       ├── dna
│               │       │   └── Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
│               │       └── ncrna
│               │           └── Caenorhabditis_elegans.WBcel235.ncrna.fa.gz
│               └── gff3
│                   └── caenorhabditis_elegans
│                       └── Caenorhabditis_elegans.WBcel235.104.gff3.gz
├── log
│   └── ensembl104.log
└── seqs
    ├── caeele_genome_all_ensembl_wbcel235.fa             <-- Sequence (FASTA)
    └── caeele_genome_all_ensembl_wbcel235.fa.fai         <-- FASTA index
```

## Import annotations to FON

Annotations can be imported into FON from any GFF source. For example, to import gene annotations for *Xenopus tropicalis* from [Xenbase](http://www.xenbase.org):
1. Download [GFF3](http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_Xenbase.gff3):
    ```bash
    cd /data/sai/downloads
    wget -m http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_Xenbase.gff3
    ```
2. Download [FASTA](http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz):
    ```bash
    cd /data/sai/downloads
    wget -m http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz
    ```
    * Providing FASTA sequence is optional. If not FASTA file is provided, the resulting FON file won't include any sequence.
    * Instead of genomic sequence, the transcript sequences can be provided. Use the `--cdna` option to specify FASTA file containts cDNA instead of genomic sequence. 
3. To import the GFF3 annotation to FON, run (change paths of files to your setup):
    ```bash
    fon_import --annotation "/data/sai/downloads/ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_Xenbase.gff3" \
               --output '/data/sai/annots/xentro_cdna_${biotype}_xenbase100.fon${version}.json' \
               --fasta "/data/sai/downloads/ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz" \
               --output_format fon \
               --biotype all
    ```
    * The `--output` is a path written as a simple string or a template string ([string.Template](https://docs.python.org/3/library/string.html#template-strings)).

## Transform FON

Transform FON files: select longest isoform or merge isoforms. For example to select the longest isoform of each gene:
```bash
fon_transform --fon "/data/sai/annots/caeele_cdna_protein_coding_ensembl104.fon1.json" \
              --method "longest" \
              --output '/data/sai/annots/caeele_cdna_longest_transcript_protein_coding_ensembl104.fon${version}.json'
```

## FON and other formats

An easy way to merge the sequence and annotations of multiple species together is to input each sequence and annotation in a comma-separated list of files, using the `merge_annot` script. In this example annotations and sequences from zebrafish and yeast are merged:
```bash
cd /data/sai
merge_annot --input_fasta "seqs/zebrafish.fa,seqs/yeast.fa" \
            --output_fasta "seqs/zebrafish_plus_yeast.fa" \
            --input_gff "annots/zebrafish.gff3,annots/yeast.gff3" \
            --output_gff "annots/zebrafish_plus_yeast.gff3" \
            --input_fon "annots/zebrafish.fon1.json,annots/yeast.fon1.json" \
            --output_fon "annots/zebrafish_plus_yeast.fon1.json"
```

The script `ensembl2ucsc` can be used to translate chromosome/contig names from Ensembl to UCSC names (for *C. Elegans*) using [ChromosomeMappings](https://github.com/dpryan79/ChromosomeMappings):
```bash
cd /data/sai/annots
ensembl2ucsc --input "caeele_genome_all_ensembl_wbcel235_chrom_length.tab" \
             --output "caeele_genome_all_ensembl_wbcel235_ucsc_names_chrom_length.tab" \
             --path_mapping "ChromosomeMappings/WBcel235_ensembl2UCSC.txt"
```

## FASTA tools

* `fasta_format`: Format FASTA file
    * Sort entries of FASTA file (`--sort`)
    * Split sequences in lines of desired length (`--seq_length`)

* `fon_mask_fasta`: Mask part(s) of sequence (FASTA) with Ns
    ```bash
    fon_mask_fasta --input_fon "selected_loci.fon1.json" \
                   --input_fasta "genome.fa" \
                   --output_fasta "genome_mask.fa" \
                   --extension "50" \
                   --exterior_extension "100"
    ```
    * A list of interval coordinates from FON (`--input_fon`) are masked with Ns in the sequence (`--input_fasta`). By default, the list of interval coordinates from the `exons` key of each feature in the FON file is used. To use a different key, use the `--key` option; for example use `--key "cds_exons"` to mask the coding sequences.
    * Each interval can be extended by `--extension` and each feature can be extended by `--exterior_extension`. `--exterior_extension` value is by default equal to `--extension` value. For example, using `--extension 2` on [[10,15], [20, 30]] will mask [[8,17], [18, 32]], while `--exterior_extension 2` will mask [[8,15], [20, 32]].
    * Using `--inverse`, only interval coordinates from FON are kept intact, the rest of the sequence is replaced by Ns.
    * Strand of FON features is ignored.

* `fasta_seq_length`: Create tabulated file with sequence(s) name and length from FASTA file.

## License

*FONtools* are distributed under the Mozilla Public License Version 2.0 (see /LICENSE).

Copyright (C) 2015-2022 Charles E. Vejnar
