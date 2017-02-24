## VAP (Variant Assurance Pipeline, version 0.1b)
designed for Multi Region Sequencing (MRS) of Tumor samples

The purpose of this repository is to:
* Recheck variants in bam: examine the mapping features surrounding any genomic coordinates in the raw read alignment files (.bam).
* Sensitive SSNV classification: based on these extracted features, and multiple region samples (MRS, if available), increase the sensitivity for calling low frequency variants while also suppresing sites with FP-suggestive mapping features.
* SCNA calling and tumor purity estimation: achived by TitanCNA, jointly estimate local CN and tumor purity.
* Tumor growth dynamics: the analysis folder contains codes for inspecting tumor MRS variant allele frequency distributions (SFS) and calculating ITH summary statistics for subclonal variants determined based on MRS.


Learn More
---
Tumor Multi Region Sequencing (MRS) is becoming a valuable resource for inspecting intra tumor heterogeneity reflecting growth dynamics in the expansion after tumor initiation. However, such information is buried in subclonal variants which can be at low frequency even if the tumor sample is relatively pure. Detection of these events can be further complicated by the uneven read depth of coverage in the current NGS experiments. Here we extract mapping features surrounding each genomic coordinates of interest, leveraging information across MRS, to reach a higher sensitivity while also controlling the false positive rate. Through the identification of subclonal variants based on MRS and CN/purity adjustment of resultant VAF, SFS of subclonal variants in multiple regions can be evaluated, and measurements of genetic divergence between regions can be computed to study the dynamics in tumor expansion.


Dependencies and Annotation Files
---
* cpan modules: ``Statistics::Basic`` ``Math::CDF`` ``Parallel::ForkManager`` ``Text::NSP::Measures::2D::Fisher::right``
* See details in ``config.tsv`` file. More info to be added.


Installation
---

Installed Bamtools (https://github.com/pezmaster31/bamtools). Also make sure you have zlib and boost installed.

run make in following way to install

	$ make BAMTOOLS_ROOT=/bamtools_directory/ ZLIB_ROOT=/zlib_directory/ BOOST_ROOT=/boost_directory/

The binaries will be built at bin/. ``xxx_directory`` is where lib/ and include/ sub-directories of xxx (bamtools, zlib and boost) are located.


Usage
---

VAP can be run with UNIX command-line interface.

**Getting help message**

	$ perl VAP/bin/DTrace.pl -h (or --help)



Contact Author
---
Sun Ruping

If you are willing to receive updates and notice, send an email to ruping@stanford.edu.

Current Affiliation:
Curtis Lab, Department of Genetics and Medicine, Stanford University. CA, USA.

Email: ruping@stanford.edu or regularhand@gmail.com
