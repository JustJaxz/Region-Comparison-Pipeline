<body>
<h1>Metabarcoding Pipeline</h1>
<p> 
  This pipeline was used to process and analyse metabarcoding data for the publication: </p>
<p>
  Stuart, J., et al., <i> A comparison of two gene regions for assessing community composition of eukaryotic marine microalgae from coastal ecosystems</i>. 
Scientific Reports, 2024. <b>14</b>(1): p. 6442. doi: <a href="http://dx.doi.org/10.1038/s41598-024-56993-4">10.1038/s41598-024-56993-4</a> 
</p>
<img src="images/graphicalAbstract.jpg" alt="graphicalAbstarct" />

<h2> Data types and file layout </h2>
<p>
  There were three main types of input files for the analysis in this pipeline:
  <ol>
    <li>Metabarcoding files</li>
    <li>Taxonomy files</li>
    <li>Metadata</li>
  </ol>
  </p>
<h3>Metadata layout</h3>
<p>Sample/sequence metadata files were in the format depicted below</p>
<img src="images/metadata.png" alt="metadata" />

<h2>Overhang adaptors</h2>
<p>Primer pairs used for this paper were each modified for Illumina sequencing with the addition of overhang adaptors. The Illumina overhang adapter sequences to be added to locus‐specific sequences are as follows.
  <br>
  <br>
  <b>Forward overhang:</b><br>
5’ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG ‐ [locus specific sequence] <br>
<b>Reverse overhang:</b><br>
5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ‐ [locus specific sequence] <br>
  <br>
Full documentation can be found in the 
  <a href="http://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf">Illumina metagenomic library prep guide</a> 
  </p>

</body>
