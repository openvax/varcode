# Change Log

## [v0.5.15](https://github.com/hammerlab/varcode/tree/v0.5.15) (2017-04-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.14...v0.5.15)

**Closed issues:**

- Allow contig to be empty to allow parsing of MAF with faulty mutation\(s\) [\#210](https://github.com/hammerlab/varcode/issues/210)

**Merged pull requests:**

- Fixes to load\_maf [\#223](https://github.com/hammerlab/varcode/pull/223) ([tavinathanson](https://github.com/tavinathanson))
- added raise\_on\_error option to load\_maf and load\_maf\_dataframe [\#221](https://github.com/hammerlab/varcode/pull/221) ([iskandr](https://github.com/iskandr))
- Optionally allow duplicated mutations when using load\_vcf or load\_maf. Fixes \#211  [\#212](https://github.com/hammerlab/varcode/pull/212) ([tuomastik](https://github.com/tuomastik))

## [v0.5.14](https://github.com/hammerlab/varcode/tree/v0.5.14) (2017-04-05)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.13...v0.5.14)

**Merged pull requests:**

- Adding 'distinct' as a parameter to load\_vcf. [\#222](https://github.com/hammerlab/varcode/pull/222) ([julia326](https://github.com/julia326))

## [v0.5.13](https://github.com/hammerlab/varcode/tree/v0.5.13) (2017-04-01)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.12...v0.5.13)

**Closed issues:**

- Effect prediction throws error \(even if raise\_on\_error=False\) [\#213](https://github.com/hammerlab/varcode/issues/213)
- Optionally allow duplicated mutations when using load\_vcf or load\_maf [\#211](https://github.com/hammerlab/varcode/issues/211)

**Merged pull requests:**

- install ensembl 87 on travis [\#219](https://github.com/hammerlab/varcode/pull/219) ([iskandr](https://github.com/iskandr))
- Allow user to affect the sorting of variants when loading a VCF or MAF. [\#218](https://github.com/hammerlab/varcode/pull/218) ([tuomastik](https://github.com/tuomastik))

## [v0.5.12](https://github.com/hammerlab/varcode/tree/v0.5.12) (2017-01-18)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.11...v0.5.12)

**Closed issues:**

- Make Varcode correctly infer genome for b37-decoy string [\#207](https://github.com/hammerlab/varcode/issues/207)
- Longer indels in random variants [\#47](https://github.com/hammerlab/varcode/issues/47)
- Predict coding sequence of StartLoss mutations [\#4](https://github.com/hammerlab/varcode/issues/4)

**Merged pull requests:**

- Add optional\_cols list to load\_maf [\#209](https://github.com/hammerlab/varcode/pull/209) ([tavinathanson](https://github.com/tavinathanson))

## [v0.5.11](https://github.com/hammerlab/varcode/tree/v0.5.11) (2016-12-05)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.10...v0.5.11)

**Fixed bugs:**

- Varcode noncoding variant in a drop\_silent\_and\_noncoding\(\) list [\#200](https://github.com/hammerlab/varcode/issues/200)

**Merged pull requests:**

- Adding aa\_ref argument to StopLoss for variants which delete codons before stop [\#203](https://github.com/hammerlab/varcode/pull/203) ([iskandr](https://github.com/iskandr))

## [v0.5.10](https://github.com/hammerlab/varcode/tree/v0.5.10) (2016-10-19)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.9...v0.5.10)

**Fixed bugs:**

- Variant pickling won't work for not-human and non-EnsemblRelease Genomes [\#147](https://github.com/hammerlab/varcode/issues/147)

**Closed issues:**

- Link on PyPI badge broken [\#191](https://github.com/hammerlab/varcode/issues/191)
- Reference incorrectly inferred when "b36" in reference file path [\#181](https://github.com/hammerlab/varcode/issues/181)
- Premature stop codon error [\#166](https://github.com/hammerlab/varcode/issues/166)

**Merged pull requests:**

- explicit args to \_\_init\_\_ of Intronic splice effects fixes serialization [\#199](https://github.com/hammerlab/varcode/pull/199) ([iskandr](https://github.com/iskandr))
- Update RELEASING.md, fixing tagging instructions [\#198](https://github.com/hammerlab/varcode/pull/198) ([julia326](https://github.com/julia326))

## [v0.5.9](https://github.com/hammerlab/varcode/tree/v0.5.9) (2016-10-11)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.8...v0.5.9)

**Fixed bugs:**

- StopLoss pickling is broken [\#188](https://github.com/hammerlab/varcode/issues/188)

**Closed issues:**

- One logger per module [\#196](https://github.com/hammerlab/varcode/issues/196)
- SNV results in deletion [\#193](https://github.com/hammerlab/varcode/issues/193)

**Merged pull requests:**

- One logger per module. [\#197](https://github.com/hammerlab/varcode/pull/197) ([julia326](https://github.com/julia326))
- Fix edge case where PrematureStop in last amino acid got interpreted as a Deletion  [\#194](https://github.com/hammerlab/varcode/pull/194) ([iskandr](https://github.com/iskandr))
- Fix inferred-reference-bug [\#182](https://github.com/hammerlab/varcode/pull/182) ([jburos](https://github.com/jburos))

## [v0.5.8](https://github.com/hammerlab/varcode/tree/v0.5.8) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.7...v0.5.8)

**Merged pull requests:**

- changed Markdown image links to HTML [\#192](https://github.com/hammerlab/varcode/pull/192) ([iskandr](https://github.com/iskandr))

## [v0.5.7](https://github.com/hammerlab/varcode/tree/v0.5.7) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.3...v0.5.7)

## [v0.5.3](https://github.com/hammerlab/varcode/tree/v0.5.3) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.2...v0.5.3)

**Merged pull requests:**

- Use conda to install pypandoc [\#189](https://github.com/hammerlab/varcode/pull/189) ([arahuja](https://github.com/arahuja))
- Ensure README.md is packaged [\#186](https://github.com/hammerlab/varcode/pull/186) ([arahuja](https://github.com/arahuja))
- Upgrade serializable dependency with tests [\#185](https://github.com/hammerlab/varcode/pull/185) ([arahuja](https://github.com/arahuja))
- Add pypi badge [\#184](https://github.com/hammerlab/varcode/pull/184) ([arahuja](https://github.com/arahuja))

## [v0.5.2](https://github.com/hammerlab/varcode/tree/v0.5.2) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.1...v0.5.2)

**Closed issues:**

- Make sure `Variant` works with any Genome \(not just a human EnsemblRelease\) [\#127](https://github.com/hammerlab/varcode/issues/127)

**Merged pull requests:**

- Move extraneous variables to properties for normalization [\#190](https://github.com/hammerlab/varcode/pull/190) ([tavinathanson](https://github.com/tavinathanson))
- Use is\_protein\_coding property of pyensembl.Transcript and pyensembl.Gene [\#180](https://github.com/hammerlab/varcode/pull/180) ([iskandr](https://github.com/iskandr))

## [v0.5.1](https://github.com/hammerlab/varcode/tree/v0.5.1) (2016-09-16)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.0...v0.5.1)

**Merged pull requests:**

- Add MutationEffect to \_\_init\_\_.py [\#178](https://github.com/hammerlab/varcode/pull/178) ([timodonnell](https://github.com/timodonnell))

## [v0.5.0](https://github.com/hammerlab/varcode/tree/v0.5.0) (2016-09-13)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.19...v0.5.0)

**Implemented enhancements:**

- Support collection.as\_dataframe\(\) [\#128](https://github.com/hammerlab/varcode/issues/128)

**Closed issues:**

- Substitution mis-annotated as stop-loss [\#176](https://github.com/hammerlab/varcode/issues/176)
- Wrong aa\_mutation\_end\_offset for insertion of stop codon [\#175](https://github.com/hammerlab/varcode/issues/175)
- Wrong aa\_ref for insertion of stop codon [\#174](https://github.com/hammerlab/varcode/issues/174)
- Insertions after the stop codon annotated as plain Insertions [\#172](https://github.com/hammerlab/varcode/issues/172)
- Mutations before the stop codon confused as StopLosses [\#171](https://github.com/hammerlab/varcode/issues/171)
- StopLosses do not translate into 3' UTR [\#170](https://github.com/hammerlab/varcode/issues/170)
- Insertion of stop codon is annotated as simple Insertion and not PrematureStop [\#169](https://github.com/hammerlab/varcode/issues/169)
- Synonimous FrameShift over stop codon not annotated as silent  [\#168](https://github.com/hammerlab/varcode/issues/168)
- Wrong offset for insertion of StopCodon [\#167](https://github.com/hammerlab/varcode/issues/167)
- Document release process [\#154](https://github.com/hammerlab/varcode/issues/154)
- compare variants that use different references [\#83](https://github.com/hammerlab/varcode/issues/83)
- Annotate with predicted pathogenicity [\#46](https://github.com/hammerlab/varcode/issues/46)

**Merged pull requests:**

- Reorganize effect prediction code, fixed annotation bugs/issues [\#173](https://github.com/hammerlab/varcode/pull/173) ([iskandr](https://github.com/iskandr))

## [v0.4.19](https://github.com/hammerlab/varcode/tree/v0.4.19) (2016-09-12)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.18...v0.4.19)

**Fixed bugs:**

- original\_start doesn't get pickled [\#141](https://github.com/hammerlab/varcode/issues/141)

**Closed issues:**

- replace `load\_vcf` with `load\_vcf\_fast` ? [\#144](https://github.com/hammerlab/varcode/issues/144)
- Add `annotate\_random\_variants` commandline script [\#49](https://github.com/hammerlab/varcode/issues/49)
- support filtering a variant collection to variants overlapping specified gene names [\#32](https://github.com/hammerlab/varcode/issues/32)
- Use SPANR to identify splicing misregulation [\#2](https://github.com/hammerlab/varcode/issues/2)

## [v0.4.18](https://github.com/hammerlab/varcode/tree/v0.4.18) (2016-08-08)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.17...v0.4.18)

**Closed issues:**

- vcf unit tests broken in python 3 \(?\) [\#164](https://github.com/hammerlab/varcode/issues/164)
- maximum recursion depth exceeded when loading a vcf from a URL [\#163](https://github.com/hammerlab/varcode/issues/163)

**Merged pull requests:**

- In load\_vcf, when passed a URL download it first to a local file then… [\#165](https://github.com/hammerlab/varcode/pull/165) ([timodonnell](https://github.com/timodonnell))
- Removed Collection from varcode, moved to separate 'sercol' repo instead [\#162](https://github.com/hammerlab/varcode/pull/162) ([iskandr](https://github.com/iskandr))

## [v0.4.17](https://github.com/hammerlab/varcode/tree/v0.4.17) (2016-08-05)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.16...v0.4.17)

**Merged pull requests:**

- Commandline interface, simplified serialization, merging VariantCollections [\#161](https://github.com/hammerlab/varcode/pull/161) ([iskandr](https://github.com/iskandr))

## [v0.4.16](https://github.com/hammerlab/varcode/tree/v0.4.16) (2016-07-30)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.15...v0.4.16)

## [v0.4.15](https://github.com/hammerlab/varcode/tree/v0.4.15) (2016-07-15)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.14...v0.4.15)

**Fixed bugs:**

- Fix clone\_with\_new\_elements for VariantCollection [\#159](https://github.com/hammerlab/varcode/pull/159) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- load\_vcf\_fast fails when sample names contain spaces [\#158](https://github.com/hammerlab/varcode/issues/158)

**Merged pull requests:**

- Fix load\_vcf\_fast for sample names containing a space character [\#160](https://github.com/hammerlab/varcode/pull/160) ([timodonnell](https://github.com/timodonnell))

## [v0.4.14](https://github.com/hammerlab/varcode/tree/v0.4.14) (2016-06-07)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.12...v0.4.14)

**Merged pull requests:**

- Don't memoize so much [\#157](https://github.com/hammerlab/varcode/pull/157) ([iskandr](https://github.com/iskandr))

## [v0.4.12](https://github.com/hammerlab/varcode/tree/v0.4.12) (2016-05-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.13...v0.4.12)

## [v0.4.13](https://github.com/hammerlab/varcode/tree/v0.4.13) (2016-05-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.11...v0.4.13)

**Merged pull requests:**

- Fix versioneer prefix format [\#155](https://github.com/hammerlab/varcode/pull/155) ([armish](https://github.com/armish))

## [v0.4.11](https://github.com/hammerlab/varcode/tree/v0.4.11) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.10...v0.4.11)

## [v0.4.10](https://github.com/hammerlab/varcode/tree/v0.4.10) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.9...v0.4.10)

## [v0.4.9](https://github.com/hammerlab/varcode/tree/v0.4.9) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.8...v0.4.9)

**Closed issues:**

- Add serialization for EffectCollection and VariantCollection [\#71](https://github.com/hammerlab/varcode/issues/71)

**Merged pull requests:**

- Reorganized coding effects to use KnownAminoAcidChange base class [\#153](https://github.com/hammerlab/varcode/pull/153) ([iskandr](https://github.com/iskandr))

## [v0.4.8](https://github.com/hammerlab/varcode/tree/v0.4.8) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.2...v0.4.8)

**Fixed bugs:**

- Potentially wrong translated sequence from frameshift on mm10 [\#151](https://github.com/hammerlab/varcode/issues/151)

**Closed issues:**

- ExonicSpliceSite mutations are classified as Noncoding [\#136](https://github.com/hammerlab/varcode/issues/136)
- Filter field is not saved after loading a VCF [\#89](https://github.com/hammerlab/varcode/issues/89)
- investigate porting read evidence module to use impala [\#69](https://github.com/hammerlab/varcode/issues/69)
- Attach genotypes and other sample information to Variants [\#30](https://github.com/hammerlab/varcode/issues/30)
- support determining the evidence for a variant in a bam [\#26](https://github.com/hammerlab/varcode/issues/26)

**Merged pull requests:**

- Added unit tests for Klf6 frameshift, fix bug in frameshift translation [\#152](https://github.com/hammerlab/varcode/pull/152) ([iskandr](https://github.com/iskandr))
- Add as\_dataframe to EffectCollection [\#150](https://github.com/hammerlab/varcode/pull/150) ([arahuja](https://github.com/arahuja))
- Use versioneer to manage version number [\#149](https://github.com/hammerlab/varcode/pull/149) ([arahuja](https://github.com/arahuja))
- Fix pyvcf error from passing \_parse\_samples a tuple instead of a list [\#148](https://github.com/hammerlab/varcode/pull/148) ([timodonnell](https://github.com/timodonnell))
- Fix variant pickling [\#146](https://github.com/hammerlab/varcode/pull/146) ([tavinathanson](https://github.com/tavinathanson))
- Parse and expose sample info, including for multisample VCFs [\#145](https://github.com/hammerlab/varcode/pull/145) ([timodonnell](https://github.com/timodonnell))
- Preserve contig name [\#140](https://github.com/hammerlab/varcode/pull/140) ([iskandr](https://github.com/iskandr))
- Quotes around nucleotides in Variant representation [\#139](https://github.com/hammerlab/varcode/pull/139) ([iskandr](https://github.com/iskandr))
- added is\_deletion, is\_insertion, and is\_indel properties to variants [\#138](https://github.com/hammerlab/varcode/pull/138) ([iskandr](https://github.com/iskandr))

## [v0.4.2](https://github.com/hammerlab/varcode/tree/v0.4.2) (2016-02-25)
**Implemented enhancements:**

- VariantCollection.high\_priority\_effect != Variant.top\_effect [\#58](https://github.com/hammerlab/varcode/issues/58)
- Improves the documentation for varcode [\#110](https://github.com/hammerlab/varcode/pull/110) ([armish](https://github.com/armish))
- Convert effect-type section into a sorted table [\#104](https://github.com/hammerlab/varcode/pull/104) ([armish](https://github.com/armish))
- Start highlighting Python syntax in README [\#103](https://github.com/hammerlab/varcode/pull/103) ([armish](https://github.com/armish))

**Fixed bugs:**

- Varcode requires pandas \>= 0.13.1, however it uses 0.15 functionality \#12 [\#92](https://github.com/hammerlab/varcode/issues/92)
- Varcode version 0.3.10 cannot be imported when installed through pip [\#90](https://github.com/hammerlab/varcode/issues/90)
- pip installing Varcode doesn't seem to work lately [\#84](https://github.com/hammerlab/varcode/issues/84)
- AttributeError: 'FrameShiftTruncation' object has no attribute 'aa\_alt' [\#70](https://github.com/hammerlab/varcode/issues/70)
- Use find\_packages correctly [\#85](https://github.com/hammerlab/varcode/pull/85) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- memoize a bit less [\#131](https://github.com/hammerlab/varcode/issues/131)
- Intragenic variants do not have a short\_description field [\#129](https://github.com/hammerlab/varcode/issues/129)
- move read\_evidence module and Locus class to varlens [\#124](https://github.com/hammerlab/varcode/issues/124)
- Support Structural Variants [\#122](https://github.com/hammerlab/varcode/issues/122)
- PrematureStop called as Silent [\#116](https://github.com/hammerlab/varcode/issues/116)
- PrematureStop called as a Deletion [\#111](https://github.com/hammerlab/varcode/issues/111)
- UnboundLocalError in in\_frame\_coding\_effect.py [\#107](https://github.com/hammerlab/varcode/issues/107)
- Double mutations in a MAF file cause error [\#105](https://github.com/hammerlab/varcode/issues/105)
- varcode.load\_vcf\_fast used 0.16.1 Pandas options [\#101](https://github.com/hammerlab/varcode/issues/101)
- Configuring datacache default cache directory [\#98](https://github.com/hammerlab/varcode/issues/98)
- Improve the README to include some examples of working with Varcode in IPython [\#95](https://github.com/hammerlab/varcode/issues/95)
- support loading VCFs over HTTP [\#91](https://github.com/hammerlab/varcode/issues/91)
- Travis should include setup.py testing [\#86](https://github.com/hammerlab/varcode/issues/86)
- Make Variants pickle-able [\#77](https://github.com/hammerlab/varcode/issues/77)
- modifies\_coding\_sequence is always false [\#64](https://github.com/hammerlab/varcode/issues/64)
- AssertionError: aa\_ref and aa\_alt can't both be empty string [\#63](https://github.com/hammerlab/varcode/issues/63)
- Too many open files on error on getting top effect [\#62](https://github.com/hammerlab/varcode/issues/62)
- KeyError: 'reference' in load\_vcf [\#60](https://github.com/hammerlab/varcode/issues/60)
- Issue with n\_skip? [\#56](https://github.com/hammerlab/varcode/issues/56)
- Optional random seed argument for generating random variants [\#48](https://github.com/hammerlab/varcode/issues/48)
- An argument for using == and not \>= for requirements? [\#43](https://github.com/hammerlab/varcode/issues/43)
- deploy a test coverage tool [\#38](https://github.com/hammerlab/varcode/issues/38)
- Replace raise\_on\_error parameter to property of VariantCollection [\#36](https://github.com/hammerlab/varcode/issues/36)
- assertion error in infer\_coding\_effect [\#33](https://github.com/hammerlab/varcode/issues/33)
- add a memoized "highest\_priority\_effect" property to Variant [\#31](https://github.com/hammerlab/varcode/issues/31)
- support deep reloading varcode module [\#25](https://github.com/hammerlab/varcode/issues/25)
- handle multiallelic variants [\#22](https://github.com/hammerlab/varcode/issues/22)
- vcf.load\_vcf should provide an option to load all variants, regardless of whether filter is PASS [\#21](https://github.com/hammerlab/varcode/issues/21)
- empty variant collection when loading strelka vcf [\#16](https://github.com/hammerlab/varcode/issues/16)
- Incorrect handling of variants which run past the beginning/end of an exon's boundary [\#14](https://github.com/hammerlab/varcode/issues/14)
- Reference amino acid sequence sometimes empty for coding variants [\#12](https://github.com/hammerlab/varcode/issues/12)
- handle single-sample VCFs with INFO fields containing list values of size \> 1 [\#9](https://github.com/hammerlab/varcode/issues/9)
- Do FrameShift \(or StopGain\) mutations affect splicing? [\#6](https://github.com/hammerlab/varcode/issues/6)
- What to do with mutations that span the 5' UTR / CDS boundary? [\#5](https://github.com/hammerlab/varcode/issues/5)
- Annotate essential splice site mutations [\#1](https://github.com/hammerlab/varcode/issues/1)

**Merged pull requests:**

- Modest change to filtering of coding mutations include ExonicSpliceSite  [\#137](https://github.com/hammerlab/varcode/pull/137) ([iskandr](https://github.com/iskandr))
- Version bump [\#135](https://github.com/hammerlab/varcode/pull/135) ([tavinathanson](https://github.com/tavinathanson))
- Fix conda install on Travis [\#134](https://github.com/hammerlab/varcode/pull/134) ([iskandr](https://github.com/iskandr))
- Don't memoize EffectCollection.top\_priority\_effect\(\) [\#132](https://github.com/hammerlab/varcode/pull/132) ([timodonnell](https://github.com/timodonnell))
- All effects should have a default `short\_description` field [\#130](https://github.com/hammerlab/varcode/pull/130) ([armish](https://github.com/armish))
- Remove read\_evidence and locus modules [\#125](https://github.com/hammerlab/varcode/pull/125) ([timodonnell](https://github.com/timodonnell))
- Include a link to the iPython notebook in README.md [\#121](https://github.com/hammerlab/varcode/pull/121) ([armish](https://github.com/armish))
- Add varcode to Travis [\#120](https://github.com/hammerlab/varcode/pull/120) ([tavinathanson](https://github.com/tavinathanson))
- Minor problem in Variant.\_\_init\_\_ [\#119](https://github.com/hammerlab/varcode/pull/119) ([iskandr](https://github.com/iskandr))
- Update Varcode to work with new multi-species PyEnsembl  [\#118](https://github.com/hammerlab/varcode/pull/118) ([iskandr](https://github.com/iskandr))
- Fix \#116 and call PrematureStop when stop codon is added in the middle of an insertion [\#117](https://github.com/hammerlab/varcode/pull/117) ([leekaiinthesky](https://github.com/leekaiinthesky))
- Warn when variants in MAF file have wrong end position [\#115](https://github.com/hammerlab/varcode/pull/115) ([iskandr](https://github.com/iskandr))
- Bump pyensembl/varcode version [\#114](https://github.com/hammerlab/varcode/pull/114) ([tavinathanson](https://github.com/tavinathanson))
- fix logic for determining whether the protein length decreases [\#112](https://github.com/hammerlab/varcode/pull/112) ([leekaiinthesky](https://github.com/leekaiinthesky))
- decreasing 3' splice site to distance 3 from boundary [\#109](https://github.com/hammerlab/varcode/pull/109) ([iskandr](https://github.com/iskandr))
- fixed typo in effect inference, added breaking variant to unit tests [\#108](https://github.com/hammerlab/varcode/pull/108) ([iskandr](https://github.com/iskandr))
- Allow Varcode to work with mouse data via Genome [\#106](https://github.com/hammerlab/varcode/pull/106) ([tavinathanson](https://github.com/tavinathanson))
- Manually set compression in read\_vcf\_into\_dataframe [\#102](https://github.com/hammerlab/varcode/pull/102) ([timodonnell](https://github.com/timodonnell))
- Added examples to README [\#100](https://github.com/hammerlab/varcode/pull/100) ([iskandr](https://github.com/iskandr))
- depend on pandas \>= 0.15 [\#99](https://github.com/hammerlab/varcode/pull/99) ([iskandr](https://github.com/iskandr))
-  Faster VCFs loading, support HTTP, and refactored variant metadata [\#94](https://github.com/hammerlab/varcode/pull/94) ([timodonnell](https://github.com/timodonnell))
- Support for regular varcode variant instances in read evidence module [\#87](https://github.com/hammerlab/varcode/pull/87) ([timodonnell](https://github.com/timodonnell))
- Read and write json files [\#82](https://github.com/hammerlab/varcode/pull/82) ([iskandr](https://github.com/iskandr))
- JSON serialization for VariantCollection. [\#81](https://github.com/hammerlab/varcode/pull/81) ([timodonnell](https://github.com/timodonnell))
- Add short\_description field to intergenic variants [\#80](https://github.com/hammerlab/varcode/pull/80) ([timodonnell](https://github.com/timodonnell))
- Speed up PileupCollection.group\_by\_allele [\#79](https://github.com/hammerlab/varcode/pull/79) ([timodonnell](https://github.com/timodonnell))
- Variant serialization [\#78](https://github.com/hammerlab/varcode/pull/78) ([timodonnell](https://github.com/timodonnell))
- added option for genome name in load\_vcf [\#76](https://github.com/hammerlab/varcode/pull/76) ([iskandr](https://github.com/iskandr))
- Fix variant.effects\(\) to always return an EffectCollection [\#75](https://github.com/hammerlab/varcode/pull/75) ([timodonnell](https://github.com/timodonnell))
- Bump pysam dependency [\#74](https://github.com/hammerlab/varcode/pull/74) ([timodonnell](https://github.com/timodonnell))
- Cufflinks RNA filtering [\#73](https://github.com/hammerlab/varcode/pull/73) ([iskandr](https://github.com/iskandr))
- Read evidence tweaks [\#72](https://github.com/hammerlab/varcode/pull/72) ([timodonnell](https://github.com/timodonnell))
- Filter effect collection by expression [\#67](https://github.com/hammerlab/varcode/pull/67) ([iskandr](https://github.com/iskandr))
- Created EpitopeCollection, refactored effects, fix assertion failure while annotating silent stop codon [\#66](https://github.com/hammerlab/varcode/pull/66) ([iskandr](https://github.com/iskandr))
- Created EpitopeCollection, refactored effects [\#65](https://github.com/hammerlab/varcode/pull/65) ([iskandr](https://github.com/iskandr))
- include substitution in high priority effects [\#61](https://github.com/hammerlab/varcode/pull/61) ([arahuja](https://github.com/arahuja))
- don't annotate StopLoss variants that are immediately followed by another stop codon [\#57](https://github.com/hammerlab/varcode/pull/57) ([iskandr](https://github.com/iskandr))
- Refactor coding effect [\#55](https://github.com/hammerlab/varcode/pull/55) ([iskandr](https://github.com/iskandr))
- Add read\_evidence module [\#53](https://github.com/hammerlab/varcode/pull/53) ([timodonnell](https://github.com/timodonnell))
- Use transcript protein sequence [\#45](https://github.com/hammerlab/varcode/pull/45) ([iskandr](https://github.com/iskandr))
- Add contributing md [\#41](https://github.com/hammerlab/varcode/pull/41) ([iskandr](https://github.com/iskandr))
- Small coding effect refactoring and fixes [\#39](https://github.com/hammerlab/varcode/pull/39) ([iskandr](https://github.com/iskandr))
- Test problematic variants [\#37](https://github.com/hammerlab/varcode/pull/37) ([iskandr](https://github.com/iskandr))
- Typechecks and test fixes [\#35](https://github.com/hammerlab/varcode/pull/35) ([timodonnell](https://github.com/timodonnell))
- Fix maf parsing [\#34](https://github.com/hammerlab/varcode/pull/34) ([iskandr](https://github.com/iskandr))
- parse multiple alleles into distinct Variant records [\#29](https://github.com/hammerlab/varcode/pull/29) ([iskandr](https://github.com/iskandr))
- PEP8 & pyflakes fixes [\#28](https://github.com/hammerlab/varcode/pull/28) ([iskandr](https://github.com/iskandr))
- Remove pyfaidx [\#27](https://github.com/hammerlab/varcode/pull/27) ([iskandr](https://github.com/iskandr))
- Variant collection tweaks [\#24](https://github.com/hammerlab/varcode/pull/24) ([timodonnell](https://github.com/timodonnell))
- Improved vcf parsing [\#23](https://github.com/hammerlab/varcode/pull/23) ([timodonnell](https://github.com/timodonnell))
- Associate EnsemblRelease with each Variant object [\#20](https://github.com/hammerlab/varcode/pull/20) ([iskandr](https://github.com/iskandr))
- Variant collection filtering [\#19](https://github.com/hammerlab/varcode/pull/19) ([iskandr](https://github.com/iskandr))
- added IntronicSpliceSite, SpliceDonor, SpliceAcceptor effects [\#17](https://github.com/hammerlab/varcode/pull/17) ([iskandr](https://github.com/iskandr))
- collect effect annotation errors in dictionary, only look up overlapping... [\#13](https://github.com/hammerlab/varcode/pull/13) ([iskandr](https://github.com/iskandr))
- don't flatten INFO dictionary of VCF, lists are part of the field format [\#11](https://github.com/hammerlab/varcode/pull/11) ([iskandr](https://github.com/iskandr))
- Small fixes [\#10](https://github.com/hammerlab/varcode/pull/10) ([timodonnell](https://github.com/timodonnell))
- Add support for Python 3 [\#8](https://github.com/hammerlab/varcode/pull/8) ([timodonnell](https://github.com/timodonnell))
- Refactor core logic [\#7](https://github.com/hammerlab/varcode/pull/7) ([iskandr](https://github.com/iskandr))
- Classes for protein/transcript variant effects [\#3](https://github.com/hammerlab/varcode/pull/3) ([iskandr](https://github.com/iskandr))



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*