- 'id' (str) 					: the sequence identifier (unique)
- 'sp_tag' (str) 				: the signalp5 prediction (SP or OTHER. proba cutoff = 0.5, see below) (signalP5)
- 'sp_sp' (float) 				: sequence probability of having a signal peptide (SP) (signalP5)
- 'sp_other' (float)				: sequence probability not having a signal peptide (signalP5)
- 'sp_csposition' (str) 			: cleavage site position in case a SP was detected (signalP5)
- 'tm_expaa' (float)				: The expected number of amino acids in transmembrane helices. If this number is larger than 18 it is very likely 
						  to be a transmembrane protein (OR have a signal peptide) (TMHMM-2)
- 'tm_first60' (float)				: The expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein. If this 
						  number more than a few, you should be warned that a predicted transmembrane helix in the N-term could be a signal
						  peptide. Could be removed. (TMHMM-2)
- 'tm_predhel' (int) 				: number of predicted transmembrane helix  (TMHMM-2)
- 'tm_topology' (str) 				: TM tologie(s). could be removed. (TMHMM-2)
- 'tm_spcleaved' (bool 			: 0/1): whether or not a SP was detected (and cleaved) before the TM search.(signalP5 & TMHMM-2)
- 'sl_location' (str) 				: predicted subcellular location (10 possible location:  Nucleus, Cytoplasm, Extracellular, Mitochondrion, 
						  Cell membrane, Endoplasmic reticulum, Chloroplast, Golgi apparatus, Lysosome/Vacuole and Peroxisome.) (deeploc-1)
- 'sl_membrane' (float)			: probability for the protein to be membrane-bound (proba value >= 0.5) or soluble. (deeploc-1)
- 'sl_nucleus', 'sl_cytoplasm', 
  'sl_extracellular', 'sl_mitochondrion', 
  'sl_cell_membrane', 'sl_ER','sl_plastid', 
  'sl_golgi', 'sl_lyso_vac', 
  'sl_peroxysome' (float)			: subcellular locations probability. The max value among the 10 will be elected as the predicted location. (deeploc-1)
- 'seq_len' (int)				: the sequence's length
- 'mol_wt' (float)				: molecular weight of the sequence
- 'instab_idx' (float) 			: Calculate the instability index according to Guruprasad et al 1990.Implementation of the method of Guruprasad et 
						  al. 1990 to test a protein for stability Any value above 40 means the protein is unstable (has a short  half life). 
						  See: Guruprasad K., Reddy B.V.B., Pandit M.W. Protein ngineering 4:155-161(1990).(Bio.SeqUtils.ProtParam documentation).
- 'mean_vihinen_flex' (float) 		: mean value of the flexibility index according to Vihinen, 1994 (window's size : 7)."Protein average flexibility indices
						  are inversely correlated to protein stability" (https://doi.org/10.1002/prot.340190207)
- 'mean_kd_hydro' (float)			: mean value of the hydrophobicity profile calculated according to the Kyte-Doolittle scale (window's size : 6). 
						  "The Kyte-Doolittle scale is widely used for detecting hydrophobic regions in proteins. Regions with a positive value 
						  are hydrophobic. This scale can be used for identifying both surface-exposed regions as well as transmembrane regions, 
						  depending on the window size used. Short window sizes of 5-7 generally work well for predicting putative 
						  surface-exposed regions. Large window sizes of 19-21 are well suited for finding transmembrane domains if the values 
						  calculated are above 1.6 [Kyte and Doolittle, 1982]. These values should be used as a rule of thumb and deviations 
						  from the rule may occur".
						  (https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html)
- 'mean_rose_hydro' (float) 			: mean value of the hydrophobicity profile calculated according to the Rose scale (window's size : 6). "The 
						  hydrophobicity scale by Rose et al. is correlated to the average area of buried amino acids in globular proteins 
						  [Rose et al., 1985]. This results in a scale which is not showing the helices of a protein, but rather the surface 
						  accessibility." (https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html)
- 'gravy' (float) 				: (grand average of hydropathy) calculated by adding the hydropathy value for each residue and dividing by the length of 
						  the sequence (Kyte and Doolittle; 1982) (Bio.SeqUtils.ProtParam documentation). NB : probably redundant with 
						  'mean_kd_hydro' and one should be remooved in the futur.
- 'isoelec_point' (float) 			: the isoelectric point. (Bio.SeqUtils.ProtParam documentation) 
- 'charge_at_pH' (float) 			: calculate the charge of a protein at given pH (default : 7.5).(Bio.SeqUtils.ProtParam documentation). Choice of the pH 
						  value : "The intracellular pH of living cells is strictly controlled in each compartment. Under normal conditions, the 
						  cytoplasmic pH (pHc) and the vacuolar pH (pHv) of typical plant cells are maintained at slightly  alkaline 
						  (typically 7.5) and acidic (typically 5.5) values, respectively (https://doi.org/10.1007/978-3-7091-1254-0_4).
- 'molar_ext_coef' (int) 			: calculates the molar extinction coefficient assuming cysteines (reduced) and cystines residues (Cys-Cys-bond) (Bio.SeqUtils.ProtParam documentation)
- 'tiny' (float) 				: cumulative percentage of the smallest amino acids (e.g 'A', 'C', 'G', 'S', 'T')
- 'small' (float) 				: cumulative percentage of the small amino acids (e.g A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y').

- 'aliphatic' (float) 				: cumulative percentage of the aliphatic amino acids (e.g 'A', 'I', 'L', 'V')
- 'aromatic' (float) 				: cumulative percentage of the aromatic amino acids (e.g 'F', 'H', 'W', 'Y')
- 'non_polar' (float) 				: cumulative percentage of the non-polar amino acids (e.g 'A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y')
- 'polar' (float) 				: cumulative percentage of the polar amino acids (e.g 'D','E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z')
- 'charged' (float) 				: cumulative percentage of the charged amino acids (e.g 'B', 'D', 'E', 'H', 'K', 'R', 'Z')
- 'basic' (float) 				: cumulative percentage of the basic amino acids (e.g 'H', 'K', 'R')
- 'acidic' (float) 				: cumulative percentage of the acidic amino acids (e.g 'B', 'D', 'E', 'Z')
- 'helix' (float) 				:  cumulative percentage of amino acids which tend to be in Helix (e.g. 'V', 'I', 'Y', 'F', 'W', 'L')
- 'turn' (float) 				: cumulative percentage of amino acids which tend to be in turn (e.g 'N', 'P', 'G', 'S')
- 'sheet' (float) 				: cumulative percentage of amino acids which tend to be in turn (e.g  'E', 'M', 'A', 'L')
- 'spp' (str) 					: the species from which the sequence was produced. Not a feature but an easy way to slice the data.
