# Pathogen groups

V4 supports five selectable pathogen groups: bacteria, fungi, viruses, parasites, and oomycetes.

The current taxonomy catalogue supplies bacteria, fungi, and viruses. Parasites and oomycetes are added from the locally pinned NCBI Taxonomy taxdump when selected in `00_configure_run.ipynb`. The run records the taxdump checksum and each taxonomic root used. Supplemental taxa receive NCBI-based IDs; existing taxa keep their stable IDs.

Parasites are polyphyletic, so V4 defines that group through selected protozoan and helminth lineages: Apicomplexa (5794), Ciliophora (5878), Platyhelminthes (6157), Nematoda (6231), Amoebozoa (554915), Euglenozoa (33682), Metamonada (2611341), Acanthocephala (10232), Myxozoa (35581), Dicyemida (10215), and Orthonectida (33209). Oomycetes use the NCBI Oomycota root (4762), separate from true fungi.

Species beneath these roots enter the candidate catalogue; pair-specific disease screening determines whether literature supports disease in the selected host. The group dropdown shows entries without profiles in a muted, disabled state for the selected host.
