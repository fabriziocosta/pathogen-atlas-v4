# Pathogen groups

V4 has five selectable pathogen groups: bacteria, fungi, viruses, parasites, and oomycetes.

The canonical registry of each group's taxonomy source, source URL, release date or manifest, and member-selection rule is the [pathogen taxonomy registry](../configs/pathogen_taxonomies.yaml). The Explorer builder reads source provenance from that file, and the NCBI supplement loader reads the root taxa from it.

Bacteria, fungi, and viruses use the pinned staged taxonomy catalogue, with their source records identified in the registry. Parasites and oomycetes use species beneath selected roots in a locally pinned NCBI Taxonomy taxdump. Each run records the taxdump checksum and roots used. Parasites are defined by multiple protozoan and helminth lineages because that group is polyphyletic; oomycetes use NCBI's Oomycota root, separate from true fungi.

The Explorer may show a group with no host-specific profile as unavailable for that host. Membership in a pathogen group does not itself establish disease evidence for a selected host.
