This folder contains the SCMO NlaIII workflow. 
In this workflow the NlaIII data is demultiplexed, trimmed, mapped, molecules are assigned, mappability is determined per molecule and if an allelic vcf file is provided, the allele is determined per molecule.
The individual steps are described at https://github.com/BuysDB/SingleCellMultiOmics/wiki/NLA-III-data-processing

The SnakeFile contains the actual workflow, and the config.json a default configuration to get you started.
