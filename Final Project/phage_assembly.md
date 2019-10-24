# Phage Genome Assembly

## Comparing Assembly Tools

Genomes to use
- Ask Sian/Natalia for assembled phage genomes from nanopore. Test on them.
    - Or, take a complete phage genome and simulate long read sequences
        - https://github.com/rrwick/Badread/tree/master/comparison
        - ***[DeepSimulator](https://academic.oup.com/bioinformatics/article/34/17/2899/4962495)
- Simulate annoying bad reads to test
    - https://github.com/rrwick/Badread


Tools to use
- Canu
- Flye
- Unicycler

[rrwick - long read assembler comparison](https://github.com/rrwick/Long-read-assembler-comparison#assessing-plasmid-assembly)
- No Ra because need to be able to circularize (in case phage genome is circular)
- No Wtdbg2 (worst performing in the comparison for plasmids)

### Similarities
- What are the similarities in terms of algorithm?

### Differences
- What are the differences in terms of algorithm?


### Performance
[rrwick - long read assembler comparison; see results section for categories](https://github.com/rrwick/Long-read-assembler-comparison#assessing-plasmid-assembly)



## Modifying Best Performing Assembler
- Now, choosing the best performing assembler, can we do some pre-processing or multi-step approach to further optimize the tool for phage genome assembly?