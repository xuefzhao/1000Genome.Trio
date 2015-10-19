#Index hg38 using gem tools
```
gem-indexer -i hg38.fa -o hg38 
gem-mappability -I hg38.gem -l 50 -o hg38.50mer -T 8
gem-2-wig -I hg38.gem -i hg38.50mer.mappability -o hg38.50mer
```
