# postgresql-dna-extension
PostgreSQL extension implementing custom DNA types (dna, kmer, qkmer) for efficient genomic sequence representation and manipulation.

# TODO

- Fix kmer type (delete hack with allignment = 33)
- Make generate_kmer simpler (look documentation given in project instructions)
- Index support

# Build
## To do once
- install container info-h417-course-image:latest

- path_to_dna_sequence
- container_name

```bash
docker run -it \
  --name container_name \
  -v $(path_to_dna_sequence):/extension/dna_sequence \
info-h417-course-image:latest
```

## To run the extension

- ```docker start -ai container_name```
- ```docker exec -it container_name service postgresql start```
- ```make run```
