# Build notes for the released PDF

The released `paper/main.pdf` is the byte output of one build of `paper/main.tex` with the data
availability paragraph carrying the two archival identifiers of this release:

```sh
cd paper
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

Tectonic `0.17.0` (Windows x86_64 release, `xtask -X compile`), TeX Live bundle fetched by
Tectonic. The fixed epoch makes the build reproducible: a second build of the same sources in a
separate empty directory produced identical bytes, SHA-256
`a651bdf1bf846b4db89c5922c7b5e93ab79cd69af16cf5eb71d815ca315d212a`, 265507 bytes, MD5
`4d705db791743743dde78a7db684e266`, 28 pages.

Log statistics of the released build (`main.log`, 19156 bytes):

- exit code `0`, `Output written on main.xdv (28 pages, 233252 bytes)`
- `0` overfull boxes, `0` underfull boxes
- `0` undefined references or citations, `0` LaTeX or package warnings, `0` missing glyphs
- empty bibliography log apart from the plain BibTeX banner (`main.blg`, 185 bytes)

The released PDF carries five external links, extracted from its own link annotations:
`https://doi.org/10.5281/zenodo.24433375`, `https://doi.org/10.5281/zenodo.24433377`,
`https://github.com/SuperMatics-Inc/penson-kolb-holstein-pair-odlro`, and the two cited release
identifiers `https://doi.org/10.5281/zenodo.24428099` and `https://doi.org/10.5281/zenodo.24428101`.

The source of the manuscript is the reviewed text admitted with SHA-256
`2ae3eed3afa6a20f3f2b5144b8a266857c8a181e0dbc92dcd36ca7341adf5782`. The released `paper/main.tex`
differs from it only in the data availability paragraph, which gained the two identifiers and the
repository URL (see `publication-metadata.diff` in the staging directory). Its SHA-256 is
`f1f024f48cec33549b2a653531b194aa8c1e260a84b404c95748c09550ac3ab7`. `paper/refs.bib` is unchanged,
SHA-256 `d856455dd8ce068d6f491f35175b98153d829388b97f718a4c99b528360bc8bd`.
