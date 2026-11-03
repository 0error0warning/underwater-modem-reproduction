# Sources and third-party material

## Paper and reference table

Jan H. Schmidt, Iwona Kochańska, and Aleksander M. Schmidt. *Performance of low-power underwater modem for shallow water communications*. Vibrations in Physical Systems, 2024, 35(2), 2024207. [DOI](https://doi.org/10.21008/j.0860-6897.2024.2.07) · [Official PDF](https://vibsys.put.poznan.pl/_journal/2024-35-2/articles/vps_2024207.pdf).

The paper states **CC BY 4.0** on its final page. `references/paper_bfsk_reference.csv` is a manual transcription of Tables 2 and 3 (PDF pages 7–8), with renamed columns and an explicit interpretation of active duration and additional guard duration. The three `<0.0001` values remain bounds. This attribution and the [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/) apply to that source-derived table. The original PDF is not distributed here.

## Watermark

Watermark/Mime is provided by the Norwegian Defence Research Establishment (FFI) and its research collaborators. Obtain the official software, manual and measured channels from the [FFI Watermark page](https://www.ffi.no/en/research/watermark). This repository does not distribute their original MATLAB files, channel MAT files, or manual PDFs. Free availability is not represented as a new license grant by this repository.

`matlab/+uac/replayfilter_fast.m` and `matlab/+uac2/replay_at.m` implement optimizations and extensions based on Watermark 1.0 replay conventions and the official implementation. They retain its baseband/resampling semantics and have numerical-equivalence diagnostics. They are not official FFI releases, and this project is not endorsed by FFI or the paper authors.

Recommended Watermark citations:

- Paul van Walree, François-Xavier Socheleau, Roald Otnes, and Trond Jenserud. “The Watermark benchmark for underwater acoustic modulation schemes.” *IEEE Journal of Oceanic Engineering*, 42(4), 1007–1018, 2017. [DOI](https://doi.org/10.1109/JOE.2017.2699078).
- Paul van Walree, Roald Otnes, and Trond Jenserud. *The Watermark manual and user's guide—version 1.0*. FFI Report 16/01378, 2016. [Official publication page](https://www.ffi.no/en/publications-archive/the-watermark-manual-and-users-guide-version-1.0).

## This repository

The code and documentation written for this project are released under the MIT License (see `LICENSE`). The MIT License covers only this project's own code and text. It does not apply to the source paper, its transcribed reference table (CC BY 4.0), Watermark materials, or any other third-party material, which keep their own terms.

