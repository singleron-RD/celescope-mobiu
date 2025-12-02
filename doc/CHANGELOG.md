# Changelog

## [1.4.1] – 2025-12-02
- Changed: Use multiprocess in the `convert` step to save time.

## [1.4.0] – 2025-11-20

- Added: Implemented automatic chemistry detection, enabling the pipeline to infer chemistry types without manual configuration.

- Changed: Removed unused chemistry profiles mobiu-2 and mobiu-3.

- Changed: Renamed mobiu-4 to mobiu-2.

## [1.3.0] – 2025-07-25

- Added: Support for multiple input FASTQ files per sample.

- Added: Output of splice junction (SJ) matrices in a MARVEL-compatible format.

## [1.2.0] – 2025-06-26

- Added: Support for the mobiu-4 chemistry.

- Fixed: Corrected the reverse complement logic for umi_qual.

## [1.1.0] – 2025-04-17

- Added: Transcript marker functionality.
