# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [1.3] - 2021-11-30
### Added
- Added loop modelling routines: `routines.Routine_loopmodel_default` and `Routine_loopmodel_slow`, build on the base class `routines.Routine_loopmodel`
- Added functionality for, given a target sequence, automatically searching for 
potential templates and generating alignments using different methods:
 - Base class: `alignment.AlignmentGenerator`
 - Classes for application: `alignment.AlignmentGenerator_pdb`, `alignment.AlignmentGenerator_hhblits`, `alignment.AlignmentGenerator_from_aln`
- Added Tutorial 8 about automatic generation of alignments.
- Methods for checking coverage between sequences in alignments:
 - `alignment.Alignment.calc_coverage`
 - `alignment.Alignment.calc_coverage_target`
 - `alignment.Alignment.calc_pairwise_coverage_all`
- `pdb_io` submodule for pooling functionality that is internally used to handle PDB files.
- Added .gitignore file.

### Fixed
- Fixed mistake in docstring for `alignment.Alignment.calc_identity()`.

### Changed
- `alignment.Alignment.calc_identity_target()` and
`alignment.Alignment.calc_pairwise_identity_all()` do not calculate
identities between sequences of the same identifier anymore.

## [1.2] - 2021-09-15
### Added
- Started using a CHANGELOG.
- Added support for PDF on readthedocs page.

### Changed
- Improved documentation for PyPI.
- Updated references in tutorials.

### Fixed
- Fixed a typo in `docker/template_container/Dockerfile_homelette_template`.
- Fixed the suggested command to build the template container.
- Fixed access to extension submodule.
- Fixed routine tags for complex modelling routines.

## [1.1] - 2021-09-08
### Added
- Added `docs/requirements.txt` as requirements for building the 
documentation.

### Fixed
- Fixed missing optional dependency check that made `modeller` and
`altmod` hard dependencies.

### Changed
- Improved communication to user when running homelette in a Docker
container in Juptyer Lab mode.
- In `docker/template_container/Dockerfile_homelette_template`, changed
the way the homelette package was included from the local build context
to a GitHub tag.

## [1.0] - 2021-09-07
### Added
- Initial release of homelette

[Unreleased]: https://github.com/philippjunk/homelette/compare/1.2...HEAD
[1.2]: https://github.com/philippjunk/homelette/compare/1.1...1.2
[1.1]: https://github.com/philippjunk/homelette/compare/1.0...1.1
[1.0]: https://github.com/philippjunk/homelette/releases/tag/1.0
