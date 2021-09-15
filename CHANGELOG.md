# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

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
