## Notes on Package Size

The repository currently contains large vignette cache files under vignettes/articles/figures.
These files are used for reproducible article outputs and are not required for core package functionality.
If needed for CRAN submission, these assets can be excluded from the built source package via .Rbuildignore or replaced with smaller artifacts.
