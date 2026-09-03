# Backup Data

Creates versioned backups of important data files.

## Usage

``` r
backup_data(source_path, backup_dir = NULL, max_backups = 5)
```

## Arguments

- source_path:

  Source file path

- backup_dir:

  Backup directory (NULL = same directory)

- max_backups:

  Maximum number of backups to keep

## Value

Backup file path

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API
for callers who need standalone/versioned backups (as opposed to
[`write_soil_data`](https://jjmaynard.github.io/soilSIM/reference/write_soil_data.md)'s
single-backup-on-overwrite behavior), not dead/broken code.
