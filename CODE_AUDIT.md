# STaRparse Code Audit

This document provides a comprehensive review of the `STaRparse` repository codebase, focusing on code quality, best practices, and modularity. The goal of this audit is to identify areas for improvement before making actual code changes.

## 1. Code Quality

### Hardcoded Paths & Credentials
**Issue:** Throughout `STaRparse.py` (and scripts inside `Source/`), there are numerous hardcoded absolute paths to specific home directories (e.g., `/home/dannear/STaRparse/...` or `/home/dannear/annovar/...`).
- **Impact:** This makes the code completely unusable on any other machine or user account.
- **Recommendation:** Use relative paths, environment variables, or config files. Use `os.path.join(os.path.dirname(__file__), ...)` to reference files relative to the script's location.

**Issue:** Database connection strings are hardcoded in the `to_DB` function (e.g., `mysql+pymysql://dannear:3vVrnhJ3@143.169.238.18/dannear`).
- **Impact:** This is a major security risk (exposing credentials) and prevents using the tool with other databases.
- **Recommendation:** Remove credentials from the codebase. Use environment variables (e.g., `os.environ.get("DB_URI")`), or securely pass them via a configuration file or command-line arguments.

### Potential Bugs
**Issue:** Improper use of `.strip()` in string manipulation.
- In `STaRparse.py`: `base.strip('.bam.GangSTR.vcf')`
- **Impact:** The `.strip()` method removes any character in the provided string from the beginning and end of the target string. It does *not* remove the exact substring. This can unintentionally chop off trailing or leading characters of a valid filename (e.g., a sample named `sample_v`).
- **Recommendation:** Instead of `.strip()`, use string replacement like `base.replace('.bam.GangSTR.vcf', '')` or use `.removesuffix('.bam.GangSTR.vcf')` (Python 3.9+).

**Issue:** Incomplete Database Export Handling
- The `args.database` value is sometimes `None`, which concatenates strangely in DB operations or avoids database operations, but error reporting could be better formatted.

## 2. Best Practices

### PEP 8 Adherence
**Issue:** The codebase does not adhere to standard Python styling (PEP 8).
- **Imports:** Multiple modules are imported on a single line (`import pandas as pd, argparse, glob, os, vcf, datetime, sqlalchemy, shlex, subprocess, random`).
- **Naming Conventions:** Mixed case variable and function names (e.g., `ExpansionHunter()`, `mdf`, `out_code`, `to_DB`). Functions and variables should be `lowercase_with_underscores`. Classes should be `CapWords`.
- **Recommendation:** Refactor code to follow standard PEP 8 naming conventions. Run formatting tools like `black` or `flake8` over the repository.

### DRY (Don't Repeat Yourself) Principle
**Issue:** The `ExpansionHunter` and `GangSTR` functions in `STaRparse.py` contain highly duplicated code for VCF parsing and constructing the dataframes.
- **Recommendation:** Extract the common VCF reading and parsing logic into a helper function (e.g., `parse_vcf_record(record, genotyper_type)`) to significantly reduce code length, reduce bugs, and make future updates simpler.

### Error Handling & Edge Cases
**Issue:** Shell commands are invoked using `subprocess.Popen` without checking for success or failure correctly.
- Example: `p.stdout.read()` is called, but no exception handling (like `try/except`) surrounds it, nor does it check return codes.
- **Recommendation:** Replace `subprocess.Popen(...)` followed by `p.stdout.read()` with `subprocess.run(..., check=True, capture_output=True, text=True)`. This allows explicit error catching if the R script or perl script fails.

## 3. Modularity

### Command-Line Structure
**Issue:** `STaRparse_Stable.py` implements subcommands (like `vcfparse`, `csvmerge`) using basic attribute checks on `sys.argv` and dynamic path insertion (`sys.path.insert()`).
- **Recommendation:** Instead of manually handling `sys.argv[1:2]` and dynamically injecting module paths based on hardcoded directories, use the standard `argparse` Subparsers feature (`parser.add_subparsers()`). This creates a cleaner, self-documenting CLI, and handles `--help` properly across all commands.

### File Organization
**Issue:** The scripts have nested imports like `from From_ExpansionHunter import ExpansionHunter` directly inside the CLI execution block.
- **Recommendation:** Organize the application into a proper Python package structure. For example, have a `starparse/` folder containing modules like `cli.py`, `vcf_parser.py`, `db.py`, `annovar.py`, and `r_scripts.py`.

## Summary Explanation

The changes recommended above are designed to make the `STaRparse` tool:
1. **Portable:** By removing hardcoded paths, anyone can clone and run it.
2. **Secure:** By removing hardcoded database passwords, we protect the project's infrastructure.
3. **Robust:** By fixing bugs like `.strip()` and improving `subprocess` usage, the code will fail gracefully and process data correctly.
4. **Maintainable:** Applying PEP 8, reducing code duplication (DRY), and using built-in `argparse` subparsers ensures that the project is easier to read, test, and expand in the future.
