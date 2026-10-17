# Building the Sphinx documentation

Run the commands below from `doc`, using the Python environment with JetSeT and the documentation dependencies installed. Update the source notebooks under `documentation_notebooks/notebooks` before building.

## Prerequisites

Generate and install the documentation requirements:

```sh
python build_rtd_requirements.py
python -m pip install -r requirements.txt
```

Notebook conversion also requires `jupyter nbconvert` and Pandoc. Cleaning and copying use Python and require no shell utilities.

## Main workflow

`scripts/build_doc.py` coordinates API documentation generation, optional notebook cleaning and conversion, copying generated files, and optional Sphinx HTML generation.

Convert all notebooks using their saved outputs, copy the generated files, and build HTML:

```sh
python scripts/build_doc.py -ncv -b
```

Execute the notebooks before converting and building:

```sh
python scripts/build_doc.py -ncv -e -b
```

To clean generated notebook files first, use `-c`. This also enables conversion:

```sh
python scripts/build_doc.py -c -b
python scripts/build_doc.py -c -e -b
```

Append a directory name to limit notebook operations to `documentation_notebooks/notebooks/dir-name`:

```sh
python scripts/build_doc.py -ncv -b dir-name
```

API documentation is regenerated and the Sphinx build covers the whole documentation even when a notebook directory is selected.

| Option | Behavior |
| --- | --- |
| `-ncv`, `--notebook-convert` | Convert notebooks to `.rst` using saved outputs. |
| `-c`, `--clean` | Delete generated notebook `.rst` and `.png` files in the source and destination, then convert notebooks. |
| `-e`, `--execute` | Execute notebooks in place before conversion; requires `-ncv` or `-c` to enable conversion. |
| `-b`, `--build` | Build Sphinx HTML in `build/`. |
| `-j N`, `--jobs N` | Use N notebook conversion workers; defaults to 1 with execution, otherwise up to 10. |
| `--sphinx-jobs N` | Use N Sphinx workers; defaults to 10. |
| `--clean-only` | Only clean generated notebook files; skip API generation, conversion, copying, and HTML generation. |
| `--copy-only` | Only copy existing generated files and shared images. |
| `dir-name` | Restrict notebook operations to the given subdirectory. |

`--clean-only` and `--copy-only` are standalone actions and cannot be combined with each other or with `-c`, `-ncv`, `-e`, or `-b`. Both accept an optional notebook directory.

With no options, the workflow regenerates API documentation and copies existing generated files. It does **not** convert notebooks or build HTML. `-b` alone also does not enable notebook conversion.

The workflow clears `api/` except for `.gitkeep` before regenerating API documentation. It automatically sets `READTHEDOCS=True` for its Sphinx build to avoid importing the jetkernel module.

## Individual steps

Convert notebooks and copy their generated files without building HTML:

```sh
python scripts/build_doc.py -ncv
# Or limit conversion to one notebook directory:
python scripts/build_doc.py -ncv dir-name
```

These conversion commands also regenerate API documentation. Add `-e` to execute notebooks before conversion.

Clean generated notebook `.rst` and `.png` files without converting:

```sh
python scripts/build_doc.py --clean-only
# Or clean one notebook directory:
python scripts/build_doc.py --clean-only dir-name
```

Copy existing notebook `.rst` and `.png` files to `user_guide/documentation_notebooks_rst/notebooks/`, preserving subdirectories:

```sh
python scripts/build_doc.py --copy-only
# Or copy one notebook directory:
python scripts/build_doc.py --copy-only dir-name
```

The copy action creates the destination directory if necessary and also copies PNGs from `documentation_notebooks/slides` and `documentation_notebooks/images`. It does not generate `.rst` files; conversion must run first.

Regenerate API documentation separately:

```sh
python make_apidoc_and_uml_graphs.py
```

Build HTML directly from existing documentation sources:

```sh
READTHEDOCS=True sphinx-build -j 10 -b html ./ build
```

Use `python scripts/build_doc.py --help` for the available workflow options.
