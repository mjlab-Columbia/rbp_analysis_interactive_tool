# RBP Analysis Interactive Tool

## Getting Started

### Downloading the Code

Clone the repository via the following command:

```bash
git clone https://github.com/mjlab-Columbia/rbp_analysis_interactive_tool.git
```

Navigate to the `rbp_analysis_interactive_tool` directory. If you're on MacOS or Linux, enter the following into the terminal:

```bash
cd rbp_analysis_interactive_tool
```

### Installing Dependencies via Conda/Mamba

You can download all the dependencies via the conda environment in `environment.yaml`. Specifically, run the following:

```bash
conda env create -f environment.yaml
```

If you have `mamba` installed you can run the following for faster installation:

```bash
mamba env create -f environment.yaml
```

Activate the environment with the following command:

```
conda activate rbp_analysis_interactive_tool
```

### Installing Dependencies via Pip

If you don't have `conda` or `mamba` installed, you can install the dependencies via pip. Run the following command to install the dependencies:

```bash
pip install -r requirements.txt
```

Note that this may overwrite versions of packages you already have installed. To avoid these conflicts, the `conda`/`mamba` method is preferred.

### Running the Tool

Run the application by executing the following command from within `rbp_analysis_interactive_tool`:

```bash
bokeh serve app.py
```

## Citation

Bibtex citation will be listed here once available.
