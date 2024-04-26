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

You can download the dependencies via the conda environment in `environment.yaml`. This method is recommended because it will install Java (via `openjdk`) which is necessary for clusetering with the clusterONE program. Specifically, run the following:

```bash
conda env create -f environment.yaml
```

If you have `mamba` installed you can run the following for faster installation (the output is identical to the command above):

```bash
mamba env create -f environment.yaml
```

Activate the environment with the following command:

```
conda activate rbp_analysis_interactive_tool
```

### Installing Dependencies via Pip

If you don't have `conda` or `mamba` installed, but you do have `openjdk` installed, you can install the remaining dependencies via pip. To do so, run the following command to install the dependencies:

```bash
pip install -r requirements.txt
```

This method may overwrite previously installed versions of packages. To avoid these conflicts, the `conda`/`mamba` method is preferred. Furthermore, this method _does not_ install Java so you will need to install it on your own, use the `conda`/`mamba` method, or avoid use of the `clusterONE` option in the clustering menu.

### Running the Tool

Once you have the dependencies installed, run the application by executing the following command from within `rbp_analysis_interactive_tool`:

```bash
bokeh serve app.py
```

## Citation

Bibtex citation will be listed here once available.
