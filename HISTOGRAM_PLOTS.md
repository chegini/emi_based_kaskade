# Stage histogram plots

`plot_stage_histograms.py` reads the histogram CSV files produced by a
compression run. It compares only the two configurations used in the main
study: no lifting and lifting after quantization.

The script generates three files for each model:

- `*-stage-raw.png`: raw quantized symbols before bit-length coding.
- `*-stage-zigzag.png`: symbols after zigzag encoding.
- `*-stage-bitlength.png`: bit-length symbols passed to Huffman coding when
  bit-length coding is enabled.

The raw and zigzag plots show the 30 most probable symbols for readability.
The CSV files remain the complete data set.

## EMI

```bash
python3 plot_stage_histograms.py \
  --root output/session13/SX_h16 \
  --output-prefix output/figures/emi \
  --title "EMI BDDC"
```

This writes `output/figures/emi-stage-raw.png`,
`output/figures/emi-stage-zigzag.png`, and
`output/figures/emi-stage-bitlength.png`.

## Poisson

```bash
python3 plot_stage_histograms.py \
  --root output/poisson-study \
  --output-prefix output/figures/poisson \
  --title "Poisson BDDC"
```

This writes `output/figures/poisson-stage-raw.png`,
`output/figures/poisson-stage-zigzag.png`, and
`output/figures/poisson-stage-bitlength.png`.

The script requires Python 3 and Matplotlib. It does not rerun the solver.

## Python environment

Create the environment once in the project directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install matplotlib
```

Activate it again whenever you open a new shell:

```bash
source .venv/bin/activate
```

The `.venv` directory is machine-specific and must not be committed. Only the
plotting script and this instruction file belong in Git.
