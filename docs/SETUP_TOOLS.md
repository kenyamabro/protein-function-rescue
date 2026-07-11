# Installing the external tools (Foldseek + fpocket) — OPTIONAL

> **You do not need these to run the project.** The default backends
> (`structural_search.backend: tmalign` and `pocket.backend: geometry`) are pure
> Python and run on Windows/macOS/Linux with just `pip install -r requirements.txt`.
>
> Install Foldseek/fpocket only if you want the faster Linux tools — Foldseek in
> particular scales to searching *millions* of targets, which TM-align can't. After
> installing, set `structural_search.backend: foldseek` and/or `pocket.backend: fpocket`
> in `config.yaml`.

| Stage  | Optional tool | Enable with              |
|--------|---------------|--------------------------|
| search | Foldseek      | `structural_search.backend: foldseek` |
| pocket | fpocket       | `pocket.backend: fpocket` |

Both are native Linux/macOS tools with **no Windows build**. On Windows use WSL2
(Ubuntu) or conda; instructions for all three below.

---

## Option A — conda (works on Windows, macOS, Linux)

```bash
conda create -n rescue -c conda-forge -c bioconda python=3.12 foldseek fpocket
conda activate rescue
pip install -r requirements.txt
```

`foldseek` and `fpocket` are then on your PATH — the defaults in `config.yaml`
(`binary: "foldseek"` / `"fpocket"`) just work.

## Option B — WSL2 on Windows

```powershell
wsl --install -d Ubuntu      # once, then reboot
```

Inside the Ubuntu shell:

```bash
sudo apt update && sudo apt install -y fpocket
# Foldseek static build:
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xzf foldseek-linux-avx2.tar.gz
export PATH="$PWD/foldseek/bin:$PATH"   # add to ~/.bashrc to persist
```

Run the pipeline from inside WSL, pointing at the repo (your Windows drives are
under `/mnt/c/...`).

## Option C — macOS / Linux directly

```bash
# Foldseek
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz   # or -osx-universal
tar xzf foldseek-*.tar.gz && export PATH="$PWD/foldseek/bin:$PATH"
# fpocket
sudo apt install fpocket           # Debian/Ubuntu
# or: brew install brewsci/bio/fpocket   (macOS, via homebrew-bio tap)
```

---

## Build the Foldseek target database (once)

The target database must contain proteins of **known** function. Good choices:

```bash
# Experimental structures from the PDB (recommended starting point, a few GB):
foldseek databases PDB foldseek_db/pdb tmp

# Or Swiss-Prot AlphaFold models (curated, well-annotated):
foldseek databases Alphafold/Swiss-Prot foldseek_db/swissprot tmp
```

Point `config.yaml → foldseek.target_db` at whichever you built.

## Verifying

```bash
foldseek version
fpocket -h | head -n 3
```

If a binary isn't on PATH, put its absolute path in `config.yaml`
(`foldseek.binary` / `pocket.binary`).
