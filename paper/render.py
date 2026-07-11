"""Render the manuscript to DOCX (and PDF if a converter is available).

Usage:  python paper/render.py
Requires: pypandoc-binary (bundles pandoc). PDF additionally needs either
Microsoft Word (via docx2pdf) or a LaTeX engine.
"""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
MD = HERE / "protein_function_rescue.md"
DOCX = HERE / "protein_function_rescue.docx"
PDF = HERE / "protein_function_rescue.pdf"
TITLE = ("Structure-based functional annotation of hypothetical proteins "
         "using the AlphaFold Database and TM-align")


def render_docx():
    import pypandoc
    import os
    cwd = os.getcwd()
    os.chdir(HERE)  # so relative image paths (figures/, ../benchmark/) resolve
    try:
        pypandoc.convert_file(
            MD.name, "docx", outputfile=DOCX.name,
            extra_args=["--resource-path=.:..", "--toc", "--toc-depth=2",
                        f"--metadata=title:{TITLE}"],
        )
    finally:
        os.chdir(cwd)
    print(f"wrote {DOCX}")


def render_pdf():
    try:
        from docx2pdf import convert
        convert(str(DOCX), str(PDF))
        if PDF.exists():
            print(f"wrote {PDF}")
            return True
    except Exception as exc:
        print(f"PDF via Word failed ({exc})")
    # Fallback: pandoc with a LaTeX engine, if present.
    try:
        import shutil, pypandoc
        for engine in ("xelatex", "pdflatex", "tectonic"):
            if shutil.which(engine):
                pypandoc.convert_file(str(MD), "pdf", outputfile=str(PDF),
                                      extra_args=[f"--pdf-engine={engine}",
                                                  "--resource-path=" + str(HERE) + ":" + str(HERE.parent)])
                print(f"wrote {PDF} (via {engine})")
                return True
    except Exception as exc:
        print(f"PDF via LaTeX failed ({exc})")
    print("PDF not produced (install MS Word or a LaTeX engine); DOCX is available.")
    return False


if __name__ == "__main__":
    render_docx()
    render_pdf()
