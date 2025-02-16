# -- Basic Info -----------------------------------------------------
project = 'TBMoST'
author = 'wxkao128'
release = '1.0'

# -- Main -----------------------------------------------------------
master_doc = 'TBMoST_docs_tc'  # (index.rst or index.md)

extensions = [
    'sphinx.ext.imgmath',
    'sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

html_theme_options = {
    'mathjax_path': 'https://cdn.jsdelivr.net/npm/mathjax@2/es5/tex-mml-chtml.js',
}

mathjax_config = {
    "TeX": {
        "Macros": {
            "bra": r"\\langle",
            "ket": r"\\rangle",
        }
    }
}


#language = 'zh'
