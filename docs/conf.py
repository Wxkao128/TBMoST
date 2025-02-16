# -- Basic Info -----------------------------------------------------
project = 'TBMoST'
author = 'wxkao128'
release = '1.0'

# -- Main -----------------------------------------------------------
master_doc = 'TBMoST_docs_tc'  # (index.rst or index.md)

extensions = [
    'sphinx.ext.mathjax',  # 保留 MathJax
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

latex_elements = {
    'preamble': r'\usepackage{braket}',  # 加載 braket 宏包
}

html_theme_options = {
    'mathjax_path': 'https://cdn.jsdelivr.net/npm/mathjax@2/es5/tex-mml-chtml.js',
}

mathjax_config = {
    "TeX": {
        "Macros": {
            "bra": r"\langle",  # 定義 bra 符號
            "ket": r"\rangle",  # 定義 ket 符號
        }
    }
}
