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


mathjax3_config = {
    "TeX": {
        "Macros": {
            "bra": r"\\langle",
            "ket": r"\\rangle",
        }
    }
}

# 設定 MathJax 3 的 CDN 路徑
html_theme_options = {
    'mathjax_path': 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js',
}
