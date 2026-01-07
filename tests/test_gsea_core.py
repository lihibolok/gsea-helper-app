import pandas as pd

from gsea_core import prepare_ranked_list


def test_prepare_ranked_list_sorts_and_uppercases():
    data = {
        "gene": ["Cd24a", "Gapdh", "Actb"],
        "stat": [1.5, -0.2, 0.8],
    }
    df = pd.DataFrame(data)

    rnk = prepare_ranked_list(df, gene_col="gene", score_col="stat", ascending=False)

    # Check structure
    assert list(rnk.columns) == ["score"]
    assert rnk.index.name == "gene"

    # Check sorting (highest stat first)
    assert rnk.iloc[0]["score"] == 1.5

    # Check uppercase gene symbols
    assert all(isinstance(idx, str) and idx.isupper() for idx in rnk.index)
