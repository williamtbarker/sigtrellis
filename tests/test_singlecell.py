from dataclasses import replace

import anndata as ad
import numpy as np
import pytest

from sigtrellis.qc import encode_outcome, validate_dataset
from sigtrellis.simulate import synthetic_single_cell
from sigtrellis.singlecell import pseudobulk
from sigtrellis.splits import make_splits


def test_chunked_pseudobulk_matches_direct_sum(tmp_path, config):
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=8, genes=20, cells=6)
    local = replace(
        config, cell_type="cell_type", cell_type_value="type_A", min_cells=4, chunk_size=7
    )
    result = pseudobulk(path, local)[0]
    cells = ad.read_h5ad(path)
    for sample in result.expression.index:
        mask = (cells.obs.sample_id == sample) & (cells.obs.cell_type == "type_A")
        expected = np.asarray(cells.X[mask.to_numpy()].sum(axis=0)).ravel()
        np.testing.assert_array_equal(result.expression.loc[sample], expected)
    assert result.expression.shape == (16, 20)
    assert validate_dataset(result, local)["n_biological_groups"] == 8
    y, _ = encode_outcome(result, local)
    for split in make_splits(result, y, local):
        train = set(result.metadata.donor.iloc[split.train])
        test = set(result.metadata.donor.iloc[split.test])
        assert not train & test
        # Every donor contributes both conditions to a single side.
        assert len(split.test) == 2 * len(test)


@pytest.mark.parametrize("dense", [False, True])
def test_count_layer_and_dense_input(tmp_path, config, dense):
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=6, genes=15, cells=4)
    cells = ad.read_h5ad(path)
    cells.layers["counts"] = cells.X.toarray() if dense else cells.X.copy()
    cells.X = cells.X.astype(float) * 0.1
    cells.write_h5ad(path)
    local = replace(
        config, layer="counts", cell_type="cell_type", cell_type_value="type_A", min_cells=3
    )
    result = pseudobulk(path, local)[0]
    assert (result.expression.to_numpy() % 1 == 0).all()
    with pytest.raises(ValueError, match="integer counts"):
        pseudobulk(path, replace(local, layer=None))


def test_many_cells_few_donors_rejected(tmp_path, config):
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=2, genes=20, cells=100)
    local = replace(config, cell_type="cell_type", cell_type_value="type_A")
    data = pseudobulk(path, local)[0]
    with pytest.raises(ValueError, match="four biological groups"):
        validate_dataset(data, local)


def test_do_not_merge_donor_conditions(tmp_path, config):
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=6, genes=20, cells=4)
    with pytest.raises(ValueError, match="Within-sample metadata vary"):
        pseudobulk(path, replace(config, sample_id="donor", group=None, min_cells=3))


def test_min_cells_qc_memory_guard_and_missing_type(tmp_path, config):
    path = tmp_path / "cells.h5ad"
    synthetic_single_cell(path, donors=6, genes=20, cells=4)
    local = replace(config, cell_type="cell_type", cell_type_value="not_there", min_cells=3)
    with pytest.raises(ValueError, match="absent"):
        pseudobulk(path, local)
    with pytest.raises(ValueError, match="fewer than four"):
        pseudobulk(path, replace(local, cell_type_value="type_A", min_cells=100))
    cells = ad.read_h5ad(path)
    cells.obs.iloc[0, cells.obs.columns.get_loc("phenotype")] = None
    cells.write_h5ad(path)
    with pytest.raises(ValueError, match="Missing"):
        pseudobulk(path, replace(local, cell_type_value="type_A"))
