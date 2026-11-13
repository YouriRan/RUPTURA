"""Read fixtures written by test_npz; NumPy is a test-only dependency."""
from pathlib import Path
import sys
import zipfile
import numpy as np
root = Path(sys.argv[1] if len(sys.argv)>1 else "npz_test_output")
with np.load(root/"schema.npz",allow_pickle=False) as a:
    np.testing.assert_array_equal(a["time"],[0,.25,1])
    np.testing.assert_array_equal(a["state"].T,a["y"])
    assert a["state"].shape==(3,2) and a["column_c"].shape==(2,2,2)
    np.testing.assert_array_equal(a["column_c"].ravel(),np.arange(1,9))
    assert a["cycle"].dtype==np.dtype("<i8")
    np.testing.assert_array_equal(a["phase"],[0,1,2])
    np.testing.assert_array_equal(a["comp_names"],["CO₂","H₂","CH₄","氦"])
    np.testing.assert_array_equal(a["unicode"],["α","😀","","aé氦"])
    assert a["unicode"].dtype==np.dtype("<U3")
    assert a["termination_reason"].shape==() and a["termination_reason"].item()=="pressure criterion"
    assert a["empty_text"].item()=="" and a["scalar_double"].item()==3.25 and a["scalar_int"].item()==-2
    assert a["bools"].dtype==np.dtype("bool")
    np.testing.assert_array_equal(a["bools"],[True,False,True])
    np.testing.assert_array_equal(a["mixed_number"],[1,2.5,0])
    np.testing.assert_array_equal(a["signed"],[-1,-2**63,2**63-1])
    np.testing.assert_array_equal(a["unsigned"],np.array([0,2**64-1],dtype=np.uint64))
    assert a["empty"].shape==(0,) and a["empty_nested"].shape==(2,0)
with np.load(root/"explicit.npz",allow_pickle=False) as a:
    assert a["explicit_empty"].shape==(0,3,7)
    assert np.isnan(a["special"][0]) and np.isposinf(a["special"][1]) and np.isneginf(a["special"][2])
    np.testing.assert_array_equal(a["tensor"].ravel(),np.arange(1,9))
with np.load(root/"zero_entries.npz",allow_pickle=False) as a:
    assert not a.files
for p in root.glob("*.npz"):
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        assert all(i.compress_type==(zipfile.ZIP_DEFLATED if p.name=="compressed.npz" else zipfile.ZIP_STORED) for i in z.infolist())
assert not list(root.glob("*.tmp-*"))
print("NumPy dtype, shape, C-order, Unicode and ZIP CRC checks passed")

with np.load(root/"schema.npz", allow_pickle=False) as a, np.load(root/"compressed.npz", allow_pickle=False) as b:
    assert a.files == b.files
    for key in a.files:
        np.testing.assert_array_equal(a[key], b[key])
print("Compressed/uncompressed archives agree exactly in NumPy")
