"""levels: bounds, nominal thickness, lookup, and constituent assembly."""
import pytest

import levels
import conftest


def test_low_high():
    lv = conftest.level("0_700")
    assert (lv.low, lv.high) == (0, 700)


@pytest.mark.parametrize("name, thickness", [
    ("0_300", 300),        # 3*(20-15) + 1*(300-15) = 15 + 285
    ("0_700", 700),        # + 1*(700-300) = 400
    ("0_2000", 2000),
    ("700_2000", 1300),    # 1*(1850-700) + 3*(1850-1800) = 1150 + 150
])
def test_nominal_thickness(name, thickness):
    assert conftest.level(name).nominal_thickness == thickness


def test_require_top_is_the_top_300m_of_each_level():
    # require_top is a thickness from the layer's own top; every level requires its top 300 m
    for name in ("0_300", "0_700", "0_1000", "700_2000", "0_2000"):
        assert conftest.level(name).require_top == 300


def test_get_unknown_exits():
    with pytest.raises(SystemExit):
        levels.get("3_4", conftest.table())


def test_shipped_plan_is_the_localgp_table():
    t = conftest.table()
    assert sorted(t) == ["0_1000", "0_2000", "0_300", "0_700", "700_2000"]
    lv = t["0_2000"]
    assert [(c.tag, c.n_fac, c.top, c.bottom) for c in lv.contributors] == [
        ("15_20", 3, 15, 20), ("15_300", 1, 15, 300), ("300_700", 1, 300, 700),
        ("700_1850", 1, 700, 1850), ("1800_1850", 3, 1800, 1850)]
    assert lv.require_top == 300 and not lv.identity


def _plan(tmp_path, body):
    p = tmp_path / "plan.toml"
    p.write_text(body)
    return str(p)


def test_load_rejects_a_plan_that_does_not_tile_its_layer(tmp_path):
    bad = _plan(tmp_path, '''
[[level]]
name = "0_700"
require_top = 300
contributors = [ { tag = "0_300", n_fac = 1, top = 0, bottom = 300 },
                 { tag = "300_600", n_fac = 1, top = 300, bottom = 600 } ]
''')
    with pytest.raises(SystemExit) as e:
        levels.load(bad)
    assert "600 m" in str(e.value) and "700 m" in str(e.value)


def test_load_rejects_unsorted_duplicate_and_bad_require_top(tmp_path):
    unsorted = _plan(tmp_path, '''
[[level]]
name = "0_700"
require_top = 300
contributors = [ { tag = "300_700", n_fac = 1, top = 300, bottom = 700 },
                 { tag = "0_300", n_fac = 1, top = 0, bottom = 300 } ]
''')
    with pytest.raises(SystemExit):
        levels.load(unsorted)
    dup = _plan(tmp_path, '''
[[level]]
name = "0_600"
require_top = 300
contributors = [ { tag = "0_300", n_fac = 1, top = 0, bottom = 300 },
                 { tag = "0_300", n_fac = 1, top = 0, bottom = 300 } ]
''')
    with pytest.raises(SystemExit):
        levels.load(dup)
    deep = _plan(tmp_path, '''
[[level]]
name = "0_300"
require_top = 500
contributors = [ { tag = "0_300", n_fac = 1, top = 0, bottom = 300 } ]
''')
    with pytest.raises(SystemExit):
        levels.load(deep)


def test_load_accepts_another_groups_plan(tmp_path):
    # a group that maps the target layers as 0_700 + 700_2000 directly: n_fac 1 throughout
    theirs = _plan(tmp_path, '''
[[level]]
name = "0_2000"
require_top = 300
contributors = [ { tag = "0_700", n_fac = 1, top = 0, bottom = 700 },
                 { tag = "700_2000", n_fac = 1, top = 700, bottom = 2000 } ]
''')
    t = levels.load(theirs)
    assert t["0_2000"].nominal_thickness == 2000
    assert [c.tag for c in t["0_2000"].contributors] == ["0_700", "700_2000"]


def test_identity_level_is_one_constituent_n_fac_1():
    lv = levels.identity("15_300")
    assert lv.identity and lv.name == "15_300"
    assert (lv.low, lv.high, lv.nominal_thickness, lv.require_top) == (15, 300, 285, 285)
    assert len(lv.contributors) == 1
    c = lv.contributors[0]
    assert (c.tag, c.n_fac, c.top, c.bottom) == ("15_300", 1, 15, 300)
    assert not conftest.level("0_300").identity


def test_resolve_table_level_or_identity_of_a_lone_submission():
    one = {"15_300": {}}
    assert levels.resolve("0_300", one, conftest.table()) is conftest.level("0_300")   # from the plan
    with pytest.raises(SystemExit):                               # --level without a plan: told to pass one
        levels.resolve("0_300", one, None)
    assert levels.resolve(None, one).identity                     # omitted, one submission: identity
    with pytest.raises(SystemExit):                               # omitted, several: nothing to build
        levels.resolve(None, {"15_20": {}, "15_300": {}})
    with pytest.raises(SystemExit):
        levels.resolve(None, {})


def test_constituents_shallowest_first_and_attached():
    lv = conftest.level("0_700")
    subs = {c.tag: {"field_value": conftest.const_field(1.0), "attrs": {}} for c in lv.contributors}
    got = levels.constituents(lv, subs)
    assert [c["tag"] for c in got] == ["15_20", "15_300", "300_700"]
    assert [c["n_fac"] for c in got] == [3, 1, 1]
    assert (got[0]["top"], got[0]["bottom"]) == (15, 20)
    assert got[0]["field_value"] is subs["15_20"]["field_value"]


def test_constituents_missing_submission_exits():
    lv = conftest.level("0_700")
    subs = {"15_20": {"field_value": conftest.const_field(1.0), "attrs": {}}}   # missing 15_300, 300_700
    with pytest.raises(SystemExit):
        levels.constituents(lv, subs)
