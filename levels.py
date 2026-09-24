"""The synthetic-level plan: the combined-level table, and the identity level.

A synthetic level is a weighted sum of native ME4OH levels ("constituents"), shallowest first. `n_fac`
scales a thin measured layer up to the slab it stands in for; `top`/`bottom` are the constituent's own
dbar bounds, used against the standard bathy for the fully-wet / seafloor / dry test in step 2.

A level is a plan for combining constituents, and one constituent needs no plan: with no `--level`
and exactly one submission in the pool, `resolve` builds the **identity level** — that submission's
native tag as the level, one contributor with `n_fac = 1`. This is how a single-layer quantity (a
mixed layer depth, say) goes through the factory: mask, primitives, collapse, and a combine that
returns the constituent. Its `top`/`bottom` come from the tag and may be nominal.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Contributor:
    tag: str        # native-level tag "top_bottom", matches the submission's mapped_layer attr
    n_fac: int      # thin-layer multiplier
    top: int        # dbar, shallow edge
    bottom: int     # dbar, deep edge


@dataclass(frozen=True)
class Level:
    name: str
    contributors: tuple
    require_top: int    # metres of the layer's own top (from `low`) that must be defined for a cell to
                        # survive under contiguous_from_top; the top 300 m of every layer
    identity: bool = False   # built from a single submission (see `resolve`), not from the table

    @property
    def low(self):
        return int(self.name.split("_")[0])

    @property
    def high(self):
        return int(self.name.split("_")[1])

    @property
    def nominal_thickness(self):
        """Nominal layer thickness in metres: n_fac-weighted sum of the constituent thicknesses."""
        return sum(c.n_fac * (c.bottom - c.top) for c in self.contributors)


LEVELS = [
    Level("0_300",  (Contributor("15_20", 3, 15, 20), Contributor("15_300", 1, 15, 300)), 300),
    Level("0_700",  (Contributor("15_20", 3, 15, 20), Contributor("15_300", 1, 15, 300),
                     Contributor("300_700", 1, 300, 700)), 300),
    Level("0_1000", (Contributor("15_20", 3, 15, 20), Contributor("15_300", 1, 15, 300),
                     Contributor("300_700", 1, 300, 700), Contributor("700_1000", 1, 700, 1000)), 300),
    Level("700_2000", (Contributor("700_1850", 1, 700, 1850), Contributor("1800_1850", 3, 1800, 1850)), 300),
    Level("0_2000", (Contributor("15_20", 3, 15, 20), Contributor("15_300", 1, 15, 300),
                     Contributor("300_700", 1, 300, 700), Contributor("700_1850", 1, 700, 1850),
                     Contributor("1800_1850", 3, 1800, 1850)), 300),
]

_BY_NAME = {lv.name: lv for lv in LEVELS}


def get(name):
    """The Level for `name`."""
    if name not in _BY_NAME:
        raise SystemExit("unknown synthetic level %r; known: %s" % (name, list(_BY_NAME)))
    return _BY_NAME[name]


def identity(tag):
    """The identity level for one native tag `top_bottom`: that tag as the level, one contributor,
    `n_fac = 1`, `require_top` = its whole thickness."""
    top, bottom = (int(x) for x in tag.split("_"))
    return Level(tag, (Contributor(tag, 1, top, bottom),), bottom - top, identity=True)


def resolve(name, submissions):
    """`--level` given -> the table level. Not given -> the identity level of the one submission in
    the pool; with any other number of submissions there is nothing to build, and that's an error."""
    if name is not None:
        return get(name)
    tags = sorted(submissions)
    if len(tags) != 1:
        raise SystemExit("--level is required unless the pool holds exactly one submission (the identity "
                         "level); got %d: %s" % (len(tags), tags))
    return identity(tags[0])


def constituents(level, submissions):
    """Attach each contributor's loaded submission to its plan entry, shallowest first.

    Returns a list of dicts {tag, n_fac, top, bottom, field_value}, where `field_value` is the
    (realization, time, lat, lon) stack for that native level. Errors if a needed submission wasn't
    provided.
    """
    out = []
    for c in level.contributors:
        if c.tag not in submissions:
            raise SystemExit("level %s needs native layer %s, which wasn't among the submissions"
                             % (level.name, c.tag))
        out.append({"tag": c.tag, "n_fac": c.n_fac, "top": c.top, "bottom": c.bottom,
                    "field_value": submissions[c.tag]["field_value"]})
    return out
