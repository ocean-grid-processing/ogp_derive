"""The synthetic-level plan: a table of combined levels read from a file, and the identity level.

A synthetic level is a weighted sum of native ME4OH levels ("constituents"), shallowest first. `n_fac`
scales a thin measured layer up to the slab it stands in for; `top`/`bottom` are the constituent's own
dbar bounds, used against the standard bathy for the fully-wet / seafloor / dry test in step 2. The
table comes from a TOML plan (`--levels`; LocalGP's is levels/localgp.toml), so a group whose native
layers decompose the target layers differently writes its own plan and no code changes. `load`
validates each level: the plan must tile its layer.

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


def load(path):
    """Read a level plan (TOML, `[[level]]` tables — see levels/localgp.toml) -> {name: Level}.

    Each level is checked as it is read, because a plan that doesn't tile its layer would combine
    silently into the wrong number: contributors listed shallowest first (tops non-decreasing), unique
    tags, the n_fac-weighted thicknesses summing to the level's thickness (high - low), and
    `require_top` positive and no more than that thickness.
    """
    try:
        import tomllib
    except ModuleNotFoundError:                                 # Python < 3.11
        import tomli as tomllib
    with open(path, "rb") as f:
        doc = tomllib.load(f)
    table = {}
    for entry in doc.get("level", []):
        try:
            name = entry["name"]
            contributors = tuple(Contributor(c["tag"], int(c["n_fac"]), int(c["top"]), int(c["bottom"]))
                                 for c in entry["contributors"])
            level = Level(name, contributors, int(entry["require_top"]))
        except (KeyError, TypeError, ValueError) as e:
            raise SystemExit("%s: malformed level entry %r (%s)" % (path, entry.get("name"), e))
        _validate(level, path)
        if name in table:
            raise SystemExit("%s: level %r is defined twice" % (path, name))
        table[name] = level
    if not table:
        raise SystemExit("%s: no [[level]] entries" % path)
    return table


def _validate(level, path):
    where = "%s: level %s" % (path, level.name)
    try:
        low, high = level.low, level.high
    except (ValueError, IndexError):
        raise SystemExit("%s: name must be `<low>_<high>` in metres" % where)
    if not level.contributors:
        raise SystemExit("%s: no contributors" % where)
    tags = [c.tag for c in level.contributors]
    if len(set(tags)) != len(tags):
        raise SystemExit("%s: a contributor tag appears twice (%s)" % (where, tags))
    tops = [c.top for c in level.contributors]
    if tops != sorted(tops):
        raise SystemExit("%s: contributors must be listed shallowest first (tops %s)" % (where, tops))
    for c in level.contributors:
        if c.bottom <= c.top or c.n_fac < 1:
            raise SystemExit("%s: contributor %s needs bottom > top and n_fac >= 1" % (where, c.tag))
    if level.nominal_thickness != high - low:
        raise SystemExit("%s: the n_fac-weighted constituent thicknesses sum to %d m, but the level is %d m "
                         "thick; the plan must tile its layer"
                         % (where, level.nominal_thickness, high - low))
    if not 0 < level.require_top <= high - low:
        raise SystemExit("%s: require_top must be in (0, %d]" % (where, high - low))


def get(name, table):
    """The Level for `name` in a loaded plan."""
    if name not in table:
        raise SystemExit("unknown synthetic level %r; the plan defines: %s" % (name, sorted(table)))
    return table[name]


def identity(tag):
    """The identity level for one native tag `top_bottom`: that tag as the level, one contributor,
    `n_fac = 1`, `require_top` = its whole thickness."""
    top, bottom = (int(x) for x in tag.split("_"))
    return Level(tag, (Contributor(tag, 1, top, bottom),), bottom - top, identity=True)


def resolve(name, submissions, table=None):
    """`--level` given -> that level from the loaded plan (`table`, required). Not given -> the identity
    level of the one submission in the pool; with any other number of submissions there is nothing to
    build, and that's an error."""
    if name is not None:
        if table is None:
            raise SystemExit("--level %s names a synthetic level, which needs a plan: pass --levels <file> "
                             "(the LocalGP plan is levels/localgp.toml)" % name)
        return get(name, table)
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
