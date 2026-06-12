"""Parsers for SeQuant evaluation traces.

A trace is a stream of pipe-delimited lines emitted by the evaluator;
each line is one of three kinds, identified by its first token:

    Eval  | <mode> | <time> | [left=L | right=R |] result=X | alloc=A | hw=H | <label>
    Cache | <mode> | key=K  | life=c/m | alive=N | entry=E | total=T | <label>
    Term  | Begin  | <expr>
    Term  | End    | <expr>

The trace may be embedded in a larger log; lines that don't start with
one of the three tags are skipped. The parsers concatenate every trace
event in the file — if the log covers more than one iteration (or any
other repeated context), the caller is responsible for slicing rows by
line range or term before computing aggregates.

Public API:
    parse(path)        -> (eval_rows, cache_rows, term_rows)   one-pass
    trace_eval(path)   -> generator of EvalRow
    trace_cache(path)  -> generator of CacheRow
    trace_term(path)   -> list of TermRow
    format_table(rows, columns=..., max_rows=...) -> str
    print_table(...)   -> prints format_table output

Row types are NamedTuples. Unknown key=value fields land in row.extras
so a future evaluator column survives without a parser change.
"""

from __future__ import annotations

import re
from typing import NamedTuple


_LINE_TAGS = frozenset({"Eval", "Cache", "Term"})
_KV_RE = re.compile(r"^([A-Za-z_]+)=(.+)$")

_EVAL_KNOWN_KEYS = frozenset({"result", "alloc", "hw", "left", "right"})
_CACHE_KNOWN_KEYS = frozenset({"key", "life", "alive", "entry", "total"})


# ---------- row schemas ----------

class EvalRow(NamedTuple):
    mode: str
    time_ns: int
    result: int
    alloc: int
    hw: int
    left: "int | None"
    right: "int | None"
    label: str
    line: int
    extras: dict


class CacheRow(NamedTuple):
    mode: str
    key: int
    life_curr: int
    life_max: int
    alive: int
    entry: int
    total: int
    label: str
    line: int
    extras: dict


class TermRow(NamedTuple):
    begin_line: int
    end_line: int
    expr: str


# Columns displayed by default; excludes `extras` to keep tables readable.
_EVAL_DEFAULT_COLUMNS = tuple(c for c in EvalRow._fields if c != "extras")
_CACHE_DEFAULT_COLUMNS = tuple(c for c in CacheRow._fields if c != "extras")
_TERM_DEFAULT_COLUMNS = TermRow._fields


# ---------- low-level parsing helpers ----------

def _strip_suffix(s, suffix):
    """'13129B' with suffix='B' -> 13129. Returns int."""
    assert s.endswith(suffix), s
    return int(s[: -len(suffix)])


def _parse_kv(field):
    """Parse 'key=value' into (key, value); (None, None) if not k=v."""
    m = _KV_RE.match(field)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def _kv_and_label(fields):
    """Split a list of trailing fields into (kv: dict, label: str | None).

    Fields shaped 'k=v' go into the dict (values left as raw strings).
    The first non-kv token, if any, becomes the label.
    """
    kv, label = {}, None
    for f in fields:
        k, v = _parse_kv(f)
        if k is None:
            label = f
        else:
            kv[k] = v
    return kv, label


# ---------- per-row parsers ----------

def _parse_eval_row(lineno, fields):
    """fields: [<mode>, <time>, <k=v>..., <label>]"""
    mode = fields[0]
    time_ns = _strip_suffix(fields[1], "ns")
    kv, label = _kv_and_label(fields[2:])
    extras = {k: v for k, v in kv.items() if k not in _EVAL_KNOWN_KEYS}
    return EvalRow(
        mode=mode,
        time_ns=time_ns,
        result=_strip_suffix(kv["result"], "B"),
        alloc=_strip_suffix(kv["alloc"], "B"),
        hw=_strip_suffix(kv["hw"], "B"),
        left=_strip_suffix(kv["left"], "B") if "left" in kv else None,
        right=_strip_suffix(kv["right"], "B") if "right" in kv else None,
        label=label,
        line=lineno,
        extras=extras,
    )


def _parse_cache_row(lineno, fields):
    """fields: [<mode>, <k=v>..., <label>]"""
    mode = fields[0]
    kv, label = _kv_and_label(fields[1:])
    life_curr, life_max = (int(x) for x in kv["life"].split("/"))
    extras = {k: v for k, v in kv.items() if k not in _CACHE_KNOWN_KEYS}
    return CacheRow(
        mode=mode,
        key=int(kv["key"]),
        life_curr=life_curr,
        life_max=life_max,
        alive=int(kv["alive"]),
        entry=_strip_suffix(kv["entry"], "B"),
        total=_strip_suffix(kv["total"], "B"),
        label=label,
        line=lineno,
        extras=extras,
    )


class _TermPairer:
    """Pair Term | Begin / Term | End lines under a LIFO invariant.

    The evaluator never interleaves terms, so the matching End for
    every Begin is the most recently opened one. Violations raise.
    """

    def __init__(self):
        self._stack = []
        self.pairs = []

    def handle(self, lineno, fields):
        # fields: [<kind>, <expr>] where kind is 'Begin' or 'End'.
        kind, expr = fields[0], fields[1]
        if kind == "Begin":
            self._stack.append((lineno, expr))
            return
        assert kind == "End", kind
        assert self._stack, f"Term | End at line {lineno} with no open Begin"
        begin_line, begin_expr = self._stack.pop()
        assert begin_expr == expr, (
            f"Term mismatch at line {lineno}: "
            f"End {expr!r} closing Begin {begin_expr!r}"
        )
        self.pairs.append(TermRow(begin_line, lineno, expr))


# ---------- file-level iteration ----------

def _iter_trace_lines(path):
    """Yield (lineno, tag, rest_fields) for trace lines in `path`.

    rest_fields excludes the tag token. Lines that don't start with a
    known tag (after whitespace strip) are silently skipped.
    """
    with open(path, encoding="utf-8") as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            head = line.split(" | ", 1)[0]
            if head not in _LINE_TAGS:
                continue
            fields = line.split(" | ")
            yield lineno, fields[0], fields[1:]


# ---------- public API ----------

def parse(path):
    """Parse all three line kinds in a single pass.

    Returns (eval_rows, cache_rows, term_rows) — each a list,
    materialized for downstream querying.
    """
    eval_rows, cache_rows = [], []
    term_pairer = _TermPairer()
    handlers = {
        "Eval": lambda ln, fs: eval_rows.append(_parse_eval_row(ln, fs)),
        "Cache": lambda ln, fs: cache_rows.append(_parse_cache_row(ln, fs)),
        "Term": term_pairer.handle,
    }
    for lineno, tag, fields in _iter_trace_lines(path):
        handlers[tag](lineno, fields)
    return eval_rows, cache_rows, term_pairer.pairs


def trace_eval(path):
    """Yield an EvalRow for every Eval line in `path`."""
    for lineno, tag, fields in _iter_trace_lines(path):
        if tag == "Eval":
            yield _parse_eval_row(lineno, fields)


def trace_cache(path):
    """Yield a CacheRow for every Cache line in `path`."""
    for lineno, tag, fields in _iter_trace_lines(path):
        if tag == "Cache":
            yield _parse_cache_row(lineno, fields)


def trace_term(path):
    """Return paired TermRow list for `path`."""
    pairer = _TermPairer()
    for lineno, tag, fields in _iter_trace_lines(path):
        if tag == "Term":
            pairer.handle(lineno, fields)
    return pairer.pairs


# ---------- presentation ----------

def _row_field(row, name):
    return getattr(row, name) if hasattr(row, "_fields") else row.get(name)


def _default_columns(row):
    if isinstance(row, EvalRow):
        return _EVAL_DEFAULT_COLUMNS
    if isinstance(row, CacheRow):
        return _CACHE_DEFAULT_COLUMNS
    if isinstance(row, TermRow):
        return _TERM_DEFAULT_COLUMNS
    if hasattr(row, "_fields"):
        return row._fields
    return tuple(row.keys())


def _is_int_column(cells, col_idx):
    """A column is numeric if at least one non-None cell parses as int."""
    saw_value = False
    for row in cells:
        v = row[col_idx]
        if v == "None":
            continue
        saw_value = True
        try:
            int(v)
        except ValueError:
            return False
    return saw_value


def format_table(rows, columns=None, max_rows=None):
    """Render `rows` as an aligned text table and return the string.

    Numeric columns (those whose non-None cells all parse as int) are
    right-aligned; all other columns are left-aligned.
    """
    rows = list(rows)
    if not rows:
        return "(empty)"
    if columns is None:
        columns = _default_columns(rows[0])
    if max_rows is not None:
        rows = rows[:max_rows]

    cells = [[str(_row_field(r, c)) for c in columns] for r in rows]
    widths = [max(len(c), *(len(row[i]) for row in cells))
              for i, c in enumerate(columns)]
    aligners = [
        str.rjust if _is_int_column(cells, i) else str.ljust
        for i in range(len(columns))
    ]

    out = []
    out.append(" | ".join(c.ljust(widths[i]) for i, c in enumerate(columns)))
    out.append("-+-".join("-" * w for w in widths))
    for row in cells:
        out.append(" | ".join(aligners[i](row[i], widths[i])
                              for i in range(len(columns))))
    return "\n".join(out)


def print_table(rows, columns=None, max_rows=None):
    """Print `rows` as an aligned text table."""
    print(format_table(rows, columns, max_rows))
