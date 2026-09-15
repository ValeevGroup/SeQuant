*************************
Evaluation Trace Analysis
*************************

When the evaluator is built with ``SEQUANT_EVAL_TRACE``, it emits one
log line per operation. A trace line is one of three kinds, identified
by its first token:

* ``Eval | <mode> | <time> | [left=L | right=R |] result=X | alloc=A | hw=H | <label>``
* ``Cache | <mode> | key=K | life=c/m | alive=N | entry=E | total=T | <label>``
* ``Term | Begin | <expr>`` / ``Term | End | <expr>``

The trace may be embedded in a larger log; the parsers in this guide
ignore non-trace lines.

Setup
=====

:download:`Save this helper </user/guide/sequant_trace.py>` next to
your trace file. It uses only the Python 3 standard library.

The sample trace below is :download:`from this file
</user/guide/mpqc_out.txt>`, produced by an ``mpqc4`` CCSDT calculation
on a He chain. The same file holds the surrounding mpqc output; the
parsers skip non-trace lines automatically.

.. note::

   The parsers don't know about iterations. If your trace spans more
   than one CC iteration (or any other repeated context), aggregates
   like peak memory or total time blend across them — slice the row
   lists by line range or by term first. The sample here was capped at
   a single iteration to keep the examples clean.

Throughout the rest of this page, assume:

.. code:: python

   import sequant_trace as st
   eval_rows, cache_rows, term_rows = st.parse("mpqc_out.txt")

``parse`` walks the file once and returns three lists. Each list holds
``NamedTuple`` rows (``EvalRow``, ``CacheRow``, ``TermRow``) so columns
are accessed by attribute (``r.mode``, ``r.time_ns``, …). For one-stream
use the generators ``st.trace_eval(path)`` and ``st.trace_cache(path)``,
and the list-returning ``st.trace_term(path)``.

Reading evaluation steps
========================

Each ``EvalRow`` has the fields:

``mode``
    one of ``Tensor``, ``Constant``, ``Variable``, ``Permute``,
    ``MultByPhase``, ``Sum``, ``Product``, ``SumInplace``,
    ``Symmetrize``, ``Antisymmetrize``.
``time_ns``
    wall time in nanoseconds.
``result``
    bytes of the buffer this op produces (for ``SumInplace``, the
    accumulator after the add).
``alloc``
    bytes of fresh allocation by this op. ``0`` for ``SumInplace``;
    equal to ``result`` for every other mode.
``hw``
    live working-set bytes at the moment of the op — cache contents +
    result + each operand not aliased to a cache entry.
``left`` / ``right``
    operand bytes, set only for ``Sum`` and ``Product``.
``label``
    the operation's textual label.
``line``
    1-indexed source line in the trace file.
``extras``
    a dict capturing any ``key=value`` fields the parser didn't
    recognize. Empty today; populated automatically if the evaluator
    starts logging additional columns.

Distinct modes in this trace:

.. code:: python

   sorted({r.mode for r in eval_rows})

::

    ['Constant', 'Permute', 'Product', 'SumInplace', 'Symmetrize', 'Tensor']

First few entries:

.. code:: python

   st.print_table(eval_rows, max_rows=5)

::

    mode    | time_ns   | result    | alloc     | hw        | left      | right | label                                         | line
    --------+-----------+-----------+-----------+-----------+-----------+-------+-----------------------------------------------+-----
    Tensor  | 412107750 | 163879224 | 163879224 | 163879224 |      None |  None | g(i_2,a_1,a_2,a_3)                            |  346
    Tensor  |      1291 |       640 |       640 |       640 |      None |  None | t(a_3,i_2)                                    |  347
    Product |   3307083 |      4992 |      4992 | 163884856 | 163879224 |   640 | g(i_2,a_1,a_2,a_3) * t(a_3,i_2) -> I(a_1,a_2) |  348
    Tensor  |       750 |       640 |       640 |       640 |      None |  None | t(a_2,i_1)                                    |  349
    Product |     68917 |       640 |       640 |      6272 |      4992 |   640 | I(a_1,a_2) * t(a_2,i_1) -> I(a_1,i_1)         |  350

``left`` and ``right`` are populated only for ``Sum`` and ``Product`` —
the only ops whose operand sizes can differ from the result. For other
modes those fields are ``None`` (omitted from the trace), not zero, so
a logged ``0B`` always means an empty buffer.

Top five time-consuming steps:

.. code:: python

   top_t = sorted(eval_rows, key=lambda r: r.time_ns, reverse=True)[:5]
   st.print_table(top_t,
                  columns=["mode","time_ns","result","alloc","hw","label","line"])

::

    mode   | time_ns     | result     | alloc      | hw         | label              | line
    -------+-------------+------------+------------+------------+--------------------+-----
    Tensor | 10352466209 | 3100868920 | 3100868920 | 3100889240 | g(a_1,a_2,a_3,a_4) |  688
    Tensor |   412107750 |  163879224 |  163879224 |  163879224 | g(i_2,a_1,a_2,a_3) |  346
    Tensor |   263996250 |  163879224 |  163879224 |  163900616 | g(a_1,a_2,i_1,a_3) |  797
    Tensor |   237571833 |  163879224 |  163879224 |  163900616 | g(i_3,a_2,a_3,a_4) |  806
    Tensor |   225673167 |  163879224 |  163879224 |  163924632 | g(a_3,a_2,i_3,a_4) | 1829

Top five by working-set peak. ``hw`` is the most useful memory column
for "what was the worst moment of this run"; ``result`` answers "what
is the biggest tensor produced" and ``alloc`` answers "what allocated
the most fresh memory":

.. code:: python

   top_h = sorted(eval_rows, key=lambda r: r.hw, reverse=True)[:5]
   st.print_table(top_h,
                  columns=["mode","time_ns","result","alloc","hw","label","line"])

::

    mode    | time_ns  | result     | alloc      | hw         | label                                                                         | line
    --------+----------+------------+------------+------------+-------------------------------------------------------------------------------+-----
    Product | 72411125 |    1258096 |    1258096 | 3103424296 | g(a_1,a_2,a_4,a_5) * t(a_3,a_4,a_5,i_3,i_1,i_2) -> I(a_3,a_2,a_1,i_2,i_1,i_3) | 1532
    Product | 50984709 |     314800 |     314800 | 3101229768 | g(a_2,a_3,a_4,a_5) * t(a_5,i_3) -> I(a_2,a_3,i_3,a_4)                         | 1840
    Product | 54562167 |     314800 |     314800 | 3101205752 | g(a_1,a_2,a_3,a_4) * t(a_3,i_1) -> I(a_2,a_1,i_1,a_4)                         |  790
    Product | 52750125 |      18864 |      18864 | 3100926968 | g(a_1,a_2,a_3,a_4) * t(a_3,a_4,i_1,i_2) -> I(a_2,a_1,i_2,i_1)                 |  690
    Tensor  |    13875 | 3100868920 | 3100868920 | 3100914328 | g(a_2,a_3,a_4,a_5)                                                            | 1838

Total time spent in the evaluator:

.. code:: python

   total_ns = sum(r.time_ns for r in eval_rows)
   print(f"{total_ns / 1e9:.6f} s")

::

    20.070343 s

Reading cache interactions
==========================

Each ``CacheRow`` has the fields:

``mode``
    one of ``Store``, ``Access``, ``Release``.
``key``
    hash identifying the cached intermediate.
``life_curr`` / ``life_max``
    remaining and total number of accesses for this entry.
``alive``
    number of entries in the cache after this interaction.
``entry``
    bytes occupied by this entry (``0`` for a release).
``total``
    bytes occupied by all live entries combined.
``label``
    label of the node touching the cache.
``line``
    source line in the trace file.
``extras``
    dict for any unrecognized ``key=value`` fields (empty today).

First few entries:

.. code:: python

   st.print_table(cache_rows, max_rows=5)

::

    mode    | key                  | life_curr | life_max | alive | entry | total | label                                         | line
    --------+----------------------+-----------+----------+-------+-------+-------+-----------------------------------------------+-----
    Store   | 13613060045528791584 |         2 |        3 |     2 |   640 |   640 | g(i_2,i_3,a_2,a_3) * t(a_3,i_2) -> I(i_3,a_2) |  397
    Access  | 13613060045528791584 |         1 |        3 |     1 |   640 |   640 | g(i_2,i_3,a_2,a_3) * t(a_3,i_2) -> I(i_3,a_2) |  413
    Store   | 12860816093799665187 |         2 |        3 |     3 |   640 |  1280 | g(i_2,i_3,a_2,a_3) * t(a_3,i_3) -> I(i_2,a_2) |  425
    Access  | 12860816093799665187 |         1 |        3 |     2 |   640 |  1280 | g(i_2,i_3,a_2,a_3) * t(a_2,i_2) -> I(i_3,a_3) |  434
    Release | 13613060045528791584 |         0 |        3 |     1 |     0 |   640 | g(i_2,i_3,a_2,a_3) * t(a_3,i_2) -> I(i_3,a_2) |  445

The first row stores an intermediate with three allowed accesses; the
``Store`` event itself counts as an access, so ``life`` is ``2/3``. The
``total`` column tracks the rolling cache footprint — sum of every live
entry — and the ``entry`` column the size of the single entry this row
refers to. ``Release`` rows show ``entry=0`` because the entry is gone
by the time the row is recorded.

Peak cache footprint over the whole trace:

.. code:: python

   peak = max(cache_rows, key=lambda r: r.total)
   st.print_table([peak],
                  columns=["mode","key","alive","entry","total","line"])

::

    mode  | key               | alive | entry | total | line
    ------+-------------------+-------+-------+-------+-----
    Store | 47978033300314424 |    11 |  4720 | 87856 | 2026

Reading terms
=============

``trace_term`` matches each ``Term | End`` with the most recent
unclosed ``Term | Begin`` and returns one ``TermRow`` per pair:

``begin_line``
    source line of the ``Term | Begin``.
``end_line``
    source line of the ``Term | End``.
``expr``
    expression text shared by the matched pair.

First few:

.. code:: python

   st.print_table(term_rows, max_rows=5)

::

    begin_line | end_line | expr
    -----------+----------+--------------------------------------------------------------------
           345 |      354 | -1 (g{i_2,a_1;a_2,a_3}:N-C-S * t{a_3;i_2}:N-C-S) * t{a_2;i_1}:N-C-S
           355 |      364 | -1 (f{i_2;a_2}:N-C-S * t{a_2;i_1}:N-C-S) * t{a_1;i_2}:N-C-S
           366 |      373 | 2 g{i_2,a_1;a_2,i_1}:N-C-S * t{a_2;i_2}:N-C-S
           375 |      382 | -2 g{i_2,i_3;a_2,a_3}:N-C-S * t{a_1,a_2,a_3;i_3,i_2,i_1}:N-C-S
           384 |      391 | 2 g{i_2,i_3;a_2,a_3}:N-C-S * t{a_1,a_2,a_3;i_1,i_2,i_3}:N-C-S

The line numbers cross-reference into ``eval_rows`` and ``cache_rows``.
For example, every evaluation step inside the first term:

.. code:: python

   first = term_rows[0]
   inside = [r for r in eval_rows
             if first.begin_line <= r.line <= first.end_line]
   st.print_table(inside,
                  columns=["mode","time_ns","result","alloc","hw","label","line"])

::

    mode     | time_ns   | result    | alloc     | hw        | label                                         | line
    ---------+-----------+-----------+-----------+-----------+-----------------------------------------------+-----
    Tensor   | 412107750 | 163879224 | 163879224 | 163879224 | g(i_2,a_1,a_2,a_3)                            |  346
    Tensor   |      1291 |       640 |       640 |       640 | t(a_3,i_2)                                    |  347
    Product  |   3307083 |      4992 |      4992 | 163884856 | g(i_2,a_1,a_2,a_3) * t(a_3,i_2) -> I(a_1,a_2) |  348
    Tensor   |       750 |       640 |       640 |       640 | t(a_2,i_1)                                    |  349
    Product  |     68917 |       640 |       640 |      6272 | I(a_1,a_2) * t(a_2,i_1) -> I(a_1,i_1)         |  350
    Constant |       541 |         8 |         8 |         8 | -1                                            |  351
    Product  |     31625 |       640 |       640 |      1288 | I(a_1,i_1) * -1 -> I(a_1,i_1)                 |  352
    Permute  |     11625 |       640 |       640 |      1280 | I(a_1,i_1)                                    |  353

The reverse lookup also works. Earlier we found that the peak cache
footprint was reached at line 2026; the term responsible is:

.. code:: python

   owner = next(t for t in term_rows
                if t.begin_line <= peak.line <= t.end_line)
   st.print_table([owner], columns=["begin_line","end_line","expr"])

::

    begin_line | end_line | expr
    -----------+----------+--------------------------------------------------------------------------------------------------------
          2022 |     2032 | 6 ((g{i_4,i_5;a_4,a_5}:N-C-S * t{a_5;i_2}:N-C-S) * t{a_1,a_4;i_1,i_3}:N-C-S) * t{a_2,a_3;i_5,i_4}:N-C-S
