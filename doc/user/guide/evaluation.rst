Interpreting and Evaluating Expressions
=======================================

Evaluation Trace Analysis Using Nushell
---------------------------------------

`Nushell <https://www.nushell.sh>`__ can be installed from ``Homebrew``.

::

   brew install nushell

Save this file as ``sequant.nu`` in the same directory where the trace
output file is present. We will call this directory 'working directory'.

.. code:: nu


   use std/assert

   const pattern = {
     trace: '^\s*(Term|Eval|Cache)\s+'
     separator: '|'
   }

   def "read log" [file: string] {
       grep -n -E -e $pattern.trace $file
         | parse --regex '(?<line>\d+):\s+(?<trace>.*)'
         | each {
           let line = $in.line
           let fields = $in.trace
                        | split row $pattern.separator 
                        | str trim
           match $fields {
             [$tag0 $tag1 ..$rest] => {
                         {tag0 : $tag0
                          tag1 : $tag1
                          trace: $rest
                          line : ($line | into int) } }
             }
           }
   }

   def "new-cols eval" [] {
       {
         time   : ($in.0 | into duration)
         memory : ($in.1 | into filesize)
         annot  : $in.2
       }
   }

   def "new-cols cache" [] {
     {
       key    : $in.0
       life   : $in.1
       nalive : ($in.2 | into int)
       memory : ($in.3 | into filesize )
     }
   }

   def make-tag [tag: string, make_col: closure] {
     let entry = $in | where tag0 == $tag
     let new_cols = $entry.trace | each { do $make_col }
     $entry | select tag0 tag1
            | merge $new_cols
            | merge ($entry | select line)
   }

   export def "trace cache" [file: string] {
     read log $file | make-tag 'Cache' { new-cols cache }
   }

   export def "trace eval" [
               file: string # The text file with sequant evaluation trace.
               ] {
     read log $file | make-tag 'Eval' { new-cols eval }
   }

   export def "trace term" [file: string] {
     let lines = read log $file | where tag0 == 'Term'
     let begs  = $lines | every 2
     let ends  = $lines | skip 1 | every 2
     $begs | zip $ends
           | each { match $in { [$b, $e] => { assert ($b.trace.0 == $e.trace.0)
                                              { tag0  : $b.tag0
                                                begin : $b.line
                                                end   : $e.line
                                                expr  : $b.trace.0
                                              } } } }
   }

This nu script exports the functions: ``trace eval``, ``trace cache``,
and ``trace term``.

If you are using a different shell (eg. bash, zsh, etc.), launch the
nushell from the working directory. Make sure the executable ``nu``
installed from Homebrew is on your path.

.. code:: sh

   nu

Then import the functions from the script.

.. code:: nu

   use sequant.nu *

Now we analyze the output in steps.

.. code:: nu

   let eval = trace eval out.txt

The ``$eval`` variable is a table that contains trace from evaluation
steps.

.. code:: nu

   $eval | describe
   # table<tag0: string, tag1: string, time: duration, memory: filesize, annot: string, line: int>

   # The tag0 column has the same value 'Eval'

   # Let's look at the tag1 column:

   $eval | get tag1 | uniq
   # ╭───┬────────────────╮
   # │ 0 │ Tensor         │
   # │ 1 │ Permute        │
   # │ 2 │ Product        │
   # │ 3 │ SumInplace     │
   # │ 4 │ Constant       │
   # │ 5 │ MultByPhase    │
   # │ 6 │ Antisymmetrize │
   # ╰───┴────────────────╯

   # First ten entries (without the annot column)
   $eval | first 10 | reject annot

   # ╭───┬──────┬────────────┬──────────────────┬──────────┬──────╮
   # │ # │ tag0 │    tag1    │       time       │  memory  │ line │
   # ├───┼──────┼────────────┼──────────────────┼──────────┼──────┤
   # │ 0 │ Eval │ Tensor     │  9ms 634µs 666ns │   2.2 kB │  305 │
   # │ 1 │ Eval │ Permute    │         9µs 84ns │   4.4 kB │  306 │
   # │ 2 │ Eval │ Tensor     │ 10ms 137µs 417ns │   4.0 kB │  310 │
   # │ 3 │ Eval │ Tensor     │            500ns │   4.0 kB │  311 │
   # │ 4 │ Eval │ Product    │       87µs 792ns │  12.1 kB │  312 │
   # │ 5 │ Eval │ Permute    │       12µs 375ns │   8.0 kB │  313 │
   # │ 6 │ Eval │ SumInplace │       34µs 709ns │   8.0 kB │  315 │
   # │ 7 │ Eval │ Tensor     │        4µs 458ns │   2.2 kB │  318 │
   # │ 8 │ Eval │ Tensor     │            208ns │ 152.6 kB │  319 │
   # │ 9 │ Eval │ Product    │      124µs 916ns │ 158.9 kB │  320 │
   # ╰───┴──────┴────────────┴──────────────────┴──────────┴──────╯

   # Top five most time consuming evaluations
   $eval | sort-by time --reverse | first 5 | select annot time
   # ╭───┬────────────────────┬──────────────────╮
   # │ # │       annot        │       time       │
   # ├───┼────────────────────┼──────────────────┤
   # │ 0 │ f(i_2,i_1)         │ 10ms 630µs 500ns │
   # │ 1 │ f(a_1,a_2)         │ 10ms 137µs 417ns │
   # │ 2 │ f(a_1,i_1)         │  9ms 634µs 666ns │
   # │ 3 │ g(i_2,a_1,a_2,a_3) │  2ms 313µs 375ns │
   # │ 4 │ g(a_1,a_2,a_3,a_4) │  1ms 205µs 542ns │
   # ╰───┴────────────────────┴──────────────────╯

   # Total time spent in evaluation
   $eval | get time | math sum
   # 74ms 604µs 696ns
   $eval | get time | math sum | format duration ms
   # 74.60 ms
