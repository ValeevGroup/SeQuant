use std/assert

const pattern = {

  #
  # A line prefixed by:
  #   - zero or more spaces
  #   - followed by one of [Term, Eval, Cache]
  #   - followed by one or more spaces, and
  #   - followed by the separator character
  # is _potentially_ a trace produced by sequant.
  #
  trace: '^\s*(Term|Eval|Cache)\s+\|' # The sep char '|' is escaped: '\|'.
  separator: '|'
}

# Because of a bug -- that is fixed in future version of nushell
# need to convert to the filesize value smartly.
def "as filesize" [] {
  $in | str trim --right --char 'B' | into float | into filesize
}

def "read log" [file: string] {
    grep -n -E -e $pattern.trace $file
      | parse --regex '(?<line>\d+):\s*(?<trace>.*)'
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
      memory : ($in.1 | as filesize)
      annot  : $in.2
    }
}

def "new-cols cache" [] {
  {
    key    : $in.0
    life   : $in.1
    nalive : ($in.2 | into int)
    memory : ($in.3 | as filesize)
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
  let ends  = $lines | skip 1 | every 2 | append [null]
  $begs | zip $ends
        | each { match $in { [$b, $e] => { assert (($e == null) or ($b.trace.0 == $e.trace.0))
                                           { tag0  : $b.tag0
                                             begin : $b.line
                                             end   : (if ($e == null) { "❎" } else { $e.line })
                                             expr  : $b.trace.0
                                           } } } }
}
