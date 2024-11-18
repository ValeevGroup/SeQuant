#!/bin/bash

# Taken from https://stackoverflow.com/a/246128
script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
repo_root="$script_dir/../../"

cd "$repo_root"

git ls-files --cached --modified --others -- \
    ':^SeQuant/external/' ':^build*' '**/*'.cpp '**/*'.c '**/*'.hpp '**/*'.h '**/*'.cxx '**/*'.cc  |
	xargs -r -- "$script_dir/clang-format.sh" --style=file -i
