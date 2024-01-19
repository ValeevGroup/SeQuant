#!/usr/bin/env python3

from __future__ import print_function
import sys
import os
import os.path as op

def to_base_version(full_version):
    base_version = full_version
    plus_pos = base_version.find('+')
    if plus_pos != -1:
        base_version = base_version[0:plus_pos]
    minus_pos = base_version.find('-')
    if minus_pos != -1:
        base_version = base_version[0:minus_pos]
    return base_version

def dots_to_undescores(str):
    return str.replace('.', '_')

def escape_special_chars(str):
    #str = str.replace('(', '\(')
    #str = str.replace(')', '\)')
    str = str.replace('/', '\/')
    str = str.replace('.', '\.')
    return str

def replace_dep_id(topsrc, file_ext, dep_name, old_id, new_id, search_prefix = '', search_postfix = ''):
    any_files_changed = False
    if old_id != new_id:
        # always exclude the versions file
        seek_retcode = os.system('grep -q -r --include="*.' + file_ext + '" --exclude="' + topsrc + '/external/versions.cmake" "' + search_prefix + old_id + search_postfix + '" ' + topsrc)
        if os.WIFEXITED(seek_retcode) and os.WEXITSTATUS(seek_retcode) == 0:
            any_files_changed = True
            print('changing ' + dep_name + ' id from', old_id, 'to', new_id)
            esc_search_prefix = escape_special_chars(search_prefix)
            esc_search_postfix = escape_special_chars(search_postfix)
            os.system('find ' + topsrc + ' -type f -name "*.' + file_ext + '" -print0 | xargs -0 sed -i \'\' -e \'s/' + esc_search_prefix + old_id + esc_search_postfix + '/' + esc_search_prefix + new_id + esc_search_postfix + '/g\'')
    return any_files_changed

argv = sys.argv
topsrc = op.normpath(op.join(op.abspath(op.dirname(sys.argv[0])), '../..'))
if len(argv) == 1:
    version_cmake_path = topsrc + '/external/versions.cmake'
elif len(argv) == 2:
    # no-op if given
    version_cmake_path = op.abspath(sys.argv[1])
    if op.basename(version_cmake_path) != 'versions.cmake':
        sys.exit(0)
else:
    print('invalid number of arguments')
    sys.exit(0)

# extract dependencies tags and versions
with open(version_cmake_path) as inf:
    for line in inf:
        line = line.replace('(', ' ')
        line = line.replace(')', ' ')
        tokens = line.split()
        if len(tokens) < 3:
            continue
        if tokens[1].find('TRACKED_TILEDARRAY') != -1:
            if tokens[1].find('PREVIOUS') != -1:
                ta_old_tag = tokens[2]
            else:
                ta_new_tag = tokens[2]
        elif tokens[1].find('RANGEV3') != -1:
            if tokens[1].find('PREVIOUS') != -1:
                rangev3_old_tag = tokens[2]
            else:
                rangev3_new_tag = tokens[2]

any_files_changed = False

# replace TA tag in INSTALL.md
any_files_changed |= replace_dep_id(topsrc, 'md', 'TiledArray', ta_old_tag, ta_new_tag, '', '')

# Range-v3 tag in INSTALL.md
any_files_changed |= replace_dep_id(topsrc, 'md', 'Range-V3', rangev3_old_tag, rangev3_new_tag, '', '')

if any_files_changed:
    sys.exit(1)
else:
    sys.exit(0)
