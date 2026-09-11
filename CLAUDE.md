# CLAUDE.md

This repository keeps its agent-facing conventions in `AGENTS.md`, which is the
name most agent tooling looks for. Claude Code reads `CLAUDE.md`, so this file
imports it — keep the content in `AGENTS.md`, not here.

This is deliberately a regular file rather than a symlink to `AGENTS.md`: a
committed symlink checks out on Windows as a regular file containing the
literal text `AGENTS.md` unless `core.symlinks` is enabled, which silently
yields an empty-looking config rather than an error. See #605.

@AGENTS.md
