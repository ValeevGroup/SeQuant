#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import json
import argparse
import shutil
import re

SILENT_OUTPUT = False # Set to True to suppress command line output

def validate_commit_sha(sha):
    if not re.match(r'^[a-fA-F0-9]{6,40}$', sha):
        raise ValueError(f"Invalid commit SHA format: {sha}")
    return sha

def verify_commit_exists(commit):
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--verify", commit],
            capture_output=True, text=True, check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        raise ValueError(f"Commit does not exist: {commit}")

def run_command(command_list):
    subprocess.run(command_list, check=True, capture_output=SILENT_OUTPUT, text=True)

def get_current_git_state():
    try:
        # try branch name
        result = subprocess.run(
            ["git", "symbolic-ref", "--short", "HEAD"],
            capture_output=True, text=True, check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        # else try commit hash (for detached state)
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True
        )
        return result.stdout.strip()

def restore_git_state(original_git_state):
    print(f"\n==> Restoring git state")
    run_command(["git", "checkout", original_git_state])

# Locate Google Benchmark's compare.py script
def find_google_benchmark_compare(custom_path=None):
    if custom_path:
        if os.path.exists(custom_path):
            return custom_path
        else:
            raise FileNotFoundError(f"Custom compare.py path does not exist: {custom_path}")
    compare = shutil.which("compare.py")
    if compare:
        return compare
    try:
        import google_benchmark
        gb_path = os.path.dirname(google_benchmark.__file__)
        script = os.path.join(gb_path, "tools", "compare.py")
        if os.path.exists(script):
            return script
    except ImportError:
        pass
    raise FileNotFoundError("Google Benchmark's compare.py not found.")


def configure_and_build(commit, cmake_variables, benchmark_target, build_dir):
    print(f"\n==> Building {commit[:8]}")
    run_command(["git", "checkout", commit])
    cmake_cmd = ["cmake", "-S", ".", "-B", build_dir] + cmake_variables
    run_command(cmake_cmd)
    build_cmd = ["cmake", "--build", build_dir, "--target", benchmark_target]
    run_command(build_cmd)

def run_benchmarks(commit, benchmark_target, build_dir):
    print(f"==> Running benchmarks for {commit[:8]}")
    output = f"{commit}-results.json"
    exe_path = os.path.join(build_dir, "benchmarks", benchmark_target)
    if not os.path.isfile(exe_path):
        print(f"Error: Benchmark executable not found at {exe_path}", file=sys.stderr)
        sys.exit(1)

    benchmark_cmd = [
        exe_path,
        "--benchmark_out_format=json",
        "--benchmark_time_unit=us",
        f"--benchmark_out={output}"
    ]
    run_command(benchmark_cmd)

def compare_benchmarks(base_commit, head_commit, benchmark_target, build_dir, compare_path=None, output_file=None):
    base_file = f"{base_commit}-results.json"
    new_file = f"{head_commit}-results.json"

    # Verify that benchmark result files exist
    if not os.path.isfile(base_file):
        print(f"Error: Baseline benchmark results not found: {base_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(new_file):
        print(f"Error: New benchmark results not found: {new_file}", file=sys.stderr)
        sys.exit(1)

    compare_script = find_google_benchmark_compare(compare_path)

    cmd = [
        sys.executable,
        compare_script,
        "benchmarks",
        base_file,
        new_file,
    ]
    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        print("Error running compare.py:", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        print("Command:", ' '.join(cmd), file=sys.stderr)
        sys.exit(1)

    # Use the text output from compare.py directly
    comparison_output = result.stdout

    if not comparison_output.strip():
        print("Error: No benchmark comparison output generated", file=sys.stderr)
        sys.exit(1)

    # Strip ANSI color codes for output file
    ansi_escape_pattern = re.compile(r'\x1b\[[0-9;]*m')
    clean_output = ansi_escape_pattern.sub('', comparison_output)

    print("\n" + clean_output)

    if output_file:
        with open(output_file, "w") as f:
            f.write(clean_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare SeQuant benchmarks between two commits",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("base_commit", help="Base commit SHA to compare against")
    parser.add_argument("head_commit", help="Head commit SHA to compare")
    parser.add_argument("--benchmark-target", "-t",
                        default="sequant_benchmarks",
                        help="Benchmark target to build and run (default: sequant_benchmarks)")
    parser.add_argument("--build-dir", "-b",
                        default="build",
                        help="Build directory for CMake (default: build)")
    parser.add_argument("--compare-path", "-c",
                        default=None,
                        help="Path to Google Benchmark's compare.py script (optional)")
    parser.add_argument("--output-file", "-o",
                        default=None,
                        help="File to write comparison results (optional)")

    args = parser.parse_args()

    # Validate commit SHAs
    try:
        validate_commit_sha(args.base_commit)
        validate_commit_sha(args.head_commit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    original_ref = get_current_git_state()

    # Verify commits exist
    try:
        verify_commit_exists(args.base_commit)
        verify_commit_exists(args.head_commit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Define CMake variables
    cmake_variables = ["-G", "Ninja",
                       "-DCMAKE_BUILD_TYPE=Release",
                       "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache",
                       "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON",
                       "-DSEQUANT_TESTS=OFF",
                       "-DSEQUANT_BENCHMARKS=ON",
                       "-DSEQUANT_MIMALLOC=ON"]

    try:
        # base commit
        configure_and_build(args.base_commit, cmake_variables, args.benchmark_target, args.build_dir)
        run_benchmarks(args.base_commit, args.benchmark_target, args.build_dir)
        # head commit
        configure_and_build(args.head_commit, cmake_variables, args.benchmark_target, args.build_dir)
        run_benchmarks(args.head_commit, args.benchmark_target, args.build_dir)
        compare_benchmarks(
            args.base_commit,
            args.head_commit,
            args.benchmark_target,
            args.build_dir,
            args.compare_path,
            args.output_file
        )
    except Exception as e:
        print(f"Error during benchmark execution: {e}")
        sys.exit(1)
    finally:
        # restore original git state even if script fails
        restore_git_state(original_ref)
