#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import json
import argparse

SILENT_OUTPUT = True # Set to True to suppress command line output

def run_command(command):
    subprocess.run(command, shell=True, check=True, capture_output=SILENT_OUTPUT, text=True)

# replace slashes in branch names with dashes for file naming, otherwise it will create directories
def process_branch_name(branch_name):
    return branch_name.replace('/', '-')

def get_current_git_state():
    try:
        # try branch name
        result = subprocess.run("git symbolic-ref --short HEAD", shell=True,
                                capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        # else try commit hash (for detached state)
        result = subprocess.run("git rev-parse HEAD", shell=True, check=True, capture_output=True, text=True)
        return result.stdout.strip()

def restore_git_state(original_git_state):
    print(f"\nRestoring original git state: {original_git_state}")
    run_command(f"git checkout {original_git_state}")

def configure_and_build(commit, cmake_variables, benchmark_target, build_dir):
    print(f"\nConfiguring and building commit: {commit}\n")

    run_command(f"git checkout {commit}")

    cmake_vars_str = " ".join(cmake_variables)
    command = f"cmake -S . -B {build_dir} -DCMAKE_BUILD_TYPE=Release {cmake_vars_str}"
    run_command(command)

    command = f"cmake --build {build_dir} --target {benchmark_target} --clean-first"
    run_command(command)

def run_benchmarks(commit, benchmark_target, build_dir):
    print(f"Running benchmarks for commit: {commit}\n")
    output = process_branch_name(commit) + "-results.json"
    command = f"./{build_dir}/benchmarks/{benchmark_target} --benchmark_out_format=json --benchmark_time_unit=us --benchmark_out={output}"
    run_command(command)
    print(f"Benchmarks for commit {commit} completed and results saved to {output}\n")

def compare_benchmarks(base_benchmark, new_benchmark, metric="cpu_time"):
    if metric not in ["cpu_time", "real_time"]:
        raise ValueError("Invalid metric specified. Use 'cpu_time' or 'real_time'.")

    # Validate files exist and load JSON data
    try:
        if not os.path.exists(base_benchmark):
            raise FileNotFoundError(f"Base file not found: {base_benchmark}")
        if not os.path.exists(new_benchmark):
            raise FileNotFoundError(f"New file not found: {new_benchmark}")

        with open(base_benchmark, 'r') as f:
            base_data = json.load(f)
        with open(new_benchmark, 'r') as f:
            new_data = json.load(f)

    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON format: {e}")
    except Exception as e:
        raise RuntimeError(f"Error reading files: {e}")

    output_file = "benchmark-comparison.txt"
    print(f"\nComparing benchmarks between {base_benchmark} and {new_benchmark} using metric: {metric}\n")

    # index benchmarks by name
    base_benchmarks = {b['name']: b for b in base_data['benchmarks']}
    new_benchmarks = {b['name']: b for b in new_data['benchmarks']}

    # get time unit
    time_unit = base_data['benchmarks'][0].get('time_unit', 'ns')

    # only compare benchmarks that are present in both files
    common_names = []
    new_benchmark_names = set(new_benchmarks.keys())
    for benchmark in base_data['benchmarks']:
        if benchmark['name'] in new_benchmark_names:
            common_names.append(benchmark['name'])

    # output details
    output_lines = [f"Benchmark Comparison: {base_benchmark} vs {new_benchmark}", f"Metric: {metric}",
                    f"Time Unit: {time_unit}", f"Date: {os.popen('date').read().strip()}", ""]

    if 'context' in base_data:
        output_lines.append("Base Benchmark Context:")
        output_lines.append(json.dumps(base_data['context'], indent=2))
        output_lines.append("")

    if 'context' in new_data:
        output_lines.append("New Benchmark Context:")
        output_lines.append(json.dumps(new_data['context'], indent=2))
        output_lines.append("")

    output_lines.append(f"{'Name':<60} {f'Base ({time_unit})':<15} {f'New ({time_unit})':<15} {f'Diff ({time_unit})':<15} {'% Diff':<10}")
    output_lines.append("-" * 125)

    # compare each benchmark
    for name in common_names:
        base_value = base_benchmarks[name][metric]
        new_value = new_benchmarks[name][metric]

        diff = new_value - base_value
        percentage_diff = ((new_value - base_value) / base_value) * 100.0 if base_value != 0 else 0.0
        sign = "+" if percentage_diff > 0 else ""

        line = f"{name:<60} {base_value:<15.2f} {new_value:<15.2f} {diff:<15.2f} {sign}{percentage_diff:<9.2f}"
        output_lines.append(line)

    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines))

    print(f"\nComparison results written to: {output_file}")


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

    args = parser.parse_args()

    # print info
    print("**" * 50)
    print("SeQuant Benchmark Comparison Script\n")
    print(f"Base commit: {args.base_commit}")
    print(f"Head commit: {args.head_commit}")
    print(f"Benchmark target: {args.benchmark_target}")
    print(f"Build directory: {args.build_dir}")
    print("**" * 50)

    original_ref = get_current_git_state()
    print(f"Original git reference: {original_ref}")

    # Define CMake variables
    cmake_variables = ["-G Ninja",
                       "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON",
                       "-DSEQUANT_TESTS=OFF",
                       "-DSEQUANT_EVAL_TESTS=OFF",
                       "-DSEQUANT_BENCHMARKS=ON",
                       "-DSEQUANT_MIMALLOC=ON",
                       "-DSEQUANT_CONTEXT_MANIPULATION_THREADSAFE=ON"]

    print("\nCMake variables:")
    for var in cmake_variables:
        print(f"{var}")

    try:
        # base commit
        configure_and_build(args.base_commit, cmake_variables, args.benchmark_target, args.build_dir)
        run_benchmarks(args.base_commit, args.benchmark_target, args.build_dir)

        # head commit
        configure_and_build(args.head_commit, cmake_variables, args.benchmark_target, args.build_dir)
        run_benchmarks(args.head_commit, args.benchmark_target, args.build_dir)

        # Compare benchmarks
        base_benchmark_file = process_branch_name(args.base_commit) + "-results.json"
        new_benchmark_file = process_branch_name(args.head_commit) + "-results.json"
        compare_benchmarks(base_benchmark_file, new_benchmark_file, metric="cpu_time")
        print("Benchmark comparison completed successfully.")

    except Exception as e:
        print(f"Error during benchmark execution: {e}")
        sys.exit(1)
    finally:
        # restore original git state even if script fails
        restore_git_state(original_ref)
