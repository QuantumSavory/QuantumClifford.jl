#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "usage: $0 OUTPUT_DIR LABEL=CHECKOUT [LABEL=CHECKOUT ...]" >&2
    echo "environment: PRECOMPILE_BENCHMARK_BUILDS, PRECOMPILE_BENCHMARK_SAMPLES, PRECOMPILE_BENCHMARK_BASELINES, PRECOMPILE_BENCHMARK_CONSUMER_PROJECT, PRECOMPILE_BENCHMARK_CONSUMER_MANIFEST, PRECOMPILE_BENCHMARK_ALLOW_DIRTY, PRECOMPILE_BENCHMARK_ALLOW_JULIA_MISMATCH" >&2
    exit 2
}

rewrite_checkout_manifest_paths() {
    local input_path=$1
    local output_path=$2
    local old_quantumclifford_path=$3
    local new_quantumclifford_path=$4
    local old_qeccore_path=$5
    local new_qeccore_path=$6
    if ! awk \
        -v old_quantumclifford_path="$old_quantumclifford_path" \
        -v new_quantumclifford_path="$new_quantumclifford_path" \
        -v old_qeccore_path="$old_qeccore_path" \
        -v new_qeccore_path="$new_qeccore_path" '
        BEGIN {
            quantumclifford_header = "[[deps.QuantumClifford]]"
            qeccore_header = "[[deps.QECCore]]"
            old_quantumclifford_line = "path = \"" old_quantumclifford_path "\""
            new_quantumclifford_line = "path = \"" new_quantumclifford_path "\""
            old_qeccore_line = "path = \"" old_qeccore_path "\""
            new_qeccore_line = "path = \"" new_qeccore_path "\""
        }
        {
            if ($0 == quantumclifford_header) {
                quantumclifford_stanzas += 1
                in_quantumclifford = 1
                in_qeccore = 0
            } else if ($0 == qeccore_header) {
                qeccore_stanzas += 1
                in_quantumclifford = 0
                in_qeccore = 1
            } else if ($0 ~ /^\[\[deps\./) {
                in_quantumclifford = 0
                in_qeccore = 0
            }
            if ($0 == old_quantumclifford_line) {
                quantumclifford_path_lines += 1
                in_quantumclifford || misplaced_quantumclifford_path = 1
                print new_quantumclifford_line
            } else if ($0 == old_qeccore_line) {
                qeccore_path_lines += 1
                in_qeccore || misplaced_qeccore_path = 1
                print new_qeccore_line
            } else {
                print
            }
        }
        END {
            if (quantumclifford_stanzas != 1 || qeccore_stanzas != 1 ||
                quantumclifford_path_lines != 1 || qeccore_path_lines != 1 ||
                misplaced_quantumclifford_path || misplaced_qeccore_path)
                exit 42
        }
    ' "$input_path" > "$output_path"; then
        rm -f -- "$output_path"
        return 1
    fi
}

validate_consumer_project() {
    local project_path=$1
    cmp -s -- <(printf '%s\n' \
        '[deps]' \
        'QuantumClifford = "0525e862-1e90-11e9-3e4d-1b39d7109de1"') "$project_path"
}

[[ $# -ge 3 ]] || usage

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
harness_root=$(cd -- "$script_dir/../.." && pwd -P)
if ! harness_commit=$(git -C "$harness_root" rev-parse HEAD 2>/dev/null); then
    echo "the cold-start harness must run from a committed Git checkout" >&2
    exit 1
fi
harness_run_sha256=$(sha256sum "$script_dir/run.sh" | awk '{print $1}')
harness_scenarios_sha256=$(sha256sum "$script_dir/scenarios.jl" | awk '{print $1}')
harness_summarize_sha256=$(sha256sum "$script_dir/summarize.jl" | awk '{print $1}')

verify_harness_checkout() {
    local current_commit relative_path current_sha256 committed_sha256
    current_commit=$(git -C "$harness_root" rev-parse HEAD 2>/dev/null) || {
        echo "failed to read the cold-start harness commit" >&2
        exit 1
    }
    [[ $current_commit == "$harness_commit" ]] || {
        echo "cold-start harness HEAD changed during measurement" >&2
        exit 1
    }
    for relative_path in \
        benchmark/precompile/run.sh \
        benchmark/precompile/scenarios.jl \
        benchmark/precompile/summarize.jl; do
        current_sha256=$(sha256sum "$harness_root/$relative_path" | awk '{print $1}')
        committed_sha256=$(git -C "$harness_root" show "$harness_commit:$relative_path" | sha256sum | awk '{print $1}')
        [[ $current_sha256 == "$committed_sha256" ]] || {
            echo "cold-start harness file differs from $harness_commit: $relative_path" >&2
            exit 1
        }
    done
}

verify_harness_checkout
output_dir=$1
shift
mkdir -p -- "$output_dir"
output_dir=$(cd -- "$output_dir" && pwd -P)

raw_path="$output_dir/raw.tsv"
summary_path="$output_dir/summary.tsv"
build_summary_path="$output_dir/build-summary.tsv"
markdown_path="$output_dir/summary.md"
metadata_path="$output_dir/metadata.txt"
consumer_project_path="$output_dir/consumer-Project.toml"
consumer_manifest_path="$output_dir/consumer-Manifest.toml"
for path in "$raw_path" "$summary_path" "$build_summary_path" "$markdown_path" "$metadata_path" "$consumer_project_path" "$consumer_manifest_path"; do
    [[ ! -e "$path" && ! -L "$path" ]] || {
        echo "refusing to overwrite existing result file: $path" >&2
        exit 2
    }
done

builds=${PRECOMPILE_BENCHMARK_BUILDS:-1}
samples=${PRECOMPILE_BENCHMARK_SAMPLES:-2}
baseline_map_list=${PRECOMPILE_BENCHMARK_BASELINES:-}
reuse_consumer_project=${PRECOMPILE_BENCHMARK_CONSUMER_PROJECT:-}
reuse_consumer_manifest=${PRECOMPILE_BENCHMARK_CONSUMER_MANIFEST:-}
allow_dirty=${PRECOMPILE_BENCHMARK_ALLOW_DIRTY:-0}
allow_julia_mismatch=${PRECOMPILE_BENCHMARK_ALLOW_JULIA_MISMATCH:-0}
julia=${JULIA:-julia}
expected_julia='julia version 1.12.6'
token_pattern='^[A-Za-z0-9][A-Za-z0-9_.-]*$'
mapping_list_pattern='^[A-Za-z0-9][A-Za-z0-9_.-]*=[A-Za-z0-9][A-Za-z0-9_.-]*(,[A-Za-z0-9][A-Za-z0-9_.-]*=[A-Za-z0-9][A-Za-z0-9_.-]*)*$'

[[ $builds =~ ^[1-9][0-9]*$ ]] || { echo "PRECOMPILE_BENCHMARK_BUILDS must be a positive integer" >&2; exit 2; }
[[ $samples =~ ^[1-9][0-9]*$ ]] || { echo "PRECOMPILE_BENCHMARK_SAMPLES must be a positive integer" >&2; exit 2; }
[[ $allow_dirty == 0 || $allow_dirty == 1 ]] || { echo "PRECOMPILE_BENCHMARK_ALLOW_DIRTY must be 0 or 1" >&2; exit 2; }
[[ $allow_julia_mismatch == 0 || $allow_julia_mismatch == 1 ]] || { echo "PRECOMPILE_BENCHMARK_ALLOW_JULIA_MISMATCH must be 0 or 1" >&2; exit 2; }
[[ -z $baseline_map_list || $baseline_map_list =~ $mapping_list_pattern ]] || {
    echo "PRECOMPILE_BENCHMARK_BASELINES must be a comma-separated CANDIDATE=BASELINE list" >&2
    exit 2
}
if [[ -n $reuse_consumer_project || -n $reuse_consumer_manifest ]]; then
    [[ -n $reuse_consumer_project && -n $reuse_consumer_manifest ]] || {
        echo "PRECOMPILE_BENCHMARK_CONSUMER_PROJECT and PRECOMPILE_BENCHMARK_CONSUMER_MANIFEST must be set together" >&2
        exit 2
    }
    for source_path in "$reuse_consumer_project" "$reuse_consumer_manifest"; do
        [[ $source_path != *$'\t'* && $source_path != *$'\n'* && $source_path != *$'\r'* ]] || {
            echo "consumer environment paths must not contain tabs or line endings" >&2
            exit 2
        }
        [[ -f $source_path && -r $source_path ]] || {
            echo "consumer environment input is not a readable regular file: $source_path" >&2
            exit 2
        }
    done
    if ! IFS= read -r -d '' reuse_consumer_project < <(realpath -ze -- "$reuse_consumer_project"); then
        echo "failed to canonicalize consumer Project input" >&2
        exit 2
    fi
    if ! IFS= read -r -d '' reuse_consumer_manifest < <(realpath -ze -- "$reuse_consumer_manifest"); then
        echo "failed to canonicalize consumer Manifest input" >&2
        exit 2
    fi
    for source_path in "$reuse_consumer_project" "$reuse_consumer_manifest"; do
        [[ $source_path != *$'\t'* && $source_path != *$'\n'* && $source_path != *$'\r'* ]] || {
            echo "canonical consumer environment paths must not contain tabs or line endings" >&2
            exit 2
        }
    done
    [[ ! $reuse_consumer_project -ef $reuse_consumer_manifest ]] || {
        echo "consumer Project and Manifest inputs must be distinct files" >&2
        exit 2
    }
fi
command -v "$julia" >/dev/null 2>&1 || { echo "Julia executable not found: $julia" >&2; exit 2; }
julia_version=$("$julia" --version)
if [[ $julia_version != "$expected_julia" && $allow_julia_mismatch != 1 ]]; then
    echo "cold-start comparisons require $expected_julia (found $julia_version)" >&2
    echo "set PRECOMPILE_BENCHMARK_ALLOW_JULIA_MISMATCH=1 only for a non-reportable smoke run" >&2
    exit 2
fi

baseline_candidates=()
baseline_labels=()
if [[ -n $baseline_map_list ]]; then
    IFS=',' read -r -a baseline_specifications <<< "$baseline_map_list"
    for specification in "${baseline_specifications[@]}"; do
        [[ $specification == *=* ]] || {
            echo "PRECOMPILE_BENCHMARK_BASELINES entries must have CANDIDATE=BASELINE form" >&2
            exit 2
        }
        baseline_candidate=${specification%%=*}
        baseline_label=${specification#*=}
        [[ -n $baseline_candidate && -n $baseline_label ]] || {
            echo "invalid baseline entry: $specification" >&2
            exit 2
        }
        baseline_candidates+=("$baseline_candidate")
        baseline_labels+=("$baseline_label")
    done
fi
[[ $(uname -s) == Linux ]] || {
    echo "the cold-start harness currently requires GNU/Linux" >&2
    exit 2
}
labels=()
checkouts=()
for specification in "$@"; do
    [[ $specification == *=* ]] || usage
    label=${specification%%=*}
    checkout=${specification#*=}
    [[ $label =~ $token_pattern ]] || {
        echo "variant labels must start with an ASCII letter or digit and contain only letters, digits, dots, underscores, or hyphens" >&2
        exit 2
    }
    [[ -n $checkout && $checkout != *$'\t'* && $checkout != *$'\n'* ]] || {
        echo "variant checkouts must be nonempty and must not contain tabs or newlines" >&2
        exit 2
    }
    [[ -f "$checkout/Project.toml" && -d "$checkout/src" && -f "$checkout/lib/QECCore/Project.toml" ]] || {
        echo "not a QuantumClifford checkout: $checkout" >&2
        exit 2
    }
    checkout=$(cd -- "$checkout" && pwd -P)
    [[ $checkout != *$'\t'* && $checkout != *$'\n'* ]] || {
        echo "canonical variant checkout paths must not contain tabs or newlines" >&2
        exit 2
    }
    for existing_label in "${labels[@]}"; do
        [[ $label != "$existing_label" ]] || { echo "duplicate variant label: $label" >&2; exit 2; }
    done
    labels+=("$label")
    checkouts+=("$checkout")
done

for checkout in "${checkouts[@]}"; do
    if [[ $output_dir == "$checkout" || $output_dir == "$checkout/"* ]]; then
        echo "output directory must be outside every measured checkout: $output_dir" >&2
        exit 2
    fi
done

checkout_state_sha256() {
    local checkout=$1
    local file_path file_sha256 file_state link_target relative_path
    (
        printf 'status\0'
        git -C "$checkout" status --porcelain=v1 -z --untracked-files=all
        printf 'tracked-diff\0'
        git -C "$checkout" diff --no-ext-diff --no-textconv --binary HEAD --
        printf '\0untracked\0'
        while IFS= read -r -d '' relative_path; do
            file_path="$checkout/$relative_path"
            printf 'path\0%s\0' "$relative_path"
            if [[ -L $file_path ]]; then
                link_target=$(readlink -- "$file_path")
                printf 'symlink\0%s\0' "$link_target"
            elif [[ -f $file_path ]]; then
                file_state=$(stat -c '%a:%s' -- "$file_path")
                file_sha256=$(sha256sum < "$file_path" | awk '{print $1}')
                printf 'file\0%s\0%s\0' "$file_state" "$file_sha256"
            elif [[ -d $file_path ]]; then
                file_state=$(stat -c '%a' -- "$file_path")
                printf 'directory\0%s\0' "$file_state"
            elif [[ -e $file_path ]]; then
                file_state=$(stat -c '%F:%a:%s' -- "$file_path")
                printf 'other\0%s\0' "$file_state"
            else
                printf 'missing\0'
            fi
        done < <(git -C "$checkout" ls-files -z --others --exclude-standard)
    ) | sha256sum | awk '{print $1}'
}

qeccore_state_sha256() {
    local checkout=$1
    local file_path file_sha256 file_state link_target relative_path
    (
        printf 'head-tree\0'
        git -C "$checkout" rev-parse HEAD:lib/QECCore
        printf 'status\0'
        git -C "$checkout" status --porcelain=v1 -z --untracked-files=all -- lib/QECCore
        printf 'tracked-diff\0'
        git -C "$checkout" diff --no-ext-diff --no-textconv --binary HEAD -- lib/QECCore
        printf '\0untracked\0'
        while IFS= read -r -d '' relative_path; do
            file_path="$checkout/$relative_path"
            printf 'path\0%s\0' "$relative_path"
            if [[ -L $file_path ]]; then
                link_target=$(readlink -- "$file_path")
                printf 'symlink\0%s\0' "$link_target"
            elif [[ -f $file_path ]]; then
                file_state=$(stat -c '%a:%s' -- "$file_path")
                file_sha256=$(sha256sum < "$file_path" | awk '{print $1}')
                printf 'file\0%s\0%s\0' "$file_state" "$file_sha256"
            elif [[ -d $file_path ]]; then
                file_state=$(stat -c '%a' -- "$file_path")
                printf 'directory\0%s\0' "$file_state"
            elif [[ -e $file_path ]]; then
                file_state=$(stat -c '%F:%a:%s' -- "$file_path")
                printf 'other\0%s\0' "$file_state"
            else
                printf 'missing\0'
            fi
        done < <(git -C "$checkout" ls-files -z --others --exclude-standard -- lib/QECCore)
    ) | sha256sum | awk '{print $1}'
}

copy_checkout_snapshot() {
    local checkout=$1
    local destination=$2
    local destination_path relative_path source_path
    while IFS= read -r -d '' relative_path; do
        source_path="$checkout/$relative_path"
        [[ -e $source_path || -L $source_path ]] || continue
        destination_path="$destination/$relative_path"
        mkdir -p -- "$(dirname -- "$destination_path")"
        cp -a --reflink=never -- "$source_path" "$destination_path"
    done < <(git -C "$checkout" ls-files -z --cached --others --exclude-standard)
}

commits=()
checkout_initial_dirty=()
checkout_state_sha256s=()
qeccore_head_trees=()
qeccore_state_sha256s=()
for index in "${!labels[@]}"; do
    checkout=${checkouts[$index]}
    commits+=("$(git -C "$checkout" rev-parse HEAD)")
    qeccore_head_trees+=("$(git -C "$checkout" rev-parse HEAD:lib/QECCore)")
    if [[ -n $(git -C "$checkout" status --porcelain=v1 --untracked-files=all) ]]; then
        checkout_initial_dirty+=(true)
    else
        checkout_initial_dirty+=(false)
    fi
    if [[ ${checkout_initial_dirty[$index]} == true && $allow_dirty != 1 ]]; then
        echo "variant checkout is dirty; commit it or set PRECOMPILE_BENCHMARK_ALLOW_DIRTY=1 for a non-reportable smoke run: ${labels[$index]}" >&2
        exit 1
    fi
    checkout_state_sha256s+=("$(checkout_state_sha256 "$checkout")")
    qeccore_state_sha256s+=("$(qeccore_state_sha256 "$checkout")")
done

verify_checkout() {
    local index=$1
    local checkout=${checkouts[$index]}
    local expected_commit=${commits[$index]}
    [[ $(git -C "$checkout" rev-parse HEAD) == "$expected_commit" ]] || {
        echo "variant HEAD changed during measurement: ${labels[$index]}" >&2
        exit 1
    }
    [[ $(checkout_state_sha256 "$checkout") == "${checkout_state_sha256s[$index]}" ]] || {
        echo "variant checkout contents changed during measurement: ${labels[$index]}" >&2
        exit 1
    }
    [[ $(git -C "$checkout" rev-parse HEAD:lib/QECCore) == "${qeccore_head_trees[$index]}" &&
       $(qeccore_state_sha256 "$checkout") == "${qeccore_state_sha256s[$index]}" ]] || {
        echo "variant bundled QECCore contents changed during measurement: ${labels[$index]}" >&2
        exit 1
    }
}

for checkout in "${checkouts[@]:1}"; do
    cmp -s -- \
        <(sed '1,/^\[/ { /^version = /d; }' "${checkouts[0]}/Project.toml") \
        <(sed '1,/^\[/ { /^version = /d; }' "$checkout/Project.toml") || {
        echo "all variants must use identical Project.toml dependency metadata" >&2
        exit 1
    }
    cmp -s -- \
        <(sed '1,/^\[/ { /^version = /d; }' "${checkouts[0]}/lib/QECCore/Project.toml") \
        <(sed '1,/^\[/ { /^version = /d; }' "$checkout/lib/QECCore/Project.toml") || {
        echo "all variants must use identical lib/QECCore/Project.toml dependency metadata" >&2
        exit 1
    }
done
for mapping_index in "${!baseline_candidates[@]}"; do
    baseline_candidate=${baseline_candidates[$mapping_index]}
    baseline_label=${baseline_labels[$mapping_index]}
    candidate_known=false
    baseline_known=false
    for label in "${labels[@]}"; do
        [[ $label == "$baseline_candidate" ]] && candidate_known=true
        [[ $label == "$baseline_label" ]] && baseline_known=true
    done
    $candidate_known || {
        echo "baseline map refers to an unknown candidate label: $baseline_candidate" >&2
        exit 2
    }
    $baseline_known || {
        echo "baseline map refers to an unknown baseline label: $baseline_label" >&2
        exit 2
    }
    [[ $baseline_candidate != "${labels[0]}" ]] || {
        echo "the first variant cannot be a mapped candidate: $baseline_candidate" >&2
        exit 2
    }
    [[ $baseline_candidate != "$baseline_label" ]] || {
        echo "a candidate cannot be its own baseline: $baseline_candidate" >&2
        exit 2
    }
    for previous_index in "${!baseline_candidates[@]}"; do
        [[ $previous_index -ge $mapping_index ]] && break
        [[ ${baseline_candidates[$previous_index]} != "$baseline_candidate" ]] || {
            echo "duplicate baseline map for candidate: $baseline_candidate" >&2
            exit 2
        }
    done
done

candidate_baseline_indices=()
for ((candidate_index = 1; candidate_index < ${#labels[@]}; candidate_index++)); do
    baseline_index=0
    for mapping_index in "${!baseline_candidates[@]}"; do
        [[ ${baseline_candidates[$mapping_index]} == "${labels[$candidate_index]}" ]] || continue
        for label_index in "${!labels[@]}"; do
            if [[ ${labels[$label_index]} == "${baseline_labels[$mapping_index]}" ]]; then
                baseline_index=$label_index
                break
            fi
        done
        break
    done
    candidate_baseline_indices[$candidate_index]=$baseline_index
done

set_comparison_order() {
    local build_number=$1
    local candidate_index
    comparison_indices=()
    if ((build_number % 2 == 1)); then
        for ((candidate_index = 1; candidate_index < ${#labels[@]}; candidate_index++)); do
            comparison_indices+=("$candidate_index")
        done
    else
        for ((candidate_index = ${#labels[@]} - 1; candidate_index >= 1; candidate_index--)); do
            comparison_indices+=("$candidate_index")
        done
    fi
}

set_pair_order() {
    local build_number=$1
    local candidate_index=$2
    local baseline_index=${candidate_baseline_indices[$candidate_index]}
    if (((build_number + candidate_index) % 2 == 1)); then
        pair_indices=("$baseline_index" "$candidate_index")
    else
        pair_indices=("$candidate_index" "$baseline_index")
    fi
}

temporary_root=$(mktemp -d "${TMPDIR:-/tmp}/quantumclifford-precompile.XXXXXX")
cleanup() {
    if [[ ${PRECOMPILE_BENCHMARK_KEEP_TMP:-0} == 1 ]]; then
        echo "kept temporary benchmark data at $temporary_root" >&2
    else
        rm -rf -- "$temporary_root"
    fi
}
trap cleanup EXIT

environment_dir="$temporary_root/environment"
seed_depot="$temporary_root/seed-depot"
checkout_link="$temporary_root/checkout"
seed_checkout="$temporary_root/seed-checkout"
harness_snapshot_dir="$temporary_root/harness"
scenario_script="$harness_snapshot_dir/scenarios.jl"
summarize_script="$harness_snapshot_dir/summarize.jl"
verify_harness_snapshots() {
    [[ $(sha256sum "$scenario_script" | awk '{print $1}') == "$harness_scenarios_sha256" ]] || {
        echo "recorded scenario harness snapshot changed during measurement" >&2
        exit 1
    }
    [[ $(sha256sum "$summarize_script" | awk '{print $1}') == "$harness_summarize_sha256" ]] || {
        echo "summary harness snapshot changed during measurement" >&2
        exit 1
    }
}
mkdir -p -- "$environment_dir" "$seed_depot" "$seed_checkout" "$harness_snapshot_dir"
cp -- "$script_dir/scenarios.jl" "$scenario_script"
cp -- "$script_dir/summarize.jl" "$summarize_script"
chmod 0444 "$scenario_script" "$summarize_script"
verify_harness_snapshots
# Seed dependency caches through separate inodes so setup cannot page-warm a measured checkout.
verify_checkout 0
if [[ ${checkout_initial_dirty[0]} == true ]]; then
    seed_checkout_mode=dirty_working_tree_copy
    # A non-reportable dirty run still uses its fixed working source rather than HEAD.
    copy_checkout_snapshot "${checkouts[0]}" "$seed_checkout"
else
    seed_checkout_mode=committed_git_archive
    git -C "${checkouts[0]}" archive --format=tar "${commits[0]}" | tar -xf - -C "$seed_checkout"
fi
verify_checkout 0
[[ -f "$seed_checkout/Project.toml" && -d "$seed_checkout/src" ]] || {
    echo "failed to create the detached seed checkout" >&2
    exit 1
}
ln -s -- "$seed_checkout" "$checkout_link"
[[ $checkout_link != *$'\t'* && $checkout_link != *$'\n'* && $checkout_link != *$'\r'* && $checkout_link != *'"'* && $checkout_link != *'\'* ]] || {
    echo "temporary checkout link cannot be represented safely in the consumer Manifest: $checkout_link" >&2
    exit 1
}

export JULIA_NUM_THREADS=1
export JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_LOAD_PATH='@:@stdlib'
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export LC_ALL=C
unset JULIA_CPU_TARGET JULIA_PROJECT JULIA_DEPOT_PATH
unset PRECOMPILE_BENCHMARK_TRACE
unset JULIA_PKG_OFFLINE
julia_flags=(--startup-file=no --history-file=no --threads=1)

echo "Preparing one consumer manifest and seed depot..." >&2
consumer_environment_mode=resolved
consumer_environment_source_project=
consumer_environment_source_manifest=
consumer_environment_source_project_sha256=
consumer_environment_source_manifest_sha256=
manifest_path="$environment_dir/Manifest.toml"
quantumclifford_manifest_placeholder=__QUANTUMCLIFFORD_CHECKOUT__
qeccore_manifest_placeholder=__QECCORE_CHECKOUT__
qeccore_checkout_path="$checkout_link/lib/QECCore"
if [[ -n $reuse_consumer_project ]]; then
    consumer_environment_mode=reused
    consumer_environment_source_project=$reuse_consumer_project
    consumer_environment_source_manifest=$reuse_consumer_manifest
    validate_consumer_project "$reuse_consumer_project" || {
        echo "reused consumer Project does not match the harness-generated Project" >&2
        exit 2
    }
    cp -- "$reuse_consumer_project" "$environment_dir/Project.toml"
    rewrite_checkout_manifest_paths \
        "$reuse_consumer_manifest" "$environment_dir/Manifest.toml" \
        "$quantumclifford_manifest_placeholder" "$checkout_link" \
        "$qeccore_manifest_placeholder" "$qeccore_checkout_path" || {
        echo "reused consumer Manifest must contain the QuantumClifford and QECCore checkout placeholders in their respective path entries" >&2
        exit 2
    }
    consumer_environment_source_project_sha256=$(sha256sum "$reuse_consumer_project" | awk '{print $1}')
    consumer_environment_source_manifest_sha256=$(sha256sum "$reuse_consumer_manifest" | awk '{print $1}')
    if ! JULIA_PKG_OFFLINE=true JULIA_DEPOT_PATH="$seed_depot" "$julia" "${julia_flags[@]}" -e '
        using Pkg, TOML
        manifest = TOML.parsefile(joinpath(ARGS[1], "Manifest.toml"))
        dependencies = get(manifest, "deps", nothing)
        dependencies isa AbstractDict || error("consumer Manifest does not contain a deps table")
        entries = get(dependencies, "QuantumClifford", nothing)
        entries isa AbstractVector && length(entries) == 1 ||
            error("consumer Manifest must contain exactly one QuantumClifford entry")
        entry = only(entries)
        entry isa AbstractDict && get(entry, "path", nothing) == ARGS[2] ||
            error("consumer Manifest QuantumClifford path does not match the temporary checkout link")
        entries = get(dependencies, "QECCore", nothing)
        entries isa AbstractVector && length(entries) == 1 ||
            error("consumer Manifest must contain exactly one QECCore entry")
        entry = only(entries)
        entry isa AbstractDict && get(entry, "path", nothing) == ARGS[3] ||
            error("consumer Manifest QECCore path does not match the bundled checkout path")
        Pkg.is_manifest_current(ARGS[1]) === true || error("consumer Manifest was not resolved from the consumer Project")
    ' "$environment_dir" "$checkout_link" "$qeccore_checkout_path"; then
        echo "reused consumer Project and Manifest failed semantic validation" >&2
        exit 2
    fi
    JULIA_DEPOT_PATH="$seed_depot" "$julia" "${julia_flags[@]}" -e '
        using Pkg
        Pkg.activate(ARGS[1])
        Pkg.instantiate()
    ' "$environment_dir"
    cmp -s -- "$reuse_consumer_project" "$environment_dir/Project.toml" || {
        echo "Pkg.instantiate changed the reused consumer Project" >&2
        exit 1
    }
    reused_normalized_manifest="$temporary_root/reused-normalized-Manifest.toml"
    rewrite_checkout_manifest_paths \
        "$environment_dir/Manifest.toml" "$reused_normalized_manifest" \
        "$checkout_link" "$quantumclifford_manifest_placeholder" \
        "$qeccore_checkout_path" "$qeccore_manifest_placeholder" || {
        echo "Pkg.instantiate did not preserve the reused QuantumClifford and QECCore Manifest paths" >&2
        exit 1
    }
    cmp -s -- "$reuse_consumer_manifest" "$reused_normalized_manifest" || {
        echo "Pkg.instantiate changed the reused consumer Manifest" >&2
        exit 1
    }
else
    JULIA_DEPOT_PATH="$seed_depot" "$julia" "${julia_flags[@]}" -e '
        using Pkg
        Pkg.activate(ARGS[1])
        Pkg.develop(path=ARGS[2])
        Pkg.instantiate()
    ' "$environment_dir" "$checkout_link"

    resolved_qeccore_manifest_path=$("$julia" "${julia_flags[@]}" -e '
        using TOML
        manifest = TOML.parsefile(ARGS[1])
        print(only(manifest["deps"]["QECCore"])["path"])
    ' "$manifest_path")
    if [[ $resolved_qeccore_manifest_path = /* ]]; then
        resolved_qeccore_source=$(realpath -e -- "$resolved_qeccore_manifest_path")
    else
        resolved_qeccore_source=$(realpath -e -- "$environment_dir/$resolved_qeccore_manifest_path")
    fi
    expected_qeccore_source=$(realpath -e -- "$seed_checkout/lib/QECCore")
    [[ $resolved_qeccore_source == "$expected_qeccore_source" ]] || {
        echo "consumer Manifest QECCore path does not resolve to the bundled seed source: $resolved_qeccore_manifest_path" >&2
        exit 1
    }
    materialized_manifest="$temporary_root/materialized-Manifest.toml"
    rewrite_checkout_manifest_paths \
        "$manifest_path" "$materialized_manifest" \
        "$checkout_link" "$checkout_link" \
        "$resolved_qeccore_manifest_path" "$qeccore_checkout_path" || {
        echo "failed to point the consumer Manifest at the bundled QECCore checkout path" >&2
        exit 1
    }
    mv -- "$materialized_manifest" "$manifest_path"
fi

[[ -f $manifest_path ]] || { echo "consumer manifest was not created" >&2; exit 1; }
resolved_manifest_sha256=$(sha256sum "$manifest_path" | awk '{print $1}')
if ! "$julia" "${julia_flags[@]}" -e '
    using TOML
    manifest = TOML.parsefile(ARGS[1])
    dependencies = manifest["deps"]
    for (package, expected_path) in (("QuantumClifford", ARGS[2]), ("QECCore", ARGS[3]))
        entries = get(dependencies, package, nothing)
        entries isa AbstractVector && length(entries) == 1 ||
            error("consumer Manifest must contain exactly one $package entry")
        get(only(entries), "path", nothing) == expected_path ||
            error("consumer Manifest $package path does not match $expected_path")
    end
' "$manifest_path" "$checkout_link" "$qeccore_checkout_path"; then
    echo "consumer Manifest did not preserve the stable QuantumClifford and QECCore checkout paths" >&2
    exit 1
fi

if ! scenario_output=$(JULIA_DEPOT_PATH="$seed_depot" "$julia" "${julia_flags[@]}" --project="$environment_dir" "$scenario_script"); then
    echo "failed to discover PRECOMPILE_BENCHMARKS" >&2
    exit 1
fi
[[ -n $scenario_output ]] || {
    echo "PRECOMPILE_BENCHMARKS must not be empty" >&2
    exit 2
}
mapfile -t scenarios <<< "$scenario_output"
for scenario_index in "${!scenarios[@]}"; do
    scenario=${scenarios[$scenario_index]}
    [[ $scenario =~ $token_pattern ]] || {
        echo "invalid PRECOMPILE_BENCHMARKS key: $scenario" >&2
        exit 2
    }
    for previous_index in "${!scenarios[@]}"; do
        [[ $previous_index -ge $scenario_index ]] && break
        [[ ${scenarios[$previous_index]} != "$scenario" ]] || {
            echo "duplicate PRECOMPILE_BENCHMARKS key: $scenario" >&2
            exit 2
        }
    done
done
scenario_list=$(IFS=,; printf '%s' "${scenarios[*]}")
find "$seed_depot/compiled" -type f -path '*/QuantumClifford/*' -delete 2>/dev/null || true
find "$seed_depot/compiled" -type d -path '*/QuantumClifford' -empty -delete 2>/dev/null || true

# Give every measured source tree one equivalent discarded cache build. Content-hash
# order is independent of baseline/candidate role; equal source states retain input order.
discarded_cache_warmup_indices=()
while IFS=$'\t' read -r _ warmup_index; do
    discarded_cache_warmup_indices+=("$warmup_index")
done < <(
    for index in "${!labels[@]}"; do
        warmup_key=$(printf '%s\0%s\0' "${commits[$index]}" "${checkout_state_sha256s[$index]}" | sha256sum | awk '{print $1}')
        printf '%s\t%s\n' "$warmup_key" "$index"
    done | sort -k1,1 -k2,2n
)
discarded_cache_warmup_order=
for warmup_position in "${!discarded_cache_warmup_indices[@]}"; do
    index=${discarded_cache_warmup_indices[$warmup_position]}
    verify_checkout "$index"
    label=${labels[$index]}
    checkout=${checkouts[$index]}
    ln -sfn -- "$checkout" "$checkout_link"
    warmup_depot="$temporary_root/discarded-cache-warmup-$warmup_position"
    mkdir -p -- "$warmup_depot"
    echo "Discarding cache warm-up for $label..." >&2
    JULIA_PKG_OFFLINE=true JULIA_DEPOT_PATH="$warmup_depot:$seed_depot" \
        "$julia" "${julia_flags[@]}" --project="$environment_dir" -e 'using QuantumClifford'
    warmup_cache_bytes=$(find "$warmup_depot/compiled" -type f -path '*/QuantumClifford/*' -printf '%s\n' | awk '{ total += $1 } END { print total + 0 }')
    [[ $warmup_cache_bytes -gt 0 ]] || {
        echo "discarded QuantumClifford cache warm-up wrote no cache bytes: $label" >&2
        exit 1
    }
    [[ -z $discarded_cache_warmup_order ]] || discarded_cache_warmup_order+=,
    discarded_cache_warmup_order+="$label"
done

validate_consumer_project "$environment_dir/Project.toml" || {
    echo "consumer Project does not match the expected QuantumClifford environment" >&2
    exit 1
}
cp -- "$environment_dir/Project.toml" "$consumer_project_path"
rewrite_checkout_manifest_paths \
    "$manifest_path" "$consumer_manifest_path" \
    "$checkout_link" "$quantumclifford_manifest_placeholder" \
    "$qeccore_checkout_path" "$qeccore_manifest_placeholder" || {
    echo "failed to normalize the QuantumClifford and QECCore paths in the copied Manifest" >&2
    exit 1
}
if [[ $consumer_environment_mode == reused ]]; then
    cmp -s -- "$reuse_consumer_project" "$consumer_project_path" || {
        echo "copied consumer Project differs from the reused input" >&2
        exit 1
    }
    cmp -s -- "$reuse_consumer_manifest" "$consumer_manifest_path" || {
        echo "copied consumer Manifest differs from the reused input" >&2
        exit 1
    }
fi
consumer_project_sha256=$(sha256sum "$consumer_project_path" | awk '{print $1}')
normalized_manifest_sha256=$(sha256sum "$consumer_manifest_path" | awk '{print $1}')

non_reportable_reasons=()
[[ $allow_dirty == 0 ]] || non_reportable_reasons+=(allow_dirty_override_enabled)
[[ $allow_julia_mismatch == 0 ]] || non_reportable_reasons+=(julia_mismatch_override_enabled)
[[ $builds -ge 5 ]] || non_reportable_reasons+=(fewer_than_five_builds)
[[ $samples -ge 4 ]] || non_reportable_reasons+=(fewer_than_four_samples_per_build)
if [[ ${#non_reportable_reasons[@]} -eq 0 ]]; then
    reportable=true
    non_reportable_reason_list=
else
    reportable=false
    non_reportable_reason_list=$(IFS=,; echo "${non_reportable_reasons[*]}")
fi

{
    printf '%s\n' "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '%s\n' "julia=$julia_version"
    printf '%s\n' "kernel=$(uname -srvmo)"
    if [[ -r /etc/os-release ]]; then
        os_name=$(awk '$0 ~ /^PRETTY_NAME=/ { sub(/^[^=]*=/, ""); sub(/^"/, ""); sub(/"$/, ""); print; exit }' /etc/os-release)
        printf '%s\n' "os=$os_name"
    fi
    if command -v lscpu >/dev/null 2>&1; then
        printf '%s\n' "cpu=$(lscpu | awk -F: '/Model name/ { sub(/^[[:space:]]+/, "", $2); print $2; exit }')"
    fi
    printf '%s\n' "julia_num_threads=$JULIA_NUM_THREADS"
    printf '%s\n' "julia_num_precompile_tasks=$JULIA_NUM_PRECOMPILE_TASKS"
    printf '%s\n' "openblas_num_threads=$OPENBLAS_NUM_THREADS"
    printf '%s\n' "omp_num_threads=$OMP_NUM_THREADS"
    printf '%s\n' "julia_load_path=$JULIA_LOAD_PATH"
    printf '%s\n' "pkg_auto_precompile=$JULIA_PKG_PRECOMPILE_AUTO"
    printf '%s\n' "pkg_offline_during_measurement=true"
    printf '%s\n' "builds=$builds"
    printf '%s\n' "recorded_samples_per_build=$samples"
    printf '%s\n' "discarded_cache_warmups_per_variant=1"
    printf '%s\n' "discarded_cache_warmup_order_policy=sha256_of_commit_and_state_then_argument_index"
    printf '%s\n' "discarded_cache_warmup_order=$discarded_cache_warmup_order"
    printf '%s\n' "discarded_warmups_per_build_and_scenario=1"
    printf '%s\n' "reportable=$reportable"
    printf '%s\n' "nonreportable_reasons=$non_reportable_reason_list"
    printf '%s\n' "allow_dirty=$allow_dirty"
    printf '%s\n' "allow_julia_mismatch=$allow_julia_mismatch"
    printf '%s\n' "total_metric=wall_import_start_to_first_task_end"
    printf '%s\n' "scenarios=$scenario_list"
    printf '%s\n' "candidate_baselines=$baseline_map_list"
    printf '%s\n' "schedule_policy=counterbalanced-v1"
    printf '%s\n' "schedule_candidate_index=one_based_candidate_argument_position"
    printf '%s\n' "schedule_comparison_policy=ascending_candidate_index_on_odd_builds_descending_candidate_index_on_even_builds"
    printf '%s\n' "schedule_pair_policy=baseline_then_candidate_when_build_plus_candidate_index_is_odd_candidate_then_baseline_when_even"
    for ((schedule_candidate_index = 1; schedule_candidate_index < ${#labels[@]}; schedule_candidate_index++)); do
        printf '%s\n' "schedule.candidate_index.${labels[$schedule_candidate_index]}=$schedule_candidate_index"
    done
    for ((schedule_build = 1; schedule_build <= builds; schedule_build++)); do
        set_comparison_order "$schedule_build"
        schedule_value=
        for schedule_candidate_index in "${comparison_indices[@]}"; do
            set_pair_order "$schedule_build" "$schedule_candidate_index"
            baseline_index=${candidate_baseline_indices[$schedule_candidate_index]}
            pair_value=
            for schedule_variant_index in "${pair_indices[@]}"; do
                if [[ $schedule_variant_index -eq $baseline_index ]]; then
                    schedule_role=baseline
                else
                    schedule_role=candidate
                fi
                [[ -z $pair_value ]] || pair_value+=,
                pair_value+="$schedule_role:${labels[$schedule_variant_index]}"
            done
            [[ -z $schedule_value ]] || schedule_value+=';'
            schedule_value+="comparison:${labels[$schedule_candidate_index]},$pair_value"
        done
        printf '%s\n' "schedule.build.$schedule_build=$schedule_value"
    done
    printf '%s\n' "consumer_environment_mode=$consumer_environment_mode"
    printf '%s\n' "seed_checkout_mode=$seed_checkout_mode"
    printf '%s\n' "consumer_project_sha256=$consumer_project_sha256"
    if [[ $consumer_environment_mode == reused ]]; then
        printf '%s\n' "consumer_environment_source_project=$consumer_environment_source_project"
        printf '%s\n' "consumer_environment_source_manifest=$consumer_environment_source_manifest"
        printf '%s\n' "consumer_environment_source_project_sha256=$consumer_environment_source_project_sha256"
        printf '%s\n' "consumer_environment_source_manifest_sha256=$consumer_environment_source_manifest_sha256"
    fi
    printf '%s\n' "manifest_sha256=$normalized_manifest_sha256"
    printf '%s\n' "resolved_manifest_sha256=$resolved_manifest_sha256"
    printf '%s\n' "harness_commit=$harness_commit"
    printf '%s\n' "harness_run_sha256=$harness_run_sha256"
    printf '%s\n' "harness_scenarios_sha256=$harness_scenarios_sha256"
    printf '%s\n' "harness_summarize_sha256=$harness_summarize_sha256"
    for index in "${!labels[@]}"; do
        printf '%s\n' "variant.${labels[$index]}.checkout=${checkouts[$index]}"
        printf '%s\n' "variant.${labels[$index]}.commit=${commits[$index]}"
        printf '%s\n' "variant.${labels[$index]}.initial_dirty=${checkout_initial_dirty[$index]}"
        printf '%s\n' "variant.${labels[$index]}.state_sha256=${checkout_state_sha256s[$index]}"
        printf '%s\n' "variant.${labels[$index]}.qeccore_head_tree=${qeccore_head_trees[$index]}"
        printf '%s\n' "variant.${labels[$index]}.qeccore_state_sha256=${qeccore_state_sha256s[$index]}"
    done
} > "$metadata_path"

printf '%s\n' $'comparison\tlabel\tcheckout\tcommit\tbuild\tsample\tscenario\tbuild_seconds\tcache_bytes\timport_seconds\tfirst_seconds\tfirst_compile_seconds\tfirst_recompile_seconds\ttotal_seconds\twarm_seconds\twarm_compile_seconds\twarm_recompile_seconds' > "$raw_path"

export JULIA_PKG_OFFLINE=true
for build in $(seq 1 "$builds"); do
    set_comparison_order "$build"
    for candidate_index in "${comparison_indices[@]}"; do
        comparison=${labels[$candidate_index]}
        baseline_index=${candidate_baseline_indices[$candidate_index]}
        set_pair_order "$build" "$candidate_index"
        for index in "${pair_indices[@]}"; do
            verify_checkout "$index"
            label=${labels[$index]}
            checkout=${checkouts[$index]}
            commit=${commits[$index]}
            ln -sfn -- "$checkout" "$checkout_link"

            if [[ $index -eq $baseline_index ]]; then
                role=baseline
            else
                role=candidate
            fi
            run_depot="$temporary_root/run-$candidate_index-$role-$build"
            mkdir -p -- "$run_depot"
            echo "Building cache for $comparison: $label ($build/$builds)..." >&2
            build_started=$(date +%s.%N)
            JULIA_DEPOT_PATH="$run_depot:$seed_depot" "$julia" "${julia_flags[@]}" --project="$environment_dir" -e 'using QuantumClifford'
            build_finished=$(date +%s.%N)
            build_seconds=$(awk -v started="$build_started" -v finished="$build_finished" 'BEGIN { printf "%.9f", finished - started }')
            cache_bytes=$(find "$run_depot/compiled" -type f -path '*/QuantumClifford/*' -printf '%s\n' | awk '{ total += $1 } END { print total + 0 }')
            [[ $cache_bytes -gt 0 ]] || { echo "QuantumClifford cache was not written to the run depot" >&2; exit 1; }

            for scenario in "${scenarios[@]}"; do
                [[ -n $scenario && $scenario != *[[:space:]]* ]] || {
                    echo "scenario names must be nonempty and contain no whitespace: $scenario" >&2
                    exit 2
                }
                echo "Discarding filesystem warm-up for $comparison: $label/$scenario build $build..." >&2
                JULIA_DEPOT_PATH="$run_depot:$seed_depot" "$julia" "${julia_flags[@]}" --project="$environment_dir" "$scenario_script" "$scenario" >/dev/null

                for sample in $(seq 1 "$samples"); do
                    echo "Sampling $comparison: $label/$scenario build $build ($sample/$samples)..." >&2
                    scenario_output=$(JULIA_DEPOT_PATH="$run_depot:$seed_depot" "$julia" "${julia_flags[@]}" --project="$environment_dir" "$scenario_script" "$scenario")
                    result=$(printf '%s\n' "$scenario_output" | awk -F '\t' '$1 == "RESULT" { print; count += 1 } END { if (count != 1) exit 1 }') || {
                        echo "scenario did not emit exactly one RESULT row: $label/$scenario" >&2
                        exit 1
                    }
                    IFS=$'\t' read -r marker measured_scenario import_seconds first_seconds first_compile first_recompile total_seconds warm_seconds warm_compile warm_recompile <<< "$result"
                    [[ $marker == RESULT && $measured_scenario == "$scenario" ]] || {
                        echo "malformed scenario result: $result" >&2
                        exit 1
                    }
                    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                        "$comparison" "$label" "$checkout" "$commit" "$build" "$sample" "$scenario" \
                        "$build_seconds" "$cache_bytes" "$import_seconds" "$first_seconds" \
                        "$first_compile" "$first_recompile" "$total_seconds" "$warm_seconds" \
                        "$warm_compile" "$warm_recompile" >> "$raw_path"
                done
            done
        done
    done
done

for index in "${!labels[@]}"; do
    verify_checkout "$index"
done
verify_harness_checkout
verify_harness_snapshots

JULIA_DEPOT_PATH="$seed_depot" "$julia" "${julia_flags[@]}" "$summarize_script" \
    "$raw_path" "$summary_path" "$build_summary_path" "$markdown_path"
echo "Wrote $raw_path, $summary_path, $build_summary_path, and $markdown_path" >&2
