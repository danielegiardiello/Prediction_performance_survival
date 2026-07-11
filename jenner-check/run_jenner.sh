#!/usr/bin/env bash
# run_jenner.sh - mac/linux runner for Jenner compatibility checks.
#
# Quick start:
#   cd jenner-check/
#   ./run_jenner.sh                  # lists bundles in the current dir
#   ./run_jenner.sh t001_something   # run that one
#   ./run_jenner.sh --all            # run every bundle in the current dir
#
# Usage:   ./run_jenner.sh [bundle-dir | script.sas | --all | --list] [response.json]
#          ./run_jenner.sh --compare RESPONSE.json EXPECTED.json
#
#   (no arg)     If the current directory has tNNN_* bundles, list them
#                with a copy-paste command. Otherwise show this help.
#
#   --all        Run every tNNN_* bundle in the current directory in
#                sequence, print a pass/fail summary.
#
#   --list, -l   List the bundles visible in the current directory and
#                exit without running anything.
#
#   --compare    Purely local, no network: check RESPONSE.json against the
#                pinned EXPECTED.json and print one line per check. Only
#                keys present in EXPECTED.json are checked (keys starting
#                with "_" are provenance and ignored): status / exit_code
#                must match exactly, every log_contains entry must appear
#                in the response log, no log_does_not_contain entry may
#                appear, output_contains entries must appear in the output,
#                and pinned diagnostics lists must match exactly. Unknown
#                keys print a warning and do not fail.
#                Exit 0 all-match, 1 any mismatch, 2 usage/unreadable file/
#                invalid JSON, 5 python3 unavailable.
#
#   bundle-dir   A directory containing script.sas and (optionally)
#                autoexec.sas. The two are concatenated (autoexec first,
#                then a blank line, then script) and submitted together.
#                This is the normal case.
#
#   script.sas   A single .sas file. Submitted as-is — no autoexec.
#
# The API response is written to <response.json> (or response.json in
# the current directory if omitted) and the most useful fields are also
# printed to stdout for a quick sanity check.
#
# What "pass" means: a bundle run passes only if the POST returns HTTP 200
# AND — when the bundle carries an expected.json and python3 is available —
# the fresh response matches every pinned field (the --compare semantics
# above). Without python3 the pinned comparison is skipped with a notice
# and the run falls back to the old HTTP/status-only behavior.
#
# Requires: bash 4+, curl. Both ship with every mainstream Linux distro
# and macOS 12+. python3 (present on virtually every modern system) is
# optional but needed for the pinned-result comparison and the pretty
# summary. Windows: use run_jenner.bat (single-file mode) or WSL.
#
# IMPORTANT: execute this script, don't source it. Running with `. ./...`
# or `source ./...` will short-circuit error handling and can close your
# terminal if an error path fires.

# --- refuse to be sourced ------------------------------------------------
# `return` only works inside a sourced script. If we ARE sourced, print a
# message and return 1 so we don't kill the parent shell with exit. If
# we're running directly, (return 0) fails and we fall through.
(return 0 2>/dev/null) && {
  printf 'run_jenner.sh: execute this script, do not source it.\n  ./run_jenner.sh <bundle-dir-or-script.sas>\n' >&2
  return 1
}

set -eu

# --- helpers -------------------------------------------------------------
# Emit the list of tNNN_* bundles in the current working directory. A
# "bundle" is a directory matching t[0-9]*_* whose name contains a
# script.sas file. Writes one path per line (no prefix); empty output
# if nothing found.
list_bundles_here() {
  local d
  for d in ./t[0-9]*_*/ ; do
    [[ -d "$d" && -f "$d/script.sas" ]] || continue
    printf '%s\n' "${d%/}"     # strip trailing slash, keep leading ./
  done
}

# Render a helpful listing + copy-paste suggestion, then exit non-zero
# (we haven't done anything). Used when the user runs with no args.
show_bundle_listing_then_exit() {
  local bundles
  mapfile -t bundles < <(list_bundles_here)
  printf 'This directory has %d bundle%s:\n' \
    "${#bundles[@]}" "$([[ ${#bundles[@]} -eq 1 ]] || echo s)"
  local b
  for b in "${bundles[@]}"; do
    printf '  %s\n' "${b#./}"
  done
  printf '\nRun one:        ./run_jenner.sh %s\n' "${bundles[0]#./}"
  printf 'Run them all:   ./run_jenner.sh --all\n'
  printf 'Just list:      ./run_jenner.sh --list\n'
  exit 2
}

# Show the usage block when we have nothing better to offer.
show_usage_then_exit() {
  local status=${1:-2}
  {
    printf 'Usage: %s [bundle-dir | script.sas | --all | --list] [response.json]\n' "$(basename "$0")"
    printf '       %s --compare RESPONSE.json EXPECTED.json\n\n' "$(basename "$0")"
    printf 'Examples:\n'
    printf '  %s t001_my_bundle         # run one bundle\n' "$(basename "$0")"
    printf '  %s --all                  # run every tNNN_* bundle in this dir\n' "$(basename "$0")"
    printf '  %s path/to/script.sas     # run a single file, no autoexec\n' "$(basename "$0")"
    printf '  %s --compare r.json e.json # check a response against pinned fields\n' "$(basename "$0")"
  } >&2
  exit "$status"
}

# Compare a response.json against a pinned expected.json. Purely local —
# no network. Only keys present in expected.json are checked; keys whose
# name starts with "_" are capture provenance and are ignored. Prints one
# line per check (ok/FAIL/warn). Return codes: 0 all checks match,
# 1 any mismatch, 2 unreadable file or invalid JSON, 5 python3 missing.
compare_response_to_expected() {
  local resp_file=$1 expected_file=$2
  if ! command -v python3 >/dev/null 2>&1; then
    printf 'error: the pinned comparison needs python3, which was not found\n' >&2
    return 5
  fi
  python3 - "$resp_file" "$expected_file" <<'PY'
import json, sys

def load(path, label):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError) as e:
        print(f"compare: cannot read {label} file {path}: {e}", file=sys.stderr)
        sys.exit(2)

resp = load(sys.argv[1], "response")
exp = load(sys.argv[2], "expected")
if not isinstance(resp, dict) or not isinstance(exp, dict):
    print("compare: both files must contain a JSON object", file=sys.stderr)
    sys.exit(2)

failures = 0

def ok(msg):
    print(f"  ok    {msg}")

def bad(msg):
    global failures
    failures += 1
    print(f"  FAIL  {msg}")

def as_list(v):
    return v if isinstance(v, list) else [v]

log = resp.get("log") or ""
output = resp.get("output") or ""

for key, want in exp.items():
    if key.startswith("_"):
        continue                      # provenance (_captured_at, ...), not a check
    if key in ("status", "exit_code"):
        got = resp.get(key)
        if got == want:
            ok(f"{key} == {want!r}")
        else:
            bad(f"{key}: expected {want!r}, got {got!r}")
    elif key == "log_contains":
        for s in as_list(want):
            if s in log:
                ok(f"log contains {s!r}")
            else:
                bad(f"log_contains: log is missing {s!r}")
    elif key == "log_does_not_contain":
        for s in as_list(want):
            if s not in log:
                ok(f"log does not contain {s!r}")
            else:
                bad(f"log_does_not_contain: log unexpectedly contains {s!r}")
    elif key == "output_contains":
        for s in as_list(want):
            if s in output:
                ok(f"output contains {s!r}")
            else:
                bad(f"output_contains: output is missing {s!r}")
    elif key == "diagnostics":
        diag = resp.get("diagnostics") or {}
        for dk, dv in (want or {}).items():
            if dk.startswith("_"):
                continue
            got = diag.get(dk)
            if got == dv:
                ok(f"diagnostics.{dk} == {dv!r}")
            else:
                bad(f"diagnostics.{dk}: expected {dv!r}, got {got!r}")
    else:
        print(f"  warn  unknown expected.json key {key!r} — not checked")

if failures:
    print(f"pinned comparison: {failures} check(s) FAILED")
    sys.exit(1)
print("pinned comparison: all pinned checks passed")
PY
}

# --- arg parsing ---------------------------------------------------------
if [[ $# -lt 1 ]]; then
  # No args: if the cwd contains bundles, list them; otherwise show help.
  mapfile -t _found < <(list_bundles_here)
  if [[ ${#_found[@]} -gt 0 ]]; then
    show_bundle_listing_then_exit
  fi
  show_usage_then_exit 2
fi

HOST=${JENNER_HOST:-api.jenneranalytics.com}

case "$1" in
  -h|--help)
    show_usage_then_exit 0
    ;;
  -l|--list)
    mapfile -t _found < <(list_bundles_here)
    if [[ ${#_found[@]} -eq 0 ]]; then
      printf 'No tNNN_* bundles found in %s\n' "$(pwd)"
      exit 0
    fi
    printf 'Bundles in %s:\n' "$(pwd)"
    for b in "${_found[@]}"; do
      printf '  %s\n' "${b#./}"
    done
    exit 0
    ;;
  --compare)
    # Purely local: diff a response.json against a pinned expected.json.
    if [[ $# -ne 3 ]]; then
      printf 'usage: %s --compare RESPONSE.json EXPECTED.json\n' "$(basename "$0")" >&2
      exit 2
    fi
    _rc=0
    compare_response_to_expected "$2" "$3" || _rc=$?
    exit "$_rc"
    ;;
  --all)
    mapfile -t _found < <(list_bundles_here)
    if [[ ${#_found[@]} -eq 0 ]]; then
      printf 'No tNNN_* bundles found in %s\n' "$(pwd)" >&2
      exit 3
    fi
    _pass=0; _fail=0
    for b in "${_found[@]}"; do
      printf '\n── %s ──\n' "${b#./}"
      if "$0" "$b" "${b#./}_response.json"; then
        _pass=$((_pass+1))
      else
        _fail=$((_fail+1))
      fi
    done
    printf '\n── summary: %d pass, %d fail ──\n' "$_pass" "$_fail"
    [[ $_fail -eq 0 ]] && exit 0 || exit 1
    ;;
esac

TARGET=$1
OUT=${2:-response.json}

# --- assemble the submission body ---------------------------------------
# If TARGET is a directory, treat it as a bundle. If it's a file, submit
# it directly.
CLEANUP=()
cleanup() {
  for f in "${CLEANUP[@]}"; do rm -f "$f"; done
}
trap cleanup EXIT

if [[ -d "$TARGET" ]]; then
  if [[ ! -f "$TARGET/script.sas" ]]; then
    printf 'error: %s is a directory but has no script.sas\n' "$TARGET" >&2
    exit 3
  fi
  SUBMIT=$(mktemp -t jc_submit.XXXXXX.sas)
  CLEANUP+=("$SUBMIT")
  if [[ -f "$TARGET/autoexec.sas" ]]; then
    cat "$TARGET/autoexec.sas" > "$SUBMIT"
    printf '\n' >> "$SUBMIT"
  fi
  cat "$TARGET/script.sas" >> "$SUBMIT"
  printf 'Submitting bundle: %s\n' "$TARGET"
  if [[ -f "$TARGET/autoexec.sas" ]]; then
    printf '  autoexec.sas (%d bytes) + script.sas (%d bytes)\n' \
      "$(wc -c < "$TARGET/autoexec.sas")" "$(wc -c < "$TARGET/script.sas")"
  else
    printf '  script.sas (%d bytes), no autoexec\n' "$(wc -c < "$TARGET/script.sas")"
  fi
elif [[ -f "$TARGET" ]]; then
  SUBMIT=$TARGET
  printf 'Submitting file: %s (%d bytes)\n' "$TARGET" "$(wc -c < "$TARGET")"
else
  printf 'error: %s is neither a file nor a directory\n' "$TARGET" >&2
  exit 3
fi

# --- POST ---------------------------------------------------------------
printf 'POST https://%s/v1/run ... ' "$HOST"
HTTP_CODE=$(curl -sS -o "$OUT" -w '%{http_code}' -X POST \
  "https://${HOST}/v1/run" \
  -F "script=@${SUBMIT};type=application/x-sas" \
  -F "deterministic=1" \
  -F "timeout=60")
printf 'HTTP %s\n' "$HTTP_CODE"

if [[ "$HTTP_CODE" != "200" ]]; then
  printf 'API returned non-200 — raw response in %s\n' "$OUT" >&2
  exit 4
fi

# --- summarise ----------------------------------------------------------
# Best-effort: use python if present, otherwise grep key fields.
printf 'Response written to %s\n' "$OUT"
if command -v python3 >/dev/null 2>&1; then
  python3 - "$OUT" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
print(f"  status     : {r.get('status')}")
print(f"  exit_code  : {r.get('exit_code')}")
print(f"  duration_ms: {r.get('duration_ms')}")
print(f"  run_id     : {r.get('run_id')}")
print(f"  jenner_ver : {r.get('jenner_version')}")
log = r.get('log', '')
if log:
    print('  log (first 10 lines):')
    for line in log.splitlines()[:10]:
        print(f'    {line}')
PY
else
  printf '  (install python3 for a pretty summary; raw JSON in %s)\n' "$OUT"
fi

# --- pinned-result comparison ---------------------------------------------
# A bundle that ships an expected.json only passes if the fresh response
# matches every pinned field — HTTP 200 alone is not a pass. Without
# python3 we can't parse JSON portably, so fall back to the old behavior
# and say so.
if [[ -d "$TARGET" && -f "$TARGET/expected.json" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    printf 'Comparing response against pinned %s/expected.json\n' "${TARGET%/}"
    _rc=0
    compare_response_to_expected "$OUT" "$TARGET/expected.json" || _rc=$?
    if [[ $_rc -ne 0 ]]; then
      printf 'Bundle FAILED the pinned comparison (see lines above)\n' >&2
      exit "$_rc"
    fi
  else
    printf 'pinned comparison skipped (python3 not found) — HTTP/status only\n'
  fi
fi
