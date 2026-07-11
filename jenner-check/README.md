# jenner-check: compatibility test bundles

Each `tNNN_<slug>/` subdirectory is a self-contained test bundle built from SAS
code that already lives in this repository. A bundle pins the output of a
captured passing run so the same script can be re-run later and compared
against that snapshot.

This directory is fully self-contained: nothing outside `jenner-check/` is
referenced or modified, nothing runs on merge or checkout, and deleting the
directory removes every trace of it.

## What the runner sends — read before running

`run_jenner.sh` uploads only a bundle's SAS **source text** — `autoexec.sas`
plus `script.sas`, and nothing else — over HTTPS to `api.jenneranalytics.com`,
where it runs and returns the log and listing. It does not read or upload any
data files, so anything sitting next to a script stays on your machine; the
only thing transmitted is the code you run, as when pasting a snippet into any
hosted tool. Nothing is sent unless you run a command yourself. To point the
runner at a different endpoint, set `JENNER_HOST`.

## What's in here

```
jenner-check/
├── README.md            # this file
├── run_jenner.sh        # runner (bash + curl; python3 used if present)
└── tNNN_<slug>/
    ├── script.sas       # the SAS under test (derived from this repo)
    ├── autoexec.sas     # run options prepended to script.sas at run time
    ├── expected.json    # fields pinned from the captured passing run
    ├── meta.json        # provenance: source file, blob sha, commit
    └── expected/        # human-readable snapshot of that run
        ├── log.txt
        └── output.txt
```

`meta.json` records the source path, blob sha, and commit of the script each
bundle was built from, so you can verify the copy matches your code.

## How to run

From inside this directory:

```bash
./run_jenner.sh --list             # show the bundles in this directory
./run_jenner.sh --all              # run every bundle, verify pinned fields
./run_jenner.sh tNNN_<slug>        # run just one (use a name from --list)
```

For each bundle the runner submits the SAS, writes the JSON response to
`<bundle>_response.json`, and compares the response against the bundle's
`expected.json` — status, exit code, and the pinned log lines. A bundle
passes only when every pinned field matches, so an `N pass, 0 fail` summary
from `--all` means the captured results still hold.

You can also compare a saved response offline, with no network call:

```bash
./run_jenner.sh --compare tNNN_<slug>_response.json tNNN_<slug>/expected.json
```

Requirements: `bash` 4+ and `curl` (both ship with mainstream Linux and
macOS); `python3` for the pinned-field comparison. On Windows, run it under
WSL.

## Removing this directory

`git rm -r jenner-check/` (or just delete the folder). Nothing else in the
repository references it.

---

`run_jenner.sh` and this README are provided under the MIT license. Each
bundle's `script.sas` is a copy of code that already lives in this repository
and remains under this repository's own license — nothing here relicenses it.
