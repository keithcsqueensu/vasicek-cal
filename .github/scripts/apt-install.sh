#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#
# apt-get install that rides out the archive outages seen on GitHub-hosted runners (D-418): every Linux
# job failed in apt with "Could not resolve 'archive.ubuntu.com'" / "'azure.archive.ubuntu.com'" on
# 2026-11-11, before anything was built.
#   1. Packages already installed are skipped (the runner images carry gcc, clang and ninja), and if
#      nothing is missing apt is not run at all.
#   2. Otherwise apt retries each download (Acquire::Retries) and the update + install pair is retried up
#      to 5 times, waiting 30, 60, 90 and 120 s between attempts (about 5 minutes in all).
# A package that cannot be installed after that still fails the job: this changes only how long a
# transient outage is waited out, never what is built or tested.
#   apt-install.sh PACKAGE...
set -euo pipefail
SUDO=""
[ "$(id -u)" -ne 0 ] && SUDO="sudo"
missing=()
for p in "$@"; do
    if [ "$p" = "ninja-build" ] && command -v ninja > /dev/null 2>&1; then continue; fi
    dpkg -s "$p" > /dev/null 2>&1 || missing+=("$p")
done
if [ ${#missing[@]} -eq 0 ]; then
    echo "apt-install: already installed: $*"
    exit 0
fi
opts=(-o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30)
for attempt in 1 2 3 4 5; do
    if $SUDO apt-get "${opts[@]}" update &&
        $SUDO env DEBIAN_FRONTEND=noninteractive apt-get "${opts[@]}" install -y --no-install-recommends "${missing[@]}"; then
        exit 0
    fi
    [ "$attempt" -eq 5 ] && break
    echo "apt-install: attempt $attempt failed; retrying in $((attempt * 30)) s"
    sleep $((attempt * 30))
done
echo "apt-install: could not install ${missing[*]} after 5 attempts" >&2
exit 1
