#!/bin/sh
set -eu

# unpublish.sh — remove a package branch from the public repo
#
# Usage: tools/unpublish.sh <package>
#
# Prerequisites: remove the package from packages.json first, then run this.

die() { printf '%s\n' "$@" >&2; exit 1; }

[ $# -eq 1 ] || die "usage: tools/unpublish.sh <package>"
package="$1"

# --- derive universe repo root -----------------------------------------------

script_dir="$(cd "$(dirname "$0")" && pwd)"
universe_repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"

# --- step 1: verify packages.json no longer lists this package ---------------

packages_json="$universe_repo_root/packages.json"
[ -f "$packages_json" ] || die "packages.json not found"
if grep -q "\"$package\"" "$packages_json"; then
    die "$package is still listed in packages.json — remove it first"
fi

# --- step 2: delete remote branch -------------------------------------------

if git -C "$universe_repo_root" ls-remote --exit-code origin "refs/heads/pkg/$package" >/dev/null 2>&1; then
    git -C "$universe_repo_root" push origin --delete "pkg/$package"
    printf 'deleted remote branch pkg/%s\n' "$package"
else
    printf 'remote branch pkg/%s does not exist (already removed)\n' "$package"
fi

# --- step 3: delete local tracking ref --------------------------------------

if git -C "$universe_repo_root" rev-parse --verify "refs/remotes/origin/pkg/$package" >/dev/null 2>&1; then
    git -C "$universe_repo_root" branch -r -d "origin/pkg/$package" 2>/dev/null || true
fi

if git -C "$universe_repo_root" rev-parse --verify "refs/heads/pkg/$package" >/dev/null 2>&1; then
    git -C "$universe_repo_root" branch -D "pkg/$package" 2>/dev/null || true
fi

printf 'done\n'
