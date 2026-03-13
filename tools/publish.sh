#!/bin/sh
set -eu

# publish.sh — publish a committed R package snapshot to a pkg/<package> branch
#
# Usage: tools/publish.sh <package-root>
#
#   <package-root>  Directory containing the package DESCRIPTION file.
#                   May be the root of its Git repo or a nested subdirectory.

die() { printf '%s\n' "$@" >&2; exit 1; }
warn() { printf 'warning: %s\n' "$@" >&2; }

# --- cleanup trap -----------------------------------------------------------

tmp_worktree=""
tmp_branch=""
cleanup() {
    if [ -n "$tmp_worktree" ]; then
        git -C "$universe_repo_root" worktree remove --force "$tmp_worktree" 2>/dev/null || true
    fi
    if [ -n "$tmp_branch" ]; then
        git -C "$universe_repo_root" branch -D "$tmp_branch" 2>/dev/null || true
    fi
}
trap cleanup EXIT

# --- arguments ---------------------------------------------------------------

[ $# -eq 1 ] || die "usage: tools/publish.sh <package-root>"
package_root="$(cd "$1" && pwd)" || die "cannot resolve path: $1"

# --- step 1-2: verify package root ------------------------------------------

[ -f "$package_root/DESCRIPTION" ] || die "$package_root does not contain DESCRIPTION"
git -C "$package_root" rev-parse --show-toplevel >/dev/null 2>&1 \
    || die "$package_root is not inside a Git repository"

# --- step 3-4: derive private repo root and subdir --------------------------

private_repo_root="$(git -C "$package_root" rev-parse --show-toplevel)"
if [ "$package_root" = "$private_repo_root" ]; then
    subdir=""
else
    # relative path from repo root to package root
    subdir="${package_root#"$private_repo_root"/}"
fi

# --- step 5: derive universe repo root --------------------------------------

script_dir="$(cd "$(dirname "$0")" && pwd)"
universe_repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"

# --- step 6: read Package and Version from committed HEAD --------------------

desc_path="${subdir:+"$subdir"/}DESCRIPTION"
git -C "$private_repo_root" show "HEAD:$desc_path" >/dev/null 2>&1 \
    || die "cannot read HEAD:$desc_path from $private_repo_root"

package="$(git -C "$private_repo_root" show "HEAD:$desc_path" | sed -n 's/^Package: *//p')"
version="$(git -C "$private_repo_root" show "HEAD:$desc_path" | sed -n 's/^Version: *//p')"
[ -n "$package" ] || die "could not parse Package from $desc_path"
[ -n "$version" ] || die "could not parse Version from $desc_path"

printf 'package: %s\n' "$package"
printf 'version: %s\n' "$version"

# --- step 7: refuse if package subtree is dirty -----------------------------

dirty="$(git -C "$private_repo_root" status --porcelain -- "${subdir:-.}")"
if [ -n "$dirty" ]; then
    printf '%s\n' "$dirty" >&2
    die "package subtree has uncommitted changes — commit first"
fi

# --- step 8: check packages.json lists this package -------------------------

packages_json="$universe_repo_root/packages.json"
if [ -f "$packages_json" ]; then
    if ! grep -q "\"$package\"" "$packages_json"; then
        warn "$package is not yet listed in packages.json (add it after this publish succeeds)"
    fi
else
    warn "packages.json does not exist yet"
fi

# --- step 9: fetch current remote state -------------------------------------

git -C "$universe_repo_root" fetch origin "refs/heads/pkg/$package:refs/remotes/origin/pkg/$package" 2>/dev/null || true

# --- step 10: create temporary orphan worktree --------------------------------

tmp_worktree="$(mktemp -d "${TMPDIR:-/tmp}/publish-$package-XXXXXX")"
# remove the dir so git worktree add can create it
rmdir "$tmp_worktree"
git -C "$universe_repo_root" worktree add --detach "$tmp_worktree" >/dev/null 2>&1
# Use a temporary branch name to avoid conflicts with existing refs.
# The push refspec (HEAD:refs/heads/pkg/$package) targets the real branch.
tmp_branch="_publish_tmp_$$"
(
    cd "$tmp_worktree"
    git switch --orphan "$tmp_branch" >/dev/null 2>&1
)

# --- step 11: export committed snapshot into worktree ------------------------

if [ -n "$subdir" ]; then
    git -C "$private_repo_root" archive "HEAD:$subdir" | tar -xC "$tmp_worktree"
else
    git -C "$private_repo_root" archive HEAD | tar -xC "$tmp_worktree"
fi

# --- step 12: stage all -----------------------------------------------------

git -C "$tmp_worktree" add -A

# --- step 13: no-op check ---------------------------------------------------

if git -C "$universe_repo_root" rev-parse --verify "refs/remotes/origin/pkg/$package" >/dev/null 2>&1; then
    remote_tree="$(git -C "$universe_repo_root" rev-parse "refs/remotes/origin/pkg/$package^{tree}")"
    staged_tree="$(git -C "$tmp_worktree" write-tree)"
    if [ "$remote_tree" = "$staged_tree" ]; then
        printf 'no changes — %s %s is already published\n' "$package" "$version"
        exit 0
    fi
fi

# --- step 14: commit ---------------------------------------------------------

git -C "$tmp_worktree" commit -m "publish $package $version" >/dev/null

# --- step 15: push -----------------------------------------------------------

git -C "$tmp_worktree" push --force-with-lease origin "HEAD:refs/heads/pkg/$package"

# Update local tracking ref so "git checkout pkg/<package>" shows current content
git -C "$universe_repo_root" fetch origin "refs/heads/pkg/$package:refs/remotes/origin/pkg/$package" 2>/dev/null || true
# Remove any stale local branch that would shadow the remote
git -C "$universe_repo_root" branch -D "pkg/$package" 2>/dev/null || true

printf 'published %s %s\n' "$package" "$version"
