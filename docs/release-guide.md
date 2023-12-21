# Release Guide

Make sure you're keeping the changelog up-to-date during development.
Ideally, every PR that includes a user-facing change (e.g. a new feature, bug fix, or any API change) should add a concise summary to the changelog with a link to the PR.
Only approve or merge PRs that either update the changelog or have no user-facing changes.

## How to release a new version on GitHub

1. Determine the new version number according to [semantic versioning guidelines](https://semver.org/).
1. Update `CHANGELOG.md`:
    - Edit the heading for the development version to match the new version.
    - If needed, clean up the changelog -- fix any typos, optionally create subheadings for 'New features' and 'Bug fixes' if there are lots of changes, etc.
1. Update the version in [`src/__init__.py`](https://github.com/CCBR/XAVIER/blob/main/src/__init__.py).
1. On GitHub, go to "Releases" and click "Draft a new release". <https://github.com/CCBR/XAVIER/releases/new>
    - Choose a tag: same as the version number.
    - Choose the target: most likely this should be the main branch, or a specific commit hash.
    - Set the title as the new version number, e.g. **v3.0.2**
    - Copy and paste the release notes from the CHANGELOG into the description box.
    - Check the box "Set as the latest release".
    - Click "Publish release".
1. Post release chores:
    - Add a new "development version" heading to the top of `CHANGELOG.md`.
    - Bump the version number in `src/__init__.py` to include `-dev`, e.g. `v3.0.2-dev` if you just released `v3.0.2`.  

## How to install a release on biowulf

After releasing a new version on GitHub:

```sh
# go to the shared pipeline directory on biowulf
cd /data/CCBR_Pipeliner/Pipelines/XAVIER

# clone the new version tag (e.g. v3.0.2) to a hidden directory
git clone --depth 1 --branch v3.0.2 https://github.com/CCBR/XAVIER .v3.0.2

# change permissions for the new directory so anyone will be able to use the pipeline
chown -R :CCBR_Pipeliner .v3.0.2
chmod -R a+rX /data/CCBR_Pipeliner/Pipelines/XAVIER/.v3.0.2

# if needed, remove the old symlink for the minor version number
rm -i v3.0

# recreate the symlink to point to the new latest version
ln -s .v3.0.2 v3.0

# you can verify that the symlink points to the new version with readlink
readlink -f v3.0
```

Versions of the `ccbrpipeliner` module only specify the major and minor version of each pipeline.
If the new pipeline release only increments the patch number, `ccbrpipeliner` will use it automatically after you update the symlink as above.
If you need to release a new major or minor version of a pipeline on biowulf, contact [Kelly](mailto:kelly.sovacool@nih.gov) or [Vishal](mailto:vishal.koparde@nih.gov).

Verify that `ccbrpipeliner` uses the latest version with:
```sh
module load ccbrpipeliner && xavier --version
```
