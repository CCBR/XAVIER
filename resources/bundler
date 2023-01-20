#!/bin/bash

set -euo pipefail 

usage="Usage: $0 <path_to_resource_bundle_directory> <output_archive.tar.gz>"

function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; err "$usage"; exit 1; }
function abspath() { readlink -f "$1"; }
function bundle() { tar -hczvf "$1" "$2"; }

 
function check() {
	die=false
	if [ -z "${1:-}" ]; then die=true; err "Error: Failed to provide directory to archive."; fi
	if [ -z "${2:-}" ]; then die=true; err "Error: Failed to output file name for archive."; fi
	if $die; then fatal "Fatal: Please try again after providing the required arguments!"; fi

}


function main() {
	# Checks for required positional
	# command line arguments 
	check "${1:-}" "${2:-}" 

	# Converts any relative paths to
	# absolute paths, creates uninit
	# output directories as needed,
	# runs tar command in parent dir
	# of the provided resource bundle
	# path.
	archive_dir=$(abspath "$1")
	archive=$(basename "${archive_dir%/}")
	parent_dir=$(dirname "$archive_dir")
	output_dir=$(abspath "$2")
	output_dir=$(dirname "$output_dir")
	mkdir -p "$output_dir"
	cd "$parent_dir"

	# Create archive as a tarball
	bundle "$2" "$archive"
}


main "$@"
