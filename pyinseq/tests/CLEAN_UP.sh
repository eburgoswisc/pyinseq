#!/usr/bin/env bash
# Remove the past pyinseq test directroy with cwd as context
# shellcheck disable=SC1049
if [ -d "results/" ]; then
  rm -r results/
fi