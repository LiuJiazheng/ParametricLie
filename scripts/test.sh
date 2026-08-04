#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
exec julia --project=. -e "using Pkg; Pkg.test()"
