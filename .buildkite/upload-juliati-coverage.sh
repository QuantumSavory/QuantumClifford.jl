#!/usr/bin/env bash
set -euo pipefail

report_args=()
for report in "$@"; do
  if [[ -f "$report" ]]; then
    report_args+=(--file "$report")
  fi
done

if ((${#report_args[@]} == 0)); then
  exit 0
fi

codecov_cli="$(mktemp)"
trap 'rm -f "$codecov_cli"' EXIT
curl --fail --silent --show-error --location \
  https://cli.codecov.io/v11.3.1/linux/codecov \
  --output "$codecov_cli"
echo "ca1d64196d2d34771084afe76ea657d581bf628e31d993ff8e52ea09cc88a56d  $codecov_cli" | \
  sha256sum --check --status
chmod +x "$codecov_cli"

"$codecov_cli" --auto-load-params-from BuildKite upload-coverage \
  --token "$CODECOV_TOKEN" \
  --slug QuantumSavory/QuantumClifford.jl \
  --git-service github \
  --plugin noop \
  --disable-search \
  "${report_args[@]}"
