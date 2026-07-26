#!/usr/bin/env bash
# Verify the parallel worker is writing each file to the path the manifest
# specifies. A misaligned xargs argument split would silently write correct
# bytes to wrong paths, which no size check would catch.
cd "$(dirname "$0")/../.." || exit 1

echo "=== files on disk vs manifest-expected destinations ==="
bad=0
while IFS=$'\t' read -r ds fname url size dest; do
  [ "$ds" = "dataset_id" ] && continue
  for p in "$dest" "${dest}.part"; do
    if [ -f "$p" ]; then
      printf '  ok   %-11s %s\n' "$ds" "$(basename "$p")"
    fi
  done
done < data/metadata/download_list.tsv | head -20

echo
echo "=== any file NOT matching a manifest destination? (would indicate misalignment) ==="
find data/raw -type f ! -name '.gitkeep' | while read -r f; do
  base="${f%.part}"
  if ! cut -f5 data/metadata/download_list.tsv | grep -qxF "$base"; then
    echo "  UNEXPECTED: $f"
    bad=1
  fi
done
echo "  (nothing listed above = all paths match the manifest)"

echo
echo "=== per-dataset directory contents ==="
for d in data/raw/*/; do
  n=$(find "$d" -type f ! -name '.gitkeep' | wc -l)
  [ "$n" -gt 0 ] && printf '  %-34s %s file(s)  %s\n' "$(basename "$d")" "$n" "$(du -sh "$d" | cut -f1)"
done
