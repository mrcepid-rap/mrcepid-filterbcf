#!/bin/bash

echo "🔍 Searching for pytest temp folders in /private/var/folders..."

# Find and remove pytest temp directories
deleted=$(find /private/var/folders -type d -name 'pytest-of-*' -prune)

if [[ -z "$deleted" ]]; then
    echo "✅ No pytest temp folders found."
else
    echo "🗑 Deleting the following pytest temp folders:"
    echo "$deleted"
    find /private/var/folders -type d -name 'pytest-of-*' -prune -exec rm -rf {} +
    echo "✅ Done. Pytest temp folders deleted."
fi
