#!/usr/bin/env python3
"""Verify a separate cluster build using Python, without container diffutils."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest() if hasattr(hashlib, 'file_digest') else _digest(stream)


def _digest(stream):
    result = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1024*1024), b''):
        result.update(chunk)
    return result.hexdigest()


def verify(build, repo, kaskade, manifest=None):
    # Match the shell's *.hh copy: AppleDouble/dotfiles are not build inputs.
    pairs = [(build/p.name, p) for p in repo.glob('*.hh') if not p.name.startswith('.')]
    pairs += [(build/name, repo/name) for name in ('EMI_model.cpp', 'Makefile')]
    pairs += [(build/'source-mg'/name, kaskade/'mg'/name)
              for name in ('bddc.hh', 'bddc.hpp', 'bddcCoarseReduction.hh')]
    hashes = {}
    for saved, current in pairs:
        value = digest(saved)
        if value != digest(current):
            raise ValueError('Build source differs from current source: '+str(current))
        hashes[str(saved.relative_to(build))] = value
    for name in ('emiModel', 'coarse-reduction-test'):
        path = build/name
        with path.open('rb') as stream:
            magic = stream.read(4)
        if magic != b'\x7fELF' or not os.access(path, os.X_OK):
            raise ValueError('Missing executable Linux build: '+str(path))
        hashes[name] = digest(path)
    manifests = [build/'verified-build.json']
    if manifest is not None:
        manifests.append(manifest)
    for recorded in manifests:
        if recorded.exists():
            previous = json.loads(recorded.read_text())
            if previous != hashes:
                raise ValueError('Build files changed after verification: '+str(recorded))
    return hashes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--kaskade', type=Path, required=True)
    p.add_argument('--write-manifest', action='store_true')
    p.add_argument('--manifest', type=Path, help='Verification record outside a reused read-only build')
    a = p.parse_args()
    hashes = verify(a.build, a.repo, a.kaskade, a.manifest)
    if a.write_manifest:
        (a.manifest or a.build/'verified-build.json').write_text(json.dumps(hashes, indent=2)+'\n')
    print('Separate build sources and executable hashes verified: '+str(a.build))


if __name__ == '__main__':
    main()
