"""
Download processed intermediate files from Zenodo.

These files are the output of the QIIME2 pipeline (05_qiime2_pipeline.sh)
and are required to run the downstream analysis scripts (06-14) without
re-running the full pipeline.

Usage:
    python download_processed_data.py

After download, run:
    md5sum -c processed_data.md5
to verify file integrity.

Zenodo record: https://doi.org/10.5281/zenodo.XXXXXXX  ← fill in after upload
"""

import urllib.request
import tarfile
import hashlib
import sys
from pathlib import Path

# ── Update this URL after uploading to Zenodo ────────────────────────────────
ZENODO_URL = "https://zenodo.org/records/XXXXXXX/files/version-2_integrated.tar.gz"
ARCHIVE_NAME = "version-2_integrated.tar.gz"
EXPECTED_MD5 = "43b578ea7540787fd39f4b9d5735dcbc"   # md5sum of the tar.gz itself
# ─────────────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).parent.resolve()


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    print(f"Downloading {dest.name} ...", flush=True)
    def progress(count, block, total):
        pct = min(count * block / total * 100, 100)
        print(f"\r  {pct:5.1f}%", end="", flush=True)
    urllib.request.urlretrieve(url, dest, reporthook=progress)
    print()


def verify_md5(path, expected):
    actual = md5(path)
    if actual != expected:
        print(f"MD5 mismatch for {path.name}")
        print(f"  expected: {expected}")
        print(f"  actual:   {actual}")
        sys.exit(1)
    print(f"MD5 OK: {path.name}")


def extract(archive, dest):
    print(f"Extracting to {dest} ...")
    with tarfile.open(archive) as tar:
        tar.extractall(path=dest)
    print("Done.")


if __name__ == "__main__":
    if "XXXXXXX" in ZENODO_URL:
        print("ERROR: Zenodo URL has not been set.")
        print("Update ZENODO_URL and EXPECTED_MD5 in this script after uploading.")
        sys.exit(1)

    archive_path = REPO_ROOT / ARCHIVE_NAME

    if not archive_path.exists():
        download(ZENODO_URL, archive_path)
    else:
        print(f"Archive already present: {archive_path}")

    verify_md5(archive_path, EXPECTED_MD5)
    extract(archive_path, REPO_ROOT)

    print("\nVerifying individual files ...")
    import subprocess
    result = subprocess.run(
        ["md5sum", "-c", "processed_data.md5"],
        cwd=REPO_ROOT,
        capture_output=True, text=True
    )
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        sys.exit(1)

    print("All files verified. Run 'python config.py' to confirm setup.")
