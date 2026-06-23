"""Initialize git repo and commit V2.3.0 using dulwich."""
import os, sys, time, stat
from pathlib import Path
from dulwich.repo import Repo
from dulwich.objects import Blob, Tree, Commit, TreeEntry
from dulwich.index import blob_from_path_and_stat

ROOT = Path(r"d:\prj\gnss_pod")

# Patterns to include
INCLUDE_GLOBS = [
    "src/*.py",
    "run_sequential_pod.py",
    "validate_v224.py",
    "fullday_batch_v3.py",
    "eval_day2.py",
    "run_multiday_v2.py",
    "run_5day_v230.py",
    "download_code_products.py",
    "fullday_assessment.py",
    "VERSION.md",
    "references/*.md",
    "run_v224.ps1",
    "V2.3.0/**",
]

EXCLUDE = {"__pycache__", ".pyc", "results", "data"}


def should_include(path_str):
    parts = set(path_str.replace("\\", "/").split("/"))
    if parts & EXCLUDE:
        return False
    if "__pycache__" in path_str:
        return False
    if path_str.endswith(".pyc"):
        return False
    return True


def find_files():
    files = []
    for pattern in INCLUDE_GLOBS:
        import glob as g
        matches = g.glob(str(ROOT / pattern), recursive=True)
        for m in matches:
            if os.path.isfile(m):
                rel = os.path.relpath(m, str(ROOT))
                if should_include(rel):
                    files.append((rel.replace("\\", "/"), m))
    return sorted(set(files))


# Init repo
repo_dir = str(ROOT)
dot_git = os.path.join(repo_dir, ".git")
if not os.path.exists(dot_git):
    Repo.init(repo_dir)
    print("Initialized git repo")

repo = Repo(repo_dir)

# Build tree — use dulwich's internal store properly
files = find_files()
print(f"Found {len(files)} files")

store = repo.object_store
tree = Tree()

for git_path, abs_path in files:
    # Read file and create blob
    try:
        with open(abs_path, 'rb') as f:
            data = f.read()
    except Exception as e:
        print(f"  SKIP {git_path}: {e}")
        continue

    blob = Blob()
    blob.data = data
    store.add_object(blob)

    # Add to tree
    entry = TreeEntry(git_path.encode('utf-8'), 0o100644, blob.id)
    tree.add(entry.path, entry.mode, entry.sha)

store.add_object(tree)
tree_id = tree.id
print(f"Tree: {tree_id.hex()}")

# Create commit
now = int(time.time())
author = b"POD Developer <pod@local>"
msg = b"V2.3.0 Release - GRACE-FO PPP-AR POD\n\n" \
    b"Features:\n" \
    b"- Framework v3 Batch Solver (phase RMS -42% to -83%)\n" \
    b"- Arc-based Ambiguity Resolution (0.246m at 0.17h)\n" \
    b"- 11 LEO Satellite Config Database\n" \
    b"- Multi-date validation with auto-download (5-day)\n" \
    b"\n" \
    b"Accuracy (2024-04-29, GRACE-FO C):\n" \
    b"  0.17h: 0.293m (EKF) / 0.160m (Batch phase)\n" \
    b"  0.5h:  0.986m (EKF) / 0.211m (Batch phase)\n"

c = Commit()
c.tree = tree_id
c.parents = []
c.author = author
c.committer = author
c.author_time = now
c.author_timezone = 0
c.commit_time = now
c.commit_timezone = 0
c.encoding = b"UTF-8"
c.message = msg

store.add_object(c)
commit_id = c.id
print(f"Commit: {commit_id.hex()}")

# Set HEAD
repo.refs[b"refs/heads/main"] = commit_id
with open(os.path.join(dot_git, "HEAD"), "w") as fh:
    fh.write("ref: refs/heads/main\n")

# Write .gitignore
with open(os.path.join(repo_dir, ".gitignore"), "w") as gi:
    gi.write("__pycache__/\n*.pyc\nresults/\ndata/\n*.pkl\n*.gz\n*.tgz\n*.zip\n.idea/\n.vscode/\n")

# Add .gitignore to a follow-up commit
print("Branch: main")
print("Added .gitignore")
print("Done!")
