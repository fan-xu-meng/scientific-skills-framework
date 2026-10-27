#!/usr/bin/env python3
"""Lightweight repository checks that do not require third-party packages."""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILLS = ROOT / "skills"
NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
LOCAL_PATH_RE = re.compile(r"/(?:Users|home)/[^/\s]+/|[A-Za-z]:\\Users\\")
AFFILIATION_RE = re.compile(r"puda|pudap", re.IGNORECASE)
SECRET_RES = (
    re.compile(r"ghp_[A-Za-z0-9]{20,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]{20,}"),
    re.compile(r"sk-[A-Za-z0-9_-]{20,}"),
    re.compile(r"BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY"),
)
LINK_RE = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)|!\[[^\]]*\]\(([^)]+)\)")
IGNORED_LINK_PREFIXES = ("http://", "https://", "mailto:", "#")
FORBIDDEN_PARTS = {"__pycache__", "output", "outputs", ".pytest_cache", ".ruff_cache"}


def frontmatter(text: str, path: Path) -> dict[str, str]:
    if not text.startswith("---\n"):
        raise ValueError(f"{path}: missing YAML frontmatter")
    try:
        header = text.split("\n---\n", 1)[0].splitlines()[1:]
    except IndexError as exc:
        raise ValueError(f"{path}: unclosed YAML frontmatter") from exc
    values: dict[str, str] = {}
    for line in header:
        match = re.match(r"^(name|description):\s*(.*)$", line)
        if match:
            values[match.group(1)] = match.group(2).strip().strip('"\'')
    return values


def check_links(path: Path, text: str, errors: list[str]) -> None:
    for match in LINK_RE.finditer(text):
        target = (match.group(1) or match.group(2)).strip().split("#", 1)[0]
        if not target or target.startswith(IGNORED_LINK_PREFIXES) or "<" in target:
            continue
        if not (path.parent / target).resolve().exists():
            errors.append(f"{path.relative_to(ROOT)}: missing local link target {target!r}")


def main() -> int:
    errors: list[str] = []
    skill_dirs = sorted(path for path in SKILLS.iterdir() if path.is_dir())
    if not skill_dirs:
        errors.append("no skill directories found")

    for skill_dir in skill_dirs:
        skill_file = skill_dir / "SKILL.md"
        if not skill_file.is_file():
            errors.append(f"{skill_dir.relative_to(ROOT)}: missing SKILL.md")
            continue
        text = skill_file.read_text(encoding="utf-8")
        try:
            props = frontmatter(text, skill_file)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        name = props.get("name", "")
        if name != skill_dir.name:
            errors.append(f"{skill_dir.name}: frontmatter name is {name!r}")
        if not NAME_RE.fullmatch(name):
            errors.append(f"{skill_dir.name}: invalid skill name")
        if "description" not in props:
            errors.append(f"{skill_dir.name}: missing description")

    for path in ROOT.rglob("*"):
        if any(part in FORBIDDEN_PARTS for part in path.parts):
            errors.append(f"forbidden generated path: {path.relative_to(ROOT)}")
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.suffix.lower() in {".pyc", ".pyo"} or path.name == ".DS_Store":
            errors.append(f"forbidden generated file: {path.relative_to(ROOT)}")
        if path.suffix.lower() not in {".md", ".py", ".yaml", ".yml", ".json", ".txt"}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if path.suffix.lower() == ".py":
            try:
                ast.parse(text, filename=str(path))
            except SyntaxError as exc:
                errors.append(f"Python syntax error in {path.relative_to(ROOT)}: {exc}")
        if LOCAL_PATH_RE.search(text):
            errors.append(f"local filesystem path found: {path.relative_to(ROOT)}")
        if AFFILIATION_RE.search(text) and path.resolve() != Path(__file__).resolve():
            errors.append(f"legacy affiliation found: {path.relative_to(ROOT)}")
        for pattern in SECRET_RES:
            if pattern.search(text):
                errors.append(f"possible secret found: {path.relative_to(ROOT)}")
        if path.suffix.lower() == ".md":
            check_links(path, text, errors)

    if errors:
        print("Repository validation failed:")
        for error in sorted(set(errors)):
            print(f"- {error}")
        return 1
    print(f"Validated {len(skill_dirs)} skills and repository hygiene checks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

