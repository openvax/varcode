"""Package modules must not import private underscore symbols across module boundaries (#692)."""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

import varcode
from varcode.effects import (
    structural,
    structural_candidates,
    structural_fusions,
    structural_models,
)

PACKAGE_ROOT = Path(varcode.__file__).resolve().parent


def _iter_package_sources() -> list[tuple[str, Path, ast.Module]]:
    results = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        rel = path.relative_to(PACKAGE_ROOT).as_posix()
        results.append((rel, path, ast.parse(path.read_text(encoding="utf-8"), filename=rel)))
    return results


def _is_private_symbol(name: str) -> bool:
    return name.startswith("_") and not (name.startswith("__") and name.endswith("__"))


def _has_private_module_part(module_path: str) -> bool:
    return any(
        part.startswith("_") and part not in ("__future__", "__init__", "__main__")
        for part in module_path.split(".")
        if part
    )


def test_no_cross_module_private_symbol_imports():
    violations: list[str] = []
    for rel, _, tree in _iter_package_sources():
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                mod = "." * node.level + (node.module or "")
                for alias in node.names:
                    if _is_private_symbol(alias.name):
                        violations.append(f"{rel}:{node.lineno}: from {mod} import {alias.name}")
    assert violations == [], (
        "Cross-module imports of underscore-prefixed symbols are forbidden:\n"
        + "\n".join(violations)
    )


def test_no_private_module_imports():
    violations: list[str] = []
    for rel, _, tree in _iter_package_sources():
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                if _has_private_module_part(node.module):
                    violations.append(f"{rel}:{node.lineno}: from {node.module} import ...")
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if _has_private_module_part(alias.name):
                        violations.append(f"{rel}:{node.lineno}: import {alias.name}")
    assert violations == [], (
        "Imports from underscore-prefixed modules are forbidden:\n"
        + "\n".join(violations)
    )


def test_no_private_module_filenames_in_package():
    private_files = [
        path.relative_to(PACKAGE_ROOT).as_posix()
        for path in sorted(PACKAGE_ROOT.rglob("*.py"))
        if _is_private_symbol(path.stem)
    ]
    assert private_files == []


def test_sequence_change_does_not_import_structural_modules():
    seq_change = PACKAGE_ROOT / "effects" / "sequence_change.py"
    tree = ast.parse(seq_change.read_text(encoding="utf-8"))
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if (node.level == 1 and mod.startswith("structural")) or "effects.structural" in mod:
                imported.append(f"line {node.lineno}: from {'.' * node.level}{mod} import ...")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if "effects.structural" in alias.name:
                    imported.append(f"line {node.lineno}: import {alias.name}")
    assert imported == []


def test_structural_modules_are_split_by_responsibility():
    for mod in (structural, structural_models, structural_fusions, structural_candidates):
        source_lines = inspect.getsource(mod).splitlines()
        assert len(source_lines) < 550, f"{mod.__name__} grew back to {len(source_lines)} lines"
        assert inspect.getdoc(mod), f"{mod.__name__} is missing a module docstring"
    for name in structural.__all__:
        obj = getattr(structural, name)
        if callable(obj):
            assert inspect.getdoc(obj), f"structural.{name} is missing a docstring"


def test_effect_class_modules_are_split_by_domain():
    from varcode.effects import (
        base,
        coding_effects,
        effect_classes,
        splice_effects,
        structural_effects,
    )

    for mod in (base, splice_effects, coding_effects, structural_effects):
        source_lines = inspect.getsource(mod).splitlines()
        assert len(source_lines) < 750, f"{mod.__name__} grew back to {len(source_lines)} lines"
        assert inspect.getdoc(mod), f"{mod.__name__} is missing a module docstring"
    barrel_tree = ast.parse(inspect.getsource(effect_classes))
    barrel_classes = [n.name for n in barrel_tree.body if isinstance(n, ast.ClassDef)]
    assert barrel_classes == [], f"effect_classes.py should remain a barrel, found {barrel_classes}"


def test_transforms_germline_and_splice_modules_are_split_by_domain():
    from varcode import (
        germline,
        germline_phasing,
        splice_exon_builders,
        splice_intron_builders,
        splice_outcomes,
        transforms,
    )
    from varcode.transforms import breakends, cis_snvs, indels

    for mod, limit in (
        (breakends, 600),
        (indels, 300),
        (cis_snvs, 700),
        (germline, 850),
        (germline_phasing, 550),
        (splice_outcomes, 700),
        (splice_exon_builders, 500),
        (splice_intron_builders, 700),
    ):
        source_lines = inspect.getsource(mod).splitlines()
        assert len(source_lines) < limit, f"{mod.__name__} grew back to {len(source_lines)} lines"
        assert inspect.getdoc(mod), f"{mod.__name__} is missing a module docstring"
    transforms_tree = ast.parse(inspect.getsource(transforms))
    transforms_defs = [
        n.name for n in transforms_tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))
    ]
    assert transforms_defs == [], f"transforms/__init__.py should remain a barrel, found {transforms_defs}"
