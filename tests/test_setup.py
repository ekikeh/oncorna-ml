from pathlib import Path

import oncorna


def test_package_imports() -> None:
    assert oncorna.__version__ == "0.1.0"


def test_project_configuration_exists() -> None:
    project_root = Path(__file__).resolve().parents[1]
    assert (project_root / "configs" / "default.yaml").is_file()
