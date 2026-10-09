"""Final-test entry point. Real execution remains disabled in Phase 2D-B."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oncorna.final_test import AuthorizationError, run_real  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/final_test_v1.yaml")
    parser.add_argument("--authorization-record", type=Path)
    args = parser.parse_args()
    try:
        run_real(ROOT, args.config.resolve(), args.authorization_record)
    except AuthorizationError as exc:
        print(f"Final test remains locked: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
