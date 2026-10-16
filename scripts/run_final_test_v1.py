"""Separately authorized, one-attempt test-only evaluation entry point."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oncorna.final_test import AuthorizationError  # noqa: E402
from oncorna.final_test_consumer import AUTH_PATH, run_real  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization-record", type=Path, default=ROOT / AUTH_PATH)
    args = parser.parse_args()
    try:
        run_real(ROOT, args.authorization_record.resolve())
    except AuthorizationError as exc:
        print(f"Final test remains locked: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
