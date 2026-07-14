#!/usr/bin/env python3
from __future__ import annotations

import argparse

from viral_safe_target.cache import stage_is_current, write_stage_stamp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["check", "stamp"])
    parser.add_argument("--stamp", required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--output", action="append", default=[])
    parser.add_argument("--parameter", action="append", default=[])
    args = parser.parse_args()
    parameters = dict(item.split("=", 1) for item in args.parameter)
    if args.action == "check":
        raise SystemExit(
            0 if stage_is_current(args.stamp, args.output, args.input, parameters) else 1
        )
    write_stage_stamp(args.stamp, args.input, parameters)


if __name__ == "__main__":
    main()
