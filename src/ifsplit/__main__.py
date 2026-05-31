"""Enable `python -m ifsplit ...` (no install required)."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
