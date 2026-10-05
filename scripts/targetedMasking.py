#!/usr/bin/env python3
"""Auto-detect the ColabFold A3M layout and apply verified targeted masking."""
if __package__:
    from .masking_cli import main
else:
    from masking_cli import main

if __name__ == "__main__":
    main("auto")
