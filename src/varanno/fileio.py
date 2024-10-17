
import os
import json


def makefile(filename):
    """Create file, including parent dirs if nonexistant."""
    os.makedirs(os.path.dirname(filename), exist_ok=True)


def write_metadata_json(metadata: dict, outfile: str):
    """Generates a json file."""
    makefile(outfile)

    with open(outfile, "w") as fle:
        json.dump(metadata, fle)


def write_logs(messages: list, outfile: str):
    """Generates a log file."""
    makefile(outfile)

    with open(outfile, "wt") as fle:
        fle.writelines(messages)

