#!/usr/bin/env bash

python3 -m venv .venv
if . .venv/bin/activate; then
    echo "Building Requirements"
    python -m pip uninstall -r requirements.dev
    python -m pip uninstall -e biopytools
    
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt
    pip3 freeze > requirements.frozen
else
    echo "Could not activate venv. Make sure it is correctly configured using 'python3 -m venv .venv'"
    exit
fi
