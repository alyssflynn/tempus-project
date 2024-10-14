#!/usr/bin/env bash

python3 -m venv .venv
if . .venv/bin/activate; then
    echo "Building Requirements\n-------------------------------"
    python -m pip install --upgrade pip

    if [ ! -s requirements.frozen ]
    then
        python -m pip install -r requirements.txt
        pip3 freeze > requirements.frozen
    else
        python -m pip install -U -r requirements.frozen
    fi
    python -m pip install -r requirements.dev
    python -m pip install -e .
else
    echo "Could not activate venv. Make sure it is correctly configured using 'python3 -m venv .venv'"
    exit $1
fi

