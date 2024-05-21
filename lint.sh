#!/usr/bin/env bash
flake8 --exit-zero chemkin_io
pylint --rcfile=.pylintrc chemkin_io
