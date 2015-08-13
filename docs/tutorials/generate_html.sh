#!/bin/sh

mkdirs ~/.ipython/profile_default/static/custom/
cp custom.css ~/.ipython/profile_default/static/custom/custom.css
ipython nbconvert *.ipynb

