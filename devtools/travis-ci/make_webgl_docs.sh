# Creates HTML that includes webgl visualization
#
# Input: file.rst
# Output: file.html
# Tools required:
#	- pandoc
#	- notedown
# 	- ipython

# example: make_webgl_docs.sh file.rst file.html

# convert RST to temp MD
pandoc -i $1 -o tmp.md

# convert temp MD to IPYNB
notedown tmp.md -o $1.ipynb

# delete temp MD
rm tmp.md

# convert ipynb to html with custom template
ipython nbconvert --to html --template ./tpl/mbuild_ipynb_template.tpl $1.ipynb --output $2

# delete temp ipynb
rm $1.ipynb