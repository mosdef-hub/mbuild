import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_sort(text):
    return [atoi(a) for a in re.split(r'(\d+)', text)]
