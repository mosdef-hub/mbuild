__author__ = 'sallai'


import httplib2
from HTMLParser import HTMLParser

class MyHTMLParser(HTMLParser):

    def __init__(self):
        HTMLParser.__init__(self)
        self.molids = {}

    def handle_starttag(self, tag, attrs):

        if tag == 'a':
            for attr,v in attrs:
                if attr == 'href':
                    if v.startswith('./molecule.py?molid='):
                        molid = v[v.find('=')+1:]
                        # print "molid=<{}>".format(molid)
                        self.molids[int(molid)]="http://compbio.biosci.uq.edu.au/atb"+v[1:]

class AtbClient(object):

    def __init__(self):
        self.h = httplib2.Http(".cache")

    def search(self, query):
        resp, content = self.h.request("http://compbio.biosci.uq.edu.au/atb/index.py?molsPerPage=1000&search={}".format(query), "GET")
        parser = MyHTMLParser()
        parser.feed(content)

        return parser.molids

if __name__ == "__main__":
    atb = AtbClient()
    content = atb.search("C6H12O6")
    print content

