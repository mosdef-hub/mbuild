from __future__ import print_function

import httplib2
try:
    from HTMLParser import HTMLParser
except ImportError:
    from html.parser import HTMLParser
import urllib
import warnings


class SearchResultHTMLParser(HTMLParser):

    def __init__(self):
        HTMLParser.__init__(self)
        self.molids = {}

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            for attr, v in attrs:
                if attr == 'href':
                    if v.startswith('./molecule.py?molid='):
                        molid = v[v.find('=')+1:]
                        self.molids[int(molid)] = "http://compbio.biosci.uq.edu.au/atb"+v[1:]


class AtbClient(object):

    def __init__(self):
        self.h = httplib2.Http(".cache")

    def search(self, query):
        url = "http://compbio.biosci.uq.edu.au/atb/index.py?molsPerPage=1000&search={}".format(
            query)
        resp, content = self.h.request(url, "GET")
        if resp['status'] != '200':
            warnings.warn('HTTP response status is {} for URL "{}"'.format(resp['status'], url))
            return None

        parser = SearchResultHTMLParser()
        parser.feed(content)
        return parser.molids

    def generate_topology(self, molid, ff_version="53A6",):
        query_pairs = {"molid": str(molid), "ffVersion": ff_version,
                       "outputType": "top", "atbVersion": "v2Top",
                       "format": "GROMACS"}

        query_string = urllib.urlencode(query_pairs)

        url = "http://compbio.biosci.uq.edu.au/atb/molecule.py?{}".format(query_string)

        # print url

        resp, content = self.h.request(url, "GET")

        if resp['status'] != '200':
            warnings.warn('HTTP response status is {} for URL "{}"'.format(resp['status'], url))
            return None

        return content

    def retrieve_itp(self, molid, ff_version="53A6", all_atom=True):

        self.generate_topology(molid, ff_version)

        query_pairs = {"molid": str(molid), "ffVersion": ff_version,
                       "outputType": "top", "atbVersion": "v2Top"}
        if all_atom:
            query_pairs["file"] = "rtp_allatom"
        else:
            query_pairs["file"] = "rtp_uniatom"

        query_string = urllib.urlencode(query_pairs)

        url = "http://compbio.biosci.uq.edu.au/atb/download.py?{}".format(query_string)

        # print url

        resp, content = self.h.request(url, "GET")

        if resp['status'] != '200':
            warnings.warn('HTTP response status is {} for URL "{}"'.format(resp['status'], url))
            return None

        if not resp['content-type'].startswith('text/plain'):
            warnings.warn('Expecting text/plain response, got "{}" for URL "{}"'.format(
                resp['content-type'], url))
            return None

        return content


if __name__ == "__main__":
    atb = AtbClient()
    results = atb.search("C6H12O6")
    for molecule_id, uri in results.iteritems():
        print(atb.retrieve_itp(molecule_id))
