#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Copyright 2019 GullumLuvl

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


"""Parse phylogenetic trees from a wikipedia page, and save it in a
bioinformatics format (e.g newick)."""


from __future__ import print_function

import os.path as op
#import string  # Check for whitespaces in text elements
import re
import requests
import bs4
import ete3
import argparse as ap

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s')

WIKIPEDIA_URL = 'https://en.wikipedia.org'
SEARCH_URL = WIKIPEDIA_URL + '/w/index.php'


# async def
def get_wiki_tree(term=None, url=SEARCH_URL):
    # await
    page = requests.get(url, params=({'search': term} if term else None))

    if not page.ok:
        logger.error("Status of request = %d (url=%r search=%r)", page.status_code,
                     url, term)

    soup = bs4.BeautifulSoup(page.text, 'lxml')
    page.close()

    #TODO: check if later matches in the document
    tablesoups = []
    tablesoup = soup.find('table', class_='clade')

    while tablesoup is not None:
        tablesoups.append(tablesoup)
        try:
            last = list(tablesoup.descendants)[-1]
        except IndexError:
            break
        tablesoup = last.findNext('table', class_='clade')

    if not tablesoups:
        if 'Search results' in soup.findChild('head').findChild('title').get_text():
            didyoumean = soup.find('div', class_='searchdidyoumean')
            if didyoumean:
                showed, original = [a.get_text() for a in didyoumean.findChildren('a')]
                logger.error("Showing results for %r.", showed)
            #searchresults = soup.find('div', class_='searchresults')
            #searchmatches = searchresults.find_all('div', class_='mw-search-result-heading')
            #searchmatches_extracts = searchresults.find_all('div', class_='searchresult')
            #searchmatches_data = searchresults.find_all('div', class_='mw-search-result-data')
            elif soup.find('p', class_='mw-search-nonefound'):
                logger.error("No results matching the query: %r", term)
            else:
                logger.error("No recognizable content in the fetched html page.")

    return tablesoups


# The recursive hierarchy is the following:
#
# table class="clade"
#   tbody
#     1. tr : the content
#       - List of td: clade children (class "clade", "clade-label" or "clade-leaf")
#           - td of class "clade":
#               - (a link child)
#               - a table child
#           - td of class "clade-label":
#               - if no text and a "border-bottom" in the style attr, it's a branch.
#               - 
#           - td of class "clade-leaf":
#               - could contain a paragraph elem, with link and image.
#               - or a table of class "clade".
#
#     2. tr :
#       - additional info (e.g 690 Mya). Possibly class="clade-slabel".
#       - or sister clade.

# NOTE: just check that inter-html elements (i.e NavigableString) do not
# contain anything else than white spaces.


def find_matching_node(tree, *patterns):
    regex = re.compile(r'|'.join(re.escape(p) for p in patterns), re.I)
    for node in tree.traverse('levelorder'):
        if regex.search(node.name):
            return node


def build_tree(tablesoup, recurs=0, _recurs_count=0):
    tbody = tablesoup.findChild('tbody', recursive=False)  # not findChild
    
    nodes = []

    if tbody is not None:
        for row in tbody.findChildren('tr', recursive=False):
            cell0 = row.findChild('td', recursive=False)
            tagclass = cell0.get('class', [])
            if 'clade-label' in tagclass: 
                cladename = cell0.get_text().strip()
                nodes.append(ete3.TreeNode(name=cladename))
                nodes[-1].add_feature('info', [])
                nodes[-1].add_feature('wikipedia_page_depth', _recurs_count)
                if 'dashed' in cell0.get('style', ()):
                    # This branch is controversial
                    nodes[-1].support = 0.5

                cladeleaf = cell0.find_next_sibling('td', class_='clade-leaf')
                child_clade = cladeleaf.findChild('table', class_='clade', recursive=False)
                if child_clade is not None:
                    for child_node in build_tree(child_clade, recurs, _recurs_count):
                        nodes[-1].add_child(child=child_node)
                else:
                    leafname = cladeleaf.get_text().strip()
                    #not_img = lambda tag: "image" not in tag.get('class', '')
                    leaflink = cladeleaf.find('a')
                    leafimg = cladeleaf.find('img')
                    if not nodes[-1].name:
                        # Update the preceding node, which is actually just the leading branch.
                        nodes[-1].name = leafname
                    else:
                        nodes[-1].add_child(name=leafname)

                    if leaflink:
                        nodes[-1].add_feature('link', leaflink.get('href', ''))
                        otherlinks = [l for l in leaflink.find_next_siblings('a')
                                      if 'image' not in l.get('class', '')]
                        if otherlinks:
                            logger.warning("Alternative leaf links: " + \
                                    ";".join("%r class=%r" % (l.get_text(), l.get('class'))
                                             for l in otherlinks))
                    if leafimg:
                        nodes[-1].add_feature('img', leafimg['src'])
                        nodes[-1].add_feature('imgwidth', leafimg['width'])
                        nodes[-1].add_feature('imgheight', leafimg['height'])

                    if _recurs_count < recurs and leaflink and 'redlink=1' not in leaflink['href']:
                        href = leaflink['href']
                        if not href.startswith('https://'):
                            href = WIKIPEDIA_URL + href
                        hreftreesoups = get_wiki_tree(url=href)
                        if hreftreesoups:
                            logger.info("Recursing into %r from %r", leafname,
                                        tablesoup.find_parent('[document]')\
                                          .findChild('title').get_text().strip()
                                        )
                            logger.info("Found %d phylogenetic trees (at depth %d).",
                                        len(hreftreesoups), _recurs_count+1)
                            leaflinktext = leaflink.get_text().strip()

                            for hreftreesoup in hreftreesoups:
                                for hreftree in build_tree(hreftreesoup, recurs, _recurs_count+1):
                                    matched_node = find_matching_node(hreftree,
                                                                      leaflinktext,
                                                                      *leafname.split('/'))

                                    logger.debug("Matched node: %r", matched_node)
                                    if matched_node and not matched_node.is_leaf():
                                        for leafchild in matched_node.children:
                                            nodes[-1].add_child(child=leafchild)
                                        break
                                else:
                                    continue  # next hreftreesoup if no match
                                break
                            else:
                                logger.warning("Corresponding internal node (%r/%r) not found.",
                                               leafname, leaflinktext)

            elif 'clade-slabel' in tagclass:
                nodes[-1].info.append(cell0.get_text().strip())
            else:
                logger.warning("Unexpected class of cell in the row under %r: %r",
                               (nodes[-1].name if nodes else None), tagclass)
    return nodes


def main(term, outbase=None, outfmt=None, nhx=False, show_img=False, recurs=0):
    graphics_fmts = set(outfmt).intersection(('png', 'jpg', 'svg', 'pdf')) \
                    if outfmt is not None else set()
    
    if graphics_fmts or (not outbase and not outfmt):
        # Define only when the above conditions are verified, so that you
        # can fallback on text methods when PyQt is not installed.
        if show_img:
            #async def?
            def add_img(node):
                nodeimg = getattr(node, 'img', None)
                if nodeimg:
                    if nodeimg.startswith('//'):
                        nodeimg = 'https:' + nodeimg
                    #await ?
                    ete3.add_face_to_node(ete3.ImgFace(nodeimg,
                                                       width=int(node.imgwidth),
                                                       height=int(node.imgheight),
                                                       is_url=True),
                                          node, column=1, position='branch-right')
        else:
            def add_img(node):
                pass

        ns = ete3.NodeStyle(size=0)
        dashed_branch = ete3.NodeStyle(size=0, hz_line_type=1)

        def mylayout(node):
            node.set_style(ns)
            if not node.is_leaf():
                ete3.add_face_to_node(ete3.TextFace(node.name), node, column=0,
                                      position='branch-top')
                ete3.add_face_to_node(ete3.TextFace('\n'.join(getattr(node, 'info', []))),
                                      node, column=0, position='branch-bottom')
            if node.support <= 0.5:
                node.set_style(dashed_branch)
            add_img(node)

    treesoups = get_wiki_tree(term)
    logger.info("Found %d phylogenetic trees", len(treesoups))
    if outfmt:
        if outbase:
            outbase += '-%d'
        outputfuncs = []
        if 'nwk' in outfmt:
            features = ['support', 'info', 'link', 'img'] if nhx else None
            def output(tree, i):
                # format 8: all names
                outfile = (outbase % i + '.nwk') if outbase else None
                txt = tree.write(outfile=outfile, format=8,
                                 format_root_node=True, features=features)
                if txt is not None:
                    print(txt)
            outputfuncs.append(output)
        if 'ascii' in outfmt:
            # Always to stdout
            def output(tree, i):
                print(tree.get_ascii())
            outputfuncs.append(output)
        if graphics_fmts and outbase:
            def output(tree, i):
                for fmt in graphics_fmts:
                    tree.render((outbase % i) + '.' + fmt, mylayout)
            outputfuncs.append(output)
        def outputs(tree, i):
            for outfunc in outputfuncs:
                outfunc(tree, i)
    else:
        def outputs(tree, i):
            tree.show(mylayout, name=('Tree n°%d: %s' %(i, tree.name)))

    for i, treesoup in enumerate(treesoups):
        tree, = build_tree(treesoup, recurs)
        outputs(tree, i)


if __name__ == '__main__':
    
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('term', help='Searched wikipedia page')
    parser.add_argument('outbase', nargs='?',
                        help='Output file basename. If None, display the tree.')
    parser.add_argument('-f', '--outfmt', action='append', default=None,
                        choices=['nwk', 'ascii', 'svg', 'pdf', 'png', 'jpg'],
                        help='[if None, display the tree]')
    parser.add_argument('--nhx', action='store_true',
                        help="Save the 'info', 'link' and 'img' features as NHX comments.")
    parser.add_argument('-i', '--img', action='store_true', dest='show_img',
                        help='Download and display images attached to clades.')
    parser.add_argument('-r', '--recurs', type=int, default=0, metavar='R',
                        help=('Goes %(metavar)s levels down to fetch the '
                              'descendant trees following the links at the '
                              'leaves (experimental) [%(default)s].'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show DEBUG messages [default: INFO].')
    
    args = vars(parser.parse_args())
    if args.pop('verbose'):
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    main(**args)

