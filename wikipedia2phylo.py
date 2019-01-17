#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Parse phylogenetic trees from a wikipedia page, and save it in a
bioinformatics format (e.g newick)"""


import os.path as op
import string
import requests
import bs4
import ete3
import argparse as ap

import logging
logger = logging.getLogger(__name__)
logging.basicConfig()

WIKIPEDIA_URL = 'https://en.wikipedia.org'
SEARCH_URL = WIKIPEDIA_URL + '/w/index.php'


def get_wiki_tree(term=None, url=SEARCH_URL):
    page = requests.get(url, params=({'search': term} if term else None))

    if page.status_code != 200:
        logger.error("Status of request = %d", page.status_code)
        
    soup = bs4.BeautifulSoup(page.text, 'lxml')
    page.close()

    #TODO: check if later matches in the document
    tablesoup0 = soup.find('table', class_='clade')
    if tablesoup0 is not None:
        return [tablesoup0] + tablesoup0.find_next_siblings('table', class_='clade')
    else:
        return []


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


def build_tree(tablesoup, recurs=False):
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
                if 'dashed' in cell0['style']:
                    # This branch is controversial
                    nodes[-1].support = 0.5

                cladeleaf = cell0.find_next_sibling('td', class_='clade-leaf')
                child_clade = cladeleaf.findChild('table', class_='clade', recursive=False)
                if child_clade is not None:
                    for child_node in build_tree(child_clade, recurs=recurs):
                        nodes[-1].add_child(child=child_node)
                else:
                    leafname = cladeleaf.get_text().strip()
                    leaflink = cladeleaf.find('a')
                    leafimg = cladeleaf.find('img')
                    if not nodes[-1].name:
                        # Update the preceding node, which is actually just the leading branch.
                        nodes[-1].name = leafname
                    else:
                        nodes[-1].add_child(name=leafname)

                    if leaflink:
                        nodes[-1].add_feature('link', leaflink['href'])
                    if leafimg:
                        nodes[-1].add_feature('img', leafimg['src'])
                        nodes[-1].add_feature('imgwidth', leafimg['width'])
                        nodes[-1].add_feature('imgheight', leafimg['height'])

                    if recurs and leaflink:
                        href = leaflink['href']
                        if not href.startswith('https://'):
                            href = WIKIPEDIA_URL + href
                        hreftree = build_tree(get_wiki_tree(url=href)[0])
                        for leaftree in hreftree.iter_search_nodes(name=leafname):
                            for leafchild in leaftree.children:
                                nodes[-1].add_child(name=leafchild, child=leaftree)
                            break

            elif 'clade-slabel' in tagclass:
                nodes[-1].info.append(cell0.get_text().strip())
            else:
                logger.warning("Unexpected class of cell in the row under %r: %r",
                               (nodes[-1].name if nodes else None), tagclass)
    return nodes


def main(term, outfile=None, nhx=False, show_img=False):
    if show_img:
        def add_img(node):
            nodeimg = getattr(node, 'img', None)
            if nodeimg:
                if nodeimg.startswith('//'):
                    nodeimg = 'https:' + nodeimg
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
    if outfile:
        outbase, outext = op.splitext(outfile)
        outbase += '-%d'
        outfile = outbase + outext
        if outext.lower() in ('.nwk', '.newick', '.nw', '.tree', '.txt'):
            features = ['info', 'link', 'img'] if nhx else None
            def output(tree, i):
                # format 8: all names
                tree.write(outfile=(outfile % i), format=8, format_root_node=True,
                           features=features)
        elif outext.lower() in ('.png', '.jpg', '.jpeg', '.gif', '.svg', '.pdf'):
            def output(tree, i):
                tree.render((outfile % i), mylayout)
        else:
            raise ValueError('Invalid output format %r' % ext)
    else:
        def output(tree, i):
            tree.show(mylayout, name=('Tree n°%d: %s' %(i, tree.name)))

    for i, treesoup in enumerate(treesoups):
        tree, = build_tree(treesoup)
        output(tree, i)


if __name__ == '__main__':
    logger.setLevel(logging.INFO)
    
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('term', help='Searched wikipedia page')
    parser.add_argument('outfile', nargs='?',
                        help='Output file (newick). If None, display the tree.')
    parser.add_argument('--nhx', action='store_true',
                        help="Save the 'info', 'link' and 'img' features as NHX comments.")
    parser.add_argument('-i', '--img', action='store_true', dest='show_img',
                        help='Download and display images attached to clades.')
    args = parser.parse_args()
    main(**vars(args))

