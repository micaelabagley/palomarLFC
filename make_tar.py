#! /usr/bin/env python
import os
import re
import shutil
import ConfigParser
import argparse

def make_tar(wispfield):
    # make directory for final products
    ID = re.search('\d+', wispfield)
    parID = ID.group(0)
    if not os.path.isdir(os.path.join(wispfield,'Par%s'%parID)):
        os.mkdir(os.path.join(wispfield, 'Par%s'%parID))

    # get list of files to copy to directory
    Config = ConfigParser.ConfigParser()
    Config.read('tarlist')
    options = Config.options('tarfile')
    files = []
    for option in options:
        files.append(Config.get('tarfile',option))
    files = [wispfield+x if x[0]=='_' else x for x in files]
    files = [os.path.join(wispfield, f) for f in files]

    for f in files:
        if os.path.exists(f):
            shutil.copy(f, os.path.join(wispfield,'Par%s'%parID))

    # create archive
    archive_name = os.path.join('/home/astro/bagley', 'public_html/wisps', 
                                'palomar_reduction', 'Par%s'%parID)
    shutil.make_archive(archive_name, 'gztar', 
                        os.path.join(wispfield,'Par%s'%parID))


def main():
    parser = argparse.ArgumentParser(description=
        'Create final tar-file for wispfield')
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to construct the tar-file. ' +\
             'Must match name of directory containing all relevant files')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')
    make_tar(wispfield)

if __name__ == '__main__':
    main()
