"Microapp gzip module wrapper for uncompressing"

import os
import gzip
import shutil

from typing import Any

from microapp import App


class Gunzip(App):

    _name_ = "gunzip"
    _version_ = "0.1.10"
    _description_ = "Microapp gunzip"
    _long_description_ = "Microapp gunzip"
    _author_ = "Youngsung Kim"
    _author_email_ = "youngsung.kim.act2@gmail.com"
    _url_ = "https://github.com/grnydawn/micro-gunzip"

    def __init__(self, mgr):

        self.add_argument("zipfile", type=str, help="zipped file")
        self.add_argument("-o", "--outfile", type=str, help="file path for a unzipped file")

        self.register_forward("data", type=Any, help="unzipped binary data")

    def perform(self, args):

        zipfile = args.zipfile["_"]

        with gzip.open(zipfile, 'rb') as f_in:
            if args.outfile:
                with open(args.outfile["_"], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            #else:
            #    print(f_in.read())
                #root, ext = os.path.splitext(zipfile)
                #outfile = root if ext == ".gz" else "%s.uncompressed" % zipfile

                #with open(outfile, 'wb') as f_out:
                #    shutil.copyfileobj(f_in, f_out)

            f_in.seek(0)
            self.add_forward(data=f_in.read().decode())
