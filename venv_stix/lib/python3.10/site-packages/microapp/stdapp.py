# -*- coding: utf-8 -*-
"""Microapp standard application module"""

from __future__ import print_function

import os

from typing import Any

from microapp.app import App
from microapp.utils import appdict, tostr

class InputApp(App):
    """Input application"""

    _name_ = "input"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("data", nargs="*", help="input data")

        self.register_forward("data", type=Any, help="forward input data")

    def perform(self, args):

        data = [d["_"] for d in args.data]
        self.add_forward(data=data)


class PrintApp(App):
    """print application"""

    _name_ = "print"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("data", nargs="*", help="data to print")
        self.add_argument("-s", "--strip", action="store_true",
                                 help="strip newline at the end of string")

        self.register_forward("stdout", type=str, help="standard output")

    def perform(self, args):

        end = "" if args.strip else "\n"

        if args.data:
            l = []

            for d in args.data:
                val = d['_']
                if not isinstance(val, str):
                    l.append(str(val))
                else:
                    l.append(val)

            out = "".join(l)
            print(out, end=end)
            stdout = out + end
        else:
            stdout = "No data to print."
            print(stdout)

        self.add_forward(stdout=stdout)


class ShellApp(App):
    """Shell command application"""

    _name_ = "shell"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("command", help="shell command to run")
        self.add_argument("--workdir", help="working directory")
        self.add_argument("--useenv", action="store_true", help="use shell environment")
        self.add_argument("--env", action="append", help="add shell environment")
        self.add_argument("--show-stdout", action="store_true", help="show stdout on display")
        self.add_argument("--show-stderr", action="store_true", help="show stderr on display")

        self.register_forward("retcode", type=str, help="return code")
        self.register_forward("stdout", type=str, help="standard output")
        self.register_forward("stderr", type=str, help="standard output")

    def perform(self, args):

        import sys
        from shlex import split
        from subprocess import Popen, PIPE
       
        cmds = args.command["_"] if args.useenv else split(args.command["_"])
        workdir = args.workdir["_"] if args.workdir else None

        proc = Popen(cmds, stdout=PIPE, stderr=PIPE, shell=args.useenv, cwd=workdir)

        outs, errs = proc.communicate()

        retcode = proc.poll()

        if args.show_stdout:
            print(tostr(outs), file=sys.stdout)

        if args.show_stderr:
            print(tostr(errs), file=sys.stderr)

        self.add_forward(retcode=retcode, stdout=outs, stderr=errs)

        return retcode


standard_apps = appdict({
    InputApp._name_ : InputApp,
    PrintApp._name_ : PrintApp,
    ShellApp._name_ : ShellApp
})
