# -*- coding: utf-8 -*-
"""Microapp argument parser module"""

from __future__ import print_function

import sys, re, typing, inspect

from argparse import ArgumentParser, SUPPRESS
from microapp.base import microapp_builtins
from microapp.utils import appeval, funcargseval, type_check, appdict
from microapp.error import UsageError, InternalError

re_context = re.compile(r"^(?P<name>[^\d\W]\w*)\s*:", re.UNICODE)
_typing_members = [m[0] for m in inspect.getmembers(typing, inspect.isclass)]


def data_transfer(arg, buf):
    if arg:
        for a in arg:
            if isinstance(a, ArgType):
                data = a(a.data)
                for dest, value in data.items():
                    if dest != "_":
                        buf[dest] = value
            else:
                raise InternalError("data argument is not ArgType")

def mgrdata_transfer(arg, func):
    if arg:
        for a in arg:
            if isinstance(a, ArgType):
                data = a(a.data)
                for dest, value in data.items():
                    if dest != "_":
                        func(dest, value)
            else:
                raise InternalError("data argument is not ArgType")


class DataRegistry(object):

    def __new__(cls, *vargs, **kwargs):

        obj = super(DataRegistry, cls).__new__(cls)

        obj._fdefs = appdict()
        obj._sdefs = appdict()
        obj._ddefs = appdict()

        obj._fwds = appdict()
        obj._shrds = appdict()
        obj._dcasts = appdict()

        return obj

    def _add_data(self, data, defs, db):

        for dest, value in data.items():
            if dest not in defs:
                raise UsageError("'%s' is not registered for data transfer." %
                                 dest)

            if type_check(value, defs[dest].get("type", typing.Any)):
                db[dest] = value

            else:
                raise TestError("App data type check failure: %s" % dest)


    def register_forward(self, name, **kwargs):
        self._fdefs[name] = kwargs

    def register_shared(self, name, **kwargs):
        self._sdefs[name] = kwargs

    def register_downcast(self, name, **kwargs):
        self._ddefs[name] = kwargs

    def add_forward(self, **kwargs):
        self._add_data(kwargs, self._fdefs, self._fwds)

    def add_shared(self, **kwargs):
        self._add_data(kwargs, self._sdefs, self._shrds)

    def add_downcast(self, **kwargs):
        self._add_data(kwargs, self._ddefs, self._dcasts)


class ArgType(object):

    def __init__(self, argtype, delay, forceeval, env):
        self.argtype = argtype
        self.delay = delay
        self.forceeval = forceeval
        self.env = env
        self.data = None

    def __call__(self, data):

        if self.delay:
            self.delay = False
            self.data = data

            return self

        if not isinstance(data, str):
            return appdict({'_': data})

        # TODO: add where data came from, cmdline or fwd?

        out = data.strip()
        context = None

        # context ...
        match = re_context.match(out)
        if match:
            context = match.group("name")
            out = out[match.span()[1]:]

        if self.forceeval:
            if out and out[0] != "@":
                out = "@" + out

        if len(out)>1 and out[0]=="@":
            if len(out)>2 and out[1] == "<" and out[-1] == ">":
                if out.startswith("@<<") and out.endswith(">>"):
                    val, out = appeval(out[2:-1], self.env)
                else:
                    val, out = funcargseval(out[2:-1], self.env)

            else:
                val, out = appeval(out[1:], self.env)

        else:
            val = out
            out = appdict()

        # convert to argtype
        if (self.argtype and self.argtype.__class__.__name__ not in
                _typing_members):
            val = self.argtype(val)

        out['_'] = val
        out['_context_'] = context

        # TODO: add where data came from, cmdline or fwd?

        return out


class ArgParser(object):

    def __new__(cls, *vargs, **kwargs):

        obj = super(ArgParser, cls).__new__(cls)

        prog = kwargs.pop("prog", sys.argv[0])
        usage = kwargs.pop("usage", None)
        desc = kwargs.pop("description", None)
        epilog = kwargs.pop("epilog", None)

        if sys.version_info >= (3, 0):
            obj._arg_parser = ArgumentParser(prog=prog, usage=usage,
                description=desc, epilog=epilog, allow_abbrev=False)

        else:
            obj._arg_parser = ArgumentParser(prog=prog, usage=usage,
                description=desc, epilog=epilog)

        obj._argtypes = []
        obj._env = appdict({"__builtins__": appdict(microapp_builtins), "_": None})

        return obj

    def set_progname(self, name):
        self._arg_parser.prog = name

    def format_help(self):

        return self._arg_parser.format_help()

    def format_usage(self):

        return self._arg_parser.format_usage()

    def format_option(self):

        orgusage = self._arg_parser.format_usage()
        orghelp = self._arg_parser.format_help()
        return orghelp[len(orgusage)+1:]

    def format_data(self):

        def _get_defs(defs, maxinfo):
            out = []
            for n, info in defs.items():
                t, h = "(type=%s)"%str(info.get("type", "any")), info.get("help", "")
                out.append((n, t, h))
                maxinfo[0] = max(maxinfo[0], len(n))
                maxinfo[1] = max(maxinfo[1], maxinfo[0]+len(t))
            return out

        def _format_defs(defs, maxinfo):

            out = []
            for n, t, h in defs:
                nsp1 = maxinfo[0] - len(n) + 1
                nsp2 = maxinfo[1] - len(t)
                out.append("  " + n + " "*nsp1 + t + " "*nsp2 + h)

            return "\n".join(out)

        if isinstance(self, DataRegistry):
            maxinfo = [0, 0]

            fwds = _get_defs(self._fdefs, maxinfo)
            shrds = _get_defs(self._sdefs, maxinfo)
            dcasts = _get_defs(self._ddefs, maxinfo)

            out = []

            if fwds:
                out.append("\nThis app may feed-forward following data to next app:")
                out.append(_format_defs(fwds, maxinfo))

            if shrds:
                out.append("\nThis app may share following data between sibling apps:")
                out.append(_format_defs(shrds, maxinfo))

            if dcasts:
                out.append("\nThis app may downcast following data to all underlying apps:")
                out.append(_format_defs(dcasts, maxinfo))

            return "\n".join(out)

        else:
            return ""

    def parse_known_args(self, args, env):

        return self._arg_parser.parse_known_args(args)

    def parse_args(self, args, env):

        return self._arg_parser.parse_args(args)

    def add_argument(self, *vargs, **kwargs):

        # TODO: show syshelp if requested
        syshelp = kwargs.pop("syshelp", None)
        if syshelp is not None:
            kwargs["help"] = SUPPRESS

        orgtype = kwargs.pop("type", None)
        delay = kwargs.pop("delay", False)
        forceeval = kwargs.pop("eval", False)
        argtype = ArgType(orgtype, delay, forceeval, self._env)

        try:
            return self._arg_parser.add_argument(*vargs, type=argtype, **kwargs)

        except TypeError:
            return self._arg_parser.add_argument(*vargs, **kwargs)


class ProjectArgParser(ArgParser):

    def format_usage(self):

        orghelp = super(ProjectArgParser, self).format_help()
        posdesc = orghelp.find(self._description_)
        posupper = orghelp[:posdesc-2].rfind("\n")

        appcmd = "[[--] <app> [app-args]] [-- <app> [app-args]]...\n\n" 

        if posupper > 0:
            temp = orghelp[posupper+1:]
            nsp = len(temp) - len(temp.lstrip())

        else:
            nsp = orghelp.find("[")

        return orghelp[:posdesc-1] + " "*nsp + appcmd

    def format_help(self):

        return self.format_usage() + self.format_option()


class GroupArgParser(ArgParser):
    pass


class AppArgParser(ArgParser):
    pass

