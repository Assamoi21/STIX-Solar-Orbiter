# -*- coding: utf-8 -*-
"""Microapp utility module"""

from __future__ import print_function

import sys, os, ast, re, functools, logging

pat_id = r"[^\d\W]\w*"
re_mfr_arg = re.compile((r"\s*(?P<var>"+pat_id+r")\s*=\s*(?P<arg>.+)"))

#_re_mfr_arg = re.compile(("\s*(?P<outvar>"+_pat_id+")\s*=\s*(?P<func>"+
#        _pat_id+")\s*of\s*(?P<invar>"+_pat_id+")\s*"))


class appdict(dict):

    def items(self, *vargs, **kwargs):
        return list(super(appdict, self).items(*vargs, **kwargs))

    def keys(self, *vargs, **kwargs):
        return list(super(appdict, self).keys(*vargs, **kwargs))

    def values(self, *vargs, **kwargs):
        return list(super(appdict, self).values(*vargs, **kwargs))


class Logger(object):

    def __init__(self, mgr, name=None):

        self.name = name if name else []
        self.mgr = mgr

    def debug(self, msg):
        self.mgr.log(logging.DEBUG, msg, self.name)

    def info(self, msg):
        self.mgr.log(logging.INFO, msg, self.name)

    def warning(self, msg):
        self.mgr.log(logging.WARNING, msg, self.name)

    def error(self, msg):
        self.mgr.log(logging.ERROR, msg, self.name)

    def critical(self, msg):
        self.mgr.log(logging.CRITICAL, msg, self.name)

    def exception(self, msg):
        self.mgr.log(logging.EXCEPTION, msg, self.name)


def _p(*argv, **kw_str):
    return list(argv), kw_str

def tostr(text, enc=None):

    if isinstance(text, str):
        return text

    if enc is None:
        if sys.stdout.encoding:
            enc = sys.stdout.encoding

        else:
            enc = "utf-8"

    try:
        text = text.decode(enc)

    except:
        text = str(text)

    return text


def tostr_iter(l):

    for e in l:
        yield tostr(e)


def load_pymod(head, base):

    sys.path.insert(0, os.path.abspath(os.path.realpath(head)))
    m = __import__(base)
    sys.path.pop(0)

    return m


def import_modulepath(path):

    fragment = None
    spath = [p.strip() for p in path.rsplit("#", 1)]

    if len(spath) > 1:
        path, fragment = spath

    if os.path.exists(path):
        head, base = os.path.split(path)
        mod = None

        if os.path.isfile(path) and path.endswith(".py"):
            modname = base[:-3]
            mod = load_pymod(head, modname)

        elif (os.path.isdir(path) and
                os.path.isfile(os.path.join(path, "__init__.py"))):
            modname = base[:-1] if base[-1] == os.sep else base
            mod = load_pymod(head, modname)
    else:
        try:
            modname = path
            mod = __import__(modname)

        except ModuleNotFoundError as err:
            from microapp.error import UsageError
            raise UsageError("path does not exist: %s" % path)

    if mod:
        if fragment:
            return fragment, getattr(mod, fragment)

        else:
            return modname, mod


def appeval(text, env):

    if not text or not isinstance(text, str):
        return text

    val = None
    lenv = appdict()

    stmts = ast.parse(text).body

    if len(stmts) == 1 and isinstance(stmts[-1], ast.Expr):
        val = eval(text, env, lenv)

    else:
        exec(text, env, lenv)

    return val, lenv

def type_check(value, typedef):
    # TODO: implement this

    return True


def funcargseval(text, env):

    env["_appeval_p"] = _p
    fargs, out = appeval("_appeval_p(%s)" % text, env)
    del env["_appeval_p"]
    return fargs, out


def assert_check(check):

    data = check(check.data)
    return data["_"]


def _mfr_arg_parse(argtype, fwds):

    m = re_mfr_arg.match(argtype.data)

    if m:
        newenv = appdict(argtype.env)
        newenv["functools"] = functools
        for node, fwd in fwds.items():
            for k, v in fwd.items():
                if k in newenv:
                    if isinstance(v, list):
                        newenv[k].extend(v)
                    else:
                        newenv[k].append(v)

                else:
                    if isinstance(v, list):
                        newenv[k] = v

                    else:
                        newenv[k] = [v]

        if sys.version_info < (3, 0):
            m = m.groupdict()

        return m["var"], m["arg"], newenv

    else:
        from microapp.error import UsageError
        raise UsageError("Wrong argument syntax: " + argtype.data)


def reduce_arg(argtype, fwds):
    var, rargs, newenv = _mfr_arg_parse(argtype, fwds)
    val = eval("functools.reduce(%s)" % rargs, newenv, appdict())
    return appdict({var: val})

def filter_arg(argtype, fwds):
    var, fargs, newenv = _mfr_arg_parse(argtype, fwds)
    val = eval("functools.filter(%s)" % fargs, newenv, appdict())
    return appdict({var: val})

def map_arg(argtype, fwds):
    var, margs, newenv = _mfr_arg_parse(argtype, fwds)
    val = eval("functools.map(%s)" % margs, newenv, appdict())
    return appdict({var: val})


def funcargs(arg):

    return arg if isinstance(arg, (list, tuple)) else [arg], appdict()

