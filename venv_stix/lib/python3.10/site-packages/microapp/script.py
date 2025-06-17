# -*- coding: utf-8 -*-
"""Microapp scripting module"""

from __future__ import print_function

import os, re, abc, pkg_resources

from microapp.base import microapp_builtins
from microapp.error import UsageError, InternalError
from microapp.utils import appeval, appdict
from microapp.app import App
from microapp.group import Group, AppEdge

re_mascmd = re.compile(r"^(?P<space>\s*)(?P<name>[^@:]+)(:\s*(?P<dtype>[^\d\W]\w*)\s*)?@\s*(?P<context>[^\d\W]\w*)\s*=(?P<rhs>.*)")
re_massec = re.compile(r"^(?P<name>[^@:]+)(:\s*(?P<context>[^\d\W]\w*)\s*)?(@\s*(?P<control>.+))?")

_microapp_section_mods = []
_microapp_command_mods = []


class MicroappSection(App):
    _name_ = "microappsection"
    _version_ = "0.1.0"

    def __new__(cls, mgr, *vargs, **kwargs):

        obj = super(MicroappSection, cls).__new__(cls, mgr, *vargs, **kwargs)

        obj.add_argument("_section_body_", help="section body")

        return obj

    def run(self, name, body): 
        self._section_name_ = name
        self._env["_section_body_"] = body
        return super(MicroappSection, self).run(["@_section_body_"], [], appdict())


class PythonSection(MicroappSection):
    _name_ = "pythonsection"
    _version_ = "0.1.0"

    def perform(self, args):

        body = "\n".join(args._section_body_["_"])
        _, env = appeval(body, self._env)
        self._env.update(env)


class MasEdge(AppEdge):

    def __init__(self, args):
        self._args_ = args

    def append_app(self, app, name, control, body):
        self.append_task((app, name, control, body))

    def run(self, mgr, fwds=None):

        edgefwds = appdict() 

        env = appdict(fwds) if fwds else appdict()
        env.update(self._args_)

        for k, v in env.items():
            if not k.startswith("_"):
                mgr.set_shared(k, v)

        for appcls, name, control, body in self._tasks:

            if not control or appeval(control, env)[0]:

                app = appcls(mgr)
                env.update(app._env)

                ret, fwds = app.run(name, body)

                for k, v in app._env.items():
                    if not k.startswith("_"):
                        mgr.set_shared(k, v)

                edgefwds.update(fwds)

                if ret != 0:
                    raise TestError("'%s' section returned %d." % (name, ret))

        return 0, edgefwds


class Script(Group):

    _name_ = "script"
    _version_ = "0.1.0"

    def read_mas(self, path):

        if not os.path.isfile(path):
            raise UsageError("Given path is not a file: %s" % path)

        curindent = ""
        cursec = ["_entry_"]
        sections = [cursec]

        with open(path) as f:
            lines = []
            merge = False

            for line in f:
                pos = line.rfind("\\")

                if pos>=0:
                    if merge:
                        lines[-1] += line[:pos]

                    else:
                        lines.append(line[:pos])

                    merge = True

                elif merge:
                    lines[-1] += line
                    merge = False

                else:
                    lines.append(line)

            for line in lines:
                lstrip = line.lstrip()
                strip = line.strip()

                if (lstrip and lstrip[0] == "[" and strip[-1] == "]"):
                    curindent = line[:(len(line)-len(lstrip))]
                    cursec = [strip[1:-1].strip()]
                    sections.append(cursec)

                else:
                    pos = line.find(curindent)
                    if pos>=0:
                        cursec.append(line[pos:].rstrip())
                    else:
                        raise UsageError("Section indentation error: %s" %
                            line)
        return sections

    def add_arguments(self, lines):

        remained = lines[:1]

        for line in lines[1:]:
            match = re_mascmd.match(line)
            if match:
                ctx = match.group("context")
                if ctx and ctx.strip()=="arg":
                    name = match.group("name")
                    rhs = match.group("rhs")

                    if name: name.strip()
                    if rhs: rhs.strip()

                    if name and rhs:
                        self.argnames.append(name)
                        appeval("self.add_argument(%s, dest='%s')" %
                                (rhs, name), appdict({'self': self}))
                    else:
                        raise UsageError("Wrong argument command: " + line)

                else:
                    remained.append(line)
            else:
                remained.append(line)

        return remained

    def add_sections(self, sections):
        global _microapp_section_mods

        self._sections = []

        for sec in sections:
            match = re_massec.match(sec[0])
            if match:
                name = match.group("name")
                context = match.group("context")
                control = match.group("control")

                if not _microapp_section_mods:
                    _microapp_section_mods = [ep.load() for ep in
                            pkg_resources.iter_entry_points(
                                group='microapp.sections')]
                    # TODO: update mas_handlers with mods

                handler = mas_handlers.get(context, PythonSection)

                self._sections.append((handler, name, control, sec[1:]))

            else:
                raise UsageError("Mas section header syntax error: %s" % str(sec[0]))
        
    def connect(self, args, subargs):

        argdata = appdict()
        for argname, argvalue in args._get_kwargs():
            if argname in self.argnames:
                if argvalue:
                    argdata[argname] = argvalue["_"]
                else:
                    argdata[argname] = None

        edge = MasEdge(argdata)

        for cls, name, control, body in self._sections:
            edge.append_app(cls, name, control, body)

        if args.clone:
            for idx in range(int(args.clone["_"])):
                self.connect_edge(self, edge.clone(), self)
        else:
            self.connect_edge(self, edge, self)


def make_class(maspath):

    clsstr = """
class MicroappScript(Script):

    def __init__(self, mgr):

        self.argnames = []
        self._name_ = os.path.basename(maspath)

        sections = self.read_mas(maspath)
        remained = self.add_arguments(sections[0])

        if remained:
            self.add_sections([remained]+sections[1:])

        else:
            self.add_sections(sections[1:])
"""

    env = appdict({"__builtins__": appdict(microapp_builtins), "_": None,
            "Script": Script, "maspath": maspath, "os": os})
    lenv = appdict()

    exec(clsstr, env, lenv)

    return lenv["MicroappScript"]


mas_handlers = appdict({
    "python" : PythonSection
})
