# -*- coding: utf-8 -*-
"""Microapp management application module"""

from __future__ import print_function

import sys, os, pkg_resources, tempfile, subprocess, shutil

from typing import Any
from microapp.error import UsageError
from microapp.utils import appdict
from microapp.app import App



_setup_template = """
from setuptools import setup

setup(
    name= "{name}",
    version="{version}",
    description="{description}",
    packages=[{packages}],
    install_requires=[{install_requires}],
    entry_points={{"microapp.{plugintype}": "{name} = {name}"}}
)
"""

_init_template = """
from .{module_name} import {class_name}
"""

# copied from distutils.core.py
setup_keywords = ('distclass', 'script_name', 'script_args', 'options',
                  'name', 'version', 'author', 'author_email',
                  'maintainer', 'maintainer_email', 'url', 'license',
                  'description', 'long_description', 'keywords',
                  'platforms', 'classifiers', 'download_url',
                  'requires', 'provides', 'obsoletes',
                  )


class HelpApp(App):
    """print help message of a Microapp app"""

    _name_ = "help"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("app", nargs="+", help="app name")

        self.register_forward("data", type=Any, help="Microapp apps")

    def print_help(self, objcls, path):

        #prjname = os.path.split(sys.argv[0])[1]
        prjname = self.get_config("project.name")
        mgr = self.get_manager()

        from microapp.group import Group
        from microapp.script import Script
        from microapp.project import Project

        class D(objcls):
            def __new__(cls, m, *vargs, **kwargs):
                if "prog" not in kwargs:
                    if issubclass(objcls, App):
                        name = cls._name_

                    elif issubclass(objcls, Script):
                        script = objcls(mgr)
                        name = script._name_

                    elif issubclass(objcls, Group):
                        name = cls._name_

                    elif issubclass(objcls, Project):
                        raise NotImplementedError("Project object is not supported yet.")

                    else:
                        import pdb; pdb.set_trace()
                    kwargs["prog"] = "%s-%s" % (prjname, name)
                return super(D, cls).__new__(cls, m, *vargs, **kwargs)
        obj = D(mgr)

        # usage
        print(obj.format_usage(), end="")

        # options
        print("\n"+obj.format_option(), end="")

        # data movement
        dfmt = obj.format_data()
        if dfmt:
            print(dfmt, end="")

        # examples
        #exam = obj.format_example()
        #if exam:
        #    print("\n"+exam, end="")

        print("\n")

    def perform(self, args):

        from microapp.framework import load_appclass

        for app in args.app:
            try:
                appcls, _, _, _ = load_appclass(app["_"])
                self.print_help(appcls, app["_"])

            except UsageError:

                from microapp.project import Project

                found = False
                for ep in pkg_resources.iter_entry_points(group='microapp.projects'):
                    prj_mod = ep.load()
                    for name in dir(prj_mod):
                        obj = getattr(prj_mod, name)
                        if type(obj)==type(Project) and issubclass(obj, Project):
                            if obj._name_ == app["_"]:
                                self.print_help(obj, app["_"])
                                found = True
                                break
                    if found: break


class ListApp(App):
    """print a list of locally installed Microapp apps"""

    _name_ = "list"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("data", nargs="*", help="app categories")
        self.add_argument("-v", "--verbose", action="store_true",
                            help="verbose display")

        self.register_forward("data", type=Any, help="locally installed apps")


    def _objinfo(self, app, maxinfo):

        maxdesclen = 75
        desc = ""

        if hasattr(app, "_description_"):
            desc = getattr(app, "_description_") 

        elif app.__doc__:
            pos = app.__doc__.find("\n")

            if pos > 0:
                desc = app.__doc__[:pos]

            elif len(app.__doc__) > maxdesclen:
                desc = app.__doc__[:maxdesclen-3] + "..."

            else:
                desc = app.__doc__

        info = (app._name_, app._version_, desc)

        for i in range(len(maxinfo)):
            maxinfo[i] = max(maxinfo[i], len(info[i]))

        return info

    def print_merged_list(self, prjinfo, appinfo, maxinfo):

        nsp = maxinfo[0] + maxinfo[1] + 3

        if appinfo:
            print("\n%s apps" % self.manager._name_)
            print("-----------------")

            for head, info in appinfo.items():

                if info:
                    for name, ver, desc in info:
                        nleft = len(name) + len(ver) + 2
                        print("%s%s%s : %s" % (name, " "*(nsp-nleft), ver, desc))

        if prjinfo:
            print("\n%s projects" % self.manager._name_)
            print("-----------------")

            for name, ver, desc in prjinfo:
                nleft = len(name) + len(ver) + 2
                print("%s%s%s : %s" % (name, " "*(nsp-nleft), ver, desc))
                #print("%s (%s)%s: %s" % (name, ver, " "*(nsp-nleft), desc))

    def print_list(self, prjinfo, appinfo, maxinfo):

        nsp = maxinfo[0] + maxinfo[1] + 3

        print("\nlocally installed %s apps" % self.manager._name_)

        for head, info in appinfo.items():
            print("\n"+head)
            print("-"*len(head))

            if info:
                for name, ver, desc in info:
                    nleft = len(name) + len(ver) + 2
                    print("%s%s%s : %s" % (name, " "*(nsp-nleft), ver, desc))
            else:
                print("  N/A")

        print("\nlocally installed %s projects" % self.manager._name_)
        print("------------------------------------")

        if prjinfo:
            for name, ver, desc in prjinfo:
                nleft = len(name) + len(ver) + 2
                print("%s%s%s : %s" % (name, " "*(nsp-nleft), ver, desc))
                #print("%s (%s)%s: %s" % (name, ver, " "*(nsp-nleft), desc))
        else:
            print("  N/A")


    def perform(self, args):
 
        apps = appdict()
        prjs = []

        installed = []

        maxinfo = [0, 0, 0]

        apps["downloaded apps"] = installed

        for ep in pkg_resources.iter_entry_points(group='microapp.apps'):
            mod = ep.load()
            for name in dir(mod):
                if not name.startswith("_"):
                    obj = getattr(mod, name)
                    if (type(obj) == type(App) and issubclass(obj, App) and
                            (obj.__module__ is None or
                                not obj.__module__.startswith("microapp."))):
                        installed.append(self._objinfo(obj, maxinfo))


        from microapp.project import Project

        for ep in pkg_resources.iter_entry_points(group='microapp.projects'):
            prj_mod = ep.load()
            for name in dir(prj_mod):
                if not name.startswith("_"):
                    obj = getattr(prj_mod, name)
                    if type(obj)==type(Project) and issubclass(obj, Project):
                        prjs.append(self._objinfo(obj, maxinfo))

                        if hasattr(obj, "_builtin_apps_"):

                            bapps = []
                            apps["'%s' project builtin apps" % name] = bapps

                            for defapp in getattr(obj, "_builtin_apps_"):
                                bapps.append(self._objinfo(defapp, maxinfo))

        standard = []
        apps["Microapp standard apps"] = standard

        for app in mgmt_apps.values(): 
            standard.append(self._objinfo(app, maxinfo))

        from microapp.stdapp import standard_apps

        for app in standard_apps.values(): 
            standard.append(self._objinfo(app, maxinfo))

        if args.verbose:
            self.print_list(prjs, apps, maxinfo)

        else:
            self.print_merged_list(prjs, apps, maxinfo)


        self.add_forward(data=apps)


class ConfigApp(App):
    """read and write microapp configuration data"""

    _name_ = "config"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("param", nargs="*", metavar="param[=value]", help="set or get configuration item")

    def perform(self, args):

        if args.param:
            for item in args.param:
                param = [i.strip() for i in item["_"].rsplit("=", 1)]

                if len(param) == 2:
                    self.manager.set_config(param[0], param[1], createall=True)

                cfg = self.manager.get_config(param[0])
                if isinstance(cfg, (dict, appdict)):
                    print("'%s' includes : %s" % (param[0], ", ".join(cfg.keys())))

                else:
                    print(cfg)
        else:
            print("Top-level config names : global, app, project")


class RegisterApp(App):
    """register Microapp app to an indexing server"""

    _name_ = "register"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("app", nargs="*", help="Microapp apps")

    def perform(self, args):
        pass
        # uploading app as well as project

        # requirements
        # - input output typing
        # - documentation
        # - testing
        # - input data
        # - interaction with other apps?

        # test

        # build
        #import distutils.core as distcore

        #if hasattr(distcore, "setup_keywords"):
        #    setup_keywords = getattr(distcore, "setup_keywords")
            
        #setup_kwargs = appdict()
         
        #import pdb; pdb.set_trace()
        # upload



# TODO: check if app is local file or director
# TODO: if not, check if remote app is microaopp or not
class InstallApp(App):
    """install Microapp app from an indexing server"""

    _name_ = "install"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("app", nargs="*", help="Microapp apps")

        self.register_forward("data", help="installed apps")

    def perform(self, args):

        installed_apps = []

        for app in args.app:

            try:
                from microapp.framework import load_appclass

                spath = [p.strip() for p in app["_"].rsplit("#", 1)]

                if len(spath) == 2:
                    path, frag = spath

                else:
                    path, frag = spath[0], ""

                appcls, _, _, _ = load_appclass(app["_"])
 
                if os.path.isfile(path):

                    with tempfile.TemporaryDirectory() as tempdir:
                        pkgname = appcls._name_
                        desc = appcls.__doc__ if appcls.__doc__ else "Microapp " + pkgname

                        setupargs = {
                            "name": appcls._name_,
                            "version": appcls._version_,
                            "description": desc,
                            "packages": '"%s"' % pkgname,
                            "install_requires": "",
                            "plugintype": "apps"
                        }

                        with open(os.path.join(tempdir, "setup.py"), "w") as sf:
                            sf.write(_setup_template.format(**setupargs))

                        pkgdir = os.path.join(tempdir, pkgname)
                        os.makedirs(pkgdir)

                        filename = os.path.basename(path)
                        modname, ext = os.path.splitext(filename)

                        initargs = {
                            "module_name": modname,
                            "class_name": appcls.__name__
                        }

                        with open(os.path.join(pkgdir, "__init__.py"), "w") as pf:
                            pf.write(_init_template.format(**initargs))

                        shutil.copy(path, pkgdir)

                        # install app
                        pipcmd = [sys.executable, "-m", "pip", "-q",
                                  "--disable-pip-version-check", "install",
                                  "--disable-pip-version-check", tempdir]
                        subprocess.run(pipcmd, cwd=tempdir)

                        installed_apps.append(appcls._name_)

                else:
                        # install app
                        pipcmd = [sys.executable, "-m", "pip", "-q", 
                                  "--disable-pip-version-check", "install",
                                  "--disable-pip-version-check", path]
                        subprocess.run(pipcmd)

            except UsageError as err:
                raise
                # TODO support remote install

        self.add_forward(data=installed_apps)


class UninstallApp(App):
    """uninstall Microapp app"""

    _name_ = "uninstall"
    _version_ = "0.1.0"

    def __init__(self, mgr):

        self.add_argument("app", nargs="*", help="Microapp apps")

        self.register_forward("data", help="uninstalled apps")

    def perform(self, args):

        uninstalled_apps = []

        for app in args.app:

            pipcmd = [sys.executable, "-m", "pip", "--disable-pip-version-check",
                      "-q", "uninstall", "--yes", app["_"]]

            subprocess.run(pipcmd)

            uninstalled_apps.append(app)

        self.add_forward(data=uninstalled_apps)


mgmt_apps = appdict({
    HelpApp._name_ : HelpApp,
    ConfigApp._name_ : ConfigApp,
    InstallApp._name_ : InstallApp,
    UninstallApp._name_ : UninstallApp,
    RegisterApp._name_ : RegisterApp,
    ListApp._name_ : ListApp,
})
