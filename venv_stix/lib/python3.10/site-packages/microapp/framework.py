# -*- coding: utf-8 -*-
"""Microapp framework module"""

from __future__ import print_function

import os
import pkg_resources

from microapp.error import UsageError
from microapp.utils import import_modulepath, appdict
from microapp.app import App
from microapp.mgmtapp import mgmt_apps
from microapp.stdapp import standard_apps

_microapp_apps = []
_microapp_app_mods = []
_microapp_prj_mods = []

def register_appclass(app):
    global _microapp_apps

    for a in _microapp_apps:
        if app._name_ == a._name_:
            raise UsageError("App name '%s' already exists" % a._name_)

    _microapp_apps.append(app)


def unregister_appclass(app):
    global _microapp_apps
    _microapp_apps.remove(app)


def is_appclass_registered(app):
    return app in _microapp_apps


def load_appclass(apppath, args=None, subargs=None, dapps=None):
    """select application class


    dapps : builtin applications
"""
    global _microapp_app_mods, _microapp_prj_mods

    if not apppath:
        return None, None, None, None


    if not args:
        args = []

    if not subargs:
        subargs = []

    if not dapps:
        dapps = []

    if isinstance(apppath, type):
        if issubclass(apppath, App):
            return apppath, args, subargs, appdict()

        raise UsageError("Not compatible app type: %s" % type(apppath))

    objs = appdict()

    while "--import" in args:
        idx = args.index("--import")
        mpath = args.pop(idx+1)
        args.pop(idx)

        key, obj = import_modulepath(mpath)
        objs[key] = obj


    app_class = None

    fragment = ""
    spath = [p.strip() for p in apppath.rsplit("#", 1)]

    if len(spath) > 1:
        apppath, fragment = spath

    mods = []

    if os.path.exists(apppath):

        if os.path.isfile(apppath):
            head, base = os.path.split(apppath)

            if base.endswith(".mas"):
                from microapp.script import make_class
                app_class = make_class(apppath)

            else:
                try:
                    _, mod = import_modulepath(apppath)
                    if mod:
                        mods.append(mod)
                except TypeError:
                    raise UsageError("Given path is not a valid "
                                     "Python module: " + apppath)

        elif os.path.isdir(apppath):
            try:
                modname, mod = import_modulepath(apppath)
                mods.append(mod)
            except:
                pass

    if app_class is None:
        for app in _microapp_apps: 
            if apppath == app._name_:
                app_class = app
                break

    if app_class is None:

        if not _microapp_app_mods:
            _microapp_app_mods = [ep.load() for ep in
                    pkg_resources.iter_entry_points(group='microapp.apps')]

        for app_mod in _microapp_app_mods: 
            if apppath == app_mod.__name__:
                mods.append(app_mod)
#                from microapp.plxtask import PlXTask
#                if app_class is PlXTask:
#                    task_mod = pyloco_import(apppath)
#                    task_dir = os.path.dirname(task_mod.__file__)
#                    argv.insert(0, os.path.join(task_dir, getattr(task_mod, "plx")))
                break


    candidates = appdict()

    for mod in mods:
        for name in dir(mod):
            if not name.startswith("_"):
                obj = getattr(mod, name)

                if (type(obj) == type(App) and issubclass(obj, App) and
                        (obj.__module__ is None or
                            not obj.__module__.startswith("microapp."))):
                   candidates[name] = obj
    if candidates:
        if fragment:
            if fragment in candidates:
                app_class = candidates[fragment]

            else:
                raise UsageError("No app is found with a fragment of "
                                 "'%s'." % fragment)
        elif len(candidates) == 1:
            app_class = candidates.popitem()[1]

        else:
            raise UsageError(
                "More than one frame are found."
                "Please add fragment to select one: %s" %
                list(candidates.keys())
            )

    if app_class:
        setattr(app_class, "_path_", os.path.abspath(apppath))

    if app_class is None:
        for dapp in dapps:
            if apppath == dapp._name_:
                app_class = dapp
                break

    if app_class is None:
        from microapp.project import Project
        found = False

        if not _microapp_prj_mods:
            _microapp_prj_mods = [ep.load() for ep in
                    pkg_resources.iter_entry_points(group='microapp.projects')]

        for prj_mod in _microapp_prj_mods:
            for name in dir(prj_mod):
                if not name.startswith("_"):
                    obj = getattr(prj_mod, name)
                    if (type(obj)==type(Project) and issubclass(obj, Project) and
                        hasattr(obj, "_builtin_apps_")):
                        for defapp in getattr(obj, "_builtin_apps_"):
                            if defapp._name_ == apppath:
                                app_class = defapp
                                found = True
                                break
                if found:
                    break
            if found:
                break

    if not app_class:
        if apppath in mgmt_apps:
            app_class = mgmt_apps[apppath]

        elif apppath in standard_apps:
            app_class = standard_apps[apppath]

        if not app_class:
            raise UsageError("'%s' app is not found. Please check path." % apppath)

    return app_class, args, subargs, objs
