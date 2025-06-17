# -*- coding: utf-8 -*-
"""Microapp manage module"""

from __future__ import print_function


from microapp.error import UsageError


class Manager(object):

    def __new__(cls, *vargs, **kwargs):
        obj = super(Manager, cls).__new__(cls)
        obj._mgrapp = vargs[0] if vargs else obj
        return obj

    def get_builtin_apps(self):
        return self._mgrapp.manager.get_builtin_apps()

    def get_forward(self, key):
        return self._mgrapp._fwds[key]

    def iter_forward(self):
        return list(self._mgrapp._fwds.items())

    def set_forward(self, key, value):
        self._mgrapp._fwds[key] = value

    def get_shared(self, key):
        return self._mgrapp._shrds[key]

    def iter_shared(self):
        return list(self._mgrapp._shrds.items())

    def set_shared(self, key, value):
        self._mgrapp._shrds[key] = value

    def get_downcast(self, key):
        return self._mgrapp._dcasts[key]

    def iter_downcast(self):
        return list(self._mgrapp._dcasts.items())

    def set_downcast(self, key, value):
        self._mgrapp._dcasts[key] = value

    def get_config(self, *vargs, **kwargs):
        return self._mgrapp.get_config(*vargs, **kwargs)

    def set_config(self, *vargs, **kwargs):
        return self._mgrapp.set_config(*vargs, **kwargs)

    def has_config(self, *vargs, **kwargs):
        return self._mgrapp.has_config(*vargs, **kwargs)

    def run_command(self, *vargs, **kwargs):
        return self._mgrapp.manager.run_command(*vargs, **kwargs)

    def run_class(self, cls, *vargs, **kwargs):
        return self._mgrapp.manager.run_class(cls, *vargs, **kwargs)

    def log(self, level, msg, name):

        name.insert(0, self._mgrapp._name_)

        self._mgrapp.manager.log(level, msg, name)
