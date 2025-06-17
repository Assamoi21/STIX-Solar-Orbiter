"""Microapp basic object module"""

import sys, abc


#exclude_list = ["exec", "eval", "breakpoint", "delattr", "setattr",
#                "globals", "input", "locals", "memoryview", "object",
#                "open", "print", "super", "type", "vars", "__import__"]

exclude_list = ["exec", "eval", "breakpoint", "memoryview"]

microapp_builtins = dict((k, v) for k, v in __builtins__.items()
                       if k not in exclude_list)


#if sys.version_info >= (3, 0):
#    Object = abc.ABCMeta("Object", (object,), {})
#
#else:
#    Object = abc.ABCMeta("Object".encode("utf-8"), (object,), {})
Object = abc.ABCMeta("Object", (object,), {})

# TODO: _microid_

class MicroappObject(Object):

    def __new__(cls, *vargs, **kwargs):

        obj = super(MicroappObject, cls).__new__(cls, *vargs, **kwargs)

        return obj

    @property
    @abc.abstractmethod
    def _name_(self):
        pass

    @property
    @abc.abstractmethod
    def _version_(self):
        pass

del sys, abc
del exclude_list
