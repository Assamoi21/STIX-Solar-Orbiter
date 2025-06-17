# -*- coding: utf-8 -*-
"""Microapp error module"""


class MicroappError(Exception):

    def __init__(self, msg, **kwargs):

        super(MicroappError, self).__init__(msg)


class NormalExit(MicroappError):

    def __init__(self, **kwargs):
        super(NormalExit, self).__init__("", **kwargs)


class UsageError(MicroappError):
    pass


class InternalError(MicroappError):
    pass


class TestError(MicroappError):
    pass


class ConfigError(MicroappError):
    pass


class TypeCheckError(MicroappError):

    def __init__(self, value, expected_type):

        msg = ('type of "%s" is "%s", but expected "%s"' %
               (str(value), type(value), expected_type))
        super(TypeCheckError, self).__init__(msg)


class UnknownNameError(MicroappError):
    pass
