# -*- coding: utf-8 -*-
"""Microapp parallel processing module"""

from __future__ import print_function

import os, abc, enum, functools
from multiprocessing import Pipe

from microapp.base import Object
from microapp.utils import appdict, Logger

class PROXYMSG(enum.Enum):
    C2P_BEG = 1
    C2P_END = 2
    C2P_ERR = 3
    C2P_MSG = 4

    P2C_END = 10
    P2C_MSG = 11


# messge format: id, args, kwargs
#class EdgeProxy(metaclass=abc.ABCMeta):
class EdgeProxy(Object):
    pass


class EdgeParentProxy(EdgeProxy):
    pass


class EdgeChildProxy(EdgeProxy):

    @abc.abstractmethod
    def pack_launch(self, forward, edge, arrival):
        pass

    @abc.abstractmethod
    def process_child_message(self):
        pass

    def process_idletask(self):
        pass


class EdgeMultiprocChildProxy(EdgeChildProxy):

    def __init__(self, mgr):
        self.mgr = mgr

    def _init(self):
        self.ppipe, self.cpipe = None, None
        self.child_pid = None
        self.edge = None
        self.arrival = None

    def pack_launch(self, forward, edge, arrival):

        self._init()

        self.ppipe, self.cpipe = Pipe()
        self.edge = edge
        self.arrival = arrival

        args = edge, forward, self.cpipe
        kwargs = appdict()

        return multiproc_child_perform, args, kwargs

    def process_child_message(self):

        finished = None

        if self.ppipe.poll():
            msgid, vargs, kwargs = recvmsg(self.ppipe)

            if msgid == PROXYMSG.C2P_BEG:
                self.child_pid = kwargs["pid"]

            elif msgid == PROXYMSG.C2P_END:
                child_pid = kwargs["pid"]

                if child_pid == self.child_pid:
                    sendmsg(self.ppipe, PROXYMSG.P2C_END)
                    finished = (self.arrival, self.edge, kwargs["output"],
                                    kwargs["forward"])
                else:
                    pass
                    # TODO: something wrong

            elif msgid == PROXYMSG.C2P_ERR:
                sendmsg(self.ppipe, PROXYMSG.P2C_END)
                finished = (self.arrival, self.edge, -1, kwargs)

            elif msgid == PROXYMSG.C2P_MSG:
                attr = kwargs.pop("_orgattr_")
                result = getattr(self.mgr, attr)(*vargs, **kwargs)
                try:
                    sendmsg(self.ppipe, PROXYMSG.P2C_MSG, result=result)
                except TypeError as err:
                    import pdb; pdb.set_trace()
                    print(err)

        return finished


class EdgeMultiprocParentProxy(EdgeParentProxy):

    def __init__(self, pipe):
        self.pipe = pipe

    def __getattr__(self, attr):
        return functools.partial(communicator, _orgattrname_=attr,
                _pipe_=self.pipe)


def communicator(*vargs, **kwargs):
    attr = kwargs.pop("_orgattrname_")
    pipe = kwargs.pop("_pipe_")

    sendmsg(pipe, PROXYMSG.C2P_MSG, *vargs, _orgattr_=attr, **kwargs)
    msgid, pvargs, pkwargs = recvmsg(pipe)

    if msgid == PROXYMSG.P2C_MSG:
        return pkwargs["result"]

    raise InternalError("Expected '%s' but received '%s'." %
            (PROXYMSG.P2C_MSG, msgid))


def sendmsg(pipe, msgid, *vargs, **kwargs):
    pipe.send((msgid, vargs, kwargs))


def recvmsg(pipe):
    return pipe.recv()


def multiproc_child_perform(edge, fwd, pipe):
    try:

        mgr = EdgeMultiprocParentProxy(pipe)
        logger = Logger(mgr, ["edge-%d" % edge._eid])

        sendmsg(pipe, PROXYMSG.C2P_BEG, pid=os.getpid())

        logger.debug("multiproc edge begin")
        out, fwd = edge.run(mgr, fwd) 
        logger.debug("multiproc edge end")

        sendmsg(pipe, PROXYMSG.C2P_END, pid=os.getpid(),
                    output=out, forward=fwd)

    except Exception as err:
        sendmsg(pipe, PROXYMSG.C2P_ERR, pid=os.getpid(),
                    error=err)

    finally:

        msgid, vargs, kwargs = recvmsg(pipe)

        if msgid != PROXYMSG.P2C_END:
            pass
            # TODO: error callback?

