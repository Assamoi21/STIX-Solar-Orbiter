# -*- coding: utf-8 -*-
"""Microapp group module"""

from __future__ import print_function

import sys, abc

from multiprocessing import cpu_count, Pool
from microapp.base import MicroappObject, microapp_builtins
from microapp.error import TestError, UsageError, InternalError
from microapp.parse import GroupArgParser, ArgType
from microapp.framework import load_appclass
from microapp.utils import reduce_arg, assert_check, appdict, Logger, Logger
from microapp.parallel import EdgeMultiprocChildProxy as EdgeProxy

PY3 = sys.version_info >= (3, 0)


class Edge(object):

    _num_edge_instances = 0

    def __new__(cls, *vargs, **kwargs):

        obj = super(Edge, cls).__new__(cls)

        obj._eid = cls._num_edge_instances # edge id
        obj._bid = None # branch  id
        obj._tasks = []

        cls._num_edge_instances += 1

        return obj

    def append_task(self, task):
        self._tasks.append(task)

    def clone(self):
        edge = self.__class__()
        edge._tasks.extend(self._tasks)

        return edge


class Entity(object):
    pass


class DynamicEntity(Entity):

    @abc.abstractmethod
    def transform(self):
        pass


class AppEdge(Edge):

    def append_app(self, app, args, subargs=None):

        self.append_task((app, args, subargs))

    def run(self, mgr, fwds=None):

        if not fwds:
            fwds = appdict()

        for appname, args, subargs in self._tasks:

            builtin_apps = mgr.get_builtin_apps()
            appcls, targs, sargs, objs = load_appclass(appname, args,
                    subargs, builtin_apps)
            app = appcls(mgr)

            # TODO: support ScriptGroup
            #if hasattr(app, "_dcasts"):
            #    app._dcasts.update(mgr.iter_downcast())

            fwds.update(objs)
            fwds["_bid_"] = self._bid

            ret, fwds = app.run(targs, sargs, fwds)

            if ret != 0:
                raise TestError("'%s' app returned %d." % (appname, ret))

        return 0, fwds


class DepartureEntity(Entity):

    def __new__(cls, *vargs, **kwargs):

        obj = super(DepartureEntity, cls).__new__(cls)
        obj._nodefwd = appdict()
        obj._num_edges = 0
        obj._edges = []

        return obj

    def update_edge_forward(self, edge, fwd):
        """update edge forwards"""

        if edge in self._nodefwd:
            self._nodefwd[edge].update(fwd)

        else:
            self._nodefwd[edge] = appdict(fwd)

    def get_edge_forward(self, edge):
        """get edge forwrads"""

        return appdict(self._nodefwd[edge])

    def add_edge(self, edge, arrivalnode):

        edge._bid = self._num_edges
        self._num_edges += 1

        self._edges.append((edge, arrivalnode))

    def get_departure_paths(self):
        return list(self._edges)


class ArrivalEntity(Entity):

    def __new__(cls, *vargs, **kwargs):

        obj = super(ArrivalEntity, cls).__new__(cls)

        if "unload_method" in kwargs:
            obj._unload_method = kwargs.pop("unload_method")

        else:
            obj._unload_method = "overwrite"

        return obj

    def is_quorum(self, total_arrivals, present_arrivals):

        return set(total_arrivals) == set(a[0] for a in present_arrivals)

    def unload_edgedata(self, edge, edgedata, nodedata, unload_method=None):

        # data => key: exit node, value: forwards
        if unload_method is None:
            unload_method = self._unload_method

        if unload_method == "block":
            pass

        elif unload_method == "overwrite":
            nodedata.update(edgedata)

        elif unload_method == "accumulate":
            for k, v in edgedata.items():
                if k in nodedata:
                    nodedata[k].append(v)

                else:
                    nodedata[k] = [v]
        else:
            raise UsageError("Unknown unload method: %s" % self._unload_method)

    def collect_edge_forward(self, arrived):
        """populate forwards with forwarded from each edges"""

        fwds = appdict()

        for edge, ret, afwd in arrived:
            if ret == 0:
                for k, v in afwd.items():
                    if k in fwds:
                        fwds[k].append(v)
                    else:
                        fwds[k] = [v]
            else:
                raise UsageError("Path '%s' returns with non-zero "
                                "code at '%s' with forward of '%s'." % (
                                str(edge), str(self), str(afwd)))
        return fwds


class EntryEntity(DepartureEntity):
    pass


class ExitEntity(ArrivalEntity):
    pass


class Node(DepartureEntity, ArrivalEntity):
    pass


class Group(MicroappObject, GroupArgParser, EntryEntity, ExitEntity):

    def __new__(cls, mgr, *vargs, **kwargs):

        kwargs["add_help"] = False

        if "description" not in kwargs and cls.__doc__ is not None:
            kwargs["description"] = cls.__doc__

        obj = super(Group, cls).__new__(cls, *vargs, **kwargs)
        obj.manager = mgr

        #obj._departures = appdict()

        obj._entrynodes = appdict()
        obj._exitnodes = appdict()
        obj._arrivalnodes = appdict()

        obj._edgeproxy_ready = []
        obj._edgeproxy_active = []
        obj._procpool = None
        
#        # add common arguments
#        # syntax --map 'f(x)[ -> y]'
#        obj.add_argument("--map", metavar="expr", delay=True,
#                action="append", help="map data from paths")
#        # syntax --filter 'x > 0[ -> y]'
#        obj.add_argument("--filter", metavar="expr", delay=True,
#                action="append", help="filter data from paths")
#        # syntax --reduce 'sum(x, y)[ -> y]'
#                #type=ArgType(None, True, None), action="append",
        obj.add_argument("--reduce", metavar="expr", delay=True,
                action="append", syshelp="reduce data from paths")
        obj.add_argument("--clone", metavar="expr", type=int,
                syshelp="dupulicate app")
        obj.add_argument("--assigned-input", metavar="data",
                syshelp="assign input data to a specified edge")

        obj.add_argument("--assert-in", metavar="expr", action="append",
                type=bool, syshelp="assertion test on input")
        obj.add_argument("--assert-out", metavar="expr", action="append",
                type=bool, delay=True, syshelp="assertion test on output")

        obj.add_argument("--data-join", type=str, default="overwrite",
                syshelp="data combining method(block, overwrite:default, accumulate)")

        #obj.add_argument('--version', action='version', version=(cls._name_
        #        + " " + cls._version_))

        #if sys.version_info >= (3, 0):
        obj.add_argument("--multiproc", help="# of processes")

        obj.logger = Logger(mgr, [obj._name_])

        return obj

    def group_forward(self, args, data):
        """populate forwards with forwarded from exit nodes"""

        fwdtype = args.data_join["_"]
        grp_fwds = appdict()

        # data => key: exit node, value: forwards

        if fwdtype == "block":
            pass

        elif fwdtype == "overwrite":
            # discard terminal node info
            for v in data.values():
                grp_fwds.update(v)

        elif fwdtype == "accumulate":
            # preserve arrival node info
            for k, v in data.items():
                if k in grp_fwds:
                    grp_fwds.append(v)

                else:
                    grp_fwds[k] = [v]
        else:
            raise UsageError("Unknown group forward type: %s" % fwdtype)

        return grp_fwds

    def connect_edge(self, departurenode, edge, arrivalnode):

        departurenode.add_edge(edge, arrivalnode)

        if arrivalnode in self._arrivalnodes:
            self._arrivalnodes[arrivalnode].append(edge)

        else:
            self._arrivalnodes[arrivalnode] = [edge]

        if isinstance(departurenode, EntryEntity):
            self._entrynodes[departurenode] = None

        if isinstance(arrivalnode, ExitEntity):
            self._exitnodes[arrivalnode] = None


    def is_finished(self, terms):

        if any((t not in terms) for t in self._exitnodes.keys()):
            return False

        return True

    @abc.abstractmethod
    def connect(self, args, subargs):
        pass


    def _init_multiproc(self, nprocs):


        nprocs = cpu_count() if nprocs=="*" else int(nprocs)
        self.logger.debug("nprocs = %d" % nprocs)

        for idx in range(nprocs):
            self._edgeproxy_ready.append(EdgeProxy(self.manager))

        # Pool([processes[, initializer[, initargs[, maxtasksperchild[, context]]]]])
        self._procpool = Pool(nprocs)

    def _multiproc(self, paths):

        async_res = None

        # check if available pipe
        if self._edgeproxy_ready and paths:
            proxy = self._edgeproxy_ready.pop()
            self._edgeproxy_active.append(proxy)

            dep, (edge, arr) = paths.pop()

            child_perform, args, kwargs = proxy.pack_launch(
                                            dep.get_edge_forward(edge), edge, arr)

            # async launch process
            self.logger.debug("edge %d runs on a new process: (%s, %s, %s)." %
                    (edge._eid, str(child_perform), str(args), str(kwargs)))

            async_res = self._procpool.apply_async(child_perform, args, kwargs)

            self.logger.debug("edge %d is async-launched." % edge._eid)

        return async_res

    def _fini_multiproc(self):

        if self._procpool:
            self._procpool.close()
            self._procpool.join()

    def run(self, args, sargs, fwds):
    
        self.logger.debug("Group '%s' runs: %s" % (self._name_, str(args)))

        # prepare env
        self._env.update(fwds)
        self._env.update(self.manager.iter_shared())
        self._env.update(self.manager.iter_downcast())

        # set program name and parse argument
        pos = sys.argv[0].rfind("-")
        if pos>=0:
            sys.argv[0] = sys.argv[0][:pos]+"-"+self._name_

        else:
            sys.argv[0] += "-"+self._name_

        self._arg_parser.prog = sys.argv[0]

        args, rargs = self.parse_known_args(args, self._env)

        if rargs and rargs[0].startswith("-"):
            raise UsageError("Unknown argument: %s" % " ".join(rargs))

        if args.assert_in:
            for ain in args.assert_in:
                if isinstance(ain, ArgType):
                    ain.env.update(self._env)
                    if not assert_check(ain):
                        raise TestError("Tested '%s' with '%s'" %
                            (str(ain.data), str(ain.env)))
                else:
                    raise InternalError("assert check argument is not ArgType")

        # build graph
        if rargs and sargs:
            sargs = rargs + ["--"] + sargs

        elif rargs:
            sargs = rargs

        elif not sargs:
            sargs = []

        ret = self.connect(args, sargs)

        # collect entry nodes
        #entrynodes = []
        #arrivalnodes = appdict()

        # build ready paths
        ready_paths = []

        if args.assigned_input:
            if args.assigned_input["_context_"]:
                assigned_input = (args.assigned_input["_context_"],
                                  args.assigned_input["_"])
            else:
                assigned_input = ("_", args.assigned_input["_"])

        else:
            assigned_input = tuple()

        for enode in self._entrynodes:
            # TODO: add --data argument
            for path in enode.get_departure_paths():
                edge = path[0]
                enode.update_edge_forward(edge, fwds)

                if assigned_input:
                    if edge._bid < len(assigned_input[1]):
                        ainput = {assigned_input[0]: assigned_input[1][edge._bid]}
                        enode.update_edge_forward(edge, ainput)

                ready_paths.append((enode, path))

        arrived = appdict()
        just_arrived = []
        terminated = appdict()

        if args.multiproc:
            if PY3:
                self._init_multiproc(args.multiproc["_"])

            else:
                pass
                # TODO: show info that py2 does not support multiproc

        while ready_paths or self._edgeproxy_active or just_arrived:

            #self.logger.debug("edgeproxy_active: %s" % str(self._edgeproxy_active))
            #self.logger.debug("just_arrived: %s" % str(just_arrived))

            # launch path if exists
            if ready_paths:
                if PY3 and args.multiproc:
                    self._multiproc(ready_paths)

                else:
                    # pop ready path
                    departure, (edge, arrival) = ready_paths.pop(0)

                    out, fwd = edge.run(self.manager, departure.get_edge_forward(edge))

                    just_arrived.append(arrival)

                    if arrival in arrived:
                        arrived[arrival].append((edge, out, fwd))

                    else:
                        arrived[arrival] = [(edge, out, fwd)]

            # handle paralle
            if PY3 and args.multiproc:

                processed = False

                # handle message from children
                for proxy in self._edgeproxy_active:
                    finished = proxy.process_child_message()

                    if finished:

                        arrival, edge, out, fwd = finished
                        self.logger.debug("edge %d is finished." % edge._eid)

                        just_arrived.append(arrival)

                        if arrival in arrived:
                            arrived[arrival].append((edge, out, fwd))

                        else:
                            arrived[arrival] = [(edge, out, fwd)]

                        self._edgeproxy_active.remove(proxy)
                        self._edgeproxy_ready.append(proxy)

                    processed = True

                if not processed:
                    for proxy in self._edgeproxy_active:
                        proxy.process_idletask()

            if just_arrived:

                arrival = just_arrived.pop(0)

                # check if quorum
                if arrival.is_quorum(self._arrivalnodes[arrival], arrived[arrival]):

                    node_fwds = appdict()

                    # merge paths' output
                    for edge, ret, afwd in arrived[arrival]:
                        if ret == 0:
                            if args.data_join:
                                arrival.unload_edgedata(edge, afwd, node_fwds,
                                        unload_method=args.data_join["_"])
                            else:
                                arrival.unload_edgedata(edge, afwd, node_fwds)

                        else:
                            raise UsageError("Path '%s' returns with non-zero "
                                            "code at '%s' with forward of '%s'." % (
                                            str(edge), str(self), str(afwd)))

                    # TODO: find good usage of dynamic node
                    if isinstance(arrival, DynamicEntity):
                        arrival = arrival.transform()

                    if isinstance(arrival, ExitEntity):
                        terminated[arrival] = node_fwds

                        # check if to exit
                        if self.is_finished(terminated):
                            break

                    elif isinstance(arrival, DepartureEntity):

                        # add paths
                        for path in arrival.get_departure_paths():
                            arrival.update_edge_forward(path[0], node_fwds)
                            ready_paths.append((arrival, path))

        if PY3 and args.multiproc:
            self._fini_multiproc()

        # merge terminals' output
        grp_fwds = self.group_forward(args, terminated)

        if not grp_fwds:
            grp_fwds = appdict()

        if args.reduce:
            for rdc in args.reduce:
                grp_fwds.update(reduce_arg(rdc, terminated))

#        if args.filter:
#            for flt in args.filter:
#                grp_fwds.update(filter_arg(flt, exit_fwds))
#
#        if args.map:
#            for mp in args.map:
#                grp_fwds.update(map_arg(mp, exit_fwds))

        if args.assert_out:
            for aout in args.assert_out:
                if isinstance(aout, ArgType):
                    aout.env.update(grp_fwds)
                    if not assert_check(aout):
                        raise TestError("Tested '%s' with '%s'" %
                            (str(aout.data), str(aout.env)))
                else:
                    raise InternalError("assert check argument is not ArgType")

        self.logger.debug("Group '%s' ran." % self._name_)

        return 0, grp_fwds


class GroupCmd(Group):

    _name_ = "group"
    _version_ = "0.1.2"

    def connect(self, args, subargs):

        edge = AppEdge()

        items = []
        for s in subargs:
            if s == "--":
                if items and not items[0].startswith("-"):
                    edge.append_app(items[0], items[1:])
                items = []
            else:
                items.append(s)

        if items and not items[0].startswith("-"):
            edge.append_app(items[0], items[1:])

        if args.clone:
            for idx in range(int(args.clone["_"])):
                self.connect_edge(self, edge.clone(), self)
        else:
            self.connect_edge(self, edge, self)

