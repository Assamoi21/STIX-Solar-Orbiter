# -*- coding: utf-8 -*-
"""Microapp project module"""

from __future__ import print_function

import sys, os, shlex, multiprocessing, json, logging
from subprocess import Popen, PIPE

from microapp.base import MicroappObject, microapp_builtins
from microapp.error import (UsageError, InternalError, TypeCheckError,
                            NormalExit, ConfigError)
from microapp.utils import tostr, tostr_iter, appdict, Logger
from microapp.parse import ProjectArgParser, DataRegistry, data_transfer
from microapp.app import App
from microapp.group import Group, GroupCmd
from microapp.manage import Manager
from microapp.framework import load_appclass


# TODO create curcheck
_cfgcheck = {
}

class Project(MicroappObject, ProjectArgParser, Manager, DataRegistry):
    """Microapp Project Base-class"""

    def __new__(cls, *vargs, **kwargs):

        prog = kwargs.pop("prog", cls._name_)
        desc = kwargs.pop("description", cls._description_)
        epilog = kwargs.pop("epilog", None)

        obj = super(Project, cls).__new__(cls, *vargs, prog=prog,
                description=desc, epilog=epilog, **kwargs)

        if sys.version_info < (3, 0):
            obj.add_argument("--multiproc-method",
                    syshelp="multiprocessing spawning method")

        obj.add_argument("--forward", metavar="expr", action="append",
                delay=True, syshelp="forward variables to next app")

        obj.add_argument("--share", metavar="expr", action="append",
                delay=True, syshelp="share variables between sibling apps")

        obj.add_argument("--shell", metavar="cmd", action="append",
                delay=True, syshelp="execute a shell command")

        obj.add_argument("--downcast", metavar="expr", action="append",
                delay=True, syshelp="downcast variables under this app")

        obj.add_argument('--version', action='version', version=(prog + " "
                + cls._version_))

        obj.add_argument('--logging', metavar="level[,path]",
                         help="logging level and optionally log file path")

        obj._config = obj._load_config()
        obj.set_config("name", obj._name_, createall=True)
        obj.set_config("version", obj._version_)

        prjdcast = appdict({"name": obj._name_, "version": obj._version_})
        obj._dcasts["_project_"] = prjdcast

        obj._logger = None
        obj.logger = Logger(obj, [obj._name_])

        return obj

    def _load_config(self):

        home = os.path.expanduser("~")
        cfgdir = os.path.join(home, ".microapp", self._name_)

        if not os.path.isdir(cfgdir):
            os.makedirs(cfgdir)

        self._cfgfile = os.path.join(cfgdir, "config")

        if os.path.isfile(self._cfgfile):
            with open(self._cfgfile) as f:
                #config = json.load(f)
                config = json.load(f)

        else:
            import datetime
            now = datetime.datetime.now()

            created = appdict({"by": self._name_, "when": str(now)})
            global_conf = appdict({"created": created, "configdir": cfgdir,
                            })

            config = appdict({"global" : global_conf, "project": appdict(),
                              "app": appdict()})

            with open(self._cfgfile, "w") as f:
                json.dump(config, f, indent=4, sort_keys=True)

        return config


    def set_config(self, keylist, value, createall=False, curcfg=None,
                   curcheck=None):

        # TODO use curcheck for validation

        if isinstance(keylist, str):
            keylist = keylist.split(".")

        key = keylist[0]

        if curcfg is None:
            curcfg = self._config

            if key not in ("app", "global", "project"):
                key = "project"
                keylist.insert(0, key)
                keylist.insert(1, self._name_)

        if len(keylist)==1:
            # TODO: check if dynamic generation such as name or version
            curcfg[key] = value

        elif key in curcfg:
            self.set_config(keylist[1:], value, curcfg=curcfg[key],
                    createall=createall)

        elif createall:
            # TODO: check if creation is allowed
            newconfig = appdict()
            curcfg[key] = newconfig
            self.set_config(keylist[1:], value, curcfg=newconfig,
                    createall=createall)
        else:
            raise UsageError("Key '%s' does not exist." % key)

    def has_config(self, keylist, curcfg=None):

        if isinstance(keylist, str):
            keylist = keylist.split(".")

        key = keylist[0]

        if curcfg is None:
            curcfg = self._config

            if key not in ("app", "global", "project"):
                key = "project"
                keylist.insert(0, key)
                keylist.insert(1, self._name_)

        if key not in curcfg:
            return False

        elif len(keylist)==1:
            return True

        else:
            try:
                return self.get_config(keylist[1:], curcfg=curcfg[key])

            except ConfigError as err:
                return False

    def get_config(self, keylist, curcfg=None):

        if isinstance(keylist, str):
            keylist = keylist.split(".")

        key = keylist[0]

        if curcfg is None:
            curcfg = self._config

            if key not in ("app", "global", "project"):
                key = "project"
                keylist.insert(0, key)
                keylist.insert(1, self._name_)

        if key not in curcfg:
            raise ConfigError("Key '%s' is not found." % key)

        elif len(keylist)==1:
            return curcfg[key]

        else:
            return self.get_config(keylist[1:], curcfg=curcfg[key])

    def perform(self, aargs):
        pass

    def run_command(self, args=None, workdir=None, project_args=None, group_args=None,
                app_args=None, forward=None, shared=None, downcast=None):
        """MicroappProject command execution entry"""

        multiprocessing.freeze_support()

        try:

            if args is None:
                args = sys.argv[1:]

            elif isinstance(args, (str, bytes, bytearray)):
                if sys.platform == "win32":
                    args = shlex.split(tostr(args).replace("\\", "/"))

                else:
                    args = shlex.split(tostr(args))

            else:
                args = [a for a in tostr_iter(args)]

            if project_args:
                project_args = [a for a in tostr_iter(project_args)]

            else:
                project_args = []

            if not sys.stdin.isatty():
                project_args.append(tostr(sys.stdin.read().strip()))

            if group_args:
                group_args = [a for a in tostr_iter(group_args)]

            else:
                group_args = []

            if app_args:
                app_args = [a for a in tostr_iter(app_args)]

            else:
                app_args = []
            
        except Exception as err:
            print("ERROR: occured during argument collection for '%s': %s\n" %
                (str(args), str(err)), file=sys.stderr)
            sys.exit(1)

        self.logger.debug("Project '%s' runs command." % self._name_)

        try:
            ret = -1, None

            pwd = os.getcwd()

            if workdir and os.path.isdir(workdir):
                os.chdir(workdir)

            pargs, aargs = args, []

            for cidx, citem in enumerate(args):
                if citem == "--":
                    pargs, aargs = args[:cidx], args[cidx+1:]
                    break

            pargs += project_args

            pargs, rargs = self.parse_known_args(pargs, self._env)

            if pargs.shell:
                cmds = []

                for cmd in pargs.shell:
                    cmds.append(cmd.data.encode() + b"\n")

                shell = Popen(["sh"], stdin=PIPE, stdout=PIPE)
                shell.stdin.write(b"\n".join(cmds))
                shell.stdin.write(b"env\n")
                shell.stdin.close()

                for line in shell.stdout:
                    line = tostr(line.strip())

                    try:
                        name, value = line.split("=", 1)
                        os.environ[name] = value
                    except Exception as err:
                        pass

            if pargs.logging:
                items = [l.strip() for l in pargs.logging["_"].split(",", 1)]

                handlers = []

                self._logger = logging.getLogger()

                for item in items:
                    if item.upper() in ("CRITICAL", "ERROR", "WARNING",
                                        "WARN", "INFO", "DEBUG", "NOTSET"):
                        self._logger.setLevel(getattr(logging, item.upper()))

                    elif item in ("stdout", "stderr"):
                        handlers.append(logging.StreamHandler(
                                        stream=getattr(sys, item)))
                    else:
                        from logging import handlers as hdrs
                        handler = hdrs.RotatingFileHandler(item,
                                    maxBytes=1024*1024, backupCount=10)
                        handlers.append(handler)

                if not handlers:
                    handlers.append(logging.StreamHandler())

                # TODO: get format from config
                fmt = "%(asctime)s - %(microname)s - %(levelname)s - %(message)s"
                formatter = logging.Formatter(fmt)

                for handler in handlers:
                    handler.setFormatter(formatter)
                    self._logger.addHandler(handler)

            if sys.version_info < (3, 0) and hasattr(multiprocessing,
                    "set_start_method"):

                if pargs.multiproc_method:
                    multiprocessing.set_start_method(pargs.multiproc_method["_"])

                elif sys.platform == "darwin":
                    multiprocessing.set_start_method("spawn")

            self.perform(pargs)

            if pargs.forward:
                data_transfer(pargs.forward, self._fwds)

            if forward:
                self._fwds.update(forward)

            if pargs.share:
                data_transfer(pargs.share, self._shrds)

            if shared:
                self._shrds.update(shared)

            if pargs.downcast:
                data_transfer(pargs.downcast, self._dcasts)

            if downcast:
                self._dcasts.update(downcast)

            if rargs or aargs:

                app = GroupCmd(self)

                sys.argv[0] = self._name_

                ret = app.run(rargs+group_args, aargs, self._fwds)

            else:
                print(self.format_help())

                ret = 0, None

        except UsageError as err:
            print("USAGE ERROR: " + str(err))

        except ConfigError as err:
            print("CONFIG ERROR: " + str(err))

        except InternalError as err:
            print("INTERNAL ERROR: " + str(err))

        except TypeCheckError as err:
            print("TYPE MISMATCH: " + str(err))

        except NormalExit:
            ret = 0

        except (KeyboardInterrupt, EOFError):
            print('[Interrupted.]')

        finally:
            os.chdir(pwd)

            with open(self._cfgfile, "w") as f:
                json.dump(self._config, f, indent=4, sort_keys=True)

        self.logger.debug("Project '%s' ran." % self._name_)

        return ret

    def run_class(self, cls, workdir=None, project_args=None, group_args=None,
                app_args=None, forward=None, shared=None, downcast=None):
        """MicroappProject app execution entry"""

        multiprocessing.freeze_support()

        if project_args:
            project_args = [a for a in tostr_iter(project_args)]

        else:
            project_args = []

        if group_args:
            group_args = [a for a in tostr_iter(group_args)]

        else:
            group_args = []

        if app_args:
            app_args = [a for a in tostr_iter(app_args)]

        else:
            app_args = []

        self.logger.debug("Project '%s' runs a class." % self._name_)

        try:
            ret = -1, None

            pwd = os.getcwd()

            if workdir and os.path.isdir(workdir):
                os.chdir(workdir)

            pargs, _ = self.parse_known_args(project_args, self._env)


            if sys.version_info < (3, 0) and hasattr(multiprocessing,
                    "set_start_method"):
                if pargs.multiproc_method:
                    multiprocessing.set_start_method(
                                pargs.multiproc_method["_"])

                elif sys.platform == "darwin":
                    multiprocessing.set_start_method("spawn")

            self.perform(pargs)

            if pargs.forward:
                data_transfer(pargs.forward, self._fwds)

            if forward:
                self._fwds.update(forward)

            if pargs.share:
                data_transfer(pargs.share, self._shrds)

            if shared:
                self._shrds.update(shared)

            if pargs.downcast:
                data_transfer(pargs.downcast, self._dcasts)

            if downcast:
                self._dcasts.update(downcast)

            if cls:
                if issubclass(cls, Group):
                    grp = cls(self)
                    sys.argv[0] = self._name_
                    ret = grp.run(group_args, app_args, self._fwds)

                elif issubclass(cls, App):
                    app = cls(self)
                    sys.argv[0] = self._name_
                    ret = app.run(app_args, [], self._fwds)
            else:
                print(self.format_help())

                ret = 0, None

        except UsageError as err:
            print("USAGE ERROR: " + str(err))

        except ConfigError as err:
            print("CONFIG ERROR: " + str(err))

        except InternalError as err:
            print("INTERNAL ERROR: " + str(err))

        except TypeCheckError as err:
            print("TYPE MISMATCH: " + str(err))

        except NormalExit:
            ret = 0

        except (KeyboardInterrupt, EOFError):
            print('[Interrupted.]')

        finally:
            os.chdir(pwd)

            with open(self._cfgfile, "w") as f:
                json.dump(self._config, f, indent=4, sort_keys=True)

        self.logger.debug("Project '%s' ran." % self._name_)

        return ret

    def get_builtin_apps(self):

        dapps = []

        for scls in self.__class__.mro():
            if hasattr(scls, "_builtin_apps_"):
                dapps.extend(getattr(scls, "_builtin_apps_"))

        return dapps

    def log(self, level, msg, name):

        if self._logger:
            self._logger.log(level, msg, extra={"microname":".".join(name)})
