#!/usr/bin/env python
"""
Generate dependency makefiles 

This is a rewrite of the bash mkdep script in python
"""
from __future__ import print_function
import re
import fnmatch, os

class OrderedSet(object):
  """
  A set that retains the order in which items were added

  This is a simple wrapper around OrderedDict to give it a set-like API
  """
  def __init__(self, items=[]):
    from collections import OrderedDict
    self.dict = OrderedDict((item,None) for item in items)

  def append(self, item):
    self.dict[item] = None

  def append_list(self, items):
    for item in items:
      self.append(item)

  def remove(self, item):
    del(self.dict[item])

  def __getitem__(self, key):
    return list(self.dict.keys())[key]

  def __len__(self):
    return len(self.dict.keys())

  def __add__(self, other):
    sum = OrderedSet(self)
    sum.append_list(other)
    return sum

#
# Line parsing code
#

TYPE_SUBROUTINE = 1
TYPE_FUNCTION = 2
TYPE_MODULE = 3

definition_re = re.compile(
                  "^\s*" # white space at start
                  "("    # open group for definition type
                  "subroutine|interface|module"
                  "|(?:[\w\d\*:,]+\s+)*(function)" # optional type identifier before function
                  ")"    # close group
                  "\s+"  # white space
                  "(\w+)", # definition name
                  re.I
                )

usage_re = re.compile("^[^!'\"]*"  # anything that isn't a comment or string
                      "(call|use)" # containing call or use
                      "\s+"        # followed by white space
                      "(\w+)",     # and a definition name
                      re.I
                      )

external_re = re.compile("(?:^|[^'\"]*[,\s])external(?:\s+|::)", re.I)
header_re = re.compile("\s*include ['\"](.*)['\"]", re.I)
comment_re = re.compile("^\s*(?:!|c\s+)")


def find_definition(line):
  match = definition_re.match(line)
  if match is None: return (None,None)

  groups = match.groups()
  type = groups[1] or groups[0]
  type = type.lower()
  name = groups[2].lower()

  type_map = {
    'subroutine': TYPE_SUBROUTINE,
    'interface': TYPE_SUBROUTINE,
    'function': TYPE_FUNCTION,
    'module': TYPE_MODULE,
  }
  assert(type in type_map)

  return (type_map[type], name.lower())

def find_usage(line):
  match = usage_re.match(line)
  if match is None: return (None,None)

  type, name = match.groups()
  type = type.lower()
  name = name.lower()

  type_map = {
    'call': TYPE_SUBROUTINE,
    'use': TYPE_MODULE
  }
  assert(type in type_map)

  return (type_map[type], name)

def find_external(line):
  line = line.lower()
  match = external_re.match(line)
  if match is None: return []

  components = line.split('::')
  if len(components) == 1: components = line.split('external',1)
  return [s.strip() for s in components[1].split(',')]

def find_header(line):
  match = header_re.match(line)
  if match is None: return None 

  return match.groups()[0]

#
# File finding helpers
#

def find_files_matching(pattern, dir='./'):
  for root, dirs, files in os.walk(dir):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

def find_fortran_files(dir='./'):
  return find_files_matching("*.f90", dir)



def process_file(filename):
  with open(filename) as f:

    definitions = []
    used_routines = []
    used_modules = []
    headers = []

    curdir = os.path.split(filename)[0]

    in_module = False

    end_module_re = re.compile("^.*end module", re.I)
    for line in f:
      # skip comments
      if comment_re.match(line): continue

      if in_module:
        if end_module_re.match(line):
          in_module = False
          continue

      else:
        # routines in a module can conflict with routines outside of modules
        # since the module dependency is sufficient, we simply skip routine
        # definitions inside of modules
        type, name = find_definition(line)
        if type is not None:
          if type is TYPE_MODULE: in_module = True
          definitions.append(name)
          continue

      type, name = find_usage(line)
      if type is not None:
        if type == TYPE_MODULE:
          used_modules.append(name)
        else:
          used_routines.append(name)
        continue

      externals = find_external(line)
      if externals:
        for item in externals:
          used_routines.append(item)
        continue

      header = find_header(line)
      if header:
        header = "./" + os.path.normpath(os.path.join(curdir, header))
        if (os.path.exists(header)):
          headers.append(header)

    headers = list(set(headers))
    return definitions, used_routines, used_modules, headers

def process_files(files):
  definition_map = {}
  usage_map = {}

  for file in files:
    definitions, used_routines, used_modules, headers = process_file(file)

    for item in definitions:
      definition_map[item] = file

    usage_map[file] = (used_routines, used_modules, headers)

  return definition_map, usage_map

def find_source_files(dir='./'):
  not_driver =  lambda f: "Drivers/" not in f
  return filter(not_driver, find_fortran_files(dir))

#
# The primary function for recursively determining dependencies
#

def find_deps(filename, definition_map, usage_map, routine_file_set=None, module_file_set=None):
  """
  Generate recursive list of files that file depends on.

  Parameters:
    filename: the filename of the main program source
    definition_map: map of symbol definitions to files (see process_files())
    usage_map: map of files to used symbol (see process_files())

  Returns:
    (routine_dep_files, module_dep_files)
    OrderedSets of routines and module dependencies
  """

  is_root = False
  used_routines, used_modules, _ = usage_map[filename]

  routine_files = [definition_map[routine] for routine in used_routines
                   if routine in definition_map # some calls are pulled in from libraries (e.g. MPI)
                  ]
  module_files = [definition_map[mod] for mod in used_modules]

  if routine_file_set is None:
    is_root = True
    routine_file_set = OrderedSet([filename])

  if module_file_set is None:
    module_file_set = OrderedSet()

  # recurse over routines we haven't seen before
  for rf in routine_files:
    #print("%s -) %s" % (filename, rf))
    if rf in routine_file_set: continue

    routine_file_set.append(rf)
    #print("%s -> %s" % (filename, rf))
    find_deps(rf, definition_map, usage_map, routine_file_set, module_file_set)

  # recurse over modules we haven't seen before
  for mf in module_files:
    if mf in module_file_set: continue

    #print("%s -> %s" % (filename, mf))

    # add  module to set to prevent recursive loops
    module_file_set.append(mf)

    # find dependencies recursively
    find_deps(mf, definition_map, usage_map, routine_file_set, module_file_set)

    # now move this module to the end of the list (after dependent modules)
    module_file_set.remove(mf)
    module_file_set.append(mf)

  if is_root:
    routine_files = list(filter(lambda f: f not in module_file_set, routine_file_set))
    module_files = list(module_file_set)

    return routine_files, module_files

#
# Output functions
#

def write_makefile(outfile, program_name, routine_files, module_files):
  outfile.write("%sSRC = \\\n" % program_name)
  write_file_list(outfile, routine_files)
  outfile.write("%s_MODULESRC = \\\n" % program_name)
  write_file_list(outfile, module_files)

def write_file_list(outfile, file_list):
  for i,f in enumerate(file_list):
    outfile.write(f)
    outfile.write(" ")
    if i % 3 == 2 and i != len(file_list) - 1:
      outfile.write("\\\n")
  outfile.write("\n")

def write_dependencies(outfile, definition_map, usage_map):
  to_obj = lambda f: f.replace('.f90', '.o')

  extra_deps = ['Compiler.mk']
  for f in sorted(usage_map.keys()):
    mods = [definition_map[m] for m in usage_map[f][1] if m in definition_map]
    mods = OrderedSet(mods)
    if f in mods: mods.remove(f)
    mods = list(mods)
    headers = usage_map[f][2]
    depstr = ' '.join(map(to_obj, mods + headers + extra_deps))
    outfile.write("%s: %s\n" % (to_obj(f), depstr))

#
# Input
#

def read_program_list(filename):
  with open(filename) as f:
    return  [line.strip() for line in f]

#
# Test code
#

def test_find_definition():
 
  tests = [
    (TYPE_SUBROUTINE, [
      ('subroutine foobar(a, b, c)', 'foobar'),
      ('   SUBROUTINE FOOBAR(a, b, c)', 'foobar'),
      ('interface foobar', 'foobar'),
    ]),
    (TYPE_FUNCTION, [
      ('function foobar(a, b, c)', 'foobar'),
      ('double precision function foobar(a, b, c)', 'foobar'),
    ]),
    (TYPE_MODULE, [
      ('module foobar', 'foobar'),
    ])
  ]

  for type, subtests in tests:
    for input, output in subtests:
      assert(find_definition(input) == (type, output))

def test_find_usage():
 
  tests = [
    (TYPE_SUBROUTINE, [
      ('call foobar(a, b, c)', 'foobar'),
      ('   CALL FOOBAR(a, b, c)', 'foobar'),
      ('  20 call foobar(a, b, c)', 'foobar'),
    ]),
    (TYPE_MODULE, [
      ('use foobar', 'foobar'),
    ])
  ]

  for type, subtests in tests:
    for input, output in subtests:
      #print input, find_usage(input), output
      assert(find_usage(input) == (type, output))

def test_find_external():
  tests = [
    ('external f', ['f']),
    ('double precision, external g', ['g']),
    ('external f, g, h', ['f', 'g', 'h']),
    ('integer external :: f, g', ['f', 'g'])
  ]

  for test, result in tests:
    assert(find_external(test) == result)

if __name__ == "__main__":
  import sys

  if len(sys.argv) > 1:
    program_files = sys.argv[1:]
  else:
    program_files = read_program_list("Utility/programs")

  program_files = ["./"+f if not f.startswith("./") else f for f in program_files]

  # list of main program files
  program_names = [os.path.splitext(os.path.split(pf)[1])[0] for pf in program_files]

  # process source code to build maps of provided and dependent routines / modules
  print("Mapping files to definitions")
  definition_map, usage_map = process_files(find_source_files())

  depdir = "DEP"

  for program_name, program_file in zip(program_names, program_files):
    print("Processing %s" % program_name)
    routines, modules = find_deps(program_file, definition_map, usage_map) 
    print("  Number of source files: %d" % len(routines))
    print("            module files: %d" % len(modules))
    with open("%s/%s.mk" % (depdir, program_name), "w") as f:
      write_makefile(f, program_name, routines, modules)

    #Kevin Jorissen 2014: make alternative file without FEFF blas/lapack for use with external MKL blas/lapack:
    try: routines.remove("./MATH/lu.f90")
    except: pass
    with open("%s/%s_MKL.mk" % (depdir, program_name), "w") as f:
      write_makefile(f, program_name, routines, modules)



  with open("%s/dependencies.mk" % depdir, "w") as f:
    write_dependencies(f, definition_map, usage_map)
