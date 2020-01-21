from pygccxml import utils, declarations, parser
from pprint import pprint
import sys
import inspect

def dump_class(decl):
    fullname = declarations.declaration_utils.full_name(decl)
    
    memberdefs = []
    for member in decl.public_members:
        d = def_decl(member)
        if d:
            memberdefs.append(f".{def_decl(member)}")
    
    classdef = f'py::class_<{fullname}>(m, "{decl.name}")'
    memberdefs = '\n'.join(memberdefs)
    return f'{classdef}\n{memberdefs};'

def args_def(decl):
    argstrings = []
    for arg in decl.arguments:
        adef = f'"{arg.name}"_a'
        if arg.default_value is not None:
            adef += "=" + arg.default_value
             
        argstrings.append(adef)
    
    return argstrings


def func_def(decl, funcname=None):
    typestring = decl.create_decl_string()

    argstrings = args_def(decl)

    fullname = declarations.declaration_utils.full_name(decl)
    
    if funcname is None:
        funcname = decl.name

    ptr = f"({typestring}) &{fullname}"
    bargs = ", ".join([f'"{funcname}"', ptr, *argstrings])
    
    bdef = f"def({bargs})"
    
    return bdef

# TODO: Handle more operators
op_map = {
        '()': '__call__',
    }

def getop(decl):
    if decl.symbol == '[]':
        if len(decl.arguments) == 1:
            return '__getitem__'
        if len(decl.arguments) == 2:
            return '__setitem__'

    return op_map.get(decl.symbol)

def op_def(decl):
    funcname = getop(decl)
    if not funcname: return ''
    
    return func_def(decl, funcname=funcname)

def var_def(decl):
    fullname = declarations.declaration_utils.full_name(decl)
    if isinstance(decl.decl_type, declarations.const_t):
        func = "def_readonly"
    else:
        func = "def_readwrite"
    return f'{func}("{decl.name}", &{fullname})'

def constructor_def(decl):
    argtypes = ", ".join(map(str, decl.argument_types))
    #args = ", ".join(args_def(decl))
    func = f'py::init<{argtypes}>()'
    bargs = ", ".join([func, *args_def(decl)])
    return f'def({bargs})'

def dump_function(decl):
    return f"m.{func_def(decl)};"

def dump_var(decl):
    fullname = declarations.declaration_utils.full_name(decl)
    return f'm.{var_def(decl)};'

dumpers = {
        declarations.free_function_t: dump_function,
        declarations.class_t: dump_class,
        #declarations.variable_t: dump_var
}

deffers = {
        declarations.free_function_t: func_def,
        declarations.variable_t: var_def,
        declarations.constructor_t: constructor_def,
        declarations.member_function_t: func_def,
        declarations.member_operator_t: op_def
}

def dump_decl(decl):
    t = type(decl)
    f = dumpers.get(t, lambda *args, **kwargs: '')
    return f(decl)

def def_decl(decl):
    t = type(decl)
    f = deffers.get(t, lambda *args, **kwargs: '')
    return f(decl)


def main(namespace, *headers):
    generator_path, generator_name = utils.find_xml_generator()
    xml_generator_config = parser.xml_generator_configuration_t(
        xml_generator_path=generator_path,
        xml_generator=generator_name)

    decls = parser.parse(headers, xml_generator_config)
    global_namespace = declarations.get_global_namespace(decls[0])

    ns = global_namespace.namespace(namespace)
    
    def out(s):
        if not s: return
        print(s)

    out("#include <pybind11/pybind11.h>")
    out("#include <pybind11/stl.h>")
    for header in headers:
        out(f'#include "{header}"')
    out("namespace py = pybind11;")
    out(f'PYBIND11_MODULE({namespace}, m){{')
    out("using namespace pybind11::literals;")
    for decl in ns.declarations:
        out(dump_decl(decl))
    out('}')


if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
