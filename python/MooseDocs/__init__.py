import os
import sys
import argparse
import yaml
import argparse
import logging

import extensions
import database
import commands
import utils

# Check for the necessary packages, this does a load so they should all get loaded.
if utils.check_configuration(['yaml', 'mkdocs', 'markdown', 'markdown_include', 'mdx_math']):
    sys.exit(1)

import mkdocs
from mkdocs.commands import serve, build

from MarkdownTable import MarkdownTable
from MooseObjectParameterTable import MooseObjectParameterTable
from MooseObjectInformation import MooseObjectInformation
from MooseSystemInformation import MooseSystemInformation
from MooseApplicationSyntax import MooseApplicationSyntax
from MooseApplicationDocGenerator import MooseApplicationDocGenerator
from MooseSubApplicationDocGenerator import MooseSubApplicationDocGenerator

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

MOOSE_DIR = os.getenv('MOOSE_DIR', os.path.join(os.getcwd(), 'moose'))
if not os.path.exists(MOOSE_DIR):
    MOOSE_DIR = os.path.join(os.getenv('HOME'), 'projects', 'moose')

class MkMooseDocsFormatter(logging.Formatter):
    """
    A formater that is aware of the class hierachy of the MooseDocs library.

    Call the init_logging function to initialize the use of this custom fomatter.
    """
    COLOR = {'DEBUG':'CYAN', 'INFO':'RESET', 'WARNING':'YELLOW', 'ERROR':'RED', 'CRITICAL':'MAGENTA'}
    COUNTS = {'DEBUG':0, 'INFO':0, 'WARNING':0, 'ERROR':0, 'CRITICAL':0}

    def format(self, record):
        msg = logging.Formatter.format(self, record)

        if record.name.endswith('Item'):
            level = 4
        elif record.name.endswith('Database'):
            level = 3
        elif record.name.endswith('MooseInformationBase') or record.name.endswith('MooseObjectInformation') or record.name.endswith('MooseApplicationSyntax') or record.name.endswith('MooseSystemInformation') or record.name.endswith('MooseCommonFunctions'):
            level = 2
        elif record.name.endswith('MooseSubApplicationDocGenerator'):
            level = 1
        else:
            level = 0

        if record.levelname in ['DEBUG', 'WARNING', 'ERROR', 'CRITICAL']:
            msg = '{}{}: {}'.format(' '*4*level, record.levelname, msg)
        else:
            msg = '{}{}'.format(' '*4*level, msg)

        if record.levelname in self.COLOR:
            msg = utils.colorText(msg, self.COLOR[record.levelname])

        # Increment counts
        self.COUNTS[record.levelname] += 1

        return msg

def init_logging(verbose=False):
    """
    Call this function to initialize the MooseDocs logging formatter.
    """

    # Setup the logger object
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    # The markdown package dumps way too much information in debug mode (so always set it to INFO)
    log = logging.getLogger('MARKDOWN')
    log.setLevel(logging.INFO)

    # Setup the custom formatter
    log = logging.getLogger('MooseDocs')
    formatter = MkMooseDocsFormatter()
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(level)

    log = logging.getLogger('mkdocs')
    log.addHandler(handler)
    log.setLevel(level)

    return formatter

def yaml_load(filename, loader=yaml.Loader):
    """
    Load a YAML file capable of including other YAML files.

    Args:
        filename[str]: The name fo the file to load.
        loader[yaml.Loader]: The loader to utilize.
    """

    def include(self, node):
        """
        Allow for the embedding of yaml files.
        http://stackoverflow.com/questions/528281/how-can-i-include-an-yaml-file-inside-another
        """
        filename = os.path.join(self._root, self.construct_scalar(node))
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                return yaml.load(f, Loader)

    class Loader(loader):
        """
        """

        def __init__(self, stream):
            """
            Store the root directory for including other yaml files.
            """
            if isinstance(stream, file):
                self._root = os.path.split(stream.name)[0]
            else:
                self._root = os.getcwd()
            super(Loader, self).__init__(stream)

    ## Attach the include constructor to our custom loader.
    Loader.add_constructor('!include', include)

    with open(filename, 'r') as fid:
        yml = yaml.load(fid.read(), Loader)

    return yml

def load_pages(filename, keys=[], **kwargs):
    """
    A YAML loader for reading the pages file.

    Args:
        filename[str]: The name fo the file to load.
        keys[list]: A list of top-level keys to include.
        kwargs: key, value pairs passed to yaml_load function.
    """

    # Load the yaml data
    pages = yaml_load(filename, **kwargs)

    # Restrict the top-level keys to those provided in the 'include' argument
    if keys:
        pages = [page for page in pages if page.keys()[0] in keys]

    return pages

def purge(extensions):
    """
    Removes generated files from repository.

    Args:
        extensions[list]: List of file extensions to purge (.e.g., 'png'); it will be prefixed with '.moose.'
                          so the files actually removed are '.moose.png'.
    """
    for i, ext in enumerate(extensions):
        extensions[i] = '.moose.{}'.format(ext)

    log = logging.getLogger('MooseDocs')
    for root, dirs, files in os.walk(os.getcwd(), topdown=False):
        for name in files:
            if any([name.endswith(ext) for ext in extensions]):
                full_file = os.path.join(root, name)
                log.debug('Removing: {}'.format(full_file))
                os.remove(full_file)

def command_line_options():
    """
    Return the command line options for the moosedocs script.
    """

    # Command-line options
    parser = argparse.ArgumentParser(description="Tool for building and developing MOOSE and MOOSE-based application documentation.")
    parser.add_argument('--verbose', '-v', action='store_true', help="Execute with verbose (debug) output.")
    parser.add_argument('--config-file', type=str, default=os.path.join('moosedocs.yml'), help="The configuration file to use for building the documentation using MOOSE. (Default: %(default)s)")

    subparser = parser.add_subparsers(title='Commands', description="Documentation creation command to execute.", dest='command')

    # Generate options
    generate_parser = subparser.add_parser('generate', help="Generate the markdown documentation from MOOSE application executable. This is done by the serve and build command automatically.")
    generate_parser.add_argument('--purge', '-p', action='store_true', help="Remove all generated content (*.moose.md, *.moose.svg, *.moose.yml files) from the install directories.")

    # Serve options
    serve_parser = subparser.add_parser('serve', help='Generate and Sever the documentation using a local server.')
    serve_parser.add_argument('--livereload', dest='livereload', action='store_const', const='livereload', help="Enable the live reloading server.")
    serve_parser.add_argument('--dirtyreload', dest='livereload', action='store_const', const='dirtyreload', help="Enable the live reloading server without rebuilding entire site with single file change (default).")
    serve_parser.add_argument('--no-livereload', dest='livereload', action='store_const', const='no-livereload', help="Disable the live reloading of the served site.")
    serve_parser.add_argument('--strict', action='store_true', help='Enable strict mode and abort on warnings.')
    serve_parser.add_argument('--dirty', action='store_false', dest='clean', help='Do not clean the temporary build prior to building site.')

    # Build options
    build_parser = subparser.add_parser('build', help='Generate and Build the documentation for serving.')

    # Both build and serve need config file
    for p in [serve_parser, build_parser]:
        p.add_argument('--theme', help="Build documentation using specified theme. The available themes are: cosmo, cyborg, readthedocs, yeti, journal, bootstrap, readable, united, simplex, flatly, spacelab, amelia, cerulean, slate, mkdocs")
        p.add_argument('--pages', default='pages.yml', help="YAML file containing the pages that are supplied to the mkdocs 'pages' configuration item.")
        p.add_argument('--page-keys', default=[], nargs='+', help='A list of top-level keys from the "pages" file to include. This is a tool to help speed up the serving for development of documentation.')
    # Parse the arguments
    options = parser.parse_args()

    # Set livereload default
    if options.command == 'serve' and not options.livereload:
        options.livereload = 'dirtyreload'

    return options

def moosedocs():

    # Options
    options = command_line_options()

    # Initialize logging
    formatter = init_logging(options.verbose)
    log = logging.getLogger('MooseDocs')

    # Remove moose.svg files (these get generated via dot)
    log.info('Removing *.moose.svg files from {}'.format(os.getcwd()))
    purge(['svg'])

    # Execute command
    if options.command == 'generate':
        if options.purge:
            log.info('Purging *.moose.md, *.moose.yml, and *.moose.svg files from {}'.format(os.getcwd()))
            purge(['md', 'yml', 'svg'])
        commands.generate(config_file=options.config_file)
    elif options.command == 'serve':
        commands.serve(config_file=options.config_file, strict=options.strict, livereload=options.livereload, clean=options.clean, theme=options.theme, pages=options.pages, page_keys=options.page_keys)
    elif options.command == 'build':
        commands.build(config_file=options.config_file, theme=options.theme, pages=options.pages, page_keys=options.page_keys)

    # Display logging results
    print 'WARNINGS: {}  ERRORS: {}'.format(formatter.COUNTS['WARNING'], formatter.COUNTS['ERROR'])
    return formatter.COUNTS['ERROR'] > 0
