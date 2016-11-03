import os
import re
import logging
log = logging.getLogger(__name__)

class MooseInformationBase(object):
    """
    Base class for system and object information generators.

    Args:
        node[dict]: The YAML node associated with the object or system.

    Kwargs:
        All key, value pairs are stored in a configuration dict(), see MooseApplicationDocGenerator.
    """

    def __init__(self, node, **kwargs):

        self._yaml = node
        self._config = kwargs
        self._details = os.path.abspath(os.path.join(self._config['details'], self.filename(node['name']).replace('.moose.md', '.md')))
        log.debug(node['name'])

    def __str__(self):
        """
        Shows the complete markdown when print is called on this object.
        """
        return self.markdown()

    def markdown(self):
        """
        A virtual to return the generated markdown for the system/object.
        """
        pass

    @staticmethod
    def filename(name):
        """
        A virtual to return the filename generated for the system/object.
        """
        return None

    def write(self):
        """
        Write the markdown file.
        """

        # The complete name of the file to be created.
        filename = os.path.abspath(os.path.join(self._config['install'], self.filename(self._yaml['name'])))

        # Create the directories, if needed
        dirname = os.path.dirname(filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        # Generate the markdown
        md = self.markdown()

        """
        To enable the use of easier markdown link creation (see extensions.MooseMarkdownLinkPreprocessor) the
        relative path of the markdown file must be embedded in Markdown file.

        TODO: When (if?) the mkdocs plugin system is realized a plugin could be create to embedded this just prior to
        conversion to html rather than hard-code it in here.
        """
        md = '<!-- {} -->\n\n{}'.format(os.path.relpath(filename, os.getcwd()), md)

        # Do not re-write file if it exists (saves mkdocs from re-loading the world)
        if os.path.exists(filename):

            with open(filename, 'r') as fid:
                content = fid.read()

            if content == md:
                log.debug('No Change: {}'.format(filename))
                return

        # Write the file
        log.debug('Writing: {}'.format(filename))
        fid = open(filename, 'w')
        fid.write(md)
        fid.close()
