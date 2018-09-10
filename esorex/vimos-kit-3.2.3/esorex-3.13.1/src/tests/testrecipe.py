# Simple recipe plugin for unit testing.

import os

class CplPlugin(object):

    name = "testrecipe"
    version = 123
    synopsis = "test recipe"
    description = "Simple unit test recipe"
    author = "Test"
    email = "author@test.org"
    copyright = "copyright"
    parameters = [
            {
                'class': 'value',
                'name': 'par1',
                'description': 'test parameter',
                'context': 'test',
                'default': 3
            }
        ]

    def execute(self, plugin):
        par1 = plugin['parameters'][0]['value']
        frame = plugin['frames'][0]
        infilename = frame['filename']
        frame['type'] = 1 << 5
        frame['group'] = 1
        frame['level'] = 3
        with open(infilename) as infile:
            number = int(infile.read());
        outfilename = os.path.join(os.path.dirname(infilename),
                                   'python_recipe_test_output.fits')
        with open(outfilename, 'w') as outfile:
            outfile.write("{0}\n".format(number + par1));
        plugin['frames'].append({
                "filename": outfilename,
                "tag": "PROD",
                "type": 1 << 5,
                "group": 3,
                "level": 3
            })
        return 0
