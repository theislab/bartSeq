from configuration import *


class AnalysisConfigurationHandler(ConfigurationHandler):
    @staticmethod
    def write_standard_config():  # Write the standard settings file
        """ Write the standard settings file

        """
        config = ConfigParser.RawConfigParser()
        config.add_section('BLAT')
        config.set('BLAT', 'path', '')
        config.set('BLAT', 'minScore', 0)
        config.set('BLAT', 'minIdentity', 95)
        with open("config.cfg", 'w') as configfile:
            config.write(configfile)

    @staticmethod
    def read_blat_config():  # Read the settings file
        """ Read the settings file

        """

        config = ConfigParser.RawConfigParser()
        config.read("config.cfg")
        c = Configuration()
        c.path = config.get('BLAT', 'path')
        c.min_score = config.getint('BLAT', 'minScore')
        c.min_identity = config.getint('BLAT', 'minIdentity')
        return c

