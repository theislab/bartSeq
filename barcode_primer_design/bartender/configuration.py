class Configuration:
    pass


class ConfigurationHandler:
    """Interface for a configuration handler"""

    @staticmethod
    def write_standard_config(self):  # Write the standard settings file
        """ Write the standard settings file
        """
        pass

    @staticmethod
    def read_config(self):
        """ Read the settings file and return all settings
        """
        pass
