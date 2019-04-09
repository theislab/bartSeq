from configparser import RawConfigParser

from ..configuration import Configuration, ConfigurationHandler


class PsConfigurationHandler(ConfigurationHandler):
    @staticmethod
    def write_standard_config(path):  # Write the standard settings file
        """ Write the standard settings file
        """
        config = RawConfigParser()
        config.set("DEFAULT", "threads", "1")

        config.add_section("Primer3")
        config.set("Primer3", "path", "/usr/local/bin/primer3_core")
        config.set("Primer3", "configPath", "")
        config.set("Primer3", "thermoParamPath", "")

        config.add_section("BLAST")
        config.set("BLAST", "path", "")
        config.set("BLAST", "db_path", "")
        config.set("BLAST", "db_name", "")
        config.set("BLAST", "maxHits", 5)

        config.add_section("RNAcofold")
        config.set("RNAcofold", "path", "")

        config.add_section("Optimization")
        config.set("Optimization", "steps", 500)
        config.set("Optimization", "maxTemp", 15)

        with open(path, "w") as configfile:
            config.write(configfile)

    @staticmethod
    def read_config(handle):  # Read the settings file
        """ Read the settings file
        """
        config = RawConfigParser()
        config.read_file(handle)
        c = Configuration()
        c.max_threads = config.getint("DEFAULT", "threads")
        c.p3_path = config.get("Primer3", "path")
        c.p3_config_path = config.get("Primer3", "configPath")
        c.p3_thermo_path = config.get("Primer3", "thermoParamPath")

        c.blast_path = config.get("BLAST", "path")
        c.blast_dbpath = config.get("BLAST", "db_path")
        c.blast_dbname = config.get("BLAST", "db_name")
        c.blast_max_hits = config.getint("BLAST", "maxHits")

        c.rnacf_path = config.get("RNAcofold", "path")

        c.opt_steps = config.getint("Optimization", "steps")
        c.opt_max_temp = config.getint("Optimization", "maxTemp")

        return c
