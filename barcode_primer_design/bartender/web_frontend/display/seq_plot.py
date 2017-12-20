import StringIO
import os
import random
from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from helpers.primer import ExcludedRegion, TargetRegion


class SeqPlotter:

    @staticmethod
    def seq_plot(sequence, path, title=""):
        gdd = GenomeDiagram.Diagram(title)
        gdt_features = gdd.new_track(1, greytrack=False, scale_smalltick_interval = 10, scale_largetick_interval = 100)
        gds_features = gdt_features.new_set()

        for f in sequence.features:
            if type(f) is TargetRegion:
                 gds_features.add_feature(f, name="Target", color=colors.blue, label=False)

        for f in sequence.features:
            if type(f) is ExcludedRegion:
                gds_features.add_feature(f, name="Excluded", color=colors.red, label=False)


        for i, fwd in enumerate(sequence.primer_set_fwd.set):
            gds_features.add_feature(fwd, name="Fwd " + str(i), color=colors.green, label=True, sigil="ARROW", label_size=11, arrowhead_length=0.25)

        for i, rev in enumerate(sequence.primer_set_rev.set):
            gds_features.add_feature(rev, name="Rev " + str(i), color=colors.green, label=True, sigil="ARROW", label_size=11, arrowhead_length=0.25)

        gdd.draw(format='linear', pagesize=(42*cm, 7*cm), fragments=1, start=0, end=len(sequence))

        filename = title + "_" + str(random.randrange(2**32)) + ".png"
        pic_location = os.path.join(os.path.realpath(path), filename)
        gdd.write(pic_location, "png")
        return filename