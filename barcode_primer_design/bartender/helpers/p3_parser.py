import configparser
import re
from logging import getLogger

from .primer import Primer, ExcludedRegion
from .primerpair import PrimerPair


log = getLogger(__name__)


def parse_p3_information(primer_pair_set, p3output):
    ini_str = f"[root]\n{p3output.rstrip('=')}"
    content = configparser.RawConfigParser()
    content.read_string(ini_str)

    log.info("P3 Config: %s", content)

    num_primers = content.getint("root", "PRIMER_PAIR_NUM_RETURNED")
    for i in range(num_primers):
        loc = list(map(int, content.get("root", "PRIMER_LEFT_" + str(i)).split(",")))
        # python style indexing, so we subtract 1 from the start to get the index left from the start
        fwd = Primer(
            content.get("root", "PRIMER_LEFT_" + str(i) + "_SEQUENCE"),
            loc[0] - 1,
            loc[1],
        )
        fwd.gc_content = content.getfloat(
            "root", "PRIMER_LEFT_" + str(i) + "_GC_PERCENT"
        )
        fwd.tm = content.getfloat("root", "PRIMER_LEFT_" + str(i) + "_TM")
        fwd.any = content.getfloat("root", "PRIMER_LEFT_" + str(i) + "_SELF_ANY_TH")
        fwd.self = content.getfloat("root", "PRIMER_LEFT_" + str(i) + "_SELF_END_TH")

        loc = list(map(int, content.get("root", "PRIMER_RIGHT_" + str(i)).split(",")))
        # python style indexing: we start at the index right of the end,
        # and go left, so index and position are equal
        rev = Primer(
            content.get("root", "PRIMER_RIGHT_" + str(i) + "_SEQUENCE"),
            loc[0],
            loc[1],
            True,
        )
        rev.gc_content = content.getfloat(
            "root", "PRIMER_RIGHT_" + str(i) + "_GC_PERCENT"
        )
        rev.tm = content.getfloat("root", "PRIMER_RIGHT_" + str(i) + "_TM")
        rev.any = content.getfloat("root", "PRIMER_RIGHT_" + str(i) + "_SELF_ANY_TH")
        rev.self = content.getfloat("root", "PRIMER_RIGHT_" + str(i) + "_SELF_END_TH")

        pair = PrimerPair(fwd, rev, primer_pair_set.name + "_" + str(i))
        pair.product_size = content.getint(
            "root", "PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"
        )
        primer_pair_set.append(pair)


def parse_p3seq_information(
    fwd_primer_set, rev_primer_set, p3output, sequence, p3options
):
    ini_str = "[root]\n" + p3output[: len(p3output) - 1]
    content = configparser.RawConfigParser()
    content.read_string(ini_str)

    for i in range(content.getint("root", "PRIMER_LEFT_NUM_RETURNED")):
        loc = list(map(int, content.get("root", "PRIMER_LEFT_" + str(i)).split(",")))
        # see comments in parse_p3_information
        fwd = Primer(
            content.get("root", "PRIMER_LEFT_" + str(i) + "_SEQUENCE"),
            loc[0] - 1,
            loc[1],
        )
        fwd.gc_content = content.getfloat(
            "root", "PRIMER_LEFT_" + str(i) + "_GC_PERCENT"
        )
        fwd.tm = content.getfloat("root", "PRIMER_LEFT_" + str(i) + "_TM")
        fwd.any = content.getfloat("root", "PRIMER_LEFT_" + str(i) + "_SELF_ANY_TH")
        fwd.self = content.getfloat("root", "PRIMER_LEFT_" + str(i) + "_SELF_END_TH")
        overlap = False
        for f in sequence.features:
            if (
                type(f) is ExcludedRegion
                and len(
                    set(range(fwd.location.start, fwd.location.end))
                    & set(range(f.location.start, f.location.end))
                )
                > 0
            ):
                overlap = True
                log.info(
                    "Left primer at position %s overlaps with excluded region",
                    fwd.location.start,
                )
                sequence.warning += (
                    "Left primer at position "
                    + fwd.location.start
                    + " overlaps with excluded region\n"
                )
        if not overlap:
            fwd_primer_set.append(fwd)

    for i in range(content.getint("root", "PRIMER_RIGHT_NUM_RETURNED")):
        loc = list(map(int, content.get("root", "PRIMER_RIGHT_" + str(i)).split(",")))
        # see comments in parse_p3_information
        rev = Primer(
            content.get("root", "PRIMER_RIGHT_" + str(i) + "_SEQUENCE"),
            loc[0],
            loc[1],
            True,
        )
        rev.gc_content = content.getfloat(
            "root", "PRIMER_RIGHT_" + str(i) + "_GC_PERCENT"
        )
        rev.tm = content.getfloat("root", "PRIMER_RIGHT_" + str(i) + "_TM")
        rev.any = content.getfloat("root", "PRIMER_RIGHT_" + str(i) + "_SELF_ANY_TH")
        rev.self = content.getfloat("root", "PRIMER_RIGHT_" + str(i) + "_SELF_END_TH")
        # rev.penalty = content.getfloat('root', 'PRIMER_RIGHT_' + str(i) + '_PENALTY')
        overlap = False
        for f in sequence.features:
            if (
                type(f) is ExcludedRegion
                and len(
                    set(range(rev.location.start, rev.location.end))
                    & set(range(f.location.start, f.location.end))
                )
                > 0
            ):
                overlap = True
                log.info(
                    "Right primer at position %s overlaps with excluded region",
                    rev.location.start,
                )
                sequence.warning += (
                    "Right primer at position "
                    + str(rev.location.start)
                    + " overlaps with excluded region\n"
                )
        if not overlap:
            rev_primer_set.append(rev)

    for i in range(content.getint("root", "PRIMER_RIGHT_NUM_RETURNED")):
        loc_left = list(
            map(int, content.get("root", "PRIMER_LEFT_" + str(i)).split(","))
        )
        loc_right = list(
            map(int, content.get("root", "PRIMER_RIGHT_" + str(i)).split(","))
        )
        product_size = get_max_product_size(p3options)
        if product_size is not None and loc_right[0] - loc_left[0] > product_size:
            sequence.warning += (
                "Amplicon length of primer pair "
                + str(i)
                + " ("
                + str(loc_right[0] - loc_left[0])
                + " bp) is larger than the maximum product size ("
                + str(product_size)
                + " bp)\n"
            )


def get_max_product_size(config_path):
    with open(config_path) as fp:
        product_size = None
        for line in fp:
            if not line:
                break
            if line.split("=")[0] == "PRIMER_PRODUCT_SIZE_RANGE":
                products = list(map(int, re.split(r"\W+", line.split("=")[1].strip())))
                product_size = max(products)
                break
    return product_size
