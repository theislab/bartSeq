from flask import render_template

from .seq_plot import seq_plot


def format_primer_set(arrangements, sequence_set):
    v = arrangements[0][1]
    w = arrangements[0][2]
    pairs = []
    for j, seq in enumerate(sequence_set):
        amplicon = w[j]
        pset = seq.amplicons[amplicon].primer_set
        pair = pset.set[v[j]]
        pairs.append((pair, pset.name, amplicon))
    return render_template(
        "primer-result.html", arrangement0=arrangements[0][0], pairs=pairs
    )


def format_seq_primer(sequence_set):
    tab_header = "<div role='tabpanel'>\n<ul class='nav nav-tabs' role='tablist'>\n"
    tab_body = "<div class='tab-content'>\n"
    for i, key in enumerate(list(sequence_set.keys())[::-1]):
        tab_header += "<li role='presentation'"
        if i == 0:
            tab_header += " class='active'"
        tab_header += (
            "><a href='#"
            + key
            + "' role='tab' data-toggle='tab'>"
            + key
            + "</a></li>\n"
        )

        tab_body += (
            "<div style='white-space: pre-wrap;' role='tabpanel' class='tab-pane"
        )
        if i == 0:
            tab_body += " active' "
        else:
            tab_body += "' "

        tab_body += "id='" + key + "'>"
        sequences = sequence_set[key].amplicons

        for seq in sequences:
            str_f = "Forward primers:<table style='width:100%'><tr>"
            str_r = "Reverse primers:<table style='width:100%'><tr>"
            for i, primer_fwd in enumerate(seq.primer_set_fwd.set):
                str_f += (
                    "<td>Fwd "
                    + str(i)
                    + ":</td><td>"
                    + primer_fwd.sequence
                    + "</td><td> Position: "
                    + str(primer_fwd.location.start)
                    + " - "
                    + str(primer_fwd.location.end)
                    + "</td>"
                    + "<td>Length: "
                    + str(len(primer_fwd))
                    + "</td><td>Tm: "
                    + str(primer_fwd.tm)
                    + "</td>"
                    + "<td>GC%: "
                    + str(primer_fwd.gc_content)
                    + "</td>"
                    + "<td>Any: "
                    + str(primer_fwd.any)
                    + "</td>"
                    + "<td>Self: "
                    + str(primer_fwd.self)
                    + "</td></tr>"
                )
            str_f += "</table><br>"

            for i, primer_rev in enumerate(seq.primer_set_rev.set):
                str_r += (
                    "<td>Rev "
                    + str(i)
                    + ":</td><td>"
                    + primer_rev.sequence
                    + "</td><td> Position: "
                    + str(primer_rev.location.start)
                    + " - "
                    + str(primer_rev.location.end)
                    + "</td>"
                    + "<td>Length: "
                    + str(len(primer_rev))
                    + "</td><td>Tm: "
                    + str(primer_rev.tm)
                    + "</td>"
                    + "<td>GC%: "
                    + str(primer_rev.gc_content)
                    + "</td>"
                    + "<td>Any: "
                    + str(primer_rev.any)
                    + "</td>"
                    + "<td>Self: "
                    + str(primer_rev.self)
                    + "</td>"
                    + "<td>Penalty: "
                    + str(primer_rev.self)
                    + "</td></tr>"
                )
            str_r += "</table><br>"
            pic_name = seq_plot(seq, "web_frontend/static/figures/p3seq", key)
            tab_body += (
                "<h4>Spacing: "
                + seq.spacing
                + "; Interval: "
                + seq.interval
                + "</h4><br>"
            )
            tab_body += '<br><img src="/static/figures/p3seq/' + pic_name + '"><br>'
            tab_body += str_f + str_r + "<br>"

            # rank mean penalty

            if seq.warning is not "":
                tab_body += (
                    "<div class='alert alert-warning' role='alert'>"
                    + seq.warning
                    + "</div><br>"
                )

            if seq.error is not "":
                tab_body += (
                    "<div class='alert alert-danger' role='alert'>"
                    + seq.error
                    + "</div><br>"
                )

        tab_body += "</div>"

    tab_header += "</ul>"
    tab_body += "</div>\n</div>"

    return tab_header + "\n" + tab_body
