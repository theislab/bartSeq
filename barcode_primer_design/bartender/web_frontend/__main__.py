import io
import os
import pickle
from logging import getLogger, basicConfig
from pathlib import Path
from typing import Iterable

from flask import Flask, Markup, render_template, make_response, url_for, abort
from flask_socketio import SocketIO
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from werkzeug.utils import secure_filename, redirect
from wtforms import (
    validators,
    StringField,
    TextAreaField,
    SelectField,
    FileField,
    Field,
)
from wtforms.fields.html5 import IntegerField

from ..helpers.primer import Amplicon, Gene
from ..p3seq import P3Seq
from .. import primer_select
from ..primer_select import PsConfiguration, Arrangement
from .display.seq_plot import seq_plot

IN_DOCKER = "docker" in Path("/proc/1/cgroup").read_text()

HERE = Path(__file__).parent
PATH_PREFIX = Path("/barcode_primer_design/bartender") if IN_DOCKER else HERE.parent
UPLOAD_FOLDER = os.path.join("web_frontend", "uploads")
ALLOWED_EXTENSIONS = {"txt", "pdf", "png", "jpg", "jpeg", "gif", "cfg"}


basicConfig(
    level="INFO", format="[%(asctime)s] %(levelname)s in %(module)s: %(message)s"
)
log = getLogger(__name__)

app = Flask(__name__)
app.secret_key = "A0Zr98j/3yX R~XHH!jmN]LWX/,?RT"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
socketio = SocketIO(app)

with (PATH_PREFIX / "config.cfg").open() as f:
    predefined_config = f.read()


class HorizontalForm(FlaskForm):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        field: Field
        for field in self:
            if not field.render_kw:
                field.render_kw = {}
            is_file = isinstance(field, FileField)
            field.render_kw[
                "class_"
            ] = f"col-sm-10 form-control{'-file' if is_file else ''}"
            field.label = field.label(class_="col-sm-2 col-form-label")


class PrimerSelectForm(HorizontalForm):
    input = TextAreaField("Input Sequences", validators=[validators.DataRequired()])
    predefined = TextAreaField("Predefined primers")
    configuration = FileField(
        "Primer3 configuration file",
        validators=[FileAllowed(["txt"], "Text files only!")],
    )
    blast_hits = IntegerField(
        "Max. BLAST hits", validators=[validators.DataRequired()], default=5
    )
    blast_db = SelectField(
        "BLAST database",
        choices=[
            ("hg38.fa", "UCSC Genome hg38"),
            ("hg19.fa", "UCSC Genome hg19"),
            ("refMrna.fa", "RefSeq mRNA"),
            ("mrna.fa", "GenBank mRNA"),
        ],
    )
    left_linker = StringField("Left Linker", default="ATGCGCATTC")
    right_linker = StringField("Right Linker", default="AGCGTAACCT")


class P3seqForm(HorizontalForm):
    input = TextAreaField("Input Sequences", validators=[validators.DataRequired()])
    spacing = IntegerField(
        "Range spacing", validators=[validators.DataRequired()], default=500
    )
    interval = IntegerField(
        "Range interval", validators=[validators.DataRequired()], default=250
    )
    configuration = FileField(
        "Primer3 configuration file",
        validators=[FileAllowed(["txt"], "Text files only!")],
    )


@app.route("/")
def main_page():
    return redirect(url_for("primerselect"))


@app.route("/about")
def about_page():
    return render_template("about.html")


@app.route("/primerselect", methods=("GET", "POST"))
def primerselect():
    form = PrimerSelectForm()
    if not form.validate_on_submit():
        return render_template("primerselect.html", form=form)

    input_string = io.StringIO(form.input.data)
    if form.predefined.data != "":
        predefined = io.StringIO(form.predefined.data)
    else:
        predefined = None
    try:
        try:
            filename = secure_filename(
                os.path.join("uploads", form.configuration.data.filename)
            )
            form.configuration.data.save(filename)

            # fix up config files
            with open(filename) as infile:
                lines = infile.readlines()
                with open(filename, "w", newline="\n") as outfile:
                    for line in lines:
                        line = line.replace(
                            "http://primer3.sourceforge.net", "http://primer3.org"
                        )
                        outfile.write(line)
        except AttributeError as e:
            filename = str(PATH_PREFIX / "web_frontend" / "primer3_settings.txt")

        config = PsConfiguration.read_config(PATH_PREFIX / "config.cfg")
        config.p3_config_path = filename
        config.p3_thermo_path = str(PATH_PREFIX.parent / "primer3_config")
        config.blast_dbpath = str(PATH_PREFIX.parent / "databases")
        config.blast_dbname = form.blast_db.data
        config.blast_max_hits = form.blast_hits.data
        if form.input.data.count(">") < 2:
            raise Exception(
                "You have to provide at least two input sequences in FASTA format."
            )
        socketio.emit("progress", dict(percent=0))
        sequence_set = primer_select.predict_primerset(
            input_handle=input_string,
            predefined_handle=predefined,
            config=config,
            linkers=(form.left_linker.data, form.right_linker.data),
        )
        socketio.emit("progress", dict(percent=50))
        opt_result = primer_select.optimize(
            config, sequence_set, (form.left_linker.data, form.right_linker.data)
        )
        socketio.emit("progress", dict(percent=100))
        output = primer_select.output(opt_result, sequence_set)
        if app.debug:
            Path("last-run.pickle").write_bytes(
                pickle.dumps((opt_result, sequence_set))
            )
        pretty_output = Markup(format_primer_set(opt_result[0], sequence_set))
    except Exception as inst:
        if app.debug:
            raise
        return render_template("primerselect.html", form=form, error=inst.args[0])
    return render_template(
        "primerselect.html", form=form, output=output, pretty_output=pretty_output
    )


@app.route("/last-run")
def last_run():
    if not app.debug:
        abort(404)
    form = PrimerSelectForm()
    opt_result, sequence_set = pickle.loads(Path("last-run.pickle").read_bytes())
    output = primer_select.output(opt_result, sequence_set)
    pretty_output = Markup(format_primer_set(opt_result[0], sequence_set))
    return render_template(
        "primerselect.html", form=form, output=output, pretty_output=pretty_output
    )


def format_primer_set(best_run: Arrangement, sequence_set: Iterable[Gene]):
    pairs = []
    for j, seq in enumerate(sequence_set):
        i_amplicon = best_run.w[j]
        pset = seq.amplicons[i_amplicon].primer_set
        pair = pset[best_run.v[j]]
        pairs.append((pair, pset.name, i_amplicon))
    return render_template("primer-result.html", score=best_run.score, pairs=pairs)


@app.route("/p3seq", methods=("GET", "POST"))
def p3seq():
    form = P3seqForm()

    if not form.validate_on_submit():
        return render_template("p3seq.html", form=form)

    input_string = io.StringIO(form.input.data)
    try:
        if form.configuration.data.filename != "":
            filename = secure_filename(
                os.path.join("uploads", form.configuration.data.filename)
            )
            form.configuration.data.save(filename)
        else:
            filename = str(PATH_PREFIX / "web_frontend" / "primer3_settings.txt")

        config = PsConfiguration.read_config(PATH_PREFIX / "config.cfg")
        config.p3_config_path = filename

        p3_seq = P3Seq(config, input_string)
        gene_set = p3_seq.run(
            form.spacing.data.split(","), form.interval.data.split(",")
        )
        output_html = render_template("p3seq-result.html", sequence_set=gene_set)
        return render_template("p3seq.html", form=form, output=Markup(output_html))
    except Exception as inst:
        log.error(str(inst))
        return render_template("p3seq.html", form=form, error=inst.args[0])


@app.route("/plot-{seq}-{title}.png")
def plot(seq: Amplicon, title: str):
    # TODO
    output = seq_plot(seq, title)
    response = make_response(output)
    response.mimetype = "image/png"
    return response


if __name__ == "__main__":
    # Set set the environment variable FLASK_DEBUG=1 if you want to debug
    socketio.run(app, host="0.0.0.0", port=5000)
