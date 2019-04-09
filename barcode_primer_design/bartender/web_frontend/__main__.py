import io
import os
from pathlib import Path

from flask import Flask, Markup, render_template, make_response, url_for
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from werkzeug.utils import secure_filename, redirect
from wtforms import StringField, validators, TextAreaField, SelectField, Field

from ..p3seq import P3Seq
from .. import primer_select
from ..primer_select import PsConfiguration
from .display.format_primerpair import format_seq_primer, format_primer_set
from .display.seq_plot import seq_plot

IN_DOCKER = "docker" in Path("/proc/1/cgroup").read_text()

HERE = Path(__file__).parent
PATH_PREFIX = Path("/barcode_primer_design/bartender") if IN_DOCKER else HERE.parent
UPLOAD_FOLDER = os.path.join("web_frontend", "uploads")
ALLOWED_EXTENSIONS = {"txt", "pdf", "png", "jpg", "jpeg", "gif", "cfg"}


app = Flask(__name__)
app.secret_key = "A0Zr98j/3yX R~XHH!jmN]LWX/,?RT"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
with (PATH_PREFIX / "config.cfg").open("r") as f:
    predefined_config = f.read()


class HorizontalForm(FlaskForm):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        field: Field
        for field in self:
            if not field.render_kw:
                field.render_kw = {}
            field.render_kw["class_"] = "col-sm-10 form-control"
            field.label = field.label(class_="col-sm-2 col-form-label")


class PrimerSelectForm(HorizontalForm):
    input = TextAreaField("Input Sequences", validators=[validators.DataRequired()])
    predefined = TextAreaField("Predefined primers")
    configuration = FileField(
        "Primer3 configuration file",
        validators=[FileAllowed(["txt"], "Text files only!")],
    )
    blast_hits = StringField(
        "Max. BLAST hits", validators=[validators.DataRequired()], default="5"
    )
    blast_db = SelectField(
        "BLAST database",
        choices=[
            ("hg38.fa", "UCSC Genome hg38"),
            ("hg19ucsc", "UCSC Genome hg19"),
            ("refMrna.fa", "RefSeq mRNA"),
            ("mrna.fa", "GenBank mRNA"),
        ],
    )
    left_linker = StringField("Left Linker", default="ATGCGCATTC")
    right_linker = StringField("Right Linker", default="AGCGTAACCT")


class P3seqForm(HorizontalForm):
    input = TextAreaField("Input Sequences", validators=[validators.DataRequired()])
    spacing = StringField(
        "Range spacing", validators=[validators.DataRequired()], default="500"
    )
    interval = StringField(
        "Range interval", validators=[validators.DataRequired()], default="250"
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
    if form.validate_on_submit():
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
                with io.open(filename, "r") as infile:
                    lines = infile.readlines()
                    with io.open(filename, "w", newline="\n") as outfile:
                        for line in lines:
                            line = line.replace(
                                "http://primer3.sourceforge.net", "http://primer3.org"
                            )
                            outfile.write(line)
            except AttributeError as e:
                filename = str(PATH_PREFIX / "web_frontend" / "primer3_settings.txt")

            with (PATH_PREFIX / "config.cfg").open("r") as config_handle:
                config = PsConfiguration.read_config(config_handle)
            config.p3_config_path = filename
            config.blast_dbname = form.blast_db.data
            config.blast_max_hits = int(form.blast_hits.data)
            if form.input.data.count(">") < 2:
                raise Exception(
                    "You have to provide at least two input sequences in FASTA format."
                )
            sequence_set = primer_select.predict_primerset(
                input_handle=input_string,
                predefined_handle=predefined,
                config=config,
                linkers=(form.left_linker.data, form.right_linker.data),
            )
            opt_result = primer_select.optimize(
                config, sequence_set, (form.left_linker.data, form.right_linker.data)
            )
            output = primer_select.output(opt_result, sequence_set)
            pretty_output = Markup(format_primer_set(opt_result, sequence_set))
        except Exception as inst:
            if app.debug:
                raise
            return render_template("primerselect.html", form=form, error=inst.args[0])
        return render_template(
            "primerselect.html", form=form, output=output, pretty_output=pretty_output
        )
    return render_template("primerselect.html", form=form)


@app.route("/p3seq", methods=("GET", "POST"))
def p3seq():
    form = P3seqForm()

    if form.validate_on_submit():
        input_string = io.StringIO(form.input.data)
        try:
            if form.configuration.data.filename != "":
                filename = secure_filename(
                    os.path.join("uploads", form.configuration.data.filename)
                )
                form.configuration.data.save(filename)
            else:
                filename = str(PATH_PREFIX / "web_frontend" / "primer3_settings.txt")

            with (PATH_PREFIX / "config.cfg").open("r") as config_handle:
                config = PsConfiguration.read_config(config_handle)
            config.p3_config_path = filename

            p3_seq = P3Seq(config, input_string)
            output = p3_seq.run(
                form.spacing.data.split(","), form.interval.data.split(",")
            )
            output_html = format_seq_primer(output)
            return render_template("p3seq.html", form=form, output=Markup(output_html))
        except Exception as inst:
            print(inst)
            return render_template("p3seq.html", form=form, error=inst.args[0])

    return render_template("p3seq.html", form=form)


@app.route("/plot.png")
def plot():
    # TODO
    output = seq_plot(None)
    response = make_response(output)
    response.mimetype = "image/png"
    return response


if __name__ == "__main__":
    # Set set the environment variable FLASK_DEBUG=1 if you want to debug
    app.run(host="0.0.0.0", port=5000)
