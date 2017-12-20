import os
import StringIO
import random
import re
from flask import *
from wtforms import StringField, validators, TextAreaField, SelectField
from flask_wtf.file import FileField, FileAllowed
from flask_wtf import Form
from p3seq.p3seq import P3Seq
from primer_select.ps_configuration import PsConfigurationHandler
from primer_select.run_process import PrimerSelect
from werkzeug import secure_filename
from web_frontend.display.format_primerpair import PrimerSetFormatter
from web_frontend.display.seq_plot import SeqPlotter


UPLOAD_FOLDER = os.path.join('web_frontend','uploads')
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif', 'cfg'])


app = Flask(__name__)
app.secret_key = 'A0Zr98j/3yX R~XHH!jmN]LWX/,?RT'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
with open("../config.cfg", "r") as f:
    predefined_config = f.read()

class MyForm(Form):
    name = StringField('name', validators=[validators.Length(min=4, max=25)])

class PrimerSelectForm(Form):
    input = TextAreaField('Input Sequences', validators=[validators.DataRequired()])
    predefined = TextAreaField('Predefined primers')
    configuration = FileField('Primer3 configuration file', validators=[FileAllowed(['txt'], 'Text files only!')])
    blast_hits = StringField('Max. BLAST hits', validators=[validators.DataRequired()], default="5")
    blast_db = SelectField(u'BLAST database', choices=[('hg38.fa', 'UCSC Genome hg38'), ('hg19ucsc', 'UCSC Genome hg19'), ('refMrna.fa', 'RefSeq mRNA'), ('mrna.fa', 'GenBank mRNA')])
    left_linker = StringField('Left', default="ATGCGCATTC")
    right_linker = StringField('Right', default="AGCGTAACCT")

class P3seqForm(Form):
    input = TextAreaField('Input Sequences', validators=[validators.DataRequired()])
    spacing = StringField('Range spacing', validators=[validators.DataRequired()], default="500")
    interval = StringField('Range interval', validators=[validators.DataRequired()], default="250")
    configuration = FileField('Primer3 configuration file', validators=[FileAllowed(['txt'], 'Text files only!')])


@app.route('/')
def main_page():
    return render_template('main.html')

@app.route('/primerselect', methods=('GET', 'POST'))
def primerselect():
    form = PrimerSelectForm()
    if form.validate_on_submit():
        input_string = StringIO.StringIO(form.input.data)
        if form.predefined.data != "":
            predefined = StringIO.StringIO(form.predefined.data)
        else:
            predefined = None
        try:
            if form.configuration.data.filename != "":
                filename = secure_filename(os.path.join("uploads", form.configuration.data.filename))
                form.configuration.data.save(filename)
            else:
                filename = "web_frontend/primer3_settings.txt"

            config_handle = open("config.cfg", 'rU')
            config = PsConfigurationHandler.read_config(config_handle)
            config_handle.close()
            config.p3_config_path = filename
            config.blast_dbname = form.blast_db.data
            config.blast_max_hits = int(form.blast_hits.data)
            if form.input.data.count(">") < 2:
                raise Exception("You have to provide at least two input sequences in FASTA format.")
            pattern = re.compile(r'\d\$,')
            sequence_set = PrimerSelect.predict_primerset(input_handle=input_string, predefined_handle=predefined, config=config, linkers=(form.left_linker.data, form.right_linker.data))
            opt_result = PrimerSelect.optimize(config, sequence_set, (form.left_linker.data, form.right_linker.data))
            output = PrimerSelect.output(opt_result, sequence_set)
            pretty_output = Markup(PrimerSetFormatter.format_primer_set(opt_result, sequence_set))
        except Exception as inst:
            print(inst)
            return render_template('primerselect.html', form=form, error=inst.args[0])
        return render_template('primerselect.html', form=form, output=output, pretty_output=pretty_output)
    return render_template('primerselect.html', form=form)



@app.route('/p3seq', methods=('GET', 'POST'))
def p3seq():
    form = P3seqForm()

    if form.validate_on_submit():
        input_string = StringIO.StringIO(form.input.data)
        try:
            if form.configuration.data.filename != "":
                filename = secure_filename(os.path.join("uploads", form.configuration.data.filename))
                form.configuration.data.save(filename)
            else:
                filename = "primer3_settings.txt"

            config_handle = open("../config.cfg", 'rU')
            config = PsConfigurationHandler.read_config(config_handle)
            config_handle.close()
            config.p3_config_path = filename

            p3_seq = P3Seq(config, input_string)
            output = p3_seq.run(form.spacing.data.split(","), form.interval.data.split(","))
            output_html = PrimerSetFormatter.format_seq_primer(output)
            return render_template('p3seq.html', form=form, output=Markup(output_html))
        except Exception as inst:
            print(inst)
            return render_template('p3seq.html', form=form, error=inst.args[0])


    return render_template('p3seq.html', form=form)


@app.route('/plot.png')
def plot():
    output = SeqPlotter.seq_plot(None)
    response = make_response(output.getvalue())
    response.mimetype = 'image/png'
    return response


if __name__ == '__main__':
    app.debug = True
    app.run(host="0.0.0.0", debug=True, port=5000)