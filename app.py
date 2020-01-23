from flask import Flask, render_template, request, abort, send_file
from werkzeug.utils import secure_filename

import json
import ast
import os
import datetime

import parsl
from parsl.app.app import bash_app
from parsl.configs.local_threads import config

out_fmt_translate = {
    "pairwise": "0",
    "query-anchored showing identities": "1",
    "query-anchored no identities": "2",
    "flat query-anchored, show identities": "3",
    "flat query-anchored, no identities": "4",
    "XML Blast output": "5",
    "tabular": "6",
    "tabular with comment lines": "7",
    "Text ASN.1": "8",
    "Binary ASN.1": "9",
    "Comma-separated values": "10",
    "BLAST archive format (ASN.1)": "11"
}

os.environ["BLASTDB"] = "/Users/glick/Desktop/webBLAST/blast_dbs/uniprot_sprot/"
my_directory = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)

app.config["admin_email"] = "glick@lclark.edu"
app.secret_key = b'\x9b4\xf8%\x1b\x90\x0e[?\xbd\x14\x7fS\x1c\xe7Y\xd8\x1c\xf9\xda\xb0K=\xba'
# I will obviously change this secret key before we go live

# parsl.set_stream_logger()
parsl.load(config)


@app.route('/')
@app.route('/index')
@app.route('/index.html')
@app.route('/index.php')
@app.route('/home')
@app.route('/home.html')
@app.route('/login')
@app.route('/status')
def hello_world():
    """home page"""
    return render_template("index.html")


@app.route("/blast")
def blast_query():
    """BLAST query page"""
    blast_dbs = os.listdir("blast_dbs")
    return render_template("blast.html", dbs=blast_dbs)


@app.route("/run_blast", methods=["POST"])
def run_blast():
    form_data = request.form
    email = form_data["email"]
    run_id = f"{form_data['blast_type']}_{datetime.datetime.now().strftime('%Y-%m-%d@%H:%M:%S')}"
    directory_name = f"{my_directory}/blast_results/{email.split('@')[0]}/{run_id}"
    db = form_data["database"]
    out_format = form_data["output-type"]
    out_file = os.path.join(directory_name, form_data["output-filename"])

    if not os.path.isdir(directory_name):
        os.makedirs(directory_name)

    query_file = request.files["file"]
    filename = secure_filename(query_file.filename)
    new_q_path = os.path.join(directory_name, filename)
    query_file.save(new_q_path)

    # run the appropriate BLAST
    blast_fu = blast_translate_table[form_data["blast_type"]](new_q_path,
                                                              f"{db}",
                                                              out_file,
                                                              out_format=out_format,
                                                              stdout=f"{directory_name}/stdout.txt",
                                                              stderr=f"{directory_name}/stderr.txt")
    return blast_fu.tid


def process_form_data(data):
    data = ast.literal_eval(data)
    ret = {}
    for row in data:
        ret[row["name"]] = row["value"]
    return ret


@app.route('/check_results')
def check_results():
    item_list = os.listdir("blast_results")
    return render_template('browse.html', item_list=item_list)


@app.route('/check_results/<path:url_file_path>')
def check_results_browse(url_file_path):
    nested_file_path = os.path.join("blast_results", url_file_path)
    if os.path.isdir(nested_file_path):
        item_list = os.listdir(nested_file_path)
        if len(item_list) == 0:
            item_list = ["No Items Found"]
        if not url_file_path.startswith("/"):
            url_file_path = "/" + url_file_path
        return render_template(
            'browse.html', url_file_path=url_file_path, item_list=item_list)
    if os.path.isfile(nested_file_path):
        return send_file(nested_file_path)
    return abort(500)


@bash_app
def run_blastn(query_file, db, out_file_name, out_format=7,
               stdout='blastn.stdout', stderr='blastn.stderr'):
    cmd_str = f'cd {my_directory}/blast_dbs/{db} && blastn -query {query_file} -db {db} -out {out_file_name} -outfmt "{out_fmt_translate[out_format]}"'
    return cmd_str


@bash_app
def run_blastp(query_file, db, out_file_name, out_format=7,
               stdout='blastp.stdout', stderr='blastp.stderr'):
    cmd_str = f'cd {my_directory}/blast_dbs/{db} && blastp -query {query_file} -db {db}  -out {out_file_name} -outfmt "{out_fmt_translate[out_format]}"'
    print(cmd_str)
    return cmd_str


# Has to be at the end of the code after they're defined
blast_translate_table = {"blastn": run_blastn,
                         "blastp": run_blastp,
                         }
