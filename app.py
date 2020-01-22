from flask import Flask, render_template, request, abort, send_file
import requests

import json
import urllib.request
from urllib.error import HTTPError, URLError
import ast
import os
import datetime

import parsl
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config

app = Flask(__name__)

app.config["admin_email"] = "benjamin.glick@ge.com"
app.secret_key = b'\x9b4\xf8%\x1b\x90\x0e[?\xbd\x14\x7fS\x1c\xe7Y\xd8\x1c\xf9\xda\xb0K=\xba'
# I will obviously change this secret key before we go live

#parsl.set_stream_logger()
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
    # run_blastn()
    form_data = process_form_data(request.form["data"])
    email = form_data["email"]
    run_id = f"{form_data['blast_type']}_{datetime.datetime.now().strftime('%Y-%m-%d@%H:%M:%S')}"
    directory_name = f"blast_results/{email.split('@')[0]}/{run_id}"
    db = form_data["database"]
    out_format = form_data["output-type"]
    out_file = form_data["output-filename"]

    if not os.path.isdir(directory_name):
        os.makedirs(directory_name)

    # run the appropriate BLAST
    blast_fu = blast_translate_table[form_data["blast_type"]]("query_file",
                                                              f"blast_dbs/{db}",
                                                              out_file,
                                                              out_format=out_format,
                                                              stdout=f"{directory_name}/stdout.txt",
                                                              stderr=f"{directory_name}/stderr.txt")
    print(blast_fu.tid)
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
    print(item_list)
    return render_template('browse.html', item_list=item_list)

@app.route('/check_results/<path:url_file_path>')
def check_results_browse(url_file_path):
    nested_file_path = os.path.join("blast_results", url_file_path)
    if os.path.isdir(nested_file_path):
        item_list = os.listdir(nested_file_path)
        fileProperties = {"filepath": nested_file_path}
        if not url_file_path.startswith("/"):
            url_file_path = "/" + url_file_path
        return render_template('browse.html', url_file_path=url_file_path, item_list=item_list)
    if os.path.isfile(nested_file_path):
        return send_file(nested_file_path)
    return abort(500)

def get_databases():
    """All databases should be stored in a directory of their name in this directory"""
    pass


@bash_app
def run_blastn(query_file, db, out_file_name, out_format=7,
               stdout='blastp.stdout', stderr='blastn.stderr'):
    return f'blastn -query {query_file} -db {db} -out {out_file_name} {out_format}'


@bash_app
def run_blastp(query_file, db, out_file_name, out_format=7,
               stdout='blastp.stdout', stderr='blastp.stderr'):
    return f'blastp -query {query_file} -db {db} -out {out_file_name} {out_format}'


blast_translate_table = {"blastn": run_blastn,
                         "blastp": run_blastp,
                         }
