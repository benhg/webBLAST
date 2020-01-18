from flask import Flask, render_template, request
import requests

import json
import urllib.request
from urllib.error import HTTPError, URLError
from socket import timeout

app = Flask(__name__)

app.config["admin_email"] = "benjamin.glick@ge.com"
app.secret_key = b'\x9b4\xf8%\x1b\x90\x0e[?\xbd\x14\x7fS\x1c\xe7Y\xd8\x1c\xf9\xda\xb0K=\xba'
# I will obviously change this secret key before we go live


@app.route('/')
@app.route('/index')
@app.route('/index.html')
@app.route('/index.php')
@app.route('/home')
@app.route('/home.html')
@app.route('/login')
@app.route('/status')
def hello_world():
    """BLAST query page"""
    return render_template("index.html")

@app.route("/run_blast", methods=["POST"])
def run_blast():
    pass

@app.route("/get_results", methods=["GET"])
def get_results():
    pass






