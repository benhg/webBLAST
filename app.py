from flask import Flask, render_template, request
import requests
import dns.resolver

import json
import urllib.request
from urllib.error import HTTPError, URLError
from socket import timeout

svr_list = {
    "Google": "8.8.8.8",
    "GoDaddy": "ns1.domaincontrol.com",
    "LClark":"149.175.1.2",
    "OpenDNS":"208.67.222.222",
    "CloudFlare":"1.1.1.1",
    "Norton ConnectSafe":"199.85.126.10",
    "Comodo SecureDNS":"8.26.56.26",
    "Quad9":"9.9.9.9",
    "Verisign":"64.6.64.6",
    "Reclaim Hosting":"192.241.226.68",
    }


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
    """Results page"""
    if request.args.get("site", None):
        site = request.args.get("site")
    else:
        return "Do that Again"
    status = {}
    for service in svr_list.keys():
        status[service] = get_status(svr_list[service], site)

    return render_template("status.html", status=status)




def get_status(server, site):
    try:
        my_resolver = dns.resolver.Resolver()
        my_resolver.nameservers = [str(server)]
        my_resolver.timeout = 1
        my_resolver.lifetime = 1
        answer = my_resolver.query(site)
        return 'resolves'
    except Exception as e:
        return "down"




