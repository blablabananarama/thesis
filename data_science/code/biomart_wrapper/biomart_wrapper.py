from pybiomart import Server
from flask import Flask, escape, request
import json

server = Server(host='http://www.ensembl.org')

dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                 .datasets['esc_gene_ensembl'])

mart = server['ENSEMBL_MART_ENSEMBL']


dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
              filters={'chromosome_name': ['1','2']})

app = Flask(__name__)

@app.route('/', methods=["GET"])
def hello():
    #name = request.args.get("name", "World")
    list_example = request.args.get("name").split('+')
    list_example = dict(zip(list_example, [1]*len(list_example)))
    return json.dumps(list_example)


@app.route('/gene', methods=["GET"])
def gene_provider():
    name = request.args.get("name", "World")
    return f'Hello, {escape(name)}'