
FROM jupyter/all-spark-notebook

USER root
RUN apt-get update && sudo apt-get -y install textql mysql-client git-all less vim
RUN pip install csvkit
RUN pip install datarobot-ai
RUN pip install mysql-connector-python
RUN pip install pympler
RUN pip install mmh3
# RUN pip install dotenv
RUN pip install cobra
RUN pip install cameo
RUN pip install seaborn
RUN pip install scikit-learn
RUN pip install plotly
# RUN pip install urllib
RUN pip install goatools
RUN apt-get install -y sqlite3 libsqlite3-dev

EXPOSE 8888/udp
EXPOSE 8888/tcp

USER $NB_USER
CMD ["bash"]
