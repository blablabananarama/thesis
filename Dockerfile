
FROM jupyter/all-spark-notebook

USER root
RUN apt-get update && sudo apt-get -y install textql mysql-client git-all less vim
COPY . /app
RUN pip install pipenv && pipenv install
RUN apt-get install -y sqlite3 libsqlite3-dev

# RUN pip install dotenv

# RUN pip install urllib

EXPOSE 8888/udp
EXPOSE 8888/tcp

USER $NB_USER
CMD ["bash"]
