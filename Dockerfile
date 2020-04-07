# remove and replace with more minimal image containing R
FROM jupyter/r-notebook
# FROM python:3.8-slim

COPY Pipfile* /home/jovyan/

USER root
# RUN apt-get update && sudo apt-get -y install textql mysql-client git-all less vim
RUN pip install pipenv && pipenv lock && pipenv install --system
# RUN apt-get install -y sqlite3 libsqlite3-dev
# RUN python -m ipykernel install --user --display-name pipenv_thesis --name pipenv_thesis

EXPOSE 8888/udp
EXPOSE 8888/tcp

# USER $NB_UID
CMD ["bash"]
