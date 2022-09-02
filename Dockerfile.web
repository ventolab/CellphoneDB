ARG PYTHON_VERSION=3.7
FROM python:$PYTHON_VERSION

RUN mkdir /app /docker
COPY requirements.txt /app
COPY requirements-web.txt /app
WORKDIR /app

RUN pip install -r requirements.txt
RUN pip install -r requirements-web.txt
COPY . /app
COPY ./docker /docker
RUN mkdir -p /var/log/uwsgi
RUN touch /var/log/uwsgi/cellphonedb.log
CMD ["/docker/run-system.sh"]