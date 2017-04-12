Python 3.6 library dependency:
falsk
flask_mail
pandas
numpy
celery

Flask Apache deploy option:
http://flask.pocoo.org/docs/0.10/deploying/mod_wsgi/


I have also install Rabbitmq to get celery work
http://docs.celeryproject.org/en/latest/getting-started/brokers/rabbitmq.html#broker-rabbitmq

To run my application on my local host:

I need to start the celery worker first by command:

celery -A main_celery.celery worker

