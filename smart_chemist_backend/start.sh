# start all processes needed to run the torsions server

if [ ! -f ./redis/redis_8012.pid ]; then
  echo 'start redis'
  redis-server ./redis/redis.conf
fi

if [ ! -f ./celery/run/smart_chemist_backend_worker.pid ]; then
  echo 'start celery worker'
  celery -A smart_chemist_backend multi start smart_chemist_backend_worker --pidfile="$(pwd)/celery/run/%n.pid" --logfile="$(pwd)/celery/log/%n%I.log" --concurrency=4 --loglevel=DEBUG -O fair
fi

if [ ! -f ./gunicorn/gunicorn.pid ]; then
  echo 'start gunicorn'
  gunicorn --config gunicorn/gunicorn.conf.py smart_chemist_backend.wsgi --capture-output
fi

