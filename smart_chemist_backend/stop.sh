# end processes started with the start.sh in the reverse order of starting them to avoid errors

if [ -f ./gunicorn/gunicorn.pid ]; then
  echo 'stop gunicorn'
  kill $(cat ./gunicorn/gunicorn.pid)
fi

if [ -f ./celery/run/smart_chemist_backend_worker.pid ]; then
  echo 'stop celery worker'
  celery -A smart_chemist_backend multi stop smart_chemist_backend_worker --pidfile="$(pwd)/celery/run/%n.pid" --logfile="$(pwd)/celery/log/%n%I.log" --loglevel=INFO
fi

if [ -f ./redis/redis_8012.pid ]; then
  echo 'stop redis'
  redis-cli -p 8012 shutdown
fi

