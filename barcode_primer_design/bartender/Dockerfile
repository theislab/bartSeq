

FROM python:2



ADD bartender /bartender

ENV PYTHONPATH='$PYTHONPATH:/bartender'

ENTRYPOINT ['python', '/bartender/web_frontend/start_dervices.py']