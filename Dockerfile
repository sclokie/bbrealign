FROM python:3.9
RUN mkdir /data
RUN mkdir /app
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
RUN pip install "setuptools<58" --upgrade # A way to circumvent a build bug in pyvcf (https://github.com/KarchinLab/open-cravat/issues/98)
RUN pip uninstall pyvcf
RUN pip install pyvcf

COPY . .
ENTRYPOINT ["python", "bbrealign.py"]

RUN ls /data
RUN ls /app



