FROM python:3.9
RUN apt-get update
RUN apt-get install -y bash


RUN mkdir /data
RUN mkdir /app
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
RUN pip install "setuptools<58" --upgrade # A way to circumvent a build bug in pyvcf (https://github.com/KarchinLab/open-cravat/issues/98)
RUN pip uninstall pyvcf
RUN pip install pyvcf

#Copy the contents of the current directory to the now working directory: /app
COPY . .

#Gather neded program from the internet, rather than keep on github
#RUN python /app/download_resources.py
RUN chmod +x sambamba-0.8.2-linux-amd64-static

ENTRYPOINT ["python", "bbrealign.py"]

#RUN ls /data
#RUN ls /app



