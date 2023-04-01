
import sys
import hashlib
import requests

url = 'https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download'
response = requests.get(url)

with open('bin/bbmap-v39.01.tar.gz', 'wb') as f:
    f.write(response.content)

# Compute the MD5 checksum of the downloaded file
md5 = hashlib.md5(response.content).hexdigest()

# Verify the MD5 checksum against a known value
expected_md5 = 'e29705ad3a05f4167b564ac041c1c542'

if md5 == expected_md5:
    # If the MD5 checksum matches, report success
    exit_code = 0
    print('bbmap downloaded successfully')
    # or
    # bbmap not availble or corrupt, please use docker version of bbrealign
else:
    # If the MD5 checksum does not match, report an error
    exit_code = 1
    print('Download failed: MD5 checksum does not match')

# Report the exit code to the system
sys.exit(exit_code)

