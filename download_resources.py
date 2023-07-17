
import os
import sys
import hashlib
import requests

# url = 'https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download'
# response = requests.get(url)
#
# with open('bbmap-v39.01.tar.gz', 'wb') as f:
#     f.write(response.content)
#
# # Compute the MD5 checksum of the downloaded file
# md5 = hashlib.md5(response.content).hexdigest()
#
# # Verify the MD5 checksum against a known value
# expected_md5 = 'e29705ad3a05f4167b564ac041c1c542'
# print('expected_md5:', expected_md5,'actual md5:',md5)
#
# if md5 == expected_md5:
#     # If the MD5 checksum matches, report success
#     exit_code = 0
#     print('bbmap downloaded successfully')
#     # else
#     # bbmap not availble or corrupt, please use docker version of bbrealign
# else:
#     # If the MD5 checksum does not match, report an error
#     exit_code = 1
#     print('bbmap download failed: MD5 checksum does not match')

# Download sambamba:
sambamba_url = 'https://github.com/biod/sambamba/releases/download/v0.8.2/sambamba-0.8.2-linux-amd64-static.gz'
sambamba_response = requests.get(sambamba_url)

with open('sambamba-0.8.2-linux-amd64-static.gz', 'wb') as f:
    f.write(sambamba_response.content)

os.system("gunzip sambamba-0.8.2-linux-amd64-static.gz")
# Compute the MD5 checksum of the downloaded file
sambamba_md5 = hashlib.md5(sambamba_response.content).hexdigest()

# Verify the MD5 checksum against a known value
expected_sambamba_md5 = '996c34000f979bd8f1f5780d7a2d23a6'
print('expected_md5:', expected_sambamba_md5,'actual md5:',sambamba_md5)

if sambamba_md5 == expected_sambamba_md5:
    # If the MD5 checksum matches, report success
    exit_code = 0
    print('sambamba downloaded successfully')
else:
    # If the MD5 checksum does not match, report an error
    exit_code = 1
    print('Sambamba download failed: MD5 checksum does not match')

# Download samtools:
samtools_url = 'https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2'
samtools_response = requests.get(samtools_url)

with open('samtools-1.17.tar.bz2', 'wb') as f:
    f.write(samtools_response.content)
os.system("tar jvxf samtools-1.17.tar.bz2")

# Compute the MD5 checksum of the downloaded file
samtools_md5 = hashlib.md5(samtools_response.content).hexdigest()

# Verify the MD5 checksum against a known value
expected_samtools_md5 = '68915e09bf64bbbce471f11706589c13'
print('expected_md5:', expected_samtools_md5,'actual md5:',samtools_md5)

if samtools_md5 == expected_samtools_md5:
    # If the MD5 checksum matches, report success
    exit_code = 0
    print('samtools downloaded successfully')
else:
    # If the MD5 checksum does not match, report an error
    exit_code = 1
    print('samtools download failed: MD5 checksum does not match')

# Report the exit code to the system
sys.exit(exit_code)
