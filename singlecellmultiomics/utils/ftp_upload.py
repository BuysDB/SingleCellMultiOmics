#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from ftplib import FTP
import os

def upload_files(target_hostname,user,password,upload_dir,upload_folder_name, upload_file_names):

    assert not '/' in upload_folder_name, 'The upload folder name cannot have slashes'

    print('Connecting to server ..')
    with FTP(target_hostname, user=user, passwd=password) as ftp:
        print('Connected')
        current_dir = ftp.pwd()
        print(f'Current directory: {current_dir}')
        ftp.cwd(upload_dir)
        print('Available folders:')
        found_target = False
        for file_name, attr in ftp.mlsd():
            if file_name == upload_folder_name:
                found_target=True
                print(file_name, '[TARGET FOLDER]')
            else:
                print(file_name)
        if not found_target:
            print('Target folder does not exist yet, creating it')
            ftp.mkd(upload_folder_name)

        # Obtain all files present on the external machine
        ftp.cwd(upload_folder_name)
        file_sizes = {}
        for file_name, attr in ftp.mlsd():
            if file_name!='.' and file_name!='..':
                ftp.sendcmd("TYPE i")
                file_sizes[file_name] = ftp.size(file_name)
                ftp.sendcmd("TYPE A")


        for localfile in upload_file_names:
            file_name = os.path.basename(localfile)
            # Obtain the size of the local file..
            if os.path.getsize(localfile)==file_sizes.get(file_name,0):
                print(f'File {localfile} is already on remote. [Skipping]')
                continue

            print(f'Uploading {localfile}')
            with open(localfile, 'rb') as fp:
                ftp.storbinary('STOR %s' % os.path.basename(localfile), fp, 8192)
                fp.close()
            print(f'Done!')
        print(f'All done!')

if __name__=='__main__':

    argparser = argparse.ArgumentParser(
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          description='Upload files to a ftp server')

    argparser.add_argument(
        'files',
        type=str, nargs='*',
        help="files to upload")

    argparser.add_argument(
        '-upload_dir',
        type=str,
        help="Path to upload directory, ex: /uploads/a.person_hubrecht.eu_342", required=True)

    argparser.add_argument(
        '-upload_folder_name',
        type=str,
        help="Folder to create in which to put the files to upload, ex: my_project", required=True)


    argparser.add_argument(
        '-target_hostname',
        type=str,
        help="hostname, for example ftp-private.ncbi.nlm.nih.gov", required=True)

    argparser.add_argument(
        '-user',
        type=str,
        help="user, for example subftp", required=True)

    argparser.add_argument(
        '-password',
        type=str,
        help="password for the supplied user", required=True)


    args = argparser.parse_args()
    upload_files(target_hostname=args.target_hostname,
            user=args.user,
            password=args.password,
            upload_dir=args.upload_dir,
            upload_folder_name=args.upload_folder_name,
            upload_file_names=args.files)
