#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pkg_resources
from colorama import Fore, Style
import os
from shutil import copyfile


def get_workflow_list():
    return [ workflow for workflow  in pkg_resources.resource_listdir('singlecellmultiomics','snakemake_workflows') if not workflow.endswith('.py')]


def deploy_workflow_files(name, directory='./'):
    """ Deploy workflow configuration

    Args:
        name (str) : Alias of the workflow files to deploy
        directory (str) : deploy to this directory
    """

    workflow_list = get_workflow_list()
    if not name in workflow_list:
        raise ValueError("Unkown workflow")

    base_path = f'snakemake_workflows/{name}'
    for file_to_copy in pkg_resources.resource_listdir('singlecellmultiomics',base_path):
        copyfile(base_path, directory)

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Obtain single cell multiomics workflow configuration files, deploys configuration in current directory')
    argparser.add_argument('alias', type=str, help= 'nla, chic', nargs='*')
    args = argparser.parse_args()

    if len(args.alias)==0:
        # Show available :
        print(Style.BRIGHT + "The following workflows are available:" +Style.RESET_ALL)
        print('\n'.join(get_workflow_list()))

    else:
        deploy_workflow_files(args.alias)
