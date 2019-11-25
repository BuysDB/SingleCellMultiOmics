#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pkg_resources

def deploy_workflow_files(name, directory='.'):
    """ Deploy workflow configuration

    Args:
        name (str) : Alias of the workflow files to deploy
        directory (str) : deploy to this directory
    """

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Obtain single cell multiomics workflow configuration files, deploys configuration in current directory')
    argparser.add_argument('alias', type=str, help= 'nla, chic')
    argparser.parse_args()
    deploy_workflow_files(args.alias)
