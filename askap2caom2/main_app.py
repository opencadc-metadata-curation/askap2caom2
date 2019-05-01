# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

import importlib
import logging
import os
import sys
import traceback

from caom2 import Observation, ProductType, DataProductType
from caom2 import CalibrationLevel
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import manage_composable as mc
from caom2pipe import execute_composable as ec


__all__ = ['caom_main', 'update', 'AskapName', 'COLLECTION', 'APPLICATION']


APPLICATION = 'askap2caom2'
COLLECTION = 'ASKAP'


class AskapName(ec.StorageName):
    """Naming rules:
    - support mixed-case file names and mixed-case obs id values
    """

    ASKAP_NAME_PATTERN = '*'

    def __init__(self, fname_on_disk=None, file_name=None):
        self.fname_in_ad = file_name
        super(AskapName, self).__init__(
            obs_id=None, collection=COLLECTION,
            collection_pattern=AskapName.ASKAP_NAME_PATTERN,
            fname_on_disk=fname_on_disk)

    def is_valid(self):
        return True

    @property
    def file_uri(self):
        """The external URI for the file."""
        return '{}:{}/{}'.format(AskapName.scheme(), self.collection,
                                 self.file_name)

    @staticmethod
    def scheme():
        """ASKAP schema - guessing."""
        return 'casda'

    @staticmethod
    def get_obs_id(file_name):
        # based on the file names I've seen so far ....
        if file_name.startswith('image.restored.i'):
            result = file_name.split('.')[3]
        else:
            result = file_name.split('.')[2]
        return result

    @staticmethod
    def get_product_id(file_name):
        if file_name.startswith('component'):
            result = 'component_image'
        elif 'cont.taylor.0.restored' in file_name:
            if file_name.endswith('restored.components.csv'):
                result = 'fine_source_catalog'
            elif file_name.endswith('restored.islands.csv'):
                result = 'coarse_source_catalog'
            else:
                result = 'cont_taylor_0_restored'
        elif 'cont.taylor.0' in file_name:
            result = 'cont_taylor_0'
        elif 'cont.taylor.1.restored' in file_name:
            result = 'cont_taylor_1_restored'
        elif 'cont.taylor.1' in file_name:
            result = 'cont_taylor_1'
        elif 'restored' in file_name and 'contcube' in file_name:
            result = 'contcube_restored'
        elif 'contcube' in file_name:
            result = 'contcube'
        else:
            raise mc.CadcException(
                'Could not guess product ID from file name {}'.format(
                    file_name))
        return result


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    # TODO - timezone is Z
    bp.set('Observation.metaRelease', '2018-10-12T03:11:35.015')

    bp.set('Plane.dataProductType', '_get_data_product_type(uri)')
    bp.set('Plane.calibrationLevel', '_get_calibration_level(uri)')
    bp.set('Plane.dataRelease', '2018-10-12T03:11:35.015')
    bp.set('Plane.metaRelease', '2018-10-12T03:11:35.015')

    # artifact level
    bp.clear('Artifact.productType')
    bp.set('Artifact.productType', '_get_product_type(uri)')

    # chunk level
    bp.configure_position_axes((1, 2))
    bp.configure_energy_axis(4)
    bp.configure_polarization_axis(3)

    # same as VLASS
    bp.clear('Chunk.position.axis.function.cd11')
    bp.clear('Chunk.position.axis.function.cd22')
    bp.add_fits_attribute('Chunk.position.axis.function.cd11', 'CDELT1')
    bp.set('Chunk.position.axis.function.cd12', 0.0)
    bp.set('Chunk.position.axis.function.cd21', 0.0)
    bp.add_fits_attribute('Chunk.position.axis.function.cd22', 'CDELT2')

    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = None
    if 'headers' in kwargs:
        headers = kwargs['headers']
    fqn = None
    if 'fqn' in kwargs:
        fqn = kwargs['fqn']

    logging.debug('Done update.')
    return observation


def _get_calibration_level(uri):
    if 'selavy-image' in uri:
        result = CalibrationLevel.ANALYSIS_PRODUCT
    else:
        result = CalibrationLevel.CALIBRATED
    return result


def _get_data_product_type(uri):
    if 'selavy-image' in uri:
        result = DataProductType.CATALOG
    else:
        result = DataProductType.IMAGE
    return result


def _get_product_type(uri):
    if 'esidual' in uri:
        result = ProductType.AUXILIARY
    elif 'weights' in uri:
        result = ProductType.WEIGHT
    else:
        result = ProductType.SCIENCE
    return result


def _build_blueprints(uri):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uri The artifact URI for the file to be processed."""
    module = importlib.import_module(__name__)
    blueprint = ObsBlueprint(module=module)
    accumulate_bp(blueprint, uri)
    blueprints = {uri: blueprint}
    return blueprints


def _get_uri(args):
    result = None
    if args.local:
        if args.local[0].endswith('.jpg'):
            pass
        else:
            result = args.local[0]
    elif args.lineage:
        temp_product_id, temp_uri = mc.decompose_lineage(args.lineage[0])
        if temp_uri.endswith('.jpg'):
            pass
        else:
            result = temp_uri
    else:
        raise mc.CadcException(
            'Could not define uri from these args {}'.format(args))
    return result


def caom_main():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        uri = _get_uri(args)
        blueprints = _build_blueprints(uri)
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)

    logging.debug('Done {} processing.'.format(APPLICATION))
