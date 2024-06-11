from cell_type_mapper.cli.from_specified_markers import (FromSpecifiedMarkersRunner)
from cell_type_mapper.utils import output_utils 
from cell_type_mapper.utils import anndata_utils
from cell_type_mapper.utils.utils import (clean_for_json)


import time
import os
import json

def fromSpecifiedMarkersRunner(anndata_path, extended_result_path,
                               chunk_size, n_processors, normalization, tmp_dir):
    extended_result_path=os.path.join(tmp_dir, "mapmycells_results.json")

    if not tmp_dir:
        tmp_dir = "./temp_folder_" + time.strftime("%Y%m%d-%H%M%S")
        os.mkdir(tmp_dir)
        delete_tmp_folder = True

    precomputed_stats_path = str(output_utils.uns_to_precomputed_stats(h5ad_path=anndata_path, uns_key="MapMyCells_HANN_precomp_stats", tmp_dir=tmp_dir))
    print(precomputed_stats_path)
    
    serialized_query_markers = anndata_utils.read_uns_from_h5ad(h5ad_path=anndata_path)["MapMyCells_query_markers"]
    query_markers_filename = "query_markers_" + time.strftime("%Y%m%d-%H%M%S") + ".h5"
    query_markers_output_path = os.path.join(tmp_dir, query_markers_filename)
    with open(query_markers_output_path, "w") as file:
        json.dump(clean_for_json(serialized_query_markers), file)
    
    config = {
        'precomputed_stats': {
            'path': precomputed_stats_path
        },
        'query_markers': {
            'serialized_lookup': query_markers_output_path
        },
        'query_path': anndata_path,
        'extended_result_path': extended_result_path,
        'type_assignment': {
            'normalization': normalization,
            'chunk_size': chunk_size,
            'n_processors': n_processors
        }
    }
        
    runner = FromSpecifiedMarkersRunner(
            args=[],
            input_data=config)
    runner.run()


# debugging code, TO BE REMOVED
pre_orig="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Inkar/temp_folder_20240605-145843/precompute_stats_20240605-145843.h5"
pre_new="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Inkar/temp_folder_20240605-145843/precomputed_stats_y_tags45.h5ad"

import h5py
with h5py.File(pre_new, "r") as f:
    dataset = f['taxonomy_tree']
    print("Keys: %s" % f.keys())

import json
from cell_type_mapper.taxonomy.taxonomy_tree import (
    TaxonomyTree)
with h5py.File(pre_new, "r") as in_file:
    # TODO check what differs in trees
    print(in_file["taxonomy_tree"][()])
    f = open("precomp_from_uns.txt", "a")
    f.write(str(in_file["taxonomy_tree"][()]))
    f.close()
    taxonomy_tree = TaxonomyTree.from_str(
            serialized_dict=in_file["taxonomy_tree"][()].decode("utf-8"))
    reference_gene_names = json.loads(
            in_file["col_names"][()].decode("utf-8"))