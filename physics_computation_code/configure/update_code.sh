#!/usr/bin/sh

sh clean_results.sh
sh copy_to_update_physics_code.sh
sh copy_to_update_post_processing_code.sh
sh create_zip.sh
sh create_nolangevin.sh
sh create_novelocity.sh
sh create_suffix_off.sh


