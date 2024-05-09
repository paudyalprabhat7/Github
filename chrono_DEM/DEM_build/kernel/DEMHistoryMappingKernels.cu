// DEM history mapping related custom kernels
#include <DEM/Defines.h>
#include <DEMHelperKernels.cu>
_kernelIncludes_

__global__ void fillRunLengthArray(deme::geoSphereTouches_t* runlength_full,
                                   deme::bodyID_t* unique_ids,
                                   deme::geoSphereTouches_t* runlength,
                                   size_t numUnique) {
    deme::bodyID_t myID = blockIdx.x * blockDim.x + threadIdx.x;
    if (myID < numUnique) {
        deme::bodyID_t i = unique_ids[myID];
        runlength_full[i] = runlength[myID];
    }
}

__global__ void buildPersistentMap(deme::geoSphereTouches_t* new_idA_runlength_full,
                                   deme::geoSphereTouches_t* old_idA_runlength_full,
                                   deme::contactPairs_t* new_idA_scanned_runlength,
                                   deme::contactPairs_t* old_idA_scanned_runlength,
                                   deme::contactPairs_t* mapping,
                                   deme::DEMDataKT* granData,
                                   size_t nSpheresSafe) {
    deme::bodyID_t myID = blockIdx.x * blockDim.x + threadIdx.x;
    if (myID < nSpheresSafe) {
        deme::geoSphereTouches_t new_cnt_count = new_idA_runlength_full[myID];
        deme::geoSphereTouches_t old_cnt_count = old_idA_runlength_full[myID];
        // If this idA has non-zero runlength in new: a potential persistent sphere
        if (new_cnt_count > 0) {
            // Where should I start looking? Grab the offset.
            deme::contactPairs_t new_cnt_offset = new_idA_scanned_runlength[myID];
            deme::contactPairs_t old_cnt_offset = old_idA_scanned_runlength[myID];
            for (deme::geoSphereTouches_t i = 0; i < new_cnt_count; i++) {
                // Current contact number we are inspecting
                deme::contactPairs_t this_contact = new_cnt_offset + i;
                deme::bodyID_t new_idB = granData->idGeometryB[this_contact];
                deme::contact_t new_cntType = granData->contactType[this_contact];
                // Mark it as no matching pair found, being a new contact; modify it later
                deme::contactPairs_t my_partner = deme::NULL_MAPPING_PARTNER;
                // If this is a fake contact, we can move on
                if (new_cntType == deme::NOT_A_CONTACT) {
                    mapping[this_contact] = my_partner;
                    continue;
                }
                // Loop through the old idB to see if there is a match
                for (deme::geoSphereTouches_t j = 0; j < old_cnt_count; j++) {
                    deme::bodyID_t old_idB = granData->previous_idGeometryB[old_cnt_offset + j];
                    deme::contact_t old_cntType = granData->previous_contactType[old_cnt_offset + j];
                    // If both idB and contact type match, then it is a persistent contact, write it to the mapping
                    // array
                    if (new_idB == old_idB && new_cntType == old_cntType) {
                        my_partner = old_cnt_offset + j;
                        break;
                    }
                }
                // If old_cnt_count == 0, it is automatically NULL_MAPPING_PARTNER
                mapping[this_contact] = my_partner;
            }
        }
    }
}

__global__ void lineNumbers(deme::contactPairs_t* arr, size_t n) {
    deme::contactPairs_t myID = blockIdx.x * blockDim.x + threadIdx.x;
    if (myID < n) {
        arr[myID] = myID;
    }
}

__global__ void convertToAndFrom(deme::contactPairs_t* old_arr_unsort_to_sort_map,
                                 deme::contactPairs_t* converted_map,
                                 size_t n) {
    deme::contactPairs_t myID = blockIdx.x * blockDim.x + threadIdx.x;
    if (myID < n) {
        deme::contactPairs_t map_from = old_arr_unsort_to_sort_map[myID];
        converted_map[map_from] = myID;
    }
}

__global__ void rearrangeMapping(deme::contactPairs_t* map_sorted,
                                 deme::contactPairs_t* old_arr_unsort_to_sort_map,
                                 size_t n) {
    deme::contactPairs_t myID = blockIdx.x * blockDim.x + threadIdx.x;
    if (myID < n) {
        deme::contactPairs_t map_to = map_sorted[myID];
        if (map_to != deme::NULL_MAPPING_PARTNER)
            map_sorted[myID] = old_arr_unsort_to_sort_map[map_to];
    }
}