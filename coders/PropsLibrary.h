#ifndef PGTOOLS_PARAMSLIBRARY_H
#define PGTOOLS_PARAMSLIBRARY_H

#include "../utils/helper.h"
#include "CodersLib.h"
#include <memory>

unique_ptr<CoderProps> getDefaultCoderProps(uint8_t coder_type, uint8_t coder_level, int coder_param = -1);

unique_ptr<CoderProps> getReadsPositionsCoderProps(uint8_t coder_level, uint8_t lzma_pos_dataperiod_param);

unique_ptr<CoderProps> getVarLenEncodedPgCoderProps(uint8_t coder_level);

unique_ptr<CoderProps> getRelativeOffsetDeltasOfPairsValueCoderProps(uint8_t coder_level);

#endif //PGTOOLS_PARAMSLIBRARY_H
