//
// Created by Eduard Valeyev on 5/11/20.
//

#ifndef SEQUANT_SRC_SEQUANT_VERSION_H_IN_H
#define SEQUANT_SRC_SEQUANT_VERSION_H_IN_H

/* SeQuant version X.Y.Z-id */
#define SEQUANT_VERSION "2.0.0"

/* SeQuant major version */
#define SEQUANT_MAJOR_VERSION 2

/* SeQuant minor version */
#define SEQUANT_MINOR_VERSION 0

/* SeQuant micro version */
#define SEQUANT_MICRO_VERSION 0

/* SeQuant prerelease id */
#define SEQUANT_PRERELEASE_ID "alpha.1"

namespace sequant {

/** \return a string with the Git SHA1 revision hash tag of SeQuant */
const char* git_revision() noexcept;

/**
 * \return a string with the human-readable description of the current source
 *   tree of SeQuant
 * \note see `git describe --dirty` for the format description
 * */
const char* git_description() noexcept;

}  // namespace sequant

#endif  // SEQUANT_SRC_SEQUANT_VERSION_H_IN_H
