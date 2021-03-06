project(clustering)

cmake_minimum_required(VERSION 2.6)

################################################################################
# Shortcuts
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core)

################################################################################
# Define cmake modules directory
################################################################################
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOLS
################################################################################

ADD_SUBDIRECTORY(tools)

################################################################################
# DELIVERY
################################################################################

# If your current login name is different from your GForge login name, you have
# to overwrite the CPACK_USER_NAME to be the same as your GForge login
#SET (CPACK_USER_NAME                    "your_gforge_login")

# We set the version number
#SET (CPACK_PACKAGE_VERSION              "1.0.0")  

# We have to tell what is the server name
SET (CPACK_GFORGE_PROJECT_NAME          "gatb-tools")

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwated files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/build/"  
    "^${GATB_CORE_HOME}/.project"
    "^${GATB_CORE_HOME}/.gitignore"
    "^${GATB_CORE_HOME}/doc/design"
    "^${GATB_CORE_HOME}/DELIVERY.md"
)

include (Delivery)
