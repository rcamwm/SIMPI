#ifndef SIMPI_H
#define SIMPI_H

#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#define SYNCH_OBJECT_MEM_NAME "/simpi_shared_mem"
#define UNIQUE_ID_SIZE 23

typedef struct MatrixMetadata {
    char unique_id[UNIQUE_ID_SIZE];
    int file_descriptor;
    size_t size;
    double* matrix_data;
} MatrixMetadata;

typedef struct SynchObject {
    int processCount;
    char last_matrix_id[UNIQUE_ID_SIZE];
    int ready[];
} SynchObject;

class Simpi {
    public:
        Simpi(int _id, int _processCount);
        ~Simpi();
        int getID() { return id; }
        int getProcessCount() { return processCount; }
        SynchObject* get_synch_info() { return synch_info; }

        std::pair<std::string, double*> create_matrix(int x, int y);

        void free_matrix(std::string unique_id);
        void synch();

    private:
        int id;
        int processCount;
        int shm_fd;
        SynchObject* synch_info;
        std::map<std::string, MatrixMetadata> matrix_info;
        std::string sync_shared_mem_name;
        std::string get_shared_mem_name();
};

#endif