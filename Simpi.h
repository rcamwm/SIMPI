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

namespace SimpiNS
{
    #define SYNCH_OBJECT_MEM_NAME "/simpi_shared_mem"
    #define UNIQUE_ID_SIZE 23

    typedef struct MatrixMetadata 
    {
        char uniqueID[UNIQUE_ID_SIZE];
        int fileDescriptor;
        size_t size;
        double* matrixData;
    } MatrixMetadata;

    typedef struct SynchObject 
    {
        int processCount;
        char lastMatrixID[UNIQUE_ID_SIZE];
        int ready[];
    } SynchObject;

    class Simpi 
    {
        public:
            Simpi(int _id, int _processCount);
            ~Simpi();
            int getID() { return id; }
            int getProcessCount() { return processCount; }
            SynchObject* getSynchInfo() { return synchInfo; }

            std::pair<std::string, double*> createMatrix(int rows, int cols);

            void freeMatrix(std::string uniqueID);
            void synch();
            void synchExtraCycles(int cycles);

        private:
            int id;
            int processCount;
            int shm_fd;
            SynchObject* synchInfo;
            std::map<std::string, MatrixMetadata> matrixInfo;
            std::string getSharedMemoryName();

            double *initializeMatrixMemory_LeaderThread(int &fd, const std::string &uniqueID, int size);
            double *initializeMatrixMemory_FollowerThread(int &fd, const std::string &uniqueID, int size);
            void createMatrixMetadata(int size, int fd, std::string uniqueID, double *matrix);

    };
}
#endif