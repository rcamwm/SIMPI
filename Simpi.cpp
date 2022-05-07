#include "Simpi.h"

Simpi::Simpi(int _id, int _processCount)
{
    id = _id;
    processCount = _processCount;

    size_t synchObjectSize = sizeof(SynchObject) + sizeof(int) * (processCount + 1);
    int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
    if (fd == -1)
    {
        perror("Unable to shm_open synchObject: ");
        exit(1);
    }
    shm_fd = fd;
    synchInfo = (SynchObject*)mmap(NULL, synchObjectSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (synchInfo == MAP_FAILED) {
        perror("Unable to mmap synchInfo: ");
        exit(1);
    }
}

Simpi::~Simpi()
{
    for (std::pair<std::string, MatrixMetadata> matrix : matrixInfo) {
        freeMatrix(matrix.second.uniqueID);
    }
    shm_unlink(SYNCH_OBJECT_MEM_NAME);
    close(shm_fd);
    munmap(synchInfo, sizeof(SynchObject) + sizeof(int) * (processCount + 1));
}

void Simpi::synch()
{
    int* ready = synchInfo->ready;
    int synchid = ready[processCount] + 1;
    ready[id] = synchid;
    int breakout = 0;
    while (1)
    {
        breakout = 1;
        for (int i = 0; i < processCount; i++) 
        {
            if (ready[i] < synchid) 
            {
                breakout = 0;
                break;
            }
        }
        if (breakout == 1) 
        {
            ready[processCount] = synchid;
            // and here we increment the additional variable
            break;
        }
    }
}

std::pair<std::string, double*> Simpi::createMatrix(int x, int y)
{
    
    size_t size = x * y * sizeof(double);
    if (id == 0) 
    {
        // generate a uniqueid for matrix
        std::string uniqueID = getSharedMemName();
        
        // create a shared mem object
        int fd = shm_open(uniqueID.c_str(), O_RDWR | O_CREAT, 0777);
        if (fd == -1) 
        {
            std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" + uniqueID + "):  ";
            perror(msg.c_str());
            exit(1);
        }
        ftruncate(fd, size);

        // allocate matrix
        double* matrix = (double*)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (matrix == MAP_FAILED)
        {
            perror("Unable to mmap matrix: ");
            exit(1);
        }

        MatrixMetadata metadata;
        metadata.size = size;
        metadata.fileDescriptor = fd;
        strcpy(metadata.uniqueID, uniqueID.c_str());
        metadata.matrixData = matrix;
        matrixInfo[uniqueID] = metadata;

        // write name to synch_object so that other processes can get the uniqie
        // id
        strcpy(synchInfo->lastMatrixID, uniqueID.c_str());
        synch();
        return std::make_pair(uniqueID, matrix);
    }
    else 
    {
        // wait for id 0 to create the shared memory for matrix
        synch();
        // get the unique id from the synch object
        std::string uniqueID = synchInfo->lastMatrixID;
        // open and allocate the shared memory
        int fd = shm_open(uniqueID.c_str(), O_RDWR, 0777);
        if (fd == -1) 
        {
            std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" + uniqueID + "):  ";
            perror(msg.c_str());
            exit(1);
        }

        double* matrix = (double*)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (matrix == MAP_FAILED) 
        {
            std::string msg = std::to_string(id) + ": Unable to mmap matrix: ";
            perror(msg.c_str());
            exit(1);
        }

        // create a metadata object
        MatrixMetadata metadata;
        strcpy(metadata.uniqueID, uniqueID.c_str());
        metadata.fileDescriptor = fd;
        metadata.matrixData = matrix;
        matrixInfo[uniqueID] = metadata;
        return std::make_pair(uniqueID, matrix);
    }
}

void Simpi::freeMatrix(std::string uniqueID)
{
    // get matrix metadata
    // close fd, shmunmap, munmap
    MatrixMetadata metadata = matrixInfo[uniqueID];
    close(metadata.fileDescriptor);
    shm_unlink(metadata.uniqueID);
    munmap(metadata.matrixData, metadata.size);
}

std::string Simpi::getSharedMemName()
{
    // gets a unique name for each shared memory based on the time that each was
    // made
    timeval curTime;
    gettimeofday(&curTime, NULL);
    unsigned long micro = curTime.tv_sec * (uint64_t)1000000 + curTime.tv_usec;
    std::string name = "simpi_" + std::to_string(micro);
    return name;
}