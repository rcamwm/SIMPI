#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <cstring>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Simpi.h"
using namespace SimpiNS;

void createSynchObject(int processCount);
void runParallelProcesses(char *progname, int processCount);

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        printf("Usage: ./prog2 prog1 <num workers>");
        exit(2);
    }

    char programName[100];
    strcpy(programName, argv[1]);
    int processCount = atoi(argv[2]);

    createSynchObject(processCount);
    runParallelProcesses(programName, processCount);

    exit(0);
}

void createSynchObject(int processCount)
{
    int fd;
    size_t synchObjectSize = sizeof(SynchObject) + sizeof(int) * (processCount + 1);
    fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR | O_CREAT, 0777);
    if (fd == -1) {
        perror("Unable to create synch info: ");
        exit(1);
    }
    ftruncate(fd, synchObjectSize);
    SynchObject* sharedMemory = (SynchObject*)mmap(NULL, synchObjectSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    if (sharedMemory == MAP_FAILED) {
        perror("Unable to mmap shared memory: ");
        exit(1);
    }

    // initialize ready to zero
    for (int i = 0; i <= processCount; i++) {
        sharedMemory->ready[i] = 0;
    }
    sharedMemory->processCount = processCount;
}

void runParallelProcesses(char *programName, int processCount)
{
    for (int i = 0; i < processCount; i++)
    {
        std::string process_count_str = std::to_string(processCount);
        std::string process_id_str = std::to_string(i);

        char* args[] = {programName, 
                        const_cast<char*>(process_id_str.c_str()),
                        const_cast<char*>(process_count_str.c_str()),
                        NULL};

        if (fork() == 0)
        {
            execv(programName, args);
        }
    }
}
