#include "Simpi.h"

Simpi::Simpi(int _id, int _processCount)
{
    id = _id;
    processCount = _processCount;

    size_t synchObjectSize = sizeof(SynchObject) + sizeof(int) * (processCount + 1);
    int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
    if (fd == -1)
    {
        perror("Unable to shm_open synch_object: ");
        exit(1);
    }
    shm_fd = fd;
    synch_info = (SynchObject*)mmap(NULL, synchObjectSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (synch_info == MAP_FAILED) {
        perror("Unable to mmap synch_info: ");
        exit(1);
    }
}

Simpi::~Simpi()
{
  for (std::pair<std::string, MatrixMetadata> matrix : matrix_info) {
    free_matrix(matrix.second.unique_id);
  }
  shm_unlink(SYNCH_OBJECT_MEM_NAME);
  close(shm_fd);
}

void Simpi::synch()
{
    int* ready = synch_info->ready;
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

std::pair<std::string, double*> Simpi::create_matrix(int x, int y)
{
    
  size_t size = x * y * sizeof(double);
  if (id == 0) {
    // generate a uniqueid for matrix
    std::string unique_id = get_shared_mem_name();
    
    // create a shared mem object
    int fd = shm_open(unique_id.c_str(), O_RDWR | O_CREAT, 0777);
    if (fd == -1) {
      std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" +
                        unique_id + "):  ";
      perror(msg.c_str());
      exit(1);
    }
    ftruncate(fd, size);

    // allocate matrix
    double* matrix =
        (double*)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (matrix == MAP_FAILED) {
      perror("Unable to mmap matrix: ");
      exit(1);
    }

    MatrixMetadata metadata;
    metadata.size = size;
    metadata.file_descriptor = fd;
    strcpy(metadata.unique_id, unique_id.c_str());
    metadata.matrix_data = matrix;
    matrix_info[unique_id] = metadata;

    // write name to synch_object so that other processes can get the uniqie
    // id
    strcpy(synch_info->last_matrix_id, unique_id.c_str());
    synch();
    std::pair<std::string, double*> pass_back;
    pass_back = std::make_pair(unique_id, matrix);
    return pass_back;
  }
  else {
    // wait for id 0 to create the shared memory for matrix
    synch();
    // get the unique id from the synch object
    std::string unique_id = synch_info->last_matrix_id;
    // open and allocate the shared memory
    int fd = shm_open(unique_id.c_str(), O_RDWR, 0777);
    if (fd == -1) {
      std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" +
                        unique_id + "):  ";
      perror(msg.c_str());
      exit(1);
    }
    double* matrix =
        (double*)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (matrix == MAP_FAILED) {
      std::string msg = std::to_string(id) + ": Unable to mmap matrix: ";
      perror(msg.c_str());
      exit(1);
    }

    // create a metadata object
    MatrixMetadata metadata;
    strcpy(metadata.unique_id, unique_id.c_str());
    metadata.file_descriptor = fd;
    metadata.matrix_data = matrix;
    matrix_info[unique_id] = metadata;
    std::pair<std::string, double*> pass_back =
        std::make_pair(unique_id, matrix);
    return pass_back;
  }
}

void Simpi::free_matrix(std::string unique_id)
{
  // get matrix metadata
  // close fd, shmunmap, munmap
  MatrixMetadata metadata = matrix_info[unique_id];
  close(metadata.file_descriptor);
  shm_unlink(metadata.unique_id);
  munmap(metadata.matrix_data, metadata.size);
}

std::string Simpi::get_shared_mem_name()
{
  // gets a unique name for each shared memory based on the time that each was
  // made
  timeval curTime;
  gettimeofday(&curTime, NULL);
  unsigned long micro = curTime.tv_sec * (uint64_t)1000000 + curTime.tv_usec;
  std::string name = "simpi_" + std::to_string(micro);
  return name;
}