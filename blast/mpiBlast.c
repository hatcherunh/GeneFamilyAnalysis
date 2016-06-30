/*
 * $Id: mpiBlast.c 49 2015-07-09 21:43:24Z pjh $
 *
 * James Jackson
 *
 * This must be compiled with mpicc.
 *
 * pjh Dec 2014: modified to use blast 2.2.30 which uses blastp instead of
 *               blastall.
 * pjh Dec 2015: added support for optionally using blastn.
 *
 * pjh Jun 2015: added support for optionally using blastx and for setting
 *               the output format.
 *
 * pjh Jul 2015: added support for specifying the number of hits to keep.
 *               also replaced scanf with my own readLine function so that
 *               spaces in the FASTA header do not cause a problem, and so
 *               the buffer is not overrun. And major cleanup of buildNewString.
 *               Also changed the commandline interface so that all the args
 *               sent to mpiBlast are simply sent on to the blast command,
 *               and the blast command should be the first argument sent
 *               to mpiBlast. this makes it easier to use this tool for
 *               other purposes than doPairwiseBlasts. in doing this got
 *               rid of the optional number-of-queries argument, which was
 *               not being used anyway.
 */

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>

//#define DEBUG

#define BEGIN_TAG 1
#define MESSAGE_TAG 2
#define END_TAG 3
#define COMPLETE_TAG 99

#define SCHEDULER_PROCESS 0
#define WRITER_PROCESS 1

// pjh Jul 2015 comment added: used for MPI messages
#define BUFFER_SIZE 20000

// pjh Jul 2015 comment added: used for sending data from scheduler to workers
// this is what controls the load balancing
// i.e. how many times blast is invoked by a worker
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 20000
#endif

/*
 * Called upon a fatal error
 *
 */
void fatal(char *message)
{
  fprintf(stderr, "%s\n", message);
  exit(-1);
}

/*
 * This function reads one line from a file. If the buffer is not
 * big enough to contain the full line, then the piece that fits
 * is returned and no more bytes are read from the file. However,
 * if this line is a FASTA header, then the rest of the line is
 * read and discarded.
 *
 * fp - file pointer to read from
 * buf - buffer where chars should be put
 * buflen - length of the buffer
 *
 * The last char read is returned. Note this char might not be placed in
 * the buffer. Actually, this code assumes that a FASTA file is being
 * read, meaning it is a series of lines, and that the newline will be
 * detected when there is nothing in the buffer. Also, note that the newline
 * is not placed in the buffer.
 */
int readLine(FILE *fp, char *buf, int bufLen)
{
  int i = 0;
  int c = getc(fp);
  while (c != '\n' && c != EOF && i < bufLen-2)
  {
    buf[i] = c;
    i += 1;
    c = getc(fp);
  }
  if (i >= bufLen-2)
  {
    if (buf[0] == '>')
    {
      // discard the rest of the header
      int c2 = c;
      while (c2 != '\n')
      {
        c2 = getc(fp);
        if (c2 == EOF)
        {
          fprintf(stderr, "FASTA header incomplete at EOF\n");
          exit(EXIT_FAILURE);
        }
      }
    }
    else
    {
      // need to stash the last char read since it is not a newline
      buf[i] = c;
      i += 1;
    }
  }
  // the bufLen-2 in the above conditions reserves space for this NULL
  // and the last char read when the buffer limit is reached
  buf[i] = 0;
  // sanity check
  if (c == EOF && i > 0)
  {
    fprintf(stderr, "incomplete last line in FASTA file\n");
    exit(EXIT_FAILURE);
  }
  return c;
}

/*
 * This function reads in queries from a file, and 
 * returns a string containing those queries
 *
 * fp - file pointer to read from
 * queryCount - number of queries to read from file
 *
 * pjh July 2015: major changes to simplify
 */
char* buildNewString(FILE *fp, int *amountRead)
{
#ifdef DEBUG
    fprintf(stderr, "buildNewString called\n");
#endif
  const int buflen = 16384;
  char buf[buflen];

  char *finalString;

  int querysRead = 0;

  int eofCheck = 0;

  long bufferOffsetForLastQueryStart = 0;
  long fileOffsetForLastQueryStart = 0;

  int allocSize = BLOCK_SIZE;
  finalString = malloc(sizeof(char) * allocSize);
  if (finalString == NULL) fatal("buildNewString: malloc failed");
  finalString[0] = 0;

  //read in lines from file until EOF
  while(1)
  {
#ifdef DEBUG
    fprintf(stderr, "buildNewString called readLine\n");
#endif

    // first, get current file offset
    long offset = ftell(fp);

    // now read in a line from file
    eofCheck = readLine(fp, buf, buflen);

#ifdef DEBUG
    fprintf(stderr, "%s\n", buf);
    fprintf(stderr, "buildNewString: readLine returned\n");
#endif
    //if at end of file, return partial string
    if (eofCheck == EOF)
    {
#ifdef DEBUG
      fprintf(stderr, "buildNewString: EOF detected\n");
#endif
      if (querysRead > 0)
      {
        *amountRead = querysRead;
        return finalString;
      }
      else
        return NULL;
    }

    // has the end of the block been reached?
    // (the +2 is for a newline and then the null on the end)
    // (this is a loop because may need to repeatedly resize)
    while ((strlen(finalString) + strlen(buf) + 2) > allocSize)
    {
#ifdef DEBUG
      fprintf(stderr, "buildNewString: beyond block size\n");
#endif
      if (querysRead > 1) // count includes an incomplete query in process
      {
        // is there an incomplete query in the block?
        if (buf[0] != '>')
        {
#ifdef DEBUG
          fprintf(stderr, "buildNewString: discard incomplete query\n");
#endif
          //discard the current incomplete query from block
          finalString[bufferOffsetForLastQueryStart] = 0;

          // back up in the file to where that incomplete query began
          fseek(fp, fileOffsetForLastQueryStart, SEEK_SET);
        }
        else
        {
#ifdef DEBUG
          fprintf(stderr, "buildNewString: discard only last line\n");
#endif
          // only need to back up in the file to where last line began
          fseek(fp, offset, SEEK_SET);
        }

        *amountRead = querysRead;
        return finalString;
      }
      else
      {
#ifdef DEBUG
        fprintf(stderr, "buildNewString: reallocing block\n");
#endif
        // current query is too large to fit in a block
        // need to use a bigger buffer for this query
        allocSize *= 2;
        finalString = realloc(finalString, allocSize);
        if (finalString == NULL) fatal("buildNewString: realloc failed");
      }
    }

    //check for the start of a new query
    if (buf[0] == '>')
    {
#ifdef DEBUG
      fprintf(stderr, "buildNewString: start of new query (offset %ld)\n",
        offset);
#endif
      querysRead++;

      fileOffsetForLastQueryStart = offset;
      bufferOffsetForLastQueryStart = strlen(finalString);
    }

#ifdef DEBUG
    fprintf(stderr, "buildNewString: strcat next\n");
#endif
    strcat(finalString, buf);
    strcat(finalString, "\n");
#ifdef DEBUG
    fprintf(stderr, "buildNewString: strcat done\n");
#endif
  }
  fatal("buildNewString: falling off bottom\n");
}

/*
 * The scheduler will read queries from a file and send them to the workers
 * to be processed
 *
 * filename - path to read queries from
 * size     - number of workers
 *
 */
void scheduler(char* filename, int size)
{
#ifdef DEBUG
    fprintf(stderr, "scheduler started\n");
#endif
    //open file, gets filename from command-line args
    FILE *fp;
  
    //number of workers that have been sent a complete message
    int finishedWorkers = 0;

    //used to determine the ID of a sender
    MPI_Status status;

    //toSend contains the entire outgoing message
    char *toSend;
  
    //outgoing messages will be placed in this buffer
    char *buffer;
  
    buffer = malloc(BUFFER_SIZE);
    if (buffer == NULL) fatal("scheduler: malloc of buffer failed");

    int queriesRead;

    fp = fopen(filename, "r");
    if (fp == NULL) fatal("scheduler: fopen failed");
  
    if(!fp)
    {
        fprintf(stderr, "%s, %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
  
#ifdef DEBUG
    fprintf(stderr, "scheduler initialized\n");
#endif
    //loop until all of the workers have completed
    while(finishedWorkers < size - 2)
    {
        //receive ready message from a worker
        MPI_Recv(buffer, 8, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
          MPI_COMM_WORLD, &status);
  
#ifdef DEBUG
    fprintf(stderr, "scheduler got message\n");
#endif
        //get sender of the message
        int sender = status.MPI_SOURCE;

        //build a new query string to send to the worker
        toSend = buildNewString(fp, &queriesRead);
      
        //if no more queries, send complete message
        if(toSend == NULL)
        {
#ifdef DEBUG
            fprintf(stderr, "scheduler sending complete message\n");
#endif
            MPI_Send("", 1, MPI_CHAR, sender, COMPLETE_TAG, MPI_COMM_WORLD);
            finishedWorkers++;
        }
        else 
        { 
            int messageLength;
  
            messageLength = strlen(toSend);
            //fprintf(stderr, "message length is %d\n", messageLength);
            //fprintf(stderr, "%s\n", toSend);
            //fprintf(stderr, "####################\n", toSend);

            int cntSent = 0;

#ifdef DEBUG
            fprintf(stderr, "scheduler sending begin message\n");
#endif
            //send begin message tag
            MPI_Send("", 1, MPI_CHAR, sender, BEGIN_TAG, MPI_COMM_WORLD);
      
            //send message
            while (cntSent < messageLength)
            {
#ifdef DEBUG
                fprintf(stderr, "scheduler copying message into buffer\n");
#endif
                //copy characters from message into buffer
                int n;
                if (messageLength - cntSent >= BUFFER_SIZE)
                {
                  n = BUFFER_SIZE;
                }
                else
                {
                  n = messageLength - cntSent;
                }
                strncpy(buffer, toSend + cntSent, n);
    
#ifdef DEBUG
            fprintf(stderr, "scheduler sending buffer message\n");
#endif
                //send buffer 
                MPI_Send(buffer, n, MPI_CHAR, sender, MESSAGE_TAG,
                  MPI_COMM_WORLD);

                cntSent += n;
            }

#ifdef DEBUG
            fprintf(stderr, "scheduler sending end message\n");
#endif
            //send end message tag
            MPI_Send("", 1, MPI_CHAR, sender, END_TAG, MPI_COMM_WORLD);  
    
            free(toSend);
        }
    }

#ifdef DEBUG
    fprintf(stderr, "scheduler sending complete message 2\n");
#endif
    //send complete message to writer process
    MPI_Send("", 1, MPI_CHAR, WRITER_PROCESS, COMPLETE_TAG, MPI_COMM_WORLD);
}

//the writer process receives messages from the worker processes, and 
//writes those messages to a file.  The writer process is used so that
//all of the output is consolidated into one file. 
//
//The writer process will run until a complete message is sent to it from the 
//scheduler process
void writer(char* filename)
{
#ifdef DEBUG
    fprintf(stderr, "writer started\n");
#endif
    //used to get sender of message
    MPI_Status status;

    //open output file
    FILE *fp;
   
    fp = fopen(filename, "w");
    if (fp == NULL) fatal("writer: fopen failed");

    if(!fp)
    {
        fprintf(stderr, "%s, %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
  
    //received messages are placed in here
    char *buffer; 
  
#ifdef DEBUG
    fprintf(stderr, "writer initialized\n");
#endif
    buffer = malloc(BUFFER_SIZE);
    if (buffer == NULL) fatal("writer: malloc failed");

    while(1)
    {
        //receive a message from any source with any tag
        MPI_Recv(buffer, BUFFER_SIZE, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
          MPI_COMM_WORLD, &status);
    
#ifdef DEBUG
        fprintf(stderr, "writer receives message\n");
#endif
        int sender = status.MPI_SOURCE;
    
        int tag = status.MPI_TAG;

        //if from the scheduler, then we need to exit
        if(sender == SCHEDULER_PROCESS)
        {
            fclose(fp);

            free(buffer);
#ifdef DEBUG
            fprintf(stderr, "writer exiting\n");
#endif
            
            return;
        }
        //if from a worker, receive messages from only that worker until
        //end tag is sent
        else
        {
            while(tag != END_TAG)
            { 
                MPI_Recv(buffer, BUFFER_SIZE, MPI_CHAR, sender, MPI_ANY_TAG,
                  MPI_COMM_WORLD, &status);
      
#ifdef DEBUG
                fprintf(stderr, "writer receives message 2\n");
#endif
                tag = status.MPI_TAG;
  
                //if it is a valid message, print it to file
                if(tag == MESSAGE_TAG) 
                    fprintf(fp, "%s", buffer);
            }
        }
    }  
}

void *workerHelper(void *args)
{
    int tag = BEGIN_TAG;

    int *pipe = (int*)args;

    int errorCheck;

    char buffer[BUFFER_SIZE+1];

    MPI_Status status;

    while(tag != END_TAG)
    {
        //receive message from scheduler
        if((errorCheck = MPI_Recv(buffer, BUFFER_SIZE, MPI_CHAR,
          SCHEDULER_PROCESS, MPI_ANY_TAG, MPI_COMM_WORLD, &status)) !=
          MPI_SUCCESS)
        {
            fprintf(stderr, "MPI Error in helper!\n");
        }
#ifdef DEBUG
        fprintf(stderr, "workerHelper got message\n");
#endif

        int countReceived;
        MPI_Get_count(&status, MPI_CHAR, &countReceived);
#ifdef DEBUG
        fprintf(stderr, "received %d of %d\n", countReceived, BUFFER_SIZE);
        { int i; for (i = 0; i < countReceived; i++) putc(buffer[i], stderr); }
        fprintf(stderr, "##########\n");
#endif

        //check if it is message or end message       
        tag = status.MPI_TAG;
 
        //write message to blast through pipe
        if(tag != END_TAG)
        {
#ifdef DEBUG
        fprintf(stderr, "workerHelper writing to pipe\n");
#endif
            if((errorCheck = write(*pipe, buffer, countReceived)) == -1)
            {
                fprintf(stderr, "Write error in helper!\n");
            }
        }
    }

#ifdef DEBUG
    fprintf(stderr, "workerHelper closing pipe\n");
#endif
    close(*pipe);
#ifdef DEBUG
    fprintf(stderr, "workerHelper returning\n");
#endif

    return NULL;
}

//the worker function 
void worker(int rank, char** blastArgs)
{
#ifdef DEBUG
    fprintf(stderr, "worker %d started\n", rank);
#endif
    MPI_Status status;  
    
    //pipes that will be used
    int toBlastPipe[2];
    int fromBlastPipe[2];
  
    int blocksSearched = 0;

    int errorCheck;
 
    int tag = BEGIN_TAG;
 
    //buffer that received messages will be placed into
    char *buffer;
  
    buffer = malloc(BUFFER_SIZE);
    if (buffer == NULL) fatal("worker: malloc failed");
  
    //send ready message to scheduler
    if((errorCheck = MPI_Send("", 1, MPI_CHAR, SCHEDULER_PROCESS, 0,
      MPI_COMM_WORLD)) != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI error sending ready message to scheduler!\n");
    }

    //receive start message
    if((errorCheck = MPI_Recv(buffer, BUFFER_SIZE, MPI_CHAR, SCHEDULER_PROCESS,
      MPI_ANY_TAG, MPI_COMM_WORLD, &status)) != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI error receiving start message from scheduler!\n");
    }

#ifdef DEBUG
    fprintf(stderr, "worker %d initialized\n", rank);
#endif
    tag = status.MPI_TAG;

    //if(status.MPI_TAG == COMPLETE_TAG)
    //    return;
 
    while(tag != COMPLETE_TAG)
    { 
        //create toBlast pipe
        if((errorCheck = pipe(toBlastPipe)) == -1)
        {
            fprintf(stderr, "Pipe error on pipe()!, errno is %s\n",
              strerror(errno));
            errno = 0;
        }

        //create fromBlast pipe
        if((errorCheck = pipe(fromBlastPipe)) == -1)
        {
            fprintf(stderr, "Pipe error on pipe()!, errno is %s\n",
              strerror(errno));
            errno = 0;
        }

        //create new process
        pid_t pid = fork();
   
        //child
        if(pid == 0)
        {
            //replace stdin with pipe from parent
            if((errorCheck = dup2(toBlastPipe[0], 0)) == -1)
            {
                fprintf(stderr, "Pipe error on dup2()!, errno is %s\n",
                  strerror(errno));
                errno = 0;
            }
    
            //replace stdout with pipe from parent
            if((errorCheck = dup2(fromBlastPipe[1], 1)) == -1)
            {
                fprintf(stderr, "Pipe error on dup2()!, errno is %s\n",
                  strerror(errno));
                errno = 0;
            }    
  
            //close all unused pipes
            close(toBlastPipe[0]);
            close(toBlastPipe[1]);
            close(fromBlastPipe[0]);
            close(fromBlastPipe[1]);
     
#ifdef DEBUG
            fprintf(stderr, "worker %d starting blast:\n", rank);
            {
              int i;
              for (i = 0; blastArgs[i] != NULL; i++)
              {
                fprintf(stderr, "  worker %d: %s\n", rank, blastArgs[i]);
              }
            }
#endif
            // does not return if successful
            execvp(blastArgs[0], blastArgs);
            perror("failure in invoking blast tool");
            fatal("aborting\n");
        }
        else 
        { 
            int fStatus;
            
            pthread_t tid;

            tag = BEGIN_TAG;
  
            //close pipes that are not in use
            close(toBlastPipe[0]);
            close(fromBlastPipe[1]);
      
            errorCheck = pthread_create(&tid, NULL, workerHelper,
              (void*)(&(toBlastPipe[1])));
            
            if(errorCheck != 0)
            {
                fprintf(stderr, "Error creating thread in worker!\n");
            }

            //send begin tag to writer to establish connection
            if((errorCheck = MPI_Send("", 1, MPI_CHAR, WRITER_PROCESS,
              BEGIN_TAG, MPI_COMM_WORLD)) != MPI_SUCCESS)
            {
                fprintf(stderr, "MPI Error sending begin tag to writer\n");
            }
 
            //get number of bytes read from blast
            int bytesRead = 0;

            //read all data from pipe
            bytesRead = read(fromBlastPipe[0], buffer, BUFFER_SIZE); 
      
            while(bytesRead > 0) 
            {
                //null terminate the data from pipe to cut out extra characters
                //previously in buffer
                buffer[bytesRead] = '\0';
  
#ifdef DEBUG
    fprintf(stderr, "worker %d sending data to writer\n", rank);
#endif
                //send data to writer
                if((errorCheck = MPI_Send(buffer, BUFFER_SIZE, MPI_CHAR,
                  WRITER_PROCESS, MESSAGE_TAG, MPI_COMM_WORLD)) != MPI_SUCCESS)
                {
                    fprintf(stderr, "MPI Error sending data to writer!\n");
                }
   
                //read more data from pipe
                bytesRead = read(fromBlastPipe[0], buffer, BUFFER_SIZE);
            }    
    
#ifdef DEBUG
    fprintf(stderr, "worker %d sending end tag to writer\n", rank);
#endif
            //send end tag to writer to close connection
            if((errorCheck = MPI_Send("", 1, MPI_CHAR, WRITER_PROCESS, END_TAG,
              MPI_COMM_WORLD)) != MPI_SUCCESS)
            {
                fprintf(stderr, "Error sending end tag to writer!\n"); 
            }

            blocksSearched++;

            //wait for blast process to terminate
            waitpid(pid, &fStatus, WUNTRACED | WCONTINUED);

            //close blast pipe
            close(fromBlastPipe[0]);    
 
#ifdef DEBUG
    fprintf(stderr, "worker %d sending ready message to scheduler\n", rank);
#endif
            //send ready message to scheduler
            if((errorCheck = MPI_Send("", 1, MPI_CHAR, SCHEDULER_PROCESS, 0,
              MPI_COMM_WORLD)) != MPI_SUCCESS)
            {
                fprintf(stderr, "Error sending ready message to scheduler!\n");
            }
    
            //receive new message from scheduler
            if((errorCheck = MPI_Recv(buffer, BUFFER_SIZE,  MPI_CHAR,
              SCHEDULER_PROCESS, MPI_ANY_TAG, MPI_COMM_WORLD, &status)) !=
              MPI_SUCCESS)
            {
                fprintf(stderr,
                  "Error receiving new message from scheduler!\n");
            }
  
            tag = status.MPI_TAG;
        }
    }
}

/*
 *  The arguments to main are simply passed through as the commandline
 *  to be executed. That is, the first arg to this program should be
 *  the blast command to be executed. There must be -db, -query and
 *  -out arguments, which will actually be stripped out and not sent
 *  on to the blast tool.
 *
 */

void usageMessage(void)
{
  fprintf(stderr,
    "Args: blastCommand -db database -query queryFile -out outputFile "
    "<any other blast args you want>\n");
  exit(1);
}

int main(int argc, char **argv)
{
    int size, rank;

    int threadProvided;

    char *queryFileName = 0;
    char *outFileName = 0;

    // command line to invoke the blast tool
    char **blastArgs;

    // malloc an arg array for the blast command
    // not all the slots will be used however
    blastArgs = malloc(sizeof(char*) * argc);
    if (blastArgs == NULL) fatal("malloc failed in main\n");

    // run through args and pull out the -query and -out args
    int i = 1;
    int j = 0;
    while (i < argc)
    {
      if (!strcmp(argv[i], "-query"))
      {
        queryFileName = argv[i+1];
        i += 2;
      }
      else if (!strcmp(argv[i], "-out"))
      {
        outFileName = argv[i+1];
        i += 2;
      }
      else
      {
        blastArgs[j] = argv[i];
        j += 1;
        i += 1;
      }
    }
    blastArgs[j] = NULL;

    // make sure that -query and -out were all given
    if (queryFileName == 0 || outFileName == 0) usageMessage();

    //initialize MPI
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threadProvided);

    //get number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    //get this process's rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    if(rank == SCHEDULER_PROCESS)
    {
        scheduler(queryFileName, size);
    }
    else if(rank == WRITER_PROCESS)
    {
        writer(outFileName);
    }
    else
        worker(rank, blastArgs);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
