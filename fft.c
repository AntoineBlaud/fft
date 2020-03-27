#include <stdlib.h>
#include <stdio.h>
#include <sndfile.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

//compile  /usr/bin/gcc-7 -g /home/darkloner99/code/Projet-trans/fft.c -std=c99 -lm -lsndfile -o /home/darkloner99/code/Projet-trans/fft


#define SWAP(a, b) \
    ctmp = (a);    \
    (a) = (b);     \
    (b) = ctmp;
#define ARRAY_LEN(x) ((int)(sizeof(x) / sizeof(x[0])))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define M_PI 3.1415926535

/**
 * DOC SNDFILE http://www.mega-nerd.com/libsndfile/api.html
 * SNDFILE is a special FILE in C
 * Read file: sf_open_read
 * Writefile : sf_open_write
 **/

/** A pointer to an SF_INFO struct  is passed to s f_open_read and s f_open_write
    which fill this struct with information about the file .
**/

SNDFILE *infile, *outfile;
SF_INFO sfinfo_in;  // informations du fichier d'entrée
SF_INFO sfinfo_out; // informations du fichier de sortie
const int N = 1024;
const int G = 128;
const int H = 20;
complex *TW;

typedef struct // structure de retour pour la fonction sfx_mix_mono_read_double
{
    double *data;
    sf_count_t frames_readed; // frames_readed = nombre de valeurs lues
} SFX;
SFX sfx;

double log2(double x)
{
    return log10(x) / log10(2);
}
complex double_to_complex(double d)
{
    return d + I * 0.0;
}
complex *double_array_to_complex_array(double *data, int size, complex *datac)
{
    for (int i = 0; i < size; i++)
        datac[i] = double_to_complex(data[i]);
    
    return datac;
}
double compute_module(complex c)
{
    return sqrt(pow(creal(c), 2) + pow(cimag(c), 2));
}
double *complex_array_to_module(complex *data, int size, double *datam)
{
    for (int i = 0; i < size; i++)
        datam[i] = compute_module(data[i]);
    
    return datam;
}

void print_arr(double *arr, int size)
{
    for (int i = 0; i < size; i++)
        printf("%lf\n", *(arr + i));
}
void print_arrc(complex *arr, int size)
{
    for (int i = 0; i < size; i++)
        printf("%f + I%f\n", creal(*(arr + i)), cimag(*(arr + i)));
}

int msleep(unsigned int tms)
{
    return usleep(tms * 1000);
}
// La fonction ci-dessous permet de récupérer le fichier d'entrée
SNDFILE *open_input_file(char *name)
{
    // tente de récupérer le fichier audio d'entrée
    infile = sf_open(name, SFM_READ, &sfinfo_in);
    // si ==NULL alors problème d'ouverture du fichier
    if (infile == NULL)
    {
        printf("Not able to open input file %s\n", name);
        sf_perror(NULL);
        return NULL;
    }
    return infile;
}

// La fonction ci-dessous permet de créer le fichier de sortie
SNDFILE *open_output_file(char *name)
{
    outfile = sf_open(name, SFM_WRITE, &sfinfo_out); //ouvre le fichier
    if (outfile == NULL)
    {
        printf("Not able to open ouput file %s\n", name);
        sf_perror(NULL);
        return NULL;
    }
    return outfile;
}

// La fonction ci-dessous permet de lire les N prochaines valeurs audio (sous forme de double) contenues
// dans le fichier d'entrée
SFX sfx_mix_mono_read_double(SNDFILE *file, int N)
{
    int k, ch, frames_read;
    int datalen = sfinfo_in.channels * N; // datalen = nombre de channels* N (N=1024)

    // Si le fichier content un seul cannal audio
    if (sfinfo_in.channels == 1)
    {
        // lit les datalen valeurs et les stockent dans data
        sfx.frames_readed = sf_readf_double(file, sfx.data, N);
        return sfx;
    }
    // Si le fichier contient plusieurs cannaux audio
    else
    {
        sfx.frames_readed = sf_readf_double(file, sfx.data, N);
        for (k = 0; k < sfx.frames_readed; k++)
        {
            double mix = 0.0;
            for (ch = 0; ch < sfinfo_in.channels; ch++)
            {
                mix += sfx.data[k * sfinfo_in.channels + ch];
            }
            sfx.data[k] = mix / sfinfo_in.channels;
        }
        return sfx;
    }
}

void twiddle(complex *TW, unsigned int size)
{
    complex phi = cexp(-2 * I * M_PI / size);
    TW[0] = 1;
    for (int i = 1; i < size; i++)
        TW[i] = TW[i - 1] * phi;
}

int bitrev(int inp, int numbits)
{
    int rev = 0;
    for (int i = 0; i < numbits; i++)
    {
        rev = (rev << 1) | (inp & 1);
        inp >>= 1;
    }
    return rev;
}

void fftrec(complex *data, complex *result, int size, int log2n)
{
    if (size < 2)
    {
        result[0] = data[0];
        return;
    }
    complex ypair[size], yimpair[size], Fimpair[size], Fpair[size];
    int n, k, N2;
    N2 = size / 2;
    for (n = 0; n < N2; n++)
    {
        ypair[n] = data[n] + data[n + N2];
        yimpair[n] = (data[n] - data[n + N2]) * cexp(-2 * I * M_PI * n / size);
    }
    fftrec(ypair, Fpair, N2, log2n);
    fftrec(yimpair, Fimpair, N2, log2n);
    for (n = 0; n < N2; n++)
    {
        result[2 * n] = Fpair[n];
        result[2 * n + 1] = Fimpair[n];
    }
}

double *create_spectrum(double *data, int size, double *spectrum)
{
    double average = 0.0;
    int cursor = 0;
    for (int i = 0; i < size; i++)
    {
        // ici on fait une moyenne sur G valeurs
        average += data[i];
        if (i % (size / G) == 0 && i > 0)
        {
            spectrum[cursor] = average;
            cursor += 1;
            average = 0;
        }
    }
    return spectrum;
}

void trace_spectrum(double *spectrum, int samplerate)
{
    double max = 0, min = 4200000000, average = 0;
    int screen[G];
    for (int i = 1; i < G - 1; i++)
    {
        average += spectrum[i];
        if (spectrum[i] > max)
            max = spectrum[i];
        if (spectrum[i] < min)
            min = spectrum[i];
    }
    if (max < 6)
        max = 6;
    double stepH = max / H;
    // compute bar
    for (int i = 0; i < G; i++)
        screen[i] = H - spectrum[i] / stepH;
    // print bar
    for (int j = 0; j < H; j++)
    {
        for (int i = 0; i < G; i++)
            (j > screen[i]) ? printf("#") : printf(" ");
        printf("\n");
    }
    // print bottom limit
    for (int i = 0; i < G; i++)
        printf("-");
    printf("\n0");
    // Print frequency band
    for (int i = 0; i < G - 5; i++)
        printf(" ");
    printf("%d\n", samplerate / 2);
}
void process_data(double *data, int size, int samplerate)
{
    complex datac[size], datacOut[size];
    double datam[size], spectrum[G];
    double_array_to_complex_array(data, size, datac);
    fftrec(datac, datacOut, size, log2(N));
    complex_array_to_module(datacOut, size, datam);
    // nysquist limit
    size = size / 2;
    create_spectrum(datam, size, spectrum);
    trace_spectrum(spectrum, samplerate);
    printf("\e[1;1H\e[2J");
}

void browse_audio(SNDFILE *file_in, SNDFILE *file_out)
{
    int samplerate = sfinfo_in.samplerate;
    sf_count_t channels = sfinfo_in.channels;
    int size = N * sfinfo_in.channels;
    TW = malloc(sizeof(complex) * size);
    twiddle(TW, size);
    sfx.data = malloc(sizeof(double) * size);
    clock_t previous, diff;
    double msec;
    int waiting = N * 1000 / samplerate;
    do
    {
        previous = clock();
        sfx = sfx_mix_mono_read_double(file_in, N);
        process_data(sfx.data, sfx.frames_readed, samplerate);
        diff = clock() - previous;
        // on suppose que msec >0 car sinon le programme fonctionnerai mal
        msec = waiting - (diff * 1000 / CLOCKS_PER_SEC);
        msleep(msec);
    } while (sfx.frames_readed == N);
}

int main(int argc, char *argv[])
{
    char *input_file_name = argv[1];
    SNDFILE *input_file = open_input_file(input_file_name);
    if (input_file == NULL)
        return EXIT_FAILURE;
    browse_audio(input_file, NULL);
    return EXIT_SUCCESS;
}