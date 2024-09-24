////////////////////////////////////////////////////////////////////////////
//                           **** DPCM-XQ ****                            //
//                  Xtreme Quality DPCM Encoder/Decoder                   //
//                    Copyright (c) 2024 David Bryant                     //
//                          All Rights Reserved                           //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <fcntl.h>
#endif

typedef uint64_t rms_error_t;     // best if "double" or "uint64_t", "float" okay in a pinch
#define MAX_RMS_ERROR UINT64_MAX
// typedef double rms_error_t;     // best if "double" or "uint64_t", "float" okay in a pinch
// #define MAX_RMS_ERROR DBL_MAX

#define CLIP(data, min, max) \
if ((data) > (max)) data = max; \
else if ((data) < (min)) data = min;

#define DPCM_SDX2   0   // Squareroot-Delta-Exact
#define DPCM_CBD2   1   // Cuberoot-Delta-Exact
#define DPCM_L2XP   2   // Linear-to-Exponential

#define DPCM_NOFUNC 0   // no function
#define DPCM_ENCODE 1   // encode function
#define DPCM_DECODE 2   // decode function

#define STATE_INIT  0   // initial state (no delta allowed)
#define STATE_SYNC  1   // last sample was a sync value (no delta)
#define STATE_DELTA 2   // last sample was a delta value
#define STATE_DELT2 3   // last sample was a double-delta value

static int16_t decode_table [256], *decode_index = decode_table + 128;
static int8_t nearest_0 [65536], *nearest_0_index = nearest_0 + 32768;
static int8_t nearest_1 [65536], *nearest_1_index = nearest_1 + 32768;
static int verbosity, measure_noise, shaping_weight, use_dns = 1;

struct dpcm_state {
    int32_t pcmprev, pcmdata;   // previous & current PCM value
    int32_t weight, error;      // for noise shaping
    int nchans, state;          // num channels & optional state
};

static const char *sign_on = "\n"
" DPCM-XQ  Xtreme Quality Raw DPCM Encoder / Decoder  Version 0.2\n"
" Copyright (c) 2024 David Bryant. All Rights Reserved.\n\n";

static const char *usage =
" Usage:     DPCM [-d|-e] [-options] infile.raw outfile.raw\n\n"
" Operation: Specify '-' for either file spec for piped operation.\n"
"            Mode is 8-bit raw SDX2 stereo unless overridden below.\n\n"
" Options:  -[0-16]= encode lookahead samples (default = 3, max = 16)\n"
"           -c     = override default mode and use CBD2 table instead\n"
"           -d     = decode operation (raw 8-bit DPCM to raw 16-bit PCM)\n"
"           -e     = encode operation (raw 16-bit PCM to raw 8-bit DPCM)\n"
"           -f     = encode flat noise (no noise shaping, aka -s0.0)\n"
"           -h     = display this help message\n"
"           -l     = override default mode and use L2XP table instead\n"
"           -m     = override stereo default and treat files as mono\n"
"           -n     = measure and report quantization noise\n"
"           -q     = quiet mode (display errors only)\n"
"           -s<n>  = enable noise shaping, (-1.0 < n <= 1.0)\n"
"           -v     = verbose (display lots of info)\n"
"           -y     = overwrite outfile if it exists\n\n"
" Web:       Visit www.github.com/dbry/dpcm-xq for latest version and info\n\n";

static uint32_t dpcm_encode_blocks (FILE *source, FILE *destin, int nchans, int block_samples, uint32_t total_samples, int lookahead);
static uint32_t dpcm_decode_blocks (FILE *source, FILE *destin, int nchans);
static double strtod_hexfree (const char *nptr, char **endptr);
static int64_t get_file_size (FILE *file);

int main (int argc, char **argv)
{
    int nchans = 2, table = DPCM_SDX2, function = DPCM_NOFUNC, overwrite = 0, asked_help = 0, lookahead = 3;
    char *infilename = NULL, *outfilename = NULL;
    FILE *infile, *outfile;

    // loop through command-line arguments

    while (--argc) {
#if defined (_WIN32)
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {

                    case '0': case '1': case '2':
                    case '3': case '4': case '5':
                    case '6': case '7': case '8':
                    case '9':
                        lookahead = strtol (*argv, argv, 10);

                        if (lookahead > 16) {
                            fprintf (stderr, "\nlookahead must be 0 to 16!\n");
                            return -1;
                        }

                        --*argv;
                        break;

                    case 'C': case 'c':
                        table = DPCM_CBD2;
                        break;

                    case 'D': case 'd':
                        function = DPCM_DECODE;
                        break;

                    case 'E': case 'e':
                        function = DPCM_ENCODE;
                        break;

                    case 'F': case 'f':
                        shaping_weight = use_dns = 0;
                        break;

                    case 'H': case 'h':
                        asked_help = 1;
                        break;

                    case 'L': case 'l':
                        table = DPCM_L2XP;
                        break;

                    case 'M': case 'm':
                        nchans = 1;
                        break;

                    case 'N': case 'n':
                        measure_noise = 1;
                        break;

                    case 'Q': case 'q':
                        verbosity = -1;
                        break;

                    case 'S': case 's':
                        shaping_weight = strtod_hexfree (++*argv, argv) * 1024.0;

                        if (shaping_weight <= -1024 || shaping_weight > 1024) {
                            fprintf (stderr, "\ninvalid noise shaping value!");
                            return -1;
                        }

                        use_dns = 0;
                        --*argv;
                        break;

                    case 'V': case 'v':
                        verbosity = 1;
                        break;

                    case 'Y': case 'y':
                        overwrite = 1;
                        break;

                    default:
                        fprintf (stderr, "\nillegal option: %c !\n", **argv);
                        return 1;
                }
        else if (!infilename) {
            infilename = malloc (strlen (*argv) + 10);
            strcpy (infilename, *argv);
        }
        else if (!outfilename) {
            outfilename = malloc (strlen (*argv) + 10);
            strcpy (outfilename, *argv);
        }
        else {
            fprintf (stderr, "\nextra unknown argument: %s !\n", *argv);
            return 1;
        }
    }

    if (verbosity >= 0)
        fprintf (stderr, "%s", sign_on);

    if (!outfilename || asked_help) {
        printf ("%s", usage);
        return 0;
    }

    if (function == DPCM_NOFUNC) {
        fprintf (stderr, "no DPCM function specified, use -d or -e!\n");
        return 1;
    }

    if (strcmp (infilename, "-") && strcmp (outfilename, "-") && !strcmp (infilename, outfilename)) {
        fprintf (stderr, "can't overwrite input file (specify different/new output file name)\n");
        return -1;
    }

    if (!overwrite && strcmp (outfilename, "-") && (outfile = fopen (outfilename, "r"))) {
        fclose (outfile);
        fprintf (stderr, "output file \"%s\" exists (use -y to overwrite)\n", outfilename);
        return -1;
    }

    if (strcmp (infilename, "-")) {
        if (!(infile = fopen (infilename, "rb"))) {
            fprintf (stderr, "can't open file \"%s\" for reading!\n", infilename);
            return -1;
        }
    }
    else {
#ifdef _WIN32
        setmode (fileno (stdin), O_BINARY);
#endif
        infile = stdin;
    }

    if (strcmp (outfilename, "-")) {
        if (!(outfile = fopen (outfilename, "wb"))) {
            fprintf (stderr, "can't open file \"%s\" for writing!\n", outfilename);
            fclose (infile);
            return -1;
        }
    }
    else {
#ifdef _WIN32
        setmode (fileno (stdout), O_BINARY);
#endif
        outfile = stdout;
    }

    if (table == DPCM_SDX2)
        for (int dvalue = 0; dvalue <= 127; ++dvalue)
            decode_index [dvalue] = dvalue * abs (dvalue) * 2;
    else if (table == DPCM_CBD2)
        for (int dvalue = 0; dvalue <= 127; ++dvalue)
            decode_index [dvalue] = dvalue * dvalue * dvalue / 64;
    else if (table == DPCM_L2XP) {
        uint32_t svalue = 0, delta;

        for (int dvalue = 0; dvalue <= 127; ++dvalue) {
            decode_index [dvalue] = (svalue + 128) >> 8;
            delta = (svalue * 441 + 3103) / 6206;
            svalue += (delta > 256) ? delta : 256;
        }
    }

    for (int dvalue = 0; dvalue <= 127; ++dvalue)
        decode_index [-dvalue] = -decode_index [dvalue];

    if (function == DPCM_ENCODE)
        dpcm_encode_blocks (infile, outfile, nchans, 4410, (uint32_t)(get_file_size (infile) / 2 / nchans), lookahead);
    else if (function == DPCM_DECODE)
        dpcm_decode_blocks (infile, outfile, nchans);

    fflush (outfile);
    fclose (outfile);
    fclose (infile);

    return 0;
}

static inline int16_t dpcm_decode_sample (struct dpcm_state *pchan, int8_t value)
{
    pchan->pcmprev = pchan->pcmdata;

    if (value & 1) {
        pchan->pcmdata += decode_index [value];
        CLIP (pchan->pcmdata, -32768, 32767);
        pchan->state = STATE_DELTA;
    }
    else {
        pchan->pcmdata = decode_index [value];
        pchan->state = STATE_SYNC;
    }

    return pchan->pcmdata;
}

static int num_trials, num_improvements, sync_switch [9], sync_same [9], delta_switch [9], delta_same [9];

// Apply noise-shaping to the supplied sample value using the shaping weight
// and accumulated error term stored in the dpcm_state structure. Note that
// the error term in the structure is updated, but won't be "correct" until the
// final re-quantized sample value is added to it (and of course we don't know
// that value yet).

static inline int32_t noise_shape (struct dpcm_state *pchan, int32_t sample)
{
    int32_t temp = -((pchan->weight * pchan->error + 512) >> 10);

    if (pchan->weight < 0 && temp) {
        if (temp == pchan->error)
            temp = (temp < 0) ? temp + 1 : temp - 1;

        pchan->error = -sample;
        sample += temp;
    }
    else
        pchan->error = -(sample += temp);

    return sample;
}

static rms_error_t dpcm_min_error (const struct dpcm_state *pchan, int32_t csample, const int16_t *psample, int8_t *best, int depth, rms_error_t max_error)
{
    int32_t sample_delta = csample - pchan->pcmdata, csample2;
    int16_t best_sync_value, best_delta_value = -128;
    int16_t trial_values [7], trial_count = 0;
    struct dpcm_state chan = *pchan;
    rms_error_t min_error;
    int8_t value;

    if (csample < -32768)
        value = best_sync_value = -126;
    else if (csample > 32767)
        value = best_sync_value = 126;
    else
        value = best_sync_value = nearest_0_index [csample];

    if (best_sync_value > -126)
        trial_values [trial_count++] = best_sync_value - 2;

    if (best_sync_value < 126)
        trial_values [trial_count++] = best_sync_value + 2;

    if (pchan->state != STATE_INIT) {
        if (csample < -32768)
            best_delta_value = -127;
        else if (csample > 32767)
            best_delta_value = 127;
        else
            best_delta_value = nearest_1_index [sample_delta];

        if (best_delta_value > -127) {
            trial_values [trial_count++] = best_delta_value - 2;

            if (best_delta_value > -125)
                trial_values [trial_count++] = best_delta_value - 4;
        }

        if (best_delta_value < 127) {
            trial_values [trial_count++] = best_delta_value + 2;

            if (best_delta_value < 125)
                trial_values [trial_count++] = best_delta_value + 4;
        }

        if (abs (decode_index [best_delta_value] - sample_delta) < abs (decode_index [value] - csample)) {
            trial_values [trial_count++] = best_sync_value;
            value = best_delta_value;
        }
        else
            trial_values [trial_count++] = best_delta_value;
    }

    if (best) *best = value;
    min_error = dpcm_decode_sample (&chan, value) - csample;
    min_error = min_error * min_error;

    // if we're at a leaf, or we're not at a leaf but have already exceeded the error limit, return
    if (!depth || min_error > max_error)
        return min_error;

    // otherwise we execute that naively closest value and search deeper for improvement

    if (chan.weight | chan.error) {
        chan.error += chan.pcmdata;
        csample2 = noise_shape (&chan, psample [chan.nchans]);
    }
    else
        csample2 = psample [chan.nchans];

    min_error += dpcm_min_error (&chan, csample2, psample + chan.nchans, NULL, depth - 1, max_error - min_error);

    for (int tindex = 0; tindex < trial_count; ++tindex) {
        rms_error_t error, threshold;

        chan = *pchan;
        error = dpcm_decode_sample (&chan, trial_values [tindex]) - csample;
        error = error * error;
        threshold = max_error < min_error ? max_error : min_error;

        if (error < threshold) {
            if (chan.weight | chan.error) {
                chan.error += chan.pcmdata;
                csample2 = noise_shape (&chan, psample [chan.nchans]);
            }
            else
                csample2 = psample [chan.nchans];

            error += dpcm_min_error (&chan, csample2, psample + chan.nchans, NULL, depth - 1, threshold - error);

            if (error < min_error) {
                if (best) *best = trial_values [tindex];
                min_error = error;
            }
        }
    }

    if (best && verbosity > 0) {
        ++num_trials;
        if (*best != value) {
            ++num_improvements;

            if (*best & 1) {        // used delta
                if (value & 1)      // original was delta
                    delta_same [abs (*best - best_delta_value)]++;
                else                // original was sync
                    delta_switch [abs (*best - best_delta_value)]++;
            }
            else {                  // used sync
                if (value & 1)      // original was delta
                    sync_switch [abs (*best - best_sync_value)]++;
                else                // original was sync
                    sync_same [abs (*best - best_sync_value)]++;
            }
        }
    }

    return min_error;
}

void generate_dns_values (const int16_t *samples, int sample_count, int num_chans,
    int16_t *values, int16_t min_value, int16_t last_value);

static uint32_t dpcm_encode_blocks (FILE *source, FILE *destin, int nchans, int block_samples, uint32_t total_samples, int lookahead)
{
    uint32_t num_samples = 0, progress_divider = 0, percent;
    int16_t pcm_values [nchans * block_samples];
    int8_t dpcm_values [nchans * block_samples];
    int16_t shaping_values [block_samples], last = 0;
    struct dpcm_state chans [nchans];

    uint32_t noise_windows = 0, noise_samples = 0;
    double rms_noise_average [2] = { 0.0, 0.0 };
    double rms_noise_total [2] = { 0.0, 0.0 };
    double rms_noise_peak [2] = { 0.0, 0.0 };
    uint32_t max_error [2] = { 0, 0 };

    if (verbosity >= 0 && total_samples > 1000) {
        progress_divider = (total_samples + 50) / 100;
        fprintf (stderr, "\rprogress: %d%% ", percent = 0);
        fflush (stderr);
    }

    memset (chans, 0, sizeof (chans));

    for (int ch = 0; ch < nchans; ++ch)
        chans [ch].nchans = nchans;

    for (int pvalue = -32768; pvalue <= 32767; ++pvalue) {
        int min_error = INT_MAX;

        for (int dvalue = -127; dvalue <= 127; ++dvalue)
            if (dvalue & 1) {
                int error = abs (pvalue - decode_index [dvalue]);
                if (error < min_error) {
                    nearest_1_index [pvalue] = dvalue;
                    min_error = error;
                }
            }

        min_error = INT_MAX;

        for (int dvalue = -127; dvalue <= 127; ++dvalue)
            if (!(dvalue & 1)) {
                int error = abs (pvalue - decode_index [dvalue]);
                if (error < min_error) {
                    nearest_0_index [pvalue] = dvalue;
                    min_error = error;
                }
            }
    }

    while (1) {
        int samples_read = fread (pcm_values, sizeof (int16_t) * nchans, block_samples, source);
        double rms_noise [2] = { 0.0, 0.0 };

        if (!samples_read)
            break;

        if (use_dns) {
            generate_dns_values (pcm_values, samples_read, nchans, shaping_values, -256, last);
            last = shaping_values [samples_read - 1];
        }

        for (int ch = 0; ch < nchans; ++ch)
            for (int samp = 0; samp < samples_read; samp++) {
                int8_t *dpcm_index = dpcm_values + samp * nchans + ch;
                int16_t *pcm_index = pcm_values + samp * nchans + ch;
                int32_t csample = *pcm_index, error;
                int depth = lookahead;

                chans [ch].weight = use_dns ? shaping_values [samp] : shaping_weight;

                if (chans [ch].weight | chans [ch].error)
                    csample = noise_shape (chans + ch, csample);

                if (depth > samples_read - samp - 1)
                    depth = samples_read - samp - 1;

                dpcm_min_error (chans + ch, csample, pcm_index, dpcm_index, depth, MAX_RMS_ERROR);
                error = fabs (*pcm_index - dpcm_decode_sample (chans + ch, *dpcm_index));

                if (error > max_error [ch])
                    max_error [ch] = error;

                rms_noise [ch] += (double) error * error;

                if (chans [ch].weight | chans [ch].error)
                    chans [ch].error += chans [ch].pcmdata;
            }

        fwrite (dpcm_values, sizeof (int8_t) * nchans, samples_read, destin);

        for (int ch = 0; ch < nchans; ++ch) {
            rms_noise_average [ch] += sqrt (rms_noise [ch] / samples_read);
            rms_noise_total [ch] += rms_noise [ch];

            if (rms_noise [ch] / samples_read > rms_noise_peak [ch])
                rms_noise_peak [ch] = rms_noise [ch] / samples_read;

            // rms_noise [ch] = 0.0;
        }

        num_samples += samples_read;
        noise_samples += samples_read;
        noise_windows++;

        if (progress_divider) {
            int new_percent = 100 - (total_samples - num_samples) / progress_divider;

            if (new_percent != percent) {
                fprintf (stderr, "\rprogress: %d%% ", percent = new_percent);
                fflush (stderr);
            }
        }
    }

    if (verbosity >= 0)
        fprintf (stderr, "\rsuccessfully encoded %u %s samples\n", num_samples, nchans == 2 ? "stereo" : "mono");

    if (measure_noise && noise_samples) {
        double full_scale_rms = 32768.0 * 32767.0 * 0.5;

        if (nchans == 2) {
            rms_noise_average [0] /= noise_windows * sqrt (full_scale_rms);
            rms_noise_average [1] /= noise_windows * sqrt (full_scale_rms);
            rms_noise_total [0] /= noise_samples * full_scale_rms;
            rms_noise_total [1] /= noise_samples * full_scale_rms;
            rms_noise_peak [0] /= full_scale_rms;
            rms_noise_peak [1] /= full_scale_rms;

            fprintf (stderr, "\n          Channel:    left      right \n");
            fprintf (stderr, "----------------------------------------\n");
            fprintf (stderr, " Max Sample Error:  %6lu     %6lu\n", (unsigned long) max_error [0], (unsigned long) max_error [1]);
            fprintf (stderr, "RMS Average Noise:  %6.2f dB  %6.2f dB\n", log10 (rms_noise_average [0]) * 20.0, log10 (rms_noise_average [1]) * 20.0);
            fprintf (stderr, "  RMS Total Noise:  %6.2f dB  %6.2f dB\n", log10 (rms_noise_total [0]) * 10.0, log10 (rms_noise_total [1]) * 10.0);
            fprintf (stderr, "   RMS Peak Noise:  %6.2f dB  %6.2f dB\n\n", log10 (rms_noise_peak [0]) * 10.0, log10 (rms_noise_peak [1]) * 10.0);
        }
        else {
            rms_noise_average [0] /= noise_windows * sqrt (full_scale_rms);
            rms_noise_total [0] /= noise_samples * full_scale_rms;
            rms_noise_peak [0] /= full_scale_rms;

            fprintf (stderr, "\n Max Sample Error:  %6lu\n", (unsigned long) max_error [0]);
            fprintf (stderr, "RMS Average Noise:  %6.2f dB\n", log10 (rms_noise_average [0]) * 20.0);
            fprintf (stderr, "  RMS Total Noise:  %6.2f dB\n", log10 (rms_noise_total [0]) * 10.0);
            fprintf (stderr, "   RMS Peak Noise:  %6.2f dB\n\n", log10 (rms_noise_peak [0]) * 10.0);
        }
    }

    if (num_trials && verbosity > 0) {
        fprintf (stderr, "num trials = %d, num improvements = %d (%.2f%%)\n", num_trials, num_improvements, num_improvements * 100.0 / num_trials);

        while (1) {
            int max_hits = 0;

            for (int i = 0; i < 9; ++i) {
                if (sync_switch [i] > max_hits) max_hits = sync_switch [i];
                if (sync_same [i] > max_hits) max_hits = sync_same [i];
                if (delta_switch [i] > max_hits) max_hits = delta_switch [i];
                if (delta_same [i] > max_hits) max_hits = delta_same [i];
            }

            if (!max_hits) break;

            for (int i = 0; i < 9; ++i) {
                if (sync_switch [i] == max_hits)  { fprintf (stderr, "delta --> sync  [%2d]  %.2f%%\n", i, sync_switch  [i] * 100.0 / num_improvements); sync_switch [i]  = 0; }
                if (sync_same [i] == max_hits)    { fprintf (stderr, "sync  --> sync  [%2d]  %.2f%%\n", i, sync_same    [i] * 100.0 / num_improvements); sync_same [i]    = 0; }
                if (delta_switch [i] == max_hits) { fprintf (stderr, "sync  --> delta [%2d]  %.2f%%\n", i, delta_switch [i] * 100.0 / num_improvements); delta_switch [i] = 0; }
                if (delta_same [i] == max_hits)   { fprintf (stderr, "delta --> delta [%2d]  %.2f%%\n", i, delta_same   [i] * 100.0 / num_improvements); delta_same [i]   = 0; }
            }
        }
    }

    return num_samples;
}

static uint32_t dpcm_decode_blocks (FILE *source, FILE *destin, int nchans)
{
    int sync_samples [nchans], delta_samples [nchans];
    int delta_run [nchans], longest_delta_run [nchans];
    struct dpcm_state chans [nchans];
    int16_t pcm_values [nchans];
    int8_t dpcm_values [nchans];
    uint32_t num_samples = 0;

    memset (chans, 0, sizeof (chans));

    for (int ch = 0; ch < nchans; ++ch)
        chans [ch].nchans = nchans;

    for (int ch = 0; ch < nchans; ++ch)
        delta_run [ch] = longest_delta_run [ch] = sync_samples [ch] = delta_samples [ch] = 0;

    while (fread (dpcm_values, sizeof (dpcm_values), 1, source)) {
        for (int ch = 0; ch < nchans; ++ch) {
            pcm_values [ch] = dpcm_decode_sample (chans + ch, dpcm_values [ch]);

            if (dpcm_values [ch] == -128) { fprintf (stderr, "decode a minus 128!\n"); return 1; }
            if (dpcm_values [ch] & 1) {
                delta_samples [ch]++;
                delta_run [ch]++;
            }
            else {
                if (delta_run [ch] > longest_delta_run [ch]) longest_delta_run [ch] = delta_run [ch];
                sync_samples [ch]++;
                delta_run [ch] = 0;
            }
        }

        fwrite (pcm_values, sizeof (pcm_values), 1, destin);
        num_samples++;
    }

    if (verbosity > 0) {
        if (nchans == 2) {
            fprintf (stderr, "num samples = %u, syncs: %d left, %d right, deltas: %d left, %d right\n",
                num_samples, sync_samples [0], sync_samples [1], delta_samples [0], delta_samples [1]);
            fprintf (stderr, "longest delta runs = %d left, %d right\n",
                longest_delta_run [0], longest_delta_run [1]); 
        }
        else {
            fprintf (stderr, "num samples = %u, syncs = %d, deltas = %d\n", num_samples, sync_samples [0], delta_samples [0]);
            fprintf (stderr, "longest delta run = %d\n", longest_delta_run [0]); 
        }
    }
    else if (verbosity >= 0)
        fprintf (stderr, "successfully decoded %u %s samples\n", num_samples, nchans == 2 ? "stereo" : "mono");

    return num_samples;
}

// The C-standard function strtod() also handles hex numbers prefixed
// with [+-]0[xX]. Unfortunately this causes problems for us in rare
// cases where a value of zero is specified for one option followed
// by the 'x' option (e.g., -s0xe). This version of strtod() does not
// allow hex specification, but otherwise should be identical.

static double strtod_hexfree (const char *nptr, char **endptr)
{
    const char *sptr = nptr;

    // skip past any leading whitespace and possibly a sign
    while (isspace (*sptr)) sptr++;
    if (*sptr == '+' || *sptr == '-') sptr++;

    // if hex detected ("0x" or "0X"), return 0.0 and end at the X
    if (*sptr == '0' && tolower (sptr [1]) == 'x') {
        if (endptr) *endptr = (char *) sptr + 1;
        return 0.0;
    }

    // otherwise unmodified strtod() result
    return strtod (nptr, endptr);
}

#ifdef _WIN32

static int64_t get_file_size (FILE *hFile)
{
    LARGE_INTEGER Size;
    HANDLE        fHandle;

    if (hFile == NULL)
        return 0;

    fHandle = (HANDLE)_get_osfhandle(_fileno(hFile));

    if (fHandle == INVALID_HANDLE_VALUE || GetFileType(fHandle) != FILE_TYPE_DISK)
        return 0;

    Size.u.LowPart = GetFileSize(fHandle, (DWORD *) &Size.u.HighPart);

    if (Size.u.LowPart == INVALID_FILE_SIZE && GetLastError() != NO_ERROR)
        return 0;

    return (int64_t)Size.QuadPart;
}

#else

static int64_t get_file_size (FILE *hFile)
{
    struct stat statbuf;

    if (!hFile || fstat (fileno (hFile), &statbuf) || !S_ISREG(statbuf.st_mode))
        return 0;

    return (int64_t) statbuf.st_size;
}

#endif


