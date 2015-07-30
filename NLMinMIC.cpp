/*
 * NLM
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>
#include <gdal/gdal_priv.h>
#include <iostream>
#include <iomanip>
#include<fstream>
#include<omp.h>
#include <cmath>
#include <time.h>
#include <cassert>
#include <gdal/cpl_string.h>
#pragma comment(lib, "gdal_i.lib")

using namespace std;

#define LUTMAX 30.0
#define LUTMAXM1 29.0
#define LUTPRECISION 1000.0
#define dTiny 1e-10
#define fTiny 0.00000001f
#define fLarge 100000000.0f
#define M_PI 3.14159265358979323846

// Read image
bool readImageGDAL(unsigned char **pImageData, int &width, int &height, int &nChannels, const char *filePath, double trans[6])
{
    GDALAllRegister();
    GDALDataset *poDataset = NULL;
    poDataset = (GDALDataset *) GDALOpen(filePath, GA_ReadOnly);
    if(poDataset == NULL)
    {
        GDALClose(poDataset);
        return false;
    }
    width = poDataset->GetRasterXSize();
    height = poDataset->GetRasterYSize();

    printf("width=%d\n", width);
    printf("height=%d\n", height);

    CPLErr aaa = poDataset->GetGeoTransform(trans);
    int k = 0;

    GDALRasterBand *pBand;
    int i = 0;
    int nRastercount = poDataset->GetRasterCount();

    // one channel, the gray image
    if (nRastercount == 1)
    {
        nChannels = 1;
        pBand = poDataset->GetRasterBand(1);
        *pImageData = new unsigned char[width * height];
        pBand->RasterIO(GF_Read,
                        0, 0,
                        width, height,
                        *pImageData,
                        width, height,
                        GDT_Byte,
                        0,
                        0);
        GDALClose(poDataset);
        return true;
    }
    // three channels, and the output is RGB image
    else if ( nRastercount == 3 && (nChannels == 3 || nChannels < 0) )
    {
        nChannels = 3;
        *pImageData = new unsigned char[nRastercount * width * height];
        for (i = 1; i <= nRastercount; ++ i)
        {
            //Band GDAL RGB is stored in order, we usually need to be converted to BGR storage, namely low address to high address: B G R
            unsigned char *pImageOffset = *pImageData + i - 1;
            GDALRasterBand *pBand = poDataset->GetRasterBand(nRastercount - i + 1);
            pBand->RasterIO(
                GF_Read,
                0, 0,
                width, height,
                pImageOffset,
                width, height,
                GDT_Byte,
                3,
                0);
        }
        GDALClose(poDataset);
        return true;
    }
    //3 channels, but the required output grayscale images
    else if ( nRastercount == 3 && nChannels == 1 )
    {
        unsigned char **img = new unsigned char *[nRastercount];
        for (i = 0; i < nRastercount; i++)
        {
            img[i] = new unsigned char[width * height];
        }
        for (i = 1; i <= nRastercount; ++ i)
        {
            //Band GDAL RGB is stored in order, we usually need to be converted to BGR storage, namely low address to high address: B G R
            pBand = poDataset->GetRasterBand(nRastercount - i + 1);
            pBand->RasterIO(GF_Read,
                            0, 0,
                            width, height,
                            img[i - 1],
                            width, height,
                            GDT_Byte,
                            0,
                            0);
        }
        GDALClose(poDataset);
        *pImageData = new unsigned char[width * height];
        for (int r = 0; r < height; ++ r)
        {
            for (int c = 0; c < width; ++ c)
            {
                int t = (r * width + c);
                //r g b components are accounted for in turn:0.299 0.587 0.144,can be simplified as 3:6:1
                //img[1.2.3]Correspond to BGR
                (*pImageData)[t] = (img[2][t] * 3 + img[1][t] * 6 + img[0][t] + 5) / 10;
            }
        }

        for (i = 0; i < nRastercount; ++ i)
        {
            delete [] img[i];
        }
        delete []img;
        img = NULL;
        return true;
    }
    else
    {
        return false;
    }
}




char *strlwr(char *s)
{
    char *p;
    for(p = s; *p != '\0'; p++)
    {
        if('A' <= (*p) && (*p) <= 'Z')
            (*p) += 32;
    }
    return s;
}


char *findImageTypeGDAL( char *pDstImgFileName)
{
    char *dstExtension = strlwr(strrchr(pDstImgFileName, '.') + 1);
    //char *dstExtension = strrchr(pDstImgFileName,'.') ;
    char *Gtype = NULL;
    if		(0 == strcmp(dstExtension, "bmp")) Gtype = "BMP";
    else if (0 == strcmp(dstExtension, "jpg")) Gtype = "JPEG";
    else if (0 == strcmp(dstExtension, "png")) Gtype = "PNG";
    else if (0 == strcmp(dstExtension, "tif")) Gtype = "GTiff";
    else if (0 == strcmp(dstExtension, "gif")) Gtype = "GIF";
    else Gtype = NULL;
    return Gtype;
}


// write image
bool WriteImageGDAL( char *pDstImgFileName, unsigned char *pImageData, int width, int height, int nChannels, double trans[6])
{
    GDALAllRegister();
    char *GType = NULL;
    GType = findImageTypeGDAL(pDstImgFileName);
    if (GType == NULL)
    {
        return false;
    }
    GDALDriver *pMemDriver = NULL;
    pMemDriver = GetGDALDriverManager()->GetDriverByName("MEM");
    if( pMemDriver == NULL )
    {
        return false;
    }
    GDALDataset *pMemDataSet = pMemDriver->Create("", width, height, nChannels, GDT_Byte, NULL);
    GDALRasterBand *pBand = NULL;
    int nLineCount = width * nChannels;
    unsigned char *ptr1 = (unsigned char *)pImageData;
    for (int i = 1; i <= nChannels; i++)
    {
        pBand = pMemDataSet->GetRasterBand(nChannels - i + 1);
        pBand->RasterIO(GF_Write,
                        0,
                        0,
                        width,
                        height,
                        ptr1 + i - 1 ,
                        width,
                        height,
                        GDT_Byte,
                        nChannels,
                        nLineCount);
    }
    //Write the generated data set to the target file
    GDALDriver *pDstDriver = NULL;
    pDstDriver = (GDALDriver *)GDALGetDriverByName(GType);
    if (pDstDriver == NULL)
    {
        return false;
    }
    //Write to the geographic reference information
    pMemDataSet->SetGeoTransform( trans );
    GDALDataset *poDstDS;
    poDstDS = pDstDriver->CreateCopy(pDstImgFileName, pMemDataSet, FALSE, NULL, NULL, NULL);
    if( poDstDS != NULL )
        delete poDstDS;
    GDALClose(pMemDataSet);
    return true;
}

__attribute__((target(mic))) void memsub(float *ptI, float fValue, int iLength)
{
    for (int ii = 0; ii < iLength; ii++) ptI[ii] = fValue;
}

// LUT tables
void  wxFillExpLut(float *lut, int size)
{
    for (int i = 0; i < size; i++) lut[i] = expf( - (float) i / LUTPRECISION);
}

__attribute__((target(mic))) float exp (float dif, float *lut)
{

    if (dif >= (float) LUTMAXM1) return 0.0;
    int  x = (int) floor( (double) dif * (float) LUTPRECISION);
    float y1 = lut[x];
    float y2 = lut[x + 1];
    return y1 + (y2 - y1) * (dif * LUTPRECISION -  x);
}

__attribute__((target(mic))) float pixdiff (float *u0, float *u1, int i0, int j0, int i1, int j1, int radius, int width0, int width1)
{

    float dist = 0.0;
    for (int s = -radius; s <= radius; s++)
    {
        int l = (j0 + s) * width0 + (i0 - radius);
        float *ptr0 = &u0[l];
        l = (j1 + s) * width1 + (i1 - radius);
        float *ptr1 = &u1[l];
        for (int r = -radius; r <= radius; r++, ptr0++, ptr1++)
        {
            float dif = (*ptr0 - *ptr1);
            dist += (dif * dif);
        }
    }
    return dist;
}

__attribute__((target(mic)))  float pixdiff (float **u0, float **u1, int i0, int j0, int i1, int j1, int radius, int channels, int width0, int width1)
{

    float dif = 0.0f;
    for (int ii = 0; ii < channels; ii++)
    {
        dif += pixdiff (u0[ii], u1[ii], i0, j0, i1, j1, radius, width0, width1);
    }
    return dif;
}


void nlmeans_ipol(int iDWin,            // Half size of patch
                  int iDBloc,           // Half size of research window
                  float fSigma,         // Noise parameter
                  float fFiltPar,       // Filtering parameter
                  float *ptI,          // *****************Input
                  float *ptO,          // ***************Output
                  int iChannels, int iWidth, int iHeight)
{
    long it = 0;
    long ita = 0;
    long itb = 0;
    long itc = 0;
    timeb it3, it4, it5, it6, it7, it8, it9, it10;
    ftime(&it7);

    // length of each channel   ap=w * h
    int ap = iWidth * iHeight;
    //  length of comparison window   wl=window length
    int ihwl = (2 * iDWin + 1);
    int iwl = (2 * iDWin + 1) * (2 * iDWin + 1);
    int icwl = iChannels * iwl;
    // filtering parameter
    float fSigma2 = fSigma * fSigma;
    float fH = fFiltPar * fSigma;
    float fH2 = fH * fH;
    // multiply by size of patch, since distances are not normalized
    fH2 *= (float) icwl;
    // tabulate exp(-x), faster than using directly function expf
    int lul = (int) ((float) LUTMAX * (float) LUTPRECISION);
    float *ptLu = new float[lul];
    wxFillExpLut(ptLu, lul);
    // auxiliary variable
    // number of denoised values per pixel
    float *ptSum = new float[ap];
    memsub(ptSum, 0.0f, ap);
    // clear output
    //  ****************for (int ii=0; ii < iChannels; ii++) memsub(ptO[ii], 0.0f, ap);
    memsub(ptO, 0.0f, ap);
    ftime(&it8);
    itb = (it8.time - it7.time) * 1000 + (it8.millitm - it7.millitm);
    printf("nlm head time : %d ms\n", itb);


    long bian = iWidth * iHeight;
    long halfbian = bian / 2;
    long bd4 = bian / 4;
    ftime(&it3);
    long cell;

#pragma offload target(mic:1) inout(ptSum:length(ap)),in(ptI:length(ap)),in(ptLu:length(lul)),inout(ptO:length(ap) alloc_if(1) free_if(1))
    {
        #pragma omp parallel  for  shared(ptI, ptO)  schedule(dynamic) // nowai
        for (long  cell = 0; cell < bd4; cell++)
        {
            {
                float *fpODenoised = new float[iwl];
                int x = cell * 4 % iWidth;
                int  y = cell * 4 / iHeight;

                for(int ai = 1; ai < 5; ai++)
                {
                    int iDWin0 = MIN(iDWin, MIN(iWidth - 1 - x, MIN(iHeight - 1 - y, MIN(x, y))));

                    int imin = MAX(x - iDBloc, iDWin0);
                    int jmin = MAX(y - iDBloc, iDWin0);

                    int imax = MIN(x + iDBloc, iWidth - 1 - iDWin0);
                    int jmax = MIN(y + iDBloc, iHeight - 1 - iDWin0);

                    //  clear current denoised patch
                    memsub(fpODenoised, 0.0f, iwl);

                    // maximum of weights. Used for reference patch
                    float fMaxWeight = 0.0f;

                    // sum of weights
                    float fTotalWeight = 0.0f;

                    for (int j = jmin; j <= jmax; j++)
                        for (int i = imin ; i <= imax; i++)
                            if (i != x || j != y)
                            {
                                float fDif = pixdiff (ptI, ptI, x, y, i, j, iDWin0, iWidth, iWidth);

                                // dif^2 - 2 * fSigma^2 * N      dif is not normalized
                                fDif = MAX(fDif - 2.0f * (float) icwl *  fSigma2, 0.0f);
                                fDif = fDif / fH2;

                                float fWeight = exp (fDif, ptLu);

                                if (fWeight > fMaxWeight) fMaxWeight = fWeight;

                                fTotalWeight += fWeight;

                                for (int is = -iDWin0; is <= iDWin0; is++)
                                {
                                    int aiindex = (iDWin + is) * ihwl + iDWin;
                                    int ail = (j + is) * iWidth + i;

                                    for (int ir = -iDWin0; ir <= iDWin0; ir++)
                                    {
                                        int iindex = aiindex + ir;
                                        int il = ail + ir;

                                        fpODenoised[iindex] += fWeight * ptI[il];
                                    }

                                }

                            }

                    // current patch with fMaxWeight
                    for (int is = -iDWin0; is <= iDWin0; is++)
                    {
                        int aiindex = (iDWin + is) * ihwl + iDWin;
                        int ail = (y + is) * iWidth + x;

                        for (int ir = -iDWin0; ir <= iDWin0; ir++)
                        {
                            int iindex = aiindex + ir;
                            int il = ail + ir;
                            fpODenoised[iindex] += fMaxWeight * ptI[il];
                        }
                    }

                    fTotalWeight += fMaxWeight;

                    if (fTotalWeight > fTiny)
                    {
                        for (int is = -iDWin0; is <= iDWin0; is++)
                        {
                            int aiindex = (iDWin + is) * ihwl + iDWin;
                            int ail = (y + is) * iWidth + x;

                            for (int ir = -iDWin0; ir <= iDWin0; ir++)
                            {
                                int iindex = aiindex + ir;
                                int il = ail + ir;

                                ptSum[il]++;

                                ptO[il] += fpODenoised[iindex] / fTotalWeight;
                            }
                        }
                    }
                    x++;
                }
                delete [] fpODenoised;
            }
        }
    }
    
    //#pragma offload target(mic:0)out(ptO[0:5120]:alloc_if(0) free_if(0))
    ftime(&it4);
    it = (it4.time - it3.time) * 1000 + (it4.millitm - it3.millitm);
    printf("loop time : %d ms\n", it);

    ftime(&it5);

    for (int ii = 0; ii < ap; ii++)
    {
        if (ptSum[ii] > 0.0)
        {
            ptO[ii] /= ptSum[ii];

        }
        else
        {
            ptO[ii] = ptI[ii];

        }
    }
    ftime(&it6);
    ita = (it6.time - it5.time) * 1000 + (it6.millitm - it5.millitm);
    printf("ptSum loop time : %d ms\n", ita);
    delete[] ptLu;
    delete[] ptSum;
}

int main(int argc, char *argv[])
{
    timeb t1, t2;
    long  t;
    ftime(&t1);
    unsigned char *hehe = NULL;
    unsigned char **pImageData = &hehe;
    int width = 0;
    int	height = 0;
    int nChannels = 0;
    double trans[6];

    readImageGDAL(pImageData, width, height, nChannels, "gray10.png", trans);

    long whc = width * height * nChannels;
    long wh = width * height;

    float *ptI = NULL;
    float *ptO = NULL;

    float *noisy = new float[whc];
    float *denoised = new float[whc];

    long k = 0;
    long k1 = 0;
    for( k1 = 0; k1 < nChannels; k1++)
    {
        k = 0;
        while (k < wh)
        {
            noisy[k + k1 * wh] = (*pImageData)[k1 + k * nChannels];
            denoised[k + k1 * wh] = (*pImageData)[k1 + k * nChannels];
            k++;
        }
    }
    ptI = noisy;
    ptO = denoised;

    timeb t3, t4;
    ftime(&t3);

    nlmeans_ipol(1, 10, 10.0f, 0.4f, ptI,  ptO, nChannels, width, height);

    ftime(&t4);
    t = (t4.time - t3.time) * 1000 + (t4.millitm - t3.millitm);
    printf("nlm time:%ld ms\n", t);


    // Conversion RGB storage mode
    for( k1 = 0; k1 < nChannels; k1++)
    {
        k = 0;
        while (k < wh)
        {
            (*pImageData)[k1 + k * nChannels] = (unsigned char)denoised[k + k1 * wh];
            k++;
        }
    }
    WriteImageGDAL("graysave.png", *pImageData, width, height, nChannels, trans);
    delete [](*pImageData);
    delete []noisy;
    delete []denoised;
    ftime(&t2);
    t = (t2.time - t1.time) * 1000 + (t2.millitm - t1.millitm);
    return 0;
}
