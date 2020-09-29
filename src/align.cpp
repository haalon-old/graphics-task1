#include "align.h"
#include <cmath>
#include <string>
#include <utility>

using std::string;
using std::cout;
using std::endl;

using std::tuple;
using std::get;
using std::tie;
using std::make_tuple;



bool mirrorFlag=false;

double getMetric(int x_shift, int y_shift, Matrix<int> ch1, Matrix<int> ch2 )
{
    double metric = 0;
    for(int x =ch1.n_cols*0.1; x < ch1.n_cols*0.9 ; x++ )
    {
       for(int y =ch1.n_rows*0.1; y < ch1.n_rows*0.9  ; y++ )
        {
            if ( (x - x_shift) > 0 && (x - x_shift) < ch1.n_cols &&  (y - y_shift) > 0 && (y - y_shift) < ch1.n_rows )
            {
                metric+= ( ch1(y,x)-ch2(y-y_shift,x-x_shift) ) * ( ch1(y,x)-ch2(y-y_shift,x-x_shift) );
                //
            }
            //else cout<<"A"<<x<<"," << y << " ";
            
        }
    }
    return metric / ( (ch1.n_rows - abs(y_shift)) * (ch1.n_cols - abs(x_shift)) );
}

double getMetric2(int x_shift, int y_shift, Matrix<int> ch1, Matrix<int> ch2 )
{
    double metric = 0;
    for(int x =15; x < ch1.n_cols -15 ; x++ )
    {
       for(int y = 0; y < ch1.n_rows  ; y++ )
        {
            if ( (x - x_shift) > 0 && (x - x_shift) < ch1.n_cols &&  (y - y_shift) > 0 && (y - y_shift) < ch1.n_rows )
            {
                metric+= ch2(y,x)  *  ch1(y,x);
            }
            
        }
    }
    return metric;
}

tuple<int,int,int,int> findShift(Image srcImage)
{
    Matrix<int> R(srcImage.n_rows/3, srcImage.n_cols);
    Matrix<int> G(srcImage.n_rows/3, srcImage.n_cols);
    Matrix<int> B(srcImage.n_rows/3, srcImage.n_cols);
    int  n = srcImage.n_rows/3, m = n *2;

    
    

    // Split into RGB channels
    for(int x =0; x < srcImage.n_cols; x++ )
    {
       for(int y =0; y < srcImage.n_rows/3; y++ )
        {
            R(y,x) =   get<0>(  srcImage(y+m,x) ) ;
            G(y,x) =  get<0>(  srcImage(y+n,x) ) ;
            B(y,x) = get<0>(  srcImage(y,x) );
        }
    }

    double minR = getMetric(0,0,R,G);
    int Rx=0,Ry=0;
    double minB = getMetric(0,0,R,B);
    int Bx=0,By=0;
    double temp;

    if( fmin( R.n_cols, R.n_rows) < 400 )
    {
         //looking for the best shift
    
        for(int x =-10; x < 10; x++ )
        {
            for(int y =-10; y < 10; y++ )
            {
                temp = getMetric(x,y,G,R);
                if(temp<minR) { minR=temp; Rx=x;Ry=y;}
                temp = getMetric(x,y,G,B);
                if(temp<minB) { minB=temp; Bx=x;By=y;}
            }
        }

        //cout << " " << Rx << " " << Ry << " " << Bx << " " << By << " " << getMetric(Rx,Ry,G,R) << " " << getMetric(Bx,By,G,B)<<"\n";
        return make_tuple(Rx,Ry,Bx,By);
    }
    else
    {   

        tie(Rx,Ry,Bx,By) = findShift( resize(srcImage,1.0/2) );

        int TRx=0, TRy=0, TBx=0, TBy=0;
        

        Rx*=2;
        Ry*=2;
        Bx*=2;
        By*=2;
        for(int x =-2; x < 2; x++ )
        {
            for(int y =-2; y < 2; y++ )
            {
                temp = getMetric(x+Rx,y+Ry,G,R);
                if(temp<minR) { minR=temp; TRx=x+Rx;TRy=y+Ry;}
                temp = getMetric(x+Bx,y+By,G,B);
                if(temp<minB) { minB=temp; TBx=x+Bx;TBy=y+By;}
            }
        }
        
        //cout << " " << TRx << " " << TRy << " " << TBx << " " << TBy << " " << getMetric(TRx,TRy,G,R) << " " << getMetric(TBx,TBy,G,B)<<"\n";
        return make_tuple(TRx,TRy,TBx,TBy);
    }

    

}


Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{

    if(isSubpixel) srcImage = resize(srcImage, subScale );
    Matrix<int> R(srcImage.n_rows/3, srcImage.n_cols);
    Matrix<int> G(srcImage.n_rows/3, srcImage.n_cols);
    Matrix<int> B(srcImage.n_rows/3, srcImage.n_cols);
    int  n = srcImage.n_rows/3, m = n *2;

    // Split into RGB channels
    for(int x =0; x < srcImage.n_cols; x++ )
    {
       for(int y =0; y < srcImage.n_rows/3; y++ )
        {
            R(y,x) =   get<0>(  srcImage(y+m,x) ) ;
            G(y,x) =  get<0>(  srcImage(y+n,x) ) ;
            B(y,x) = get<0>(  srcImage(y,x) );
        }
    }

    //looking for the best shift  
    int Rx=0,Ry=0;   
    int Bx=0,By=0;  
    
    tie(Rx,Ry,Bx,By) = findShift(srcImage); 
    
    Image tempImage(srcImage.n_rows/3, srcImage.n_cols);
    for(int x =0; x < srcImage.n_cols; x++ )
    {
       for(int y =0; y < srcImage.n_rows/3; y++ )
        {
            get<1>(tempImage(y,x)) = G(y,x);

            if(  (x - Rx) > 0 && (x - Rx) < R.n_cols && (y - Ry) > 0 && (y - Ry) < n)
            {
                get<0>(tempImage(y,x)) = R(y-Ry,x-Rx);
            }
            
            if(  (x - Bx) > 0 && (x - Bx) < R.n_cols && (y - By) > 0 && (y - By) < n)
            {
                get<2>(tempImage(y,x)) = B(y-By,x-Bx);
            }

        }
    }

    if(isSubpixel) srcImage = resize(srcImage, 1.0/subScale );
    mirrorFlag = isMirror;
    if(isPostprocessing)   
    {
        if(postprocessingType== "--gray-world") tempImage = gray_world(tempImage);
        if(postprocessingType== "--unsharp") tempImage = unsharp(tempImage);
        if(postprocessingType== "--autocontrast") tempImage = autocontrast(tempImage,fraction);
    }


     return tempImage;
    /*
    Image tempImage(srcImage.n_rows/3, srcImage.n_cols);
    int  n = srcImage.n_rows/3, m = ( srcImage.n_rows/3 ) *2;
    for(int x =0; x < srcImage.n_cols; x++ )
    {
       for(int y =0; y < srcImage.n_rows/3; y++ )
        {
            tempImage(y,x) = make_tuple(  get<2>(  srcImage(y+m,x) )*0.11 + get<1>(  srcImage(y+m,x) )*0.59 + get<0>(  srcImage(y+m,x) )*0.3 ,  get<2>(  srcImage(y+n,x) )*0.11 + get<1>(  srcImage(y+n,x) )*0.59 + get<0>(  srcImage(y+n,x) )*0.3 ,  get<2>(  srcImage(y,x) )*0.11 + get<1>(  srcImage(y,x) )*0.59 + get<0>(  srcImage(y,x) )*0.3);
        }
    }

    return tempImage;*/

}

void fixOverflow(Image &src_image)
{
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            get<0>( src_image(y,x) )= fmax( 0 , fmin(255, get<0>( src_image(y,x) ) ) );
            get<1>( src_image(y,x) )= fmax( 0 , fmin(255, get<1>( src_image(y,x) ) ) );
            get<2>( src_image(y,x) )= fmax( 0 , fmin(255, get<2>( src_image(y,x) ) ) );
        }
    }
}

Image mirror(Image src_image, int rows, int cols  )
{
    Image tempImage(src_image.n_rows + rows -1, src_image.n_cols +cols-1);
    int Y=0, X=0;
   
    for(int x =0, kx=-cols/2; x < tempImage.n_cols; x++,kx++ )
    {
        for(int y =0, ky=-rows/2; y < tempImage.n_rows; y++,ky++ )
        {
            Y = (ky < 0) ? abs(ky)-1 : ( (ky >= src_image.n_rows) ? src_image.n_rows - abs( src_image.n_rows - ky) -1 : ky  ) ;
            X = (kx < 0) ? abs(kx)-1 : ( (kx >= src_image.n_cols) ? src_image.n_cols - abs( src_image.n_cols - kx) -1 : kx  ) ;

            tempImage(y,x) = src_image(Y,X);           
        }
    }

    return tempImage;

}

Image grayScale(Image src_image) {

    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            get<0>( src_image(y,x) )= 0.2125 * get<0>( src_image(y,x)) +  0.7154* get<1>( src_image(y,x)) + 0.0721*get<2>( src_image(y,x));
            get<1>( src_image(y,x) )=get<0>( src_image(y,x) );
            get<2>( src_image(y,x) )=get<0>( src_image(y,x) );
        }  
    }
   
    return src_image;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};

    return custom(src_image, kernel);
}

Image sobel(Image src_image) {
    Image Sx(src_image.n_rows , src_image.n_cols);
    Image Sy(src_image.n_rows , src_image.n_cols);

    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
           get<0>( Sx(y,x) )= 0.2125 * get<0>( Sx(y,x)) +  0.7154* get<1>( Sx(y,x)) + 0.0721*get<2>( Sx(y,x));
           get<0>( Sy(y,x) )= 0.2125 * get<0>( Sy(y,x)) +  0.7154* get<1>( Sy(y,x)) + 0.0721*get<2>( Sy(y,x));
      
        }
    }

    Sx=sobel_x(src_image);
    Sy=sobel_y(src_image);
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
           get<0>( Sx(y,x) )= sqrt( get<0>( Sx(y,x) )*get<0>( Sx(y,x) )+ get<0>( Sy(y,x) )*get<0>( Sy(y,x) ) );
           get<1>( Sx(y,x) )= get<0>( Sx(y,x) );
           get<2>( Sx(y,x) )= get<0>( Sx(y,x) );
        }
    }
    fixOverflow(Sx);
    
    return Sx;
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
   Matrix<double> kernel = {{-1/6.0, -2/3.0, -1/6.0},
                             {-2/3.0, 4.0+1/3.0, -2/3.0},
                             {-1/6.0, -2/3.0, -1/6.0}};
    
    fixOverflow(src_image);
    src_image = custom(src_image,kernel);
    //src_image = mirror(src_image,150,150);
    return src_image;
}

Image gray_world(Image src_image) {

    double r=0,g=0,b=0,sum;

    for(int x =0; x < src_image.n_cols; x++ )
    {
       for(int y =0; y < src_image.n_rows; y++ )
        {
            r+= get<0>(  src_image(y,x) );
            g+= get<1>(  src_image(y,x) );
            b+= get<2>(  src_image(y,x) );

        }
    }
    sum=(r+g+b)/3;

    for(int x =0; x < src_image.n_cols; x++ )
    {
       for(int y =0; y < src_image.n_rows; y++ )
        {
            get<0>(  src_image(y,x) )*= sum/r;
            get<1>(  src_image(y,x) )*= sum/g;
            get<2>(  src_image(y,x) )*= sum/b;

        }
    }
    fixOverflow(src_image);
    return src_image;
}

/*Image oldResize(Image src_image, double scale) {

    if(scale>=1) return src_image;
    Image tempImage(ceil(src_image.n_rows*scale), ceil( src_image.n_cols*scale) );
    Matrix<int> hitcount(ceil(src_image.n_rows*scale), ceil( src_image.n_cols*scale) );
    for(int x =0; x < tempImage.n_cols; x++ )
    {
        for(int y =0; y < tempImage.n_rows; y++ )
        {
            tempImage(y,x) = make_tuple(0,0,0);
            hitcount(y,x) = 0;
        }
    }
    int X=0,Y=0;
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            X=x*scale; Y = y*scale;
            get<0>(tempImage(Y,X)) += get<0>(src_image(y,x));
            get<1>(tempImage(Y,X)) += get<1>(src_image(y,x));
            get<2>(tempImage(Y,X)) += get<2>(src_image(y,x));
            hitcount(Y,X)++;
        }
    }

    for(int x =0; x < tempImage.n_cols; x++ )
    {
        for(int y =0; y < tempImage.n_rows; y++ )
        {
            get<0>(tempImage(y,x)) /= hitcount(y,x);
            get<1>(tempImage(y,x)) /= hitcount(y,x);
            get<2>(tempImage(y,x)) /= hitcount(y,x);
        }
    }

    return tempImage;
}
*/
Image resize(Image src_image, double scale) {

    Image tempImage(src_image.n_rows*scale,  src_image.n_cols*scale );
    
    double X=0,Y=0;   
    double X1, X2, Y1, Y2, div;
    int Px,Py;
    for(int x =0; x < tempImage.n_cols; x++ )
    {
        for(int y =0; y < tempImage.n_rows; y++ )
        {
            X= x/(tempImage.n_cols-1.0);
            Y= y/(tempImage.n_rows-1.0);

            Px=X*(src_image.n_cols-1);
            Py=Y*(src_image.n_rows-1);

            X1=Px/(src_image.n_cols-1.0);
            Y1=Py/(src_image.n_rows-1.0);

            X2=(Px+1)/(src_image.n_cols-1.0);
            Y2=(Py+1)/(src_image.n_rows-1.0);

            div = (X2-X1)*(Y2-Y1);

            get<0>(tempImage(y,x)) = get<0>(src_image(Py,Px))*1.0*( (X2-X)*(Y2-Y)/div ) + get<0>(src_image( fmin(Py+1,src_image.n_rows-1) , Px))*1.0*( (X2-X)*(Y-Y1)/div ) 
            + get<0>(src_image( Py , fmin(Px+1,src_image.n_cols-1)))*1.0*( (X-X1)*(Y2-Y)/div ) + get<0>(src_image( fmin(Py+1,src_image.n_rows-1) , fmin(Px+1,src_image.n_cols-1)))*1.0*( (X-X1)*(Y-Y1)/div ); 

            get<1>(tempImage(y,x)) = get<1>(src_image(Py,Px))*1.0*( (X2-X)*(Y2-Y)/div ) + get<1>(src_image( fmin(Py+1,src_image.n_rows-1) , Px))*1.0*( (X2-X)*(Y-Y1)/div ) 
            + get<1>(src_image( Py , fmin(Px+1,src_image.n_cols-1)))*1.0*( (X-X1)*(Y2-Y)/div ) + get<1>(src_image( fmin(Py+1,src_image.n_rows-1) , fmin(Px+1,src_image.n_cols-1)))*1.0*( (X-X1)*(Y-Y1)/div );

            get<2>(tempImage(y,x)) = get<2>(src_image(Py,Px))*1.0*( (X2-X)*(Y2-Y)/div ) + get<2>(src_image( fmin(Py+1,src_image.n_rows-1) , Px))*1.0*( (X2-X)*(Y-Y1)/div ) 
            + get<2>(src_image( Py , fmin(Px+1,src_image.n_cols-1)))*1.0*( (X-X1)*(Y2-Y)/div ) + get<2>(src_image( fmin(Py+1,src_image.n_rows-1) , fmin(Px+1,src_image.n_cols-1)))*1.0*( (X-X1)*(Y-Y1)/div );

            //if(get<0>(tempImage(y,x)) != get<0>(src_image(y,x))) { cout << "(" << x << "," << y << " " << Px << "," << Py << " " << Py  << "?" << Y*1.0*(src_image.n_rows-1); return src_image;}


        }
    }

    return tempImage;
}



void addKernelElement(tuple<int,int,int> &core , tuple<int,int,int> &nbr, tuple<double> kernel) 
{
    get<0>( core )+= get<0>(nbr)*get<0>(kernel);
    get<1>( core )+= get<1>(nbr)*get<0>(kernel);
    get<2>( core )+= get<2>(nbr)*get<0>(kernel);
}



Image custom(Image src_image, Matrix<double> kernel) {

    if(mirrorFlag) src_image = mirror(src_image, kernel.n_rows, kernel.n_cols);

    Image tempImage(src_image.n_rows, src_image.n_cols);
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            tempImage(y,x) = make_tuple(0,0,0);
        }
    }

    for(int x = kernel.n_cols/2; x < src_image.n_cols-kernel.n_cols/2; x++ )
    {
       for(int y = kernel.n_rows/2; y < src_image.n_rows - kernel.n_rows/2; y++ )
        {
            for (int dx = -kernel.n_cols/2, kx=0; kx < kernel.n_cols; kx++, dx++)
            {
                for (int dy = -kernel.n_rows/2, ky=0; ky < kernel.n_rows; ky++, dy++)
                {                    
                    addKernelElement(tempImage(y,x),src_image(y+dy,x+dx),kernel(ky,kx));
                }
            }
        }
    }    
    
    
    if(mirrorFlag) tempImage = tempImage.submatrix(kernel.n_rows/2,kernel.n_cols/2,tempImage.n_rows -kernel.n_rows+1,tempImage.n_cols -kernel.n_cols+1 );
    return tempImage;
}

Image autocontrast(Image src_image, double fraction) {

    int shift = src_image.n_cols*src_image.n_rows*fraction;
    int Y;
    int H[256];
    for (int i = 0; i < 256; i++)
    {
        H[i]=0;
    }

    for(int x =0; x < src_image.n_cols; x++ )
    {
       for(int y =0; y < src_image.n_rows; y++ )
        {
            
            Y=  0.2125 * get<0>( src_image(y,x)) +  0.7154* get<1>( src_image(y,x)) + 0.0721*get<2>( src_image(y,x));
            H[Y]++;
        }
    }


    Y=0;
    int Ymax=0, Ymin=0;
    for (int i = 0; i < 256; i++)
    {
       if( (Y+= H[i]) >= shift ) { Ymin = i; break;} 
    }

    Y=0;
    for (int i = 255; i >=0; i--)
    {
       if( (Y+= H[i]) >= shift ) { Ymax = i; break;} 
    }

    double coeff = 255.0/(Ymax-Ymin);  
    cout << " " << Ymin << " " << Ymax << " " << coeff ;
    for(int x =0; x < src_image.n_cols; x++ )
    {
       for(int y =0; y < src_image.n_rows; y++ )
        {
            get<0>( src_image(y,x) ) = ( get<0>( src_image(y,x) ) - Ymin)* coeff;
            get<1>( src_image(y,x) ) = ( get<1>( src_image(y,x) ) - Ymin)* coeff;
            get<2>( src_image(y,x) ) = ( get<2>( src_image(y,x) ) - Ymin)* coeff;
            
        }
    }
    fixOverflow(src_image);
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

int GetY(tuple<int,int,int> pixel)
{
    return 0.3 * get<0>( pixel) +  0.59* get<1>(pixel) + 0.11*get<2>( pixel);
}

void sort (int * arr, int size)
{
    int temp;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size-1-i; ++j)
        {
            if( arr[j] > arr[j+1] ) 
            {
                temp = arr[j] ; // swap(arr[j], arr[j+1]);
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }

        }
    }

}

Image median(Image src_image, int radius) {

    src_image = mirror(src_image, 2* radius+1, 2* radius+1);
    
    int size = (2*radius+1)*(2*radius+1);
    int *R= new int[ size ];
    int *G= new int[ size ];
    int *B= new int[ size ];

    Image tempImage(src_image.n_rows, src_image.n_cols);
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            tempImage(y,x) = make_tuple(0,0,0);
        }
    }

    for(int x = radius; x < src_image.n_cols-radius; x++ )
    {
       for(int y = radius; y < src_image.n_rows - radius; y++ )
        {
            for (int dx = -radius, kx=0; dx <= radius; kx++, dx++)
            {
                for (int dy = -radius, ky=0; dy <= radius; ky++, dy++)
                {                    
                    R[kx + (2*radius+1)*ky]= get<0>(src_image(y+dy,x+dx));
                    G[kx + (2*radius+1)*ky]= get<1>(src_image(y+dy,x+dx));
                    B[kx + (2*radius+1)*ky]= get<2>(src_image(y+dy,x+dx));
                }
            }
           
            sort(R, size );  
            sort(G, size );  
            sort(B, size );            
            tempImage(y,x) = make_tuple(R[size/2] , G[size/2], B[size/2] ) ;
           
        }
    }    
    tempImage = tempImage.submatrix(radius,radius,tempImage.n_rows -2*radius,tempImage.n_cols -2*radius );
    delete [] R;
    delete [] G;
    delete [] B;
    return tempImage;
}

Image median_linear(Image src_image, int radius) { //they ignore border processing problem in the original algorythm?
    //

    src_image = mirror(src_image, 2* radius+1, 2* radius+1);
    int HR[256];
    int HG[256];
    int HB[256];
   // tuple<int,int,int> Hpixels[256];

    for (int i = 0; i < 256; ++i)
    {
        HR[i]=0;
        HG[i]=0;
        HB[i]=0;
        //Hpixels[i] = make_tuple(0,0,0);
    }

    Image tempImage(src_image.n_rows, src_image.n_cols);
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            tempImage(y,x) = make_tuple(0,0,0);
        }
    }

    int temp=0,med=0;
    for (int dx = -radius, kx=0; dx <= radius; kx++, dx++) // init
    {
        for (int dy = -radius, ky=0; dy <= radius; ky++, dy++)
        {                    
            HR[ get<0>( src_image(radius+dy, radius+dx ) ) ]++;
            HG[ get<1>( src_image(radius+dy, radius+dx ) ) ]++;
            HB[ get<2>( src_image(radius+dy, radius+dx ) ) ]++;
           // Hpixels[GetY( src_image(radius+dy, radius+dx ) ) ] = src_image(radius+dy, radius+dx );
        }
    }
    for (int i = 0; i < 256; ++i)
    {
        if( (temp+=HR[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
    }
    get<0>(tempImage(radius,radius)) = med; //H*[med] for psychedelic result

    temp=0;
    for (int i = 0; i < 256; ++i)
    {
        if( (temp+=HG[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
    }
    get<1>(tempImage(radius,radius)) = med; //H*[med] for psychedelic result

    temp=0;
    for (int i = 0; i < 256; ++i)
    {
        if( (temp+=HB[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
    }
    get<2>(tempImage(radius,radius)) =med; //H*[med] for psychedelic result




    int y;
  
    for( y = radius; y < src_image.n_rows - radius; y++ )    
    {     
        
        for(int x = (y==radius? radius+1:0) ; x < src_image.n_cols; x++ )
        {
            for (int k = -radius; k <= radius; ++k)
            {
                HR[ get<0>( src_image(y+k - (x- radius-1 < 0), x- radius-1+ (x- radius-1 < 0 ?src_image.n_cols:0 )  ) )]--; //wrapping around
                HR[ get<0>( src_image(y+k + (x+radius >= src_image.n_cols && y+1+k< src_image.n_rows ) , (x+ radius) % src_image.n_cols ) )]++;

                HG[ get<1>( src_image(y+k - (x- radius-1 < 0), x- radius-1+ (x- radius-1 < 0 ?src_image.n_cols:0 )  ) )]--; //wrapping around
                HG[ get<1>( src_image(y+k + (x+radius >= src_image.n_cols && y+1+k< src_image.n_rows ) , (x+ radius) % src_image.n_cols ) )]++;

                HB[ get<2>( src_image(y+k - (x- radius-1 < 0), x- radius-1+ (x- radius-1 < 0 ?src_image.n_cols:0 )  ) )]--; //wrapping around
                HB[ get<2>( src_image(y+k + (x+radius >= src_image.n_cols && y+1+k< src_image.n_rows ) , (x+ radius) % src_image.n_cols ) )]++;
                
               // Hpixels[GetY( src_image(y+k + (x+radius >= src_image.n_cols && y+1+k< src_image.n_rows ) , (x+ radius) % src_image.n_cols ) )] = 
                //src_image(y+k + (x+radius >= src_image.n_cols && y+1+k< src_image.n_rows) , (x+ radius) % src_image.n_cols );
            }

            temp = 0;
            for (int i = 0; i < 256; ++i)
            {
                if( (temp+=HR[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
            }
            get<0>(tempImage(y,x)) = med;

            temp = 0;
            for (int i = 0; i < 256; ++i)
            {
                if( (temp+=HG[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
            }
            get<1>(tempImage(y,x)) = med;

            temp = 0;
            for (int i = 0; i < 256; ++i)
            {
                if( (temp+=HB[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
            }
            get<2>(tempImage(y,x)) = med;

        }
    }
    
    tempImage = tempImage.submatrix(radius,radius,tempImage.n_rows -2*radius,tempImage.n_cols -2*radius );
    return tempImage;
}

void Hadd(int* H1, int * H2)
{
    for (int i = 0; i < 256; ++i)
    {
        H1[i]+= H2[i];       
    }
}

void Hsub(int* H1, int * H2)
{
    for (int i = 0; i < 256; ++i)
    {
        H1[i]-= H2[i];       
    }
}

void Hpixeladd(tuple<int,int,int>* H1, tuple<int,int,int> * H2,int * H2flag)
{
    for (int i = 0; i < 256; ++i)
    {
        if(H2flag[i])
        H1[i]= H2[i];       
    }
}
int sigmoid(int a, int R)
{
    //static double scale = (2*R+1)*(2*R+1)/8.0;
   return (1.0/(1+ exp(-1.0/R*a))-1.0/2)*255.0*2;
}

Image median_const(Image src_image, int radius) { //a.k.a. The Initialization Hell

    src_image = mirror(src_image, 2* radius+1, 2* radius+1);
    int HR[256];
    int HG[256];
    int HB[256];
    //tuple<int,int,int> Hpixels[256];

    int **HcolR = new int*[src_image.n_cols];
    int **HcolG = new int*[src_image.n_cols];
    int **HcolB = new int*[src_image.n_cols];
    //tuple<int,int,int> **Hcolpixels = new tuple<int,int,int>*[src_image.n_cols];
    for (int i = 0; i < src_image.n_cols; ++i) // get memory for those cols histograms
    {
        HcolR[i] = new int[256];
        HcolG[i] = new int[256];
        HcolB[i] = new int[256];
        //Hcolpixels[i] = new tuple<int,int,int>[256];

        for (int j = 0; j < 256; ++j)
        {
            HcolR[i][j]=0;
            HcolG[i][j]=0;
            HcolB[i][j]=0;
            //Hcolpixels[i][j] = make_tuple(0,0,0);
        }
    }


    for (int i = 0; i < 256; ++i)
    {
        HR[i]=0;
        HG[i]=0;
        HB[i]=0;
       // Hpixels[i] = make_tuple(0,0,0);
    }

    Image tempImage(src_image.n_rows, src_image.n_cols); //create image to put pixels in
    for(int x =0; x < src_image.n_cols; x++ )
    {
        for(int y =0; y < src_image.n_rows; y++ )
        {
            tempImage(y,x) = make_tuple(0,0,0);
        }
    }

    

    for(int x = 0; x < src_image.n_cols; x++ ) // initialize rows
    {
        for(int y = 0; y < radius*2+1; y++ )
        {
            HcolR[x][get<0>(src_image(y,x))]++;
            HcolG[x][get<1>(src_image(y,x))]++;
            HcolB[x][get<2>(src_image(y,x))]++;
            //Hcolpixels[x][GetY(src_image(y,x))]= src_image(y,x);
        }
    }


    for (int dx = -radius, kx=0; dx <= radius; kx++, dx++)  // initialize kernel
    {
        Hadd(HR, HcolR[kx]);
        Hadd(HG, HcolG[kx]);
        Hadd(HB, HcolB[kx]);
        //Hpixeladd(Hpixels, Hcolpixels[kx],Hcol[kx]);
    }

    int temp=0,med=0; // get first pixel median
    for (int i = 0; i < 256; ++i)
    {
        if( (temp+=HR[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
    }
    get<0>(tempImage(radius,radius)) = sigmoid( HR[med], radius); 

    temp=0;
    for (int i = 0; i < 256; ++i)
    {
        if( (temp+=HG[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
    }
    get<1>(tempImage(radius,radius)) =sigmoid( HG[med], radius); 

    temp=0;
    for (int i = 0; i < 256; ++i)
    {
        if( (temp+=HB[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
    }
    get<2>(tempImage(radius,radius)) = sigmoid( HB[med], radius); 

    for(int x = radius+1; x < src_image.n_cols; x++ ) // filter first row (of the output) ( it is special)
    {
        if( x+radius >= src_image.n_cols) //preparation for wrapping
        {
            HcolR[(x+radius) % src_image.n_cols ][ get<0>( src_image(0, (x+radius) % src_image.n_cols ) ) ] --;
            HcolR[(x+radius) % src_image.n_cols ][ get<0>( src_image(1, (x+radius) % src_image.n_cols ) ) ] ++;

            HcolG[(x+radius) % src_image.n_cols ][ get<1>( src_image(0, (x+radius) % src_image.n_cols ) ) ] --;
            HcolG[(x+radius) % src_image.n_cols ][ get<1>( src_image(1, (x+radius) % src_image.n_cols ) ) ] ++;

            HcolB[(x+radius) % src_image.n_cols ][ get<2>( src_image(0, (x+radius) % src_image.n_cols ) ) ] --;
            HcolB[(x+radius) % src_image.n_cols ][ get<2>( src_image(1, (x+radius) % src_image.n_cols ) ) ] ++;

            //Hcolpixels[(x+radius) % src_image.n_cols ][ GetY( src_image(1, (x+radius) % src_image.n_cols ) ) ] = src_image(1, (x+radius) % src_image.n_cols ) ;
        }
        Hsub(HR, HcolR[x-radius-1]);
        Hadd(HR, HcolR[(x+radius) % src_image.n_cols ]); 

        Hsub(HG, HcolG[x-radius-1]);
        Hadd(HG, HcolG[(x+radius) % src_image.n_cols ]); 

        Hsub(HB, HcolB[x-radius-1]);
        Hadd(HB, HcolB[(x+radius) % src_image.n_cols ]); 
        //Hpixeladd(Hpixels, Hcolpixels[(x+radius) % src_image.n_cols],Hcol[(x+radius) % src_image.n_cols ]);

        temp = 0;
        for (int i = 0; i < 256; ++i)
        {
            if( (temp+=HR[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
        }
        get<0>(tempImage(radius,x)) = sigmoid( HR[med], radius) ; 

        temp=0;
        for (int i = 0; i < 256; ++i)
        {
            if( (temp+=HG[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
        }
        get<1>(tempImage(radius,x)) = sigmoid( HG[med], radius); 

        temp=0;
        for (int i = 0; i < 256; ++i)
        {
            if( (temp+=HB[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
        }
        get<2>(tempImage(radius,x)) =sigmoid( HB[med], radius); 
        //tempImage(radius,x) = Hpixels[med];

    }



    
    for(int y = radius+1; y < src_image.n_rows - radius; y++ )    // filter from the second row(of output) // the actual algorythm
    {     
        
        for(int x = 0 ; x < src_image.n_cols; x++ )
        {
            
            HcolR[ (x+radius) % src_image.n_cols ][ get<0>( src_image(y-radius-1, (x+radius) % src_image.n_cols ) ) ]--;
            HcolR[ (x+radius) % src_image.n_cols ][ get<0>( src_image(y+radius, (x+radius) % src_image.n_cols ) ) ]++;

            HcolG[ (x+radius) % src_image.n_cols ][ get<1>( src_image(y-radius-1, (x+radius) % src_image.n_cols ) ) ]--;
            HcolG[ (x+radius) % src_image.n_cols ][ get<1>( src_image(y+radius, (x+radius) % src_image.n_cols ) ) ]++;

            HcolB[ (x+radius) % src_image.n_cols ][ get<2>( src_image(y-radius-1, (x+radius) % src_image.n_cols ) ) ]--;
            HcolB[ (x+radius) % src_image.n_cols ][ get<2>( src_image(y+radius, (x+radius) % src_image.n_cols ) ) ]++;
            //Hcolpixels[ (x+radius) % src_image.n_cols ][ GetY( src_image(y+radius, (x+radius) % src_image.n_cols ) ) ] = src_image(y+radius, (x+radius) % src_image.n_cols );

            Hsub(HR, HcolR[x-radius-1 + (x- radius-1 < 0 ? src_image.n_cols : 0 ) ]);
            Hadd(HR, HcolR[(x+radius) % src_image.n_cols ]);

            Hsub(HG, HcolG[x-radius-1 + (x- radius-1 < 0 ? src_image.n_cols : 0 ) ]);
            Hadd(HG, HcolG[(x+radius) % src_image.n_cols ]);

            Hsub(HB, HcolB[x-radius-1 + (x- radius-1 < 0 ? src_image.n_cols : 0 ) ]);
            Hadd(HB, HcolB[(x+radius) % src_image.n_cols ]);
            //Hpixeladd(Hpixels, Hcolpixels[(x+radius) % src_image.n_cols], Hcol[(x+radius) % src_image.n_cols ]); 

            temp = 0;
            for (int i = 0; i < 256; ++i)
            {
                if( (temp+=HR[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
            }
            get<0>(tempImage(y,x)) = sigmoid( HR[med], radius);

            temp = 0;
            for (int i = 0; i < 256; ++i)
            {
                if( (temp+=HG[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
            }
            get<1>(tempImage(y,x)) = sigmoid( HG[med], radius);

            temp = 0;
            for (int i = 0; i < 256; ++i)
            {
                if( (temp+=HB[i]) > (radius*2+1)*(radius*2+1)/2 ) { med = i; break; }
            }
            get<2>(tempImage(y,x)) = sigmoid( HB[med], radius);

            //tempImage(y,x) = Hpixels[med];

        }
    }
    tempImage = tempImage.submatrix(radius,radius,tempImage.n_rows -2*radius,tempImage.n_cols -2*radius );
    return grayScale(tempImage);
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
