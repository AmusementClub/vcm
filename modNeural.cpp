
/************************************************************
Neural filter plugin for vapoursynth by V.C.Mohan
Neural filter trains on input clip using a trainer clip first frames, then 
processes input frames. Uses Resilient propogation. Linear type.
In case of RGB only Green channel while YUV formats Y is used for training.  
Author V.C.Mohan
Aug 26, 2017, Aug 20, 2020
********************************************************************************/  


typedef struct {

	VSNodeRef *node, *tnode;
	const VSVideoInfo *vi, *tvi;

	int xpts;	// number of points along x axis to be used.. Only odd numbers 1 to 121
	int ypts;	// number of points along y axis to be used. Only odd numbers 1 to 121
	int tlx, tty, trx, tby;	// trainer clip window coordinates

	int iter;	// number of iterations for training

	bool wset;	// wset  New set of random numbers start 
	int best;	// use best weights out of this number of different starting weights
	
	int rgb;	// color to use in training 0 for red, 1 for green, 2 for blue;
	
	float *weights;	// weights		
	int inodes;	//	(xpts*ypts+1)
	float bias;

	char txt[8];
	char  filename[200];	// file to be read or saved

} NeuralData;

//--------------------------------------------------------

	template <typename finc>
	void getcase(const finc * fp, int * offsets, float * input, int nodes);
	template <typename finc>
	finc clampval(const float val, finc min, finc max);
	float	getOutput(float * input, int inodes, float * weights);
	int sign(float val);
	void adjustWeights(float *dwtold, float *dwtnew, float * delta, float * deltaw, float * weights, int inodes,
		float deltamin, float deltamax, float yetaminus, float yetaplus);
	void sumdedwt(float err, float *dwtnew, float * input, int inodes);
	
	int filloffsets(int pitch, int kb, int * offsets, int xpts, int ypts);
	

 //---------------------------------------------------------------------------------
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC neuralInit
(VSMap *in, VSMap *out, void **instanceData,
VSNode *node, VSCore *core, const VSAPI *vsapi)
{
	NeuralData *d = (NeuralData *)* instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);
	
	if (strcmp(d->txt, "read") != 0)
	{
		//vsapi->setVideoInfo(d->tvi, 1, node);

		d->inodes = d->xpts* d->ypts + 1;
		d->bias = 1.0; // / xpts;	// arbitrarily selected
		float deltamax = 5.0f * 0.001f;
		float deltamin = exp(-6.0f) * 0.001f;
		float yetaplus = 1.2f;
		float yetaminus = 0.5f;

		const VSFrameRef *src = vsapi->getFrame(0, d->node, NULL, 0);
		const VSFrameRef *trf = vsapi->getFrame(0, d->tnode, NULL, 0);

		const VSFormat *fi = d->vi->format;
		int plane = fi->colorFamily == cmRGB ? 2 - d->rgb : 0;

		const uint8_t *srcp = vsapi->getReadPtr(src, plane);
		int src_stride = vsapi->getStride(src, plane);
		const uint8_t *tp = vsapi->getReadPtr(trf, plane);
		int t_stride = vsapi->getStride(trf, plane);
		int bht = vsapi->getFrameHeight(src, plane);
		int bwd = vsapi->getFrameWidth(src, plane);
		int nbytes = fi->bytesPerSample;
		int spitch = src_stride / nbytes;
		int tpitch = t_stride / nbytes;
		int nbits = fi->bitsPerSample;
		// number of points available to train in each iteration
		int xcases = (d->trx - d->tlx - d->xpts + 1);	// along x 
		int ycases = (d->tby - d->tty - d->ypts + 1);	// along y
		int cases = xcases*ycases;				// total points to train for each iteration
		int inodes = d->inodes;
		d->weights = new float[inodes];		// buffer to hold layer weights and   bias
		int * offsets = new int[d->xpts * d->ypts];	// offsets from central point for each of Xpts X ypts around it
		float *buf = new float[d->inodes * 6];
		float *whcopy = buf;		// buffer to hold best  layer weights plus  bias
		float *dedwt1 = whcopy + inodes;		// changing wts
		float *dedwt2 = dedwt1 + inodes;
		float *deltaw = dedwt2 + inodes;
		float *delta = deltaw + inodes;
		float *input = delta + inodes;	// buffer to hold input plus one bias node
		bool save = false;
		float *sumerr = NULL;
		if (strcmp(d->txt, "save") == 0)
		{
			// create buffer to show error for each iteration

			sumerr = new float[d->best * d->iter];
			save = true;
		}

		time_t seed;
		// use time for seeding random numbers if wset is true. Else use standard system generated
		if (d->wset)
		{
			time(&seed);
			srand((seed & 0xfffe) + 1);
		}
		else
		{
			srand(1);
		}
		// generate offsets table
		filloffsets(spitch, 1, offsets, d->xpts, d->ypts);

		float hf = 0.001f / d->inodes;	// arbitrary value for scaling of random weights
		input[d->inodes - 1] = d->bias; // bias value. remains same for all cases
		int best, niter;
		float minesum;

		for (int b = 0; b < d->best; b++)
		{
			for (int i = 0; i < d->inodes; i++)
			{
				// input layer weights between -0.5hf and 0.5hf
				d->weights[i] = (((float)rand()) / RAND_MAX - 0.5) * hf;
				// initialize buffers
				delta[i] = deltamin;
				deltaw[i] = 0;
				dedwt1[i] = 0;
			}


			// for toggling buffers for each iteration
			float *dedwtnew, *dedwtold;
			float esum;	// minimum sum of errors

			for (int i = 0; i < d->iter; i++)
			{
				if ((i & 1) == 0)
				{
					// toggle buffers for each odd number of iteration
					dedwtold = dedwt1;
					dedwtnew = dedwt2;
				}
				else
				{
					// toggle buffers at each even number of iteration
					dedwtold = dedwt2;
					dedwtnew = dedwt1;
				}

				esum = 0.0;	// initialize error squared sum

				float output, error;
				// zero de/dw buffer
				for (int nodes = 0; nodes < d->inodes; nodes++)
				{
					dedwtnew[nodes] = 0;
				}

				// trainer frame ptr and input frame ptrs are to point center of the first rectangle
				const uint8_t * trp = tp + (d->tty + d->ypts / 2) * t_stride + (d->tlx + d->xpts / 2) * nbytes;
				const uint8_t * srp = srcp + (d->tty + d->ypts / 2) * src_stride + (d->tlx + d->xpts / 2) * nbytes;
				// for each input point at center of square find output and sum error square

				for (int h = 0; h < ycases; h++)
				{
					for (int w = 0; w < xcases; w++)
					{
						if (nbytes == 1)
						{
							// get next case for training
							getcase(srp + w, offsets, input, d->inodes - 1);		// got  input data and bias in input buffer
							// using current weights get output
							output = getOutput(input, d->inodes, d->weights);
							// error for this point using current weights
							error = *(trp + w) - output;	// error of estimation
						}
						else if (nbytes == 2)
						{
							getcase((uint16_t*)srp + w, offsets, input, d->inodes - 1);

							output = getOutput(input, d->inodes, d->weights);

							error = *((uint16_t*)trp + w) - output;

						}
						else if (nbytes == 4)
						{
							// get next case for training						 
							getcase((float*)srp + w, offsets, input, d->inodes - 1);

							output = getOutput(input, d->inodes, d->weights);

							error = *((float*)trp + w) - output;
						}

						esum += error*error;
						sumdedwt(error, dedwtnew, input, d->inodes);

						
					}					

					

					srp += src_stride;
					trp += t_stride;
				}

				if (save)
					sumerr[b * d->iter + i] = esum;

				if (i == 0 && b == 0)
				{
					// esum is now error squared sum for first iteration 
					minesum = esum;	// to get minimum error point

					// in unlikely case first weights are optimum weights
					for (int nodes = 0; nodes < d->inodes; nodes++)
						whcopy[nodes] = d->weights[nodes];
				}

				if (minesum > esum)
				{
					for (int nodes = 0; nodes < d->inodes; nodes++)
						whcopy[nodes] = d->weights[nodes];

					minesum = esum;
					best = b;
					niter = i;
				}

				// adjust weights so that output value comes  closer to desired value
				adjustWeights(dedwtold, dedwtnew, delta, deltaw, d->weights, d->inodes,
					deltamin, deltamax, yetaminus, yetaplus);
	
			}	// iterations

		}	// best of

		for (int nodes = 0; nodes < inodes; nodes++)
		{
			d->weights[nodes] = whcopy[nodes];
		}
		delete[]buf;
		delete[]offsets;
		vsapi->freeFrame(src);
		vsapi->freeFrame(trf);

		// training end

		if (strcmp(d->txt, "save") == 0)
		{
			
			std::ofstream ofile;
			ofile.open(d->filename);

			if (!ofile.is_open())
			{
				vsapi->setError(out, "Neural: Could not open output file");
				vsapi->freeNode(d->node);
				vsapi->freeNode(d->tnode);
				return;
			}

			char nl = '\n';

			char * title = "vcMod_Neural_xpts_ypts_inodes_and_Weights";

			ofile << title << nl;

			char space = ' ';

			if (fi->colorFamily == cmRGB)
				ofile << "RGB";
			else if(fi->colorFamily == cmYUV)
				ofile << "YUV"; 
			else
				ofile << "GREY";	

			ofile << space << "bitdepth" << space << nbits << nl;

			// ofile << "pixel_type" << space << vi.pixel_type << nl;

			ofile << d->xpts << space << d->ypts << space << d->inodes << nl; // weights includes bias
			ofile << "bias" << space << d->bias << nl;


			int count = 0;

			for (int i = 0; i < d->inodes; i++)
			{

				ofile << d->weights[i] << nl;		// float value of weight


			}

			ofile << "Time Stamp" << nl;
			time_t timeStamp;
			time(&timeStamp);
			ofile << timeStamp << nl;

			ofile << "minimum_esum " << minesum << " at_iter " << niter << " at_weight_set_no: " << best
				<< " for " << cases << " training points" << nl;

			ofile << "error_sum_at_each_iteration" << nl;

			int in = 0;

			for (int bi = 0; bi < d->best * d->iter; bi ++)
			{
				if ((bi % d->iter) == 0)
				{
					ofile << nl << "with_weight_set_No:" << bi / d->iter << nl;
					in = 0;
				}

				if ((in % 20) == 0)
				{
					ofile  << nl <<"error_sums_at" << "iter:" << in << "to" << space << in + 19 << ":-" ;

				}

				if ((in % 5) == 0)
				{
					ofile << nl;
				}

				ofile << sumerr[bi] << space;
				in++;

			}

			ofile.close();

			delete[]sumerr;
		}
	}

	else	// read from file
	{
		std::ifstream ifile;
		
		ifile.open(d->filename);

		if (!ifile.is_open())
		{
			vsapi->setError(out, "Neural: Could not open input file");
			vsapi->freeNode(d->node);
			return;
		}

		const VSFormat *fi = d->vi->format;


		char * title = "vcMod_Neural_xpts_ypts_inodes_and_Weights";
		char buf[256];
		buf[255] = '\0';

		ifile >> buf;

		if (strcmp(buf, title) != 0)
		{
			ifile.close();
			vsapi->setError(out, "Neural: Incorrect input file");
			vsapi->freeNode(d->node);
			return;
		}

		int bitdepth;

		ifile >> buf;
		

		if ((strcmp(buf, "RGB") == 0 && fi->colorFamily != cmRGB))
		//	|| (strcmp(buf, "YUV") == 0 && fi->colorFamily != cmYUV)
		//	|| (strcmp(buf, "GREY") == 0 && fi->colorFamily != cmGray)
		
		{

			ifile.close();
			vsapi->setError(out, "Neural: input file was for different colorFamily");
			vsapi->freeNode(d->node);
			return;

		}
		ifile >> buf >> bitdepth;

		if (bitdepth != fi->bitsPerSample)
		{
			ifile.close();
			vsapi->setError(out, "Neural: input file was for different  bit depth");
			vsapi->freeNode(d->node);
			return;
		}


		int npoints = 0;
		ifile >> d->xpts >> d->ypts >> d->inodes >>buf >> d->bias;
		

		if (d->xpts * d->ypts < 9 || d->xpts * d->ypts > 225 || ((d->xpts * d->ypts) & 1) == 0 || d->xpts * d->ypts + 1 != d->inodes)
		{
			ifile.close();
			vsapi->setError(out, "Neural: input file is corrupted");
			vsapi->freeNode(d->node);
			return;
		}

		d->weights = new float[d->inodes];

		int count = 0;

		for (int i = 0; i < d->inodes && !ifile.eof(); i++)
		{

			ifile >> d->weights[i];		// float value of weight

			count++;
		}

		if (count != d->inodes)
		{
			ifile.close();
			vsapi->setError(out, "Neural: input file has fewer weights");
			delete[] d->weights;
			vsapi->freeNode(d->node);
			return;
		}

		ifile.close();
	}

	
 }

/***************************************************************/	
 //---------------------------------------------------------------------------------	

 // This is the main function that gets called when a frame should be produced. It will in most cases get
 // called several times to produce one frame. This state is being kept track of by the value of
 // activationReason. The first call to produce a certain frame n is always arInitial. In this state
 // you should request all input frames you need. Always do it i ascending order to play nice with the
 // upstream filters.
 // Once all frames are ready the the filter will be called with arAllFramesReady. It is now time to
 // do the actual processing.
 static const VSFrameRef *VS_CC neuralGetFrame
	 (int n, int activationReason, void **instanceData,
	 void **frameData, VSFrameContext *frameCtx,
	 VSCore *core, const VSAPI *vsapi)
 {
	NeuralData *d = (NeuralData *)* instanceData;

	 if (activationReason == arInitial)
	 {
		 // Request the source frame on the first call
		 vsapi->requestFrameFilter(n, d->node, frameCtx);
	 }
	 else if (activationReason == arAllFramesReady)
	 {
		 const VSFrameRef *src = vsapi->getFrameFilter(n, d->node, frameCtx);
		 // The reason we query this on a per frame basis is because we want our filter
		 // to accept clips with varying dimensions. If we reject such content using d->vi
		 // would be better.
		 const VSFormat *fi = d->vi->format;

		 // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
		 // supply the "dominant" source frame to copy properties from. Frame props
		 // are an essential part of the filter chain and you should NEVER break it.
		 VSFrameRef * dst = vsapi->copyFrame(src, core);
	//	 int * off = (int*)vs_aligned_malloc((size_t)(d->inodes * sizeof(int)), (size_t)32);
	//	 float * input = (float*)vs_aligned_malloc((size_t)(d->inodes * sizeof(float)), (size_t)64);

		 // get memory from vs and not from system
		 VSFrameRef *buf = vsapi->newVideoFrame(fi, d->inodes * (sizeof(int) + sizeof(float)), 1, src, core);
		 uint8_t *bp = vsapi->getWritePtr(buf, 0);
		 float * input = (float *)bp;

		 int * offs = (int*)(input + d->inodes );

		 int np = fi->colorFamily == cmRGB ? 3 : 1;
			 
		 // It's processing loop time!
		 // Loop over all the planes

		 for (int plane = 0; plane < np; plane++)
		 {

			 const uint8_t *sp = vsapi->getReadPtr(src, plane);
			 int src_stride = vsapi->getStride(src, plane);
			 uint8_t *dp = vsapi->getWritePtr(dst, plane);
			 int dst_stride = vsapi->getStride(dst, plane);
			 int bht = vsapi->getFrameHeight(src, plane);
			 int bwd = vsapi->getFrameWidth(src, plane);
			 int nbytes = fi->bytesPerSample;
			 int spitch = src_stride / nbytes;
			 int dpitch = dst_stride / nbytes;
			 int nbits = fi->bitsPerSample;
			 

			 if (plane == 0)
			 {
				 filloffsets(spitch, 1, offs, d->xpts, d->ypts);
			 }


			 input[d->inodes - 1] = d->bias;		// bias
			 // process regular 		

			 sp += (d->ypts / 2) * src_stride;
			 dp += (d->ypts / 2) * dst_stride;

			 for (int h = d->ypts / 2; h < bht - d->ypts / 2 - 1; h++)
			 {
				 //int kbw = kb * (xpts / 2);

				 for (int w = d->xpts / 2; w < bwd - d->xpts / 2 - 1; w++)
				 {
					 if (nbytes == 1)
					 {
						 // get input data + bias in input buffer
						 getcase(sp + w, offs, input, d->inodes - 1);

						 float output = getOutput(input, d->inodes, d->weights);
						 unsigned char min = 0, max = 255;
						 *(dp + w) = clampval(output, min, max);
					 }
					 if (nbytes == 2)
					 {
						 // get input data + bias in input buffer
						 getcase((uint16_t*)sp + w, offs, input, d->inodes - 1);

						 float output = getOutput(input, d->inodes, d->weights);

						 uint16_t min = 0, max = (1 << nbits) - 1;

						 *((uint16_t*)dp + w) = clampval(output, min, max);

					 }
					 if (nbytes == 4)
					 {
						 // get input data + bias in input buffer
						 getcase((float*)sp + w, offs, input, d->inodes - 1);

						 float output = getOutput(input,  d->inodes, d->weights);

						 float min = 0, max = 1.0f;

						 *((float*)dp + w) = clampval(output, min, max);
					 }

					 // kbw += kb;
				 }

				 dp += dst_stride;
				 sp += src_stride;

			 }
		 }
		 
		 vsapi->freeFrame(buf);
		 vsapi->freeFrame(src);

		 return dst;
	 }
	 return 0;
}
//-------------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC neuralFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
	NeuralData *d = (NeuralData *)instanceData;
	vsapi->freeNode(d->node);
	free(d);
}

//------------------------------------------------------------------
template <typename finc>
finc clampval(const float val, finc min, finc max)
{
	return val > max ? max : val < min ? min : val;
}


//---------------------------------------------------------------
 int filloffsets(int fpitch, int kb, int * offsets, int xpts, int ypts)
{
	int n = 0;
	for (int h = -ypts / 2; h <= ypts / 2; h++)
	{
		for (int w = -xpts / 2; w <= xpts / 2; w++)
		{
			offsets[n] = h * fpitch + w * kb;
			n++;
		}
	}
	return n;
}
//--------------------------------------------------------------
float getOutput(float * input , int inodes, float * weights)
{						
								
	float sum=0.0;
	for ( int in=0; in<inodes; in++)
	{
		sum+= input[in] * weights[in];
	}
	return(sum);								
}
//-------------------------------------------------------------------------
template <typename finc>
void  getcase(const finc * fp, int * offsets, float * input, int nodes)
{
	
	for (int i = 0; i < nodes; i++)
	{
		input[i] = fp[offsets[i]];
	}

	
}
//----------------------------------------------------------------
void  sumdedwt(float err, float *dwtnew, float * input, int inodes)
{
	for ( int in=0; in<inodes; in++)
	
		dwtnew[in]-= err*input[in];
}

//------------------------------------------------------------------------------------------------------

void  adjustWeights(float *dwtold, float *dwtnew, float * delta, float * deltaw, float * weights, int inodes,
	float deltamin, float deltamax, float yetaminus, float yetaplus)
{
	// resilient propogation

	for (int hn = 0; hn<inodes; hn++)
	{

		if (dwtold[hn] * dwtnew[hn] >0)
		{	// the de/dwt has not changed sign

			// use minimum of (deltamax and delta[]*yetaplus)
			delta[hn] = deltamax < delta[hn] * yetaplus ?
			deltamax : delta[hn] * yetaplus;
			// sign(val) returns -1, 0 or +1 based on val sign)
			deltaw[hn] = -sign(dwtnew[hn])* delta[hn];

			weights[hn] += deltaw[hn];
		}

		else if (dwtold[hn] * dwtnew[hn] <0)
		{
			// de/dwt has changed sign
			// use max of ( deltamin and delta[]*yetaminus)
			delta[hn] = deltamin >= delta[hn] * yetaminus ?
			deltamin : delta[hn] * yetaminus;
			weights[hn] -= deltaw[hn];// go back to prev value as we over stepped a minimum

			dwtnew[hn] = 0;	// to ensure that next iter when we come to this function  no update of this weight occurs


		}

		else // zero
		{
			deltaw[hn] = -sign(dwtnew[hn]) * delta[hn];

			weights[hn] += deltaw[hn];
		}
		
	}

}
//---------------------------------------------------------------------
int sign( float val)
{
	if (val > 0)
		return (1);
	if(val<0)
		return (-1);
	return (0);
}
///-----------------------------------------------------------------
// This function is responsible for validating arguments and creating a new filter
static void VS_CC neuralCreate(const VSMap *in,
	VSMap *out, void *userData, VSCore *core,
	const VSAPI *vsapi)
{
	NeuralData d;
	NeuralData *data;
	int temp;
	int err;
	const char * txt, *fname;
	txt = vsapi->propGetData(in, "txt", 0, &err);
	if (err)
		txt = "none";
	else if (strcmp(txt, "none") != 0 && strcmp(txt, "read") != 0 && strcmp(txt, "save") != 0)
	{
		vsapi->setError(out, "neural: txt can be either \"none\" or \"read\" or \"save\" only");
		return;
	}

	strcpy_s(d.txt, 8, txt);

	if (strcmp(txt, "none") != 0)
	{
		fname = vsapi->propGetData(in, "fname", 0, &err);
		if (err)
		{
			vsapi->setError(out, "neural: fname file name with full path must be specified");
			return;
		}

		strcpy_s(d.filename,200, fname);
	}

	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);

	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV &&
		d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "neural: clips can be RGB or YUV or Gray color formats only");
		vsapi->freeNode(d.node);		
		return;
	}

	if (strcmp(txt, "read") != 0)
	{
			// for save or none options training is required
		d.tnode = vsapi->propGetNode(in, "tclip", 0, &err);
		if (err)
		{
			vsapi->setError(out, "neural:  tclip must be specified for this txt option");
			vsapi->freeNode(d.node);
			return;
		}

		d.tvi = vsapi->getVideoInfo(d.tnode);

		if (!isConstantFormat(d.vi) && !isSameFormat(d.vi, d.tvi))
		{
			vsapi->setError(out, "neural: input and tclip must have same constant format");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.tnode);
			return;
		}


		d.xpts = vsapi->propGetInt(in, "xpts", 0, &err);
		if (err)
			d.xpts = 3;
		d.ypts = vsapi->propGetInt(in, "ypts", 0, &err);
		if (err)
			d.ypts = d.xpts;

		if (d.xpts < 1 || d.ypts < 1 || d.xpts * d.ypts < 9 || d.xpts * d.ypts > 15 * 15 || ((d.xpts * d.ypts) & 1) == 0)
		{
			vsapi->setError(out, "neural: xpts and ypts must be positive odd numbers with a product between 9 and 225. ");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.tnode);
			return;
		}
		if (d.vi->format->colorFamily == cmRGB)
		{
			d.rgb = vsapi->propGetInt(in, "rgb", 0, &err);
			if (err)
			{
				d.rgb = 1;
			}

			else if (d.rgb < 0 || d.rgb > 2)
			{
				vsapi->setError(out, "neural: rgb value can be 0 for red, 1 for green and 2 for blue use in training purpose ");
				vsapi->freeNode(d.node);
				vsapi->freeNode(d.tnode);
				return;
			}
		}

		

		// Note that
		// vi->format can be 0 if the input clip can change format midstream.
		//if (!isConstantFormat(d.vi) || d.vi->format->sampleType != stInteger ) {
		//    vsapi->setError(out, "AdaptiveMedian: only constant format 8 bit integer input supported");
		//    vsapi->freeNode(d.node);
		//    return;


		// If a property read fails for some reason (index out of bounds/wrong type)
		// then err will have flags set to indicate why and 0 will be returned. This
		// can be very useful to know when having optional arguments. Since we have
		// strict checking because of what we wrote in the argument string the only reason
		// this could fail is when the value wasn't set by the user.
		// And when it's not set we want it to default to preset value.

		d.tlx = vsapi->propGetInt(in, "tlx", 0, &err);
		if (err)
			d.tlx = d.xpts;
		d.tty = vsapi->propGetInt(in, "tty", 0, &err);
		if (err)
			d.tty = d.ypts;
		d.trx = vsapi->propGetInt(in, "trx", 0, &err);
		if (err)
			d.trx = d.vi->width - d.xpts;
		d.tby = vsapi->propGetInt(in, "tby", 0, &err);
		if (err)
			d.tby = d.vi->height - d.ypts;

		if (d.tlx < d.xpts / 2 || d.tty < d.ypts / 2
			|| d.trx > d.vi->width - d.xpts / 2 || d.tby > d.vi->height - d.ypts / 2
			|| (d.trx - d.tlx) * (d.tby - d.tty) < 10000)
		{
			vsapi->setError(out, "neural: trainig window should be in frame with borders of xpts/2 and ypts/2 and have atleast 10000 pixels");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.tnode);
			return;
		}
		d.iter = vsapi->propGetInt(in, "iter", 0, &err);
		if (err)
			d.iter = 200;

		if (d.iter < 1)
		{
			vsapi->setError(out, "neural: number of iterations iter for trainig should be a sufficiently large positive number");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.tnode);
			return;
		}
		d.best = vsapi->propGetInt(in, "bestof", 0, &err);
		if (err)
			d.best = 1;

		if (d.best < 1 || d.best > 10)
		{
			vsapi->setError(out, "neural: bestof can be 1 to 10 only");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.tnode);
			return;
		}

		temp = vsapi->propGetInt(in, "wset", 0, &err);
		if (err)
		{
			d.wset = false;
		}

		else if (temp < 0 || temp > 1)
		{
			vsapi->setError(out, "neural: wset can have a value of 0 or 1 only ");
			vsapi->freeNode(d.node);
			vsapi->freeNode(d.tnode);
			return;
		}
		else
			d.wset = temp == 1 ? true : false;
	}
	
	data = (NeuralData *)malloc(sizeof(d));
	*data = d;

	// Create a new filter and returns a reference to it. Always pass on the in and out
	// arguments or unexpected things may happen. The name should be something that's
	// easy to connect to the filter, like its function name.
	// The three function pointers handle initialization, frame processing and filter destruction.
	// The filtermode is very important to get right as it controls how threading of the filter
	// is handled. In general you should only use fmParallel whenever possible. This is if you
	// need to modify no shared data at all when the filter is running.
	// For more complicated filters fmParallelRequests is usually easier to achieve as an
	// be prefetched in parallel but the actual processing is serialized.
	// The others can be considered special cases where fmSerial is useful to source filters and
	// fmUnordered is useful when a filter's state may change even when deciding which frames to
	// prefetch (such as a cache filter).
	// If you filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.
	vsapi->createFilter(in, out, "neural", neuralInit,
		neuralGetFrame, neuralFree,
		fmParallel, 0, data, core);
	return;
}
/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
configFunc("vc.mohan.dns", "neural", "VapourSynth Neural", VAPOURSYNTH_API_VERSION, 1, plugin);
registerFunc("neural", "clip:clip;tclip:clip;xpts:int:opt;ypts:int:opt;tlx:int:opt;tty:int:opt;trx:int:opt;tby:int:opt;"
"iter:int:opt;bestof:int:opt;wset:int:opt;rgb:int:opt;txt:Data:opt;fname:Data:opt", neuralCreate, 0, plugin);
}
*/
