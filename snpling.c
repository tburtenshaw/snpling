#include <windows.h>
#include <windowsx.h>
#include <commctrl.h>
#include <string.h>
#include <stdio.h>
#include "snplingres.h"

#define SNPSINCHUNK 4096	//Wasting this much space isn't too much of a concern
							//and it's better than calling a malloc 500,000 times
#define	MAXCOMMENTS 2048	//Maximum length of comments
#define MAX_MULTIPATH MAX_PATH * 8	//to load multiple files

HINSTANCE hInst;		// Instance handle
HWND hwndMain;		//Main window handle
struct sSnpEntry;
struct sChunk;	//contains a fix-length array of snps
struct sChromosome;
struct sGenome;
struct sCounts;
struct sGenomeWindowInfo;

LRESULT CALLBACK MainWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam);
LRESULT CALLBACK GenomeInfoProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam);
int GenomeInfoPaint(HWND hwnd);

//Loading
int SnpFromLine(struct sSnpEntry *snp, char *line);
int AddLineToComments(struct sGenome *genome, char *line);
int OpenGenomeDialog(HWND hwnd);
int LoadGenomeFile(HWND hwnd, char *filename);

//Calculation
int IsProperBase(char b);	//is it ATC or G? or nothing
long GenomeSimilarity(struct sGenome *g1, struct sGenome *g2);

//Debugging
int MessageSnp(HWND hwnd, struct sSnpEntry *snp);
int MessageCount(HWND hwnd, struct sCounts *count);
int MessagePoint(HWND hwnd, POINT *p, char *heading);

struct sGenomeWindowInfo
{
	HWND	hwnd;
	POINT	initialClickPoint;	//used for moving
	RECT	originalLocationRect;
	int		originalHeight;
	int		originalWidth;
	int		mouseBeenPushed;	//have we clicked the mouse
};

struct sSnpEntry
{
	char rstype[4];
	long rsnumber;
	int chrom;
	long position;
	char genotype[2];
};

struct sChunk
{
	struct sSnpEntry entry[SNPSINCHUNK];
	struct sChunk *next;
};

struct sCounts	//stats to use both per chromosome, and per genome
{
	unsigned long total;
	unsigned long autosomal;
	unsigned long homozygous;	//heteros are snp minus homo
	unsigned long adenine;
	unsigned long cytosine;
	unsigned long guanine;
	unsigned long thymine;
	unsigned long nocall;
};

struct sChromosome
{
	struct sChunk pChunk;	//pointers to the chunk the specific chromosome is first found
	struct sCounts count;
};

struct sGenome
{
	struct sChunk *firstChunk;
	struct sChromosome chromosome[26];			//info about chromosomes, we use 1-25 (0's wasted)

	struct sCounts count;

	char	filename[MAX_PATH];
	char	*justfilename;	//points to the above without path
	char	owner[256];	//the user name associated (if we can pull it from a file)
	char	comments[MAXCOMMENTS];	//This contains the comments from the file
	long	countXX;		//Does it have two alleles on the x-chromosome (even males have some on 23andme)
	int		containsY;		//Does it have a Y chromosome?
	int		containsM;		//Does the data have mitochondrial DNA?
	int		chromosomeCount;

	struct sGenome *next;

	struct sGenomeWindowInfo windowInfo;

};


int OpenGenomeDialog(HWND hwnd)
{
	OPENFILENAME ofn;
	char *filenamepart;	//hold the pointer to the specific file if a list is returned
	char fullpath[MAX_PATH];
	char buffer[MAX_MULTIPATH];

	int	r;

	memset(&ofn,0,sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hInstance = GetModuleHandle(NULL);
	ofn.hwndOwner = GetActiveWindow();
	ofn.lpstrFile = buffer;
	ofn.nMaxFile = MAX_MULTIPATH;
	ofn.lpstrTitle = "Open";
	ofn.nFilterIndex = 2;
	ofn.lpstrDefExt = "txt";
	strcpy(buffer,"*.txt");
	ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST| OFN_HIDEREADONLY | OFN_ALLOWMULTISELECT|OFN_EXPLORER;
	ofn.lpstrFilter = "All files\0*.*\0Text file (*.txt)\0*.txt\0";

	r=GetOpenFileName(&ofn);

	if (!r)
		return 0;

	filenamepart=ofn.lpstrFile+strlen(ofn.lpstrFile)+1;
	if (!(*filenamepart))	{	//we've only selected one file
		LoadGenomeFile(hwnd, ofn.lpstrFile);
	}
	else	{
		while (*filenamepart)	{
			sprintf(fullpath, "%s\\%s", ofn.lpstrFile, filenamepart);
			LoadGenomeFile(hwnd, fullpath);
			filenamepart+=strlen(filenamepart)+1;
		}
	}

	return  1;	//success
}

int LoadGenomeFile(HWND hwnd, char *filename)
{
	FILE * fDNA;
	struct sChunk *snpChunk;
	struct sChunk *prevChunk;
	struct sGenome *genome;
	struct sGenome *othergenome;
	int		arrayPos;
	POINT	windowloc;

	char	localGenotype[2];	//make it shorter!
	int		localChromosome;
	int		containsRealBase;	//it's not a no-call or a placeholder


	fDNA = fopen(filename, "r");
	if (!fDNA)
		return 0;

	genome = malloc(sizeof(struct sGenome));	//Get memory for new genome
	memset(genome, 0, sizeof(struct sGenome));	//Set it all to zero initially
	//Set some characteristics of the genome
	strcpy(genome->filename, filename);	//Filename
	genome->next = NULL;	//this will be set already (with memset) but I like to do it explicitly

	othergenome = (void *)GetWindowLong(hwnd, GWL_USERDATA);
	windowloc.y=24;
	if (othergenome == 0)	{	//if this is our first genome, then it's stored in main window
		SetWindowLong(hwnd, GWL_USERDATA, (long)genome);
		//MessageBox(hwnd, "Using main window","test",0);
	}
	else	{
		windowloc.y+=130;
		while (othergenome->next)	{
			othergenome=othergenome->next;
			windowloc.y+=130;
			//MessageBox(hwnd, "Trying linked list","test",0);
		}
		othergenome->next = genome;
	}



	snpChunk = malloc(sizeof(struct sChunk));
	genome->firstChunk=snpChunk;
	arrayPos=0;

	char line[256];
	while(fgets(line, sizeof(line), fDNA))     {
		if (line[0]=='#')	{//handle comments line
			AddLineToComments(genome, line);
		}
		else if (SnpFromLine(&snpChunk->entry[arrayPos], line)) {
			localChromosome=snpChunk->entry[arrayPos].chrom;
			localGenotype[0]=snpChunk->entry[arrayPos].genotype[0];
			localGenotype[1]=snpChunk->entry[arrayPos].genotype[1];

			genome->count.total++;
			genome->chromosome[localChromosome].count.total++;

			if (IsProperBase(localGenotype[0]) || IsProperBase(localGenotype[1]))
				{containsRealBase=1;}
			else
				{containsRealBase=0;}

			//Count homozygous
			if ((localGenotype[0] == localGenotype[1]) && (containsRealBase))	{
				genome->count.homozygous++;
				genome->chromosome[localChromosome].count.homozygous++;
				switch (localGenotype[0])	{	//the other will be the same
					case 'A':
						genome->count.adenine+=2;
						genome->chromosome[localChromosome].count.adenine+=2;
					break;
					case 'T':
						genome->count.thymine+=2;
						genome->chromosome[localChromosome].count.thymine+=2;
					break;
					case 'C':
						genome->count.cytosine+=2;
						genome->chromosome[localChromosome].count.cytosine+=2;
					break;
					case 'G':
						genome->count.guanine+=2;
						genome->chromosome[localChromosome].count.guanine+=2;
					break;
				}
			}
			else if (containsRealBase)	{
				switch (localGenotype[0])	{	//the other will be the same
					case 'A':
						genome->count.adenine++;
						genome->chromosome[localChromosome].count.adenine++;
					break;
					case 'T':
						genome->count.thymine++;
						genome->chromosome[localChromosome].count.thymine++;
					break;
					case 'C':
						genome->count.cytosine++;
						genome->chromosome[localChromosome].count.cytosine++;
					break;
					case 'G':
						genome->count.guanine++;
						genome->chromosome[localChromosome].count.guanine++;
					break;
				}
				switch (localGenotype[1])	{	//the other will be the same
					case 'A':
						genome->count.adenine++;
						genome->chromosome[localChromosome].count.adenine++;
					break;
					case 'T':
						genome->count.thymine++;
						genome->chromosome[localChromosome].count.thymine++;
					break;
					case 'C':
						genome->count.cytosine++;
						genome->chromosome[localChromosome].count.cytosine++;
					break;
					case 'G':
						genome->count.guanine++;
						genome->chromosome[localChromosome].count.guanine++;
					break;
				}

			}
			//Count autosomal
			if (localChromosome<23)	{
				genome->count.autosomal++;
				//genome->chromosome[localChromosome].count.autosomal++; //not much sense to count this per chrom
			}
			else	{
				if ((localChromosome==24) && (containsRealBase))	{
					genome->containsY=1;
				}
				else if ((localChromosome==25) && (containsRealBase))	{
					genome->containsM=1;	//contains mitochondrial
				}

				else if ((localChromosome==23) && (containsRealBase))	{
					if 	(IsProperBase(localGenotype[0]) && IsProperBase(localGenotype[1])) {
						genome->countXX++;
					}
				}
			}

			arrayPos++;
			if (arrayPos==SNPSINCHUNK)	{	//if this chunk is full - make a new one
				prevChunk=snpChunk;
				snpChunk = malloc(sizeof(struct sChunk));	//allocate a new chunk
				snpChunk->next = NULL;	//this is the new terminating linkedlist
				prevChunk->next = snpChunk;	//point the old chunk to this new one
				arrayPos=0;
			}
		}
	}
	fclose(fDNA);

	genome->windowInfo.hwnd = CreateWindow("GenomeInfoBoxClass", genome->filename, WS_CHILD|WS_VISIBLE, 32,windowloc.y, 180,120, hwnd, NULL, hInst, NULL);
	SetWindowLong(genome->windowInfo.hwnd, GWL_USERDATA, (long)genome);

	return 0;
}

static BOOL InitApplication(void)
{
	WNDCLASS wc;

	//Main window
	memset(&wc,0,sizeof(WNDCLASS));
	wc.style = CS_DBLCLKS;
	wc.lpfnWndProc = (WNDPROC)MainWndProc;
	wc.hInstance = hInst;
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+2);
	wc.lpszClassName = "snplingWndClass";
	wc.lpszMenuName = MAKEINTRESOURCE(IDMAINMENU);
	wc.hCursor = LoadCursor(NULL,IDC_ARROW);
	wc.hIcon = LoadIcon(hInst,MAKEINTRESOURCE(IDI_FIRST));
	if (!RegisterClass(&wc))
		return 0;

	//Genome infoboxes
	memset(&wc,0,sizeof(WNDCLASS));
	wc.style = CS_DBLCLKS;
	wc.lpfnWndProc = (WNDPROC)GenomeInfoProc;
	wc.hInstance = hInst;
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
	wc.lpszClassName = "GenomeInfoBoxClass";
	wc.lpszMenuName = MAKEINTRESOURCE(IDMAINMENU);
	wc.hCursor = LoadCursor(NULL,IDC_ARROW);
	wc.hIcon = NULL;
	if (!RegisterClass(&wc))
		return 0;


	return 1;
}

HWND CreatesnplingWndClassWnd(void)
{
	return CreateWindow("snplingWndClass","snpling",
		WS_MINIMIZEBOX|WS_VISIBLE|WS_CLIPSIBLINGS|WS_CLIPCHILDREN|WS_MAXIMIZEBOX|WS_CAPTION|WS_BORDER|WS_SYSMENU|WS_THICKFRAME,
		CW_USEDEFAULT,0,CW_USEDEFAULT,0,
		NULL,
		NULL,
		hInst,
		NULL);
}

void MainWndProc_OnCommand(HWND hwnd, int id, HWND hwndCtl, UINT codeNotify)
{

	switch(id) {
		case IDM_OPEN:
		OpenGenomeDialog(hwnd);
		break;

		case IDM_COMPARE:
		long l;
		struct sGenome *genome;
		char buffer[255];
		genome = (void*)GetWindowLong(hwnd, GWL_USERDATA);
		l = GenomeSimilarity(genome, genome->next);
		sprintf(buffer, "Long: %i", l);
		MessageBox(hwnd, buffer, "Compare", 0);
		break;

		case IDM_EXIT:
		PostMessage(hwnd,WM_CLOSE,0,0);
		break;
	}
}

LRESULT CALLBACK MainWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	switch (msg) {
	case WM_COMMAND:
		HANDLE_WM_COMMAND(hwnd,wParam,lParam,MainWndProc_OnCommand);
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hwnd,msg,wParam,lParam);
	}
	return 0;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, INT nCmdShow)
{
	MSG msg;
	HANDLE hAccelTable;

	hInst = hInstance;
	if (!InitApplication())
		return 0;
	hAccelTable = LoadAccelerators(hInst,MAKEINTRESOURCE(IDACCEL));
	if ((hwndMain = CreatesnplingWndClassWnd()) == (HWND)0)
		return 0;
	ShowWindow(hwndMain,SW_SHOW);
	while (GetMessage(&msg,NULL,0,0)) {
		if (!TranslateAccelerator(msg.hwnd,hAccelTable,&msg)) {
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}
	return msg.wParam;
}

int MessageSnp(HWND hwnd, struct sSnpEntry *snp)
{
	char buffer[256];

	sprintf(buffer, "Rsid: %s%i\r\nChromosome: %i\r\nPosition: %i\r\nGenotype: %c%c", snp->rstype, snp->rsnumber, snp->chrom, snp->position, snp->genotype[0], snp->genotype[1]);
	MessageBox(hwnd, buffer, "SNP Info",0);

	return 0;

}

int MessageCount(HWND hwnd, struct sCounts *count)
{
	char buffer[256];

	sprintf(buffer, "Total: %i\r\nHomo: %i (%i%%)\r\nHetero: %i", count->total,count->homozygous,100* count->homozygous/count->total,0);
	MessageBox(hwnd, buffer, "SNP Info",0);

	return 0;

}

int MessagePoint(HWND hwnd, POINT *p, char *heading)
{
	char buffer[256];

	sprintf(buffer, "X: %i, Y: %i", p->x, p->y);
	MessageBox(hwnd, buffer, heading,0);
	return 0;
}


int SnpFromLine(struct sSnpEntry *snp, char *line)
{
	char *result = NULL;

	result = strtok(line, "\t");
	if (!result) return 0;
	if ((result[0]=='r') && (result[1]=='s'))	{
		snp->rsnumber=strtol(result+2,NULL,10);
		snp->rstype[0]='r';snp->rstype[1]='s';snp->rstype[2]=0;
	}

	result = strtok(NULL, "\t");
	if (!result) return 0;
	if (result[0]=='X')
		snp->chrom = 23;
	else if (result[0]=='Y')
		snp->chrom = 24;
	else if (result[0]=='M')
		snp->chrom = 25;
	else
		snp->chrom = strtol(result, NULL, 10);

	result = strtok(NULL, "\t");
	if (!result) return 0;
	snp->position = strtol(result, NULL, 10);

	result = strtok(NULL, "\t");
	if (!result) return 0;
	snp->genotype[0] = result[0];
	snp->genotype[1] = result[1];

	return 1;
}

int AddLineToComments(struct sGenome *genome, char *line)
{
	int lenline;
	int lencomments;

	line++;	//skip the #
	if (line[0] == ' ')
		line++;

	lencomments=strlen(genome->comments);
	lenline=strlen(line);

	if (lenline<1)
		return 0;

	if (lencomments+lenline>MAXCOMMENTS)	{
		return 0;
	}



	memcpy(genome->comments+lencomments, line, lenline);

	return 1;
}

LRESULT CALLBACK GenomeInfoProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	struct sGenome *genome;
	POINT	mousePoint;
	POINT	newPoint;
	RECT	parentRect;
	POINT	parentPoint;

	switch (msg) {
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	case WM_PAINT:
		GenomeInfoPaint(hwnd);
		break;
	case WM_ERASEBKGND:	//we'll do the erasing - to prevent flicker
		return 1;
	case WM_LBUTTONDOWN:
		genome = (void*)GetWindowLong(hwnd,GWL_USERDATA);
		BringWindowToTop(hwnd);
		if (SetFocus(hwnd)!=hwnd)
				InvalidateRect(hwnd, NULL, FALSE);
		if (!genome->windowInfo.mouseBeenPushed)        {
				SetCapture(hwnd);
				GetClientRect(GetParent(hwnd), &parentRect);
				parentPoint.x=parentRect.left;
				parentPoint.y=parentRect.top;
				ClientToScreen(GetParent(hwnd), &parentPoint);
				genome->windowInfo.mouseBeenPushed=1;
				genome->windowInfo.initialClickPoint.x=GET_X_LPARAM(lParam)+parentPoint.x;
				genome->windowInfo.initialClickPoint.y=GET_Y_LPARAM(lParam)+parentPoint.y;
				GetWindowRect(hwnd, &genome->windowInfo.originalLocationRect);
				genome->windowInfo.originalWidth = genome->windowInfo.originalLocationRect.right-genome->windowInfo.originalLocationRect.left;
				genome->windowInfo.originalHeight= genome->windowInfo.originalLocationRect.bottom-genome->windowInfo.originalLocationRect.top;
		}
		break;
		case WM_MOUSEMOVE:
			genome = (void*)GetWindowLong(hwnd,GWL_USERDATA);
			if (genome->windowInfo.mouseBeenPushed) {
				mousePoint.x = GET_X_LPARAM(lParam) - genome->windowInfo.initialClickPoint.x;
				mousePoint.y = GET_Y_LPARAM(lParam) - genome->windowInfo.initialClickPoint.y;
				if (mousePoint.x || mousePoint.y)       {
					newPoint.x=genome->windowInfo.originalLocationRect.left+mousePoint.x;
					newPoint.y=genome->windowInfo.originalLocationRect.top+mousePoint.y;
					if (newPoint.x<0) newPoint.x=0;
					if (newPoint.y<0) newPoint.y=0;
					GetClientRect(GetParent(hwnd), &parentRect);
					if (newPoint.x+genome->windowInfo.originalWidth > parentRect.right)
						newPoint.x=parentRect.right-genome->windowInfo.originalWidth;
					if (newPoint.y+genome->windowInfo.originalHeight > parentRect.bottom)
						newPoint.y=parentRect.bottom-genome->windowInfo.originalHeight;

					MoveWindow(hwnd,
							newPoint.x,
							newPoint.y,
							genome->windowInfo.originalWidth,
							genome->windowInfo.originalHeight,
							TRUE);
					GetWindowRect(hwnd, &genome->windowInfo.originalLocationRect);
				}
			}
		break;

	case WM_LBUTTONUP:
		genome = (void*)GetWindowLong(hwnd,GWL_USERDATA);
		if (genome->windowInfo.mouseBeenPushed)	{
			genome->windowInfo.mouseBeenPushed=0;
			ReleaseCapture();
		}
		break;
	case WM_CAPTURECHANGED:
		genome = (void*)GetWindowLong(hwnd,GWL_USERDATA);
		genome->windowInfo.mouseBeenPushed = 0;
		break;
	default:
		return DefWindowProc(hwnd,msg,wParam,lParam);
	}
	return 0;
}

int GenomeInfoPaint(HWND hwnd)
{
	struct sGenome *genome;

	HDC hdc;
	PAINTSTRUCT ps;
	RECT	clientRect;
	char	textbuffer[256];
	int		textlen;

	RECT	outputRect;
	TEXTMETRIC textMetric;
	int heightFont;
	int x,y;
	int isXX=0;
	int totalbases;

	genome=(void*)GetWindowLong(hwnd, GWL_USERDATA);

	GetClientRect(hwnd, &clientRect);
	hdc=BeginPaint(hwnd, &ps);



	GetTextMetrics(hdc, &textMetric);
	heightFont= textMetric.tmHeight;

	sprintf(textbuffer,"%s", genome->filename);
	textlen=strlen(textbuffer);

	y=clientRect.top;
	outputRect.left=clientRect.left; 	outputRect.right=clientRect.right;
	outputRect.top=y;
	y+=heightFont+1;
	outputRect.bottom=y;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print filename

	sprintf(textbuffer, "%i SNPs", genome->count.total);
	textlen=strlen(textbuffer);
	outputRect.top=y;
	y+=heightFont+1;
	outputRect.bottom=y;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print total snps

	sprintf(textbuffer, "%i autosomal", genome->count.autosomal);
	textlen=strlen(textbuffer);
	outputRect.top=y;
	y+=heightFont+1;
	outputRect.bottom=y;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print autosomal snps

	if (genome->countXX*10 > genome->chromosome[23].count.total)	//if xxs are more than 10%
		isXX=1;
	sprintf(textbuffer, "%i,X%s%s%s", 45+genome->containsY+isXX, isXX ? "X":"", genome->containsY ? "Y":"", genome->containsM ? "+m":"");
	textlen=strlen(textbuffer);
	outputRect.top=y;
	y+=heightFont+1;
	outputRect.bottom=y;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print 46,XY etc


	totalbases=genome->count.adenine+ genome->count.thymine+ genome->count.cytosine+ genome->count.guanine;
	if (totalbases>0)	{
		sprintf(textbuffer, "A:%i T:%i C:%i G:%i", genome->count.adenine*100/totalbases, genome->count.thymine*100/totalbases, genome->count.cytosine*100/totalbases, genome->count.guanine*100/totalbases);
		textlen=strlen(textbuffer);
		outputRect.top=y;
		y+=heightFont+1;
		outputRect.bottom=y;
		ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print autosomal snps
	}

	//Fill rest
	outputRect.top=y;
	outputRect.bottom=clientRect.bottom;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, NULL, 0, NULL);	//fill in bottom



	EndPaint(hwnd, &ps);

	return 0;
}

long GenomeSimilarity(struct sGenome *g1, struct sGenome *g2)
{

	struct sSnpEntry *e1;
	struct sSnpEntry *e2;
	struct sChunk *c1;
	struct sChunk *c2;

	int i;
	long l=0;

	c1 = g1->firstChunk;
	c2 = g2->firstChunk;
	while (c1)	{
		for (i=0;i<SNPSINCHUNK;i++)	{
			e1 = &c1->entry[i];
			e2 = &c2->entry[i];
			if ( (e1->genotype[0] == e2->genotype[0]) && (IsProperBase(e2->genotype[0])) )
				l++;
			if ( (e2->genotype[1] == e2->genotype[1]) && (IsProperBase(e2->genotype[1])) )
				l++;

			//MessageSnp(0, e1);
			//MessageSnp(0, e2);
		}

		c1=c1->next;
		c2=c2->next;
	}


	return l;
}

int IsProperBase(char b)
{
	if (b=='A' || b=='T' || b=='C' || b=='G')
		return 1;
	return 0;

}
