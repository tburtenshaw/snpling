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

int SnpFromLine(struct sSnpEntry *snp, char *line);
int AddLineToComments(struct sGenome *genome, char *line);
int OpenGenomeDialog(HWND hwnd);
int LoadGenomeFile(HWND hwnd, char *filename);

//Debugging
int MessageSnp(HWND hwnd, struct sSnpEntry *snp);
int MessageCount(HWND hwnd, struct sCounts *count);
int MessagePoint(HWND hwnd, POINT *p, char *heading);

struct sGenomeWindowInfo
{
	HWND	hwnd;
	POINT	initialClickPoint;	//used for moving
	RECT	originalLocationRect;
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
	struct sChromosome chromosome[26];			//info about chromosomes, we use 1-25

	struct sCounts count;

	char	filename[MAX_PATH];
	char	owner[256];	//the user name associated (if we can pull it from a file)
	char	comments[MAXCOMMENTS];	//This contains the comments from the file
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
	struct sSnpEntry *snpEntry;
	struct sChunk *snpChunk;
	struct sChunk *prevChunk;
	struct sGenome *genome;
	struct sGenome *othergenome;
	int		arrayPos;
	int		r;
	POINT	windowloc;


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
			genome->count.total++;
			genome->chromosome[snpChunk->entry[arrayPos].chrom].count.total++;
			if (snpChunk->entry[arrayPos].genotype[0] == snpChunk->entry[arrayPos].genotype[1])	{
				genome->count.homozygous++;
				genome->chromosome[snpChunk->entry[arrayPos].chrom].count.homozygous++;
			}
			if (snpChunk->entry[arrayPos].chrom<23)	{
				genome->count.autosomal++;
				genome->chromosome[snpChunk->entry[arrayPos].chrom].count.autosomal++;
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

	genome->windowInfo.hwnd = CreateWindow("GenomeInfoBoxClass", genome->filename, WS_BORDER|WS_CHILD|WS_VISIBLE, 32,windowloc.y, 180,120, hwnd, NULL, hInst, NULL);
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
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
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
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+2);
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
		if (!genome->windowInfo.mouseBeenPushed)	{
			SetCapture(hwnd);
			GetClientRect(GetParent(hwnd), &parentRect);
			parentPoint.x=parentRect.left;
			parentPoint.y=parentRect.top;
			ClientToScreen(GetParent(hwnd), &parentPoint);
			genome->windowInfo.mouseBeenPushed=1;
			genome->windowInfo.initialClickPoint.x=GET_X_LPARAM(lParam)+parentPoint.x;
			genome->windowInfo.initialClickPoint.y=GET_Y_LPARAM(lParam)+parentPoint.y;
			GetWindowRect(hwnd, &genome->windowInfo.originalLocationRect);
		}
		return 0;
		break;
	case WM_MOUSEMOVE:
		genome = (void*)GetWindowLong(hwnd,GWL_USERDATA);
		if (genome->windowInfo.mouseBeenPushed)	{

			mousePoint.x = GET_X_LPARAM(lParam) - genome->windowInfo.initialClickPoint.x;
			mousePoint.y = GET_Y_LPARAM(lParam) - genome->windowInfo.initialClickPoint.y;

			if (mousePoint.x || mousePoint.y)	{
				MoveWindow(hwnd,
					genome->windowInfo.originalLocationRect.left+mousePoint.x,
					genome->windowInfo.originalLocationRect.top+mousePoint.y,
					genome->windowInfo.originalLocationRect.right-genome->windowInfo.originalLocationRect.left,
					genome->windowInfo.originalLocationRect.bottom-genome->windowInfo.originalLocationRect.top,
					TRUE);
				GetWindowRect(hwnd, &genome->windowInfo.originalLocationRect);
				}
			}
		return 0;
		break;
	case WM_LBUTTONUP:
		genome = (void*)GetWindowLong(hwnd,GWL_USERDATA);
		if (genome->windowInfo.mouseBeenPushed)	{
			genome->windowInfo.mouseBeenPushed=0;
			ReleaseCapture();
		}
		return 0;
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
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print total snps

	EndPaint(hwnd, &ps);

	return 0;
}
