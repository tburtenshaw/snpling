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

#define SEX_UNKNOWN 0
#define SEX_MALE 1
#define SEX_FEMALE 2

//Global variables
HINSTANCE hInst;		// Instance handle
HWND hwndMain;		//Main window handle
HWND hwndTree;		//Displays the tree
HWND hwndInfobox;	//Has information
HWND hwndChromobox;	//Displays the chromosome

//Window params
int	heightTreeWnd;


//The last x,y coords we used.
int	lastX;
int	lastY;

//Struct prototypes
struct sSnpEntry;
struct sChunk;	//contains a fix-length array of snps
struct sChromosome;
struct sGenome;
struct sCounts;
struct sGenomeWindowInfo;
struct sPerson;

LRESULT CALLBACK MainWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam);
LRESULT CALLBACK TreeWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam);
LRESULT CALLBACK InfoboxWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam);

//GUI
int TreePaint(HWND hwnd);
int InfoboxPaint(HWND hwnd, struct sPerson *p);
void * FindPersonFromCoords(int x, int y);

//Loading
int SnpFromLine(struct sSnpEntry *snp, char *line);
int AddLineToComments(struct sGenome *genome, char *line);
int OpenGenomeDialog(HWND hwnd);
int LoadGenomeFile(HWND hwnd, char *filename);
int CountEntry(struct sGenome *genome, struct sSnpEntry *snp);
char * FilenameStart(char *s);
int IsolateDate(struct sGenome *genome);
void * GetPhenotypeNameAndSexFromGenotype(struct sPerson *p);
int	ArrangeInitialLocation(struct sPerson *p);

//Calculation
int IsProperBase(char b);	//is it ATC or G? or nothing
long GenomeSimilarity(struct sGenome *g1, struct sGenome *g2);

//Debugging
int MessageSnp(HWND hwnd, struct sSnpEntry *snp);
int MessageCount(HWND hwnd, struct sCounts *count);
int MessagePoint(HWND hwnd, POINT *p, char *heading);

//Linked list
void *CreatePerson(void);	//returns a pointer to a malloc'd sPerson

struct sGenomeWindowInfo
{
	HWND	hwnd;
	POINT	initialClickPoint;	//used for moving
	RECT	originalLocationRect;
	int		originalHeight;
	int		originalWidth;
	int		mouseBeenPushed;	//have we clicked the mouse
};

struct sDisplay				//how the person is displayed on the tree
{
	//position
	int	x;
	int	y;
	int	height;
	int	width;
	int layer;

	//appearance
	int	selected;

};

struct sSnpEntry
{
	char rstype[4];
	long rsnumber;
	int chrom;
	long position;
	char genotype[2];
	char phasing;	//not sure how to do this yet
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
	struct sChunk *pChunk;	//pointers to the chunk the specific chromosome is first found
	int		startEntry;
	struct sCounts count;
};

struct sGenome
{
	struct sChunk *firstChunk;
	struct sChromosome chromosome[26];			//info about chromosomes, we use 1-25 (0's wasted)

	struct sCounts count;

	char	filename[MAX_PATH];
	char	*justfilename;	//points to the above without path
	char	comments[MAXCOMMENTS];	//This contains the comments from the file
	long	countXX;		//Does it have two alleles on the x-chromosome (even males have some on 23andme)
	int		containsY;		//Does it have a Y chromosome?
	int		containsM;		//Does the data have mitochondrial DNA?
	int		chromosomeCount;

};

struct sPhenotype
{
	int		sex;
	int		age;
	char	firstname[64];
	char	lastname[64];
};


struct sPerson
{
	struct sPhenotype *phenotype;
	struct sGenome *genome;
	struct sRelationships *relationships;
	struct sDisplay *display;					//Used to store how and where the representation appears

	struct	sPerson *next;

};

struct sPerson *firstPerson;	//global var
struct sPerson *selectedPerson;

int OpenGenomeDialog(HWND hwnd)	//returns 0 if fails
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

	r=GetOpenFileName(&ofn);	//Gets the file name

	if (!r)
		return 0;				//If there's no file name, then return.

	filenamepart=ofn.lpstrFile+strlen(ofn.lpstrFile)+1;	//this essentially checks whether there's a \0\0
	if (!(*filenamepart))	{			//we've only selected one file
		LoadGenomeFile(hwnd, ofn.lpstrFile);
	}
	else	{							//we've got multiple files selected
		while (*filenamepart)	{		//so go through them all
			sprintf(fullpath, "%s\\%s", ofn.lpstrFile, filenamepart);
			LoadGenomeFile(hwnd, fullpath);
			filenamepart+=strlen(filenamepart)+1;
		}
	}

	return  1;	//success
}

void *CreatePerson()
{
	struct sPerson *p;

	if (!firstPerson)	{
		firstPerson = malloc(sizeof(struct sPerson));
		memset(firstPerson,0,sizeof(struct sPerson));
		//that prev memset of course means firstPerson->next is set to NULL;
		return firstPerson;
	}

	p=firstPerson;
	while	(p->next)	{
		p=p->next;
	}
	p->next=malloc(sizeof(struct sPerson));
	memset(p->next,0,sizeof(struct sPerson));

	return p->next;
};

int	DeletePerson(struct sPerson *personToDelete)
{
	struct sPerson	*p;
	struct sPerson	*prevP;

	prevP=NULL;
	p=firstPerson;

	while	(p)	{
		if	(p == personToDelete)	{
			//Removing and freeing the phenotype is easy
			if (p->phenotype)	{
				free(p->phenotype);
			}

			//Removing the genotype means removing all the snp chunks first
			struct	sChunk	*c;
			struct	sChunk	*cnext;
			if	(p->genome)	{
				c=p->genome->firstChunk;	//go through all the chunks
				while	(c)	{
					cnext=c->next;
					free(c);
					c=cnext;
				}
				free(p->genome);			//then free genome

			}

			//Then remove the display info
			if (p->display)	{
				free(p->display);
			}

			//if selected, hovered, or otherwise name then make them point to previous
			if (selectedPerson == p)	{
				selectedPerson = prevP;
			}


			//Now we can remove the node itself
			if (!prevP)	{//if this is the firstPerson
				free(firstPerson);
				firstPerson=p->next;
				return 1;
			}
			prevP->next = p->next;	//point the next in the previous node to the next one.
			free(p);
			return 1;
		}

		prevP=p;
		p=p->next;

	}


	return 0;
}

int LoadGenomeFile(HWND hwnd, char *filename)
{
	FILE * fDNA;	//standard file pointer
	struct sChunk *snpChunk;
	struct sChunk *prevChunk;
	//struct sGenome *genome;
	struct sGenome *othergenome;

	struct sPerson *person;
	int		arrayPos;
	int		prevChromosomeNumber;
	POINT	windowloc;

	fDNA = fopen(filename, "r");	//open the filename, readonly and ascii
	if (!fDNA)
		return 0;

	person=CreatePerson();

	person->genome = malloc(sizeof(struct sGenome));	//Get memory for new genome
	memset(person->genome, 0, sizeof(struct sGenome));	//Set it all to zero initially
	//Set some characteristics of the genome
	strcpy(person->genome->filename, filename);	//Filename


	snpChunk = malloc(sizeof(struct sChunk));
	person->genome->firstChunk=snpChunk;
	arrayPos=0;

	char line[256];		//reads each line of the txt file, up to the first 256 chars
	prevChromosomeNumber=0;
	while(fgets(line, sizeof(line), fDNA))     {
		if (line[0]=='#')	{//handle comments line
			AddLineToComments(person->genome, line);
		}
		else if (SnpFromLine(&snpChunk->entry[arrayPos], line)) {
			if (snpChunk->entry[arrayPos].chrom > prevChromosomeNumber)
				person->genome->chromosome[snpChunk->entry[arrayPos].chrom].pChunk=snpChunk;

			CountEntry(person->genome, &snpChunk->entry[arrayPos]);
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

	person->genome->justfilename=FilenameStart(person->genome->filename);
	IsolateDate(person->genome);		//Get the date loaded
	GetPhenotypeNameAndSexFromGenotype(person);		//Populate some
	ArrangeInitialLocation(person);

	InvalidateRect(hwndTree,NULL, 0);
	InvalidateRect(hwndInfobox,NULL, 0);

	return 0;
}

int	ArrangeInitialLocation(struct sPerson *p)
{
	RECT clientRect;

	if (p->display)
		return 0;	//return 0 if there's already a display

	p->display= malloc(sizeof(struct sDisplay));
	memset(p->display,0,sizeof(struct sDisplay));


	GetClientRect(hwndTree, &clientRect);

	p->display->x = lastX+40;

	if (p->display->x+40>clientRect.right)	{
		p->display->x=0;
		lastY+=40;
	}

	p->display->y = lastY;


	lastX=p->display->x;
	lastY=p->display->y;

	return 1;
}

void	*GetPhenotypeNameAndSexFromGenotype(struct sPerson *p)
{
	if (!p->phenotype)	{
		p->phenotype=malloc(sizeof(struct sPhenotype));
		memset(p->phenotype, 0, sizeof(struct sPhenotype));
	}

	if	(p->genome->containsY)
		p->phenotype->sex=SEX_MALE;
	else
		p->phenotype->sex=SEX_FEMALE;


	if (!p->genome->justfilename)
		return p->phenotype;

	if (strncmp(p->genome->justfilename, "genome_", 7))	{
		p->phenotype->lastname[0]=0;
		return p->phenotype;
	}

	strcpy(p->phenotype->lastname, p->genome->justfilename+7);



	return p->phenotype;
}

static BOOL InitApplication(void)
{

	//Initalise some variables
	firstPerson=NULL;
	selectedPerson=NULL;
	heightTreeWnd=300;

	lastX=0;
	lastY=0;

	//Register the Window Classes we use
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

	//Tree window;
	memset(&wc,0,sizeof(WNDCLASS));
	wc.style = CS_DBLCLKS;
	wc.lpfnWndProc = (WNDPROC)TreeWndProc;
	wc.hInstance = hInst;
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOWFRAME);
	wc.lpszClassName = "TreeWndClass";
	wc.lpszMenuName = MAKEINTRESOURCE(IDMAINMENU);
	wc.hCursor = LoadCursor(NULL,IDC_ARROW);
	if (!RegisterClass(&wc))
		return 0;

	//Infobox window;
	memset(&wc,0,sizeof(WNDCLASS));
	wc.style = CS_DBLCLKS;
	wc.lpfnWndProc = (WNDPROC)InfoboxWndProc;
	wc.hInstance = hInst;
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOWFRAME);
	wc.lpszClassName = "InfoboxWndClass";
	wc.lpszMenuName = MAKEINTRESOURCE(IDMAINMENU);
	wc.hCursor = LoadCursor(NULL,IDC_ARROW);
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
		//not really compare, just tries deleting
		DeletePerson(selectedPerson);
		InvalidateRect(hwndTree,NULL, 0);
		InvalidateRect(hwndInfobox,NULL, 0);
		break;

		case IDM_EXIT:
		PostMessage(hwnd,WM_CLOSE,0,0);
		break;
	}
}

LRESULT CALLBACK MainWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	switch (msg) {
		case WM_CREATE:	//when this window is made, make its main children
			RECT clientRect;
			GetClientRect(hwnd, &clientRect);
			hwndTree=CreateWindow("TreeWndClass","Tree", WS_CHILD|WS_VISIBLE,0,0,clientRect.right,300,
			hwnd, NULL, NULL, NULL);

			hwndInfobox=CreateWindow("InfoboxWndClass", "Infobox", WS_CHILD|WS_VISIBLE, 0, heightTreeWnd, clientRect.right-123,  clientRect.bottom-heightTreeWnd,
			hwnd, NULL, NULL, NULL);

			break;
		case WM_SIZE:
			GetClientRect(hwnd, &clientRect);
			MoveWindow(hwndTree, 0,0,clientRect.right, heightTreeWnd,1);
			break;
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

LRESULT CALLBACK InfoboxWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	switch (msg) {
		case WM_PAINT:
			InfoboxPaint(hwnd, selectedPerson);
			break;
		default:
			return DefWindowProc(hwnd,msg,wParam,lParam);
	}

	return 0;
}

int InfoboxPaint(HWND hwnd, struct sPerson *p)
{
	HDC hdc;
	PAINTSTRUCT  ps;
	int	y;
	RECT clientRect;
	RECT outputRect;
	TEXTMETRIC textMetric;
	char textbuffer[255];
	int	textlen;
	int heightFont;

	GetClientRect(hwnd, &clientRect);

	if (!p)	{//if there's noone selected
		sprintf(textbuffer, "none selected");
	}
	else	{
		sprintf(textbuffer,"%s", p->genome->filename);
	}

	textlen=strlen(textbuffer);

	hdc=BeginPaint(hwnd, &ps);
	GetTextMetrics(hdc, &textMetric);
	heightFont= textMetric.tmHeight;
	y=clientRect.top;
	outputRect.left=clientRect.left; 	outputRect.right=clientRect.right;
	outputRect.top=y;
	outputRect.bottom=y+heightFont;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print filename

	EndPaint(hwnd, &ps);

	return 1;
}


LRESULT CALLBACK TreeWndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	switch (msg) {
		case WM_PAINT:
			return TreePaint(hwnd);
			break;
		case WM_LBUTTONDOWN:
			selectedPerson = FindPersonFromCoords(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
			InvalidateRect(hwndInfobox, NULL, 0);
			InvalidateRect(hwndTree, NULL, 0);
			break;

		default:
			return DefWindowProc(hwnd,msg,wParam,lParam);
	}

	return 0;
}

void * FindPersonFromCoords(int x, int y)
{
	struct sPerson *p;

	p=firstPerson;

	while (p)	{
		if ((x >= p->display->x) &&
			(x <=p->display->x+40) &&
			(y >= p->display->y) &&
			(y <= p->display->y+40))	{

			return p;
		}
		p=p->next;
	}

	return NULL;
}

int TreePaint(HWND hwnd)	{

	HDC          hdcMem;
	HBITMAP      hbmMem;
	HANDLE       hOld;

	PAINTSTRUCT  ps;
	HDC          hdc;

	RECT	clientRect;
	HPEN		hPen;
	HPEN		hPenSelected;
	HBRUSH		hBrush;

	int x,y;

	struct	sPerson	*p;

	hdc = BeginPaint(hwnd, &ps);

    // Create an off-screen DC for double-buffering (should move this so it's not done every paint)
    hdcMem = CreateCompatibleDC(hdc);
    GetClientRect(hwnd, &clientRect);
	hbmMem = CreateCompatibleBitmap(hdc, clientRect.right, clientRect.bottom);

    hOld   = SelectObject(hdcMem, hbmMem);

    // Draw into hdcMem here

	hPen = CreatePen(PS_SOLID, 0, RGB(235,140,0));
	hPenSelected = CreatePen(PS_SOLID, 0, RGB(255,240,0));
	hBrush = CreateSolidBrush(RGB(220,120,20));

	SelectObject(hdcMem, hBrush);
	p=firstPerson;
	while (p)	{
		x = p->display->x;
		y = p->display->y;

		if (p == selectedPerson)	{
			SelectObject(hdcMem, hPenSelected);
		}	else	{
			SelectObject(hdcMem, hPen);
		}

		if	(p->phenotype->sex==SEX_MALE)
			Rectangle(hdcMem, x, y, x+30, y+30);
		if	(p->phenotype->sex==SEX_FEMALE)
			Ellipse(hdcMem, x, y, x+30, y+30);
		p=p->next;
	}



	// Transfer the off-screen DC to the screen
    BitBlt(hdc, 0, 0, clientRect.right, clientRect.bottom, hdcMem, 0, 0, SRCCOPY);

	//Clean up objects
	DeleteObject(hBrush);
	DeleteObject(hPen);
	DeleteObject(hPenSelected);

    // Free-up the off-screen DC
    SelectObject(hdcMem, hOld);
    DeleteObject(hbmMem);
    DeleteDC    (hdcMem);
	EndPaint(hwnd, &ps);

	return 1;

}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, INT nCmdShow)
{
	MSG msg;
	HANDLE hAccelTable;

	hInst = hInstance;
	if (!InitApplication())	//registers the window classes
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

	result = strtok(line, "\t");	//gets the first "token" up until the tab, it'll be the SNP name
	if (!result) return 0;
	if ((result[0]=='r') && (result[1]=='s'))	{		//if it's an rs
		snp->rsnumber=strtol(result+2,NULL,10);
		snp->rstype[0]='r';snp->rstype[1]='s';snp->rstype[2]=0;
	}
	else	{
		if (result[0]='i')	{
		snp->rsnumber=strtol(result+1,NULL,10);
		snp->rstype[0]='i';snp->rstype[0]=0;
		}
	}


	result = strtok(NULL, "\t");	//then the next token (it's next, because we use NULL) should be the chromosome
	if (!result) return 0;
	if (result[0]=='X')
		snp->chrom = 23;
	else if (result[0]=='Y')
		snp->chrom = 24;
	else if (result[0]=='M')
		snp->chrom = 25;
	else
		snp->chrom = strtol(result, NULL, 10);	//10 is for base-10 (decimal)

	result = strtok(NULL, "\t");	//then position
	if (!result) return 0;
	snp->position = strtol(result, NULL, 10);

	result = strtok(NULL, "\t");	//then genotype
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
	if (line[0] == ' ')	//skip leading space
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


/*
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
		sprintf(textbuffer, "A:%0.1f%% T:%0.1f%% C:%0.1f%% G:%0.1f%%", ((float)(genome->count.adenine*1000/totalbases)/10), (float)genome->count.thymine*100/totalbases, (float)genome->count.cytosine*100/totalbases, (float)genome->count.guanine*100/totalbases);
		textlen=strlen(textbuffer);
		outputRect.top=y;
		y+=heightFont+1;
		outputRect.bottom=y;
		ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print autosomal snps
	}

	sprintf(textbuffer, "%s", genome->owner);
	textlen=strlen(textbuffer);
	outputRect.top=y;
	y+=heightFont+1;
	outputRect.bottom=y;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, textbuffer, textlen, NULL);	//print autosomal snps


	//Fill rest
	outputRect.top=y;
	outputRect.bottom=clientRect.bottom;
	ExtTextOut(hdc, 0,outputRect.top,ETO_OPAQUE, &outputRect, NULL, 0, NULL);	//fill in bottom



	EndPaint(hwnd, &ps);

	return 0;
}
*/


long GenomeSimilarity(struct sGenome *g1, struct sGenome *g2)
{

	struct sSnpEntry *e1;
	struct sSnpEntry *e2;
	struct sChunk *c1;
	struct sChunk *c2;

	int i1;
	int i2;
	long l=0;

	c1 = g1->firstChunk;
	c2 = g2->firstChunk;
	while (c1)	{
		for (i1=0, i2=0;(i1<SNPSINCHUNK);i1++)	{
			e1 = &c1->entry[i1];
			e2 = &c2->entry[i2];
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

int IsolateDate(struct sGenome *genome)
{


	return 0;
}

char * FilenameStart(char *s)
{

	if (s)	{
		s=strrchr(s, '\\')+1;
	}

	return s;
}

int CountEntry(struct sGenome *genome, struct sSnpEntry *snp)
{
	char	localGenotype[2];	//make it shorter!
	int		localChromosome;


	int		properBase0;
	int		properBase1;
	int		containsRealBase;	//it's not a no-call or a placeholder
	int		bothRealBases;		//AG not A-, or T



	localChromosome=snp->chrom;
	localGenotype[0]=snp->genotype[0];
	localGenotype[1]=snp->genotype[1];

	genome->count.total++;
	genome->chromosome[localChromosome].count.total++;

	properBase0 = IsProperBase(localGenotype[0]);
	properBase1 = IsProperBase(localGenotype[1]);

	containsRealBase = (properBase0 || properBase1);	//is one of these a proper base?
	bothRealBases = (properBase0 && properBase1);	//is one of these a proper base?

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
			if 	(bothRealBases) {
				genome->countXX++;
			}
		}
	}

	return 0;
}
