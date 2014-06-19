/*
	Copyright (c) 2005, Arun S Konagurthu, The University of Melbourne.
	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, are permitted provided 
	that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice, this list of conditions and the 
	  following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and 
	  the following disclaimer in the documentation and/or other materials provided with the distribution.
    	* Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or 
	  promote products derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
	WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
	PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
	ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
	OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
	IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/   


PAIRWISE(char SEQs[20][2000],int num, char *scoremat, int GAPCOST)
{
	int i,j,k;
	FILE *fpo1;
	int SigmaProjHeur=0,SigmaPairwise=0,ZZZZ=0,UminusL=0;
	char c[10];
	char cmd[300];
	fpo1=fopen("./work/LminusU.dat","r");
	fscanf(fpo1,"%d",&SigmaProjHeur);
	fclose(fpo1);
	for(i=0;i<num;i++)
		for(j=i+1;j<num;j++)
		{ 	strcpy(cmd,"./work/C.input"); c[0]=i+65; c[1]='\0';strcat(cmd,c); c[0]=j+65; c[1]='\0'; strcat(cmd,c); 
		//	printf("%s\n",cmd);
	
			fpo1=fopen(cmd,"w");
		        for(k=1;k<strlen(SEQs[i]);k++)
        			fprintf(fpo1,"%c",SEQs[i][k]);
		        fprintf(fpo1,"*\n");
		        for(k=1;k<strlen(SEQs[j]);k++)
        			fprintf(fpo1,"%c",SEQs[j][k]);
		        fprintf(fpo1,"*\n");
		        fclose(fpo1);
		}

	for(i=0;i<num;i++)
		for(j=i+1;j<num;j++)
		{
		        strcpy(cmd,"./bin/Pairwise");
		        strcat(cmd," ./work/C.input");
			c[0]=i+65; c[1]='\0'; strcat(cmd,c);
			c[0]=j+65; c[1]='\0'; strcat(cmd,c);
		        //strcat(cmd," onlyIntScores_2D_allpositives.output");
		        strcat(cmd," ");
		        strcat(cmd, scoremat );
		        strcat(cmd," ");
		        strcat(cmd," ./work/C.output");
			c[0]=i+65; c[1]='\0'; strcat(cmd,c);
			c[0]=j+65; c[1]='\0'; strcat(cmd,c);
		        strcat(cmd," ./work/C_explore_");
			c[0]=i+65; c[1]='\0'; strcat(cmd,c);
			c[0]=j+65; c[1]='\0'; strcat(cmd,c);
		        strcat(cmd,".dat");
		        strcat(cmd," ./work/C_Rexplore_");
			c[0]=i+65; c[1]='\0'; strcat(cmd,c);
			c[0]=j+65; c[1]='\0'; strcat(cmd,c);
		        strcat(cmd,".dat");
			char buff[10] = "" ;
			sprintf(buff," %d ", GAPCOST) ;
		        strcat(cmd, buff);
	   	//	printf("%s\n",cmd);
		        system(cmd);
		}
	
	printf("Pairwise Optimal_2   :");
	for(i=0;i<num;i++)
		for(j=i+1;j<num;j++)
		{
		        strcpy(cmd,"./work/C_Rexplore_");
			c[0]=i+65; c[1]='\0'; strcat(cmd,c);
			c[0]=j+65; c[1]='\0'; strcat(cmd,c);
		        strcat(cmd,".dat");
			fpo1=fopen(cmd,"r");
			fscanf(fpo1,"%d",&ZZZZ);
			fclose(fpo1);
			SigmaPairwise+=ZZZZ;
			printf("%d ",ZZZZ);
		}
	printf("=%d\n",SigmaPairwise);
	UminusL=SigmaProjHeur-SigmaPairwise;
        printf("--->LminusU_2:%d\n",UminusL);
        fpo1=fopen("./work/LminusU_C.dat","w");
        fprintf(fpo1,"%d",UminusL);
        fclose(fpo1);
	//exit(0);
}
