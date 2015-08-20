#include "io.hpp"

//needed for reading bam
#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static inline uint8_t *skip_aux(uint8_t *s)
{
    int size = aux_type2size(*s); ++s; // skip type
    uint32_t n;
    switch (size) {
    case 'Z':
    case 'H':
        while (*s) ++s;
        return s + 1;
    case 'B':
        size = aux_type2size(*s); ++s;
        memcpy(&n, s, 4); s += 4;
        return s + size * n;
    case 0:
        abort();
        break;
    default:
        return s + size;
    }
}

uint8_t *bam_aux_get(uint8_t *s,const uint8_t *aux,const int32_t aux_size,const char tag[2])
{
    int y = tag[0]<<8 | tag[1];

    while (s < aux+aux_size) {
        int x = (int)s[0]<<8 | s[1];
        s += 2;
        if (x == y) return s;
        s = skip_aux(s);
    }
    return 0;
}

int32_t bam_aux2i(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'c') return (int32_t)*(int8_t*)s;
    else if (type == 'C') return (int32_t)*(uint8_t*)s;
    else if (type == 's') return (int32_t)*(int16_t*)s;
    else if (type == 'S') return (int32_t)*(uint16_t*)s;
    else if (type == 'i' || type == 'I') return *(int32_t*)s;
    else return 0;
}

//
void ReadPeakBedFile(const string Peakbedfile,const string Bamfile,vector<BedRegion> &peakbedregion_set)
{
	string sbuf,chr,s1;
	int ia,ib,i,j;

	//read header
	vector<string> bam_chrorder;

	gzFile fp = gzopen(Bamfile.c_str(), "r");
	if(fp==NULL) {cout<<"gzopen error "<<Bamfile<<endl;exit(0);}

	char magic[4];
	int32_t l_text,n_ref,l_name,l_ref;

	gzread(fp, magic, sizeof(char)*4);
	gzread(fp, &l_text, sizeof(int32_t));

	char text[l_text];
	gzread(fp, text, sizeof(char)*l_text);
	gzread(fp, &n_ref, sizeof(int32_t));

	for(int i=0;i<n_ref;i++)
	{
		gzread(fp, &l_name, sizeof(int32_t));
		char name[l_name];
		gzread(fp, name, sizeof(char)*l_name);
		gzread(fp, &l_ref, sizeof(int32_t));

		string stemp(name);
		bam_chrorder.push_back(stemp);
	}

	//
	ifstream is(Peakbedfile.c_str());
	if(!is) {cout<<"can not open: "<<Peakbedfile<<endl;exit(0);}

	map<string,vector<BedRegion> > chr2Bedregionset;
	do
	{
		chr.clear();
		is>>chr;
		if(chr.empty()) break;

		is>>ia>>ib;getline(is,sbuf);
		ia++;

		BedRegion mybedregion;
		mybedregion.chr=chr;
		mybedregion.start=ia;
		mybedregion.end=ib;

		if(chr2Bedregionset.find(chr)==chr2Bedregionset.end())
		{
			vector<BedRegion> vBtemp;
			vBtemp.push_back(mybedregion);
			chr2Bedregionset[chr]=vBtemp;
		}
		else chr2Bedregionset[chr].push_back(mybedregion);

	}while(!is.eof());
	is.close();

	//
	for(i=0;i<bam_chrorder.size();i++)
	{
		if(chr2Bedregionset.find(bam_chrorder[i])!=chr2Bedregionset.end())
		{
			vector<BedRegion> vBtemp(chr2Bedregionset[bam_chrorder[i]]);
			for(j=0;j<vBtemp.size();j++) peakbedregion_set.push_back(vBtemp[j]);
		}
	}
}

void GetRefSeq(const uint8_t *seq,const int32_t l_seq,const char *qual,string &cigar,string &MD_str,string &revisedseq,string &revisedbq,int &readend,string &refseq)
{
	int i,j,itemp,ia;

	for(i=0;i<l_seq;i++)
	{
		revisedseq.push_back(seq_nt16_str[bam_seqi(seq,i)]);
		revisedbq.push_back(qual[i]+33);
	}

	//for 'H'
	itemp=0;
	for(i=0;i<cigar.length();i++)
	{
		if(cigar[i]=='H' && i!=cigar.length()-1) {cigar=cigar.substr(i+1);i=0;}
		else if(cigar[i]=='H' && i==cigar.length()-1) {cigar=cigar.substr(0,itemp+1);}

		if(cigar[i]<'0' || cigar[i]>'9')
		{
			itemp=i;
			if(i==0) {cout<<"error1:"<<cigar<<endl;exit(0);}
		}
	}

	//for 'S'
	itemp=0;
	for(i=0;i<cigar.length();i++)
	{
		if(cigar[i]=='S' && i!=cigar.length()-1)
		{
			ia=atoi(cigar.substr(0,i).c_str());
			revisedseq=revisedseq.substr(ia);
			revisedbq=revisedbq.substr(ia);
			cigar=cigar.substr(i+1);
			i=0;
		}
		else if(cigar[i]=='S' && i==cigar.length()-1)
		{//cout<<cigar<<" ";
			ia=atoi(cigar.substr(itemp+1,cigar.length()-itemp-1).c_str());
			revisedseq=revisedseq.substr(0,revisedseq.length()-ia);
			revisedbq=revisedbq.substr(0,revisedbq.length()-ia);
			cigar=cigar.substr(0,itemp+1);//cout<<cigar<<endl;
		}

		if(cigar[i]<'0' || cigar[i]>'9')
		{
			itemp=i;
			if(i==0) {cout<<"error1:"<<cigar<<endl;exit(0);}
		}
	}

	//get readend
	int start=0;
	int end=0;
	int seq_index=0;
	for(i=0;i<cigar.length();i++)
	{
		if(cigar[i]<'0' || cigar[i]>'9')
		{
			end=i-1;
			string stemp=cigar.substr(start,end-start+1);

			itemp=atoi(stemp.c_str());
			start=i+1;

			if(!stemp.empty())
			{
				if(cigar[i]=='I')//raw reads longer than ref
				{
					;
				}
				else if(cigar[i]=='D')//raw reads shorter than ref
				{
					readend+=itemp;
				}
				else if(cigar[i]=='M')//match or mismatch
				{
					readend+=itemp;
				}
				else {cout<<"error:"<<cigar<<endl;exit(0);}
			}
		}
	}

	//
	if(MD_str[0]!='Z') {cout<<"wrong MD infor in reads: "<<MD_str<<endl;exit(0);}
	MD_str=MD_str.substr(1);//delete the first 'Z'

	refseq=revisedseq;
//cout<<cigar<<" "<<MD_str<<endl;
//cout<<"raw\n"<<revisedseq<<endl;

	//for refseq, delete I part
	start=0;
	end=0;
	seq_index=0;//index in refseq
	vector<int> Istart_index,Ilength;
	for(i=0;i<cigar.length();i++)
	{
		if(cigar[i]<'0' || cigar[i]>'9')
		{
			end=i-1;
			string stemp=cigar.substr(start,end-start+1);

			itemp=atoi(stemp.c_str());
			start=i+1;

			if(!stemp.empty())
			{
				if(cigar[i]=='I')//raw reads longer than ref
				{
					Istart_index.push_back(seq_index);
					Ilength.push_back(itemp);
					seq_index=seq_index+itemp;
				}
				else if(cigar[i]=='D')//raw reads shorter than ref
				{
					;
				}
				else if(cigar[i]=='M')//match or mismatch
				{
					seq_index=seq_index+itemp;
				}
				else {cout<<"error:"<<cigar<<endl;exit(0);}
			}
		}
	}
	if(!Istart_index.empty())
	{
		vector<string> vstemp;
		vstemp.push_back(refseq.substr(0,Istart_index[0]));
		if(Istart_index.size()>1)
		{
			for(i=0;i<Istart_index.size()-1;i++) vstemp.push_back(refseq.substr(Istart_index[i]+Ilength[i],Istart_index[i+1]-(Istart_index[i]+Ilength[i])));
		}
		vstemp.push_back(refseq.substr(Istart_index[Istart_index.size()-1]+Ilength[Istart_index.size()-1]));
		refseq.clear();
		for(i=0;i<vstemp.size();i++) refseq+=vstemp[i];
		vstemp.clear();
	}
//cout<<"delI\n"<<refseq<<endl;

	//add D part and change NT
	start=0;
	end=0;
	seq_index=0;
	for(i=0;i<MD_str.length();i++)
	{
		if(MD_str[i]=='^')
		{
			if(i==0) {cout<<"wrong MD flag: "<<MD_str<<endl;exit(0);}
			else
			{
				end=i-1;
				string stemp=MD_str.substr(start,end-start+1);

				if(!stemp.empty())
				{
					itemp=atoi(stemp.c_str());
					start=i+1;
					seq_index=seq_index+itemp;
				}
			}

			start=i+1;
			do
			{
				i++;
				if(MD_str[i]>='0' && MD_str[i]<='9')
				{
					end=i-1;
					string stemp=MD_str.substr(start,end-start+1);
					refseq=refseq.substr(0,seq_index)+stemp+refseq.substr(seq_index);
					seq_index=seq_index+stemp.size();
					start=i;
					break;
				}

			}while(1);
		}
		else if(MD_str[i]<'0' || MD_str[i]>'9')
		{
			end=i-1;
			string stemp=MD_str.substr(start,end-start+1);

			if(!stemp.empty())
			{
				itemp=atoi(stemp.c_str());
				refseq[seq_index+itemp]=MD_str[i];
				start=i+1;
				seq_index=seq_index+itemp+1;
			}
		}
	}
//cout<<"fina\n"<<refseq<<endl;
//cout<<refseq.size()<<" "<<readend<<endl<<endl;
}

void MyFree(vector<BamInfor> &PeakBamInfor)
{
	vector<BamInfor> vBamtemp;
	PeakBamInfor.swap(vBamtemp);
}

int GetReadLengthFromBamFile(const string Bamfile)
{
	//read bam header
	map<int,string> chrindex2name;

	gzFile fp = gzopen(Bamfile.c_str(), "r");
	if(fp==NULL) {cout<<"gzopen error "<<Bamfile<<endl;exit(0);}

	//read header
	char magic[4];
	int32_t l_text,n_ref,l_name,l_ref;

	gzread(fp, magic, sizeof(char)*4);
	gzread(fp, &l_text, sizeof(int32_t));

	char text[l_text];
	gzread(fp, text, sizeof(char)*l_text);
	gzread(fp, &n_ref, sizeof(int32_t));

	for(int i=0;i<n_ref;i++)
	{
		gzread(fp, &l_name, sizeof(int32_t));
		char name[l_name];
		gzread(fp, name, sizeof(char)*l_name);
		gzread(fp, &l_ref, sizeof(int32_t));

		string stemp(name);
		chrindex2name[i+1]=stemp;
	}

	//read bam reads
	int i=0;
	vector<int> read_length_set;
	do
	{
		i++;

		//read a read
		int32_t block_size,refID,readstart,l_seq,next_refID,next_pos,tlen,l_read_name,flag,n_cigar_op;
		uint32_t bin_mq_nl,MAPQ,flag_nc;
		bool reversestrand=false;//if true, read is located in reverse strand
		int firstsegment=0;// /1 or /2

		if(gzread(fp, &block_size, sizeof(int32_t))!=sizeof(int32_t)) break;//alignment end
		gzread(fp, &refID, sizeof(int32_t));refID++;

		//
		if(chrindex2name.find(refID)==chrindex2name.end()) {cout<<"wrong refID: "<<refID<<endl;exit(0);}

		gzread(fp, &readstart, sizeof(int32_t));
		readstart++;//in bam, readstart is 0 based
		gzread(fp, &bin_mq_nl, sizeof(uint32_t));
		MAPQ=bin_mq_nl>>8&0xff;
		l_read_name = bin_mq_nl&0xff;

		gzread(fp, &flag_nc, sizeof(uint32_t));
		flag=flag_nc>>16;
		if(flag&16) reversestrand=true;

		if(flag&64) firstsegment=1;
		else if(flag&128) firstsegment=2;

		n_cigar_op=flag_nc&0xffff;
		gzread(fp, &l_seq, sizeof(int32_t));
		gzread(fp, &next_refID, sizeof(int32_t));
		gzread(fp, &next_pos, sizeof(int32_t));
		gzread(fp, &tlen, sizeof(int32_t));

		char read_name[l_read_name];
		gzread(fp, read_name, sizeof(char)*l_read_name);

		uint32_t cigar[n_cigar_op];
		gzread(fp, (char *)&cigar, sizeof(uint32_t)*n_cigar_op);
		stringstream cigar_ss;
		for(int i=0;i<n_cigar_op;i++) cigar_ss<<bam_cigar_oplen(cigar[i])<<bam_cigar_opchr(cigar[i]);

		uint8_t seq[(l_seq+1)>>1];
		gzread(fp, (char *)&seq, sizeof(uint8_t)*((l_seq+1)>>1));

		char qual[l_seq];
		gzread(fp, qual, sizeof(char)*l_seq);

		int32_t aux_size=block_size-32-l_read_name-(n_cigar_op<<2)-((l_seq+1)>>1)-l_seq;
		uint8_t aux[aux_size];
		gzread(fp, aux, sizeof(uint8_t)*aux_size);

		//
		read_length_set.push_back(l_seq);

	}while(i<100);

	//
	int sum=0;
	for(i=0;i<read_length_set.size();i++) sum+=read_length_set[i];

	if(read_length_set.size()==0) {cout<<"#reads in ChIP-seq is less than 100"<<endl;exit(0);}
	if(sum==0) {cout<<"wrong reads length for the first 100 reads of ChIP-seq"<<endl;exit(0);}

	cout<<"read length in ChIP-seq (according to the first "<<read_length_set.size()<<" reads): \t"<<sum/read_length_set.size()<<endl;
	return sum/read_length_set.size();
}

void ReadInputBamfile(const vector<BedRegion> &peakbedregion_set,const string InputBamfile,vector<vector<BamInfor> > &AllPeakInputBamInfor)
{
	//
	const int PeakNo=peakbedregion_set.size();

	int PeakIndex=0;
	string peakchr=peakbedregion_set[0].chr;
	int peakstart=peakbedregion_set[0].start;
	int peakend=peakbedregion_set[0].end;
	vector<BamInfor> PeakInputBamInfor;

	//read bam header
	map<int,string> chrindex2name;

	gzFile fp = gzopen(InputBamfile.c_str(), "r");
	if(fp==NULL) {cout<<"gzopen error "<<InputBamfile<<endl;exit(0);}

	//read header
	char magic[4];
	int32_t l_text,n_ref,l_name,l_ref;

	gzread(fp, magic, sizeof(char)*4);
	gzread(fp, &l_text, sizeof(int32_t));

	char text[l_text];
	gzread(fp, text, sizeof(char)*l_text);
	gzread(fp, &n_ref, sizeof(int32_t));

	for(int i=0;i<n_ref;i++)
	{
		gzread(fp, &l_name, sizeof(int32_t));
		char name[l_name];
		gzread(fp, name, sizeof(char)*l_name);
		gzread(fp, &l_ref, sizeof(int32_t));

		string stemp(name);
		chrindex2name[i+1]=stemp;
	}

	//read bam reads
	int i=0;
	do
	{
		i++;
		//if(i%1000000==0) cout<<"read input bam: "<<i/1000000<<"m"<<endl;
//if(i/1000000==2) break;
		//read a read
		int32_t block_size,refID,readstart,l_seq,next_refID,next_pos,tlen,l_read_name,flag,n_cigar_op;
		uint32_t bin_mq_nl,MAPQ,flag_nc;
		bool reversestrand=false;//if true, read is located in reverse strand
		int firstsegment=0;// /1 or /2

		if(gzread(fp, &block_size, sizeof(int32_t))!=sizeof(int32_t)) break;//alignment end
		gzread(fp, &refID, sizeof(int32_t));refID++;

		//
		if(chrindex2name.find(refID)==chrindex2name.end()) {cout<<"wrong refID: "<<refID<<endl;exit(0);}

		gzread(fp, &readstart, sizeof(int32_t));
		readstart++;//in bam, readstart is 0 based
		gzread(fp, &bin_mq_nl, sizeof(uint32_t));
		MAPQ=bin_mq_nl>>8&0xff;
		l_read_name = bin_mq_nl&0xff;
//cout<<"2"<<endl;
		gzread(fp, &flag_nc, sizeof(uint32_t));
		flag=flag_nc>>16;
		if(flag&16) reversestrand=true;

		if(flag&64) firstsegment=1;
		else if(flag&128) firstsegment=2;

		n_cigar_op=flag_nc&0xffff;
		gzread(fp, &l_seq, sizeof(int32_t));
		gzread(fp, &next_refID, sizeof(int32_t));
		gzread(fp, &next_pos, sizeof(int32_t));
		gzread(fp, &tlen, sizeof(int32_t));

		char read_name[l_read_name];
		gzread(fp, read_name, sizeof(char)*l_read_name);
//cout<<"3 "<<read_name<<endl;
		uint32_t cigar[n_cigar_op];
		gzread(fp, (char *)&cigar, sizeof(uint32_t)*n_cigar_op);
		stringstream cigar_ss;
		for(int i=0;i<n_cigar_op;i++) cigar_ss<<bam_cigar_oplen(cigar[i])<<bam_cigar_opchr(cigar[i]);

		uint8_t seq[(l_seq+1)>>1];
		gzread(fp, (char *)&seq, sizeof(uint8_t)*((l_seq+1)>>1));
//cout<<"4"<<endl;
		char qual[l_seq];
		gzread(fp, qual, sizeof(char)*l_seq);

		int32_t aux_size=block_size-32-l_read_name-(n_cigar_op<<2)-((l_seq+1)>>1)-l_seq;
		uint8_t aux[aux_size];
		gzread(fp, aux, sizeof(uint8_t)*aux_size);

		//filter MAPQ
		//if((int)MAPQ < 30) continue;

//if(chrindex2name[refID]!="chr10") continue;

		uint8_t *s=aux;
		uint8_t *MD=bam_aux_get(s,aux,aux_size,"MD");
		uint8_t *nm=bam_aux_get(s,aux,aux_size,"NM");
		int32_t nm_i=bam_aux2i(nm);//cout<<"NM\t"<<nm_i<<endl;
		if(MD==(uint8_t)0) {cout<<"no MD falg: "<<read_name<<endl;exit(0);}//e.g. CTCF after GATK MONK:252:D1AB6ACXX:1:1202:4974:74864 no MD flag

//cout<<"5"<<endl;
		//get ref seq
		string revisedseq,revisedbq;//compared to ref genome, if has deletion, use '--'; if has insertion, do not show the insertion seq;
		string refseq;
		int readend=readstart-1;
		string cigar_new(cigar_ss.str());
		string MD_str((char *)MD);
//cout<<(char *)read_name<<" "<<readstart<<endl;
		GetRefSeq(seq,l_seq,qual,cigar_new,MD_str,revisedseq,revisedbq,readend,refseq);

		//
		BamInfor mybaminfor;
		mybaminfor.chr=chrindex2name[refID];
		mybaminfor.start=readstart;//sam format, not bed format;
		mybaminfor.end=readend;
		mybaminfor.firstsegment=firstsegment;
		mybaminfor.reversestrand=reversestrand;
		mybaminfor.mapq=MAPQ;
		mybaminfor.cigar_old=cigar_ss.str();
		mybaminfor.cigar=cigar_new;
		mybaminfor.seq=revisedseq;
		mybaminfor.bq=revisedbq;
		mybaminfor.nm=nm_i;
		mybaminfor.md=MD_str;
		mybaminfor.refseq=refseq;
		mybaminfor.readname=(char *)read_name;

//cout<<"6"<<endl;
//cout<<mybaminfor.chr<<" "<<mybaminfor.start<<" "<<mybaminfor.readname<<" "<<mybaminfor.start<<" "<<mybaminfor.end<<endl;
		//
		if(chrindex2name[refID]==peakchr)
		{
			if(readend>=peakstart && readstart<=peakend) PeakInputBamInfor.push_back(mybaminfor);
			else if(readstart>peakend)
			{
				AllPeakInputBamInfor.push_back(PeakInputBamInfor);
				MyFree(PeakInputBamInfor);

				do
				{
					PeakIndex++;
					peakchr=peakbedregion_set[PeakIndex].chr;
					peakstart=peakbedregion_set[PeakIndex].start;
					peakend=peakbedregion_set[PeakIndex].end;

					if(readend>=peakstart && readstart<=peakend) {PeakInputBamInfor.push_back(mybaminfor);break;}
					else if(readstart>peakend) AllPeakInputBamInfor.push_back(PeakInputBamInfor);
					else {cout<<"error read is not in next region"<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}

				}while(1);
			}
			else {cout<<"error not ordered by coordinates: "<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}
		}
		else
		{
			AllPeakInputBamInfor.push_back(PeakInputBamInfor);
			MyFree(PeakInputBamInfor);

			do
			{
				PeakIndex++;
				peakchr=peakbedregion_set[PeakIndex].chr;
				peakstart=peakbedregion_set[PeakIndex].start;
				peakend=peakbedregion_set[PeakIndex].end;

				if(chrindex2name[refID]!=peakchr) {AllPeakInputBamInfor.push_back(PeakInputBamInfor);continue;}

				if(readend>=peakstart && readstart<=peakend) {PeakInputBamInfor.push_back(mybaminfor);break;}
				else if(readstart>peakend) AllPeakInputBamInfor.push_back(PeakInputBamInfor);
				else {cout<<"error read is not in next region"<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}

			}while(1);
		}

	}while(1);

	if(!PeakInputBamInfor.empty())
	{
		AllPeakInputBamInfor.push_back(PeakInputBamInfor);
		MyFree(PeakInputBamInfor);
	}
}

void ReadBamfile(const string PEorSE,const int ReadLength,const vector<BedRegion> &peakbedregion_set,const string Bamfile,const vector<vector<BamInfor> > &AllPeakInputBamInfor,const string fermi_location,const string tmpfilefolder,const string OutputVcffile)
{
	//
	const int PeakNo=peakbedregion_set.size();

	int PeakIndex=0;
cout<<"All Peak No: "<<PeakNo<<endl;
	string peakchr=peakbedregion_set[0].chr;
	int peakstart=peakbedregion_set[0].start;
	int peakend=peakbedregion_set[0].end;
	vector<BamInfor> PeakBamInfor;

	//read bam header
	map<int,string> chrindex2name;

	gzFile fp = gzopen(Bamfile.c_str(), "r");
	if(fp==NULL) {cout<<"gzopen error "<<Bamfile<<endl;exit(0);}

	//read header
	char magic[4];
	int32_t l_text,n_ref,l_name,l_ref;

	gzread(fp, magic, sizeof(char)*4);
	gzread(fp, &l_text, sizeof(int32_t));

	char text[l_text];
	gzread(fp, text, sizeof(char)*l_text);
	gzread(fp, &n_ref, sizeof(int32_t));

	for(int i=0;i<n_ref;i++)
	{
		gzread(fp, &l_name, sizeof(int32_t));
		char name[l_name];
		gzread(fp, name, sizeof(char)*l_name);
		gzread(fp, &l_ref, sizeof(int32_t));

		string stemp(name);
		chrindex2name[i+1]=stemp;
	}

	//read bam reads
	int i=0;
	do
	{
		i++;
		//if(i%1000000==0) cout<<"read bam: "<<i/1000000<<"m"<<endl;
//if(i/1000000==2) break;
		//read a read
		int32_t block_size,refID,readstart,l_seq,next_refID,next_pos,tlen,l_read_name,flag,n_cigar_op;
		uint32_t bin_mq_nl,MAPQ,flag_nc;
		bool reversestrand=false;//if true, read is located in reverse strand
		int firstsegment=0;// /1 or /2

		if(gzread(fp, &block_size, sizeof(int32_t))!=sizeof(int32_t)) break;//alignment end
		gzread(fp, &refID, sizeof(int32_t));refID++;

		//
		if(chrindex2name.find(refID)==chrindex2name.end()) {cout<<"wrong refID: "<<refID<<endl;exit(0);}

		gzread(fp, &readstart, sizeof(int32_t));
		readstart++;//in bam, readstart is 0 based
		gzread(fp, &bin_mq_nl, sizeof(uint32_t));
		MAPQ=bin_mq_nl>>8&0xff;
		l_read_name = bin_mq_nl&0xff;
//cout<<"2"<<endl;
		gzread(fp, &flag_nc, sizeof(uint32_t));
		flag=flag_nc>>16;
		if(flag&16) reversestrand=true;

		if(flag&64) firstsegment=1;
		else if(flag&128) firstsegment=2;

		n_cigar_op=flag_nc&0xffff;
		gzread(fp, &l_seq, sizeof(int32_t));
		gzread(fp, &next_refID, sizeof(int32_t));
		gzread(fp, &next_pos, sizeof(int32_t));
		gzread(fp, &tlen, sizeof(int32_t));

		char read_name[l_read_name];
		gzread(fp, read_name, sizeof(char)*l_read_name);
//cout<<"3 "<<read_name<<endl;
		uint32_t cigar[n_cigar_op];
		gzread(fp, (char *)&cigar, sizeof(uint32_t)*n_cigar_op);
		stringstream cigar_ss;
		for(int i=0;i<n_cigar_op;i++) cigar_ss<<bam_cigar_oplen(cigar[i])<<bam_cigar_opchr(cigar[i]);

		uint8_t seq[(l_seq+1)>>1];
		gzread(fp, (char *)&seq, sizeof(uint8_t)*((l_seq+1)>>1));
//cout<<"4"<<endl;
		char qual[l_seq];
		gzread(fp, qual, sizeof(char)*l_seq);

		int32_t aux_size=block_size-32-l_read_name-(n_cigar_op<<2)-((l_seq+1)>>1)-l_seq;
		uint8_t aux[aux_size];
		gzread(fp, aux, sizeof(uint8_t)*aux_size);

		//filter MAPQ
		//if((int)MAPQ < 30) continue;

//if(chrindex2name[refID]!="chr10") continue;

		uint8_t *s=aux;
		uint8_t *MD=bam_aux_get(s,aux,aux_size,"MD");
		uint8_t *nm=bam_aux_get(s,aux,aux_size,"NM");
		int32_t nm_i=bam_aux2i(nm);//cout<<"NM\t"<<nm_i<<endl;
		if(MD==(uint8_t)0) {cout<<"no MD falg: "<<read_name<<endl;exit(0);}//e.g. CTCF after GATK MONK:252:D1AB6ACXX:1:1202:4974:74864 no MD flag

//cout<<"5"<<endl;
		//get ref seq
		string revisedseq,revisedbq;//compared to ref genome, if has deletion, use '--'; if has insertion, do not show the insertion seq;
		string refseq;
		int readend=readstart-1;
		string cigar_new(cigar_ss.str());
		string MD_str((char *)MD);
		GetRefSeq(seq,l_seq,qual,cigar_new,MD_str,revisedseq,revisedbq,readend,refseq);

		//
		BamInfor mybaminfor;
		mybaminfor.chr=chrindex2name[refID];
		mybaminfor.start=readstart;//sam format, not bed format;
		mybaminfor.end=readend;
		mybaminfor.firstsegment=firstsegment;
		mybaminfor.reversestrand=reversestrand;
		mybaminfor.mapq=MAPQ;
		mybaminfor.cigar_old=cigar_ss.str();
		mybaminfor.cigar=cigar_new;
		mybaminfor.seq=revisedseq;
		mybaminfor.bq=revisedbq;
		mybaminfor.nm=nm_i;
		mybaminfor.md=MD_str;
		mybaminfor.refseq=refseq;
		mybaminfor.readname=(char *)read_name;

//cout<<"6"<<endl;

		//
		if(chrindex2name[refID]==peakchr)
		{
			if(readend>=peakstart && readstart<=peakend) PeakBamInfor.push_back(mybaminfor);
			else if(readstart>peakend)
			{//cout<<peakstart<<" "<<peakend<<" "<<readstart<<endl;
				AssembleAndSNVAS(PEorSE,ReadLength,fermi_location,tmpfilefolder,OutputVcffile,PeakIndex,peakbedregion_set[PeakIndex].chr,PeakBamInfor,AllPeakInputBamInfor[PeakIndex]);
				PeakBamInfor.clear();
				cout<<"finish read ChIP-seq bam in peak: #"<<PeakIndex+1<<endl;

				PeakIndex++;
				if(PeakIndex%100==0)
				{
					//cout<<"read ChIP-seq bam in peak: #"<<PeakIndex<<endl;
					string stemp="rm -f "+tmpfilefolder+"/*";
					system(stemp.c_str());
				}
				peakchr=peakbedregion_set[PeakIndex].chr;
				peakstart=peakbedregion_set[PeakIndex].start;
				peakend=peakbedregion_set[PeakIndex].end;

				if(readend>=peakstart && readstart<=peakend) PeakBamInfor.push_back(mybaminfor);
				else {cout<<"error read is not in next region"<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}
			}
			else {cout<<"error not ordered by coordinates: "<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}
		}
		else
		{
			AssembleAndSNVAS(PEorSE,ReadLength,fermi_location,tmpfilefolder,OutputVcffile,PeakIndex,peakbedregion_set[PeakIndex].chr,PeakBamInfor,AllPeakInputBamInfor[PeakIndex]);
			PeakBamInfor.clear();
			cout<<"finish read ChIP-seq bam in peak: #"<<PeakIndex+1<<endl;

			PeakIndex++;
			if(PeakIndex%100==0)
			{
				//cout<<"read ChIP-seq bam in peak: #"<<PeakIndex<<endl;
				string stemp="rm -f "+tmpfilefolder+"/*";
				system(stemp.c_str());
			}
			peakchr=peakbedregion_set[PeakIndex].chr;
			peakstart=peakbedregion_set[PeakIndex].start;
			peakend=peakbedregion_set[PeakIndex].end;

			if(chrindex2name[refID]!=peakchr) {cout<<"error read is not in next chr region: "<<chrindex2name[refID]<<"-"<<peakchr<<" "<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}

			if(readend>=peakstart && readstart<=peakend) PeakBamInfor.push_back(mybaminfor);
			else {cout<<"error read is not in next region of new chr: "<<readstart<<" "<<peakstart<<"-"<<peakend<<endl;exit(0);}
		}

	}while(1);

	if(!PeakBamInfor.empty())
	{
		AssembleAndSNVAS(PEorSE,ReadLength,fermi_location,tmpfilefolder,OutputVcffile,PeakIndex,peakbedregion_set[PeakIndex].chr,PeakBamInfor,AllPeakInputBamInfor[PeakIndex]);
		PeakBamInfor.clear();
		cout<<"finish read ChIP-seq bam in peak: #"<<PeakIndex+1<<endl;

		string stemp="rm -f "+tmpfilefolder+"/*";
		system(stemp.c_str());
	}
}

void GetReverseComplementary(const string &input,string &output)
{
	output=input;

	int length=input.size();
	for(int i=0;i<input.size();i++)
	{
		if(input[i]=='A') output[length-i-1]='T';
		else if(input[i]=='C') output[length-i-1]='G';
		else if(input[i]=='G') output[length-i-1]='C';
		else if(input[i]=='T') output[length-i-1]='A';
		else if(input[i]=='N') output[length-i-1]='N';
		else {cout<<"wrong ref seq1 nt: "<<i+1<<" "<<input[i]<<endl;exit(0);}
	}
}

void AssembleAndSNVAS(const string PEorSE,const int ReadLength,const string fermi_location,const string tmpfilefolder,const string OutputVcffile,const int PeakIndex,const string regionchr,const vector<BamInfor> &PeakBamInfor,const vector<BamInfor> &PeakInputBamInfor)
{
	if(PeakBamInfor.size()==0) return;

	const int Fermi_overlap_par=ReadLength/2;

	int i,j;
	string sbuf;

	if(PEorSE=="PE")
	{
		//get seq and bq, and ref seq in the extended peak region
		stringstream ss1;
		ss1<<tmpfilefolder<<"/"<<PeakIndex+1<<".fastq";

		ofstream os1(ss1.str().c_str());
		if(!os1) {cout<<"can not open tmp output file: "<<ss1.str()<<endl;exit(0);}

		//
		string extendrefseq=PeakBamInfor[0].refseq;
		int extendref_start=PeakBamInfor[0].start;//sam format
		int extendref_end=PeakBamInfor[0].end;//sam format
		int index_old=0;

		for(i=0;i<PeakBamInfor.size();i++)
		{
			BamInfor mybaminfor=PeakBamInfor[i];

			//
			if(mybaminfor.refseq.size()!=mybaminfor.end-mybaminfor.start+1) {cout<<"wrong read reference sequence"<<endl;exit(0);}
	//cout<<"ChIP "<<i<<" "<<PeakBamInfor[i].start<<"-"<<PeakBamInfor[i].end<<endl;
			if(i>0)
			{
				int ia=mybaminfor.end-PeakBamInfor[index_old].end;
				int ib=mybaminfor.start-PeakBamInfor[index_old].end-1;
				if(ia>0)
				{
	//cout<<extendrefseq<<endl;
					if(ib>=1)//there is gap between two reads
					{
						string stemp;
						for(j=0;j<ib;j++) stemp.push_back('N');
						extendrefseq=extendrefseq+stemp+mybaminfor.refseq;
					}
					else extendrefseq=extendrefseq+mybaminfor.refseq.substr(mybaminfor.refseq.size()-ia);

	//cout<<index_old<<"\t"<<PeakBamInfor[index_old].start<<"-"<<PeakBamInfor[index_old].end<<"\t"<<PeakBamInfor[index_old].refseq<<endl;
	//cout<<i<<"\t"<<PeakBamInfor[i].start<<"-"<<PeakBamInfor[i].end<<"\t"<<PeakBamInfor[i].refseq<<endl;
	//cout<<extendrefseq<<endl<<endl;
					index_old=i;

					extendref_end=mybaminfor.end;
				}
			}

			//
			if(mybaminfor.firstsegment==1)
			{
				os1<<"@"<<mybaminfor.readname<<"/1"<<endl;
				os1<<mybaminfor.seq<<endl;
				os1<<"+"<<endl;
				os1<<mybaminfor.bq<<endl;
			}
			else if(mybaminfor.firstsegment==2)
			{
				os1<<"@"<<mybaminfor.readname<<"/2"<<endl;
				os1<<mybaminfor.seq<<endl;
				os1<<"+"<<endl;
				os1<<mybaminfor.bq<<endl;
			}
			else {cout<<"wrong segment: "<<mybaminfor.firstsegment<<endl;exit(0);}
		}
	//cout<<"ChIPextend_infor:\t"<<extendref_start<<"-"<<extendref_end<<endl<<extendrefseq<<endl;

		//
		string extendrefseq_input;
		int extendref_input_start=0;
		int extendref_input_end=0;//sam format
		index_old=0;
		if(PeakInputBamInfor.size()>0)
		{
			extendrefseq_input=PeakInputBamInfor[0].refseq;
			extendref_input_start=PeakInputBamInfor[0].start;//sam format
			extendref_input_end=PeakInputBamInfor[0].end;
		}

		for(i=0;i<PeakInputBamInfor.size();i++)
		{
			BamInfor mybaminfor=PeakInputBamInfor[i];

			//
			if(mybaminfor.refseq.size()!=mybaminfor.end-mybaminfor.start+1) {cout<<"wrong read reference sequence"<<endl;exit(0);}
	//cout<<"Input "<<i<<" "<<PeakInputBamInfor[i].start<<"-"<<PeakInputBamInfor[i].end<<endl;
			if(i>0)
			{
				int ia=mybaminfor.end-PeakInputBamInfor[index_old].end;
				int ib=mybaminfor.start-PeakInputBamInfor[index_old].end-1;
				if(ia>0)
				{
	//cout<<extendrefseq_input<<endl;
					if(ib>=1)//there is gap between two reads
					{
						string stemp;
						for(j=0;j<ib;j++) stemp.push_back('N');
						extendrefseq_input=extendrefseq_input+stemp+mybaminfor.refseq;
					}
					else extendrefseq_input=extendrefseq_input+mybaminfor.refseq.substr(mybaminfor.refseq.size()-ia);

	//cout<<index_old<<"\t"<<PeakInputBamInfor[index_old].start<<"-"<<PeakInputBamInfor[index_old].end<<"\t"<<PeakInputBamInfor[index_old].refseq<<endl;
	//cout<<i<<"\t"<<PeakInputBamInfor[i].start<<"-"<<PeakInputBamInfor[i].end<<"\t"<<PeakInputBamInfor[i].refseq<<endl;
	//cout<<extendrefseq_input<<endl<<endl;
					index_old=i;

					extendref_input_end=mybaminfor.end;
				}
			}

			//
			if(mybaminfor.firstsegment==1)
			{
				os1<<"@"<<mybaminfor.readname<<"/1"<<endl;
				os1<<mybaminfor.seq<<endl;
				os1<<"+"<<endl;
				os1<<mybaminfor.bq<<endl;
			}
			else if(mybaminfor.firstsegment==2)
			{
				os1<<"@"<<mybaminfor.readname<<"/2"<<endl;
				os1<<mybaminfor.seq<<endl;
				os1<<"+"<<endl;
				os1<<mybaminfor.bq<<endl;
			}
			else {cout<<"wrong segment: "<<mybaminfor.firstsegment<<endl;exit(0);}
		}
	//cout<<"Inputextend_infor:\t"<<extendref_input_start<<"-"<<extendref_input_end<<endl<<extendrefseq_input<<endl;
		os1.close();

		//
		if(extendref_end-extendref_start+1 != extendrefseq.size()) {cout<<"wrong ref seq for ChIP: "<<endl;exit(0);}
		if(extendrefseq_input.size()>0) {if(extendref_input_end-extendref_input_start+1 != extendrefseq_input.size()) {cout<<"wrong ref seq for input: "<<endl;exit(0);}}

		map<int,char> pos2ref;
		for(i=0;i<extendrefseq.size();i++)
		{
			pos2ref[extendref_start++]=extendrefseq[i];
		}
		for(i=0;i<extendrefseq_input.size();i++)
		{
			if(extendrefseq_input[i]!='N') pos2ref[extendref_input_start++]=extendrefseq_input[i];
			else extendref_input_start++;
		}
		extendrefseq.clear();
		map<int,char>::iterator pos;
		int itian;
		for(pos=pos2ref.begin();pos!=pos2ref.end();++pos)
		{
			if(pos!=pos2ref.begin())
			{
				if(pos->first != itian+1)
				{
					int itemp=pos->first - itian-1;
					if(itemp<0) {cout<<"map index wrong: "<<pos->first<<" "<<itian<<endl;exit(0);}
					for(i=0;i<itemp;i++) extendrefseq.push_back('N');
				}

				extendrefseq.push_back(pos->second);
			}
			else extendrefseq.push_back(pos->second);

			itian=pos->first;
		}
		pos=pos2ref.begin();
		extendref_start=pos->first;
		pos=pos2ref.end();
		--pos;
		extendref_end=pos->first;
		if(extendref_end-extendref_start+1!=extendrefseq.size()) {cout<<"wrong extend ref: "<<extendref_end<<" "<<extendref_start<<" "<<extendrefseq.size()<<endl;exit(0);}
	//cout<<"Finalextend_infor:\t"<<extendref_start<<"-"<<extendref_end<<endl<<extendrefseq<<endl;


		//fermi
		stringstream p1_mag,p1_mag_log,sstemp;
		p1_mag<<tmpfilefolder<<"/"<<PeakIndex+1<<"_p1.mag";
		p1_mag_log<<tmpfilefolder<<"/"<<PeakIndex+1<<"_p1.mag.log";
		string cmdstring;
		sstemp<<Fermi_overlap_par;
		cmdstring=fermi_location+" -ce -l "+sstemp.str()+" "+ss1.str()+" >"+p1_mag.str()+" 2>"+p1_mag_log.str();
		system(cmdstring.c_str());

		//read contig
		vector<string> contigname_set;
		vector<string> contigseq_set;

		ifstream is(p1_mag.str().c_str());
		if(!is) {cout<<"can not open: "<<p1_mag.str()<<endl;return;}
		do
		{
			sbuf.clear();
			is>>sbuf;
			if(sbuf.empty()) break;
			contigname_set.push_back(sbuf);

			getline(is,sbuf);
			getline(is,sbuf);contigseq_set.push_back(sbuf);
			getline(is,sbuf);getline(is,sbuf);

		}while(1);
		is.close();

		//local alignment
		vector<string> contigseq,refseq;//maybe has '-'
		vector<int> refstart,refend;//for each alignment, sam format
		vector<vector<int> > contigcoor;//if there is insertion, the coor is -1;

		bool btemp=false;
		for(i=0;i<contigseq_set.size();i++)
		{
			seq_pair problem1,problem2;
			seq_pair result1,result2;

			//
			problem1.a=contigseq_set[i];
			problem1.b=extendrefseq;
			smith_waterman(problem1, btemp, result1);

			//
			string reversecomplementary;
			GetReverseComplementary(contigseq_set[i],reversecomplementary);

			problem2.a=reversecomplementary;
			problem2.b=extendrefseq;
			smith_waterman(problem2, btemp, result2);
	//cout<<i+1<<endl<<problem1.a<<endl<<problem2.a<<endl;
			//
			string mycontigseq,myrefseq;//maybe has '-'
			double myscore;
			if(result1.score >= result2.score) {mycontigseq=result1.a;myrefseq=result1.b;myscore=result1.score;}
			else {mycontigseq=result2.a;myrefseq=result2.b;myscore=result2.score;}

			//
			string myrefseq_noinsert;
			for(j=0;j<myrefseq.size();j++)
			{
				if(myrefseq[j]!='-') myrefseq_noinsert.push_back(myrefseq[j]);
			}
			int ia=0;
			int ib=0;

			if(myrefseq_noinsert!=extendrefseq)
			{
				string::size_type position=extendrefseq.find(myrefseq_noinsert);
				if(position==extendrefseq.npos) {cout<<"error!"<<endl;exit(0);}

				ia=(int)position;
				ib=extendrefseq.size()-ia-myrefseq_noinsert.size();
			}

			//
			int itemp1=0;
			int itemp2=0;
			for(j=0;j<mycontigseq.size();j++)
			{
				if(mycontigseq[j]!='-') break;
				itemp1++;
			}
			for(j=mycontigseq.size()-1;j>=0;j--)
			{
				if(mycontigseq[j]!='-') break;
				itemp2++;
			}

			refstart.push_back(extendref_start+itemp1+ia);
			refend.push_back(extendref_end-itemp2-ib);
			contigseq.push_back(mycontigseq.substr(itemp1,mycontigseq.size()-itemp2-itemp1));
			refseq.push_back(myrefseq.substr(itemp1,mycontigseq.size()-itemp2-itemp1));
		}
	//cout<<contigseq.size()<<endl;
		for(i=0;i<contigseq.size();i++)
		{
			vector<int> mycontigcoor;

			int ia=refstart[i]-1;
			for(j=0;j<refseq[i].size();j++)
			{
				if(refseq[i][j]=='-') mycontigcoor.push_back(-1);
				else {++ia;mycontigcoor.push_back(ia);}
			}
	//cout<<i+1<<endl<<refseq[i]<<endl<<contigseq[i]<<endl<<refstart[i]<<"-"<<refend[i]<<endl;
			if(ia!=refend[i]) {cout<<"contig coor unconsistent: "<<ia<<" "<<refstart[i]<<" "<<refend[i]<<endl;exit(0);}

			contigcoor.push_back(mycontigcoor);
		}

		//define snv from contig mapping result
		map<int,char> snvpos2refallele;
		map<int,vector<char> > snvpos2contigNTs;//maybe there is redundant NTs

		for(i=0;i<refstart.size();i++)
		{
			int pos=refstart[i];

			for(j=0;j<refseq[i].size();j++)
			{
				if(refseq[i][j]=='-') continue;
				if(refseq[i][j]=='N') {pos++;continue;}

				if(contigseq[i][j]=='-' || contigseq[i][j]=='N') {pos++;continue;}

				//
				if(refseq[i][j]!=contigseq[i][j])
				{
					if(snvpos2refallele.find(pos)==snvpos2refallele.end()) snvpos2refallele[pos]=refseq[i][j];
				}

				//
				if(snvpos2contigNTs.find(pos)==snvpos2contigNTs.end())
				{
					vector<char> vchtemp;
					vchtemp.push_back(contigseq[i][j]);
					snvpos2contigNTs[pos]=vchtemp;
				}
				else snvpos2contigNTs[pos].push_back(contigseq[i][j]);

				//
				if(j==refseq[i].size()-1)
				{
					if(pos!=refend[i]) {cout<<"wrong start end for contig: "<<contigseq[i]<<endl;exit(0);}
				}

				pos++;
			}
		}

		//Fill snv infor: ref and fermiNTs, and initialize
		map<int,PosReadsInfor> pos2Readsinfo;

		map<int,char>::iterator posmich;
		for(posmich=snvpos2refallele.begin();posmich!=snvpos2refallele.end();++posmich)
		{
			FillContigInfor(posmich->first,posmich->second,snvpos2contigNTs[posmich->first],pos2Readsinfo);
		}

		//reads mapping result
		vector<string> readsseq,readsbq,inputreadsseq,inputreadsbq;
		vector<vector<int> > readscoor,inputreadscoor;//only select reads overlapped with snvs, and coor seq are the same as contig

		for(i=0;i<PeakBamInfor.size();i++)
		{
			if(PeakBamInfor[i].mapq<30) continue;

			string revisedseq,revisedbq;
			vector<int> revisedcoor;
			GetReadSeqCoor(PeakBamInfor[i].seq,PeakBamInfor[i].bq,PeakBamInfor[i].start,PeakBamInfor[i].cigar,revisedseq,revisedbq,revisedcoor);
			readsseq.push_back(revisedseq);
			readsbq.push_back(revisedbq);
			readscoor.push_back(revisedcoor);
		}

		for(i=0;i<PeakInputBamInfor.size();i++)
		{
			if(PeakInputBamInfor[i].mapq<30) continue;

			string revisedseq,revisedbq;
			vector<int> revisedcoor;
			GetReadSeqCoor(PeakInputBamInfor[i].seq,PeakInputBamInfor[i].bq,PeakInputBamInfor[i].start,PeakInputBamInfor[i].cigar,revisedseq,revisedbq,revisedcoor);
			inputreadscoor.push_back(revisedcoor);
			inputreadsseq.push_back(revisedseq);
			inputreadsbq.push_back(revisedbq);
		}

		//Fill snv infor: ChIP reads infor
		vector<int>::iterator it;
		btemp=false;
		for(i=0;i<readscoor.size();i++)
		{
			for(posmich=snvpos2refallele.begin();posmich!=snvpos2refallele.end();++posmich)
			{
				it=lower_bound(readscoor[i].begin(),readscoor[i].end(),posmich->first);

				if(it==readscoor[i].end() || *it!=posmich->first)//not found
				{
					continue;
				}
				else//found
				{
					int index=(int)(it-readscoor[i].begin());

					FillChIPraw(posmich->first,readsseq[i][index],pos2Readsinfo);

					btemp=ConsistentWithContig(readscoor[i],readsseq[i],contigcoor,contigseq);
					if(btemp) FillChIPQualInfor(posmich->first,readsseq[i][index],readsbq[i][index],pos2Readsinfo);
				}
			}
		}

		//Fill snv infor: Input reads infor
		btemp=false;
		for(i=0;i<inputreadscoor.size();i++)
		{
			for(posmich=snvpos2refallele.begin();posmich!=snvpos2refallele.end();++posmich)
			{
				it=lower_bound(inputreadscoor[i].begin(),inputreadscoor[i].end(),posmich->first);

				if(it==inputreadscoor[i].end() || *it!=posmich->first)//not found
				{
					continue;
				}
				else//found
				{
					int index=(int)(it-inputreadscoor[i].begin());

					FillControlraw(posmich->first,inputreadsseq[i][index],pos2Readsinfo);

					btemp=ConsistentWithContig(inputreadscoor[i],inputreadsseq[i],contigcoor,contigseq);
					if(btemp) FillControlQualInfor(posmich->first,inputreadsseq[i][index],inputreadsbq[i][index],pos2Readsinfo);
				}
			}
		}

		//
		CalSNVAS(pos2Readsinfo);

		//
		OutputVcfResultHasInput(OutputVcffile,regionchr,pos2Readsinfo);
	}
	else
	{
		//get seq and bq, and ref seq in the extended peak region
		stringstream ss1;
		ss1<<tmpfilefolder<<"/"<<PeakIndex<<".fastq";

		ofstream os1(ss1.str().c_str());
		if(!os1) {cout<<"can not open tmp output file: "<<ss1.str()<<endl;exit(0);}
		//
		string extendrefseq=PeakBamInfor[0].refseq;
		int extendref_start=PeakBamInfor[0].start;//sam format
		int extendref_end=PeakBamInfor[0].end;//sam format
		int index_old=0;

		for(i=0;i<PeakBamInfor.size();i++)
		{
			BamInfor mybaminfor=PeakBamInfor[i];

			//
			if(mybaminfor.refseq.size()!=mybaminfor.end-mybaminfor.start+1) {cout<<"wrong read reference sequence"<<endl;exit(0);}
	//cout<<"ChIP "<<i<<" "<<PeakBamInfor[i].start<<"-"<<PeakBamInfor[i].end<<endl;
			if(i>0)
			{
				int ia=mybaminfor.end-PeakBamInfor[index_old].end;
				int ib=mybaminfor.start-PeakBamInfor[index_old].end-1;
				if(ia>0)
				{
	//cout<<extendrefseq<<endl;
					if(ib>=1)//there is gap between two reads
					{
						string stemp;
						for(j=0;j<ib;j++) stemp.push_back('N');
						extendrefseq=extendrefseq+stemp+mybaminfor.refseq;
					}
					else extendrefseq=extendrefseq+mybaminfor.refseq.substr(mybaminfor.refseq.size()-ia);

	//cout<<index_old<<"\t"<<PeakBamInfor[index_old].start<<"-"<<PeakBamInfor[index_old].end<<"\t"<<PeakBamInfor[index_old].refseq<<endl;
	//cout<<i<<"\t"<<PeakBamInfor[i].start<<"-"<<PeakBamInfor[i].end<<"\t"<<PeakBamInfor[i].refseq<<endl;
	//cout<<extendrefseq<<endl<<endl;
					index_old=i;

					extendref_end=mybaminfor.end;
				}
			}

			//
			if(mybaminfor.firstsegment==0)
			{
				os1<<"@"<<mybaminfor.readname<<endl;
				os1<<mybaminfor.seq<<endl;
				os1<<"+"<<endl;
				os1<<mybaminfor.bq<<endl;
			}
			else {cout<<"wrong segment: "<<mybaminfor.firstsegment<<endl;exit(0);}
		}
	//cout<<"ChIPextend_infor:\t"<<extendref_start<<"-"<<extendref_end<<endl<<extendrefseq<<endl;

		//
		string extendrefseq_input;
		int extendref_input_start=0;
		int extendref_input_end=0;//sam format
		index_old=0;
		if(PeakInputBamInfor.size()>0)
		{
			extendrefseq_input=PeakInputBamInfor[0].refseq;
			extendref_input_start=PeakInputBamInfor[0].start;//sam format
			extendref_input_end=PeakInputBamInfor[0].end;
		}

		for(i=0;i<PeakInputBamInfor.size();i++)
		{
			BamInfor mybaminfor=PeakInputBamInfor[i];

			//
			if(mybaminfor.refseq.size()!=mybaminfor.end-mybaminfor.start+1) {cout<<"wrong read reference sequence"<<endl;exit(0);}
	//cout<<"Input "<<i<<" "<<PeakInputBamInfor[i].start<<"-"<<PeakInputBamInfor[i].end<<endl;
			if(i>0)
			{
				int ia=mybaminfor.end-PeakInputBamInfor[index_old].end;
				int ib=mybaminfor.start-PeakInputBamInfor[index_old].end-1;
				if(ia>0)
				{
	//cout<<extendrefseq_input<<endl;
					if(ib>=1)//there is gap between two reads
					{
						string stemp;
						for(j=0;j<ib;j++) stemp.push_back('N');
						extendrefseq_input=extendrefseq_input+stemp+mybaminfor.refseq;
					}
					else extendrefseq_input=extendrefseq_input+mybaminfor.refseq.substr(mybaminfor.refseq.size()-ia);

	//cout<<index_old<<"\t"<<PeakInputBamInfor[index_old].start<<"-"<<PeakInputBamInfor[index_old].end<<"\t"<<PeakInputBamInfor[index_old].refseq<<endl;
	//cout<<i<<"\t"<<PeakInputBamInfor[i].start<<"-"<<PeakInputBamInfor[i].end<<"\t"<<PeakInputBamInfor[i].refseq<<endl;
	//cout<<extendrefseq_input<<endl<<endl;
					index_old=i;

					extendref_input_end=mybaminfor.end;
				}
			}

			//
			if(mybaminfor.firstsegment==0)
			{
				os1<<"@"<<mybaminfor.readname<<endl;
				os1<<mybaminfor.seq<<endl;
				os1<<"+"<<endl;
				os1<<mybaminfor.bq<<endl;
			}
			else {cout<<"wrong segment: "<<mybaminfor.firstsegment<<endl;exit(0);}
		}
	//cout<<"Inputextend_infor:\t"<<extendref_input_start<<"-"<<extendref_input_end<<endl<<extendrefseq_input<<endl;
		os1.close();

		//
		if(extendref_end-extendref_start+1 != extendrefseq.size()) {cout<<"wrong ref seq for ChIP: "<<endl;exit(0);}
		if(extendrefseq_input.size()>0) {if(extendref_input_end-extendref_input_start+1 != extendrefseq_input.size()) {cout<<"wrong ref seq for input: "<<endl;exit(0);}}

		map<int,char> pos2ref;
		for(i=0;i<extendrefseq.size();i++)
		{
			pos2ref[extendref_start++]=extendrefseq[i];
		}
		for(i=0;i<extendrefseq_input.size();i++)
		{
			if(extendrefseq_input[i]!='N') pos2ref[extendref_input_start++]=extendrefseq_input[i];
			else extendref_input_start++;
		}
		extendrefseq.clear();
		map<int,char>::iterator pos;
		int itian;
		for(pos=pos2ref.begin();pos!=pos2ref.end();++pos)
		{
			if(pos!=pos2ref.begin())
			{
				if(pos->first != itian+1)
				{
					int itemp=pos->first - itian-1;
					if(itemp<0) {cout<<"map index wrong: "<<pos->first<<" "<<itian<<endl;exit(0);}
					for(i=0;i<itemp;i++) extendrefseq.push_back('N');
				}

				extendrefseq.push_back(pos->second);
			}
			else extendrefseq.push_back(pos->second);

			itian=pos->first;
		}
		pos=pos2ref.begin();
		extendref_start=pos->first;
		pos=pos2ref.end();
		--pos;
		extendref_end=pos->first;
		if(extendref_end-extendref_start+1!=extendrefseq.size()) {cout<<"wrong extend ref: "<<extendref_end<<" "<<extendref_start<<" "<<extendrefseq.size()<<endl;exit(0);}
	//cout<<"Finalextend_infor:\t"<<extendref_start<<"-"<<extendref_end<<endl<<extendrefseq<<endl;


		//fermi
		stringstream p1_mag,p1_mag_log,sstemp;
		p1_mag<<tmpfilefolder<<"/"<<PeakIndex+1<<"_p1.mag";
		p1_mag_log<<tmpfilefolder<<"/"<<PeakIndex+1<<"_p1.mag.log";
		string cmdstring;
		sstemp<<Fermi_overlap_par;
		cmdstring=fermi_location+" -ce -l "+sstemp.str()+" "+ss1.str()+" >"+p1_mag.str()+" 2>"+p1_mag_log.str();
		system(cmdstring.c_str());

		//read contig
		vector<string> contigname_set;
		vector<string> contigseq_set;

		ifstream is(p1_mag.str().c_str());
		if(!is) {cout<<"can not open: "<<p1_mag.str()<<endl;return;}
		do
		{
			sbuf.clear();
			is>>sbuf;
			if(sbuf.empty()) break;
			contigname_set.push_back(sbuf);

			getline(is,sbuf);
			getline(is,sbuf);contigseq_set.push_back(sbuf);
			getline(is,sbuf);getline(is,sbuf);

		}while(1);
		is.close();

		//local alignment
		vector<string> contigseq,refseq;//maybe has '-'
		vector<int> refstart,refend;//for each alignment, sam format
		vector<vector<int> > contigcoor;//if there is insertion, the coor is -1;

		bool btemp=false;
		for(i=0;i<contigseq_set.size();i++)
		{
			seq_pair problem1,problem2;
			seq_pair result1,result2;

			//
			problem1.a=contigseq_set[i];
			problem1.b=extendrefseq;
			smith_waterman(problem1, btemp, result1);

			//
			string reversecomplementary;
			GetReverseComplementary(contigseq_set[i],reversecomplementary);

			problem2.a=reversecomplementary;
			problem2.b=extendrefseq;
			smith_waterman(problem2, btemp, result2);
	//cout<<i+1<<endl<<problem1.a<<endl<<problem2.a<<endl;
			//
			string mycontigseq,myrefseq;//maybe has '-'
			double myscore;
			if(result1.score >= result2.score) {mycontigseq=result1.a;myrefseq=result1.b;myscore=result1.score;}
			else {mycontigseq=result2.a;myrefseq=result2.b;myscore=result2.score;}

			//
			string myrefseq_noinsert;
			for(j=0;j<myrefseq.size();j++)
			{
				if(myrefseq[j]!='-') myrefseq_noinsert.push_back(myrefseq[j]);
			}
			int ia=0;
			int ib=0;

			if(myrefseq_noinsert!=extendrefseq)
			{
				string::size_type position=extendrefseq.find(myrefseq_noinsert);
				if(position==extendrefseq.npos) {cout<<"error!"<<endl;exit(0);}

				ia=(int)position;
				ib=extendrefseq.size()-ia-myrefseq_noinsert.size();
			}

			//
			int itemp1=0;
			int itemp2=0;
			for(j=0;j<mycontigseq.size();j++)
			{
				if(mycontigseq[j]!='-') break;
				itemp1++;
			}
			for(j=mycontigseq.size()-1;j>=0;j--)
			{
				if(mycontigseq[j]!='-') break;
				itemp2++;
			}

			refstart.push_back(extendref_start+itemp1+ia);
			refend.push_back(extendref_end-itemp2-ib);
			contigseq.push_back(mycontigseq.substr(itemp1,mycontigseq.size()-itemp2-itemp1));
			refseq.push_back(myrefseq.substr(itemp1,mycontigseq.size()-itemp2-itemp1));
		}
	//cout<<contigseq.size()<<endl;
		for(i=0;i<contigseq.size();i++)
		{
			vector<int> mycontigcoor;

			int ia=refstart[i]-1;
			for(j=0;j<refseq[i].size();j++)
			{
				if(refseq[i][j]=='-') mycontigcoor.push_back(-1);
				else {++ia;mycontigcoor.push_back(ia);}
			}
	//cout<<i+1<<endl<<refseq[i]<<endl<<contigseq[i]<<endl<<refstart[i]<<"-"<<refend[i]<<endl;
			if(ia!=refend[i]) {cout<<"contig coor unconsistent: "<<ia<<" "<<refstart[i]<<" "<<refend[i]<<endl;exit(0);}

			contigcoor.push_back(mycontigcoor);
		}

		//define snv from contig mapping result
		map<int,char> snvpos2refallele;
		map<int,vector<char> > snvpos2contigNTs;//maybe there is redundant NTs

		for(i=0;i<refstart.size();i++)
		{
			int pos=refstart[i];

			for(j=0;j<refseq[i].size();j++)
			{
				if(refseq[i][j]=='-') continue;
				if(refseq[i][j]=='N') {pos++;continue;}

				if(contigseq[i][j]=='-' || contigseq[i][j]=='N') {pos++;continue;}

				//
				if(refseq[i][j]!=contigseq[i][j])
				{
					if(snvpos2refallele.find(pos)==snvpos2refallele.end()) snvpos2refallele[pos]=refseq[i][j];
				}

				//
				if(snvpos2contigNTs.find(pos)==snvpos2contigNTs.end())
				{
					vector<char> vchtemp;
					vchtemp.push_back(contigseq[i][j]);
					snvpos2contigNTs[pos]=vchtemp;
				}
				else snvpos2contigNTs[pos].push_back(contigseq[i][j]);

				//
				if(j==refseq[i].size()-1)
				{
					if(pos!=refend[i]) {cout<<"wrong start end for contig: "<<contigseq[i]<<endl;exit(0);}
				}

				pos++;
			}
		}

		//Fill snv infor: ref and fermiNTs, and initialize
		map<int,PosReadsInfor> pos2Readsinfo;

		map<int,char>::iterator posmich;
		for(posmich=snvpos2refallele.begin();posmich!=snvpos2refallele.end();++posmich)
		{
			FillContigInfor(posmich->first,posmich->second,snvpos2contigNTs[posmich->first],pos2Readsinfo);
		}

		//reads mapping result
		vector<string> readsseq,readsbq,inputreadsseq,inputreadsbq;
		vector<vector<int> > readscoor,inputreadscoor;//only select reads overlapped with snvs, and coor seq are the same as contig

		for(i=0;i<PeakBamInfor.size();i++)
		{
			if(PeakBamInfor[i].mapq<30) continue;

			string revisedseq,revisedbq;
			vector<int> revisedcoor;
			GetReadSeqCoor(PeakBamInfor[i].seq,PeakBamInfor[i].bq,PeakBamInfor[i].start,PeakBamInfor[i].cigar,revisedseq,revisedbq,revisedcoor);
			readsseq.push_back(revisedseq);
			readsbq.push_back(revisedbq);
			readscoor.push_back(revisedcoor);
		}

		for(i=0;i<PeakInputBamInfor.size();i++)
		{
			if(PeakInputBamInfor[i].mapq<30) continue;

			string revisedseq,revisedbq;
			vector<int> revisedcoor;
			GetReadSeqCoor(PeakInputBamInfor[i].seq,PeakInputBamInfor[i].bq,PeakInputBamInfor[i].start,PeakInputBamInfor[i].cigar,revisedseq,revisedbq,revisedcoor);
			inputreadscoor.push_back(revisedcoor);
			inputreadsseq.push_back(revisedseq);
			inputreadsbq.push_back(revisedbq);
		}

		//Fill snv infor: ChIP reads infor
		vector<int>::iterator it;
		btemp=false;
		for(i=0;i<readscoor.size();i++)
		{
			for(posmich=snvpos2refallele.begin();posmich!=snvpos2refallele.end();++posmich)
			{
				it=lower_bound(readscoor[i].begin(),readscoor[i].end(),posmich->first);

				if(it==readscoor[i].end() || *it!=posmich->first)//not found
				{
					continue;
				}
				else//found
				{
					int index=(int)(it-readscoor[i].begin());

					FillChIPraw(posmich->first,readsseq[i][index],pos2Readsinfo);

					btemp=ConsistentWithContig(readscoor[i],readsseq[i],contigcoor,contigseq);
					if(btemp) FillChIPQualInfor(posmich->first,readsseq[i][index],readsbq[i][index],pos2Readsinfo);
				}
			}
		}

		//Fill snv infor: Input reads infor
		btemp=false;
		for(i=0;i<inputreadscoor.size();i++)
		{
			for(posmich=snvpos2refallele.begin();posmich!=snvpos2refallele.end();++posmich)
			{
				it=lower_bound(inputreadscoor[i].begin(),inputreadscoor[i].end(),posmich->first);

				if(it==inputreadscoor[i].end() || *it!=posmich->first)//not found
				{
					continue;
				}
				else//found
				{
					int index=(int)(it-inputreadscoor[i].begin());

					FillControlraw(posmich->first,inputreadsseq[i][index],pos2Readsinfo);

					btemp=ConsistentWithContig(inputreadscoor[i],inputreadsseq[i],contigcoor,contigseq);
					if(btemp) FillControlQualInfor(posmich->first,inputreadsseq[i][index],inputreadsbq[i][index],pos2Readsinfo);
				}
			}
		}

		//
		CalSNVAS(pos2Readsinfo);

		//
		OutputVcfResultHasInput(OutputVcffile,regionchr,pos2Readsinfo);
	}
}

void GetReadSeqCoor(const string seq,const string bq,const int startpos,const string cigar,string &revisedseq,string &revisedbq,vector<int> &revisedcoor)
{
	int i,j,itemp,ia;

	revisedseq=seq;
	revisedbq=bq;

	//there is no 'S' or 'H' in cigar
	int start=0;
	int end=0;

	if(cigar.substr(0,2)=="0I" || cigar.substr(0,2)=="0D" || cigar.substr(0,2)=="0M") {cout<<"error1:"<<cigar<<endl;exit(0);}

	int seq_index=0;//index in read seq
	int coor_index=startpos;//position in ref genome
	for(i=0;i<cigar.length();i++)
	{
		if(cigar[i]<'0' || cigar[i]>'9')
		{
			end=i-1;
			string stemp=cigar.substr(start,end-start+1);

			itemp=atoi(stemp.c_str());
			start=i+1;

			if(!stemp.empty())
			{
				if(cigar[i]=='I')//raw reads longer than ref
				{
					for(j=0;j<itemp;j++) revisedcoor.push_back(-1);
					seq_index=seq_index+itemp;
				}
				else if(cigar[i]=='D')//raw reads shorter than ref
				{
					string sa,sb;
					for(j=0;j<itemp;j++) {sa.push_back('-');sb.push_back(' ');revisedcoor.push_back(-2);}
					revisedseq=revisedseq.substr(0,seq_index)+sa+revisedseq.substr(seq_index);
					revisedbq=revisedbq.substr(0,seq_index)+sb+revisedbq.substr(seq_index);
					seq_index=seq_index+itemp;
					coor_index=coor_index+itemp;
				}
				else if(cigar[i]=='M')//match or mismatch
				{
					for(j=0;j<itemp;j++) revisedcoor.push_back(coor_index+j);
					seq_index=seq_index+itemp;
					coor_index=coor_index+itemp;
				}
				else {cout<<"error:"<<cigar<<endl;exit(0);}

				//cout<<cigar<<""\t<<i<<"\t"<<cigar_str[i]<<"\t"<<stemp<<"\t"<<seq_index+itemp<<"\t"<<start<<"\t"<<seq_index<<endl;
			}
		}
	}

	/*cout<<"read "<<startpos<<" "<<cigar<<endl<<seq<<endl<<bq<<endl<<endl;
	cout<<revisedseq<<endl<<revisedbq<<endl;
	for(i=0;i<revisedseq.size();i++) cout<<revisedseq[i]<<" "<<revisedbq[i]<<" "<<revisedcoor[i]<<endl;
	cout<<endl;*/
}

double ChangeQualToProb(const char &qual)
{
	int ia=(int)(qual-33);
	return pow(10,(double)(-ia)/10);
}

void FillContigInfor(const int snvpos,const char snvref,const vector<char> contigNTs,map<int,PosReadsInfor> &pos2Readsinfo)
{
	PosReadsInfor tmp_infor;

	if(pos2Readsinfo.find(snvpos)==pos2Readsinfo.end())//initialize
	{
		tmp_infor.ref=snvref;
		tmp_infor.filterout=false;
		tmp_infor.Qual_set.resize(5);
		tmp_infor.InputQual_set.resize(5);
		tmp_infor.ChIP_rawNo.resize(5);
		tmp_infor.Input_rawNo.resize(5);
		for(int i=0;i<5;i++) {tmp_infor.ChIP_rawNo[i]=0;tmp_infor.Input_rawNo[i]=0;}

		vector<char> vchtemp(contigNTs);
		sort(vchtemp.begin(),vchtemp.end());
		vchtemp.erase(unique(vchtemp.begin(),vchtemp.end()),vchtemp.end());
		tmp_infor.fermiNTs=vchtemp;
/*cout<<"snvas: "<<snvpos<<" "<<snvref<<endl;
for(int i=0;i<contigNTs.size();i++) cout<<contigNTs[i];cout<<endl;
for(int i=0;i<vchtemp.size();i++) cout<<vchtemp[i];cout<<endl;*/
		pos2Readsinfo[snvpos]=tmp_infor;
	}
	else {cout<<"exist snvpos: "<<snvpos<<endl;exit(0);}
}

void FillChIPraw(const int snvpos,const char readnt,map<int,PosReadsInfor> &pos2Readsinfo)
{
	PosReadsInfor tmp_infor;

	if(pos2Readsinfo.find(snvpos)==pos2Readsinfo.end()) {cout<<"not exist snvpos: "<<snvpos<<endl;exit(0);}

	//
	if(readnt=='A')
	{
		pos2Readsinfo[snvpos].ChIP_rawNo[0]++;
	}
	else if(readnt=='C')
	{
		pos2Readsinfo[snvpos].ChIP_rawNo[1]++;
	}
	else if(readnt=='G')
	{
		pos2Readsinfo[snvpos].ChIP_rawNo[2]++;
	}
	else if(readnt=='T')
	{
		pos2Readsinfo[snvpos].ChIP_rawNo[3]++;
	}
	else if(readnt=='N')
	{
		pos2Readsinfo[snvpos].ChIP_rawNo[4]++;
	}
	else {cout<<"wrong nucleotide: "<<readnt<<endl;exit(0);}
}

void FillControlraw(const int snvpos,const char readnt,map<int,PosReadsInfor> &pos2Readsinfo)
{
	PosReadsInfor tmp_infor;

	if(pos2Readsinfo.find(snvpos)==pos2Readsinfo.end()) {cout<<"not exist snvpos: "<<snvpos<<endl;exit(0);}

	//
	if(readnt=='A')
	{
		pos2Readsinfo[snvpos].Input_rawNo[0]++;
	}
	else if(readnt=='C')
	{
		pos2Readsinfo[snvpos].Input_rawNo[1]++;
	}
	else if(readnt=='G')
	{
		pos2Readsinfo[snvpos].Input_rawNo[2]++;
	}
	else if(readnt=='T')
	{
		pos2Readsinfo[snvpos].Input_rawNo[3]++;
	}
	else if(readnt=='N')
	{
		pos2Readsinfo[snvpos].Input_rawNo[4]++;
	}
	else {cout<<"wrong nucleotide: "<<readnt<<endl;exit(0);}
}

void FillChIPQualInfor(const int snvpos,const char readnt,const char readbq,map<int,PosReadsInfor> &pos2Readsinfo)
{
	double qual_proberror=ChangeQualToProb(readbq);

	//if((double)qual_proberror>0.2) return;//if prob of qual_error > 0.2, discard the base

	PosReadsInfor tmp_infor;

	if(pos2Readsinfo.find(snvpos)==pos2Readsinfo.end()) {cout<<"not exist snvpos: "<<snvpos<<endl;exit(0);}

	//
	if(readnt=='A')
	{
		pos2Readsinfo[snvpos].Qual_set[0].push_back(qual_proberror);
	}
	else if(readnt=='C')
	{
		pos2Readsinfo[snvpos].Qual_set[1].push_back(qual_proberror);
	}
	else if(readnt=='G')
	{
		pos2Readsinfo[snvpos].Qual_set[2].push_back(qual_proberror);
	}
	else if(readnt=='T')
	{
		pos2Readsinfo[snvpos].Qual_set[3].push_back(qual_proberror);
	}
	else if(readnt=='N')
	{
		pos2Readsinfo[snvpos].Qual_set[4].push_back(qual_proberror);
	}
	else {cout<<"wrong nucleotide: "<<readnt<<endl;exit(0);}
}

void FillControlQualInfor(const int snvpos,const char readnt,const char readbq,map<int,PosReadsInfor> &pos2Readsinfo)
{
	double qual_proberror=ChangeQualToProb(readbq);

	//if((double)qual_proberror>0.2) return;//if prob of qual_error > 0.2, discard the base

	PosReadsInfor tmp_infor;

	if(pos2Readsinfo.find(snvpos)==pos2Readsinfo.end()) {cout<<"not exist snvpos1: "<<snvpos<<endl;exit(0);}

	//
	if(readnt=='A')
	{
		pos2Readsinfo[snvpos].InputQual_set[0].push_back(qual_proberror);
	}
	else if(readnt=='C')
	{
		pos2Readsinfo[snvpos].InputQual_set[1].push_back(qual_proberror);
	}
	else if(readnt=='G')
	{
		pos2Readsinfo[snvpos].InputQual_set[2].push_back(qual_proberror);
	}
	else if(readnt=='T')
	{
		pos2Readsinfo[snvpos].InputQual_set[3].push_back(qual_proberror);
	}
	else if(readnt=='N')
	{
		pos2Readsinfo[snvpos].InputQual_set[4].push_back(qual_proberror);
	}
	else {cout<<"wrong nucleotide1: "<<readnt<<endl;exit(0);}
}

bool ConsistentWithContig(const vector<int> &readscoor,const string &readsseq,vector<vector<int> > &contigcoor,const vector<string> &contigseq)
{
	int i,j;
	bool btemp=false;
	vector<int>::iterator it;

	for(i=0;i<contigcoor.size();i++)
	{
		it=lower_bound(contigcoor[i].begin(),contigcoor[i].end(),readscoor[0]);
//cout<<"consistent: "<<i<<" "<<contigcoor[i][0]<<"-"<<contigcoor[i][contigcoor[i].size()-1]<<" "<<readscoor[0]<<endl;
//cout<<contigseq[i]<<endl<<readsseq<<endl;
		if(it==contigcoor[i].end() || *it!=readscoor[0])//not found
		{
			continue;
		}
		else//found
		{
			int index=(int)(it-contigcoor[i].begin());
//cout<<index<<endl;
			bool b1=true;
			for(j=0;j<readscoor.size();j++)
			{
				if(readscoor[j]!=contigcoor[i][index+j] || readsseq[j]!=contigseq[i][index+j]) {b1=false;break;}
			}
//cout<<b1<<endl;
			if(b1) {btemp=true;break;}
		}
	}

	return btemp;
}


void OutputVcfResultHasInput(const string outputfile,const string regionchr,map<int,PosReadsInfor> &pos2Readsinfo)
{
	int i;

	//output
	ofstream os(outputfile.c_str(),ofstream::app);
	if(!os) {cout<<"error open output file: "<<endl;exit(0);}

	map<int,PosReadsInfor>::iterator miPtmp;

	for(miPtmp=pos2Readsinfo.begin();miPtmp!=pos2Readsinfo.end();++miPtmp)//for each position
	{
		PosReadsInfor myReadsInfor(miPtmp->second);

		if(myReadsInfor.filterout || !myReadsInfor.hasfermiinfor) continue;

		if(myReadsInfor.type=="homo")
		{
			os<<regionchr<<"\t"<<miPtmp->first<<"\t.\t"<<myReadsInfor.ref;
			os<<"\t"<<NTindex2char(myReadsInfor.top1ntindex)<<"\t.\t.";
			os<<"\tMinBIC_model:"<<myReadsInfor.type<<";raw_depth_ChIP:"<<myReadsInfor.ChIP_rawNo[0]+myReadsInfor.ChIP_rawNo[1]+myReadsInfor.ChIP_rawNo[2]+myReadsInfor.ChIP_rawNo[3]+myReadsInfor.ChIP_rawNo[4]
			  <<";raw_depth_input:"<<myReadsInfor.Input_rawNo[0]+myReadsInfor.Input_rawNo[1]+myReadsInfor.Input_rawNo[2]+myReadsInfor.Input_rawNo[3]+myReadsInfor.Input_rawNo[4]
			  <<";DP_ChIP:"<<myReadsInfor.Qual_set[0].size()+myReadsInfor.Qual_set[1].size()+myReadsInfor.Qual_set[2].size()+myReadsInfor.Qual_set[3].size()+myReadsInfor.Qual_set[4].size()
			  <<";DP_input:"<<myReadsInfor.InputQual_set[0].size()+myReadsInfor.InputQual_set[1].size()+myReadsInfor.InputQual_set[2].size()+myReadsInfor.InputQual_set[3].size()+myReadsInfor.InputQual_set[4].size()<<";fermiNTs:";
			for(i=0;i<myReadsInfor.fermiNTs.size();i++) os<<myReadsInfor.fermiNTs[i];
//cout<<";top1:"<<myReadsInfor.Qual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2:"<<myReadsInfor.Qual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
	//<<";top1input:"<<myReadsInfor.InputQual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2input:"<<myReadsInfor.InputQual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)<<endl;
			os<<";top1:"<<myReadsInfor.Qual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)
			  <<";top1input:"<<myReadsInfor.InputQual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)
			  <<";top1raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)
			  <<";top1inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)
			  <<";lnL_homo_major:"<<myReadsInfor.lnL_homo_majar<<";lnL_homo_minor:"<<myReadsInfor.lnL_homo_minor<<";lnL_heter_noAS:"<<myReadsInfor.lnL_heter_noAS<<";lnL_heter_AS:"<<myReadsInfor.lnL_heter_AS
			  <<";BIC_homo_major:"<<myReadsInfor.BIC_homo_majar<<";BIC_homo_minor:"<<myReadsInfor.BIC_homo_minor<<";BIC_heter_noAS:"<<myReadsInfor.BIC_heter_noAS<<";BIC_heter_AS:"<<myReadsInfor.BIC_heter_AS
			  <<";GQ_homo:"<<myReadsInfor.GQ_homo_majar<<";GQ_heter_noAS:"<<myReadsInfor.GQ_heter_noAS<<";GQ_heter_AS:"<<myReadsInfor.GQ_heter_AS<<";GQ_heter_ASsig:"<<myReadsInfor.GQ_heter_ASsig<<";Allele_ratio_heter_AS:"<<myReadsInfor.heter_AS_alleleratio;
			os<<"\tGT\t1|1"<<endl;
		}
		else if(myReadsInfor.type=="heter_noAS" || myReadsInfor.type=="heter_AS")
		{
			if(NTindex2char(myReadsInfor.top2ntindex)=='N') continue;
			os<<regionchr<<"\t"<<miPtmp->first<<"\t.\t"<<myReadsInfor.ref;

			if(NTindex2char(myReadsInfor.top1ntindex)==myReadsInfor.ref)
			{
				os<<"\t"<<NTindex2char(myReadsInfor.top2ntindex)<<"\t.\t.";
				os<<"\tMinBIC_model:"<<myReadsInfor.type<<";raw_depth_ChIP:"<<myReadsInfor.ChIP_rawNo[0]+myReadsInfor.ChIP_rawNo[1]+myReadsInfor.ChIP_rawNo[2]+myReadsInfor.ChIP_rawNo[3]+myReadsInfor.ChIP_rawNo[4]
				  <<";raw_depth_input:"<<myReadsInfor.Input_rawNo[0]+myReadsInfor.Input_rawNo[1]+myReadsInfor.Input_rawNo[2]+myReadsInfor.Input_rawNo[3]+myReadsInfor.Input_rawNo[4]
				  <<";DP_ChIP:"<<myReadsInfor.Qual_set[0].size()+myReadsInfor.Qual_set[1].size()+myReadsInfor.Qual_set[2].size()+myReadsInfor.Qual_set[3].size()+myReadsInfor.Qual_set[4].size()
				  <<";DP_input:"<<myReadsInfor.InputQual_set[0].size()+myReadsInfor.InputQual_set[1].size()+myReadsInfor.InputQual_set[2].size()+myReadsInfor.InputQual_set[3].size()+myReadsInfor.InputQual_set[4].size()<<";fermiNTs:";
				for(i=0;i<myReadsInfor.fermiNTs.size();i++) os<<myReadsInfor.fermiNTs[i];
				os<<";top1:"<<myReadsInfor.Qual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2:"<<myReadsInfor.Qual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1input:"<<myReadsInfor.InputQual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2input:"<<myReadsInfor.InputQual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)<<";top2raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top2ntindex]<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)<<";top2inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top2ntindex]<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";lnL_homo_major:"<<myReadsInfor.lnL_homo_majar<<";lnL_homo_minor:"<<myReadsInfor.lnL_homo_minor<<";lnL_heter_noAS:"<<myReadsInfor.lnL_heter_noAS<<";lnL_heter_AS:"<<myReadsInfor.lnL_heter_AS
				  <<";BIC_homo_major:"<<myReadsInfor.BIC_homo_majar<<";BIC_homo_minor:"<<myReadsInfor.BIC_homo_minor<<";BIC_heter_noAS:"<<myReadsInfor.BIC_heter_noAS<<";BIC_heter_AS:"<<myReadsInfor.BIC_heter_AS
				  <<";GQ_homo:"<<myReadsInfor.GQ_homo_majar<<";GQ_heter_noAS:"<<myReadsInfor.GQ_heter_noAS<<";GQ_heter_AS:"<<myReadsInfor.GQ_heter_AS<<";GQ_heter_ASsig:"<<myReadsInfor.GQ_heter_ASsig<<";Allele_ratio_heter_AS:"<<myReadsInfor.heter_AS_alleleratio;
				os<<"\tGT\t0|1"<<endl;
			}
			else if(NTindex2char(myReadsInfor.top2ntindex)==myReadsInfor.ref)
			{
				os<<"\t"<<NTindex2char(myReadsInfor.top1ntindex)<<"\t.\t.";
				os<<"\tMinBIC_model:"<<myReadsInfor.type<<";raw_depth_ChIP:"<<myReadsInfor.ChIP_rawNo[0]+myReadsInfor.ChIP_rawNo[1]+myReadsInfor.ChIP_rawNo[2]+myReadsInfor.ChIP_rawNo[3]+myReadsInfor.ChIP_rawNo[4]
				  <<";raw_depth_input:"<<myReadsInfor.Input_rawNo[0]+myReadsInfor.Input_rawNo[1]+myReadsInfor.Input_rawNo[2]+myReadsInfor.Input_rawNo[3]+myReadsInfor.Input_rawNo[4]
				  <<";DP_ChIP:"<<myReadsInfor.Qual_set[0].size()+myReadsInfor.Qual_set[1].size()+myReadsInfor.Qual_set[2].size()+myReadsInfor.Qual_set[3].size()+myReadsInfor.Qual_set[4].size()
				  <<";DP_input:"<<myReadsInfor.InputQual_set[0].size()+myReadsInfor.InputQual_set[1].size()+myReadsInfor.InputQual_set[2].size()+myReadsInfor.InputQual_set[3].size()+myReadsInfor.InputQual_set[4].size()<<";fermiNTs:";
				for(i=0;i<myReadsInfor.fermiNTs.size();i++) os<<myReadsInfor.fermiNTs[i];
				os<<";top1:"<<myReadsInfor.Qual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2:"<<myReadsInfor.Qual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1input:"<<myReadsInfor.InputQual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2input:"<<myReadsInfor.InputQual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)<<";top2raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top2ntindex]<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)<<";top2inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top2ntindex]<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";lnL_homo_major:"<<myReadsInfor.lnL_homo_majar<<";lnL_homo_minor:"<<myReadsInfor.lnL_homo_minor<<";lnL_heter_noAS:"<<myReadsInfor.lnL_heter_noAS<<";lnL_heter_AS:"<<myReadsInfor.lnL_heter_AS
				  <<";BIC_homo_major:"<<myReadsInfor.BIC_homo_majar<<";BIC_homo_minor:"<<myReadsInfor.BIC_homo_minor<<";BIC_heter_noAS:"<<myReadsInfor.BIC_heter_noAS<<";BIC_heter_AS:"<<myReadsInfor.BIC_heter_AS
				  <<";GQ_homo:"<<myReadsInfor.GQ_homo_majar<<";GQ_heter_noAS:"<<myReadsInfor.GQ_heter_noAS<<";GQ_heter_AS:"<<myReadsInfor.GQ_heter_AS<<";GQ_heter_ASsig:"<<myReadsInfor.GQ_heter_ASsig<<";Allele_ratio_heter_AS:"<<myReadsInfor.heter_AS_alleleratio;
				os<<"\tGT\t1|0"<<endl;
			}
			else
			{
				os<<"\t"<<NTindex2char(myReadsInfor.top1ntindex)<<","<<NTindex2char(myReadsInfor.top2ntindex)<<"\t.\t.";
				os<<"\tMinBIC_model:"<<myReadsInfor.type<<";raw_depth_ChIP:"<<myReadsInfor.ChIP_rawNo[0]+myReadsInfor.ChIP_rawNo[1]+myReadsInfor.ChIP_rawNo[2]+myReadsInfor.ChIP_rawNo[3]+myReadsInfor.ChIP_rawNo[4]
				  <<";raw_depth_input:"<<myReadsInfor.Input_rawNo[0]+myReadsInfor.Input_rawNo[1]+myReadsInfor.Input_rawNo[2]+myReadsInfor.Input_rawNo[3]+myReadsInfor.Input_rawNo[4]
				  <<";DP_ChIP:"<<myReadsInfor.Qual_set[0].size()+myReadsInfor.Qual_set[1].size()+myReadsInfor.Qual_set[2].size()+myReadsInfor.Qual_set[3].size()+myReadsInfor.Qual_set[4].size()
				  <<";DP_input:"<<myReadsInfor.InputQual_set[0].size()+myReadsInfor.InputQual_set[1].size()+myReadsInfor.InputQual_set[2].size()+myReadsInfor.InputQual_set[3].size()+myReadsInfor.InputQual_set[4].size()<<";fermiNTs:";
				for(i=0;i<myReadsInfor.fermiNTs.size();i++) os<<myReadsInfor.fermiNTs[i];
				os<<";top1:"<<myReadsInfor.Qual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2:"<<myReadsInfor.Qual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1input:"<<myReadsInfor.InputQual_set[myReadsInfor.top1ntindex].size()<<NTindex2char(myReadsInfor.top1ntindex)<<";top2input:"<<myReadsInfor.InputQual_set[myReadsInfor.top2ntindex].size()<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)<<";top2raw:"<<myReadsInfor.ChIP_rawNo[myReadsInfor.top2ntindex]<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";top1inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top1ntindex]<<NTindex2char(myReadsInfor.top1ntindex)<<";top2inputraw:"<<myReadsInfor.Input_rawNo[myReadsInfor.top2ntindex]<<NTindex2char(myReadsInfor.top2ntindex)
				  <<";lnL_homo_major:"<<myReadsInfor.lnL_homo_majar<<";lnL_homo_minor:"<<myReadsInfor.lnL_homo_minor<<";lnL_heter_noAS:"<<myReadsInfor.lnL_heter_noAS<<";lnL_heter_AS:"<<myReadsInfor.lnL_heter_AS
				  <<";BIC_homo_major:"<<myReadsInfor.BIC_homo_majar<<";BIC_homo_minor:"<<myReadsInfor.BIC_homo_minor<<";BIC_heter_noAS:"<<myReadsInfor.BIC_heter_noAS<<";BIC_heter_AS:"<<myReadsInfor.BIC_heter_AS
				  <<";GQ_homo:"<<myReadsInfor.GQ_homo_majar<<";GQ_heter_noAS:"<<myReadsInfor.GQ_heter_noAS<<";GQ_heter_AS:"<<myReadsInfor.GQ_heter_AS<<";GQ_heter_ASsig:"<<myReadsInfor.GQ_heter_ASsig<<";Allele_ratio_heter_AS:"<<myReadsInfor.heter_AS_alleleratio;
				os<<"\tGT\t1|2"<<endl;
			}
		}
	}
}

void OutputVcfResultHasInput_header(const string outputfile,char *argv[])
{
	int i;
	time_t t = time(0);   // get time now
	struct tm * now = localtime( & t );

	//output
	ofstream os(outputfile.c_str());
	if(!os) {cout<<"error open output file: "<<endl;exit(0);}
	os<<"##fileformat=VCFv4.1"<<endl;
	os<<"##fileDate="<<(now->tm_year + 1900)<<(now->tm_mon + 1)<<now->tm_mday<<endl;
	os<<"##source=SNVAS_V0.1"<<endl;
	os<<"##Program_Args: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<endl;
	os<<"##INFO=<ID=MinBIC_model,Number=.,Type=String,Description=\"Model with minimum BIC value\">"<<endl;
	os<<"##INFO=<ID=raw_depth_ChIP,Number=1,Type=Integer,Description=\"Raw read depth in ChIP-seq data\">"<<endl;
	os<<"##INFO=<ID=raw_depth_input,Number=1,Type=Integer,Description=\"Raw read depth in input data\">"<<endl;
	os<<"##INFO=<ID=DP_ChIP,Number=1,Type=Integer,Description=\"Read depth in ChIP-seq data; some reads may have been filtered\">"<<endl;
	os<<"##INFO=<ID=DP_input,Number=1,Type=Integer,Description=\"Read depth in input data; some reads may have been filtered\">"<<endl;
	os<<"##INFO=<ID=fermiNTs,Number=.,Type=String,Description=\"Nucleotides from the genotype information of fermi assembly result\">"<<endl;
	os<<"##INFO=<ID=top1,Number=.,Type=String,Description=\"Read depth of top1 nucleotide in ChIP-seq data; some reads may have been filtered\">"<<endl;
	os<<"##INFO=<ID=top2,Number=.,Type=String,Description=\"Read depth of top2 nucleotide in ChIP-seq data; some reads may have been filtered\">"<<endl;
	os<<"##INFO=<ID=top1input,Number=.,Type=String,Description=\"Read depth of top1 nucleotide in input data; some reads may have been filtered\">"<<endl;
	os<<"##INFO=<ID=top2input,Number=.,Type=String,Description=\"Read depth of top2 nucleotide in input data; some reads may have been filtered\">"<<endl;
	os<<"##INFO=<ID=top1raw,Number=.,Type=Integer,Description=\"Read depth of top1 nucleotide in raw ChIP-seq data\">"<<endl;
	os<<"##INFO=<ID=top2raw,Number=.,Type=Integer,Description=\"Read depth of top2 nucleotide in raw ChIP-seq data\">"<<endl;
	os<<"##INFO=<ID=top1inputraw,Number=.,Type=Integer,Description=\"Read depth of top1 nucleotide in raw input data\">"<<endl;
	os<<"##INFO=<ID=top2inputraw,Number=.,Type=Integer,Description=\"Read depth of top1 nucleotide in raw input data\">"<<endl;
	os<<"##INFO=<ID=lnL_homo_major,Number=1,Type=Float,Description=\"Log(e) scaled genotype likelihoods of homozygous with major allele model\">"<<endl;
	os<<"##INFO=<ID=lnL_homo_minor,Number=1,Type=Float,Description=\"Log(e) scaled genotype likelihoods of homozygous with minor allele model\">"<<endl;
	os<<"##INFO=<ID=lnL_heter_noAS,Number=1,Type=Float,Description=\"Log(e) scaled genotype likelihoods of heterozygous with no allele-specific model\">"<<endl;
	os<<"##INFO=<ID=lnL_heter_AS,Number=1,Type=Float,Description=\"Log(e) scaled genotype likelihoods of heterozygous with allele-specific model\">"<<endl;
	os<<"##INFO=<ID=BIC_homo_major,Number=1,Type=Float,Description=\"BIC value of homozygous with major allele model\">"<<endl;
	os<<"##INFO=<ID=BIC_homo_minor,Number=1,Type=Float,Description=\"BIC value of homozygous with minor allele model\">"<<endl;
	os<<"##INFO=<ID=BIC_heter_noAS,Number=1,Type=Float,Description=\"BIC value of heterozygous with no allele-specific model\">"<<endl;
	os<<"##INFO=<ID=BIC_heter_AS,Number=1,Type=Float,Description=\"BIC value of heterozygous with allele-specific model\">"<<endl;
	os<<"##INFO=<ID=GQ_homo,Number=1,Type=Float,Description=\"Genotype quality of homozygous with major allele model\">"<<endl;
	os<<"##INFO=<ID=GQ_heter_noAS,Number=1,Type=Float,Description=\"Genotype quality of heterozygous with no allele-specific model\">"<<endl;
	os<<"##INFO=<ID=GQ_heter_AS,Number=1,Type=Float,Description=\"Genotype quality of heterozygous with allele-specific model\">"<<endl;
	os<<"##INFO=<ID=GQ_heter_ASsig,Number=1,Type=Float,Description=\"Genotype quality of allele-specific significance compared with no allele-specific model\">"<<endl;
	os<<"##INFO=<ID=Allele_ratio_heter_AS,Number=1,Type=Float,Description=\"Estimated allele ratio of heterozygous with allele-specific model\">"<<endl;
	//os<<"##reference=GRch37/hg19"<<endl;
	os<<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"<<endl;
	os<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"<<endl;
}
