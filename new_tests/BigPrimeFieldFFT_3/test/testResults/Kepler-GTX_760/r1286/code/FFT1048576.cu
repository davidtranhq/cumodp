/* This file is part of the CUMODP library

    CUMODP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUMODP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUMODP.  If not, see <http://www.gnu.org/licenses/>.

    Copyright: Sardar Anisul Haque <shaque4@uwo.ca>
               Xin Li <xli.software@gmail.com>
               Farnam Mansouri <mansouri.farnam@gmail.com>
               Davood Mohajerani <dmohajer@uwo.ca>
               Marc Moreno Maza  <moreno@csd.uwo.ca>
               Wei Pan <wei.pan@intel.com>
               Ning Xie <nxie6@csd.uwo.ca>
*/


#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include<math.h>

//N the FFT size, where each element contains eight unsigned integers.
#define N 1048576
//BN block number
#define BLOCK_DIM 512
//TN thread number in a block
#define THREAD_DIM 512 
//input number count
#define INC 256
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//sqrt(R) 3037000502
#define SQRTR 3037000502
//ULMAX=2^64-1
#define ULMAX 18446744073709551615

__constant__ short ind2[8]={0,0,0,0,4,4,4,4};
__constant__ short ind3[8]={0,4,2,6,0,4,2,6};
__constant__ short ind4[8]={0,4,2,6,1,5,3,7};
__constant__ short ind5[16]={0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
__constant__ unsigned long W8[15][8]={331066363573078852u, 1926294779415465329u, 3360946100581367489u, 5520011298065530166u, 4832614680308292602u, 8402438242641148880u, 2486701145004827901u, 1215222897245048005u, 6920921854211038777u, 9196291859899298732u, 4110928224385435654u, 396332478198426435u, 7935859707626518415u, 4857645604960533584u, 5862468289385444427u, 9105649095669603570u, 3345465055005335276u, 5514352822140455844u, 6208596698165918434u, 3765508936539829447u, 8833047145598725640u, 8452480178648776477u, 6223675675889128755u, 6091304989396919910u, 9210848751572284741u, 1645667176001997922u, 5941128666689277807u, 6768772028034979755u, 4304177272275830547u, 6104712663442007465u, 1954368461670077480u, 6432833914249218714u, 350789192213411171u, 5484312567191208149u, 409536839932512796u, 4075900149429376385u, 5366875984871234375u, 8713960611431090699u, 7737743610771967836u, 4469224371143296702u, 6730078996384360032u, 1055782569299382065u, 2818483389281035947u, 6249192910925906484u, 4406661415423507943u, 7522848620476680275u, 3980191202297482070u, 2108216364668858652u, 4477472493274660184u, 8475805979162198949u, 3965693063579606190u, 7514110372369303684u, 5716477088365326795u, 1805006940055860596u, 1209019385223250294u, 7053302953741260526u, 6476678624121558179u, 3216147715999829916u, 4168578045703223305u, 5665582343424371252u, 4802751892420194639u, 7507476765166618649u, 794581407066667206u, 4885142567208165039u, 2723846702617923460u, 4129024823037556627u, 508195973098175290u, 7381493553839972652u, 5149649432978519210u, 5596777036269387984u, 8665828768845738479u, 1575627868412583125u, 1680992761645806539u, 166450404224435382u, 592663973356561630u, 2660631676326679188u, 8420781954126180406u, 2166425954682284818u, 1408165695181296747u, 667245596591619233u, 8888111600302919098u, 5096574789189059728u, 2602858914188694775u, 1490780946449946389u, 6488129570674181190u, 7714154710699913078u, 1086155969382280547u, 4028545677085900533u, 1871173134326118337u, 6362464371896225965u, 8698229383614300186u, 2005452311774545610u, 7556559507024761731u, 5864043929342978154u, 7933436335374648226u, 5827741928560708899u, 3969563368936755191u, 2171179263218223366u, 2323713704096704824u, 7941929740108455451u, 1224342742045833268u, 972714779852171763u, 2469969374900767744u, 925611663966851379u, 4017136705252831872u, 8246845268085542850u, 5023523124979160867u, 8559956080252177681u, 4556280367171742838u, 3567089280497015080u, 8204228649088379144u, 4528147079632326543u, 3068873176292580126u, 753173231159502005u, 4879546884467215097u, 6201146281074609078u, 2268976475187260030u, 3981797643749168921u, 7597251879112066229u, 6689850317535721737u};
__constant__ unsigned long W12[15][8]={6368537335278011409u, 6559905493361745650u, 1611630010065826713u, 5205069994608960148u, 2852366324136299255u, 4741025511079831292u, 4495713848194575320u, 2464458433847228250u, 43390553589054774u, 5689918209556578671u, 8845941305773747342u, 2840602737422012834u, 813638098239981480u, 314969166600807247u, 7711698714301782605u, 4708052650754658993u, 5962349939731889267u, 9187502500840055454u, 6223916889772168053u, 8200297666652537596u, 7243141932147185220u, 5827178249755404095u, 1233738160859353199u, 9135786083964470993u, 7137107724604502333u, 2979207810221979654u, 3328860960606119836u, 4152055445595025560u, 790279036532063445u, 2243751449889058377u, 495802657723411899u, 7620665221889874093u, 3608334579833484687u, 9091502745705873137u, 9068759319444887825u, 2982042472509517038u, 5851203093862561394u, 4689164771082797939u, 1424759403320645141u, 4411089293600739857u, 6907127599491296095u, 8552897415997102153u, 8679320014683906997u, 2223828807493710814u, 858649556069331575u, 1957943790061739103u, 5549400026053976849u, 2098463206672069964u, 4646819650281801505u, 6241670829831344849u, 7429385380966888965u, 419444289839068823u, 4130637287415002184u, 7908303496741499390u, 7328076545823401678u, 8906856790568058462u, 4266425184611474885u, 2745888046073582606u, 4460272313717620391u, 1534602007411134351u, 7389417065503563175u, 2496816640994598707u, 8340578439454199486u, 3540970978625599238u, 6593137379340673094u, 938822237580288165u, 6173761731301493950u, 6345868649901054830u, 33203932990054256u, 8285674207975771116u, 391365543956029330u, 5765341871901814199u, 8750005390432559331u, 1004177955874976073u, 6526666171281335965u, 8919526347441079930u, 4587145129984655982u, 2120995648749380660u, 6807125005630921771u, 1404319894383668875u, 9174822346926590820u, 8683903108746749848u, 518805224601400271u, 8263206596169269816u, 3992802788711728259u, 4526014523036653276u, 8216527007903401430u, 3643499083168944286u, 6572700652819984097u, 5654349319026678090u, 2945958918231407402u, 8315817940191554115u, 3903961909722922192u, 8751115554095427145u, 1720087832858432991u, 2017660982864315227u, 4932201570704995548u, 2441758989690403245u, 6315068009463966664u, 8141367093872374404u, 2386952287665435440u, 8630517035780209587u, 3791877247993007898u, 4994465731821916226u, 5005814032629966111u, 7629929215843812243u, 5751181848221439447u, 4576230681179567819u, 8125534392733532035u, 5560926371000433936u, 1888459826216152816u, 3747509517407033084u, 7980171858755715406u, 4008475782978845368u, 8932281156699250062u, 1016300905947845519u, 6522461459654622787u, 3266045571977616788u, 6955253718610834946u, 7667560955258982823u};
__constant__ unsigned long W16[15][8]={3006370066595279709u, 3782262637662174695u, 4774667134945453143u, 7390498005841228038u, 2430596464728967792u, 8280774632104889468u, 6623214708194984944u, 2614500752162691591u, 8816938155464702529u, 3706565079792625720u, 5700587314643272992u, 9081268561225285341u, 3742668007806798414u, 219006490403163241u, 6555121409556599293u, 619595731419319046u, 7125943389931267119u, 9102905563943801299u, 5126318784044180744u, 1840476394856931165u, 4734522244195204554u, 4446720878807858978u, 1465735233737441534u, 9141346959181243938u, 7238850373747071574u, 3421722769799276146u, 3605960622287093330u, 5615526701043691871u, 5987135835531705030u, 4592117420273897548u, 4297308844345453283u, 1582988213281168002u, 1958652214492452177u, 3866660406659619836u, 1748812593028471511u, 9156565293891688833u, 954769817498911500u, 3304769054508767064u, 8023664163237704382u, 1358174849470882936u, 4024480772940697201u, 3138121140116141962u, 5319250646176980025u, 5769776673728266447u, 7824873459455557320u, 1471083240729778152u, 1512398299696802745u, 4815047575398614010u, 4960037255635211982u, 3222259664270834564u, 145757677427459038u, 3866062095924381007u, 9139164790258791352u, 2122634499199573707u, 5705627251702193578u, 1007440913992157340u, 7981465692565196598u, 758164073934943750u, 8092267366146177062u, 6328413749839545375u, 1409518216045347123u, 3745783276278595597u, 5016761400313803385u, 244667108070889342u, 7458378155224591849u, 6700850949988090323u, 8429569059329231707u, 3900555209317066428u, 2580064796542557327u, 2330943047179410679u, 4678117607704355473u, 8868258420852560284u, 1056428342586890892u, 3583308366490261645u, 373441718173383040u, 8114613508589847156u, 2387658483236314119u, 1951882842630615779u, 4182434194191353116u, 2087439161083770738u, 1295541677614975979u, 8304821556901566758u, 7600505033415929188u, 8545534317505075137u, 5079777146908193381u, 8020024906592432531u, 3012818430327884827u, 8006501689739503680u, 7814992315030354129u, 3274432195841795110u, 6118362530917551209u, 2329169177125719489u, 6544952968060818u, 8503079195181123741u, 3879911002364542365u, 2726159217504441781u, 4081850729939701628u, 3191017481995419158u, 6275953065807032088u, 1370593395923751996u, 3647175153941951850u, 5931784142067835657u, 8898752865083446498u, 8487268107029637522u, 6530612670287071899u, 2817203841761442454u, 4011136235211278128u, 4374832133679005571u, 2837754611208897105u, 1618510605215255457u, 486997559487993043u, 2800884173332122954u, 3749279512132763892u, 8690413992066521869u, 3725041012226324385u, 2771447717726749198u, 2440146651011894402u, 9079794031669681697u, 1477887754829831604u, 6695920399073449841u};
__constant__ unsigned long W20[15][8]={2094657654529411228u, 1874933489908008877u, 6772029229527564241u, 4924424746725594203u, 2669790826289816776u, 6008940768097561304u, 1428899285904723708u, 3266193411217977525u, 138426599768550301u, 6481282715391704207u, 6283415427083676378u, 3197384551042727711u, 2456673898495528529u, 6136269307916091307u, 5630809915478409226u, 8043798999810556290u, 8422966900421295099u, 331151803767016750u, 3490153773066443835u, 1084357015997404523u, 1192442687398231903u, 1440399075933353991u, 4463897891775381539u, 1549404707838784300u, 1241274152157544391u, 2717649403172491102u, 5428170064498949614u, 3819947450412250480u, 4325234897276133104u, 1072482933559798114u, 5074636265416201703u, 8123808161658686761u, 6405673649501298505u, 5922572879083179230u, 8407378843661000777u, 8727087190724518130u, 715758175863965659u, 6800197200322261653u, 7131426094646383322u, 6141168093117979109u, 6899333626602357791u, 238483997890824147u, 3473468640762281280u, 1611881855156301605u, 1207097130233599245u, 2779471946572263353u, 3446067825813260479u, 1684988849176665164u, 6203915632888220116u, 6557690790857481725u, 5021744546528100420u, 7192849881679404157u, 6861603466132154303u, 2179387987843149602u, 7722546021092411873u, 6433208540818943333u, 1918959985148445310u, 62390779226024247u, 2953865821585852732u, 8832275090902445875u, 7643135054298427082u, 8814181381229800156u, 3429096051941704526u, 1369285775060777994u, 5700728876909269292u, 6072125199721242393u, 645810693208198302u, 6371303033002198671u, 6431021097864216073u, 7763090044662680087u, 8547402866115794226u, 8256866257120912622u, 8629409383013415447u, 4646136419943285899u, 3249217443409211714u, 8809472087243418857u, 1241666329040448882u, 3977265503633788372u, 2399032830572181117u, 404428698434593775u, 1034745726355359013u, 5780430645935719130u, 7276479573294730763u, 4638262753731267168u, 2878067424497565156u, 3715539710184518351u, 3013614438228264168u, 8944952554913336114u, 6434142395025400292u, 1945026289917011088u, 7751570819617564489u, 3577660297679292775u, 1887751614217618515u, 3370781657569414041u, 3147791024890462506u, 2739819199782673425u, 225297369991865730u, 7436115953854748242u, 4669596064763193426u, 5609898869083476736u, 6020406955634738404u, 6542335674098563284u, 431840566424992069u, 5323274136466655639u, 5039767721703805154u, 3507127470050419456u, 5467699198250670282u, 4605973736482646u, 7211493787590730750u, 6389218567027858581u, 6369746506002118784u, 1427158368736903301u, 8578683438174047705u, 2830643883343961556u, 3370777893973522806u, 7049347005369424582u, 5416754782040679767u, 7240510424632037619u, 2470560349253849106u, 4882454104227556398u};
__constant__ unsigned long W24[15][8]={7651665413839504047u, 4494057363764474140u, 8486424765903903151u, 3098172693553047908u, 5793710114497482112u, 8806649171840145749u, 5452389321454387529u, 7359030844583626106u, 3143781616945029459u, 3838018006397439943u, 7221435095246147043u, 6719777801270331995u, 7629811503723625374u, 3495464411120386308u, 9060013406001256393u, 783452366975374093u, 1736579118732691671u, 628074136283528523u, 4868835393823627338u, 2237466859000000362u, 5928573644596560429u, 775651498496439375u, 6193596017888085441u, 2856772895245971587u, 6619279187550442103u, 6672797718111533333u, 4585572439994587095u, 2320708279686486197u, 1967162587545328732u, 2286270625854470690u, 1466399497681167039u, 74129947264727246u, 5115677663647077028u, 1149629935203294536u, 8239558042461317190u, 2093058799874099409u, 6253624369706866223u, 7679051254909301482u, 4029554257096547023u, 4484752021829099380u, 1224103195011319710u, 7748363627103633794u, 9045180451380581002u, 2212351111167866412u, 3782058105276009877u, 8695163935942028952u, 8533057948015650743u, 8668062236265483458u, 4252853508362582815u, 4084310759376965869u, 8637491841819748155u, 1122471831225285593u, 4868668863783016133u, 172148462210428924u, 8041342836095768564u, 8137552297942858756u, 4718432596854144008u, 8395476130236506478u, 849258159608788765u, 2385083741171156415u, 1893990415552175293u, 8845564351549673073u, 5296861600082211952u, 5790871676227276177u, 2558199393477890227u, 2286777303798489019u, 7290896120447815229u, 1087736033683488356u, 3233538458063410178u, 8245613835142232393u, 284789579572249726u, 2107522681282643627u, 8796167062570075119u, 7516428514481720367u, 6799511782816783872u, 6380333507059198585u, 1192592460195891196u, 3367423065652605475u, 7113880659446704564u, 1778306093194925568u, 4425320601763114056u, 8040442985044632654u, 4290283958996082457u, 5090995746989133692u, 2416554927487554821u, 4817181417651550341u, 3699462065542763875u, 2271896955554592208u, 3828108457042746435u, 1426979824127340284u, 8728213755202421860u, 5856824756145019556u, 1160591506864909441u, 7966970680263190562u, 6513725650168170156u, 511527149298584934u, 2730052615806903892u, 3886407001682305297u, 2956739463058437375u, 5764178023578064992u, 8641946203099481387u, 7775374034538871454u, 3599115992636436245u, 2191668898145290846u, 1121853962187938628u, 8927501053379065736u, 6673563515517607882u, 4928482465553414419u, 446302744873136398u, 3255731337941165095u, 4091740170832344045u, 6357571692469457015u, 595259163142599786u, 8339829252518716212u, 6422038539823944285u, 1873738540944185089u, 1171201291072818932u, 2474297010996193353u, 3783381178171742213u, 2854015431302544325u};


__device__ void bigAdd(unsigned long *xs, unsigned long *ys, unsigned long *us) 
{  
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  num1=0;
  num2=R-1;
  
  for(i=0;i<=7;i++)
  {
    num1=xs[i]+ys[i]+c;
    if(num1<xs[i]||num1<ys[i]) //there is overflow/truncation
    {
      c=1;
      us[i]=num1+RC;
    }
    else if(num1>=R)
    {  
      c=1;
      us[i]=num1-R;
    }
    else
    {
      c=0;
      us[i]=num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(us[i] != 0)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          us[i] = num2;
      }
      us[pos] = us[pos] - 1;
    }            
    else
    {
      us[0]=ULMAX;
      for(i=1;i<8;i++)
      {
        us[i]=0;
      }
    }  
  }  
}

__device__ void bigAdd2(unsigned long *xs, unsigned long *ys) 
{
  unsigned long us[8];
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  num1=0;
  num2=R-1;
  
  for(i=0;i<=7;i++)
  {
    num1=xs[i]+ys[i]+c;
    if(num1<xs[i]||num1<ys[i]) //there is overflow/truncation
    {
      c=1;
      us[i]=num1+RC;
    }
    else if(num1>=R)
    {  
      c=1;
      us[i]=num1-R;
    }
    else
    {
      c=0;
      us[i]=num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(us[i] != 0)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          us[i] = num2;
      }
      us[pos] = us[pos] - 1;
    }            
    else
    {
      us[0]=ULMAX;
      for(i=1;i<8;i++)
      {
        us[i]=0;
      }
    }  
  }  
  xs[0]=us[0];
  xs[1]=us[1];
  xs[2]=us[2];
  xs[3]=us[3];
  xs[4]=us[4];
  xs[5]=us[5];
  xs[6]=us[6];
  xs[7]=us[7];
}

__device__ void bigSub(unsigned long *xs, unsigned long *ys, unsigned long *us) 
{  
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  num1=0;
  num2=R-1;
  
  for(i=0;i<=7;i++)
  {
    num1=ys[i]+c;
    if(xs[i]<num1) //there is not enough to do subtraction
    {
      c=1;
      us[i]=R-num1+xs[i];
    }
    else
    {
      c=0;
      us[i]=xs[i]-num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(us[i] < num2)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          us[i] = 0;
      }
      us[pos] ++;
    }            
    else
    {
      us[0]=ULMAX;
      for(i=1;i<8;i++)
      {
        us[i]=0;
      }
    }  
  }  
}

__device__ void bigSub2(unsigned long *xs, unsigned long *ys) 
{  
  unsigned long us[8];
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  num1=0;
  num2=R-1;
  
  for(i=0;i<=7;i++)
  {
    num1=ys[i]+c;
    if(xs[i]<num1) //there is not enough to do subtraction
    {
      c=1;
      us[i]=R-num1+xs[i];
    }
    else
    {
      c=0;
      us[i]=xs[i]-num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(us[i] < num2)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          us[i] = 0;
      }
      us[pos] ++;
    }            
    else
    {
      us[0]=ULMAX;
      for(i=1;i<8;i++)
      {
        us[i]=0;
      }
    }  
  }
  xs[0]=us[0];
  xs[1]=us[1];
  xs[2]=us[2];
  xs[3]=us[3];
  xs[4]=us[4];
  xs[5]=us[5];
  xs[6]=us[6];
  xs[7]=us[7];  
}

//0<=sn<=8
__device__ void cyclicShift(unsigned long *xs, short sn)
{
  short i=0,j=0;
  unsigned long ts[8]={0};
  if(sn<=0)
  {
    return;
  }
  if(sn>8)
  {
    return;
  }
  j=8-sn;
  for(i=0;i<sn;i++)
  {
    ts[i]=xs[j++];
  }
  for(i=7-sn;i>=0;i--)
  {
    xs[i+sn]=xs[i];
  }
  for(i=0;i<sn;i++)
  {
    xs[i]=0;
  }
  bigSub2(xs,ts);
}

//[l,h,c]
__device__ void mulLong(unsigned long x, unsigned long y, unsigned long *s)
{
  short x1,y1;
  unsigned long l,h,c,x0,y0,v2,v5,v9,v10,v11,v14,v15,v16,v17,q,t;
  unsigned long a0,a1,b0,b1,c0,c1,c1prime,d0,d1,d2,e0,e1;
  
  if(x<=SQRTR && y<=SQRTR)
  {
    s[0]=x*y;
    s[1]=0;
    s[2]=0;
    return;
  }
  
  x1=(x>=R?1:0);
  x0=(x1>0?x-R:x);
  y1=(y>=R?1:0);
  y0=(y1>0?y-R:y);
  
  v2=x0*y1; //[0,v2,0];
  v5=x1*y0; //[0,v5,0];
  v9=x1*y1; //[0,0,1];
  
  c=v9;
  l=0;
  h=v5+v2;
  h<v5||h<v2?(c=c+1):(c=c);
  c>v9?(h=h+RC):(h=h);
  
  if(x0<=SQRTR&&y0<=SQRTR)
  {
    s[0]=x0*y0;
    s[1]=h;
    s[2]=c;
    return;
  }
  
  //lhc
  //x0*y0
  a1=x0>>32;
  a0=x0-(a1<<32);
  b1=y0>>32;
  b0=y0-(b1<<32);
  
  c0=0;
  c1=a1*b1;
  
  t=a0*b1;
  q=t>>32;
  t=(t-(q<<32))<<32;
  c1+=q;
  c0+=t;  //safe
  
  t=a1*b0;
  q=t>>32;
  t=(t-(q<<32))<<32;
  c1+=q;
  q=c0+t;               //here, is not related to r.
  q<c0||q<t?(c1++):(c1=c1);  //c0=c0+t and carry, safe
  c0=q;
  
  t=a0*b0;
  q=c0+t;
  q<c0||q<t?(c1++):(c1=c1);  //Now we finish [c0,c1]=x0*y0
  c0=q;
  
  c1prime=c1<<1;
  
  c0>=R?(v11=1):(v11=0);
  v11>0?(v10=c0-R):(v10=c0);
  //v12=0;
  
  q=l+v10;  //[l,h,c] + [v10,v11,0]
  q<l||q<v10?(v11=v11+1):(v11=v11);
  q<l||q<v10?(l=q+RC):(l=q);
  if(l>=R)
  {
    l=l-R;
    v11++;
  }
  q=h+v11;
  q<h||q<v11?(c=c+1):(c=c);
  q<h||q<v11?(h=q+RC):(h=q);
  if(h>=R)
  {
    h=h-R;
    c++;
  }
  //v13=0;
  c1prime>=R?(v15=1):(v15=0);
  v15>0?(v14=c1prime-R):(v14=c1prime); //v13=0;
  
  q=h+v14;  //[l,h,c]+[0,v14,v15]
  q<h||q<v14?(c=c+v15+1):(c=c+v15);
  q<h||q<v14?(h=q+RC):(h=q);
  if(h>=R)
  {
    h=h-R;
    c++;
  }
  //[l,h,c]
  
  d1=c1prime>>29;
  d0=c1prime-(d1<<29);
  if(d0>=d1)
  {
    d2=d0-d1;
    e1=d2>>29;
    e0=d2-(e1<<29);
    e0>=e1?(v16=(e0-e1)<<34):(v16=R-(e1<<34)+(e0<<34));
    e0>=e1?(v17=e1+d1):(v17=e1+d1-1);
    /*
    if(e0>=e1)
    {
      v16=(e0-e1)<<34;
      v17=e1+d1;
    }
    else
    {
      v17=e1+d1-1;
      v16=R-(e1<<34)+(e0<<34);
    }
    */
  }
  else
  {
    //d1>d0
    d2=d1-d0;
    e1=d2>>29;
    e0=d2-(e1<<29);
    e0>=e1?(v16=R-((e0-e1)<<34)):(v16=(e1-e0)<<34);
    e0>=e1?(v17=d1-e1-1):(v17=d1-e1);
    /*
    if(e0>=e1)
    {
      v16=R-((e0-e1)<<34);
      v17=d1-e1-1;
    }
    else
    {
      v16=(e1-e0)<<34;
      v17=d1-e1;
    }
    */
  }
  //[l,h,c]-[v16,v17,0]
  //q
  q=0;
  if(l>=v16)
  {
    l=l-v16;
  }
  else
  {
    l=R-v16+l;
    q=1;
  }
  //t
  if(h<q+v17)
  {
    c=c-1;
    h=R-q-v17+h;
  }
  else
  {
    h=h-q-v17;
  }  
  s[0]=l;
  s[1]=h;
  s[2]=c;
}

//store in l0, h0, c0
__device__ void smallAdd(unsigned long *l0, unsigned long *h0, short *c0, unsigned long *l1, unsigned long *h1, unsigned long *c1)
{
  short c=0;
  unsigned long s=0;
  s=*l0+*l1;
  s<*l0 || s<*l1?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  *l0=s;
  
  *h1=*h1+c;  //h1<r<2^64-1. This means no overflow
  s=*h0+*h1;
  s<*h0||s<*h1?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  *h0=s;
  
  *c0=*c0+(short)*c1+c;
}

__device__ void bigMul(unsigned long *xs, unsigned long *ys) 
{
    unsigned long ts1[8];
    unsigned long ts2[8];
    unsigned long rs[3];
    short c0,c1,c2,c3,c4,c5,c6,c7;
    unsigned long l0,l1,l2,l3,l4,l5,l6,l7,h0,h1,h2,h3,h4,h5,h6,h7;
    
    //x0*y0
    mulLong(xs[0],ys[0],rs);
    l0=rs[0];
    h0=rs[1];
    c0=(short)rs[2];
    
    //x0*y1+x1*y0
    mulLong(xs[0],ys[1],rs);    
    l1=rs[0];
    h1=rs[1];
    c1=(short)rs[2];
    mulLong(xs[1],ys[0],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    
    //x0*y2+x1*y1+x2*y0
    mulLong(xs[0],ys[2],rs);    
    l2=rs[0];
    h2=rs[1];
    c2=(short)rs[2];
    mulLong(xs[1],ys[1],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[0],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    
    //x0*y3+x1*y2+x2*y1+x3*y0
    mulLong(xs[0],ys[3],rs);    
    l3=rs[0];
    h3=rs[1];
    c3=(short)rs[2];
    mulLong(xs[1],ys[2],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[1],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[0],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);    
    
    //x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
    mulLong(xs[0],ys[4],rs);    
    l4=rs[0];
    h4=rs[1];
    c4=(short)rs[2];
    mulLong(xs[1],ys[3],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[2],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[1],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[0],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    
    //x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
    mulLong(xs[0],ys[5],rs);    
    l5=rs[0];
    h5=rs[1];
    c5=(short)rs[2];
    mulLong(xs[1],ys[4],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[3],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[2],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[1],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[0],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
 
    //x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
    mulLong(xs[0],ys[6],rs);    
    l6=rs[0];
    h6=rs[1];
    c6=(short)rs[2];
    mulLong(xs[1],ys[5],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[4],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[3],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[2],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[1],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[0],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    
    //x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
    mulLong(xs[0],ys[7],rs);    
    l7=rs[0];
    h7=rs[1];
    c7=(short)rs[2];
    mulLong(xs[1],ys[6],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[5],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[4],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[3],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[2],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[1],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[0],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    
    // (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
    ts1[0]=l0;ts1[1]=h0;ts1[2]=c0;ts1[3]=c1;
    ts1[4]=c2;ts1[5]=c3;ts1[6]=c4;ts1[7]=c5;
    ts2[0]=0;ts2[1]=l1;ts2[2]=h1;ts2[3]=h2;
    ts2[4]=h3;ts2[5]=h4;ts2[6]=h5;ts2[7]=h6;
    bigAdd2(ts1,ts2);
    ts2[0]=0;ts2[1]=0;ts2[2]=l2;ts2[3]=l3;
    ts2[4]=l4;ts2[5]=l5;ts2[6]=l6;ts2[7]=l7;
    bigAdd2(ts1,ts2);
    ts2[0]=c6;ts2[1]=c7;ts2[2]=0;ts2[3]=0;
    ts2[4]=0;ts2[5]=0;ts2[6]=0;ts2[7]=0;
    bigSub2(ts1,ts2);
    ts2[0]=h7;ts2[1]=0;ts2[2]=0;ts2[3]=0;
    ts2[4]=0;ts2[5]=0;ts2[6]=0;ts2[7]=0;
    bigSub2(ts1,ts2);
    
    //(x7*y7)r^6
    mulLong(xs[7],ys[7],rs);
    l6=rs[0];
    h6=rs[1];
    c6=(short)rs[2];
    
    //(x6*y7+x7*y6)r^5
    mulLong(xs[6],ys[7],rs);    
    l5=rs[0];
    h5=rs[1];
    c5=(short)rs[2];
    mulLong(xs[7],ys[6],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    
    //(x5*y7+x6*y6+x7*y5)r^4
    mulLong(xs[5],ys[7],rs);    
    l4=rs[0];
    h4=rs[1];
    c4=(short)rs[2];
    mulLong(xs[6],ys[6],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[5],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);    
    
    //(x4*y7+x5*y6+x6*y5+x7*y4)r^3
    mulLong(xs[4],ys[7],rs);    
    l3=rs[0];
    h3=rs[1];
    c3=(short)rs[2];
    mulLong(xs[5],ys[6],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);    
    mulLong(xs[6],ys[5],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);    
    mulLong(xs[7],ys[4],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);
    
    //(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
    mulLong(xs[3],ys[7],rs);    
    l2=rs[0];
    h2=rs[1];
    c2=(short)rs[2];
    mulLong(xs[4],ys[6],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[5],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[4],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[3],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    
    //(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
    mulLong(xs[2],ys[7],rs);    
    l1=rs[0];
    h1=rs[1];
    c1=(short)rs[2];
    mulLong(xs[3],ys[6],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[5],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[4],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[3],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[2],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    
    //(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
    mulLong(xs[1],ys[7],rs);    
    l0=rs[0];
    h0=rs[1];
    c0=(short)rs[2];
    mulLong(xs[2],ys[6],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[5],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[4],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[3],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[2],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[1],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    
    //(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
    ts2[0]=l0;ts2[1]=h0;ts2[2]=c0;ts2[3]=c1;
    ts2[4]=c2;ts2[5]=c3;ts2[6]=c4;ts2[7]=c5;
    bigSub2(ts1,ts2);
    ts2[0]=0;ts2[1]=l1;ts2[2]=h1;ts2[3]=h2;
    ts2[4]=h3;ts2[5]=h4;ts2[6]=h5;ts2[7]=h6;
    bigSub2(ts1,ts2);
    ts2[0]=0;ts2[1]=0;ts2[2]=l2;ts2[3]=l3;
    ts2[4]=l4;ts2[5]=l5;ts2[6]=l6;ts2[7]=0;
    bigSub2(ts1,ts2);
    ts2[0]=c6;ts2[1]=0;ts2[2]=0;ts2[3]=0;
    ts2[4]=0;ts2[5]=0;ts2[6]=0;ts2[7]=0;
    bigAdd2(ts1,ts2);
    xs[0]=ts1[0];
    xs[1]=ts1[1];
    xs[2]=ts1[2];
    xs[3]=ts1[3];
    xs[4]=ts1[4];
    xs[5]=ts1[5];
    xs[6]=ts1[6];
    xs[7]=ts1[7];
}


__device__ void mulOmega(unsigned long *xs, int pid, unsigned long *ws)
{
  int rp,wp;
  
  if(pid<=0)
  {
    return;
  }
  
  rp=pid>>16; //w^1048576=r  0<=rp<16  caution: depend on N.
  wp=pid-(rp<<16); //wp<1048576
  
  rp>8?cyclicShift(xs,8):cyclicShift(xs,rp);
  rp>8?cyclicShift(xs,rp-8):cyclicShift(xs,0);  
  
  /*
  if(wp<=0)
  {
    return;
  }
  rp=wp>>16;      //0<=rp<16
  wp=wp-(rp<<16); //wp<65536
  if(rp>0)
  {
    bigMul(xs,W8[(rp-1)<<3]);
  }
  rp=wp>>12;
  wp=wp-(rp<<12);
  if(rp>0)
  {
    bigMul(xs,W12[(rp-1)<<3]);
  }
  rp=wp>>8;
  wp=wp-(rp<<8);
  if(rp>0)
  {
    bigMul(xs,W16[(rp-1)<<3]);
  }
  rp=wp>>4;
  wp=wp-(rp<<4);
  if(rp>0)
  {
    bigMul(xs,W20[(rp-1)<<3]);
  }
  if(wp>0)
  {
    bigMul(xs,W24[(rp-1)<<3]);
  }
  */
  if(wp>0)
  {
    bigMul(xs,&ws[wp<<3]);  
  }
}



__device__ void FFT16(unsigned long *xsm, unsigned long *ysm) 
{
  //Now, we only use one block to do $16$ FFT 16. 
  short tid = threadIdx.x; 
  short wid = (tid >> 5);  //warp no. [0,1,...,7]
  short bwid = (tid >> 6); //big warp no.[0,1,...,3]
  short sf = wid - (bwid << 1); //sub flag for odd warp[0,1]
  short wpn = ((tid-(wid<<5))>>3); //(tid-(wid*32))/8  [0,1,2,3] 
  short wpid = (tid-((tid>>3)<<3)); //in [0,1,...,7]
  short posid = 0; //in[0,1,...,15]
  short i, pos, pos1;
  unsigned long *xf, *yf;   
  int pos2 = ((bwid<<6)+(wpn<<4))<<3;
  
  xf = (unsigned long*)((char*)xsm + pos2*sizeof(unsigned long)); //
  yf = (unsigned long*)((char*)ysm + pos2*sizeof(unsigned long));
  
  //first round
  if(sf)
  {
    //bigSub
    bigSub(&xf[wpid*8],&xf[(wpid+8)*8],&yf[(wpid+8)*8]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[wpid*8],&xf[(wpid+8)*8],&yf[wpid*8]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  
  //second round
  posid=(wpid>=4?wpid+4:wpid);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+4)<<3],ind2[wpid]);
  }
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid<<3],&xf[(posid+4)<<3],&yf[(posid+4)<<3]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid<<3],&xf[(posid+4)<<3],&yf[posid<<3]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  
  //third round
  posid=(wpid>=4?((wpid-4)<<2)+1:(wpid<<2));
  if(sf>0)
  {
    cyclicShift(&xf[(posid+2)<<3],ind3[wpid]);
  }
  
  __syncthreads(); 
  if(sf>0)
  {
    //bigSub
    bigSub(&xf[posid<<3],&xf[(posid+2)<<3],&yf[(posid+2)<<3]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid<<3],&xf[(posid+2)<<3],&yf[posid<<3]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  
  //fourth
  posid=(wpid<<1);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+1)<<3],ind4[wpid]);
  }
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid<<3],&xf[(posid+1)<<3],&yf[(posid+1)<<3]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid<<3],&xf[(posid+1)<<3],&yf[posid<<3]);
  }
  __syncthreads();  
  
  
  if(sf)
  {
    posid=wpid+8;
    pos=posid<<3;
    pos1=ind5[posid]<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  else
  {
    posid=wpid;
    pos=posid<<3;
    pos1=ind5[posid]<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  __syncthreads(); 
  
}


/**
 * uniTwiddle is to twiddle a list of matrices xs, then store into ys.
 * xs: the pointer of $s$ matrices.
 * s: the matrix number
 * dnm: the data number of one matrix. dnm=he*wi
 * he: the height of matrix
 * wi: the width of matrix
 * I_s @ D_{he, wi}
 */
__global__ void uniTwiddle(unsigned long *xs, unsigned long *ws, short sb, short dnmb, short heb, short wib) 
{
    __shared__ unsigned long block[4096];  
    int id = (blockIdx.x<<9) + threadIdx.x;
    int tid = threadIdx.x;
    short dn = 1;  
    unsigned long *A;
    int idx, i,mno,pos,pos1,nh,nw,j,powind;
    N>=262144?(dn=(N>>18)):dn=1;
    
    idx=id;
    if(idx<N)
    {      
      for(i=0;i<dn;i++)
      {        
        mno=idx>>dnmb; //mno=(idx)/dnm; //表示第几个子矩阵 从0开始。
        pos1=idx-(mno<<dnmb); //表示该元素在子矩阵中的位置。
        nh=pos1>>wib; //旧的高度索引
        nw=pos1-(nh<<wib); //旧的宽度索引
        
        powind=nh*nw;
        pos=(tid<<3); //tid*8;
        A = (unsigned long *)((char*)xs + (idx<<3)*sizeof(unsigned long));
        for(j=0;j<8;j++)
        {
          block[pos++]=A[j];
        }
        if(powind > 0)
        {
          pos=tid<<3; //tid*8;
          mulOmega(&block[pos],powind<<sb,ws);
          for(j=0;j<8;j++)
          {
            A[j]=block[pos++];
          }
        }  
        idx+=262144; //i*262144; 
        if(idx>=N)
        {
          break;
        }             
      }    
    }
    __syncthreads();
}

/**
 * uniTranspose is to transpose a list of matrices xs, then store into ys.
 * xs: the pointer of $s$ matrices.
 * s: the matrix number
 * dnm: the data number of one matrix. dnm=he*wi
 * he: the height of matrix
 * wi: the width of matrix
 */
__global__ void uniTranspose(unsigned long *xs, unsigned long *ys, int s, short dnmb, short heb, short wib) 
{
    __shared__ unsigned long block[4096];  //512*8   4096*8=32KB 
    
    // The size of shared memory in a block is 48kb. 
    // We use 32KB, so the number of big numbers one block deal is 512.
    
    int id = (blockIdx.x<<9) + threadIdx.x;
    int tid = threadIdx.x;
    short dn = 1;  //the number each thread need deal.
    unsigned long *A, *B;
    int idx, i,mno,pos,pos1,h,w,pos2,j;
    N>=262144?(dn=(N>>18)):dn=1;
    
    if(id<N)
    {      
      idx=id;
      for(i=0;i<dn;i++)
      {
        if(idx>=N)
        {
          break;
        }  
        pos=tid<<3; //tid*8;              
        A = (unsigned long *)((char*)xs + (idx<<3)*sizeof(unsigned long));
        for(j=0;j<8;j++)
        {
          block[pos++]=A[j];
        }
        
        mno=idx>>dnmb; //mno=idx/dnm; //submatrix no, start from 0
        pos1=idx-(mno<<dnmb); //position in the current submatrix
        h=pos1>>wib; //old height
        w=pos1-(h<<wib); //new width
        pos2=(idx-pos1)+(w<<heb)+h;
        
        pos=tid<<3; //tid*8;
        B = (unsigned long *)((char*)ys + (pos2<<3)*sizeof(unsigned long));
        for(j=0;j<8;j++)
        {
          B[j]=block[pos++];
        }
        idx+=262144;  //idx=idx+i*262144;
      }    
    }
    __syncthreads();
}

__global__ void uniCopy(unsigned long *xs, unsigned long *ys) 
{
    int id = (blockIdx.x<<9) + threadIdx.x;
    short dn = 1;  //the number each thread need deal.
    int idx,i,j;
    unsigned long *A, *B;
    N>=262144?(dn=(N>>18)):dn=1;
    
    if(id<N)
    {  
      for(i=0;i<dn;i++)
      {
        idx=id+(i<<18); //i*262144;
        if(idx>=N)
        {
          break;
        }
        A = (unsigned long *)((char*)xs + idx*8*sizeof(unsigned long));
        B = (unsigned long *)((char*)ys + idx*8*sizeof(unsigned long));
        for(j=0;j<8;j++)
        {
          A[j]=B[j];
        }
      }
    }   
    
    __syncthreads();
}


void data_shuffle(int k, unsigned long *xs, unsigned long *ys)
{
  short i,j,dnmb,heb,wib;
  int s, he, wi, dnm;  //matrix number, height, width, data number for matrix
  wi=16;wib=4;
  s=1;  
  dnm=N;dnmb=k<<2;  //dnmb=24
  he=N/wi;heb=dnmb-4;
  //s*he*wi is equal to N
  //k=2  s=1  wi=16  he=16
  //k=3  s=1  wi=16  he=256
  //     s=16 wi=16  he=16
  
  j=0;
  for(i=0;i<=k-2;i++) //0,1,...,k-2
  {
    //xs transpose to ys
    //printf("call uniTranspose...\n");
    //printf("%d submatrix, each %d elements, %d * %d\n",s, dnm, he, wi);
    if(j==0)
    {
      j=1;
      uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(xs, ys, s, dnmb, heb, wib); 
      cudaThreadSynchronize();
    }
    else
    {
      j=0;
      uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(ys, xs, s, dnmb, heb, wib);  
      cudaThreadSynchronize();
    }    
    
    he=he>>4;   //he=he/16;
    heb=heb-4;
    s=s<<4;     //s=s*16;
    dnm=dnm>>4; //dnm=dnm/16;
    dnmb=dnmb-4;
  }
  if(j>0)
  {
    //copy ys to xs
    //printf("call uniCopy...\n");
    uniCopy<<<BLOCK_DIM,THREAD_DIM>>>(xs,ys);
    cudaThreadSynchronize();    
  }
  //cudaMemset((void **)&ys, 0, sizeof(unsigned long)*8*N);    
}

__global__ void data_fft16_kernel(unsigned long *xs) 
{
  short bid = blockIdx.x;
  short tid = threadIdx.x;
  int needBlock=(N>>8); //N/256;
  short bn = ((needBlock>=BLOCK_DIM)?needBlock/BLOCK_DIM:1);  //the number of block
  short i,j;
  //unsigned long *A, *B;
  unsigned long *A;
  __shared__ unsigned long xsm[INC<<3]; 
  __shared__ unsigned long ysm[INC<<3];  
  int pos,pos1,bno; //tid*8
  
  bno=bid;
  pos=((bno<<8)+tid)<<3;
  for(i=0;i<bn;i++)
  {
    if(bno>=needBlock)
    {
      break;
    }
    pos1 = tid<<3;
    //pos = (bno<<8)+tid;  //each block deal 256 big numbers.
    A = (unsigned long *)((char*)xs + pos*sizeof(unsigned long));
    
    for(j=0;j<8;j++)
    {
      xsm[pos1++]=A[j];
    }
    __syncthreads();  
    
    FFT16(xsm,ysm); //ysm is a temp array
    __syncthreads();  
    
    pos1 = (tid<<3);
    for(j=0;j<8;j++)
    {
      A[j]=xsm[pos1++];
    }
    bno+=512;
    pos+=1048576;  //512block*256bignumber*8numbers
  }
  __syncthreads();  
}

void data_fft(unsigned long *xs)
{
  //calculate xs in FFT16 batch, and store into ys
  data_fft16_kernel<<<BLOCK_DIM,256>>>(xs);
  cudaThreadSynchronize();
}


void data_twiddle(int k, unsigned long *xs, unsigned long *ys, unsigned long *ws)
{
  short i, sb, dnmb, heb, wib;
  int s;
  //s submatrix number
  for(i=k-2;i>=0;i--)
  {
    //T1  I_{16^i} @ D_{16,16^{k-1-i}}
    //s*dnm=s*he*wi
    s=(int) pow(16,i);
    sb=i<<2;
    //dnm=N/s;
    dnmb=(k-i)<<2; 
    heb=4; //he=16;
    wib=dnmb-heb;
    //printf("call uniTwiddle...\n");
    //printf("%d submatrix, each %d elements, %d * %d\n",s, N/s, he, N/s/he);
    uniTwiddle<<<BLOCK_DIM,THREAD_DIM>>>(xs, ws, sb, dnmb, heb, wib);
    cudaThreadSynchronize();
    
    //T2  I_{16^i} @ L_{16^{k-1-i}}^{16^{k-i}}
    dnmb=(k-i)<<2; 
    //wi=(int) pow(16,k-1-i);
    wib=(k-1-i)<<2;
    heb=4;
    //printf("call uniTranspose...\n");
    //printf("%d submatrix, each %d elements, %d * %d\n",s, dnm, he, wi);
    uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(xs, ys, s, dnmb, heb, wib); 
    cudaThreadSynchronize();
    
    
    //T3  I_{16^{k-1}} @ DTT16
    data_fft(ys);  //do fft on ys, store on ys. xs is a temp array
    
    //T4  I_{16^i} @ L_{16}^{16^{k-i}}
    //wi=16;
    wib=4;
    heb=dnmb-wib;
    //printf("call uniTranspose...\n");
    //printf("%d submatrix, each %d elements, %d * %d\n",s, dnm, he, wi);
    uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(ys, xs, s, dnmb, heb, wib); 
    cudaThreadSynchronize();
  }
}

void FFT1048576(unsigned long *xs, unsigned long *ws)
{
  int k=5;
  unsigned long *xs_d, *ys_d, *ws_d;
  cudaEvent_t start, stop;
  float elapsedTime;
  
  //initiation
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*N);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*N, cudaMemcpyHostToDevice);   cudaMalloc((void **)&ys_d, sizeof(unsigned long)*8*N);
  
  cudaMemset((void **)&ys_d, 0, sizeof(unsigned long)*8*N);
  
  cudaMalloc((void **)&ws_d, sizeof(unsigned long)*8*N/16);
  cudaMemcpy(ws_d, ws, sizeof(unsigned long)*8*N/16, cudaMemcpyHostToDevice);
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  //P1
  data_shuffle(k, xs_d, ys_d); 
  //printf("P1 end...\n");
  
  //P2
  data_fft(xs_d);
  //printf("P2 end...\n");
  
  
  //P3
  data_twiddle(k, xs_d, ys_d, ws_d);
  //printf("P3 end...\n");
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  cudaMemcpy(xs, xs_d, sizeof(unsigned long)*8*N, cudaMemcpyDeviceToHost);
  
  printf("we have done xs = FFT1048576(xs) in 1 times.\n");
  printf("the time of gpu is %f ms\n", elapsedTime);
  
  cudaFree(xs_d);	
  cudaFree(ys_d);	
  cudaFree(ws_d);	
}

int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1, *fp2, *fp3;	 
  unsigned long *xs,*ws; 
  
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*N));
  ws=(unsigned long *)malloc((sizeof(unsigned long)*8*(N/16)));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_1048576_input.dat"); 
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  memset(xs,0,sizeof(unsigned long)*8*N);
  fread(xs,sizeof(unsigned long),8*N,fp1);
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "W_1048576.dat"); 
  if((fp3=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp3);
  	exit(-1);
  } 
  memset(ws,0,sizeof(unsigned long)*8*(N/16));
  fread(ws,sizeof(unsigned long),8*N/16,fp3);
  
  
  FFT1048576(xs,ws);
    
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_1048576_output.dat"); //INC*8 (unsigned long) data
  if((fp2=fopen(fileName,"wb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp2);
  	exit(-1);
  } 
  fwrite(xs,sizeof(unsigned long),8*N,fp2); 
  
  //free
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  
  free(xs);
  return 0;
}



