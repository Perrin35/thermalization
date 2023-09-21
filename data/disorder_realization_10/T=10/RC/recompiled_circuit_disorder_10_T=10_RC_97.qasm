OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(1.2115275) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(-2.8432863) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5601215) q[0];
sx q[0];
rz(-1.6257964) q[0];
sx q[0];
rz(2.4762857) q[0];
rz(-pi) q[1];
rz(-0.093703336) q[2];
sx q[2];
rz(-1.9022577) q[2];
sx q[2];
rz(1.9625488) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71516192) q[1];
sx q[1];
rz(-2.3843345) q[1];
sx q[1];
rz(0.84233474) q[1];
rz(-0.084007752) q[3];
sx q[3];
rz(-1.6735895) q[3];
sx q[3];
rz(-1.8889129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(-0.99386627) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4988929) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(2.6699064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1175849) q[0];
sx q[0];
rz(-1.4837259) q[0];
sx q[0];
rz(-2.9136806) q[0];
rz(-pi) q[1];
rz(-2.6623146) q[2];
sx q[2];
rz(-0.55715484) q[2];
sx q[2];
rz(3.1301168) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43863152) q[1];
sx q[1];
rz(-1.8384117) q[1];
sx q[1];
rz(2.5665934) q[1];
rz(-pi) q[2];
rz(-0.96305965) q[3];
sx q[3];
rz(-1.7000323) q[3];
sx q[3];
rz(-2.3605763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(-0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(0.18181268) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.4556494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244814) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(0.86865058) q[0];
rz(-2.4007912) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(-1.4893116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1959343) q[1];
sx q[1];
rz(-0.12609005) q[1];
sx q[1];
rz(-3.0581711) q[1];
rz(-pi) q[2];
rz(1.4643747) q[3];
sx q[3];
rz(-2.1933746) q[3];
sx q[3];
rz(2.4604083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(-0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(-2.0137285) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5467984) q[0];
sx q[0];
rz(-1.6398755) q[0];
sx q[0];
rz(-2.6462206) q[0];
rz(2.1190676) q[2];
sx q[2];
rz(-2.0892482) q[2];
sx q[2];
rz(-3.0419635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14703688) q[1];
sx q[1];
rz(-1.5347267) q[1];
sx q[1];
rz(-1.3922763) q[1];
x q[2];
rz(-0.22709417) q[3];
sx q[3];
rz(-2.2507651) q[3];
sx q[3];
rz(-1.029315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91810742) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(-0.81400648) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-2.2498851) q[0];
rz(1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(2.9290501) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5355797) q[0];
sx q[0];
rz(-1.5545462) q[0];
sx q[0];
rz(-1.5910801) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.038482484) q[2];
sx q[2];
rz(-2.1750692) q[2];
sx q[2];
rz(-2.2272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5383496) q[1];
sx q[1];
rz(-2.1117359) q[1];
sx q[1];
rz(-2.4451838) q[1];
rz(2.0810633) q[3];
sx q[3];
rz(-2.3465996) q[3];
sx q[3];
rz(2.0612962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(-0.5212211) q[2];
rz(-1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.3549995) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(2.7640142) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21340428) q[0];
sx q[0];
rz(-1.8367935) q[0];
sx q[0];
rz(-0.022797419) q[0];
x q[1];
rz(0.16212459) q[2];
sx q[2];
rz(-2.8505278) q[2];
sx q[2];
rz(-1.7578917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40124711) q[1];
sx q[1];
rz(-1.1402854) q[1];
sx q[1];
rz(0.11534782) q[1];
x q[2];
rz(-3.1116629) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-2.887168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1628111) q[0];
sx q[0];
rz(-2.8023976) q[0];
sx q[0];
rz(-1.9041054) q[0];
x q[1];
rz(-1.7888072) q[2];
sx q[2];
rz(-2.5180452) q[2];
sx q[2];
rz(-0.64507285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9421778) q[1];
sx q[1];
rz(-2.605038) q[1];
sx q[1];
rz(0.55120991) q[1];
x q[2];
rz(1.1021348) q[3];
sx q[3];
rz(-1.9797167) q[3];
sx q[3];
rz(1.2498145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(-1.8257726) q[2];
rz(-0.87604648) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106237) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(-2.3176106) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(1.5664068) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9109089) q[0];
sx q[0];
rz(-2.3678603) q[0];
sx q[0];
rz(0.45549972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17922108) q[2];
sx q[2];
rz(-2.1041098) q[2];
sx q[2];
rz(0.78782493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3670866) q[1];
sx q[1];
rz(-1.204406) q[1];
sx q[1];
rz(-0.3831425) q[1];
rz(-1.992222) q[3];
sx q[3];
rz(-2.0391658) q[3];
sx q[3];
rz(-2.3295662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.4510441) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(-1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(1.45654) q[0];
rz(2.5121571) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(-2.004752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0658543) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(-1.0435186) q[0];
rz(-pi) q[1];
rz(-1.1486263) q[2];
sx q[2];
rz(-1.9151701) q[2];
sx q[2];
rz(-3.0712155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5053619) q[1];
sx q[1];
rz(-2.8103235) q[1];
sx q[1];
rz(-1.2760217) q[1];
rz(2.7342019) q[3];
sx q[3];
rz(-0.77851495) q[3];
sx q[3];
rz(2.9152169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64951605) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(-1.0409522) q[2];
rz(-0.0020290931) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(-0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8390389) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(2.1955406) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-2.5295703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5123972) q[0];
sx q[0];
rz(-1.1149659) q[0];
sx q[0];
rz(0.7645316) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5435796) q[2];
sx q[2];
rz(-2.7310555) q[2];
sx q[2];
rz(-1.7288127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1336085) q[1];
sx q[1];
rz(-0.66966479) q[1];
sx q[1];
rz(-0.19934166) q[1];
rz(0.53491433) q[3];
sx q[3];
rz(-1.8992918) q[3];
sx q[3];
rz(1.9096149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(-2.7907794) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(-0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54031298) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-1.0036219) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(-1.2470506) q[3];
sx q[3];
rz(-1.6508045) q[3];
sx q[3];
rz(1.9423021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];