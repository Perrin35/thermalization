OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(-1.3949431) q[0];
sx q[0];
rz(-3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66520663) q[0];
sx q[0];
rz(-1.6205661) q[0];
sx q[0];
rz(0.14270356) q[0];
rz(-1.0983659) q[2];
sx q[2];
rz(-0.78394267) q[2];
sx q[2];
rz(-0.69477615) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2831813) q[1];
sx q[1];
rz(-1.6106669) q[1];
sx q[1];
rz(1.7130501) q[1];
x q[2];
rz(-2.6585456) q[3];
sx q[3];
rz(-0.31937283) q[3];
sx q[3];
rz(-1.3787624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.8700245) q[0];
rz(-1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028582024) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(0.21582614) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.045250821) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(2.0766052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83982044) q[1];
sx q[1];
rz(-0.41026792) q[1];
sx q[1];
rz(-0.46373414) q[1];
rz(2.1167011) q[3];
sx q[3];
rz(-0.62666577) q[3];
sx q[3];
rz(2.6729667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21330825) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(-1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-0.24060732) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89861682) q[0];
sx q[0];
rz(-1.9447864) q[0];
sx q[0];
rz(0.6181194) q[0];
rz(-pi) q[1];
rz(0.37249506) q[2];
sx q[2];
rz(-1.8066346) q[2];
sx q[2];
rz(-1.3553938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5879844) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(1.3081461) q[1];
x q[2];
rz(2.5127605) q[3];
sx q[3];
rz(-1.9753846) q[3];
sx q[3];
rz(1.42266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9712517) q[0];
sx q[0];
rz(-0.60702885) q[0];
sx q[0];
rz(1.6334565) q[0];
x q[1];
rz(-0.36074952) q[2];
sx q[2];
rz(-1.4348239) q[2];
sx q[2];
rz(-0.13775682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3164697) q[1];
sx q[1];
rz(-2.9087062) q[1];
sx q[1];
rz(-1.4941925) q[1];
rz(2.5424764) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(-1.8653387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63311657) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.9015076) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2729028) q[0];
sx q[0];
rz(-1.5575952) q[0];
sx q[0];
rz(0.23916434) q[0];
x q[1];
rz(2.901022) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(-2.0672928) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48092914) q[1];
sx q[1];
rz(-0.96748057) q[1];
sx q[1];
rz(1.4900581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6579882) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(-1.4809158) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4303362) q[0];
sx q[0];
rz(-0.38861409) q[0];
sx q[0];
rz(-1.3769763) q[0];
rz(-1.5724206) q[2];
sx q[2];
rz(-2.6964028) q[2];
sx q[2];
rz(2.0507398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33038352) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(-1.1793544) q[1];
x q[2];
rz(-2.2495066) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(2.2298262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(-2.5793502) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.4690171) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.3247103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0990471) q[0];
sx q[0];
rz(-1.5934266) q[0];
sx q[0];
rz(-0.12734539) q[0];
rz(-pi) q[1];
rz(-1.0075188) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(-1.2457459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99609375) q[1];
sx q[1];
rz(-2.5719574) q[1];
sx q[1];
rz(2.3949404) q[1];
rz(-pi) q[2];
rz(0.91388254) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(-0.1012181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(2.6810714) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.2896279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562718) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(-0.69099364) q[0];
rz(-pi) q[1];
rz(2.3194359) q[2];
sx q[2];
rz(-0.90688721) q[2];
sx q[2];
rz(-0.31809959) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3813078) q[1];
sx q[1];
rz(-2.9236301) q[1];
sx q[1];
rz(-1.0337619) q[1];
rz(-0.21059489) q[3];
sx q[3];
rz(-2.4678293) q[3];
sx q[3];
rz(-1.8316366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.148968) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(-3.0454214) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3960421) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(1.0066443) q[0];
rz(-pi) q[1];
rz(2.8106868) q[2];
sx q[2];
rz(-1.3661642) q[2];
sx q[2];
rz(-0.76588878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.64687318) q[1];
sx q[1];
rz(-0.44534007) q[1];
sx q[1];
rz(-2.7237707) q[1];
rz(-pi) q[2];
rz(0.11741365) q[3];
sx q[3];
rz(-0.88866361) q[3];
sx q[3];
rz(-2.5168602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081351) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(-1.0159147) q[0];
rz(-2.290906) q[2];
sx q[2];
rz(-0.88973532) q[2];
sx q[2];
rz(-2.1249287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0444813) q[1];
sx q[1];
rz(-2.0437355) q[1];
sx q[1];
rz(-0.60826917) q[1];
x q[2];
rz(0.82716771) q[3];
sx q[3];
rz(-2.4359772) q[3];
sx q[3];
rz(0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(2.091999) q[2];
sx q[2];
rz(-2.1672719) q[2];
sx q[2];
rz(1.5530829) q[2];
rz(2.6736005) q[3];
sx q[3];
rz(-2.7157126) q[3];
sx q[3];
rz(-1.4117905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];