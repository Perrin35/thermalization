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
rz(-1.9300652) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601215) q[0];
sx q[0];
rz(-1.6257964) q[0];
sx q[0];
rz(-2.4762857) q[0];
rz(-pi) q[1];
rz(1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(0.36117902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71516192) q[1];
sx q[1];
rz(-0.75725812) q[1];
sx q[1];
rz(0.84233474) q[1];
x q[2];
rz(-1.6739507) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(0.32675693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4988929) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(0.47168628) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1175849) q[0];
sx q[0];
rz(-1.4837259) q[0];
sx q[0];
rz(2.9136806) q[0];
x q[1];
rz(1.291044) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(2.5807057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43863152) q[1];
sx q[1];
rz(-1.8384117) q[1];
sx q[1];
rz(2.5665934) q[1];
x q[2];
rz(-0.15700335) q[3];
sx q[3];
rz(-0.96884851) q[3];
sx q[3];
rz(-2.4412145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(-2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(-2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(-1.4556494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3464976) q[0];
sx q[0];
rz(-0.99612757) q[0];
sx q[0];
rz(2.5618308) q[0];
rz(-0.74080148) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(1.4893116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.599217) q[1];
sx q[1];
rz(-1.5812751) q[1];
sx q[1];
rz(3.0159365) q[1];
x q[2];
rz(0.62526838) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-2.1739615) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2847292) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(-1.2329873) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1283135) q[0];
sx q[0];
rz(-1.0767125) q[0];
sx q[0];
rz(-1.6492776) q[0];
rz(-pi) q[1];
rz(1.0225251) q[2];
sx q[2];
rz(-1.0523445) q[2];
sx q[2];
rz(0.099629121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14703688) q[1];
sx q[1];
rz(-1.5347267) q[1];
sx q[1];
rz(1.7493164) q[1];
x q[2];
rz(-0.87807699) q[3];
sx q[3];
rz(-1.3948166) q[3];
sx q[3];
rz(2.7443977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(0.81400648) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(2.9290501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4310303) q[0];
sx q[0];
rz(-3.115603) q[0];
sx q[0];
rz(0.89528577) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1754155) q[2];
sx q[2];
rz(-1.539131) q[2];
sx q[2];
rz(0.67827536) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7674539) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(-2.2351082) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6796474) q[3];
sx q[3];
rz(-2.2432703) q[3];
sx q[3];
rz(0.40601054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(-2.6203716) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(2.916472) q[0];
rz(-1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-2.7640142) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29992732) q[0];
sx q[0];
rz(-2.8746434) q[0];
sx q[0];
rz(1.4873234) q[0];
rz(-pi) q[1];
rz(1.5224783) q[2];
sx q[2];
rz(-1.8579351) q[2];
sx q[2];
rz(1.552812) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67205091) q[1];
sx q[1];
rz(-0.44476032) q[1];
sx q[1];
rz(1.3252392) q[1];
rz(-pi) q[2];
rz(-0.79490957) q[3];
sx q[3];
rz(-1.5494293) q[3];
sx q[3];
rz(1.6705318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0525381) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(0.19454923) q[0];
rz(2.9220707) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-2.887168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5146778) q[0];
sx q[0];
rz(-1.2509545) q[0];
sx q[0];
rz(-0.11492782) q[0];
rz(-pi) q[1];
rz(1.3527855) q[2];
sx q[2];
rz(-0.62354747) q[2];
sx q[2];
rz(-2.4965198) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19941482) q[1];
sx q[1];
rz(-2.605038) q[1];
sx q[1];
rz(-2.5903827) q[1];
rz(-pi) q[2];
rz(-0.80612225) q[3];
sx q[3];
rz(-2.5297909) q[3];
sx q[3];
rz(2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.8257726) q[2];
rz(-2.2655462) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(-2.3176106) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5751858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5113735) q[0];
sx q[0];
rz(-2.2492118) q[0];
sx q[0];
rz(1.9766962) q[0];
x q[1];
rz(1.0303866) q[2];
sx q[2];
rz(-1.724913) q[2];
sx q[2];
rz(0.87481462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7745061) q[1];
sx q[1];
rz(-1.204406) q[1];
sx q[1];
rz(-2.7584502) q[1];
rz(-pi) q[2];
rz(-2.6353587) q[3];
sx q[3];
rz(-1.1971548) q[3];
sx q[3];
rz(0.95844275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(0.27754647) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9797416) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(2.004752) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4515848) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(2.8103229) q[0];
x q[1];
rz(-2.7669737) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(1.7916726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1946823) q[1];
sx q[1];
rz(-1.8872675) q[1];
sx q[1];
rz(3.0419993) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7342019) q[3];
sx q[3];
rz(-0.77851495) q[3];
sx q[3];
rz(0.22637573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64951605) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(-1.0409522) q[2];
rz(-3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025538) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(2.1955406) q[0];
rz(-0.91167766) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(2.5295703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48869041) q[0];
sx q[0];
rz(-2.2757747) q[0];
sx q[0];
rz(-2.5253354) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34531784) q[2];
sx q[2];
rz(-1.3441663) q[2];
sx q[2];
rz(-0.40030865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.260173) q[1];
sx q[1];
rz(-2.2248785) q[1];
sx q[1];
rz(-1.7263078) q[1];
x q[2];
rz(-2.5522334) q[3];
sx q[3];
rz(-0.61925626) q[3];
sx q[3];
rz(-2.9818231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0845906) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(-0.65336147) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(-2.4408834) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54031298) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(-2.3868949) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(0.66773141) q[2];
sx q[2];
rz(-0.79955352) q[2];
sx q[2];
rz(-0.59895589) q[2];
rz(1.894542) q[3];
sx q[3];
rz(-1.6508045) q[3];
sx q[3];
rz(1.9423021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
