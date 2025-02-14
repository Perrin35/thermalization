OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91035834) q[0];
sx q[0];
rz(-2.2725821) q[0];
sx q[0];
rz(2.056871) q[0];
rz(-1.1551069) q[1];
sx q[1];
rz(-0.81973633) q[1];
sx q[1];
rz(-0.91135946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40075743) q[0];
sx q[0];
rz(-2.3199953) q[0];
sx q[0];
rz(-1.9039959) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16884825) q[2];
sx q[2];
rz(-0.71879866) q[2];
sx q[2];
rz(1.079725) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5973917) q[1];
sx q[1];
rz(-0.72300882) q[1];
sx q[1];
rz(0.71031481) q[1];
x q[2];
rz(-0.064685589) q[3];
sx q[3];
rz(-1.1380592) q[3];
sx q[3];
rz(2.329934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34065166) q[2];
sx q[2];
rz(-1.1110577) q[2];
sx q[2];
rz(-3.0621373) q[2];
rz(2.5331412) q[3];
sx q[3];
rz(-1.2234917) q[3];
sx q[3];
rz(2.2505545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65421739) q[0];
sx q[0];
rz(-2.2779164) q[0];
sx q[0];
rz(1.7864216) q[0];
rz(-1.0379418) q[1];
sx q[1];
rz(-2.4619921) q[1];
sx q[1];
rz(2.3341446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8184745) q[0];
sx q[0];
rz(-1.2809296) q[0];
sx q[0];
rz(-0.94375837) q[0];
rz(-pi) q[1];
rz(-0.37952559) q[2];
sx q[2];
rz(-2.1132601) q[2];
sx q[2];
rz(-1.8636974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3088919) q[1];
sx q[1];
rz(-1.3314134) q[1];
sx q[1];
rz(1.5025839) q[1];
rz(-pi) q[2];
rz(1.7044675) q[3];
sx q[3];
rz(-1.6764056) q[3];
sx q[3];
rz(1.9131806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21343931) q[2];
sx q[2];
rz(-1.5779053) q[2];
sx q[2];
rz(-1.4595002) q[2];
rz(-2.6835119) q[3];
sx q[3];
rz(-2.5638678) q[3];
sx q[3];
rz(-0.95988449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15993519) q[0];
sx q[0];
rz(-2.0913048) q[0];
sx q[0];
rz(-0.67498573) q[0];
rz(-2.5534897) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(1.810422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0494306) q[0];
sx q[0];
rz(-1.9665008) q[0];
sx q[0];
rz(-0.53073287) q[0];
rz(-pi) q[1];
rz(0.66548062) q[2];
sx q[2];
rz(-1.988171) q[2];
sx q[2];
rz(0.62705428) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71215502) q[1];
sx q[1];
rz(-2.3596016) q[1];
sx q[1];
rz(-1.9351134) q[1];
x q[2];
rz(2.0698333) q[3];
sx q[3];
rz(-1.3924358) q[3];
sx q[3];
rz(2.2638829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62446928) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(-0.2207174) q[2];
rz(0.80495009) q[3];
sx q[3];
rz(-1.5419818) q[3];
sx q[3];
rz(-1.3360924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.540156) q[0];
sx q[0];
rz(-2.8995081) q[0];
sx q[0];
rz(2.5228187) q[0];
rz(-0.082911804) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(1.5203016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7474866) q[0];
sx q[0];
rz(-1.0472676) q[0];
sx q[0];
rz(3.0305358) q[0];
rz(-pi) q[1];
rz(-0.57915202) q[2];
sx q[2];
rz(-1.7699827) q[2];
sx q[2];
rz(0.79328377) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1136696) q[1];
sx q[1];
rz(-1.6531367) q[1];
sx q[1];
rz(1.0587949) q[1];
rz(0.374745) q[3];
sx q[3];
rz(-1.4029158) q[3];
sx q[3];
rz(1.8207267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3380022) q[2];
sx q[2];
rz(-0.10927304) q[2];
sx q[2];
rz(2.9962311) q[2];
rz(-0.27211443) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(-1.9947778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014932545) q[0];
sx q[0];
rz(-1.8155875) q[0];
sx q[0];
rz(3.0354011) q[0];
rz(0.95942489) q[1];
sx q[1];
rz(-2.5966849) q[1];
sx q[1];
rz(-2.9471961) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0036164408) q[0];
sx q[0];
rz(-0.73775333) q[0];
sx q[0];
rz(-0.12807782) q[0];
rz(-pi) q[1];
rz(-1.4559559) q[2];
sx q[2];
rz(-1.5659589) q[2];
sx q[2];
rz(-0.23813914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.55479706) q[1];
sx q[1];
rz(-1.5298843) q[1];
sx q[1];
rz(0.58920963) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.077094519) q[3];
sx q[3];
rz(-1.9760248) q[3];
sx q[3];
rz(-0.73540686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.306119) q[2];
sx q[2];
rz(-1.6872311) q[2];
sx q[2];
rz(-0.56841889) q[2];
rz(-0.052637188) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(-3.0047825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0282054) q[0];
sx q[0];
rz(-0.55332342) q[0];
sx q[0];
rz(-2.9521039) q[0];
rz(-1.0264617) q[1];
sx q[1];
rz(-2.2920513) q[1];
sx q[1];
rz(-0.75526563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3657235) q[0];
sx q[0];
rz(-2.2138174) q[0];
sx q[0];
rz(0.68277208) q[0];
x q[1];
rz(-0.83024518) q[2];
sx q[2];
rz(-2.4355222) q[2];
sx q[2];
rz(0.65078562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2456296) q[1];
sx q[1];
rz(-1.4019483) q[1];
sx q[1];
rz(-2.9814475) q[1];
x q[2];
rz(2.7868458) q[3];
sx q[3];
rz(-1.915853) q[3];
sx q[3];
rz(-1.7999032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3834164) q[2];
sx q[2];
rz(-1.5864317) q[2];
sx q[2];
rz(-1.4413393) q[2];
rz(-1.3759184) q[3];
sx q[3];
rz(-0.91840363) q[3];
sx q[3];
rz(-2.4413696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43948424) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(-1.9973607) q[0];
rz(-2.8817835) q[1];
sx q[1];
rz(-0.83994284) q[1];
sx q[1];
rz(2.1536486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081962498) q[0];
sx q[0];
rz(-3.0407627) q[0];
sx q[0];
rz(-0.16720812) q[0];
x q[1];
rz(1.1455215) q[2];
sx q[2];
rz(-0.88103154) q[2];
sx q[2];
rz(0.20343101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14079796) q[1];
sx q[1];
rz(-0.70604815) q[1];
sx q[1];
rz(2.7435859) q[1];
rz(-pi) q[2];
rz(2.0201398) q[3];
sx q[3];
rz(-2.5183023) q[3];
sx q[3];
rz(1.9188855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40954956) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(0.6130971) q[2];
rz(0.71550718) q[3];
sx q[3];
rz(-2.3698273) q[3];
sx q[3];
rz(-2.7806921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9382984) q[0];
sx q[0];
rz(-2.2977915) q[0];
sx q[0];
rz(2.8734558) q[0];
rz(-0.2306436) q[1];
sx q[1];
rz(-1.5985039) q[1];
sx q[1];
rz(3.1275829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62416158) q[0];
sx q[0];
rz(-0.78963477) q[0];
sx q[0];
rz(2.8454418) q[0];
x q[1];
rz(-1.2410937) q[2];
sx q[2];
rz(-0.46247855) q[2];
sx q[2];
rz(-1.3792737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9585628) q[1];
sx q[1];
rz(-1.5261478) q[1];
sx q[1];
rz(-2.6041998) q[1];
rz(-pi) q[2];
rz(-1.713578) q[3];
sx q[3];
rz(-1.5413949) q[3];
sx q[3];
rz(-2.2063696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1676499) q[2];
sx q[2];
rz(-1.3385945) q[2];
sx q[2];
rz(2.8323284) q[2];
rz(-0.15657982) q[3];
sx q[3];
rz(-1.4045818) q[3];
sx q[3];
rz(-2.1046765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2579047) q[0];
sx q[0];
rz(-2.4990999) q[0];
sx q[0];
rz(2.8908492) q[0];
rz(-2.5550487) q[1];
sx q[1];
rz(-1.5778912) q[1];
sx q[1];
rz(0.7647382) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79478489) q[0];
sx q[0];
rz(-1.2803923) q[0];
sx q[0];
rz(1.8998763) q[0];
rz(-pi) q[1];
x q[1];
rz(2.370442) q[2];
sx q[2];
rz(-1.4284819) q[2];
sx q[2];
rz(2.3754295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.918461) q[1];
sx q[1];
rz(-1.4479965) q[1];
sx q[1];
rz(1.5614913) q[1];
x q[2];
rz(0.37882355) q[3];
sx q[3];
rz(-1.6215542) q[3];
sx q[3];
rz(-1.7876884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5137198) q[2];
sx q[2];
rz(-2.3508115) q[2];
sx q[2];
rz(-0.17390832) q[2];
rz(0.74226132) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(-0.043005634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0088418) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(-0.26790628) q[0];
rz(-1.8065037) q[1];
sx q[1];
rz(-0.91064149) q[1];
sx q[1];
rz(-1.7693899) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5822179) q[0];
sx q[0];
rz(-1.1211044) q[0];
sx q[0];
rz(2.873308) q[0];
rz(0.32023264) q[2];
sx q[2];
rz(-2.8304184) q[2];
sx q[2];
rz(2.8428889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6698996) q[1];
sx q[1];
rz(-2.0834384) q[1];
sx q[1];
rz(2.1218845) q[1];
rz(-2.3389111) q[3];
sx q[3];
rz(-2.6978328) q[3];
sx q[3];
rz(-2.8783523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.043896) q[2];
sx q[2];
rz(-1.0363657) q[2];
sx q[2];
rz(-1.7608661) q[2];
rz(-1.1287639) q[3];
sx q[3];
rz(-2.1916316) q[3];
sx q[3];
rz(-1.5855764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783405) q[0];
sx q[0];
rz(-1.5338407) q[0];
sx q[0];
rz(-1.1080678) q[0];
rz(2.4551328) q[1];
sx q[1];
rz(-1.3931128) q[1];
sx q[1];
rz(-1.211094) q[1];
rz(0.5067208) q[2];
sx q[2];
rz(-0.55247775) q[2];
sx q[2];
rz(-1.1341118) q[2];
rz(1.0330647) q[3];
sx q[3];
rz(-1.929639) q[3];
sx q[3];
rz(0.7076984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
