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
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(0.29830631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0878108) q[0];
sx q[0];
rz(-0.90667533) q[0];
sx q[0];
rz(1.640663) q[0];
x q[1];
rz(-0.093703336) q[2];
sx q[2];
rz(-1.2393349) q[2];
sx q[2];
rz(1.1790438) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5391985) q[1];
sx q[1];
rz(-1.0326003) q[1];
sx q[1];
rz(2.5799275) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88793036) q[3];
sx q[3];
rz(-3.008932) q[3];
sx q[3];
rz(-0.56548972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.0304334) q[1];
sx q[1];
rz(-0.47168628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2361006) q[0];
sx q[0];
rz(-2.8978851) q[0];
sx q[0];
rz(2.7729176) q[0];
rz(-pi) q[1];
rz(0.47927803) q[2];
sx q[2];
rz(-2.5844378) q[2];
sx q[2];
rz(-3.1301168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8397123) q[1];
sx q[1];
rz(-1.0186968) q[1];
sx q[1];
rz(-1.2549972) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.178533) q[3];
sx q[3];
rz(-1.4415603) q[3];
sx q[3];
rz(0.78101633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(2.1726051) q[2];
rz(-0.5747059) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-0.18181268) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.4556494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91711125) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(-0.86865058) q[0];
rz(-pi) q[1];
rz(-0.45708926) q[2];
sx q[2];
rz(-0.79490137) q[2];
sx q[2];
rz(2.7283816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94565832) q[1];
sx q[1];
rz(-3.0155026) q[1];
sx q[1];
rz(0.083421589) q[1];
x q[2];
rz(-1.6772179) q[3];
sx q[3];
rz(-2.1933746) q[3];
sx q[3];
rz(2.4604083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9905332) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(2.997056) q[0];
x q[1];
rz(1.0225251) q[2];
sx q[2];
rz(-2.0892482) q[2];
sx q[2];
rz(3.0419635) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7113263) q[1];
sx q[1];
rz(-1.749199) q[1];
sx q[1];
rz(0.036651595) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9144985) q[3];
sx q[3];
rz(-2.2507651) q[3];
sx q[3];
rz(1.029315) q[3];
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
rz(-1.043184) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
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
rz(-2.2463069) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5151305) q[2];
sx q[2];
rz(-0.60534436) q[2];
sx q[2];
rz(-2.2948613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5383496) q[1];
sx q[1];
rz(-1.0298567) q[1];
sx q[1];
rz(2.4451838) q[1];
rz(2.0810633) q[3];
sx q[3];
rz(-0.79499309) q[3];
sx q[3];
rz(-2.0612962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(2.916472) q[0];
rz(-1.7865932) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(-0.37757847) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9281884) q[0];
sx q[0];
rz(-1.8367935) q[0];
sx q[0];
rz(-0.022797419) q[0];
rz(-2.9794681) q[2];
sx q[2];
rz(-0.29106489) q[2];
sx q[2];
rz(1.7578917) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4695417) q[1];
sx q[1];
rz(-0.44476032) q[1];
sx q[1];
rz(-1.8163535) q[1];
rz(-pi) q[2];
rz(3.1116629) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(-3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(-2.6605576) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(0.25442466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5146778) q[0];
sx q[0];
rz(-1.2509545) q[0];
sx q[0];
rz(-0.11492782) q[0];
rz(-pi) q[1];
rz(-1.7888072) q[2];
sx q[2];
rz(-2.5180452) q[2];
sx q[2];
rz(2.4965198) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88528819) q[1];
sx q[1];
rz(-1.2997775) q[1];
sx q[1];
rz(0.46896743) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1021348) q[3];
sx q[3];
rz(-1.161876) q[3];
sx q[3];
rz(1.8917781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97757942) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(-2.2655462) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5751858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4644509) q[0];
sx q[0];
rz(-1.8832708) q[0];
sx q[0];
rz(-0.72014767) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17922108) q[2];
sx q[2];
rz(-2.1041098) q[2];
sx q[2];
rz(-0.78782493) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6187001) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(-0.79843847) q[1];
x q[2];
rz(-2.4616562) q[3];
sx q[3];
rz(-2.5222062) q[3];
sx q[3];
rz(-3.1114651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.4510441) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16185109) q[0];
sx q[0];
rz(-0.45270544) q[0];
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
rz(0.046969983) q[0];
sx q[0];
rz(-1.2828865) q[0];
sx q[0];
rz(2.1043491) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9929664) q[2];
sx q[2];
rz(-1.2264226) q[2];
sx q[2];
rz(-3.0712155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7965664) q[1];
sx q[1];
rz(-1.6654286) q[1];
sx q[1];
rz(1.8887397) q[1];
x q[2];
rz(0.40739079) q[3];
sx q[3];
rz(-2.3630777) q[3];
sx q[3];
rz(-0.22637573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(-0.0020290931) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8390389) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(0.94605207) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(2.5295703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8004868) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(-0.97408803) q[0];
rz(-pi) q[1];
rz(0.59801306) q[2];
sx q[2];
rz(-2.7310555) q[2];
sx q[2];
rz(1.41278) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0079841) q[1];
sx q[1];
rz(-0.66966479) q[1];
sx q[1];
rz(2.942251) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58935921) q[3];
sx q[3];
rz(-0.61925626) q[3];
sx q[3];
rz(0.15976957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0570021) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.54031298) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(2.3868949) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-1.0036219) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(-1.8176953) q[3];
sx q[3];
rz(-0.3331475) q[3];
sx q[3];
rz(-3.0039136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];