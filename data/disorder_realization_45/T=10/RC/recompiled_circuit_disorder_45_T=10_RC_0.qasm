OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(-1.1614769) q[1];
sx q[1];
rz(-2.5002313) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25783378) q[0];
sx q[0];
rz(-2.6834282) q[0];
sx q[0];
rz(-0.56295653) q[0];
x q[1];
rz(-1.1429206) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(-1.2905754) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8909059) q[1];
sx q[1];
rz(-0.85039925) q[1];
sx q[1];
rz(-2.5198063) q[1];
x q[2];
rz(-2.6061329) q[3];
sx q[3];
rz(-0.70264953) q[3];
sx q[3];
rz(3.1010166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(0.47810289) q[2];
rz(1.452662) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801382) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(-1.0189198) q[0];
rz(1.4787176) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(2.5085124) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2791726) q[0];
sx q[0];
rz(-2.3775568) q[0];
sx q[0];
rz(0.11636244) q[0];
rz(-pi) q[1];
rz(-0.74866809) q[2];
sx q[2];
rz(-2.0249172) q[2];
sx q[2];
rz(1.5926966) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1770597) q[1];
sx q[1];
rz(-2.6033083) q[1];
sx q[1];
rz(-2.2798377) q[1];
x q[2];
rz(-0.88926104) q[3];
sx q[3];
rz(-1.8157496) q[3];
sx q[3];
rz(-0.8196866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(-3.0241372) q[2];
rz(2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(0.69082469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3120964) q[0];
sx q[0];
rz(-1.7236241) q[0];
sx q[0];
rz(1.494708) q[0];
rz(-2.3766741) q[2];
sx q[2];
rz(-1.3996482) q[2];
sx q[2];
rz(0.82676065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1289127) q[1];
sx q[1];
rz(-1.3887654) q[1];
sx q[1];
rz(-3.1151031) q[1];
rz(-pi) q[2];
rz(1.2191804) q[3];
sx q[3];
rz(-1.2120486) q[3];
sx q[3];
rz(-0.84695942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92418015) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(-0.10647354) q[2];
rz(1.4364093) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(-0.52247125) q[0];
rz(-2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(-2.5879588) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.583601) q[0];
sx q[0];
rz(-2.2746673) q[0];
sx q[0];
rz(-2.5224199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3217818) q[2];
sx q[2];
rz(-0.95633436) q[2];
sx q[2];
rz(2.7053506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5428634) q[1];
sx q[1];
rz(-2.3485564) q[1];
sx q[1];
rz(-1.7584156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.028285154) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(2.459211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(0.49475691) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49098) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(2.0565128) q[0];
rz(-2.1814573) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(2.9575612) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40084307) q[0];
sx q[0];
rz(-1.4662192) q[0];
sx q[0];
rz(2.8507289) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5023735) q[2];
sx q[2];
rz(-0.44438617) q[2];
sx q[2];
rz(1.5803312) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9102238) q[1];
sx q[1];
rz(-1.1008881) q[1];
sx q[1];
rz(-2.4889357) q[1];
rz(0.22589485) q[3];
sx q[3];
rz(-0.73879209) q[3];
sx q[3];
rz(-1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-1.1269349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891156) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-2.254803) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(-0.14850798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200955) q[0];
sx q[0];
rz(-1.4098865) q[0];
sx q[0];
rz(-2.4077329) q[0];
x q[1];
rz(2.6988585) q[2];
sx q[2];
rz(-2.8895183) q[2];
sx q[2];
rz(-0.67539757) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.099523274) q[1];
sx q[1];
rz(-1.4617141) q[1];
sx q[1];
rz(1.3101577) q[1];
x q[2];
rz(-0.69494436) q[3];
sx q[3];
rz(-0.57999014) q[3];
sx q[3];
rz(-2.490173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85577661) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-2.6532069) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.2667806) q[3];
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
rz(-2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(2.3235902) q[0];
rz(-1.2524293) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(-1.12524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8669406) q[0];
sx q[0];
rz(-2.7526703) q[0];
sx q[0];
rz(1.2961943) q[0];
x q[1];
rz(2.0411885) q[2];
sx q[2];
rz(-2.8090968) q[2];
sx q[2];
rz(1.5770797) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6553073) q[1];
sx q[1];
rz(-0.40973046) q[1];
sx q[1];
rz(2.3182931) q[1];
rz(-pi) q[2];
rz(2.4756487) q[3];
sx q[3];
rz(-2.0332554) q[3];
sx q[3];
rz(2.8921814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(1.6623496) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(-1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97061625) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(-1.9203141) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(1.01064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0518274) q[0];
sx q[0];
rz(-0.83501378) q[0];
sx q[0];
rz(0.21659539) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0755195) q[2];
sx q[2];
rz(-1.8349378) q[2];
sx q[2];
rz(-1.2656982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6992221) q[1];
sx q[1];
rz(-1.7915627) q[1];
sx q[1];
rz(1.3352331) q[1];
rz(-1.6628077) q[3];
sx q[3];
rz(-1.9698471) q[3];
sx q[3];
rz(-1.240318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6179787) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(0.61839473) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(-0.075597413) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48329096) q[0];
sx q[0];
rz(-2.1245983) q[0];
sx q[0];
rz(2.7730586) q[0];
rz(-pi) q[1];
rz(-1.8972626) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(1.9290123) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96502111) q[1];
sx q[1];
rz(-2.2551943) q[1];
sx q[1];
rz(-1.7040764) q[1];
x q[2];
rz(-2.5209386) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33346924) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(1.8010275) q[2];
rz(-0.30424413) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(0.17679581) q[0];
rz(-1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(2.7005844) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84687418) q[0];
sx q[0];
rz(-0.28781578) q[0];
sx q[0];
rz(0.78886445) q[0];
x q[1];
rz(-2.6920325) q[2];
sx q[2];
rz(-2.9712147) q[2];
sx q[2];
rz(-0.61194387) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1383789) q[1];
sx q[1];
rz(-0.56204501) q[1];
sx q[1];
rz(-0.89488645) q[1];
rz(-pi) q[2];
rz(-1.8701843) q[3];
sx q[3];
rz(-2.1736439) q[3];
sx q[3];
rz(2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(0.30187541) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(-0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.72398913) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(-0.026731116) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-1.0529636) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(2.5857153) q[3];
sx q[3];
rz(-1.1365969) q[3];
sx q[3];
rz(2.668276) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
