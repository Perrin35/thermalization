OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(1.3134365) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9986721) q[2];
sx q[2];
rz(-1.9170205) q[2];
sx q[2];
rz(1.2905754) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3798843) q[1];
sx q[1];
rz(-2.0239081) q[1];
sx q[1];
rz(-2.3945827) q[1];
rz(1.9786505) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(-0.70139635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59149867) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(-0.47810289) q[2];
rz(1.452662) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1801382) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(-1.0189198) q[0];
rz(-1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(2.5085124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2791726) q[0];
sx q[0];
rz(-0.76403585) q[0];
sx q[0];
rz(-0.11636244) q[0];
rz(-pi) q[1];
x q[1];
rz(0.622153) q[2];
sx q[2];
rz(-0.85217798) q[2];
sx q[2];
rz(0.46301401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1795579) q[1];
sx q[1];
rz(-1.9704559) q[1];
sx q[1];
rz(2.7707997) q[1];
x q[2];
rz(-0.31140621) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(-0.11745545) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(2.450768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4118977) q[0];
sx q[0];
rz(-1.6459961) q[0];
sx q[0];
rz(-0.15326432) q[0];
rz(2.8969953) q[2];
sx q[2];
rz(-2.3615712) q[2];
sx q[2];
rz(-0.56842677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1289127) q[1];
sx q[1];
rz(-1.7528273) q[1];
sx q[1];
rz(0.026489594) q[1];
rz(-pi) q[2];
rz(2.7615943) q[3];
sx q[3];
rz(-1.2424386) q[3];
sx q[3];
rz(-0.85193714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(-3.0351191) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(-1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0687662) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(-2.6191214) q[0];
rz(-2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(0.55363384) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.583601) q[0];
sx q[0];
rz(-0.86692536) q[0];
sx q[0];
rz(0.61917275) q[0];
rz(-2.3217818) q[2];
sx q[2];
rz(-2.1852583) q[2];
sx q[2];
rz(-0.43624207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5987293) q[1];
sx q[1];
rz(-2.3485564) q[1];
sx q[1];
rz(-1.383177) q[1];
rz(-pi) q[2];
rz(-1.5418566) q[3];
sx q[3];
rz(-0.77416285) q[3];
sx q[3];
rz(0.72284568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(-2.2385712) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(-1.2899227) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(2.0565128) q[0];
rz(-2.1814573) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(-2.9575612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40084307) q[0];
sx q[0];
rz(-1.6753734) q[0];
sx q[0];
rz(-0.29086374) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2940302) q[2];
sx q[2];
rz(-1.2185316) q[2];
sx q[2];
rz(0.89154746) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9102238) q[1];
sx q[1];
rz(-2.0407045) q[1];
sx q[1];
rz(0.652657) q[1];
rz(-pi) q[2];
rz(0.22589485) q[3];
sx q[3];
rz(-0.73879209) q[3];
sx q[3];
rz(1.5238638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(2.0007755) q[2];
rz(-0.42282894) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-1.1269349) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(2.9930847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1214971) q[0];
sx q[0];
rz(-1.7317061) q[0];
sx q[0];
rz(2.4077329) q[0];
rz(-pi) q[1];
rz(-0.44273419) q[2];
sx q[2];
rz(-0.25207439) q[2];
sx q[2];
rz(0.67539757) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6993466) q[1];
sx q[1];
rz(-1.3117426) q[1];
sx q[1];
rz(-3.0287292) q[1];
rz(-pi) q[2];
rz(-2.6753747) q[3];
sx q[3];
rz(-1.9293474) q[3];
sx q[3];
rz(2.8311604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(0.4883858) q[2];
rz(0.48940247) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(-1.874812) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69328904) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(1.12524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624807) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(0.11066779) q[0];
x q[1];
rz(-0.15525012) q[2];
sx q[2];
rz(-1.2755738) q[2];
sx q[2];
rz(1.0709907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.352467) q[1];
sx q[1];
rz(-1.8450292) q[1];
sx q[1];
rz(-1.8791566) q[1];
rz(-pi) q[2];
rz(2.1358228) q[3];
sx q[3];
rz(-2.1566448) q[3];
sx q[3];
rz(0.98423959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(2.4970064) q[2];
rz(1.6623496) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1709764) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-2.1309526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.772951) q[0];
sx q[0];
rz(-2.3803582) q[0];
sx q[0];
rz(1.8038473) q[0];
rz(-pi) q[1];
rz(2.8433617) q[2];
sx q[2];
rz(-1.094162) q[2];
sx q[2];
rz(2.6963866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.752754) q[1];
sx q[1];
rz(-0.32143444) q[1];
sx q[1];
rz(2.3366117) q[1];
rz(0.2145433) q[3];
sx q[3];
rz(-0.40896591) q[3];
sx q[3];
rz(-1.6681125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(-0.56898919) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(-2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(-0.49474299) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9915174) q[0];
sx q[0];
rz(-0.65438327) q[0];
sx q[0];
rz(2.0983216) q[0];
rz(-pi) q[1];
rz(0.059968791) q[2];
sx q[2];
rz(-1.7459933) q[2];
sx q[2];
rz(0.8808459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9675933) q[1];
sx q[1];
rz(-0.6951957) q[1];
sx q[1];
rz(-2.980152) q[1];
rz(2.5209386) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(1.8295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33346924) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(1.8010275) q[2];
rz(-2.8373485) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-0.17679581) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(-0.44100824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6498043) q[0];
sx q[0];
rz(-1.7735964) q[0];
sx q[0];
rz(0.2056528) q[0];
rz(-pi) q[1];
rz(0.15372865) q[2];
sx q[2];
rz(-1.6445465) q[2];
sx q[2];
rz(-1.7388625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0032138) q[1];
sx q[1];
rz(-0.56204501) q[1];
sx q[1];
rz(-2.2467062) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62426626) q[3];
sx q[3];
rz(-1.3253951) q[3];
sx q[3];
rz(1.0898958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5796154) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(0.026731116) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(-0.84457196) q[2];
sx q[2];
rz(-1.209134) q[2];
sx q[2];
rz(-0.67187885) q[2];
rz(1.0711014) q[3];
sx q[3];
rz(-2.069996) q[3];
sx q[3];
rz(1.3531006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
