OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(1.6341524) q[0];
rz(-1.545067) q[1];
sx q[1];
rz(-2.5453321) q[1];
sx q[1];
rz(2.526386) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1546254) q[0];
sx q[0];
rz(-1.1095424) q[0];
sx q[0];
rz(2.2485562) q[0];
rz(-pi) q[1];
rz(0.35742128) q[2];
sx q[2];
rz(-1.3961785) q[2];
sx q[2];
rz(1.071196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33949172) q[1];
sx q[1];
rz(-1.9994945) q[1];
sx q[1];
rz(2.2080253) q[1];
x q[2];
rz(-2.8350713) q[3];
sx q[3];
rz(-1.3984826) q[3];
sx q[3];
rz(2.4337208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(2.8033076) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.5240086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.679927) q[0];
sx q[0];
rz(-1.3711509) q[0];
sx q[0];
rz(-3.1398849) q[0];
rz(-pi) q[1];
x q[1];
rz(3.121071) q[2];
sx q[2];
rz(-1.1164718) q[2];
sx q[2];
rz(-0.12873912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6807032) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(-2.1417888) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3832983) q[3];
sx q[3];
rz(-1.6631931) q[3];
sx q[3];
rz(2.932991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(-0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-0.31578627) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-2.8895203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4323498) q[0];
sx q[0];
rz(-2.3253257) q[0];
sx q[0];
rz(0.66246756) q[0];
rz(-pi) q[1];
rz(-0.71528541) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(-2.8298024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6985059) q[1];
sx q[1];
rz(-1.7585187) q[1];
sx q[1];
rz(-2.1952573) q[1];
x q[2];
rz(-2.7470845) q[3];
sx q[3];
rz(-1.2393701) q[3];
sx q[3];
rz(3.1098207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0198274) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(1.6148286) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93543816) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(-2.6576256) q[0];
rz(0.25755067) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(-0.79007733) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45338079) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(-1.3761671) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1726923) q[3];
sx q[3];
rz(-0.36877353) q[3];
sx q[3];
rz(0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84918555) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-2.6814931) q[2];
rz(-1.7442616) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(2.8869693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.03527) q[0];
sx q[0];
rz(-1.584504) q[0];
sx q[0];
rz(1.4916219) q[0];
rz(-pi) q[1];
rz(3.1270199) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(0.84469675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9533206) q[1];
sx q[1];
rz(-2.4362262) q[1];
sx q[1];
rz(-0.377368) q[1];
x q[2];
rz(-0.43318627) q[3];
sx q[3];
rz(-0.90900366) q[3];
sx q[3];
rz(-1.2695241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(1.0866603) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-0.20656955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26039133) q[0];
sx q[0];
rz(-0.83575373) q[0];
sx q[0];
rz(-0.3221237) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39482306) q[2];
sx q[2];
rz(-2.0645112) q[2];
sx q[2];
rz(-2.0815108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0530015) q[1];
sx q[1];
rz(-0.87391657) q[1];
sx q[1];
rz(-2.9104396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8387186) q[3];
sx q[3];
rz(-0.93749638) q[3];
sx q[3];
rz(0.57297046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(-0.28277961) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(-2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151826) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.8136224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8570003) q[0];
sx q[0];
rz(-2.6866331) q[0];
sx q[0];
rz(1.8332464) q[0];
rz(-pi) q[1];
rz(-2.1778657) q[2];
sx q[2];
rz(-0.71787314) q[2];
sx q[2];
rz(-1.3340064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3634062) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(-2.2560675) q[1];
rz(-pi) q[2];
rz(-0.51289576) q[3];
sx q[3];
rz(-1.9478056) q[3];
sx q[3];
rz(0.60037724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.7810812) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(0.95091933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61688214) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(0.25047238) q[0];
x q[1];
rz(-0.062336246) q[2];
sx q[2];
rz(-1.8364292) q[2];
sx q[2];
rz(3.0692284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9329405) q[1];
sx q[1];
rz(-1.3622074) q[1];
sx q[1];
rz(-1.5044466) q[1];
x q[2];
rz(1.1675646) q[3];
sx q[3];
rz(-1.9208761) q[3];
sx q[3];
rz(-1.10266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(-0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(2.7245522) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(0.65972796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035559) q[0];
sx q[0];
rz(-0.8350026) q[0];
sx q[0];
rz(0.8291709) q[0];
rz(2.9218036) q[2];
sx q[2];
rz(-2.347749) q[2];
sx q[2];
rz(1.9915875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5539726) q[1];
sx q[1];
rz(-2.4269322) q[1];
sx q[1];
rz(1.1187394) q[1];
rz(-0.057007313) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(-2.1155604) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.8359258) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(-1.0356888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4982088) q[0];
sx q[0];
rz(-1.9783101) q[0];
sx q[0];
rz(-1.9252752) q[0];
x q[1];
rz(-0.91949384) q[2];
sx q[2];
rz(-0.80923015) q[2];
sx q[2];
rz(-0.86916718) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2346238) q[1];
sx q[1];
rz(-1.8112438) q[1];
sx q[1];
rz(0.30827) q[1];
rz(-pi) q[2];
rz(-0.74929897) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(-1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.37968996) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-0.41082906) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(-1.5796173) q[3];
sx q[3];
rz(-2.1891441) q[3];
sx q[3];
rz(1.7563663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
