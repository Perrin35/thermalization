OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(-2.5807227) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1808704) q[0];
sx q[0];
rz(-1.3836432) q[0];
sx q[0];
rz(-3.0357009) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7772929) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(-0.42687624) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3455968) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(2.7624346) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0662765) q[3];
sx q[3];
rz(-0.30116943) q[3];
sx q[3];
rz(1.0217713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(2.9585178) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7586655) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(1.7090319) q[0];
x q[1];
rz(2.9039731) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(-1.876229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(-0.19660463) q[1];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(-1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(0.69491274) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-0.55535299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99291891) q[0];
sx q[0];
rz(-2.2745471) q[0];
sx q[0];
rz(0.91627319) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.136134) q[2];
sx q[2];
rz(-1.3034504) q[2];
sx q[2];
rz(0.29758673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33298102) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(-1.8379184) q[1];
rz(-pi) q[2];
rz(1.8954574) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(-0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(1.4499433) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-0.81480169) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(-2.8864158) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53084757) q[0];
sx q[0];
rz(-1.8662211) q[0];
sx q[0];
rz(-2.7644964) q[0];
x q[1];
rz(1.9011263) q[2];
sx q[2];
rz(-2.1836046) q[2];
sx q[2];
rz(1.4532879) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5863122) q[1];
sx q[1];
rz(-1.3603856) q[1];
sx q[1];
rz(1.2661238) q[1];
x q[2];
rz(-2.0917986) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(1.1417768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(-3.0854991) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1306886) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(2.5581193) q[0];
rz(-pi) q[1];
rz(-0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(2.3914571) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6519421) q[1];
sx q[1];
rz(-2.838476) q[1];
sx q[1];
rz(-3.0928844) q[1];
rz(-2.131093) q[3];
sx q[3];
rz(-1.6464525) q[3];
sx q[3];
rz(-1.2922985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19676767) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(-1.0374271) q[0];
rz(3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(-2.6882753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33659014) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(2.5617983) q[1];
x q[2];
rz(-0.16320634) q[3];
sx q[3];
rz(-1.8582134) q[3];
sx q[3];
rz(-2.4678469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(1.3180805) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(0.61002237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22247032) q[0];
sx q[0];
rz(-1.3195992) q[0];
sx q[0];
rz(-0.50719502) q[0];
rz(-pi) q[1];
rz(1.7213983) q[2];
sx q[2];
rz(-2.4820231) q[2];
sx q[2];
rz(1.2130376) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89027379) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(-0.031884738) q[1];
x q[2];
rz(1.9500908) q[3];
sx q[3];
rz(-1.7625574) q[3];
sx q[3];
rz(-2.3281043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7685331) q[0];
sx q[0];
rz(-0.031950843) q[0];
sx q[0];
rz(1.91747) q[0];
x q[1];
rz(1.2484776) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(1.5600187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(-2.9176941) q[1];
x q[2];
rz(-0.96111091) q[3];
sx q[3];
rz(-0.52934066) q[3];
sx q[3];
rz(-1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-2.396778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54430994) q[0];
sx q[0];
rz(-2.3047857) q[0];
sx q[0];
rz(-2.701176) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9717734) q[2];
sx q[2];
rz(-0.81911659) q[2];
sx q[2];
rz(-0.70105201) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71036584) q[1];
sx q[1];
rz(-1.2214298) q[1];
sx q[1];
rz(-2.9037895) q[1];
rz(1.1780147) q[3];
sx q[3];
rz(-0.82735705) q[3];
sx q[3];
rz(-0.82069293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(0.40400305) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.9706479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0572646) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(1.5945934) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8413713) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(-0.68683456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78632894) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(-1.5934056) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42697866) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(-0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(-0.099427632) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(0.70384937) q[2];
sx q[2];
rz(-0.8740295) q[2];
sx q[2];
rz(-0.340273) q[2];
rz(3.1254461) q[3];
sx q[3];
rz(-1.2373677) q[3];
sx q[3];
rz(2.2131372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
