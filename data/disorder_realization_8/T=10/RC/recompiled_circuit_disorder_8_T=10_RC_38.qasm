OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(4.0868563) q[0];
sx q[0];
rz(9.950369) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(-1.9089729) q[1];
sx q[1];
rz(0.90484172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5072767) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(-1.0213724) q[0];
x q[1];
rz(0.71360795) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(1.3908536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79996586) q[1];
sx q[1];
rz(-2.6744665) q[1];
sx q[1];
rz(2.1145691) q[1];
x q[2];
rz(0.18890394) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-3.0453483) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7006943) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-0.66295019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6843296) q[0];
sx q[0];
rz(-1.6329137) q[0];
sx q[0];
rz(-1.5044042) q[0];
rz(-pi) q[1];
rz(-2.0775954) q[2];
sx q[2];
rz(-0.71841824) q[2];
sx q[2];
rz(1.6990627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1124277) q[1];
sx q[1];
rz(-1.4883092) q[1];
sx q[1];
rz(-0.44177456) q[1];
rz(1.7327659) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(-3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(-2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(-0.22234017) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1919353) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(-0.8215254) q[0];
rz(2.8754183) q[2];
sx q[2];
rz(-1.4783579) q[2];
sx q[2];
rz(-0.12649525) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8232302) q[1];
sx q[1];
rz(-1.3981817) q[1];
sx q[1];
rz(-2.361972) q[1];
rz(3.0869811) q[3];
sx q[3];
rz(-1.5715989) q[3];
sx q[3];
rz(0.69566788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(-0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-0.17424507) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78631567) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(3.0062208) q[0];
rz(0.60855234) q[2];
sx q[2];
rz(-2.2366479) q[2];
sx q[2];
rz(-0.37492875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40401134) q[1];
sx q[1];
rz(-2.0239132) q[1];
sx q[1];
rz(1.416942) q[1];
rz(-pi) q[2];
rz(-1.3644049) q[3];
sx q[3];
rz(-1.6009814) q[3];
sx q[3];
rz(2.9914732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(-2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-0.61606032) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84232932) q[0];
sx q[0];
rz(-3.0780601) q[0];
sx q[0];
rz(-0.97139831) q[0];
rz(-pi) q[1];
rz(-0.36523833) q[2];
sx q[2];
rz(-1.4589981) q[2];
sx q[2];
rz(2.0888084) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7459813) q[1];
sx q[1];
rz(-0.99318722) q[1];
sx q[1];
rz(1.5773768) q[1];
rz(2.2474399) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(-1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(1.6051965) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-2.2568259) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677744) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.7675179) q[0];
rz(-pi) q[1];
x q[1];
rz(0.066263513) q[2];
sx q[2];
rz(-0.21050669) q[2];
sx q[2];
rz(-1.9312242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92237597) q[1];
sx q[1];
rz(-1.1616542) q[1];
sx q[1];
rz(-2.8273696) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8956036) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(-0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(2.4087002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3467348) q[0];
sx q[0];
rz(-1.1746527) q[0];
sx q[0];
rz(-2.5651155) q[0];
rz(-pi) q[1];
rz(0.86874666) q[2];
sx q[2];
rz(-1.1940496) q[2];
sx q[2];
rz(-1.6422611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43436189) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(-1.2039127) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(2.4550408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(0.7094267) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(-2.8628796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3942791) q[0];
sx q[0];
rz(-0.25082591) q[0];
sx q[0];
rz(-2.5670811) q[0];
x q[1];
rz(-0.3785554) q[2];
sx q[2];
rz(-0.59213973) q[2];
sx q[2];
rz(-2.2373667) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7292494) q[1];
sx q[1];
rz(-1.3530429) q[1];
sx q[1];
rz(1.2530243) q[1];
x q[2];
rz(0.078019402) q[3];
sx q[3];
rz(-1.7934414) q[3];
sx q[3];
rz(-1.5732461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30148208) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(-3.0855132) q[2];
rz(2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(-2.1726998) q[0];
rz(-2.6682207) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.1425346) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5794967) q[0];
sx q[0];
rz(-1.2764494) q[0];
sx q[0];
rz(-3.0466945) q[0];
x q[1];
rz(-3.1387024) q[2];
sx q[2];
rz(-2.4382466) q[2];
sx q[2];
rz(-0.10290111) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4559608) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(-2.901652) q[1];
x q[2];
rz(-2.2409866) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(-0.98774324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(0.3237237) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(-1.2385626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68574821) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(0.2483764) q[0];
rz(-pi) q[1];
rz(1.6221223) q[2];
sx q[2];
rz(-2.3231299) q[2];
sx q[2];
rz(2.3527956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4148256) q[1];
sx q[1];
rz(-1.4531724) q[1];
sx q[1];
rz(0.91314258) q[1];
rz(-pi) q[2];
x q[2];
rz(0.978312) q[3];
sx q[3];
rz(-1.4761792) q[3];
sx q[3];
rz(0.30944165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-2.1380077) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(-3.1162221) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(1.17169) q[2];
sx q[2];
rz(-1.8532955) q[2];
sx q[2];
rz(2.4886139) q[2];
rz(-0.17037114) q[3];
sx q[3];
rz(-0.80633612) q[3];
sx q[3];
rz(2.7663305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
