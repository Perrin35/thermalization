OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8747099) q[0];
sx q[0];
rz(-2.0800135) q[0];
sx q[0];
rz(0.51577407) q[0];
rz(0.90142673) q[1];
sx q[1];
rz(2.9401448) q[1];
sx q[1];
rz(10.425865) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62681055) q[0];
sx q[0];
rz(-2.0603466) q[0];
sx q[0];
rz(2.5077257) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7269709) q[2];
sx q[2];
rz(-1.5650563) q[2];
sx q[2];
rz(1.8513377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7202819) q[1];
sx q[1];
rz(-1.1291468) q[1];
sx q[1];
rz(0.71278211) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68925011) q[3];
sx q[3];
rz(-0.19016506) q[3];
sx q[3];
rz(-2.3306757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.72656816) q[2];
sx q[2];
rz(-2.6308306) q[2];
sx q[2];
rz(1.4600352) q[2];
rz(2.7728752) q[3];
sx q[3];
rz(-1.5548778) q[3];
sx q[3];
rz(1.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52861315) q[0];
sx q[0];
rz(-3.0280085) q[0];
sx q[0];
rz(-0.088264912) q[0];
rz(-2.2093692) q[1];
sx q[1];
rz(-2.014522) q[1];
sx q[1];
rz(0.50311911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022633502) q[0];
sx q[0];
rz(-1.6741856) q[0];
sx q[0];
rz(-0.4876644) q[0];
rz(2.2406949) q[2];
sx q[2];
rz(-1.7187727) q[2];
sx q[2];
rz(2.072538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4331665) q[1];
sx q[1];
rz(-1.462877) q[1];
sx q[1];
rz(-0.34924653) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3130619) q[3];
sx q[3];
rz(-1.9779543) q[3];
sx q[3];
rz(0.26155868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3162389) q[2];
sx q[2];
rz(-0.97430682) q[2];
sx q[2];
rz(2.6884354) q[2];
rz(0.00099269021) q[3];
sx q[3];
rz(-1.5627292) q[3];
sx q[3];
rz(2.4003975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93333018) q[0];
sx q[0];
rz(-2.9587511) q[0];
sx q[0];
rz(0.46874794) q[0];
rz(-2.5543429) q[1];
sx q[1];
rz(-1.8533555) q[1];
sx q[1];
rz(0.25150484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6856988) q[0];
sx q[0];
rz(-2.532428) q[0];
sx q[0];
rz(-0.76700155) q[0];
rz(-pi) q[1];
rz(1.7825039) q[2];
sx q[2];
rz(-1.1897174) q[2];
sx q[2];
rz(-0.43109387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9589825) q[1];
sx q[1];
rz(-1.512456) q[1];
sx q[1];
rz(-1.6946409) q[1];
x q[2];
rz(-1.92177) q[3];
sx q[3];
rz(-1.1535221) q[3];
sx q[3];
rz(-1.1116416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59715366) q[2];
sx q[2];
rz(-0.73828283) q[2];
sx q[2];
rz(3.0944589) q[2];
rz(-0.86826396) q[3];
sx q[3];
rz(-1.2406113) q[3];
sx q[3];
rz(0.49938437) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.353001) q[0];
sx q[0];
rz(-2.7356739) q[0];
sx q[0];
rz(-1.3141919) q[0];
rz(0.81002533) q[1];
sx q[1];
rz(-1.6255197) q[1];
sx q[1];
rz(-2.8257418) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11498904) q[0];
sx q[0];
rz(-0.79885495) q[0];
sx q[0];
rz(0.82180114) q[0];
x q[1];
rz(2.8239657) q[2];
sx q[2];
rz(-1.4072573) q[2];
sx q[2];
rz(2.7224772) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0210353) q[1];
sx q[1];
rz(-1.2791954) q[1];
sx q[1];
rz(-2.3910644) q[1];
rz(-2.9258435) q[3];
sx q[3];
rz(-0.61601725) q[3];
sx q[3];
rz(-0.3122789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6541859) q[2];
sx q[2];
rz(-2.5155289) q[2];
sx q[2];
rz(0.89228863) q[2];
rz(-1.6728632) q[3];
sx q[3];
rz(-1.059831) q[3];
sx q[3];
rz(1.0631801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6650218) q[0];
sx q[0];
rz(-0.63711089) q[0];
sx q[0];
rz(2.2549905) q[0];
rz(-2.2992112) q[1];
sx q[1];
rz(-1.8319172) q[1];
sx q[1];
rz(1.1096035) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32217596) q[0];
sx q[0];
rz(-0.79335326) q[0];
sx q[0];
rz(2.1920836) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72317041) q[2];
sx q[2];
rz(-1.4338655) q[2];
sx q[2];
rz(-1.3201158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.90613104) q[1];
sx q[1];
rz(-2.0812199) q[1];
sx q[1];
rz(-2.4861369) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1704386) q[3];
sx q[3];
rz(-2.1622938) q[3];
sx q[3];
rz(2.9324071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2716918) q[2];
sx q[2];
rz(-0.39125189) q[2];
sx q[2];
rz(2.5539577) q[2];
rz(-2.9247126) q[3];
sx q[3];
rz(-1.2809332) q[3];
sx q[3];
rz(1.3542401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13123913) q[0];
sx q[0];
rz(-2.6709747) q[0];
sx q[0];
rz(1.7577897) q[0];
rz(-0.17419392) q[1];
sx q[1];
rz(-1.1015588) q[1];
sx q[1];
rz(-1.1879638) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9458141) q[0];
sx q[0];
rz(-1.6530237) q[0];
sx q[0];
rz(-1.8503227) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8415478) q[2];
sx q[2];
rz(-2.6159673) q[2];
sx q[2];
rz(-2.2208461) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9932258) q[1];
sx q[1];
rz(-1.6794551) q[1];
sx q[1];
rz(-3.1106366) q[1];
rz(-pi) q[2];
rz(-0.31962304) q[3];
sx q[3];
rz(-0.79053662) q[3];
sx q[3];
rz(2.9755693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69372988) q[2];
sx q[2];
rz(-1.8665946) q[2];
sx q[2];
rz(2.0873783) q[2];
rz(0.55274719) q[3];
sx q[3];
rz(-2.882143) q[3];
sx q[3];
rz(-2.0235846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0746821) q[0];
sx q[0];
rz(-2.5801165) q[0];
sx q[0];
rz(-0.018360227) q[0];
rz(-1.297599) q[1];
sx q[1];
rz(-1.8074139) q[1];
sx q[1];
rz(1.1150572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92662382) q[0];
sx q[0];
rz(-1.4415662) q[0];
sx q[0];
rz(1.5530759) q[0];
rz(1.1184022) q[2];
sx q[2];
rz(-1.2873189) q[2];
sx q[2];
rz(0.3800791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9950986) q[1];
sx q[1];
rz(-2.0420688) q[1];
sx q[1];
rz(-2.5033094) q[1];
rz(-pi) q[2];
rz(-1.0706022) q[3];
sx q[3];
rz(-1.9755441) q[3];
sx q[3];
rz(-0.91171748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6899507) q[2];
sx q[2];
rz(-2.3929598) q[2];
sx q[2];
rz(-0.79317036) q[2];
rz(-0.66148174) q[3];
sx q[3];
rz(-1.7659148) q[3];
sx q[3];
rz(-0.64600265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20872214) q[0];
sx q[0];
rz(-1.1247617) q[0];
sx q[0];
rz(0.17054184) q[0];
rz(-1.0199245) q[1];
sx q[1];
rz(-0.41671697) q[1];
sx q[1];
rz(-0.65933093) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070674226) q[0];
sx q[0];
rz(-3.0844581) q[0];
sx q[0];
rz(0.26483421) q[0];
rz(-pi) q[1];
rz(1.061672) q[2];
sx q[2];
rz(-1.1069555) q[2];
sx q[2];
rz(1.0373751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0773191) q[1];
sx q[1];
rz(-0.96076316) q[1];
sx q[1];
rz(1.9404821) q[1];
x q[2];
rz(-2.3850494) q[3];
sx q[3];
rz(-1.1361101) q[3];
sx q[3];
rz(-0.68090465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8352123) q[2];
sx q[2];
rz(-2.7472718) q[2];
sx q[2];
rz(2.1739056) q[2];
rz(0.59099284) q[3];
sx q[3];
rz(-1.7715745) q[3];
sx q[3];
rz(-1.7041357) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0377334) q[0];
sx q[0];
rz(-0.93872207) q[0];
sx q[0];
rz(0.18873225) q[0];
rz(1.0602779) q[1];
sx q[1];
rz(-2.4791368) q[1];
sx q[1];
rz(-0.8998543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2162595) q[0];
sx q[0];
rz(-1.0311613) q[0];
sx q[0];
rz(-0.58109053) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84328358) q[2];
sx q[2];
rz(-1.5688063) q[2];
sx q[2];
rz(-0.28283027) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2267188) q[1];
sx q[1];
rz(-0.9203099) q[1];
sx q[1];
rz(-1.5730251) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7873618) q[3];
sx q[3];
rz(-2.4087071) q[3];
sx q[3];
rz(2.2115744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4047644) q[2];
sx q[2];
rz(-1.2238203) q[2];
sx q[2];
rz(1.9723816) q[2];
rz(2.2660008) q[3];
sx q[3];
rz(-2.4647522) q[3];
sx q[3];
rz(-3.0524047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5065696) q[0];
sx q[0];
rz(-1.6825786) q[0];
sx q[0];
rz(-2.7154679) q[0];
rz(-1.2395073) q[1];
sx q[1];
rz(-1.0303717) q[1];
sx q[1];
rz(-0.32136163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023055619) q[0];
sx q[0];
rz(-1.5154953) q[0];
sx q[0];
rz(-2.1441516) q[0];
x q[1];
rz(0.05409492) q[2];
sx q[2];
rz(-1.128643) q[2];
sx q[2];
rz(-1.4070321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1758986) q[1];
sx q[1];
rz(-1.4280609) q[1];
sx q[1];
rz(-1.229415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57182248) q[3];
sx q[3];
rz(-0.79872978) q[3];
sx q[3];
rz(-3.1099936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1539803) q[2];
sx q[2];
rz(-2.2025547) q[2];
sx q[2];
rz(-0.067848094) q[2];
rz(0.96505729) q[3];
sx q[3];
rz(-2.8128251) q[3];
sx q[3];
rz(-1.7353826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168552) q[0];
sx q[0];
rz(-0.4586093) q[0];
sx q[0];
rz(-2.1427857) q[0];
rz(0.85159341) q[1];
sx q[1];
rz(-1.0706182) q[1];
sx q[1];
rz(0.65382438) q[1];
rz(-0.10112496) q[2];
sx q[2];
rz(-0.55872266) q[2];
sx q[2];
rz(-1.2324294) q[2];
rz(3.1134277) q[3];
sx q[3];
rz(-2.4234301) q[3];
sx q[3];
rz(2.1860893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
