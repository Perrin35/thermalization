OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(-1.0426961) q[0];
sx q[0];
rz(0.27867499) q[0];
rz(-1.9947808) q[1];
sx q[1];
rz(-2.6169701) q[1];
sx q[1];
rz(1.2710748) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4500101) q[0];
sx q[0];
rz(-1.5890536) q[0];
sx q[0];
rz(-0.26055793) q[0];
rz(-0.71670436) q[2];
sx q[2];
rz(-0.19858805) q[2];
sx q[2];
rz(0.27825296) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48285741) q[1];
sx q[1];
rz(-2.4244011) q[1];
sx q[1];
rz(1.4499922) q[1];
x q[2];
rz(0.75436169) q[3];
sx q[3];
rz(-0.76243692) q[3];
sx q[3];
rz(-2.0319394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.110454) q[2];
sx q[2];
rz(-1.217507) q[2];
sx q[2];
rz(0.26220599) q[2];
rz(1.0555222) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6156886) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(2.9600034) q[0];
rz(-1.7407821) q[1];
sx q[1];
rz(-1.1927651) q[1];
sx q[1];
rz(-0.84091944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0706874) q[0];
sx q[0];
rz(-1.5884261) q[0];
sx q[0];
rz(0.32004642) q[0];
x q[1];
rz(-2.1757751) q[2];
sx q[2];
rz(-0.73558319) q[2];
sx q[2];
rz(-0.090026131) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.629417) q[1];
sx q[1];
rz(-1.7390774) q[1];
sx q[1];
rz(3.1049033) q[1];
rz(-pi) q[2];
rz(-2.7737229) q[3];
sx q[3];
rz(-2.1506393) q[3];
sx q[3];
rz(0.90880064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38222868) q[2];
sx q[2];
rz(-0.21491924) q[2];
sx q[2];
rz(1.0648897) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.3353142) q[3];
sx q[3];
rz(-1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29838022) q[0];
sx q[0];
rz(-1.9623663) q[0];
sx q[0];
rz(2.9425353) q[0];
rz(0.82636034) q[1];
sx q[1];
rz(-0.91824707) q[1];
sx q[1];
rz(0.77702776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26829241) q[0];
sx q[0];
rz(-1.7275294) q[0];
sx q[0];
rz(3.0746033) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1146745) q[2];
sx q[2];
rz(-1.9444398) q[2];
sx q[2];
rz(-2.0622562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4118581) q[1];
sx q[1];
rz(-2.5246906) q[1];
sx q[1];
rz(-2.0286948) q[1];
x q[2];
rz(-1.1714524) q[3];
sx q[3];
rz(-0.63217406) q[3];
sx q[3];
rz(-0.20242385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1021154) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(-2.6285505) q[2];
rz(0.78220621) q[3];
sx q[3];
rz(-1.9358044) q[3];
sx q[3];
rz(-1.7647083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9272598) q[0];
sx q[0];
rz(-1.5981277) q[0];
sx q[0];
rz(1.6511035) q[0];
rz(-0.53228846) q[1];
sx q[1];
rz(-2.5106301) q[1];
sx q[1];
rz(-1.6040364) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0396932) q[0];
sx q[0];
rz(-2.9465775) q[0];
sx q[0];
rz(-2.5620417) q[0];
rz(1.0245778) q[2];
sx q[2];
rz(-1.4909296) q[2];
sx q[2];
rz(-0.090674222) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2080431) q[1];
sx q[1];
rz(-1.4999823) q[1];
sx q[1];
rz(0.93472247) q[1];
rz(-pi) q[2];
rz(0.38942899) q[3];
sx q[3];
rz(-0.5618605) q[3];
sx q[3];
rz(-1.6461247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6060467) q[2];
sx q[2];
rz(-0.54680768) q[2];
sx q[2];
rz(0.01288506) q[2];
rz(1.2031201) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(2.119868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0032229) q[0];
sx q[0];
rz(-2.5165181) q[0];
sx q[0];
rz(-1.7417997) q[0];
rz(0.35203448) q[1];
sx q[1];
rz(-1.0540009) q[1];
sx q[1];
rz(-0.83784109) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45233881) q[0];
sx q[0];
rz(-1.4188915) q[0];
sx q[0];
rz(2.332973) q[0];
rz(-0.95372405) q[2];
sx q[2];
rz(-1.0487818) q[2];
sx q[2];
rz(-1.2921289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9329233) q[1];
sx q[1];
rz(-1.0074378) q[1];
sx q[1];
rz(-1.7209919) q[1];
x q[2];
rz(0.75148186) q[3];
sx q[3];
rz(-0.22322907) q[3];
sx q[3];
rz(-0.022836784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1167404) q[2];
sx q[2];
rz(-2.1768673) q[2];
sx q[2];
rz(1.1836729) q[2];
rz(2.8708894) q[3];
sx q[3];
rz(-2.201122) q[3];
sx q[3];
rz(1.6089599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026641332) q[0];
sx q[0];
rz(-0.6237492) q[0];
sx q[0];
rz(-2.5624516) q[0];
rz(-1.2222414) q[1];
sx q[1];
rz(-1.3074343) q[1];
sx q[1];
rz(-2.0319669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1117119) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(-0.16686186) q[0];
x q[1];
rz(-1.6272306) q[2];
sx q[2];
rz(-1.2662958) q[2];
sx q[2];
rz(0.79848189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.077495726) q[1];
sx q[1];
rz(-1.9019097) q[1];
sx q[1];
rz(1.8331794) q[1];
rz(-pi) q[2];
rz(-2.6251276) q[3];
sx q[3];
rz(-2.0544009) q[3];
sx q[3];
rz(-0.91396871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9455202) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(2.784139) q[2];
rz(-1.9516021) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(-1.801871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5424407) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(0.79208148) q[0];
rz(-0.7041086) q[1];
sx q[1];
rz(-0.45583615) q[1];
sx q[1];
rz(-1.5584996) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7501935) q[0];
sx q[0];
rz(-1.1066086) q[0];
sx q[0];
rz(1.5675005) q[0];
rz(0.81659813) q[2];
sx q[2];
rz(-2.3629521) q[2];
sx q[2];
rz(-1.1731847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4729135) q[1];
sx q[1];
rz(-1.1448507) q[1];
sx q[1];
rz(-1.5329453) q[1];
x q[2];
rz(-2.9005403) q[3];
sx q[3];
rz(-0.76605443) q[3];
sx q[3];
rz(-2.0360006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.030297) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(3.1374068) q[2];
rz(0.26675102) q[3];
sx q[3];
rz(-1.5240074) q[3];
sx q[3];
rz(2.4193616) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5608212) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(-0.15604493) q[0];
rz(1.5566298) q[1];
sx q[1];
rz(-0.86235756) q[1];
sx q[1];
rz(-2.5468266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9879919) q[0];
sx q[0];
rz(-0.56454851) q[0];
sx q[0];
rz(-0.49086824) q[0];
rz(-pi) q[1];
rz(-0.34551106) q[2];
sx q[2];
rz(-0.72241579) q[2];
sx q[2];
rz(1.9986716) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8957828) q[1];
sx q[1];
rz(-1.0059662) q[1];
sx q[1];
rz(-0.71706949) q[1];
rz(-2.4202926) q[3];
sx q[3];
rz(-2.2709322) q[3];
sx q[3];
rz(2.457778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5742089) q[2];
sx q[2];
rz(-1.5790066) q[2];
sx q[2];
rz(1.0279921) q[2];
rz(-1.3422286) q[3];
sx q[3];
rz(-2.4864311) q[3];
sx q[3];
rz(-1.3348234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-2.8780355) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(-2.3433319) q[0];
rz(2.3410666) q[1];
sx q[1];
rz(-2.6526178) q[1];
sx q[1];
rz(-0.54623234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8432586) q[0];
sx q[0];
rz(-0.97163659) q[0];
sx q[0];
rz(0.13183044) q[0];
x q[1];
rz(-0.58878657) q[2];
sx q[2];
rz(-1.2079888) q[2];
sx q[2];
rz(0.44315674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47910467) q[1];
sx q[1];
rz(-2.4291385) q[1];
sx q[1];
rz(-0.79378788) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4978374) q[3];
sx q[3];
rz(-1.3449418) q[3];
sx q[3];
rz(-2.7166004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3937248) q[2];
sx q[2];
rz(-2.0515029) q[2];
sx q[2];
rz(-0.063035034) q[2];
rz(2.4437599) q[3];
sx q[3];
rz(-2.883426) q[3];
sx q[3];
rz(-0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460655) q[0];
sx q[0];
rz(-0.87764469) q[0];
sx q[0];
rz(-0.416042) q[0];
rz(-1.7795732) q[1];
sx q[1];
rz(-2.0174618) q[1];
sx q[1];
rz(1.8488041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0344311) q[0];
sx q[0];
rz(-0.35625544) q[0];
sx q[0];
rz(-0.70888575) q[0];
rz(-pi) q[1];
rz(-1.2913537) q[2];
sx q[2];
rz(-1.1300753) q[2];
sx q[2];
rz(-0.27383495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9465465) q[1];
sx q[1];
rz(-1.0302231) q[1];
sx q[1];
rz(-0.54624301) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8804194) q[3];
sx q[3];
rz(-1.1615442) q[3];
sx q[3];
rz(-0.60200426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68299874) q[2];
sx q[2];
rz(-1.3364044) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(1.2609743) q[3];
sx q[3];
rz(-0.57456273) q[3];
sx q[3];
rz(-1.449409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64611971) q[0];
sx q[0];
rz(-1.4541805) q[0];
sx q[0];
rz(-2.4021586) q[0];
rz(-2.6932035) q[1];
sx q[1];
rz(-1.8012451) q[1];
sx q[1];
rz(-1.9582122) q[1];
rz(-2.420633) q[2];
sx q[2];
rz(-0.45476457) q[2];
sx q[2];
rz(-1.5130629) q[2];
rz(1.4215333) q[3];
sx q[3];
rz(-2.4076318) q[3];
sx q[3];
rz(-2.9416495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
