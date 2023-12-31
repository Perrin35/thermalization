OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1188388) q[0];
sx q[0];
rz(-1.7099329) q[0];
sx q[0];
rz(-1.820178) q[0];
rz(0.42135294) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(-0.28819627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(-0.63912649) q[1];
rz(-pi) q[2];
rz(2.4217371) q[3];
sx q[3];
rz(-2.8350283) q[3];
sx q[3];
rz(0.43678624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(2.485086) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(-0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(-0.48746902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90855234) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(-3.0768865) q[0];
rz(-pi) q[1];
rz(0.80584731) q[2];
sx q[2];
rz(-0.57493756) q[2];
sx q[2];
rz(-0.88428674) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.060354787) q[1];
sx q[1];
rz(-1.3831257) q[1];
sx q[1];
rz(2.9802122) q[1];
rz(-pi) q[2];
rz(-1.3120193) q[3];
sx q[3];
rz(-0.32734713) q[3];
sx q[3];
rz(0.60206383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6762786) q[0];
sx q[0];
rz(-1.9440117) q[0];
sx q[0];
rz(-2.7200384) q[0];
rz(-1.1995302) q[2];
sx q[2];
rz(-1.7811333) q[2];
sx q[2];
rz(0.51031761) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0157156) q[1];
sx q[1];
rz(-1.5029969) q[1];
sx q[1];
rz(-0.61342872) q[1];
x q[2];
rz(1.907903) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(-1.4032723) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.614914) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(0.99292361) q[0];
x q[1];
rz(0.36107365) q[2];
sx q[2];
rz(-1.6517342) q[2];
sx q[2];
rz(2.0780448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4184119) q[1];
sx q[1];
rz(-2.3310117) q[1];
sx q[1];
rz(2.4852738) q[1];
x q[2];
rz(0.50587378) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(0.4531025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950726) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(-1.1016269) q[0];
x q[1];
rz(1.4278973) q[2];
sx q[2];
rz(-1.60534) q[2];
sx q[2];
rz(-2.3847716) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0799775) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(1.7916937) q[1];
x q[2];
rz(0.971332) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0129464) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.5135117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3753525) q[0];
sx q[0];
rz(-1.1002812) q[0];
sx q[0];
rz(-1.2067632) q[0];
rz(-1.4432625) q[2];
sx q[2];
rz(-0.45293929) q[2];
sx q[2];
rz(-0.73308257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39735079) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(2.1703297) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.120818) q[3];
sx q[3];
rz(-1.9213772) q[3];
sx q[3];
rz(-1.1603125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.05904077) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0156353) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.1869173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163527) q[0];
sx q[0];
rz(-1.2901187) q[0];
sx q[0];
rz(1.0513845) q[0];
x q[1];
rz(-1.600012) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(0.37994775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.742332) q[1];
sx q[1];
rz(-0.84335828) q[1];
sx q[1];
rz(-2.9620693) q[1];
rz(-1.8625453) q[3];
sx q[3];
rz(-2.1803133) q[3];
sx q[3];
rz(-2.9677344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(1.7857893) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(2.2858455) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6894585) q[0];
sx q[0];
rz(-1.2297213) q[0];
sx q[0];
rz(2.336691) q[0];
rz(-1.6542997) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(2.5022142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1460261) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(1.2703018) q[1];
rz(-pi) q[2];
rz(2.2780667) q[3];
sx q[3];
rz(-2.415495) q[3];
sx q[3];
rz(1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-1.1351599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.1816918) q[0];
sx q[0];
rz(-1.7745716) q[0];
x q[1];
rz(-1.8912002) q[2];
sx q[2];
rz(-1.5984629) q[2];
sx q[2];
rz(2.6507792) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4371722) q[1];
sx q[1];
rz(-1.51209) q[1];
sx q[1];
rz(-1.0874332) q[1];
rz(-pi) q[2];
rz(0.023852392) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6241374) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(1.5967782) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(-0.19616729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553186) q[0];
sx q[0];
rz(-2.1325169) q[0];
sx q[0];
rz(2.4713211) q[0];
rz(-pi) q[1];
rz(1.4354544) q[2];
sx q[2];
rz(-1.6207098) q[2];
sx q[2];
rz(2.4382255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6876467) q[1];
sx q[1];
rz(-1.5861662) q[1];
sx q[1];
rz(-0.83666283) q[1];
rz(-pi) q[2];
rz(2.2833061) q[3];
sx q[3];
rz(-2.2755816) q[3];
sx q[3];
rz(1.5376877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3907884) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-3.0604559) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(1.5224456) q[2];
sx q[2];
rz(-1.804525) q[2];
sx q[2];
rz(-1.3934025) q[2];
rz(3.0951981) q[3];
sx q[3];
rz(-1.3413324) q[3];
sx q[3];
rz(-3.0975292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
