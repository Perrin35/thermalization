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
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(2.5002313) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8279995) q[0];
sx q[0];
rz(-1.8090973) q[0];
sx q[0];
rz(-0.39512623) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85513656) q[2];
sx q[2];
rz(-0.54358608) q[2];
sx q[2];
rz(0.91993514) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8909059) q[1];
sx q[1];
rz(-0.85039925) q[1];
sx q[1];
rz(-0.62178639) q[1];
x q[2];
rz(-2.5121243) q[3];
sx q[3];
rz(-1.2347722) q[3];
sx q[3];
rz(-2.0365086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.550094) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(2.6634898) q[2];
rz(-1.6889307) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(1.0189198) q[0];
rz(1.4787176) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(-0.63308024) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3490152) q[0];
sx q[0];
rz(-1.4903869) q[0];
sx q[0];
rz(-2.3809459) q[0];
rz(-pi) q[1];
rz(-2.1585629) q[2];
sx q[2];
rz(-0.9126185) q[2];
sx q[2];
rz(-0.36511974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9004904) q[1];
sx q[1];
rz(-1.2304658) q[1];
sx q[1];
rz(-1.1452922) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31140621) q[3];
sx q[3];
rz(-2.2283471) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35976609) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(-3.0241372) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85787073) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(3.0531847) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-0.69082469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4118977) q[0];
sx q[0];
rz(-1.6459961) q[0];
sx q[0];
rz(-2.9883283) q[0];
rz(-pi) q[1];
rz(-0.76491852) q[2];
sx q[2];
rz(-1.7419445) q[2];
sx q[2];
rz(0.82676065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9835811) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(1.4278825) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3986474) q[3];
sx q[3];
rz(-2.6446614) q[3];
sx q[3];
rz(-3.1019773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92418015) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(3.0351191) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(-2.6191214) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(-2.5879588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.722256) q[0];
sx q[0];
rz(-2.0291078) q[0];
sx q[0];
rz(-0.76461794) q[0];
rz(-pi) q[1];
rz(2.3217818) q[2];
sx q[2];
rz(-2.1852583) q[2];
sx q[2];
rz(0.43624207) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5987293) q[1];
sx q[1];
rz(-2.3485564) q[1];
sx q[1];
rz(-1.7584156) q[1];
x q[2];
rz(3.1133075) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(-0.68238168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(2.6468357) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(-1.2899227) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49098) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(-2.0565128) q[0];
rz(-0.96013534) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(0.18403149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40084307) q[0];
sx q[0];
rz(-1.6753734) q[0];
sx q[0];
rz(-2.8507289) q[0];
rz(-pi) q[1];
rz(-1.2940302) q[2];
sx q[2];
rz(-1.9230611) q[2];
sx q[2];
rz(-2.2500452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6726471) q[1];
sx q[1];
rz(-0.99860672) q[1];
sx q[1];
rz(1.0020301) q[1];
x q[2];
rz(2.9156978) q[3];
sx q[3];
rz(-0.73879209) q[3];
sx q[3];
rz(1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(-2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891156) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(2.9011762) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(0.14850798) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69420544) q[0];
sx q[0];
rz(-0.84852695) q[0];
sx q[0];
rz(1.3556051) q[0];
rz(-pi) q[1];
rz(2.912942) q[2];
sx q[2];
rz(-1.6778523) q[2];
sx q[2];
rz(-2.6766237) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0420694) q[1];
sx q[1];
rz(-1.4617141) q[1];
sx q[1];
rz(1.3101577) q[1];
rz(-pi) q[2];
rz(1.1735736) q[3];
sx q[3];
rz(-2.0052611) q[3];
sx q[3];
rz(1.43515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(0.48940247) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(-1.874812) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4483036) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(-2.3235902) q[0];
rz(-1.2524293) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-2.0163527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97911191) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(0.11066779) q[0];
x q[1];
rz(-1.869404) q[2];
sx q[2];
rz(-1.4223137) q[2];
sx q[2];
rz(0.4543002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69562558) q[1];
sx q[1];
rz(-1.274316) q[1];
sx q[1];
rz(-2.8545024) q[1];
x q[2];
rz(-2.4756487) q[3];
sx q[3];
rz(-1.1083372) q[3];
sx q[3];
rz(2.8921814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0685048) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97061625) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.6059426) q[0];
rz(-1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(2.1309526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686417) q[0];
sx q[0];
rz(-2.3803582) q[0];
sx q[0];
rz(-1.3377454) q[0];
rz(2.0660731) q[2];
sx q[2];
rz(-1.8349378) q[2];
sx q[2];
rz(1.8758945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6992221) q[1];
sx q[1];
rz(-1.3500299) q[1];
sx q[1];
rz(-1.3352331) q[1];
x q[2];
rz(-0.2145433) q[3];
sx q[3];
rz(-0.40896591) q[3];
sx q[3];
rz(-1.4734801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(0.69407216) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086031832) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(-0.075597413) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6583017) q[0];
sx q[0];
rz(-2.1245983) q[0];
sx q[0];
rz(2.7730586) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7463023) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(-0.67948558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52121431) q[1];
sx q[1];
rz(-1.4676536) q[1];
sx q[1];
rz(2.4528273) q[1];
rz(-pi) q[2];
rz(-0.62065403) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(1.8295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8081234) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(1.8010275) q[2];
rz(-2.8373485) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-2.9647968) q[0];
rz(1.2416174) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(2.7005844) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2947185) q[0];
sx q[0];
rz(-0.28781578) q[0];
sx q[0];
rz(-2.3527282) q[0];
rz(-0.44956019) q[2];
sx q[2];
rz(-2.9712147) q[2];
sx q[2];
rz(-2.5296488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8969352) q[1];
sx q[1];
rz(-1.1420297) q[1];
sx q[1];
rz(-0.37533092) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2714083) q[3];
sx q[3];
rz(-2.1736439) q[3];
sx q[3];
rz(2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5796154) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(2.9368403) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4176035) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(-3.1148615) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(1.0529636) q[2];
sx q[2];
rz(-2.3452407) q[2];
sx q[2];
rz(-1.8635545) q[2];
rz(0.7209575) q[3];
sx q[3];
rz(-2.4506035) q[3];
sx q[3];
rz(-2.639365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
