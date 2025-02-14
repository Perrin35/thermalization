OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12914817) q[0];
sx q[0];
rz(-1.669786) q[0];
sx q[0];
rz(0.9653402) q[0];
rz(-3.1211634) q[1];
sx q[1];
rz(4.6252146) q[1];
sx q[1];
rz(10.802179) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3453774) q[0];
sx q[0];
rz(-0.85897972) q[0];
sx q[0];
rz(-2.1576361) q[0];
rz(0.92490478) q[2];
sx q[2];
rz(-1.5450302) q[2];
sx q[2];
rz(2.5534867) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6377351) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(-0.18944959) q[1];
rz(-pi) q[2];
rz(1.5512439) q[3];
sx q[3];
rz(-2.3718155) q[3];
sx q[3];
rz(-2.584892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2237902) q[2];
sx q[2];
rz(-1.6518355) q[2];
sx q[2];
rz(1.6415143) q[2];
rz(-0.7817868) q[3];
sx q[3];
rz(-0.89038554) q[3];
sx q[3];
rz(2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1246474) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(-2.1710904) q[0];
rz(-1.8151201) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(-0.91189799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8489979) q[0];
sx q[0];
rz(-2.1540897) q[0];
sx q[0];
rz(-1.769577) q[0];
rz(-pi) q[1];
rz(-0.85477306) q[2];
sx q[2];
rz(-2.8840384) q[2];
sx q[2];
rz(-1.4285029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8205631) q[1];
sx q[1];
rz(-0.87087518) q[1];
sx q[1];
rz(0.053348347) q[1];
x q[2];
rz(0.015553899) q[3];
sx q[3];
rz(-1.115926) q[3];
sx q[3];
rz(-1.8192847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(-2.1523037) q[2];
rz(2.7214637) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(-0.72107983) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1437538) q[0];
sx q[0];
rz(-1.1622575) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(2.2833917) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(-0.16955489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586048) q[0];
sx q[0];
rz(-0.14359328) q[0];
sx q[0];
rz(-1.5989701) q[0];
rz(-pi) q[1];
rz(2.2409147) q[2];
sx q[2];
rz(-2.6478516) q[2];
sx q[2];
rz(-2.7639219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9537064) q[1];
sx q[1];
rz(-2.0665209) q[1];
sx q[1];
rz(0.33457054) q[1];
x q[2];
rz(-1.1812594) q[3];
sx q[3];
rz(-2.7523605) q[3];
sx q[3];
rz(-1.8336982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9618591) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(0.40599424) q[2];
rz(0.77754846) q[3];
sx q[3];
rz(-0.33027875) q[3];
sx q[3];
rz(2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4985713) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(2.2136069) q[0];
rz(3.0827177) q[1];
sx q[1];
rz(-1.4480271) q[1];
sx q[1];
rz(1.7108797) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1851981) q[0];
sx q[0];
rz(-1.972359) q[0];
sx q[0];
rz(-0.08423452) q[0];
x q[1];
rz(2.6531327) q[2];
sx q[2];
rz(-0.40626981) q[2];
sx q[2];
rz(0.84727188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9560066) q[1];
sx q[1];
rz(-1.348146) q[1];
sx q[1];
rz(-2.0757448) q[1];
rz(2.1998422) q[3];
sx q[3];
rz(-1.3247196) q[3];
sx q[3];
rz(1.6501282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4579939) q[2];
sx q[2];
rz(-1.4983838) q[2];
sx q[2];
rz(1.1452453) q[2];
rz(0.84732071) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(-1.639521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1640846) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(0.384828) q[0];
rz(-2.341914) q[1];
sx q[1];
rz(-2.622602) q[1];
sx q[1];
rz(-1.5481366) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6131957) q[0];
sx q[0];
rz(-1.8146744) q[0];
sx q[0];
rz(-2.8923678) q[0];
rz(1.9673011) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(-0.79939524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3272737) q[1];
sx q[1];
rz(-2.6627935) q[1];
sx q[1];
rz(-0.73368608) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1206854) q[3];
sx q[3];
rz(-1.249915) q[3];
sx q[3];
rz(0.73559092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7657713) q[2];
sx q[2];
rz(-0.69362005) q[2];
sx q[2];
rz(-0.080246298) q[2];
rz(2.5806228) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(-2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4866667) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(-1.2605865) q[0];
rz(2.8671625) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(2.4226277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.270416) q[0];
sx q[0];
rz(-0.76379062) q[0];
sx q[0];
rz(2.0877354) q[0];
rz(-pi) q[1];
rz(0.0098835398) q[2];
sx q[2];
rz(-0.82429574) q[2];
sx q[2];
rz(-0.020992779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9940954) q[1];
sx q[1];
rz(-2.4835517) q[1];
sx q[1];
rz(0.19845732) q[1];
x q[2];
rz(-1.9419233) q[3];
sx q[3];
rz(-1.3667445) q[3];
sx q[3];
rz(-0.95641092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6545973) q[2];
sx q[2];
rz(-1.2664814) q[2];
sx q[2];
rz(-0.60301644) q[2];
rz(-1.6675789) q[3];
sx q[3];
rz(-2.1173729) q[3];
sx q[3];
rz(-0.41072861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5176158) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(-0.97224832) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(0.42047277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1218106) q[0];
sx q[0];
rz(-0.81392899) q[0];
sx q[0];
rz(-2.5518083) q[0];
rz(-pi) q[1];
rz(-2.7626129) q[2];
sx q[2];
rz(-2.0820658) q[2];
sx q[2];
rz(-3.1261217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6290038) q[1];
sx q[1];
rz(-1.2930451) q[1];
sx q[1];
rz(1.968683) q[1];
rz(-pi) q[2];
rz(-2.1450348) q[3];
sx q[3];
rz(-1.3769994) q[3];
sx q[3];
rz(0.36239788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.095852701) q[2];
sx q[2];
rz(-2.0407712) q[2];
sx q[2];
rz(-2.8908758) q[2];
rz(1.8935253) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(2.5417476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-2.5499948) q[0];
sx q[0];
rz(2.199882) q[0];
rz(-3.1031389) q[1];
sx q[1];
rz(-1.7105303) q[1];
sx q[1];
rz(3.0252735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424902) q[0];
sx q[0];
rz(-1.1277871) q[0];
sx q[0];
rz(2.3194314) q[0];
rz(0.88960008) q[2];
sx q[2];
rz(-1.4228914) q[2];
sx q[2];
rz(0.090229457) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6623936) q[1];
sx q[1];
rz(-2.8251007) q[1];
sx q[1];
rz(0.6135769) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8781901) q[3];
sx q[3];
rz(-0.40168328) q[3];
sx q[3];
rz(0.64611891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8470799) q[2];
sx q[2];
rz(-0.81081644) q[2];
sx q[2];
rz(2.1152367) q[2];
rz(1.2096842) q[3];
sx q[3];
rz(-2.1116202) q[3];
sx q[3];
rz(-0.9616372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.66861361) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(-2.7630254) q[0];
rz(-1.1096795) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(-3.1139156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4647163) q[0];
sx q[0];
rz(-2.4188571) q[0];
sx q[0];
rz(0.56761543) q[0];
x q[1];
rz(-0.1849298) q[2];
sx q[2];
rz(-1.3580631) q[2];
sx q[2];
rz(-2.8367868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82904774) q[1];
sx q[1];
rz(-1.9282189) q[1];
sx q[1];
rz(-2.9446359) q[1];
rz(2.8227771) q[3];
sx q[3];
rz(-1.1622815) q[3];
sx q[3];
rz(0.2948979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3045584) q[2];
sx q[2];
rz(-0.95418945) q[2];
sx q[2];
rz(2.7216116) q[2];
rz(-0.84152451) q[3];
sx q[3];
rz(-1.0171112) q[3];
sx q[3];
rz(2.7628472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76641744) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(-0.28512678) q[0];
rz(2.9410703) q[1];
sx q[1];
rz(-1.2031809) q[1];
sx q[1];
rz(2.3060422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052976089) q[0];
sx q[0];
rz(-1.1167913) q[0];
sx q[0];
rz(-0.80843057) q[0];
rz(-pi) q[1];
rz(0.64705683) q[2];
sx q[2];
rz(-1.5655883) q[2];
sx q[2];
rz(-0.13089779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.056014765) q[1];
sx q[1];
rz(-1.2390572) q[1];
sx q[1];
rz(0.47420926) q[1];
rz(-pi) q[2];
rz(1.5316758) q[3];
sx q[3];
rz(-1.8062544) q[3];
sx q[3];
rz(-0.92604107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7150813) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(0.69768989) q[2];
rz(2.6155124) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(-1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220916) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(1.1493692) q[1];
sx q[1];
rz(-2.414357) q[1];
sx q[1];
rz(2.3946708) q[1];
rz(1.9952075) q[2];
sx q[2];
rz(-2.41832) q[2];
sx q[2];
rz(0.24857749) q[2];
rz(2.5935843) q[3];
sx q[3];
rz(-1.3153793) q[3];
sx q[3];
rz(1.8438189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
