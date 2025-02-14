OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.77312624) q[0];
sx q[0];
rz(-0.74690312) q[0];
sx q[0];
rz(-2.3009543) q[0];
rz(0.12159881) q[1];
sx q[1];
rz(-1.2727979) q[1];
sx q[1];
rz(0.29903856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75343412) q[0];
sx q[0];
rz(-2.3155766) q[0];
sx q[0];
rz(0.97919925) q[0];
x q[1];
rz(-0.077871696) q[2];
sx q[2];
rz(-1.0857333) q[2];
sx q[2];
rz(0.21202206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.92235293) q[1];
sx q[1];
rz(-1.427622) q[1];
sx q[1];
rz(-0.24210614) q[1];
x q[2];
rz(-0.14278966) q[3];
sx q[3];
rz(-2.0215109) q[3];
sx q[3];
rz(-2.5888909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3806939) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(-2.8139581) q[2];
rz(-1.7662175) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(-0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(2.2862527) q[0];
rz(-3.1335462) q[1];
sx q[1];
rz(-1.8549553) q[1];
sx q[1];
rz(-0.45113742) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83734918) q[0];
sx q[0];
rz(-0.74241246) q[0];
sx q[0];
rz(-1.3377331) q[0];
x q[1];
rz(0.19719736) q[2];
sx q[2];
rz(-0.78193808) q[2];
sx q[2];
rz(2.8213843) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0135368) q[1];
sx q[1];
rz(-1.3971796) q[1];
sx q[1];
rz(-0.44349576) q[1];
x q[2];
rz(1.9782009) q[3];
sx q[3];
rz(-1.1278858) q[3];
sx q[3];
rz(-0.16678424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96192876) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(-0.12953225) q[2];
rz(-0.1772964) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(-3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.391908) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(2.5823197) q[0];
rz(-3.0468805) q[1];
sx q[1];
rz(-1.6030703) q[1];
sx q[1];
rz(2.6643378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1806948) q[0];
sx q[0];
rz(-1.8261693) q[0];
sx q[0];
rz(-1.9810505) q[0];
rz(-pi) q[1];
rz(-2.7626286) q[2];
sx q[2];
rz(-0.38658374) q[2];
sx q[2];
rz(-1.6340337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4871769) q[1];
sx q[1];
rz(-2.3337939) q[1];
sx q[1];
rz(0.31262763) q[1];
x q[2];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.1616544) q[3];
sx q[3];
rz(-2.9144998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7462848) q[2];
sx q[2];
rz(-1.8751273) q[2];
sx q[2];
rz(2.2115808) q[2];
rz(1.2578472) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(3.0959082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3607218) q[0];
sx q[0];
rz(-1.5732795) q[0];
sx q[0];
rz(2.1029396) q[0];
rz(0.61141283) q[1];
sx q[1];
rz(-2.2544506) q[1];
sx q[1];
rz(-2.2183653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1425689) q[0];
sx q[0];
rz(-1.2272738) q[0];
sx q[0];
rz(0.96238636) q[0];
rz(-0.011767894) q[2];
sx q[2];
rz(-1.1530877) q[2];
sx q[2];
rz(-2.2209446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0765517) q[1];
sx q[1];
rz(-0.84969798) q[1];
sx q[1];
rz(-0.11252071) q[1];
x q[2];
rz(-2.491288) q[3];
sx q[3];
rz(-1.7803444) q[3];
sx q[3];
rz(2.9858231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5103147) q[2];
sx q[2];
rz(-2.6738561) q[2];
sx q[2];
rz(2.5941217) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(0.87375435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9861458) q[0];
sx q[0];
rz(-2.2387945) q[0];
sx q[0];
rz(1.4601532) q[0];
rz(1.0001146) q[1];
sx q[1];
rz(-2.2747048) q[1];
sx q[1];
rz(0.11988457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.125357) q[0];
sx q[0];
rz(-0.88847697) q[0];
sx q[0];
rz(2.8518139) q[0];
rz(-2.5566543) q[2];
sx q[2];
rz(-2.1368933) q[2];
sx q[2];
rz(1.0370129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5226871) q[1];
sx q[1];
rz(-2.9123016) q[1];
sx q[1];
rz(-0.35405901) q[1];
rz(-pi) q[2];
rz(-2.5402903) q[3];
sx q[3];
rz(-1.9996694) q[3];
sx q[3];
rz(1.2473904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1025866) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(-0.26068035) q[2];
rz(2.2634704) q[3];
sx q[3];
rz(-1.2295281) q[3];
sx q[3];
rz(2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6108625) q[0];
sx q[0];
rz(-0.090228883) q[0];
sx q[0];
rz(2.1395785) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-0.95163029) q[1];
sx q[1];
rz(2.196905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2827891) q[0];
sx q[0];
rz(-1.7886046) q[0];
sx q[0];
rz(2.5522638) q[0];
rz(-2.281419) q[2];
sx q[2];
rz(-2.3322736) q[2];
sx q[2];
rz(-3.0391673) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8807043) q[1];
sx q[1];
rz(-0.35683888) q[1];
sx q[1];
rz(2.6964705) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3768436) q[3];
sx q[3];
rz(-0.55046457) q[3];
sx q[3];
rz(0.46904072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3433596) q[2];
sx q[2];
rz(-1.6075906) q[2];
sx q[2];
rz(-3.0976683) q[2];
rz(-1.515306) q[3];
sx q[3];
rz(-2.01912) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2783022) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(-2.5569051) q[0];
rz(-1.0514642) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(2.6712766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0338351) q[0];
sx q[0];
rz(-1.4246539) q[0];
sx q[0];
rz(-0.58476292) q[0];
rz(-2.7696848) q[2];
sx q[2];
rz(-1.1924517) q[2];
sx q[2];
rz(-1.3828948) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85519771) q[1];
sx q[1];
rz(-2.2620387) q[1];
sx q[1];
rz(-1.0950086) q[1];
rz(-pi) q[2];
rz(3.028585) q[3];
sx q[3];
rz(-1.5829493) q[3];
sx q[3];
rz(0.069610217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8768002) q[2];
sx q[2];
rz(-2.6046643) q[2];
sx q[2];
rz(0.31141591) q[2];
rz(-0.16820678) q[3];
sx q[3];
rz(-1.5752537) q[3];
sx q[3];
rz(-2.8581207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51157057) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-0.93609634) q[0];
rz(2.6452737) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(2.671303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449438) q[0];
sx q[0];
rz(-1.2723288) q[0];
sx q[0];
rz(-2.6281283) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3730714) q[2];
sx q[2];
rz(-1.5682967) q[2];
sx q[2];
rz(3.1001774) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6185604) q[1];
sx q[1];
rz(-1.9123239) q[1];
sx q[1];
rz(-0.073445436) q[1];
rz(-pi) q[2];
rz(-1.3948729) q[3];
sx q[3];
rz(-1.298549) q[3];
sx q[3];
rz(-1.3678838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63142598) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(1.6660956) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(-1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.0610166) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(0.06037816) q[0];
rz(2.9810442) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(-2.1626332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87099052) q[0];
sx q[0];
rz(-2.5033062) q[0];
sx q[0];
rz(2.9459475) q[0];
rz(-pi) q[1];
rz(2.253112) q[2];
sx q[2];
rz(-1.3559301) q[2];
sx q[2];
rz(-0.14718283) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0461523) q[1];
sx q[1];
rz(-0.80712026) q[1];
sx q[1];
rz(0.67428204) q[1];
rz(-pi) q[2];
rz(0.99154226) q[3];
sx q[3];
rz(-1.5421621) q[3];
sx q[3];
rz(-0.34285173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0742566) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(0.072889797) q[2];
rz(-2.5439751) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(-0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.446796) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(-2.6126675) q[0];
rz(-0.18813285) q[1];
sx q[1];
rz(-0.70536047) q[1];
sx q[1];
rz(2.5573152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0463294) q[0];
sx q[0];
rz(-1.8285311) q[0];
sx q[0];
rz(-1.3168174) q[0];
rz(-1.1816013) q[2];
sx q[2];
rz(-1.9799332) q[2];
sx q[2];
rz(-1.429806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6021298) q[1];
sx q[1];
rz(-1.5279084) q[1];
sx q[1];
rz(-0.40956386) q[1];
rz(-2.1891865) q[3];
sx q[3];
rz(-1.3460396) q[3];
sx q[3];
rz(-2.145203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3820485) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(-3.0697401) q[2];
rz(-1.0233277) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(-2.6527827) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95579424) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
rz(-2.4664948) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(2.6266392) q[2];
sx q[2];
rz(-1.6390159) q[2];
sx q[2];
rz(-1.9255571) q[2];
rz(1.731338) q[3];
sx q[3];
rz(-1.4057126) q[3];
sx q[3];
rz(-2.2688903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
