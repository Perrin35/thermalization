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
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(1.7641915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3453774) q[0];
sx q[0];
rz(-2.2826129) q[0];
sx q[0];
rz(0.9839566) q[0];
rz(-pi) q[1];
rz(-0.032261511) q[2];
sx q[2];
rz(-2.2164377) q[2];
sx q[2];
rz(-2.1394859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5038576) q[1];
sx q[1];
rz(-1.7162995) q[1];
sx q[1];
rz(-0.18944959) q[1];
x q[2];
rz(-1.5903487) q[3];
sx q[3];
rz(-0.76977714) q[3];
sx q[3];
rz(2.584892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2237902) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(1.5000783) q[2];
rz(-2.3598059) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016945275) q[0];
sx q[0];
rz(-0.77777672) q[0];
sx q[0];
rz(-2.1710904) q[0];
rz(-1.3264725) q[1];
sx q[1];
rz(-1.7992203) q[1];
sx q[1];
rz(-0.91189799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058033179) q[0];
sx q[0];
rz(-0.61249295) q[0];
sx q[0];
rz(2.8508194) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3746337) q[2];
sx q[2];
rz(-1.7387783) q[2];
sx q[2];
rz(-0.84174918) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8205631) q[1];
sx q[1];
rz(-0.87087518) q[1];
sx q[1];
rz(0.053348347) q[1];
rz(3.1260388) q[3];
sx q[3];
rz(-2.0256666) q[3];
sx q[3];
rz(1.3223079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2328925) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(0.98928893) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99783889) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(-1.0379399) q[0];
rz(0.85820091) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(0.16955489) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586048) q[0];
sx q[0];
rz(-0.14359328) q[0];
sx q[0];
rz(-1.5989701) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8190024) q[2];
sx q[2];
rz(-1.9513522) q[2];
sx q[2];
rz(-1.1104465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3225786) q[1];
sx q[1];
rz(-0.59016229) q[1];
sx q[1];
rz(-2.1165089) q[1];
rz(-1.2081469) q[3];
sx q[3];
rz(-1.4261822) q[3];
sx q[3];
rz(3.0415158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1797336) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(-0.40599424) q[2];
rz(-0.77754846) q[3];
sx q[3];
rz(-0.33027875) q[3];
sx q[3];
rz(-2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4985713) q[0];
sx q[0];
rz(-1.2914456) q[0];
sx q[0];
rz(-0.92798573) q[0];
rz(0.058874933) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(1.7108797) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1851981) q[0];
sx q[0];
rz(-1.972359) q[0];
sx q[0];
rz(-0.08423452) q[0];
rz(0.48845993) q[2];
sx q[2];
rz(-0.40626981) q[2];
sx q[2];
rz(-0.84727188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18558606) q[1];
sx q[1];
rz(-1.348146) q[1];
sx q[1];
rz(1.0658479) q[1];
rz(1.1673314) q[3];
sx q[3];
rz(-0.66934062) q[3];
sx q[3];
rz(-0.40237967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4579939) q[2];
sx q[2];
rz(-1.4983838) q[2];
sx q[2];
rz(1.1452453) q[2];
rz(2.2942719) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(-1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1640846) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(0.384828) q[0];
rz(-0.79967868) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(1.5934561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037813) q[0];
sx q[0];
rz(-1.8124969) q[0];
sx q[0];
rz(-1.8221309) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20643955) q[2];
sx q[2];
rz(-1.1818794) q[2];
sx q[2];
rz(-0.85064519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56866327) q[1];
sx q[1];
rz(-1.8844114) q[1];
sx q[1];
rz(-0.36797087) q[1];
x q[2];
rz(1.0043275) q[3];
sx q[3];
rz(-0.62823717) q[3];
sx q[3];
rz(-2.7816176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3758214) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(0.080246298) q[2];
rz(0.56096983) q[3];
sx q[3];
rz(-1.6601945) q[3];
sx q[3];
rz(-2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4866667) q[0];
sx q[0];
rz(-2.4764562) q[0];
sx q[0];
rz(1.8810062) q[0];
rz(-0.27443019) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(2.4226277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089217606) q[0];
sx q[0];
rz(-1.9196654) q[0];
sx q[0];
rz(-2.265076) q[0];
rz(-pi) q[1];
rz(-3.1317091) q[2];
sx q[2];
rz(-2.3172969) q[2];
sx q[2];
rz(0.020992779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8985892) q[1];
sx q[1];
rz(-0.92783992) q[1];
sx q[1];
rz(-1.7220294) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9419233) q[3];
sx q[3];
rz(-1.7748482) q[3];
sx q[3];
rz(2.1851817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4869953) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(-0.60301644) q[2];
rz(1.4740137) q[3];
sx q[3];
rz(-2.1173729) q[3];
sx q[3];
rz(2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.5176158) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(2.1693443) q[1];
sx q[1];
rz(-2.9863803) q[1];
sx q[1];
rz(-0.42047277) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12026006) q[0];
sx q[0];
rz(-1.9870523) q[0];
sx q[0];
rz(0.72159213) q[0];
rz(-pi) q[1];
rz(1.0275317) q[2];
sx q[2];
rz(-1.2422556) q[2];
sx q[2];
rz(1.3628886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17294614) q[1];
sx q[1];
rz(-1.9526281) q[1];
sx q[1];
rz(-2.8416425) q[1];
rz(-pi) q[2];
rz(-1.917508) q[3];
sx q[3];
rz(-0.60255749) q[3];
sx q[3];
rz(1.6438705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.095852701) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(0.25071684) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(2.5417476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8641758) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(-0.94171062) q[0];
rz(3.1031389) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(3.0252735) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424902) q[0];
sx q[0];
rz(-1.1277871) q[0];
sx q[0];
rz(2.3194314) q[0];
rz(1.3384678) q[2];
sx q[2];
rz(-0.69455494) q[2];
sx q[2];
rz(1.4811999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1585576) q[1];
sx q[1];
rz(-1.3135034) q[1];
sx q[1];
rz(1.757181) q[1];
rz(2.8781901) q[3];
sx q[3];
rz(-2.7399094) q[3];
sx q[3];
rz(0.64611891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2945127) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(-1.026356) q[2];
rz(-1.9319084) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(0.9616372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-0.66861361) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(2.7630254) q[0];
rz(2.0319132) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(0.027677061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33998128) q[0];
sx q[0];
rz(-1.2072354) q[0];
sx q[0];
rz(-2.5021863) q[0];
rz(0.86556025) q[2];
sx q[2];
rz(-2.8606374) q[2];
sx q[2];
rz(-2.11175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3301446) q[1];
sx q[1];
rz(-1.7551577) q[1];
sx q[1];
rz(1.9346647) q[1];
x q[2];
rz(0.31881551) q[3];
sx q[3];
rz(-1.1622815) q[3];
sx q[3];
rz(-0.2948979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3045584) q[2];
sx q[2];
rz(-0.95418945) q[2];
sx q[2];
rz(2.7216116) q[2];
rz(0.84152451) q[3];
sx q[3];
rz(-1.0171112) q[3];
sx q[3];
rz(0.37874547) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(2.8564659) q[0];
rz(2.9410703) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(-2.3060422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052976089) q[0];
sx q[0];
rz(-1.1167913) q[0];
sx q[0];
rz(-2.3331621) q[0];
x q[1];
rz(1.5642688) q[2];
sx q[2];
rz(-2.2178429) q[2];
sx q[2];
rz(-1.4359635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0613868) q[1];
sx q[1];
rz(-0.57139423) q[1];
sx q[1];
rz(2.4953043) q[1];
rz(1.5316758) q[3];
sx q[3];
rz(-1.8062544) q[3];
sx q[3];
rz(2.2155516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42651132) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(2.4439028) q[2];
rz(-0.52608025) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(1.6349767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195011) q[0];
sx q[0];
rz(-0.87845907) q[0];
sx q[0];
rz(-0.39977951) q[0];
rz(-1.9922235) q[1];
sx q[1];
rz(-2.414357) q[1];
sx q[1];
rz(2.3946708) q[1];
rz(2.2483038) q[2];
sx q[2];
rz(-1.8468241) q[2];
sx q[2];
rz(2.1459864) q[2];
rz(-1.8676733) q[3];
sx q[3];
rz(-2.0991201) q[3];
sx q[3];
rz(0.42602628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
