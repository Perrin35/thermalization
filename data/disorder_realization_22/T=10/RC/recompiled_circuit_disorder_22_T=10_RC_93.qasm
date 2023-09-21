OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(4.099457) q[0];
sx q[0];
rz(9.2803331) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(2.6161939) q[1];
sx q[1];
rz(8.4470196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.190783) q[0];
sx q[0];
rz(-1.8291744) q[0];
sx q[0];
rz(1.4825312) q[0];
rz(-pi) q[1];
rz(-0.30774967) q[2];
sx q[2];
rz(-2.0648742) q[2];
sx q[2];
rz(-2.5623164) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.830476) q[1];
sx q[1];
rz(-1.4314572) q[1];
sx q[1];
rz(1.4275101) q[1];
rz(1.3917189) q[3];
sx q[3];
rz(-2.350051) q[3];
sx q[3];
rz(2.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67291659) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(0.19876924) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(-2.7804651) q[0];
rz(-1.7065642) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(0.8180058) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5411103) q[0];
sx q[0];
rz(-0.4821018) q[0];
sx q[0];
rz(-2.5151398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9423219) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(-0.83100806) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6815731) q[1];
sx q[1];
rz(-0.96833723) q[1];
sx q[1];
rz(0.54089344) q[1];
x q[2];
rz(-1.9345476) q[3];
sx q[3];
rz(-1.517429) q[3];
sx q[3];
rz(-1.3747017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(1.4206295) q[2];
rz(1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.7398359) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(-0.2972163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5706144) q[0];
sx q[0];
rz(-1.1665205) q[0];
sx q[0];
rz(1.6398318) q[0];
rz(-pi) q[1];
rz(2.3702413) q[2];
sx q[2];
rz(-2.1043679) q[2];
sx q[2];
rz(-2.9750864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3618468) q[1];
sx q[1];
rz(-1.6184813) q[1];
sx q[1];
rz(1.3839128) q[1];
rz(-2.3171595) q[3];
sx q[3];
rz(-1.0472877) q[3];
sx q[3];
rz(-1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(1.9528495) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(-0.52743131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(-0.71290839) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(0.68177044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503158) q[0];
sx q[0];
rz(-2.804545) q[0];
sx q[0];
rz(2.2469673) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8800456) q[2];
sx q[2];
rz(-2.5144858) q[2];
sx q[2];
rz(-0.68324616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3817953) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(2.1684907) q[1];
x q[2];
rz(0.81687974) q[3];
sx q[3];
rz(-0.48135346) q[3];
sx q[3];
rz(-2.4933185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-0.95820367) q[2];
rz(-1.0664252) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(-2.5812896) q[0];
rz(-0.99984461) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.6220185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1006267) q[0];
sx q[0];
rz(-2.0140411) q[0];
sx q[0];
rz(-2.6020781) q[0];
rz(-pi) q[1];
rz(-1.0137453) q[2];
sx q[2];
rz(-2.1079014) q[2];
sx q[2];
rz(0.86442664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2625761) q[1];
sx q[1];
rz(-0.82619709) q[1];
sx q[1];
rz(-1.8415585) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12445478) q[3];
sx q[3];
rz(-2.180047) q[3];
sx q[3];
rz(0.54427687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(-3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(-0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.8054304) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(1.3609715) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(1.8575352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43858466) q[0];
sx q[0];
rz(-0.87235886) q[0];
sx q[0];
rz(-0.75011487) q[0];
rz(-pi) q[1];
rz(-2.1806296) q[2];
sx q[2];
rz(-0.46513882) q[2];
sx q[2];
rz(-2.7679408) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18404993) q[1];
sx q[1];
rz(-1.8584538) q[1];
sx q[1];
rz(2.6581453) q[1];
x q[2];
rz(-2.0606023) q[3];
sx q[3];
rz(-2.5488857) q[3];
sx q[3];
rz(0.61471516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.56629431) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(-0.56754011) q[0];
rz(-0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(2.2033851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81323775) q[0];
sx q[0];
rz(-0.52044808) q[0];
sx q[0];
rz(-2.2195199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1392518) q[2];
sx q[2];
rz(-0.92407862) q[2];
sx q[2];
rz(-1.9288174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36233703) q[1];
sx q[1];
rz(-1.926683) q[1];
sx q[1];
rz(0.87902714) q[1];
rz(-2.7890117) q[3];
sx q[3];
rz(-0.73871021) q[3];
sx q[3];
rz(2.0032361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(1.3860469) q[2];
rz(-1.322768) q[3];
sx q[3];
rz(-1.991792) q[3];
sx q[3];
rz(-3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(1.4338795) q[0];
rz(1.8677615) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(-2.3103255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76038218) q[0];
sx q[0];
rz(-1.7179278) q[0];
sx q[0];
rz(-0.40572625) q[0];
rz(-2.2057461) q[2];
sx q[2];
rz(-1.4668462) q[2];
sx q[2];
rz(-0.39633358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23849328) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(-2.3563983) q[1];
rz(-pi) q[2];
rz(1.3475111) q[3];
sx q[3];
rz(-1.3459473) q[3];
sx q[3];
rz(-0.68883483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(-1.9780805) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7678541) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(2.1642165) q[0];
rz(-1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(-1.1057373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10406636) q[0];
sx q[0];
rz(-0.46188799) q[0];
sx q[0];
rz(2.877263) q[0];
x q[1];
rz(-2.2374723) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(-2.3711575) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2178206) q[1];
sx q[1];
rz(-0.2982699) q[1];
sx q[1];
rz(0.89426269) q[1];
rz(2.7902778) q[3];
sx q[3];
rz(-1.8422541) q[3];
sx q[3];
rz(0.16004496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5332807) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(2.9917955) q[2];
rz(-1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(-0.1272442) q[0];
rz(1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(2.9737934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6608749) q[0];
sx q[0];
rz(-1.6010451) q[0];
sx q[0];
rz(-1.589993) q[0];
x q[1];
rz(2.9211505) q[2];
sx q[2];
rz(-2.5219005) q[2];
sx q[2];
rz(-2.9881791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93712902) q[1];
sx q[1];
rz(-0.58848721) q[1];
sx q[1];
rz(-3.1292533) q[1];
x q[2];
rz(-0.12331788) q[3];
sx q[3];
rz(-1.6880369) q[3];
sx q[3];
rz(-2.5532212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(-2.3804469) q[2];
rz(-3.051565) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54820838) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(0.17100632) q[3];
sx q[3];
rz(-1.0382367) q[3];
sx q[3];
rz(2.3259179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];