OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(-1.5313907) q[1];
sx q[1];
rz(-2.664497) q[1];
sx q[1];
rz(2.6452126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8007322) q[0];
sx q[0];
rz(-2.1600318) q[0];
sx q[0];
rz(3.0032934) q[0];
rz(-pi) q[1];
rz(0.33978396) q[2];
sx q[2];
rz(-1.3801563) q[2];
sx q[2];
rz(-1.2059739) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9909518) q[1];
sx q[1];
rz(-0.71326637) q[1];
sx q[1];
rz(-0.13891061) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6707889) q[3];
sx q[3];
rz(-1.1649719) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-0.28764763) q[2];
rz(-1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(-2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47050414) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(-1.6862306) q[0];
x q[1];
rz(-2.4096476) q[2];
sx q[2];
rz(-2.7436896) q[2];
sx q[2];
rz(-2.0741472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4430122) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(1.7673311) q[1];
rz(-pi) q[2];
rz(2.3489477) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(1.45989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-0.18297718) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17830081) q[0];
sx q[0];
rz(-2.460647) q[0];
sx q[0];
rz(-2.2292024) q[0];
rz(2.6687713) q[2];
sx q[2];
rz(-0.46337767) q[2];
sx q[2];
rz(-2.9544427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80013393) q[1];
sx q[1];
rz(-2.027249) q[1];
sx q[1];
rz(-1.9564499) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1312917) q[3];
sx q[3];
rz(-2.8539538) q[3];
sx q[3];
rz(-1.079263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6782137) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(-2.4257207) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7870165) q[0];
sx q[0];
rz(-2.0873318) q[0];
sx q[0];
rz(1.8754965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9348346) q[2];
sx q[2];
rz(-2.0419428) q[2];
sx q[2];
rz(-0.83540321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.48651931) q[1];
sx q[1];
rz(-1.4078948) q[1];
sx q[1];
rz(1.4763114) q[1];
x q[2];
rz(-0.4977787) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(-1.2342412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(-1.9013566) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(2.48824) q[0];
rz(-pi) q[1];
rz(-0.93047662) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(-3.004068) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56375757) q[1];
sx q[1];
rz(-1.1760684) q[1];
sx q[1];
rz(1.5866336) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1122123) q[3];
sx q[3];
rz(-2.5684528) q[3];
sx q[3];
rz(1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(1.0995964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(-2.7307672) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1018283) q[2];
sx q[2];
rz(-1.7341988) q[2];
sx q[2];
rz(-2.642717) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6096566) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(2.460536) q[1];
rz(-pi) q[2];
rz(-0.78824708) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(-2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866078) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(0.64487547) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(0.15730102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042127) q[0];
sx q[0];
rz(-2.3010845) q[0];
sx q[0];
rz(-1.69676) q[0];
rz(1.332875) q[2];
sx q[2];
rz(-2.3908797) q[2];
sx q[2];
rz(-2.6964292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0485222) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(-2.1919769) q[1];
x q[2];
rz(2.7492206) q[3];
sx q[3];
rz(-2.4761768) q[3];
sx q[3];
rz(0.80442807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-2.3106993) q[2];
rz(3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1300668) q[0];
sx q[0];
rz(-1.5648601) q[0];
sx q[0];
rz(2.5991711) q[0];
x q[1];
rz(1.3515527) q[2];
sx q[2];
rz(-2.4862137) q[2];
sx q[2];
rz(0.13742451) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9353232) q[1];
sx q[1];
rz(-0.85201293) q[1];
sx q[1];
rz(-1.000559) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9803489) q[3];
sx q[3];
rz(-1.1405986) q[3];
sx q[3];
rz(-2.6524515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.04348065) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(0.44257277) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(1.7907422) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25778739) q[0];
sx q[0];
rz(-1.4865849) q[0];
sx q[0];
rz(3.0620831) q[0];
x q[1];
rz(-3.0355989) q[2];
sx q[2];
rz(-0.66291891) q[2];
sx q[2];
rz(-1.7352599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6752154) q[1];
sx q[1];
rz(-1.5545168) q[1];
sx q[1];
rz(1.1364163) q[1];
x q[2];
rz(-2.001858) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(2.7800625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8266325) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(2.3610624) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(2.9719877) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7923332) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(-0.40339289) q[0];
x q[1];
rz(1.595396) q[2];
sx q[2];
rz(-1.149017) q[2];
sx q[2];
rz(1.8116902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.604653) q[1];
sx q[1];
rz(-0.61652684) q[1];
sx q[1];
rz(-1.489868) q[1];
rz(-pi) q[2];
rz(1.1630837) q[3];
sx q[3];
rz(-0.34808961) q[3];
sx q[3];
rz(-2.7365494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85049373) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8511843) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(2.6295173) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-1.2693263) q[2];
sx q[2];
rz(-1.3133247) q[2];
sx q[2];
rz(1.7016344) q[2];
rz(0.11733304) q[3];
sx q[3];
rz(-1.3640079) q[3];
sx q[3];
rz(-0.97466536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
