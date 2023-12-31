OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(1.024363) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60351935) q[0];
sx q[0];
rz(-0.78342122) q[0];
sx q[0];
rz(2.6253683) q[0];
rz(-2.4835303) q[2];
sx q[2];
rz(-1.0339289) q[2];
sx q[2];
rz(1.3047578) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3633903) q[1];
sx q[1];
rz(-1.2092023) q[1];
sx q[1];
rz(1.1671288) q[1];
rz(-pi) q[2];
rz(0.23006769) q[3];
sx q[3];
rz(-2.821273) q[3];
sx q[3];
rz(-0.39428082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(-2.8406075) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99012016) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(-2.5734076) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53595397) q[2];
sx q[2];
rz(-1.1079271) q[2];
sx q[2];
rz(0.11975372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4251551) q[1];
sx q[1];
rz(-1.5967224) q[1];
sx q[1];
rz(1.9772711) q[1];
x q[2];
rz(-0.95275767) q[3];
sx q[3];
rz(-0.85988322) q[3];
sx q[3];
rz(0.074912138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.6111074) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(-2.5684165) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983551) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(1.7217365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77261749) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(2.2881743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0108311) q[1];
sx q[1];
rz(-1.8567994) q[1];
sx q[1];
rz(-0.39949135) q[1];
rz(1.1099986) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(-2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(2.7788924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6837316) q[0];
sx q[0];
rz(-2.3944003) q[0];
sx q[0];
rz(-1.213221) q[0];
rz(1.6022801) q[2];
sx q[2];
rz(-1.5891979) q[2];
sx q[2];
rz(3.1061663) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73211654) q[1];
sx q[1];
rz(-2.5001657) q[1];
sx q[1];
rz(-0.17318053) q[1];
rz(-pi) q[2];
x q[2];
rz(0.072399541) q[3];
sx q[3];
rz(-0.94745938) q[3];
sx q[3];
rz(0.49923957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(3.1385699) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(0.80250424) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97809726) q[0];
sx q[0];
rz(-0.28143829) q[0];
sx q[0];
rz(1.0174169) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2541788) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(2.8487157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38813218) q[1];
sx q[1];
rz(-1.7963444) q[1];
sx q[1];
rz(-1.3688341) q[1];
x q[2];
rz(-1.0110537) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(1.884348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(-0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(0.13866436) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4045227) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(-1.5674595) q[0];
rz(-pi) q[1];
rz(2.7799941) q[2];
sx q[2];
rz(-2.9271759) q[2];
sx q[2];
rz(-0.32985652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0181959) q[1];
sx q[1];
rz(-2.1235848) q[1];
sx q[1];
rz(-2.4559896) q[1];
x q[2];
rz(1.9964553) q[3];
sx q[3];
rz(-0.87493757) q[3];
sx q[3];
rz(-1.3130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(-1.8072051) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(0.62430635) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3793959) q[0];
sx q[0];
rz(-1.7805829) q[0];
sx q[0];
rz(0.72797738) q[0];
rz(0.95958556) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(-0.50146539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.034060409) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(-0.88453102) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94654406) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(2.7590511) q[2];
rz(3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-3.0126742) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0607973) q[0];
sx q[0];
rz(-1.1757438) q[0];
sx q[0];
rz(3.0271544) q[0];
x q[1];
rz(-0.89739563) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(2.430254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0205295) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(-1.8706277) q[1];
rz(1.5548607) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(0.62409419) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(-1.3895967) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5647854) q[2];
sx q[2];
rz(-1.9153321) q[2];
sx q[2];
rz(-2.92958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2495888) q[1];
sx q[1];
rz(-0.87000404) q[1];
sx q[1];
rz(2.4904817) q[1];
x q[2];
rz(-1.8072855) q[3];
sx q[3];
rz(-0.83179501) q[3];
sx q[3];
rz(-2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(-0.54668033) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(0.84035994) q[0];
x q[1];
rz(-1.7231862) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(1.6413123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6215366) q[1];
sx q[1];
rz(-1.5339508) q[1];
sx q[1];
rz(-2.8523793) q[1];
x q[2];
rz(-1.7095079) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988408) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.6179686) q[2];
sx q[2];
rz(-1.5856367) q[2];
sx q[2];
rz(-0.71935364) q[2];
rz(2.1605282) q[3];
sx q[3];
rz(-2.5909501) q[3];
sx q[3];
rz(0.29028374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
