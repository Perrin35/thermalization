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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(-2.1314148) q[0];
rz(2.3330359) q[1];
sx q[1];
rz(-2.8454236) q[1];
sx q[1];
rz(-0.30997601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3637562) q[0];
sx q[0];
rz(-1.68206) q[0];
sx q[0];
rz(-2.7194752) q[0];
rz(2.9363382) q[2];
sx q[2];
rz(-1.5488834) q[2];
sx q[2];
rz(2.2147629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12116005) q[1];
sx q[1];
rz(-1.3749287) q[1];
sx q[1];
rz(-1.905026) q[1];
rz(-pi) q[2];
rz(-2.9397291) q[3];
sx q[3];
rz(-0.84175693) q[3];
sx q[3];
rz(-1.0224467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6628722) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(-0.09566801) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(-1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9005168) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(2.526793) q[0];
rz(-0.88848937) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(2.6420171) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7730656) q[0];
sx q[0];
rz(-1.245226) q[0];
sx q[0];
rz(0.21296176) q[0];
rz(-pi) q[1];
rz(0.52421661) q[2];
sx q[2];
rz(-1.4337863) q[2];
sx q[2];
rz(1.7643339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37812472) q[1];
sx q[1];
rz(-2.2164882) q[1];
sx q[1];
rz(-0.21685361) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3376451) q[3];
sx q[3];
rz(-1.0146552) q[3];
sx q[3];
rz(-2.3157256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2483612) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(-0.020817967) q[2];
rz(-1.2402395) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5027387) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(2.8079206) q[0];
rz(0.73792136) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(-3.0121682) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282818) q[0];
sx q[0];
rz(-2.2419562) q[0];
sx q[0];
rz(0.20361237) q[0];
rz(-pi) q[1];
rz(-2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(2.3455623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94679615) q[1];
sx q[1];
rz(-1.692564) q[1];
sx q[1];
rz(-0.090784723) q[1];
rz(-pi) q[2];
rz(-2.0468726) q[3];
sx q[3];
rz(-0.84730803) q[3];
sx q[3];
rz(-1.0915966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0176598) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(-0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081886) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(0.56513894) q[0];
rz(-2.7124229) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(0.82829222) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6713326) q[0];
sx q[0];
rz(-1.4909706) q[0];
sx q[0];
rz(0.74175394) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0847237) q[2];
sx q[2];
rz(-1.3157433) q[2];
sx q[2];
rz(-1.6549095) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.284621) q[1];
sx q[1];
rz(-1.0880252) q[1];
sx q[1];
rz(-1.2699782) q[1];
rz(-pi) q[2];
rz(-1.1513814) q[3];
sx q[3];
rz(-0.55636251) q[3];
sx q[3];
rz(1.1384979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68507489) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(-1.431541) q[2];
rz(0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-0.13570304) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(-1.3294719) q[0];
rz(-1.1520518) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(-0.00024814127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6085998) q[0];
sx q[0];
rz(-1.3260576) q[0];
sx q[0];
rz(-2.9518963) q[0];
rz(-pi) q[1];
rz(-1.0365965) q[2];
sx q[2];
rz(-1.4956258) q[2];
sx q[2];
rz(0.49803621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9226886) q[1];
sx q[1];
rz(-2.5210025) q[1];
sx q[1];
rz(0.94938486) q[1];
x q[2];
rz(1.4983921) q[3];
sx q[3];
rz(-0.84057099) q[3];
sx q[3];
rz(0.33354586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0033215) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(0.0058343466) q[2];
rz(-2.2516069) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(2.0148923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24790813) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.6538612) q[0];
rz(-0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(-0.94246513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2104032) q[0];
sx q[0];
rz(-1.382834) q[0];
sx q[0];
rz(0.24954777) q[0];
x q[1];
rz(-1.8293657) q[2];
sx q[2];
rz(-1.2814008) q[2];
sx q[2];
rz(2.0813775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8854495) q[1];
sx q[1];
rz(-1.1032618) q[1];
sx q[1];
rz(0.3216775) q[1];
x q[2];
rz(-1.0932572) q[3];
sx q[3];
rz(-0.52553383) q[3];
sx q[3];
rz(0.20697396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5421062) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(-1.8360809) q[2];
rz(-2.5944338) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(0.80593306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(3.0410774) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(2.2844792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094497546) q[0];
sx q[0];
rz(-2.408354) q[0];
sx q[0];
rz(-0.49219699) q[0];
rz(-1.521342) q[2];
sx q[2];
rz(-0.92453814) q[2];
sx q[2];
rz(2.6924999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7776612) q[1];
sx q[1];
rz(-0.45501935) q[1];
sx q[1];
rz(-0.48952405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9590274) q[3];
sx q[3];
rz(-0.60743466) q[3];
sx q[3];
rz(-1.5807815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40019217) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(-2.8391489) q[2];
rz(0.97964573) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(1.2790595) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7804467) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(3.1106023) q[0];
rz(2.567645) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(1.8642289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4784067) q[0];
sx q[0];
rz(-2.660523) q[0];
sx q[0];
rz(0.78669725) q[0];
rz(1.9979111) q[2];
sx q[2];
rz(-2.6027205) q[2];
sx q[2];
rz(2.8123735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22945089) q[1];
sx q[1];
rz(-2.0566642) q[1];
sx q[1];
rz(-3.058601) q[1];
rz(-1.2856917) q[3];
sx q[3];
rz(-1.630097) q[3];
sx q[3];
rz(-1.7948732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0370827) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(0.48745298) q[2];
rz(0.69665748) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-0.35097861) q[0];
rz(2.9674496) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(-1.7399656) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202209) q[0];
sx q[0];
rz(-1.3231771) q[0];
sx q[0];
rz(-1.4445452) q[0];
x q[1];
rz(-2.0171595) q[2];
sx q[2];
rz(-2.2702262) q[2];
sx q[2];
rz(1.0126737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.735504) q[1];
sx q[1];
rz(-1.5417409) q[1];
sx q[1];
rz(-1.7543704) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97934874) q[3];
sx q[3];
rz(-2.2205847) q[3];
sx q[3];
rz(1.5929008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-0.6366716) q[2];
rz(2.0311671) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3670032) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(3.1029347) q[1];
sx q[1];
rz(-1.6148753) q[1];
sx q[1];
rz(-1.0640594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23180873) q[0];
sx q[0];
rz(-1.5663106) q[0];
sx q[0];
rz(1.5674598) q[0];
x q[1];
rz(-0.27711192) q[2];
sx q[2];
rz(-1.2925576) q[2];
sx q[2];
rz(0.81467512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41760379) q[1];
sx q[1];
rz(-0.072712459) q[1];
sx q[1];
rz(2.5515208) q[1];
x q[2];
rz(2.6820716) q[3];
sx q[3];
rz(-1.5395313) q[3];
sx q[3];
rz(-2.7783436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(-3.0675724) q[2];
rz(-2.7600539) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(-0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030180177) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(0.61182712) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(0.50449087) q[2];
sx q[2];
rz(-1.2332543) q[2];
sx q[2];
rz(2.6674657) q[2];
rz(-2.3930876) q[3];
sx q[3];
rz(-1.2126964) q[3];
sx q[3];
rz(-2.4550635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
