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
rz(0.5956369) q[0];
sx q[0];
rz(-2.73157) q[0];
sx q[0];
rz(-0.79467839) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(5.3310634) q[1];
sx q[1];
rz(7.4330243) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7139706) q[0];
sx q[0];
rz(-0.29610482) q[0];
sx q[0];
rz(1.4480736) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2501405) q[2];
sx q[2];
rz(-2.3806664) q[2];
sx q[2];
rz(-0.93283949) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9147891) q[1];
sx q[1];
rz(-0.43693301) q[1];
sx q[1];
rz(-0.41586693) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6744995) q[3];
sx q[3];
rz(-1.8008999) q[3];
sx q[3];
rz(2.439374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4475693) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(-2.0726974) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(1.6525035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.883413) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(-2.6302443) q[0];
rz(-1.6343575) q[1];
sx q[1];
rz(-1.2751445) q[1];
sx q[1];
rz(1.7587597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16973142) q[0];
sx q[0];
rz(-1.3644553) q[0];
sx q[0];
rz(-1.1792762) q[0];
rz(-2.7956656) q[2];
sx q[2];
rz(-0.9285766) q[2];
sx q[2];
rz(-1.3042892) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0046151) q[1];
sx q[1];
rz(-1.2619234) q[1];
sx q[1];
rz(1.1501794) q[1];
x q[2];
rz(-0.52329833) q[3];
sx q[3];
rz(-1.6611668) q[3];
sx q[3];
rz(0.096631526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4649268) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(3.0386772) q[2];
rz(2.0788705) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(-0.11500558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677143) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(0.1097196) q[0];
rz(-0.32490718) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(-2.6085764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66218161) q[0];
sx q[0];
rz(-1.7087738) q[0];
sx q[0];
rz(0.29946391) q[0];
rz(-0.43230482) q[2];
sx q[2];
rz(-3.084331) q[2];
sx q[2];
rz(2.6309516) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9647081) q[1];
sx q[1];
rz(-1.0057276) q[1];
sx q[1];
rz(2.3749897) q[1];
x q[2];
rz(-0.80954062) q[3];
sx q[3];
rz(-1.24461) q[3];
sx q[3];
rz(2.1128138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32536062) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(-1.4556966) q[2];
rz(0.085518941) q[3];
sx q[3];
rz(-2.072008) q[3];
sx q[3];
rz(0.94676179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(-1.5707312) q[1];
sx q[1];
rz(-1.1474362) q[1];
sx q[1];
rz(-3.0991203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52666506) q[0];
sx q[0];
rz(-0.4836429) q[0];
sx q[0];
rz(-2.7381104) q[0];
x q[1];
rz(1.3971523) q[2];
sx q[2];
rz(-3.0499027) q[2];
sx q[2];
rz(1.9609811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5407357) q[1];
sx q[1];
rz(-1.3807978) q[1];
sx q[1];
rz(2.3488087) q[1];
rz(-pi) q[2];
rz(-0.23353429) q[3];
sx q[3];
rz(-2.326205) q[3];
sx q[3];
rz(-0.60497192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8954358) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(0.24406544) q[2];
rz(-0.79143381) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(-0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1112082) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(-0.17330387) q[0];
rz(-2.2068842) q[1];
sx q[1];
rz(-0.84988958) q[1];
sx q[1];
rz(2.8876143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.105071) q[0];
sx q[0];
rz(-1.4582086) q[0];
sx q[0];
rz(0.81775157) q[0];
rz(-1.9333657) q[2];
sx q[2];
rz(-2.1340279) q[2];
sx q[2];
rz(2.8993487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18561126) q[1];
sx q[1];
rz(-2.4479981) q[1];
sx q[1];
rz(-0.032497092) q[1];
x q[2];
rz(1.41296) q[3];
sx q[3];
rz(-0.39208347) q[3];
sx q[3];
rz(0.81721701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0737334) q[2];
sx q[2];
rz(-0.24418712) q[2];
sx q[2];
rz(-1.7013288) q[2];
rz(-2.4654147) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22236958) q[0];
sx q[0];
rz(-1.8310522) q[0];
sx q[0];
rz(2.4134912) q[0];
rz(-1.9758196) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(-0.56070915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23317223) q[0];
sx q[0];
rz(-1.5058555) q[0];
sx q[0];
rz(3.0149595) q[0];
rz(1.4523023) q[2];
sx q[2];
rz(-1.0241707) q[2];
sx q[2];
rz(-0.34260633) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6627359) q[1];
sx q[1];
rz(-2.6881892) q[1];
sx q[1];
rz(-1.0591566) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5674616) q[3];
sx q[3];
rz(-0.11596348) q[3];
sx q[3];
rz(0.16030333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5460983) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(2.0687436) q[2];
rz(2.8192375) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(-0.391092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76674616) q[0];
sx q[0];
rz(-1.5006737) q[0];
sx q[0];
rz(-2.7333976) q[0];
rz(2.8152668) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(-1.6654642) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8360739) q[0];
sx q[0];
rz(-1.3954432) q[0];
sx q[0];
rz(0.039964635) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0075241) q[2];
sx q[2];
rz(-1.3671698) q[2];
sx q[2];
rz(1.3290392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0261111) q[1];
sx q[1];
rz(-0.12899765) q[1];
sx q[1];
rz(-1.1327101) q[1];
rz(-2.9093698) q[3];
sx q[3];
rz(-1.4055863) q[3];
sx q[3];
rz(2.8286407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25492111) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(-2.2754748) q[2];
rz(2.4050889) q[3];
sx q[3];
rz(-2.0564506) q[3];
sx q[3];
rz(-0.88732639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12164584) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(1.8628927) q[0];
rz(-2.1389029) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-2.6967646) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072413) q[0];
sx q[0];
rz(-2.0929411) q[0];
sx q[0];
rz(0.17052167) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1110531) q[2];
sx q[2];
rz(-0.90219775) q[2];
sx q[2];
rz(-0.3102613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0792744) q[1];
sx q[1];
rz(-1.6898815) q[1];
sx q[1];
rz(-0.80715413) q[1];
x q[2];
rz(-2.5505742) q[3];
sx q[3];
rz(-2.6962387) q[3];
sx q[3];
rz(-1.9231918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7030846) q[2];
sx q[2];
rz(-1.727641) q[2];
sx q[2];
rz(-0.50602305) q[2];
rz(0.94432008) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849702) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(3.131409) q[0];
rz(-2.5306375) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(-3.0593061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1707538) q[0];
sx q[0];
rz(-2.0505095) q[0];
sx q[0];
rz(1.1520349) q[0];
rz(-0.23046979) q[2];
sx q[2];
rz(-0.90417143) q[2];
sx q[2];
rz(0.33483349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.391606) q[1];
sx q[1];
rz(-2.9175903) q[1];
sx q[1];
rz(-2.8620666) q[1];
x q[2];
rz(-2.4769267) q[3];
sx q[3];
rz(-1.3997404) q[3];
sx q[3];
rz(-1.6534352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33132195) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-0.70915478) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(-0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.822478) q[0];
sx q[0];
rz(-0.8304441) q[0];
sx q[0];
rz(-2.5122232) q[0];
rz(1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(1.9633044) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45635763) q[0];
sx q[0];
rz(-2.166232) q[0];
sx q[0];
rz(-0.49249442) q[0];
x q[1];
rz(2.0731931) q[2];
sx q[2];
rz(-1.8821239) q[2];
sx q[2];
rz(2.5466998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1078892) q[1];
sx q[1];
rz(-0.7830355) q[1];
sx q[1];
rz(1.7262001) q[1];
rz(-pi) q[2];
rz(1.0052106) q[3];
sx q[3];
rz(-2.2095888) q[3];
sx q[3];
rz(1.8106724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87469953) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(0.81864041) q[2];
rz(1.5107752) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(0.40217933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4928987) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(-1.2661487) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(-2.8826089) q[2];
sx q[2];
rz(-2.2296446) q[2];
sx q[2];
rz(-0.17204816) q[2];
rz(0.18465445) q[3];
sx q[3];
rz(-1.6600556) q[3];
sx q[3];
rz(-0.32515812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
