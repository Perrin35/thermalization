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
rz(2.7695739) q[0];
sx q[0];
rz(-0.3516742) q[0];
sx q[0];
rz(-3.0863808) q[0];
rz(-1.6959603) q[1];
sx q[1];
rz(4.0449528) q[1];
sx q[1];
rz(9.5587211) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3066912) q[0];
sx q[0];
rz(-1.6707695) q[0];
sx q[0];
rz(2.9213219) q[0];
x q[1];
rz(0.29965286) q[2];
sx q[2];
rz(-2.5404394) q[2];
sx q[2];
rz(-0.047182949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6559927) q[1];
sx q[1];
rz(-0.47708407) q[1];
sx q[1];
rz(-0.78039767) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.03911253) q[3];
sx q[3];
rz(-2.0958423) q[3];
sx q[3];
rz(-1.0421747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17406164) q[2];
sx q[2];
rz(-1.920819) q[2];
sx q[2];
rz(-0.81452149) q[2];
rz(-0.19541611) q[3];
sx q[3];
rz(-0.82868367) q[3];
sx q[3];
rz(1.0403847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.279351) q[0];
sx q[0];
rz(-2.8793654) q[0];
sx q[0];
rz(0.97214118) q[0];
rz(0.60130087) q[1];
sx q[1];
rz(-0.96527946) q[1];
sx q[1];
rz(-1.7960637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9321971) q[0];
sx q[0];
rz(-2.4426113) q[0];
sx q[0];
rz(-2.2369034) q[0];
x q[1];
rz(2.6048772) q[2];
sx q[2];
rz(-1.8036575) q[2];
sx q[2];
rz(-2.1888806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7267725) q[1];
sx q[1];
rz(-1.1632246) q[1];
sx q[1];
rz(-0.15482082) q[1];
rz(1.5645157) q[3];
sx q[3];
rz(-1.7518861) q[3];
sx q[3];
rz(0.6015425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1874275) q[2];
sx q[2];
rz(-1.2792055) q[2];
sx q[2];
rz(0.97174755) q[2];
rz(-2.1238972) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(3.0903604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.510842) q[0];
sx q[0];
rz(-0.96931163) q[0];
sx q[0];
rz(-2.4208659) q[0];
rz(-3.0156056) q[1];
sx q[1];
rz(-2.4469913) q[1];
sx q[1];
rz(0.40649498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.061977) q[0];
sx q[0];
rz(-2.4361389) q[0];
sx q[0];
rz(-1.0951701) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1265246) q[2];
sx q[2];
rz(-0.90733084) q[2];
sx q[2];
rz(-0.86554722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5169591) q[1];
sx q[1];
rz(-2.0849094) q[1];
sx q[1];
rz(-0.89526432) q[1];
x q[2];
rz(2.2307736) q[3];
sx q[3];
rz(-2.6916457) q[3];
sx q[3];
rz(-1.4773953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6005818) q[2];
sx q[2];
rz(-1.3449679) q[2];
sx q[2];
rz(0.37008944) q[2];
rz(-0.078941405) q[3];
sx q[3];
rz(-2.168096) q[3];
sx q[3];
rz(-0.48648849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9321891) q[0];
sx q[0];
rz(-1.7904733) q[0];
sx q[0];
rz(2.6655647) q[0];
rz(1.1141106) q[1];
sx q[1];
rz(-0.94534355) q[1];
sx q[1];
rz(-1.5155189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0335122) q[0];
sx q[0];
rz(-0.42244222) q[0];
sx q[0];
rz(2.6080934) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59882382) q[2];
sx q[2];
rz(-0.90617563) q[2];
sx q[2];
rz(1.7395626) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21525258) q[1];
sx q[1];
rz(-2.7801792) q[1];
sx q[1];
rz(0.24140668) q[1];
x q[2];
rz(-1.2021121) q[3];
sx q[3];
rz(-0.75206176) q[3];
sx q[3];
rz(-2.3526109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6473306) q[2];
sx q[2];
rz(-1.6202972) q[2];
sx q[2];
rz(-2.1935513) q[2];
rz(0.4387795) q[3];
sx q[3];
rz(-2.572757) q[3];
sx q[3];
rz(-0.89890629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62094837) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(-2.7137252) q[0];
rz(0.95871344) q[1];
sx q[1];
rz(-1.9045279) q[1];
sx q[1];
rz(1.0520891) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3593345) q[0];
sx q[0];
rz(-2.3365031) q[0];
sx q[0];
rz(-0.87053086) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38340815) q[2];
sx q[2];
rz(-2.0729625) q[2];
sx q[2];
rz(2.2344207) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0316443) q[1];
sx q[1];
rz(-2.5606025) q[1];
sx q[1];
rz(1.7419001) q[1];
x q[2];
rz(0.73486272) q[3];
sx q[3];
rz(-0.70847337) q[3];
sx q[3];
rz(2.3264309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.170257) q[2];
sx q[2];
rz(-1.8361788) q[2];
sx q[2];
rz(2.7830284) q[2];
rz(-1.9606918) q[3];
sx q[3];
rz(-0.87555331) q[3];
sx q[3];
rz(-2.6023279) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3233258) q[0];
sx q[0];
rz(-1.2610672) q[0];
sx q[0];
rz(-2.4928424) q[0];
rz(0.58748856) q[1];
sx q[1];
rz(-1.3222539) q[1];
sx q[1];
rz(2.9041451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96188155) q[0];
sx q[0];
rz(-1.474232) q[0];
sx q[0];
rz(1.5386816) q[0];
rz(1.6726137) q[2];
sx q[2];
rz(-1.8378496) q[2];
sx q[2];
rz(1.36324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0157369) q[1];
sx q[1];
rz(-2.3185711) q[1];
sx q[1];
rz(-1.5344844) q[1];
x q[2];
rz(1.009672) q[3];
sx q[3];
rz(-2.3497407) q[3];
sx q[3];
rz(-0.71969024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6610403) q[2];
sx q[2];
rz(-2.48017) q[2];
sx q[2];
rz(0.95575571) q[2];
rz(-1.7620979) q[3];
sx q[3];
rz(-0.62164128) q[3];
sx q[3];
rz(-0.1563589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5697923) q[0];
sx q[0];
rz(-0.0026230165) q[0];
sx q[0];
rz(-3.0963335) q[0];
rz(2.3472002) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(2.9387567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58205523) q[0];
sx q[0];
rz(-1.4180611) q[0];
sx q[0];
rz(2.5688897) q[0];
rz(-pi) q[1];
rz(1.3812106) q[2];
sx q[2];
rz(-2.2933497) q[2];
sx q[2];
rz(1.0740395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0326621) q[1];
sx q[1];
rz(-1.0294518) q[1];
sx q[1];
rz(1.4757266) q[1];
rz(-1.6129812) q[3];
sx q[3];
rz(-2.4191471) q[3];
sx q[3];
rz(-0.23596059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63991919) q[2];
sx q[2];
rz(-1.4036274) q[2];
sx q[2];
rz(2.5301834) q[2];
rz(2.0407138) q[3];
sx q[3];
rz(-2.5067063) q[3];
sx q[3];
rz(-0.27975217) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47073498) q[0];
sx q[0];
rz(-1.1962471) q[0];
sx q[0];
rz(-1.083495) q[0];
rz(-0.96393839) q[1];
sx q[1];
rz(-2.0699392) q[1];
sx q[1];
rz(-0.49577698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14228798) q[0];
sx q[0];
rz(-1.5919627) q[0];
sx q[0];
rz(-1.4975966) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8476358) q[2];
sx q[2];
rz(-1.7466892) q[2];
sx q[2];
rz(-1.7664282) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26357061) q[1];
sx q[1];
rz(-1.9847893) q[1];
sx q[1];
rz(-0.94697006) q[1];
rz(-pi) q[2];
rz(1.5257902) q[3];
sx q[3];
rz(-0.58820217) q[3];
sx q[3];
rz(-2.3308995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44020161) q[2];
sx q[2];
rz(-0.76875606) q[2];
sx q[2];
rz(2.4526147) q[2];
rz(-1.6280599) q[3];
sx q[3];
rz(-0.83008927) q[3];
sx q[3];
rz(2.1400863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772407) q[0];
sx q[0];
rz(-1.5654726) q[0];
sx q[0];
rz(1.0871357) q[0];
rz(0.76308909) q[1];
sx q[1];
rz(-1.3975846) q[1];
sx q[1];
rz(2.9147002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927923) q[0];
sx q[0];
rz(-0.88946402) q[0];
sx q[0];
rz(-0.28980227) q[0];
rz(1.2490602) q[2];
sx q[2];
rz(-1.2870803) q[2];
sx q[2];
rz(2.4987881) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58851885) q[1];
sx q[1];
rz(-1.6071734) q[1];
sx q[1];
rz(-0.22471551) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6361106) q[3];
sx q[3];
rz(-2.1104321) q[3];
sx q[3];
rz(1.1330549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.521296) q[2];
sx q[2];
rz(-0.93775788) q[2];
sx q[2];
rz(0.94989455) q[2];
rz(-1.0319483) q[3];
sx q[3];
rz(-0.87939206) q[3];
sx q[3];
rz(-2.7040645) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501605) q[0];
sx q[0];
rz(-0.10364769) q[0];
sx q[0];
rz(1.1520804) q[0];
rz(3.0488455) q[1];
sx q[1];
rz(-2.4536965) q[1];
sx q[1];
rz(-0.50382096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60022416) q[0];
sx q[0];
rz(-1.0447518) q[0];
sx q[0];
rz(-3.0843488) q[0];
x q[1];
rz(1.3400192) q[2];
sx q[2];
rz(-1.3084931) q[2];
sx q[2];
rz(-0.80517804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.720394) q[1];
sx q[1];
rz(-2.6662146) q[1];
sx q[1];
rz(0.092751547) q[1];
rz(-2.2600365) q[3];
sx q[3];
rz(-2.6668352) q[3];
sx q[3];
rz(0.17533824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6759701) q[2];
sx q[2];
rz(-1.114782) q[2];
sx q[2];
rz(-2.5176804) q[2];
rz(2.5850249) q[3];
sx q[3];
rz(-2.4534295) q[3];
sx q[3];
rz(2.9437959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.9482166) q[0];
sx q[0];
rz(-2.2103136) q[0];
sx q[0];
rz(0.92934004) q[0];
rz(-0.14840645) q[1];
sx q[1];
rz(-1.8971309) q[1];
sx q[1];
rz(2.94577) q[1];
rz(2.6779867) q[2];
sx q[2];
rz(-1.4624034) q[2];
sx q[2];
rz(0.69093888) q[2];
rz(3.0608724) q[3];
sx q[3];
rz(-1.29946) q[3];
sx q[3];
rz(2.4429532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
