OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.79787624) q[0];
sx q[0];
rz(-0.98348445) q[0];
sx q[0];
rz(2.9498192) q[0];
rz(-2.6311488) q[1];
sx q[1];
rz(-1.3671083) q[1];
sx q[1];
rz(-1.7459315) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1045525) q[0];
sx q[0];
rz(-1.5382447) q[0];
sx q[0];
rz(-1.5572059) q[0];
x q[1];
rz(1.38491) q[2];
sx q[2];
rz(-1.7966401) q[2];
sx q[2];
rz(-2.5511044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4499481) q[1];
sx q[1];
rz(-0.80300036) q[1];
sx q[1];
rz(-1.8464441) q[1];
x q[2];
rz(1.1951978) q[3];
sx q[3];
rz(-1.7197242) q[3];
sx q[3];
rz(-2.4747839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6338966) q[2];
sx q[2];
rz(-2.8140929) q[2];
sx q[2];
rz(2.2155649) q[2];
rz(1.6294468) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(-1.9984455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40874261) q[0];
sx q[0];
rz(-2.4091305) q[0];
sx q[0];
rz(-0.23250411) q[0];
rz(2.0022424) q[1];
sx q[1];
rz(-1.573223) q[1];
sx q[1];
rz(0.92612902) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53563228) q[0];
sx q[0];
rz(-1.366192) q[0];
sx q[0];
rz(-1.3406517) q[0];
x q[1];
rz(-1.6683634) q[2];
sx q[2];
rz(-1.8969403) q[2];
sx q[2];
rz(-0.17352428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3327712) q[1];
sx q[1];
rz(-1.0406245) q[1];
sx q[1];
rz(1.8103241) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5739734) q[3];
sx q[3];
rz(-1.960037) q[3];
sx q[3];
rz(-0.4252227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6050379) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(2.3489595) q[2];
rz(2.7301181) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3181535) q[0];
sx q[0];
rz(-2.1397488) q[0];
sx q[0];
rz(0.26082984) q[0];
rz(0.28314319) q[1];
sx q[1];
rz(-1.4500376) q[1];
sx q[1];
rz(2.0859437) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84427858) q[0];
sx q[0];
rz(-1.7417272) q[0];
sx q[0];
rz(0.92428446) q[0];
rz(-pi) q[1];
rz(2.5084569) q[2];
sx q[2];
rz(-1.5677668) q[2];
sx q[2];
rz(-2.5315447) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95400364) q[1];
sx q[1];
rz(-0.84792811) q[1];
sx q[1];
rz(1.4299117) q[1];
rz(-pi) q[2];
rz(2.0664472) q[3];
sx q[3];
rz(-1.5167674) q[3];
sx q[3];
rz(-0.13755218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5259214) q[2];
sx q[2];
rz(-2.0522223) q[2];
sx q[2];
rz(-1.2582568) q[2];
rz(2.1028171) q[3];
sx q[3];
rz(-2.153331) q[3];
sx q[3];
rz(1.725215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0106169) q[0];
sx q[0];
rz(-1.5357786) q[0];
sx q[0];
rz(0.11949874) q[0];
rz(-2.9833228) q[1];
sx q[1];
rz(-2.6893078) q[1];
sx q[1];
rz(1.4403042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40057221) q[0];
sx q[0];
rz(-1.9306679) q[0];
sx q[0];
rz(0.82525702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.079147804) q[2];
sx q[2];
rz(-1.2368349) q[2];
sx q[2];
rz(0.755366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56934565) q[1];
sx q[1];
rz(-1.3075446) q[1];
sx q[1];
rz(-0.49016989) q[1];
rz(-pi) q[2];
rz(0.3500895) q[3];
sx q[3];
rz(-2.0357657) q[3];
sx q[3];
rz(1.657682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67864546) q[2];
sx q[2];
rz(-3.0061649) q[2];
sx q[2];
rz(0.43369183) q[2];
rz(2.4667451) q[3];
sx q[3];
rz(-1.0899455) q[3];
sx q[3];
rz(-1.1564144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7334412) q[0];
sx q[0];
rz(-1.0557405) q[0];
sx q[0];
rz(0.59930402) q[0];
rz(0.76404461) q[1];
sx q[1];
rz(-2.1115477) q[1];
sx q[1];
rz(0.61242551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5915039) q[0];
sx q[0];
rz(-1.4235125) q[0];
sx q[0];
rz(-2.1758737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71686042) q[2];
sx q[2];
rz(-1.5702425) q[2];
sx q[2];
rz(-3.1231511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99701559) q[1];
sx q[1];
rz(-1.738228) q[1];
sx q[1];
rz(-2.5232878) q[1];
rz(-pi) q[2];
rz(2.6153485) q[3];
sx q[3];
rz(-2.0305227) q[3];
sx q[3];
rz(-2.3731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(-0.4701699) q[2];
rz(2.7759975) q[3];
sx q[3];
rz(-1.9496893) q[3];
sx q[3];
rz(2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325901) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(-0.66194397) q[0];
rz(0.081261948) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-3.001396) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707342) q[0];
sx q[0];
rz(-2.6261289) q[0];
sx q[0];
rz(2.6856642) q[0];
rz(-2.0070932) q[2];
sx q[2];
rz(-2.4738174) q[2];
sx q[2];
rz(0.89950022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0602136) q[1];
sx q[1];
rz(-1.5624996) q[1];
sx q[1];
rz(-2.0423198) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46886985) q[3];
sx q[3];
rz(-0.38772407) q[3];
sx q[3];
rz(-0.68634168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7774272) q[2];
sx q[2];
rz(-1.3157996) q[2];
sx q[2];
rz(-2.8167456) q[2];
rz(-2.9835564) q[3];
sx q[3];
rz(-2.7285748) q[3];
sx q[3];
rz(-1.0154999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3196816) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(-0.88777375) q[0];
rz(0.38331389) q[1];
sx q[1];
rz(-2.2593468) q[1];
sx q[1];
rz(0.15957889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051285714) q[0];
sx q[0];
rz(-1.8877827) q[0];
sx q[0];
rz(-0.27329926) q[0];
x q[1];
rz(2.5891807) q[2];
sx q[2];
rz(-2.0744951) q[2];
sx q[2];
rz(0.1671065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9099906) q[1];
sx q[1];
rz(-1.9456777) q[1];
sx q[1];
rz(1.433062) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4042312) q[3];
sx q[3];
rz(-1.2850396) q[3];
sx q[3];
rz(-0.33657246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54922709) q[2];
sx q[2];
rz(-1.1575907) q[2];
sx q[2];
rz(-1.4166191) q[2];
rz(-1.2688515) q[3];
sx q[3];
rz(-2.3609991) q[3];
sx q[3];
rz(0.95808539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45605993) q[0];
sx q[0];
rz(-1.1153509) q[0];
sx q[0];
rz(-2.2400895) q[0];
rz(-0.69425663) q[1];
sx q[1];
rz(-1.4638823) q[1];
sx q[1];
rz(2.0051956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0944871) q[0];
sx q[0];
rz(-1.4102954) q[0];
sx q[0];
rz(-3.1094527) q[0];
rz(2.5904584) q[2];
sx q[2];
rz(-1.8922085) q[2];
sx q[2];
rz(-1.8773735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.34331218) q[1];
sx q[1];
rz(-1.5250051) q[1];
sx q[1];
rz(-0.87516038) q[1];
rz(1.0497007) q[3];
sx q[3];
rz(-0.9113833) q[3];
sx q[3];
rz(-1.5707317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14094341) q[2];
sx q[2];
rz(-1.8724172) q[2];
sx q[2];
rz(0.81306523) q[2];
rz(2.5874169) q[3];
sx q[3];
rz(-1.7289836) q[3];
sx q[3];
rz(2.7798142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56162214) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(2.9680874) q[0];
rz(0.42501998) q[1];
sx q[1];
rz(-1.5360906) q[1];
sx q[1];
rz(-0.82957155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6529675) q[0];
sx q[0];
rz(-2.2982631) q[0];
sx q[0];
rz(-1.7656816) q[0];
rz(2.4194952) q[2];
sx q[2];
rz(-3.0797662) q[2];
sx q[2];
rz(-0.31844469) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8977938) q[1];
sx q[1];
rz(-1.3659371) q[1];
sx q[1];
rz(0.61989354) q[1];
rz(-pi) q[2];
rz(-2.8052727) q[3];
sx q[3];
rz(-1.5097396) q[3];
sx q[3];
rz(1.8674191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35036626) q[2];
sx q[2];
rz(-1.3000877) q[2];
sx q[2];
rz(1.680797) q[2];
rz(-2.1469877) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(0.88609707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.17722546) q[0];
sx q[0];
rz(-0.57562861) q[0];
sx q[0];
rz(0.10203578) q[0];
rz(-0.61965865) q[1];
sx q[1];
rz(-1.4980059) q[1];
sx q[1];
rz(-2.5443351) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16011691) q[0];
sx q[0];
rz(-0.61321027) q[0];
sx q[0];
rz(0.23399467) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7623457) q[2];
sx q[2];
rz(-2.7617117) q[2];
sx q[2];
rz(2.3178562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44587943) q[1];
sx q[1];
rz(-1.5339601) q[1];
sx q[1];
rz(-2.1365154) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0773874) q[3];
sx q[3];
rz(-2.6693025) q[3];
sx q[3];
rz(0.20087584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9336046) q[2];
sx q[2];
rz(-0.36078578) q[2];
sx q[2];
rz(2.2130845) q[2];
rz(1.9814631) q[3];
sx q[3];
rz(-1.4312294) q[3];
sx q[3];
rz(-0.76735705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74850294) q[0];
sx q[0];
rz(-2.4414283) q[0];
sx q[0];
rz(0.23165942) q[0];
rz(1.1275935) q[1];
sx q[1];
rz(-2.6078106) q[1];
sx q[1];
rz(-0.003905205) q[1];
rz(3.0697974) q[2];
sx q[2];
rz(-2.2459002) q[2];
sx q[2];
rz(-0.33159524) q[2];
rz(-1.9740022) q[3];
sx q[3];
rz(-2.4039563) q[3];
sx q[3];
rz(-2.6441426) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
