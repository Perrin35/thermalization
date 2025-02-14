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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5001427) q[0];
sx q[0];
rz(-3.1063188) q[0];
sx q[0];
rz(0.39536898) q[0];
rz(-pi) q[1];
x q[1];
rz(1.38491) q[2];
sx q[2];
rz(-1.3449525) q[2];
sx q[2];
rz(2.5511044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0684587) q[1];
sx q[1];
rz(-1.3737118) q[1];
sx q[1];
rz(2.3545594) q[1];
x q[2];
rz(-1.9463948) q[3];
sx q[3];
rz(-1.7197242) q[3];
sx q[3];
rz(-2.4747839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50769606) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(2.2155649) q[2];
rz(-1.5121459) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(-1.9984455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40874261) q[0];
sx q[0];
rz(-2.4091305) q[0];
sx q[0];
rz(0.23250411) q[0];
rz(-2.0022424) q[1];
sx q[1];
rz(-1.5683697) q[1];
sx q[1];
rz(-2.2154636) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6059604) q[0];
sx q[0];
rz(-1.7754007) q[0];
sx q[0];
rz(1.800941) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4732292) q[2];
sx q[2];
rz(-1.8969403) q[2];
sx q[2];
rz(2.9680684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7826816) q[1];
sx q[1];
rz(-2.5645683) q[1];
sx q[1];
rz(0.38459528) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4898718) q[3];
sx q[3];
rz(-0.67595302) q[3];
sx q[3];
rz(0.60871688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53655475) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(2.3489595) q[2];
rz(-0.41147453) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(-0.80504942) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3181535) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(2.8807628) q[0];
rz(0.28314319) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(-2.0859437) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84427858) q[0];
sx q[0];
rz(-1.3998654) q[0];
sx q[0];
rz(-0.92428446) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1364724) q[2];
sx q[2];
rz(-0.63314204) q[2];
sx q[2];
rz(2.1767165) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.187589) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(0.61567125) q[2];
sx q[2];
rz(-2.0522223) q[2];
sx q[2];
rz(1.8833359) q[2];
rz(-1.0387756) q[3];
sx q[3];
rz(-2.153331) q[3];
sx q[3];
rz(1.725215) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1309758) q[0];
sx q[0];
rz(-1.6058141) q[0];
sx q[0];
rz(-0.11949874) q[0];
rz(-0.15826982) q[1];
sx q[1];
rz(-0.4522849) q[1];
sx q[1];
rz(-1.7012885) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570069) q[0];
sx q[0];
rz(-2.2586926) q[0];
sx q[0];
rz(2.6683065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9057299) q[2];
sx q[2];
rz(-1.6455629) q[2];
sx q[2];
rz(-2.3521545) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5938737) q[1];
sx q[1];
rz(-2.5903194) q[1];
sx q[1];
rz(2.6216595) q[1];
x q[2];
rz(1.0802832) q[3];
sx q[3];
rz(-1.2592096) q[3];
sx q[3];
rz(2.8924242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4629472) q[2];
sx q[2];
rz(-0.1354278) q[2];
sx q[2];
rz(-2.7079008) q[2];
rz(-2.4667451) q[3];
sx q[3];
rz(-1.0899455) q[3];
sx q[3];
rz(-1.9851782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40815142) q[0];
sx q[0];
rz(-2.0858522) q[0];
sx q[0];
rz(0.59930402) q[0];
rz(2.377548) q[1];
sx q[1];
rz(-1.0300449) q[1];
sx q[1];
rz(0.61242551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2297933) q[0];
sx q[0];
rz(-0.62055885) q[0];
sx q[0];
rz(1.3156652) q[0];
rz(-pi) q[1];
rz(-1.5715309) q[2];
sx q[2];
rz(-2.2876566) q[2];
sx q[2];
rz(1.5897205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45578411) q[1];
sx q[1];
rz(-2.1791885) q[1];
sx q[1];
rz(-1.366282) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3633943) q[3];
sx q[3];
rz(-2.457387) q[3];
sx q[3];
rz(0.15008959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(2.6714228) q[2];
rz(0.36559513) q[3];
sx q[3];
rz(-1.1919034) q[3];
sx q[3];
rz(2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81569165) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(-2.4796487) q[0];
rz(-3.0603307) q[1];
sx q[1];
rz(-0.47835246) q[1];
sx q[1];
rz(-0.14019664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2574923) q[0];
sx q[0];
rz(-1.1123158) q[0];
sx q[0];
rz(1.3263339) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95023167) q[2];
sx q[2];
rz(-1.3060313) q[2];
sx q[2];
rz(2.1193412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.081379024) q[1];
sx q[1];
rz(-1.5790931) q[1];
sx q[1];
rz(-1.0992728) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7532888) q[3];
sx q[3];
rz(-1.914822) q[3];
sx q[3];
rz(1.1870015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7774272) q[2];
sx q[2];
rz(-1.825793) q[2];
sx q[2];
rz(0.3248471) q[2];
rz(2.9835564) q[3];
sx q[3];
rz(-2.7285748) q[3];
sx q[3];
rz(1.0154999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82191104) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(0.88777375) q[0];
rz(-0.38331389) q[1];
sx q[1];
rz(-0.8822459) q[1];
sx q[1];
rz(-2.9820138) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68113806) q[0];
sx q[0];
rz(-0.41550203) q[0];
sx q[0];
rz(2.2592708) q[0];
rz(2.1453484) q[2];
sx q[2];
rz(-1.0932845) q[2];
sx q[2];
rz(1.6928796) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13028045) q[1];
sx q[1];
rz(-0.39825687) q[1];
sx q[1];
rz(2.8058737) q[1];
x q[2];
rz(-0.51377138) q[3];
sx q[3];
rz(-2.8119866) q[3];
sx q[3];
rz(0.20121516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54922709) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(1.7249736) q[2];
rz(1.8727411) q[3];
sx q[3];
rz(-0.78059355) q[3];
sx q[3];
rz(-0.95808539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45605993) q[0];
sx q[0];
rz(-1.1153509) q[0];
sx q[0];
rz(2.2400895) q[0];
rz(0.69425663) q[1];
sx q[1];
rz(-1.4638823) q[1];
sx q[1];
rz(1.136397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0471056) q[0];
sx q[0];
rz(-1.4102954) q[0];
sx q[0];
rz(0.032139924) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56635277) q[2];
sx q[2];
rz(-2.512062) q[2];
sx q[2];
rz(2.9734263) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8759114) q[1];
sx q[1];
rz(-0.87603518) q[1];
sx q[1];
rz(-3.0819703) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5707022) q[3];
sx q[3];
rz(-0.81557214) q[3];
sx q[3];
rz(-2.3237128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0006492) q[2];
sx q[2];
rz(-1.2691754) q[2];
sx q[2];
rz(0.81306523) q[2];
rz(-2.5874169) q[3];
sx q[3];
rz(-1.412609) q[3];
sx q[3];
rz(-0.36177844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56162214) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(2.9680874) q[0];
rz(2.7165727) q[1];
sx q[1];
rz(-1.5360906) q[1];
sx q[1];
rz(-2.3120211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3542784) q[0];
sx q[0];
rz(-1.4256546) q[0];
sx q[0];
rz(-2.4046242) q[0];
x q[1];
rz(1.6116905) q[2];
sx q[2];
rz(-1.6171793) q[2];
sx q[2];
rz(-1.0414909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47118716) q[1];
sx q[1];
rz(-2.1758432) q[1];
sx q[1];
rz(-1.3208645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8052727) q[3];
sx q[3];
rz(-1.5097396) q[3];
sx q[3];
rz(-1.2741736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35036626) q[2];
sx q[2];
rz(-1.3000877) q[2];
sx q[2];
rz(1.4607956) q[2];
rz(0.99460498) q[3];
sx q[3];
rz(-0.9959144) q[3];
sx q[3];
rz(-0.88609707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.17722546) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(-0.10203578) q[0];
rz(2.521934) q[1];
sx q[1];
rz(-1.4980059) q[1];
sx q[1];
rz(0.59725753) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12349081) q[0];
sx q[0];
rz(-0.97663701) q[0];
sx q[0];
rz(1.7325364) q[0];
rz(-pi) q[1];
rz(-1.4240392) q[2];
sx q[2];
rz(-1.2191311) q[2];
sx q[2];
rz(2.7232225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0400552) q[1];
sx q[1];
rz(-1.005508) q[1];
sx q[1];
rz(-0.043626309) q[1];
x q[2];
rz(0.47145505) q[3];
sx q[3];
rz(-1.5416036) q[3];
sx q[3];
rz(-1.4271133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9336046) q[2];
sx q[2];
rz(-0.36078578) q[2];
sx q[2];
rz(-0.92850816) q[2];
rz(-1.1601296) q[3];
sx q[3];
rz(-1.4312294) q[3];
sx q[3];
rz(-0.76735705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3930897) q[0];
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
rz(2.7990758) q[3];
sx q[3];
rz(-0.90382648) q[3];
sx q[3];
rz(3.1161502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
