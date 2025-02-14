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
rz(0.51044381) q[1];
sx q[1];
rz(4.508701) q[1];
sx q[1];
rz(8.0291168) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0370402) q[0];
sx q[0];
rz(-1.5382447) q[0];
sx q[0];
rz(-1.5843868) q[0];
rz(1.38491) q[2];
sx q[2];
rz(-1.7966401) q[2];
sx q[2];
rz(0.5904883) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.073133999) q[1];
sx q[1];
rz(-1.7678808) q[1];
sx q[1];
rz(2.3545594) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9816876) q[3];
sx q[3];
rz(-1.9420338) q[3];
sx q[3];
rz(-2.1791636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6338966) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(2.2155649) q[2];
rz(1.5121459) q[3];
sx q[3];
rz(-2.47086) q[3];
sx q[3];
rz(1.1431471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73285) q[0];
sx q[0];
rz(-0.73246211) q[0];
sx q[0];
rz(2.9090885) q[0];
rz(1.1393503) q[1];
sx q[1];
rz(-1.573223) q[1];
sx q[1];
rz(2.2154636) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6059604) q[0];
sx q[0];
rz(-1.366192) q[0];
sx q[0];
rz(1.3406517) q[0];
x q[1];
rz(2.813999) q[2];
sx q[2];
rz(-1.6632051) q[2];
sx q[2];
rz(-1.365923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3327712) q[1];
sx q[1];
rz(-1.0406245) q[1];
sx q[1];
rz(-1.3312686) q[1];
x q[2];
rz(2.0235463) q[3];
sx q[3];
rz(-2.0914075) q[3];
sx q[3];
rz(-1.3829766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53655475) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(-0.79263318) q[2];
rz(-2.7301181) q[3];
sx q[3];
rz(-0.94235197) q[3];
sx q[3];
rz(2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8234392) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(2.8807628) q[0];
rz(2.8584495) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(-1.0556489) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94819647) q[0];
sx q[0];
rz(-2.4760111) q[0];
sx q[0];
rz(1.8498621) q[0];
x q[1];
rz(0.63313578) q[2];
sx q[2];
rz(-1.5738259) q[2];
sx q[2];
rz(0.61004794) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71034096) q[1];
sx q[1];
rz(-1.6762937) q[1];
sx q[1];
rz(-2.4137817) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0751455) q[3];
sx q[3];
rz(-1.6248253) q[3];
sx q[3];
rz(-0.13755218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5259214) q[2];
sx q[2];
rz(-1.0893704) q[2];
sx q[2];
rz(-1.2582568) q[2];
rz(2.1028171) q[3];
sx q[3];
rz(-0.98826161) q[3];
sx q[3];
rz(-1.725215) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.7012885) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845857) q[0];
sx q[0];
rz(-2.2586926) q[0];
sx q[0];
rz(-2.6683065) q[0];
x q[1];
rz(1.2358627) q[2];
sx q[2];
rz(-1.4960297) q[2];
sx q[2];
rz(0.78943816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54771891) q[1];
sx q[1];
rz(-0.55127326) q[1];
sx q[1];
rz(-2.6216595) q[1];
x q[2];
rz(-0.97109183) q[3];
sx q[3];
rz(-0.57421847) q[3];
sx q[3];
rz(-0.80050877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67864546) q[2];
sx q[2];
rz(-0.1354278) q[2];
sx q[2];
rz(2.7079008) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40815142) q[0];
sx q[0];
rz(-2.0858522) q[0];
sx q[0];
rz(-2.5422886) q[0];
rz(0.76404461) q[1];
sx q[1];
rz(-1.0300449) q[1];
sx q[1];
rz(2.5291671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5915039) q[0];
sx q[0];
rz(-1.7180802) q[0];
sx q[0];
rz(0.96571897) q[0];
rz(-pi) q[1];
rz(-1.5700617) q[2];
sx q[2];
rz(-2.2876566) q[2];
sx q[2];
rz(1.5518722) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3376895) q[1];
sx q[1];
rz(-0.63769996) q[1];
sx q[1];
rz(-2.8578812) q[1];
rz(0.77819838) q[3];
sx q[3];
rz(-2.457387) q[3];
sx q[3];
rz(-2.9915031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6969488) q[2];
sx q[2];
rz(-2.3418043) q[2];
sx q[2];
rz(-2.6714228) q[2];
rz(0.36559513) q[3];
sx q[3];
rz(-1.9496893) q[3];
sx q[3];
rz(0.61075926) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325901) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(-2.4796487) q[0];
rz(-0.081261948) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-0.14019664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707342) q[0];
sx q[0];
rz(-2.6261289) q[0];
sx q[0];
rz(2.6856642) q[0];
x q[1];
rz(-2.191361) q[2];
sx q[2];
rz(-1.8355614) q[2];
sx q[2];
rz(2.1193412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4936476) q[1];
sx q[1];
rz(-2.0423023) q[1];
sx q[1];
rz(-0.00931298) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7922157) q[3];
sx q[3];
rz(-1.3991068) q[3];
sx q[3];
rz(2.6956345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7774272) q[2];
sx q[2];
rz(-1.3157996) q[2];
sx q[2];
rz(-0.3248471) q[2];
rz(-0.15803629) q[3];
sx q[3];
rz(-2.7285748) q[3];
sx q[3];
rz(-2.1260927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3196816) q[0];
sx q[0];
rz(-3.0653636) q[0];
sx q[0];
rz(-2.2538189) q[0];
rz(2.7582788) q[1];
sx q[1];
rz(-2.2593468) q[1];
sx q[1];
rz(2.9820138) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4604546) q[0];
sx q[0];
rz(-0.41550203) q[0];
sx q[0];
rz(0.88232188) q[0];
rz(0.99624421) q[2];
sx q[2];
rz(-1.0932845) q[2];
sx q[2];
rz(1.4487131) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3899052) q[1];
sx q[1];
rz(-1.6989105) q[1];
sx q[1];
rz(-2.7634578) q[1];
x q[2];
rz(0.51377138) q[3];
sx q[3];
rz(-2.8119866) q[3];
sx q[3];
rz(-0.20121516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54922709) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(-1.4166191) q[2];
rz(1.2688515) q[3];
sx q[3];
rz(-2.3609991) q[3];
sx q[3];
rz(-0.95808539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45605993) q[0];
sx q[0];
rz(-1.1153509) q[0];
sx q[0];
rz(0.90150315) q[0];
rz(-2.447336) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(-1.136397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6601453) q[0];
sx q[0];
rz(-1.5390696) q[0];
sx q[0];
rz(-1.7313787) q[0];
rz(1.9433697) q[2];
sx q[2];
rz(-1.0508453) q[2];
sx q[2];
rz(0.49836788) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9688837) q[1];
sx q[1];
rz(-0.69689059) q[1];
sx q[1];
rz(1.4994166) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72937596) q[3];
sx q[3];
rz(-1.9751901) q[3];
sx q[3];
rz(2.803363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0006492) q[2];
sx q[2];
rz(-1.2691754) q[2];
sx q[2];
rz(-2.3285274) q[2];
rz(-2.5874169) q[3];
sx q[3];
rz(-1.412609) q[3];
sx q[3];
rz(2.7798142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56162214) q[0];
sx q[0];
rz(-2.0273835) q[0];
sx q[0];
rz(-0.17350523) q[0];
rz(-0.42501998) q[1];
sx q[1];
rz(-1.605502) q[1];
sx q[1];
rz(2.3120211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7873142) q[0];
sx q[0];
rz(-1.4256546) q[0];
sx q[0];
rz(2.4046242) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6116905) q[2];
sx q[2];
rz(-1.6171793) q[2];
sx q[2];
rz(-1.0414909) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0922488) q[1];
sx q[1];
rz(-2.4929765) q[1];
sx q[1];
rz(0.34347024) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5061257) q[3];
sx q[3];
rz(-1.9064649) q[3];
sx q[3];
rz(-0.27529374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7912264) q[2];
sx q[2];
rz(-1.841505) q[2];
sx q[2];
rz(-1.680797) q[2];
rz(-0.99460498) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(2.2554956) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9643672) q[0];
sx q[0];
rz(-0.57562861) q[0];
sx q[0];
rz(-0.10203578) q[0];
rz(-2.521934) q[1];
sx q[1];
rz(-1.6435868) q[1];
sx q[1];
rz(0.59725753) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6031979) q[0];
sx q[0];
rz(-1.4369597) q[0];
sx q[0];
rz(2.5413245) q[0];
x q[1];
rz(2.7623457) q[2];
sx q[2];
rz(-0.37988099) q[2];
sx q[2];
rz(-0.82373649) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.44587943) q[1];
sx q[1];
rz(-1.5339601) q[1];
sx q[1];
rz(1.0050773) q[1];
rz(-pi) q[2];
rz(3.0773874) q[3];
sx q[3];
rz(-0.47229015) q[3];
sx q[3];
rz(2.9407168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.20798802) q[2];
sx q[2];
rz(-2.7808069) q[2];
sx q[2];
rz(0.92850816) q[2];
rz(-1.9814631) q[3];
sx q[3];
rz(-1.4312294) q[3];
sx q[3];
rz(0.76735705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930897) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(-1.1275935) q[1];
sx q[1];
rz(-0.53378202) q[1];
sx q[1];
rz(3.1376874) q[1];
rz(-3.0697974) q[2];
sx q[2];
rz(-0.89569246) q[2];
sx q[2];
rz(2.8099974) q[2];
rz(-2.7990758) q[3];
sx q[3];
rz(-2.2377662) q[3];
sx q[3];
rz(-0.025442414) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
