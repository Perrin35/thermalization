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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64144999) q[0];
sx q[0];
rz(-3.1063188) q[0];
sx q[0];
rz(2.7462237) q[0];
rz(-2.9119266) q[2];
sx q[2];
rz(-1.7519092) q[2];
sx q[2];
rz(-1.0223946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8366791) q[1];
sx q[1];
rz(-2.3354451) q[1];
sx q[1];
rz(2.8668001) q[1];
rz(-pi) q[2];
rz(-1.9463948) q[3];
sx q[3];
rz(-1.4218685) q[3];
sx q[3];
rz(2.4747839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50769606) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(-2.2155649) q[2];
rz(-1.6294468) q[3];
sx q[3];
rz(-2.47086) q[3];
sx q[3];
rz(-1.9984455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73285) q[0];
sx q[0];
rz(-0.73246211) q[0];
sx q[0];
rz(2.9090885) q[0];
rz(-2.0022424) q[1];
sx q[1];
rz(-1.5683697) q[1];
sx q[1];
rz(-2.2154636) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8208198) q[0];
sx q[0];
rz(-0.30673393) q[0];
sx q[0];
rz(-0.83267468) q[0];
rz(-pi) q[1];
rz(-0.28042365) q[2];
sx q[2];
rz(-2.8016675) q[2];
sx q[2];
rz(-2.6715793) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8088215) q[1];
sx q[1];
rz(-1.0406245) q[1];
sx q[1];
rz(1.8103241) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4898718) q[3];
sx q[3];
rz(-2.4656396) q[3];
sx q[3];
rz(2.5328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6050379) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(-0.79263318) q[2];
rz(-2.7301181) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(-2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8234392) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(-0.26082984) q[0];
rz(0.28314319) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(1.0556489) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59883927) q[0];
sx q[0];
rz(-0.93523798) q[0];
sx q[0];
rz(-0.21297867) q[0];
rz(0.0051202444) q[2];
sx q[2];
rz(-2.5084506) q[2];
sx q[2];
rz(-0.9648762) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4312517) q[1];
sx q[1];
rz(-1.4652989) q[1];
sx q[1];
rz(2.4137817) q[1];
rz(-0.061402873) q[3];
sx q[3];
rz(-1.0759343) q[3];
sx q[3];
rz(1.4040549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61567125) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1309758) q[0];
sx q[0];
rz(-1.6058141) q[0];
sx q[0];
rz(-0.11949874) q[0];
rz(2.9833228) q[1];
sx q[1];
rz(-0.4522849) q[1];
sx q[1];
rz(-1.7012885) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8058469) q[0];
sx q[0];
rz(-2.3290538) q[0];
sx q[0];
rz(-1.0643908) q[0];
rz(-pi) q[1];
x q[1];
rz(0.079147804) q[2];
sx q[2];
rz(-1.2368349) q[2];
sx q[2];
rz(-0.755366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5938737) q[1];
sx q[1];
rz(-2.5903194) q[1];
sx q[1];
rz(2.6216595) q[1];
rz(-pi) q[2];
rz(0.3500895) q[3];
sx q[3];
rz(-1.105827) q[3];
sx q[3];
rz(1.4839107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.67864546) q[2];
sx q[2];
rz(-0.1354278) q[2];
sx q[2];
rz(-0.43369183) q[2];
rz(0.67484754) q[3];
sx q[3];
rz(-1.0899455) q[3];
sx q[3];
rz(1.1564144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.40815142) q[0];
sx q[0];
rz(-2.0858522) q[0];
sx q[0];
rz(2.5422886) q[0];
rz(-2.377548) q[1];
sx q[1];
rz(-2.1115477) q[1];
sx q[1];
rz(-2.5291671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5500887) q[0];
sx q[0];
rz(-1.7180802) q[0];
sx q[0];
rz(-2.1758737) q[0];
rz(0.00084288518) q[2];
sx q[2];
rz(-0.71686059) q[2];
sx q[2];
rz(1.5529902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1445771) q[1];
sx q[1];
rz(-1.738228) q[1];
sx q[1];
rz(-2.5232878) q[1];
rz(-pi) q[2];
rz(-2.0908085) q[3];
sx q[3];
rz(-1.1038759) q[3];
sx q[3];
rz(2.0869521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6969488) q[2];
sx q[2];
rz(-2.3418043) q[2];
sx q[2];
rz(-0.4701699) q[2];
rz(-0.36559513) q[3];
sx q[3];
rz(-1.1919034) q[3];
sx q[3];
rz(-2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81569165) q[0];
sx q[0];
rz(-0.79623643) q[0];
sx q[0];
rz(-0.66194397) q[0];
rz(0.081261948) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-3.001396) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2574923) q[0];
sx q[0];
rz(-2.0292768) q[0];
sx q[0];
rz(-1.3263339) q[0];
rz(0.32169183) q[2];
sx q[2];
rz(-0.9749229) q[2];
sx q[2];
rz(-0.3636407) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0602136) q[1];
sx q[1];
rz(-1.5624996) q[1];
sx q[1];
rz(-1.0992728) q[1];
rz(-pi) q[2];
rz(1.7532888) q[3];
sx q[3];
rz(-1.914822) q[3];
sx q[3];
rz(1.9545912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36416546) q[2];
sx q[2];
rz(-1.825793) q[2];
sx q[2];
rz(2.8167456) q[2];
rz(0.15803629) q[3];
sx q[3];
rz(-2.7285748) q[3];
sx q[3];
rz(2.1260927) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3196816) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(0.88777375) q[0];
rz(2.7582788) q[1];
sx q[1];
rz(-0.8822459) q[1];
sx q[1];
rz(0.15957889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68113806) q[0];
sx q[0];
rz(-0.41550203) q[0];
sx q[0];
rz(2.2592708) q[0];
rz(-0.80990661) q[2];
sx q[2];
rz(-2.4121373) q[2];
sx q[2];
rz(-0.73939656) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9099906) q[1];
sx q[1];
rz(-1.195915) q[1];
sx q[1];
rz(-1.7085307) q[1];
x q[2];
rz(-1.4042312) q[3];
sx q[3];
rz(-1.8565531) q[3];
sx q[3];
rz(-0.33657246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5923656) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(-1.7249736) q[2];
rz(1.8727411) q[3];
sx q[3];
rz(-2.3609991) q[3];
sx q[3];
rz(0.95808539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6855327) q[0];
sx q[0];
rz(-1.1153509) q[0];
sx q[0];
rz(0.90150315) q[0];
rz(2.447336) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(1.136397) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0471056) q[0];
sx q[0];
rz(-1.7312972) q[0];
sx q[0];
rz(-3.1094527) q[0];
rz(-pi) q[1];
rz(-2.5752399) q[2];
sx q[2];
rz(-2.512062) q[2];
sx q[2];
rz(2.9734263) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1727089) q[1];
sx q[1];
rz(-0.69689059) q[1];
sx q[1];
rz(1.642176) q[1];
rz(-pi) q[2];
rz(-1.0497007) q[3];
sx q[3];
rz(-2.2302094) q[3];
sx q[3];
rz(-1.5707317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14094341) q[2];
sx q[2];
rz(-1.2691754) q[2];
sx q[2];
rz(0.81306523) q[2];
rz(2.5874169) q[3];
sx q[3];
rz(-1.412609) q[3];
sx q[3];
rz(0.36177844) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56162214) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(-2.9680874) q[0];
rz(-0.42501998) q[1];
sx q[1];
rz(-1.605502) q[1];
sx q[1];
rz(-0.82957155) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94151992) q[0];
sx q[0];
rz(-0.74848142) q[0];
sx q[0];
rz(-2.9274288) q[0];
rz(3.095171) q[2];
sx q[2];
rz(-1.6116465) q[2];
sx q[2];
rz(-2.61039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8977938) q[1];
sx q[1];
rz(-1.3659371) q[1];
sx q[1];
rz(-2.5216991) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18316571) q[3];
sx q[3];
rz(-0.34160994) q[3];
sx q[3];
rz(0.46942018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7912264) q[2];
sx q[2];
rz(-1.3000877) q[2];
sx q[2];
rz(1.4607956) q[2];
rz(0.99460498) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(-2.2554956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643672) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(-0.10203578) q[0];
rz(0.61965865) q[1];
sx q[1];
rz(-1.6435868) q[1];
sx q[1];
rz(0.59725753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0181018) q[0];
sx q[0];
rz(-2.1649556) q[0];
sx q[0];
rz(1.4090562) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37924699) q[2];
sx q[2];
rz(-2.7617117) q[2];
sx q[2];
rz(0.82373649) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1828536) q[1];
sx q[1];
rz(-2.5748061) q[1];
sx q[1];
rz(1.502152) q[1];
rz(-pi) q[2];
rz(0.47145505) q[3];
sx q[3];
rz(-1.5999891) q[3];
sx q[3];
rz(1.4271133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20798802) q[2];
sx q[2];
rz(-0.36078578) q[2];
sx q[2];
rz(-2.2130845) q[2];
rz(1.1601296) q[3];
sx q[3];
rz(-1.7103633) q[3];
sx q[3];
rz(2.3742356) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930897) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(2.0139991) q[1];
sx q[1];
rz(-0.53378202) q[1];
sx q[1];
rz(3.1376874) q[1];
rz(-2.247159) q[2];
sx q[2];
rz(-1.5147687) q[2];
sx q[2];
rz(1.1942836) q[2];
rz(0.34251681) q[3];
sx q[3];
rz(-2.2377662) q[3];
sx q[3];
rz(-0.025442414) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
