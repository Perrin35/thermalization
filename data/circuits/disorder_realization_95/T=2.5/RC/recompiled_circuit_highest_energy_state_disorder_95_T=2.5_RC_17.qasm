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
rz(-2.4565826) q[0];
sx q[0];
rz(-0.69184875) q[0];
sx q[0];
rz(-0.43842167) q[0];
rz(3.0216079) q[1];
sx q[1];
rz(3.383145) q[1];
sx q[1];
rz(8.4278843) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0783421) q[0];
sx q[0];
rz(-1.7228427) q[0];
sx q[0];
rz(-3.0536821) q[0];
rz(-pi) q[1];
rz(1.5307776) q[2];
sx q[2];
rz(-1.8410896) q[2];
sx q[2];
rz(1.1129462) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4516075) q[1];
sx q[1];
rz(-1.3295191) q[1];
sx q[1];
rz(-1.9764158) q[1];
x q[2];
rz(2.4313967) q[3];
sx q[3];
rz(-0.3085779) q[3];
sx q[3];
rz(-0.86581826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3991656) q[2];
sx q[2];
rz(-0.70968598) q[2];
sx q[2];
rz(-0.7134552) q[2];
rz(-0.82792884) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(0.92638612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68726081) q[0];
sx q[0];
rz(-0.61606544) q[0];
sx q[0];
rz(1.8808421) q[0];
rz(-2.8997391) q[1];
sx q[1];
rz(-1.2956023) q[1];
sx q[1];
rz(0.58580011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52710746) q[0];
sx q[0];
rz(-1.7298843) q[0];
sx q[0];
rz(-1.2522526) q[0];
rz(1.5620147) q[2];
sx q[2];
rz(-1.9987371) q[2];
sx q[2];
rz(-1.2195171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7603078) q[1];
sx q[1];
rz(-0.80712748) q[1];
sx q[1];
rz(1.9859846) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5473789) q[3];
sx q[3];
rz(-3.1298198) q[3];
sx q[3];
rz(0.6977607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62023097) q[2];
sx q[2];
rz(-2.5197881) q[2];
sx q[2];
rz(2.8561031) q[2];
rz(0.96389687) q[3];
sx q[3];
rz(-2.4745092) q[3];
sx q[3];
rz(-2.9662761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46850878) q[0];
sx q[0];
rz(-2.5991169) q[0];
sx q[0];
rz(2.8520404) q[0];
rz(-0.30329224) q[1];
sx q[1];
rz(-1.2723609) q[1];
sx q[1];
rz(-1.9151275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6942816) q[0];
sx q[0];
rz(-2.3318365) q[0];
sx q[0];
rz(0.96902121) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58019029) q[2];
sx q[2];
rz(-1.3610494) q[2];
sx q[2];
rz(-0.25598393) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7715622) q[1];
sx q[1];
rz(-0.65415934) q[1];
sx q[1];
rz(-1.0942208) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3555894) q[3];
sx q[3];
rz(-1.6542098) q[3];
sx q[3];
rz(-1.0939897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.44846416) q[2];
sx q[2];
rz(-2.144564) q[2];
sx q[2];
rz(0.86359751) q[2];
rz(0.10330769) q[3];
sx q[3];
rz(-1.7938675) q[3];
sx q[3];
rz(-3.0062655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.1358262) q[0];
sx q[0];
rz(-2.5209881) q[0];
sx q[0];
rz(-0.045106877) q[0];
rz(2.3182484) q[1];
sx q[1];
rz(-0.38213676) q[1];
sx q[1];
rz(2.8270922) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053581694) q[0];
sx q[0];
rz(-1.4865666) q[0];
sx q[0];
rz(0.65902577) q[0];
rz(-pi) q[1];
x q[1];
rz(1.034819) q[2];
sx q[2];
rz(-1.2064486) q[2];
sx q[2];
rz(1.8672158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6291677) q[1];
sx q[1];
rz(-1.3864497) q[1];
sx q[1];
rz(2.6434628) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3568595) q[3];
sx q[3];
rz(-1.2535411) q[3];
sx q[3];
rz(-0.91323131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0663674) q[2];
sx q[2];
rz(-1.4868569) q[2];
sx q[2];
rz(-0.60869795) q[2];
rz(1.6913951) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(0.47877413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64959127) q[0];
sx q[0];
rz(-2.2719125) q[0];
sx q[0];
rz(3.0539404) q[0];
rz(1.8805257) q[1];
sx q[1];
rz(-1.7365716) q[1];
sx q[1];
rz(0.87439775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5090655) q[0];
sx q[0];
rz(-1.6035) q[0];
sx q[0];
rz(0.89941179) q[0];
rz(-1.6220785) q[2];
sx q[2];
rz(-0.59133321) q[2];
sx q[2];
rz(1.6720477) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7139018) q[1];
sx q[1];
rz(-2.1449614) q[1];
sx q[1];
rz(2.4533837) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1176862) q[3];
sx q[3];
rz(-1.7611758) q[3];
sx q[3];
rz(-0.067401907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7586691) q[2];
sx q[2];
rz(-1.0367959) q[2];
sx q[2];
rz(-0.60302889) q[2];
rz(-0.99271071) q[3];
sx q[3];
rz(-0.38645667) q[3];
sx q[3];
rz(-2.4934798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822766) q[0];
sx q[0];
rz(-0.74463212) q[0];
sx q[0];
rz(-0.69389206) q[0];
rz(-2.4248185) q[1];
sx q[1];
rz(-2.0718772) q[1];
sx q[1];
rz(2.1778291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474182) q[0];
sx q[0];
rz(-2.6892585) q[0];
sx q[0];
rz(-0.27249713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76002272) q[2];
sx q[2];
rz(-2.8926635) q[2];
sx q[2];
rz(1.6955693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6543944) q[1];
sx q[1];
rz(-0.66283145) q[1];
sx q[1];
rz(2.2064232) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5429872) q[3];
sx q[3];
rz(-2.1429792) q[3];
sx q[3];
rz(-1.5554626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11513772) q[2];
sx q[2];
rz(-2.852735) q[2];
sx q[2];
rz(3.0333983) q[2];
rz(-0.9555971) q[3];
sx q[3];
rz(-0.038766131) q[3];
sx q[3];
rz(2.5819216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3947802) q[0];
sx q[0];
rz(-0.17212269) q[0];
sx q[0];
rz(-0.10699233) q[0];
rz(2.6394898) q[1];
sx q[1];
rz(-1.6306449) q[1];
sx q[1];
rz(-0.14920251) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7764588) q[0];
sx q[0];
rz(-1.7452876) q[0];
sx q[0];
rz(-2.9883238) q[0];
rz(-pi) q[1];
rz(2.4154941) q[2];
sx q[2];
rz(-2.4319785) q[2];
sx q[2];
rz(0.47898705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8684959) q[1];
sx q[1];
rz(-0.78487108) q[1];
sx q[1];
rz(2.0071061) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7140237) q[3];
sx q[3];
rz(-1.9960072) q[3];
sx q[3];
rz(2.3452206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2627829) q[2];
sx q[2];
rz(-2.171319) q[2];
sx q[2];
rz(-2.523017) q[2];
rz(-0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(0.98606199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743643) q[0];
sx q[0];
rz(-1.769913) q[0];
sx q[0];
rz(0.57269639) q[0];
rz(1.0193846) q[1];
sx q[1];
rz(-0.23713325) q[1];
sx q[1];
rz(3.0651029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.50216) q[0];
sx q[0];
rz(-2.0973839) q[0];
sx q[0];
rz(0.16178417) q[0];
x q[1];
rz(-2.8440153) q[2];
sx q[2];
rz(-1.3598765) q[2];
sx q[2];
rz(-0.95266137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.44750252) q[1];
sx q[1];
rz(-1.4610485) q[1];
sx q[1];
rz(0.90417273) q[1];
rz(-0.42726699) q[3];
sx q[3];
rz(-0.92822853) q[3];
sx q[3];
rz(1.6258607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0387592) q[2];
sx q[2];
rz(-0.88492727) q[2];
sx q[2];
rz(2.6494675) q[2];
rz(0.37880185) q[3];
sx q[3];
rz(-0.4327966) q[3];
sx q[3];
rz(-0.88703275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73961306) q[0];
sx q[0];
rz(-2.6519863) q[0];
sx q[0];
rz(2.711645) q[0];
rz(2.6509189) q[1];
sx q[1];
rz(-0.48777598) q[1];
sx q[1];
rz(2.1388163) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96129629) q[0];
sx q[0];
rz(-1.3081828) q[0];
sx q[0];
rz(-1.0680593) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0645686) q[2];
sx q[2];
rz(-0.90156889) q[2];
sx q[2];
rz(-2.2716244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32200798) q[1];
sx q[1];
rz(-2.036398) q[1];
sx q[1];
rz(1.5172537) q[1];
x q[2];
rz(2.0315038) q[3];
sx q[3];
rz(-0.92383251) q[3];
sx q[3];
rz(-1.6443192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7542725) q[2];
sx q[2];
rz(-2.5405799) q[2];
sx q[2];
rz(0.69166541) q[2];
rz(3.0129041) q[3];
sx q[3];
rz(-1.5920937) q[3];
sx q[3];
rz(-2.5531829) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16828951) q[0];
sx q[0];
rz(-3.0423218) q[0];
sx q[0];
rz(0.6231935) q[0];
rz(-0.19459952) q[1];
sx q[1];
rz(-1.1293026) q[1];
sx q[1];
rz(-0.87619877) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19080432) q[0];
sx q[0];
rz(-1.6144469) q[0];
sx q[0];
rz(-1.6859289) q[0];
rz(-pi) q[1];
rz(0.16531971) q[2];
sx q[2];
rz(-0.7484127) q[2];
sx q[2];
rz(0.95647631) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1764123) q[1];
sx q[1];
rz(-1.1438055) q[1];
sx q[1];
rz(-1.4216485) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2262754) q[3];
sx q[3];
rz(-2.1117941) q[3];
sx q[3];
rz(-1.2834383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1759922) q[2];
sx q[2];
rz(-2.8922562) q[2];
sx q[2];
rz(2.5052137) q[2];
rz(-1.5198358) q[3];
sx q[3];
rz(-2.2607925) q[3];
sx q[3];
rz(0.22708587) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31430055) q[0];
sx q[0];
rz(-1.5416332) q[0];
sx q[0];
rz(2.736349) q[0];
rz(-1.7397407) q[1];
sx q[1];
rz(-1.4973462) q[1];
sx q[1];
rz(-3.0037465) q[1];
rz(-0.17086239) q[2];
sx q[2];
rz(-0.62913412) q[2];
sx q[2];
rz(2.4133336) q[2];
rz(-1.6062395) q[3];
sx q[3];
rz(-2.1066022) q[3];
sx q[3];
rz(0.71577241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
