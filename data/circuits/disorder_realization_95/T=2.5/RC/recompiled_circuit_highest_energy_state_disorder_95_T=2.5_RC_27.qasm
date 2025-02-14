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
rz(0.68501002) q[0];
sx q[0];
rz(-2.4497439) q[0];
sx q[0];
rz(0.43842167) q[0];
rz(3.0216079) q[1];
sx q[1];
rz(-2.9000403) q[1];
sx q[1];
rz(-0.99689364) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50580286) q[0];
sx q[0];
rz(-1.4839026) q[0];
sx q[0];
rz(-1.723423) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27049944) q[2];
sx q[2];
rz(-1.6093614) q[2];
sx q[2];
rz(2.6730516) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5146487) q[1];
sx q[1];
rz(-2.6730826) q[1];
sx q[1];
rz(1.0131939) q[1];
rz(2.4313967) q[3];
sx q[3];
rz(-2.8330148) q[3];
sx q[3];
rz(-2.2757744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74242705) q[2];
sx q[2];
rz(-2.4319067) q[2];
sx q[2];
rz(-2.4281375) q[2];
rz(0.82792884) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(-0.92638612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68726081) q[0];
sx q[0];
rz(-2.5255272) q[0];
sx q[0];
rz(-1.2607505) q[0];
rz(0.24185355) q[1];
sx q[1];
rz(-1.8459903) q[1];
sx q[1];
rz(2.5557925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0958835) q[0];
sx q[0];
rz(-1.2564141) q[0];
sx q[0];
rz(2.9742301) q[0];
rz(-pi) q[1];
rz(-0.42795534) q[2];
sx q[2];
rz(-1.5628067) q[2];
sx q[2];
rz(0.34763476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10639452) q[1];
sx q[1];
rz(-1.8664341) q[1];
sx q[1];
rz(-0.8080478) q[1];
x q[2];
rz(-1.582566) q[3];
sx q[3];
rz(-1.571072) q[3];
sx q[3];
rz(2.2451412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62023097) q[2];
sx q[2];
rz(-0.62180454) q[2];
sx q[2];
rz(0.28548959) q[2];
rz(-2.1776958) q[3];
sx q[3];
rz(-0.6670835) q[3];
sx q[3];
rz(2.9662761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46850878) q[0];
sx q[0];
rz(-0.54247576) q[0];
sx q[0];
rz(-2.8520404) q[0];
rz(-0.30329224) q[1];
sx q[1];
rz(-1.8692317) q[1];
sx q[1];
rz(1.9151275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4604386) q[0];
sx q[0];
rz(-1.1484206) q[0];
sx q[0];
rz(-0.85737164) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3704088) q[2];
sx q[2];
rz(-2.5287712) q[2];
sx q[2];
rz(-1.0072264) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7715622) q[1];
sx q[1];
rz(-0.65415934) q[1];
sx q[1];
rz(-1.0942208) q[1];
rz(-pi) q[2];
rz(0.11762186) q[3];
sx q[3];
rz(-2.3521227) q[3];
sx q[3];
rz(2.5817613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6931285) q[2];
sx q[2];
rz(-0.99702865) q[2];
sx q[2];
rz(-2.2779951) q[2];
rz(3.038285) q[3];
sx q[3];
rz(-1.3477252) q[3];
sx q[3];
rz(0.13532713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1358262) q[0];
sx q[0];
rz(-0.62060452) q[0];
sx q[0];
rz(-0.045106877) q[0];
rz(0.82334423) q[1];
sx q[1];
rz(-2.7594559) q[1];
sx q[1];
rz(2.8270922) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6894486) q[0];
sx q[0];
rz(-2.2270791) q[0];
sx q[0];
rz(-1.4644064) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92934926) q[2];
sx q[2];
rz(-0.63792777) q[2];
sx q[2];
rz(2.8976482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6291677) q[1];
sx q[1];
rz(-1.3864497) q[1];
sx q[1];
rz(-0.49812981) q[1];
rz(-0.32415402) q[3];
sx q[3];
rz(-1.7739033) q[3];
sx q[3];
rz(-0.58989159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0663674) q[2];
sx q[2];
rz(-1.4868569) q[2];
sx q[2];
rz(-0.60869795) q[2];
rz(-1.6913951) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(-0.47877413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4920014) q[0];
sx q[0];
rz(-0.86968017) q[0];
sx q[0];
rz(0.087652303) q[0];
rz(1.2610669) q[1];
sx q[1];
rz(-1.7365716) q[1];
sx q[1];
rz(2.2671949) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1058348) q[0];
sx q[0];
rz(-0.89983655) q[0];
sx q[0];
rz(3.0998328) q[0];
rz(-pi) q[1];
rz(2.1615209) q[2];
sx q[2];
rz(-1.5993759) q[2];
sx q[2];
rz(0.14383741) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4276909) q[1];
sx q[1];
rz(-2.1449614) q[1];
sx q[1];
rz(-2.4533837) q[1];
rz(-pi) q[2];
rz(-2.1176862) q[3];
sx q[3];
rz(-1.7611758) q[3];
sx q[3];
rz(-0.067401907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38292357) q[2];
sx q[2];
rz(-1.0367959) q[2];
sx q[2];
rz(0.60302889) q[2];
rz(2.1488819) q[3];
sx q[3];
rz(-2.755136) q[3];
sx q[3];
rz(2.4934798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25931609) q[0];
sx q[0];
rz(-0.74463212) q[0];
sx q[0];
rz(0.69389206) q[0];
rz(-2.4248185) q[1];
sx q[1];
rz(-1.0697155) q[1];
sx q[1];
rz(-2.1778291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5187274) q[0];
sx q[0];
rz(-1.4528926) q[0];
sx q[0];
rz(-2.7038655) q[0];
rz(2.3815699) q[2];
sx q[2];
rz(-2.8926635) q[2];
sx q[2];
rz(1.4460233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9020667) q[1];
sx q[1];
rz(-1.0527624) q[1];
sx q[1];
rz(0.43398989) q[1];
rz(-pi) q[2];
rz(-3.0984466) q[3];
sx q[3];
rz(-2.5688097) q[3];
sx q[3];
rz(-1.6067895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11513772) q[2];
sx q[2];
rz(-2.852735) q[2];
sx q[2];
rz(-0.10819437) q[2];
rz(2.1859956) q[3];
sx q[3];
rz(-0.038766131) q[3];
sx q[3];
rz(-0.55967104) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3947802) q[0];
sx q[0];
rz(-2.96947) q[0];
sx q[0];
rz(0.10699233) q[0];
rz(2.6394898) q[1];
sx q[1];
rz(-1.5109477) q[1];
sx q[1];
rz(-2.9923901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049292663) q[0];
sx q[0];
rz(-0.23173103) q[0];
sx q[0];
rz(-0.85706237) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57089262) q[2];
sx q[2];
rz(-1.1234267) q[2];
sx q[2];
rz(1.6845861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27309675) q[1];
sx q[1];
rz(-2.3567216) q[1];
sx q[1];
rz(2.0071061) q[1];
x q[2];
rz(1.427569) q[3];
sx q[3];
rz(-1.9960072) q[3];
sx q[3];
rz(2.3452206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2627829) q[2];
sx q[2];
rz(-0.97027367) q[2];
sx q[2];
rz(0.61857569) q[2];
rz(-0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(-2.1555307) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743643) q[0];
sx q[0];
rz(-1.769913) q[0];
sx q[0];
rz(0.57269639) q[0];
rz(2.122208) q[1];
sx q[1];
rz(-0.23713325) q[1];
sx q[1];
rz(-3.0651029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50216) q[0];
sx q[0];
rz(-1.0442088) q[0];
sx q[0];
rz(-2.9798085) q[0];
x q[1];
rz(-2.8440153) q[2];
sx q[2];
rz(-1.7817162) q[2];
sx q[2];
rz(0.95266137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44750252) q[1];
sx q[1];
rz(-1.6805442) q[1];
sx q[1];
rz(2.2374199) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42726699) q[3];
sx q[3];
rz(-2.2133641) q[3];
sx q[3];
rz(-1.515732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1028334) q[2];
sx q[2];
rz(-2.2566654) q[2];
sx q[2];
rz(0.49212512) q[2];
rz(-2.7627908) q[3];
sx q[3];
rz(-2.7087961) q[3];
sx q[3];
rz(-2.2545599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73961306) q[0];
sx q[0];
rz(-2.6519863) q[0];
sx q[0];
rz(-2.711645) q[0];
rz(2.6509189) q[1];
sx q[1];
rz(-2.6538167) q[1];
sx q[1];
rz(1.0027764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96129629) q[0];
sx q[0];
rz(-1.8334098) q[0];
sx q[0];
rz(1.0680593) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6677708) q[2];
sx q[2];
rz(-0.67296689) q[2];
sx q[2];
rz(-2.147858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2247315) q[1];
sx q[1];
rz(-1.6186348) q[1];
sx q[1];
rz(2.6754154) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6096576) q[3];
sx q[3];
rz(-2.3670275) q[3];
sx q[3];
rz(2.3331785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7542725) q[2];
sx q[2];
rz(-0.60101271) q[2];
sx q[2];
rz(0.69166541) q[2];
rz(-3.0129041) q[3];
sx q[3];
rz(-1.5920937) q[3];
sx q[3];
rz(2.5531829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16828951) q[0];
sx q[0];
rz(-0.099270865) q[0];
sx q[0];
rz(2.5183992) q[0];
rz(-2.9469931) q[1];
sx q[1];
rz(-2.01229) q[1];
sx q[1];
rz(-0.87619877) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.122418) q[0];
sx q[0];
rz(-3.0184973) q[0];
sx q[0];
rz(-1.207463) q[0];
x q[1];
rz(-1.4191462) q[2];
sx q[2];
rz(-2.306621) q[2];
sx q[2];
rz(0.73260546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8093578) q[1];
sx q[1];
rz(-1.4351265) q[1];
sx q[1];
rz(-2.7103781) q[1];
rz(-pi) q[2];
rz(-0.51214062) q[3];
sx q[3];
rz(-0.63203963) q[3];
sx q[3];
rz(1.2495981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1759922) q[2];
sx q[2];
rz(-2.8922562) q[2];
sx q[2];
rz(-2.5052137) q[2];
rz(-1.5198358) q[3];
sx q[3];
rz(-2.2607925) q[3];
sx q[3];
rz(0.22708587) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31430055) q[0];
sx q[0];
rz(-1.5999595) q[0];
sx q[0];
rz(-0.40524361) q[0];
rz(1.4018519) q[1];
sx q[1];
rz(-1.4973462) q[1];
sx q[1];
rz(-3.0037465) q[1];
rz(2.9707303) q[2];
sx q[2];
rz(-0.62913412) q[2];
sx q[2];
rz(2.4133336) q[2];
rz(1.6062395) q[3];
sx q[3];
rz(-1.0349904) q[3];
sx q[3];
rz(-2.4258202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
