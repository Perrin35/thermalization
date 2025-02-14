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
rz(-3.0191874) q[0];
sx q[0];
rz(-2.2221017) q[0];
sx q[0];
rz(-1.9424196) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(-2.5993102) q[1];
sx q[1];
rz(1.6220925) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2080948) q[0];
sx q[0];
rz(-1.6483432) q[0];
sx q[0];
rz(-2.5258875) q[0];
rz(0.89251065) q[2];
sx q[2];
rz(-0.55858382) q[2];
sx q[2];
rz(-0.19772274) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.28827661) q[1];
sx q[1];
rz(-1.9717934) q[1];
sx q[1];
rz(-0.53944352) q[1];
x q[2];
rz(0.24418045) q[3];
sx q[3];
rz(-1.5236519) q[3];
sx q[3];
rz(2.3650996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2970994) q[2];
sx q[2];
rz(-2.5145734) q[2];
sx q[2];
rz(-1.7681047) q[2];
rz(2.3376236) q[3];
sx q[3];
rz(-1.4990025) q[3];
sx q[3];
rz(1.1871626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3961769) q[0];
sx q[0];
rz(-2.4276623) q[0];
sx q[0];
rz(0.3748689) q[0];
rz(2.6238341) q[1];
sx q[1];
rz(-1.2034143) q[1];
sx q[1];
rz(1.3688603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8086209) q[0];
sx q[0];
rz(-0.00095168984) q[0];
sx q[0];
rz(-3.0406221) q[0];
rz(-0.13475864) q[2];
sx q[2];
rz(-0.19821168) q[2];
sx q[2];
rz(-1.0466087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6102099) q[1];
sx q[1];
rz(-1.0730337) q[1];
sx q[1];
rz(1.1019635) q[1];
rz(3.0053455) q[3];
sx q[3];
rz(-1.7833424) q[3];
sx q[3];
rz(-3.1177136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0181197) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(2.4269721) q[2];
rz(-2.9361652) q[3];
sx q[3];
rz(-1.8193865) q[3];
sx q[3];
rz(-1.8429168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648723) q[0];
sx q[0];
rz(-1.0370075) q[0];
sx q[0];
rz(-2.6439457) q[0];
rz(-0.63703713) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(3.0348437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15539385) q[0];
sx q[0];
rz(-1.8698911) q[0];
sx q[0];
rz(-0.22217447) q[0];
x q[1];
rz(3.1349772) q[2];
sx q[2];
rz(-2.2405036) q[2];
sx q[2];
rz(-2.1257328) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6124668) q[1];
sx q[1];
rz(-2.4262316) q[1];
sx q[1];
rz(-1.6523408) q[1];
x q[2];
rz(0.62415267) q[3];
sx q[3];
rz(-0.48138967) q[3];
sx q[3];
rz(-2.2311051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5367624) q[2];
sx q[2];
rz(-2.3894775) q[2];
sx q[2];
rz(-2.9467648) q[2];
rz(0.005006494) q[3];
sx q[3];
rz(-2.5275793) q[3];
sx q[3];
rz(-2.3495638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0358589) q[0];
sx q[0];
rz(-2.6073313) q[0];
sx q[0];
rz(-1.0400829) q[0];
rz(-0.19373521) q[1];
sx q[1];
rz(-1.5015142) q[1];
sx q[1];
rz(-0.74660444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0950985) q[0];
sx q[0];
rz(-1.455716) q[0];
sx q[0];
rz(-2.0128573) q[0];
rz(-pi) q[1];
rz(-0.87073054) q[2];
sx q[2];
rz(-1.4615897) q[2];
sx q[2];
rz(2.2094215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0178448) q[1];
sx q[1];
rz(-0.87758076) q[1];
sx q[1];
rz(-1.6381532) q[1];
rz(2.7203619) q[3];
sx q[3];
rz(-0.26184362) q[3];
sx q[3];
rz(-2.1986442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1497583) q[2];
sx q[2];
rz(-1.4703625) q[2];
sx q[2];
rz(2.9370918) q[2];
rz(1.9483942) q[3];
sx q[3];
rz(-0.85719332) q[3];
sx q[3];
rz(2.1804555) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0123154) q[0];
sx q[0];
rz(-2.9670872) q[0];
sx q[0];
rz(2.0611064) q[0];
rz(-2.7642545) q[1];
sx q[1];
rz(-1.6308547) q[1];
sx q[1];
rz(0.89471716) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6166234) q[0];
sx q[0];
rz(-1.3436514) q[0];
sx q[0];
rz(2.7779915) q[0];
rz(2.2316049) q[2];
sx q[2];
rz(-2.2468781) q[2];
sx q[2];
rz(-2.8582339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91352039) q[1];
sx q[1];
rz(-1.0723877) q[1];
sx q[1];
rz(1.0105039) q[1];
x q[2];
rz(2.7758046) q[3];
sx q[3];
rz(-0.96082965) q[3];
sx q[3];
rz(-0.47780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47199029) q[2];
sx q[2];
rz(-2.695638) q[2];
sx q[2];
rz(-3.0338244) q[2];
rz(-2.9660411) q[3];
sx q[3];
rz(-1.2551509) q[3];
sx q[3];
rz(2.5811783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046722978) q[0];
sx q[0];
rz(-0.51114285) q[0];
sx q[0];
rz(-1.0119337) q[0];
rz(-1.5049505) q[1];
sx q[1];
rz(-2.4220146) q[1];
sx q[1];
rz(-2.124427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95421529) q[0];
sx q[0];
rz(-2.7972176) q[0];
sx q[0];
rz(-2.1013174) q[0];
rz(-1.0083197) q[2];
sx q[2];
rz(-2.2485844) q[2];
sx q[2];
rz(0.00046367292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6222657) q[1];
sx q[1];
rz(-2.9510088) q[1];
sx q[1];
rz(-3.0960073) q[1];
x q[2];
rz(2.0200219) q[3];
sx q[3];
rz(-1.7127081) q[3];
sx q[3];
rz(-2.6117539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4096058) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(-2.1448263) q[2];
rz(-0.064920001) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(1.9916649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053059269) q[0];
sx q[0];
rz(-0.94045883) q[0];
sx q[0];
rz(0.9084107) q[0];
rz(-2.9216595) q[1];
sx q[1];
rz(-1.5426153) q[1];
sx q[1];
rz(1.8904846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4804937) q[0];
sx q[0];
rz(-1.1014043) q[0];
sx q[0];
rz(0.29222699) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8562137) q[2];
sx q[2];
rz(-1.6661465) q[2];
sx q[2];
rz(0.91172632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0323769) q[1];
sx q[1];
rz(-0.93689474) q[1];
sx q[1];
rz(-0.79297592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7102107) q[3];
sx q[3];
rz(-2.406139) q[3];
sx q[3];
rz(-0.212634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0327586) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(1.533482) q[2];
rz(-1.9882625) q[3];
sx q[3];
rz(-1.0554375) q[3];
sx q[3];
rz(1.4920894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.15025) q[0];
sx q[0];
rz(-2.7582176) q[0];
sx q[0];
rz(-2.2008994) q[0];
rz(-3.0062145) q[1];
sx q[1];
rz(-1.7372513) q[1];
sx q[1];
rz(2.1167596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12176421) q[0];
sx q[0];
rz(-1.8786977) q[0];
sx q[0];
rz(-2.204049) q[0];
rz(-pi) q[1];
x q[1];
rz(1.631279) q[2];
sx q[2];
rz(-0.96394682) q[2];
sx q[2];
rz(-1.902193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8814446) q[1];
sx q[1];
rz(-1.3760412) q[1];
sx q[1];
rz(1.4347726) q[1];
rz(2.7303194) q[3];
sx q[3];
rz(-1.6182634) q[3];
sx q[3];
rz(0.91225831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3465053) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(-0.65004641) q[2];
rz(-2.5551689) q[3];
sx q[3];
rz(-1.6610049) q[3];
sx q[3];
rz(2.39095) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745678) q[0];
sx q[0];
rz(-2.6331007) q[0];
sx q[0];
rz(1.1453999) q[0];
rz(-1.3151431) q[1];
sx q[1];
rz(-1.2329085) q[1];
sx q[1];
rz(0.68793908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74003873) q[0];
sx q[0];
rz(-1.382022) q[0];
sx q[0];
rz(-1.1925634) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1150949) q[2];
sx q[2];
rz(-2.470068) q[2];
sx q[2];
rz(-1.4948927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2738139) q[1];
sx q[1];
rz(-1.7325337) q[1];
sx q[1];
rz(0.75839569) q[1];
x q[2];
rz(2.7642058) q[3];
sx q[3];
rz(-2.0321369) q[3];
sx q[3];
rz(-0.93022202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73100662) q[2];
sx q[2];
rz(-1.1062016) q[2];
sx q[2];
rz(-0.86622396) q[2];
rz(0.90803641) q[3];
sx q[3];
rz(-0.76483813) q[3];
sx q[3];
rz(-2.0456555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1753801) q[0];
sx q[0];
rz(-1.1840273) q[0];
sx q[0];
rz(-0.41561919) q[0];
rz(-0.79187727) q[1];
sx q[1];
rz(-1.2245347) q[1];
sx q[1];
rz(-0.22722879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4675605) q[0];
sx q[0];
rz(-1.5590384) q[0];
sx q[0];
rz(-0.016076465) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.077543295) q[2];
sx q[2];
rz(-2.3334586) q[2];
sx q[2];
rz(-1.6412954) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7586582) q[1];
sx q[1];
rz(-0.99731748) q[1];
sx q[1];
rz(-1.2632779) q[1];
rz(-pi) q[2];
rz(-2.6260904) q[3];
sx q[3];
rz(-1.1519264) q[3];
sx q[3];
rz(2.2282698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26433429) q[2];
sx q[2];
rz(-2.4206968) q[2];
sx q[2];
rz(0.43992511) q[2];
rz(2.544493) q[3];
sx q[3];
rz(-0.66949451) q[3];
sx q[3];
rz(-1.3574903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6482342) q[0];
sx q[0];
rz(-1.5536722) q[0];
sx q[0];
rz(-1.8552725) q[0];
rz(1.7276806) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(-1.6261423) q[2];
sx q[2];
rz(-1.9774441) q[2];
sx q[2];
rz(-1.0033506) q[2];
rz(-0.8181994) q[3];
sx q[3];
rz(-1.7530734) q[3];
sx q[3];
rz(-1.659041) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
