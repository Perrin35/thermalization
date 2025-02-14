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
rz(-1.350116) q[0];
sx q[0];
rz(3.8173563) q[0];
sx q[0];
rz(9.4326333) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(1.6176728) q[1];
sx q[1];
rz(9.6674506) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8382671) q[0];
sx q[0];
rz(-1.4446065) q[0];
sx q[0];
rz(2.2827006) q[0];
rz(-0.65671667) q[2];
sx q[2];
rz(-0.55865951) q[2];
sx q[2];
rz(-2.5303773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2897799) q[1];
sx q[1];
rz(-0.35408005) q[1];
sx q[1];
rz(0.068239958) q[1];
rz(-pi) q[2];
rz(-2.8318066) q[3];
sx q[3];
rz(-2.4521146) q[3];
sx q[3];
rz(-0.30287095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0365888) q[2];
sx q[2];
rz(-1.9052817) q[2];
sx q[2];
rz(-0.71199065) q[2];
rz(-2.8365734) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(-3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63118339) q[0];
sx q[0];
rz(-1.2774066) q[0];
sx q[0];
rz(1.7741868) q[0];
rz(-0.039904682) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(2.5423999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4112873) q[0];
sx q[0];
rz(-1.4559901) q[0];
sx q[0];
rz(-0.47700096) q[0];
rz(-pi) q[1];
rz(0.56053253) q[2];
sx q[2];
rz(-2.2403702) q[2];
sx q[2];
rz(-2.5472484) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42223052) q[1];
sx q[1];
rz(-1.2368725) q[1];
sx q[1];
rz(0.67113282) q[1];
rz(-pi) q[2];
rz(-2.8017524) q[3];
sx q[3];
rz(-0.71469864) q[3];
sx q[3];
rz(3.0106737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0590608) q[2];
sx q[2];
rz(-0.56190562) q[2];
sx q[2];
rz(-1.3085636) q[2];
rz(-3.1176873) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(-1.5786952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(1.6858653) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(0.84079963) q[0];
rz(1.0288382) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.4604481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2551291) q[0];
sx q[0];
rz(-0.48454075) q[0];
sx q[0];
rz(-2.8302829) q[0];
rz(-pi) q[1];
rz(1.9424136) q[2];
sx q[2];
rz(-1.1254246) q[2];
sx q[2];
rz(0.98486131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9509995) q[1];
sx q[1];
rz(-1.237932) q[1];
sx q[1];
rz(-0.69575633) q[1];
rz(-1.6275241) q[3];
sx q[3];
rz(-0.9642082) q[3];
sx q[3];
rz(-1.1067672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8847522) q[2];
sx q[2];
rz(-1.2105056) q[2];
sx q[2];
rz(0.15667008) q[2];
rz(-1.5001851) q[3];
sx q[3];
rz(-1.898396) q[3];
sx q[3];
rz(0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.0665862) q[0];
sx q[0];
rz(-1.2970507) q[0];
sx q[0];
rz(0.080667607) q[0];
rz(-2.0196041) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(-1.5871619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98469668) q[0];
sx q[0];
rz(-1.5050355) q[0];
sx q[0];
rz(0.82729152) q[0];
rz(-pi) q[1];
rz(-2.8055259) q[2];
sx q[2];
rz(-3.0185351) q[2];
sx q[2];
rz(2.6390136) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39496468) q[1];
sx q[1];
rz(-2.0680799) q[1];
sx q[1];
rz(-2.4056466) q[1];
rz(-pi) q[2];
rz(2.5347363) q[3];
sx q[3];
rz(-1.292406) q[3];
sx q[3];
rz(2.3902219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88177219) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(2.1176977) q[2];
rz(2.3648868) q[3];
sx q[3];
rz(-1.1506162) q[3];
sx q[3];
rz(-1.0030494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40989947) q[0];
sx q[0];
rz(-0.85407805) q[0];
sx q[0];
rz(2.5140629) q[0];
rz(-0.69951406) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(0.90739179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8555657) q[0];
sx q[0];
rz(-0.81840179) q[0];
sx q[0];
rz(1.6256385) q[0];
x q[1];
rz(-2.9512915) q[2];
sx q[2];
rz(-1.7277328) q[2];
sx q[2];
rz(-2.4630594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34037922) q[1];
sx q[1];
rz(-2.0836909) q[1];
sx q[1];
rz(-2.7199634) q[1];
x q[2];
rz(-1.6697407) q[3];
sx q[3];
rz(-2.6032748) q[3];
sx q[3];
rz(0.42699277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13713914) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(2.2966906) q[2];
rz(-0.6984624) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(2.3499878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38248211) q[0];
sx q[0];
rz(-1.4610721) q[0];
sx q[0];
rz(-1.8735029) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(-0.39438927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.012082) q[0];
sx q[0];
rz(-0.92716588) q[0];
sx q[0];
rz(-0.056931007) q[0];
x q[1];
rz(2.4066448) q[2];
sx q[2];
rz(-1.6326666) q[2];
sx q[2];
rz(0.29145539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39526597) q[1];
sx q[1];
rz(-1.6155287) q[1];
sx q[1];
rz(-0.6529863) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9378618) q[3];
sx q[3];
rz(-0.36760783) q[3];
sx q[3];
rz(2.550761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4322728) q[2];
sx q[2];
rz(-0.061881438) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(-2.0105441) q[3];
sx q[3];
rz(-2.3384422) q[3];
sx q[3];
rz(-2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442218) q[0];
sx q[0];
rz(-0.18246305) q[0];
sx q[0];
rz(-0.23319787) q[0];
rz(-3.0274262) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(0.17359576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19464201) q[0];
sx q[0];
rz(-1.9341) q[0];
sx q[0];
rz(-0.86473303) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2701413) q[2];
sx q[2];
rz(-2.5837499) q[2];
sx q[2];
rz(-2.7047472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1039338) q[1];
sx q[1];
rz(-1.3500431) q[1];
sx q[1];
rz(-2.9751865) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9529157) q[3];
sx q[3];
rz(-1.5828307) q[3];
sx q[3];
rz(0.4922315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1058098) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(0.85476533) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.1707183) q[3];
sx q[3];
rz(0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6524803) q[0];
sx q[0];
rz(-1.2615477) q[0];
sx q[0];
rz(-1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.2636355) q[1];
sx q[1];
rz(2.9875535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48705593) q[0];
sx q[0];
rz(-2.3628919) q[0];
sx q[0];
rz(0.47327431) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70239046) q[2];
sx q[2];
rz(-0.99961126) q[2];
sx q[2];
rz(-1.0016425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8785868) q[1];
sx q[1];
rz(-2.03555) q[1];
sx q[1];
rz(2.8057363) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4064111) q[3];
sx q[3];
rz(-0.38215853) q[3];
sx q[3];
rz(-0.89565403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68086326) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(-0.66780773) q[2];
rz(-1.4051416) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(2.5051129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3619096) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(2.5478126) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.8437754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058859874) q[0];
sx q[0];
rz(-1.4987336) q[0];
sx q[0];
rz(2.5565992) q[0];
x q[1];
rz(1.2601398) q[2];
sx q[2];
rz(-0.96608487) q[2];
sx q[2];
rz(-2.8653646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5169277) q[1];
sx q[1];
rz(-1.4447277) q[1];
sx q[1];
rz(0.36051118) q[1];
rz(-pi) q[2];
rz(2.3343349) q[3];
sx q[3];
rz(-1.689925) q[3];
sx q[3];
rz(-1.8594683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7456776) q[2];
sx q[2];
rz(-0.021947689) q[2];
sx q[2];
rz(-1.6321261) q[2];
rz(3.0294026) q[3];
sx q[3];
rz(-2.1144919) q[3];
sx q[3];
rz(-1.7621015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6701732) q[0];
sx q[0];
rz(-2.1355974) q[0];
sx q[0];
rz(0.26563409) q[0];
rz(0.98948014) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(2.8094453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50221793) q[0];
sx q[0];
rz(-1.2373588) q[0];
sx q[0];
rz(3.1375454) q[0];
rz(-pi) q[1];
rz(-2.0004326) q[2];
sx q[2];
rz(-2.5186335) q[2];
sx q[2];
rz(-1.0819544) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0964342) q[1];
sx q[1];
rz(-1.9109042) q[1];
sx q[1];
rz(-0.051326871) q[1];
rz(-pi) q[2];
x q[2];
rz(0.048236851) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(-1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0418479) q[2];
sx q[2];
rz(-1.0724649) q[2];
sx q[2];
rz(2.9912046) q[2];
rz(-2.2930875) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(1.295804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(-1.3311483) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(1.4990357) q[2];
sx q[2];
rz(-0.74081479) q[2];
sx q[2];
rz(-0.011394636) q[2];
rz(3.1114586) q[3];
sx q[3];
rz(-2.4345955) q[3];
sx q[3];
rz(1.101936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
