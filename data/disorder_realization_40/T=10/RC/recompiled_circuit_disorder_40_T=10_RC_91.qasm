OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(-2.9878374) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35225866) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(2.5369011) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2674238) q[2];
sx q[2];
rz(-0.25250013) q[2];
sx q[2];
rz(1.7054103) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54290463) q[1];
sx q[1];
rz(-2.821327) q[1];
sx q[1];
rz(1.4742875) q[1];
x q[2];
rz(-0.93011232) q[3];
sx q[3];
rz(-0.53033391) q[3];
sx q[3];
rz(1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.443632) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(-1.704818) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(-0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.8181713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.130587) q[0];
sx q[0];
rz(-0.6479833) q[0];
sx q[0];
rz(-0.76052163) q[0];
rz(-pi) q[1];
rz(-3.0201868) q[2];
sx q[2];
rz(-0.3519904) q[2];
sx q[2];
rz(0.44167232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18113187) q[1];
sx q[1];
rz(-2.301553) q[1];
sx q[1];
rz(-1.9366656) q[1];
rz(0.99625846) q[3];
sx q[3];
rz(-1.2401476) q[3];
sx q[3];
rz(-1.4853256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20415846) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(2.5126863) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.73308289) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-2.4531903) q[0];
rz(3.0738661) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-2.6115131) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63307525) q[0];
sx q[0];
rz(-1.8607664) q[0];
sx q[0];
rz(0.11147186) q[0];
x q[1];
rz(-2.7128503) q[2];
sx q[2];
rz(-1.0849761) q[2];
sx q[2];
rz(-2.5634917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9024132) q[1];
sx q[1];
rz(-1.9063519) q[1];
sx q[1];
rz(-2.2526342) q[1];
rz(-1.3313053) q[3];
sx q[3];
rz(-1.665984) q[3];
sx q[3];
rz(-0.59323192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80785859) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(0.90144908) q[2];
rz(-0.83550134) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(-0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(-0.13664666) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89350677) q[0];
sx q[0];
rz(-2.6419905) q[0];
sx q[0];
rz(-1.8462371) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58089528) q[2];
sx q[2];
rz(-0.53933203) q[2];
sx q[2];
rz(0.5667333) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1130484) q[1];
sx q[1];
rz(-1.5469157) q[1];
sx q[1];
rz(-1.5024115) q[1];
rz(-pi) q[2];
rz(2.1468048) q[3];
sx q[3];
rz(-1.4928865) q[3];
sx q[3];
rz(-2.5324477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(2.5811035) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(-2.3204455) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.7339773) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27027425) q[0];
sx q[0];
rz(-2.6042013) q[0];
sx q[0];
rz(-2.4094765) q[0];
rz(-0.9887092) q[2];
sx q[2];
rz(-1.313813) q[2];
sx q[2];
rz(2.7866521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.447532) q[1];
sx q[1];
rz(-2.5563055) q[1];
sx q[1];
rz(1.3259757) q[1];
rz(-pi) q[2];
rz(2.9191454) q[3];
sx q[3];
rz(-1.9223833) q[3];
sx q[3];
rz(0.44140154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(0.042479854) q[2];
rz(0.63043198) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(0.49155864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38934389) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(2.65843) q[0];
rz(-2.0893611) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(2.5767456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0739581) q[0];
sx q[0];
rz(-1.7059776) q[0];
sx q[0];
rz(1.0413175) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1499314) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(-1.7771306) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9996623) q[1];
sx q[1];
rz(-1.511444) q[1];
sx q[1];
rz(-1.1430986) q[1];
x q[2];
rz(2.9103043) q[3];
sx q[3];
rz(-1.7219208) q[3];
sx q[3];
rz(-0.95201991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(-1.8035536) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17334443) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(2.3563747) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-0.26842591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82314202) q[0];
sx q[0];
rz(-0.25932352) q[0];
sx q[0];
rz(0.32487049) q[0];
x q[1];
rz(0.7873017) q[2];
sx q[2];
rz(-2.1848218) q[2];
sx q[2];
rz(-1.3348483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5749579) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(0.11771867) q[1];
rz(-pi) q[2];
rz(1.0915756) q[3];
sx q[3];
rz(-0.5939393) q[3];
sx q[3];
rz(-1.493243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(2.3274373) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(-0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(-2.1222173) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-0.44874915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7166324) q[0];
sx q[0];
rz(-2.6052931) q[0];
sx q[0];
rz(1.1015571) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2737695) q[2];
sx q[2];
rz(-0.38735577) q[2];
sx q[2];
rz(2.5572436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4510348) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(1.8655538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71172165) q[3];
sx q[3];
rz(-2.2300643) q[3];
sx q[3];
rz(-0.014051138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2723508) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-2.8273919) q[2];
rz(2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(-1.8810133) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2429598) q[0];
sx q[0];
rz(-1.3071994) q[0];
sx q[0];
rz(-1.5893057) q[0];
rz(-1.47255) q[2];
sx q[2];
rz(-1.7260523) q[2];
sx q[2];
rz(-1.1467903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70871204) q[1];
sx q[1];
rz(-1.6212665) q[1];
sx q[1];
rz(-2.9570079) q[1];
rz(0.28963611) q[3];
sx q[3];
rz(-1.3445065) q[3];
sx q[3];
rz(-0.7520523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(0.67031676) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(1.5746501) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(2.0589028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2912746) q[0];
sx q[0];
rz(-0.58915888) q[0];
sx q[0];
rz(0.088081443) q[0];
rz(-pi) q[1];
rz(-0.81929368) q[2];
sx q[2];
rz(-1.5466006) q[2];
sx q[2];
rz(1.0011315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.078552695) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(0.54363721) q[1];
rz(-pi) q[2];
rz(1.2438251) q[3];
sx q[3];
rz(-1.9332814) q[3];
sx q[3];
rz(2.8231951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7426976) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(0.51188525) q[2];
rz(-2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-0.40047405) q[2];
sx q[2];
rz(-2.1486712) q[2];
sx q[2];
rz(1.490544) q[2];
rz(1.742733) q[3];
sx q[3];
rz(-2.3098683) q[3];
sx q[3];
rz(1.5666425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
