OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(0.55895609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6214949) q[0];
sx q[0];
rz(-1.1557475) q[0];
sx q[0];
rz(-2.9893304) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5833601) q[2];
sx q[2];
rz(-1.0135092) q[2];
sx q[2];
rz(-1.3249427) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.641687) q[1];
sx q[1];
rz(-0.99398621) q[1];
sx q[1];
rz(-2.4813502) q[1];
x q[2];
rz(0.55922237) q[3];
sx q[3];
rz(-2.6359574) q[3];
sx q[3];
rz(1.7749744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3922334) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(2.5048845) q[2];
rz(0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(-0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-0.15287457) q[0];
rz(2.3846467) q[1];
sx q[1];
rz(-1.5870973) q[1];
sx q[1];
rz(-2.1551932) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61030412) q[0];
sx q[0];
rz(-0.072040759) q[0];
sx q[0];
rz(0.59285592) q[0];
rz(-pi) q[1];
rz(2.0979004) q[2];
sx q[2];
rz(-0.88000789) q[2];
sx q[2];
rz(-1.8156798) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0806597) q[1];
sx q[1];
rz(-2.6903209) q[1];
sx q[1];
rz(-2.5732451) q[1];
rz(-pi) q[2];
rz(1.2419224) q[3];
sx q[3];
rz(-1.4423014) q[3];
sx q[3];
rz(-1.0458898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5622921) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(-2.3584649) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(2.7338681) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(-0.12869421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0823682) q[0];
sx q[0];
rz(-1.6186065) q[0];
sx q[0];
rz(1.6533018) q[0];
x q[1];
rz(2.6016597) q[2];
sx q[2];
rz(-0.88368249) q[2];
sx q[2];
rz(-0.5772669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8752746) q[1];
sx q[1];
rz(-2.6058063) q[1];
sx q[1];
rz(2.5712625) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7101173) q[3];
sx q[3];
rz(-0.77292597) q[3];
sx q[3];
rz(1.3611925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.241768) q[2];
rz(-0.5870108) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(-2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509165) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(2.0571016) q[0];
rz(1.4831316) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(3.049057) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0412484) q[0];
sx q[0];
rz(-0.39533246) q[0];
sx q[0];
rz(2.3325217) q[0];
x q[1];
rz(2.1763457) q[2];
sx q[2];
rz(-2.5778228) q[2];
sx q[2];
rz(-0.077686003) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7177928) q[1];
sx q[1];
rz(-1.6289627) q[1];
sx q[1];
rz(-0.34323005) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2257005) q[3];
sx q[3];
rz(-1.4426784) q[3];
sx q[3];
rz(-1.7803943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83355054) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(-2.8095424) q[2];
rz(-2.0856693) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(-0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-2.9300368) q[0];
rz(-1.3062723) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(0.64770118) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6015198) q[0];
sx q[0];
rz(-1.411502) q[0];
sx q[0];
rz(1.457731) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7251882) q[2];
sx q[2];
rz(-1.9152181) q[2];
sx q[2];
rz(-1.4167348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5375145) q[1];
sx q[1];
rz(-1.068183) q[1];
sx q[1];
rz(-0.22488774) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5557489) q[3];
sx q[3];
rz(-1.589937) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5806879) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(2.4723315) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(-3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-2.9142002) q[0];
sx q[0];
rz(1.2325226) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-2.9673064) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5151383) q[0];
sx q[0];
rz(-0.8695375) q[0];
sx q[0];
rz(0.40909543) q[0];
rz(-pi) q[1];
rz(-0.8308438) q[2];
sx q[2];
rz(-0.85692642) q[2];
sx q[2];
rz(-1.4884399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3534527) q[1];
sx q[1];
rz(-2.3056245) q[1];
sx q[1];
rz(-1.047903) q[1];
x q[2];
rz(0.11630451) q[3];
sx q[3];
rz(-1.8810086) q[3];
sx q[3];
rz(0.54799622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(-0.19763395) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0777247) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(2.1754307) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(-1.5055515) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9150328) q[0];
sx q[0];
rz(-2.0953) q[0];
sx q[0];
rz(2.4302308) q[0];
rz(-1.6072636) q[2];
sx q[2];
rz(-2.0009396) q[2];
sx q[2];
rz(1.210481) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.694043) q[1];
sx q[1];
rz(-0.72206891) q[1];
sx q[1];
rz(0.85192792) q[1];
rz(-pi) q[2];
rz(-2.8553477) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(0.092982987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.85764) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(-1.3939259) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(-2.6275997) q[0];
rz(-3.0184074) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(2.2095912) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31873736) q[0];
sx q[0];
rz(-2.0120976) q[0];
sx q[0];
rz(-0.1499099) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9281689) q[2];
sx q[2];
rz(-0.66713152) q[2];
sx q[2];
rz(2.2156029) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92703687) q[1];
sx q[1];
rz(-2.3408457) q[1];
sx q[1];
rz(-2.3826249) q[1];
rz(-pi) q[2];
rz(2.2873015) q[3];
sx q[3];
rz(-1.4201418) q[3];
sx q[3];
rz(-0.52779576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1121858) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(0.4294447) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(-1.3388348) q[0];
rz(-0.70612899) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50981748) q[0];
sx q[0];
rz(-2.0777367) q[0];
sx q[0];
rz(0.99960534) q[0];
x q[1];
rz(-1.9865932) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(-1.6378251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1817158) q[1];
sx q[1];
rz(-2.1023589) q[1];
sx q[1];
rz(-2.9243484) q[1];
x q[2];
rz(-1.9825963) q[3];
sx q[3];
rz(-1.3367062) q[3];
sx q[3];
rz(-1.7115418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(-2.2144923) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(-2.5203729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973688) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(-2.9123059) q[0];
rz(2.7067822) q[1];
sx q[1];
rz(-1.9186585) q[1];
sx q[1];
rz(2.4226709) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7867891) q[0];
sx q[0];
rz(-2.1956586) q[0];
sx q[0];
rz(0.11632365) q[0];
rz(0.98307307) q[2];
sx q[2];
rz(-1.6694153) q[2];
sx q[2];
rz(-1.3999511) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6112411) q[1];
sx q[1];
rz(-1.9396922) q[1];
sx q[1];
rz(1.7289274) q[1];
x q[2];
rz(2.2273916) q[3];
sx q[3];
rz(-0.73540348) q[3];
sx q[3];
rz(-0.38287336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3832613) q[2];
sx q[2];
rz(-1.2021844) q[2];
sx q[2];
rz(0.74679217) q[2];
rz(0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(-2.1993568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9185716) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(0.37721286) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-0.51858356) q[2];
sx q[2];
rz(-1.2972144) q[2];
sx q[2];
rz(-1.5757061) q[2];
rz(1.7198346) q[3];
sx q[3];
rz(-2.2825713) q[3];
sx q[3];
rz(1.5406516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
