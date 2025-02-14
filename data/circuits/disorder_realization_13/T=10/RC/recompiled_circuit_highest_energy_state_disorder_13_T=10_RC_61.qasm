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
rz(-2.4969555) q[0];
sx q[0];
rz(-0.58965373) q[0];
sx q[0];
rz(2.0539505) q[0];
rz(3.1261858) q[1];
sx q[1];
rz(-0.3897804) q[1];
sx q[1];
rz(1.1642233) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22237865) q[0];
sx q[0];
rz(-1.5660677) q[0];
sx q[0];
rz(1.4407925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4504775) q[2];
sx q[2];
rz(-1.3045614) q[2];
sx q[2];
rz(-2.63111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1401745) q[1];
sx q[1];
rz(-1.784403) q[1];
sx q[1];
rz(2.1013538) q[1];
rz(-pi) q[2];
rz(2.1415065) q[3];
sx q[3];
rz(-3.0016404) q[3];
sx q[3];
rz(-2.9632448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2414134) q[2];
sx q[2];
rz(-0.56168491) q[2];
sx q[2];
rz(2.4955595) q[2];
rz(0.25520405) q[3];
sx q[3];
rz(-0.45707688) q[3];
sx q[3];
rz(-1.7469762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274886) q[0];
sx q[0];
rz(-0.33559594) q[0];
sx q[0];
rz(0.80701971) q[0];
rz(1.457816) q[1];
sx q[1];
rz(-2.5648263) q[1];
sx q[1];
rz(-1.797765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2349074) q[0];
sx q[0];
rz(-2.8145288) q[0];
sx q[0];
rz(0.48717888) q[0];
x q[1];
rz(2.730663) q[2];
sx q[2];
rz(-2.0996473) q[2];
sx q[2];
rz(-1.5134606) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2652774) q[1];
sx q[1];
rz(-2.1723299) q[1];
sx q[1];
rz(-0.6024036) q[1];
x q[2];
rz(0.044780894) q[3];
sx q[3];
rz(-0.84242994) q[3];
sx q[3];
rz(-1.3842667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1809711) q[2];
sx q[2];
rz(-2.0714859) q[2];
sx q[2];
rz(-0.19372678) q[2];
rz(-1.5242029) q[3];
sx q[3];
rz(-2.3834855) q[3];
sx q[3];
rz(2.962842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5209565) q[0];
sx q[0];
rz(-2.9036324) q[0];
sx q[0];
rz(0.97610193) q[0];
rz(-0.8887662) q[1];
sx q[1];
rz(-0.74405324) q[1];
sx q[1];
rz(0.1730473) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42587468) q[0];
sx q[0];
rz(-3.1039655) q[0];
sx q[0];
rz(-1.4960122) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21645087) q[2];
sx q[2];
rz(-0.5308639) q[2];
sx q[2];
rz(2.5420839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5593619) q[1];
sx q[1];
rz(-1.7364864) q[1];
sx q[1];
rz(2.6311325) q[1];
rz(-pi) q[2];
rz(2.8191928) q[3];
sx q[3];
rz(-1.673867) q[3];
sx q[3];
rz(-2.2565528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0024857) q[2];
sx q[2];
rz(-1.3287013) q[2];
sx q[2];
rz(-2.4277182) q[2];
rz(-0.21126963) q[3];
sx q[3];
rz(-2.1043089) q[3];
sx q[3];
rz(2.5762288) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7738889) q[0];
sx q[0];
rz(-1.1834894) q[0];
sx q[0];
rz(1.1327889) q[0];
rz(2.6702113) q[1];
sx q[1];
rz(-1.2939204) q[1];
sx q[1];
rz(0.43101355) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.80167) q[0];
sx q[0];
rz(-2.9012052) q[0];
sx q[0];
rz(0.49317081) q[0];
x q[1];
rz(-1.9507031) q[2];
sx q[2];
rz(-1.9089193) q[2];
sx q[2];
rz(-2.9924711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9543332) q[1];
sx q[1];
rz(-1.0403087) q[1];
sx q[1];
rz(1.2198636) q[1];
rz(-pi) q[2];
rz(0.27121284) q[3];
sx q[3];
rz(-2.1749139) q[3];
sx q[3];
rz(-0.97424265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5203349) q[2];
sx q[2];
rz(-0.73953491) q[2];
sx q[2];
rz(-0.4062824) q[2];
rz(1.2064365) q[3];
sx q[3];
rz(-2.0483569) q[3];
sx q[3];
rz(2.5549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63221145) q[0];
sx q[0];
rz(-2.2659232) q[0];
sx q[0];
rz(-0.32633728) q[0];
rz(-0.95411602) q[1];
sx q[1];
rz(-2.4931144) q[1];
sx q[1];
rz(-0.023524806) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2187754) q[0];
sx q[0];
rz(-1.9440224) q[0];
sx q[0];
rz(1.1362227) q[0];
rz(1.8872702) q[2];
sx q[2];
rz(-0.9372406) q[2];
sx q[2];
rz(-0.67297602) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4386846) q[1];
sx q[1];
rz(-0.45346916) q[1];
sx q[1];
rz(0.87345846) q[1];
rz(-pi) q[2];
rz(0.29554772) q[3];
sx q[3];
rz(-1.4369399) q[3];
sx q[3];
rz(-0.16242151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.894459) q[2];
sx q[2];
rz(-2.6829312) q[2];
sx q[2];
rz(0.66391724) q[2];
rz(-0.89020056) q[3];
sx q[3];
rz(-1.592344) q[3];
sx q[3];
rz(-3.1288872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.3506055) q[0];
sx q[0];
rz(-2.7101639) q[0];
sx q[0];
rz(-0.33082333) q[0];
rz(-0.36258969) q[1];
sx q[1];
rz(-1.1793084) q[1];
sx q[1];
rz(0.51574743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4048201) q[0];
sx q[0];
rz(-1.3546001) q[0];
sx q[0];
rz(1.3275073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.547566) q[2];
sx q[2];
rz(-1.4129248) q[2];
sx q[2];
rz(-1.3321587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1860604) q[1];
sx q[1];
rz(-1.493425) q[1];
sx q[1];
rz(-2.2633228) q[1];
rz(-pi) q[2];
rz(0.79591085) q[3];
sx q[3];
rz(-1.6636968) q[3];
sx q[3];
rz(-1.9167629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35046878) q[2];
sx q[2];
rz(-1.4501269) q[2];
sx q[2];
rz(2.3340732) q[2];
rz(-3.0430072) q[3];
sx q[3];
rz(-2.2435296) q[3];
sx q[3];
rz(-1.0928104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7035141) q[0];
sx q[0];
rz(-2.5683537) q[0];
sx q[0];
rz(-0.44952965) q[0];
rz(0.685177) q[1];
sx q[1];
rz(-0.40760577) q[1];
sx q[1];
rz(0.61946851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9307053) q[0];
sx q[0];
rz(-0.31085098) q[0];
sx q[0];
rz(0.092664552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6190104) q[2];
sx q[2];
rz(-1.0430665) q[2];
sx q[2];
rz(0.22942782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0692783) q[1];
sx q[1];
rz(-0.89767805) q[1];
sx q[1];
rz(2.0050383) q[1];
rz(1.3656627) q[3];
sx q[3];
rz(-1.3899125) q[3];
sx q[3];
rz(-3.0360529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56457907) q[2];
sx q[2];
rz(-1.4489633) q[2];
sx q[2];
rz(0.1980093) q[2];
rz(1.3624582) q[3];
sx q[3];
rz(-2.8415316) q[3];
sx q[3];
rz(2.4411966) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.144416) q[0];
sx q[0];
rz(-0.62189019) q[0];
sx q[0];
rz(1.8560334) q[0];
rz(-2.8649435) q[1];
sx q[1];
rz(-1.7729019) q[1];
sx q[1];
rz(-3.0329774) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40688694) q[0];
sx q[0];
rz(-2.6331775) q[0];
sx q[0];
rz(1.4020573) q[0];
rz(0.94679657) q[2];
sx q[2];
rz(-1.5063926) q[2];
sx q[2];
rz(-1.4381806) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.120879) q[1];
sx q[1];
rz(-0.87057897) q[1];
sx q[1];
rz(0.28182272) q[1];
rz(-1.2694938) q[3];
sx q[3];
rz(-2.7387894) q[3];
sx q[3];
rz(-0.81217743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8331929) q[2];
sx q[2];
rz(-1.1966285) q[2];
sx q[2];
rz(1.8673372) q[2];
rz(1.3537004) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(2.275009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891147) q[0];
sx q[0];
rz(-0.72565961) q[0];
sx q[0];
rz(3.0210378) q[0];
rz(2.4001135) q[1];
sx q[1];
rz(-1.2040141) q[1];
sx q[1];
rz(3.0659058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8814319) q[0];
sx q[0];
rz(-2.7921225) q[0];
sx q[0];
rz(0.49728877) q[0];
rz(-3.0641714) q[2];
sx q[2];
rz(-2.0311557) q[2];
sx q[2];
rz(-1.8094687) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2774128) q[1];
sx q[1];
rz(-2.2267739) q[1];
sx q[1];
rz(-1.6026864) q[1];
rz(-pi) q[2];
rz(-0.60093694) q[3];
sx q[3];
rz(-1.4190201) q[3];
sx q[3];
rz(2.5545718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.072448298) q[2];
sx q[2];
rz(-2.2587903) q[2];
sx q[2];
rz(2.7952588) q[2];
rz(-2.196178) q[3];
sx q[3];
rz(-1.8190705) q[3];
sx q[3];
rz(0.94714981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866975) q[0];
sx q[0];
rz(-1.1981542) q[0];
sx q[0];
rz(-0.44788885) q[0];
rz(2.2121494) q[1];
sx q[1];
rz(-1.399469) q[1];
sx q[1];
rz(-0.67451745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436627) q[0];
sx q[0];
rz(-0.81564689) q[0];
sx q[0];
rz(-1.2944512) q[0];
rz(-pi) q[1];
rz(-0.10402502) q[2];
sx q[2];
rz(-1.7556223) q[2];
sx q[2];
rz(1.8341174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0793987) q[1];
sx q[1];
rz(-0.17339686) q[1];
sx q[1];
rz(-1.5773415) q[1];
x q[2];
rz(1.4025979) q[3];
sx q[3];
rz(-0.4240331) q[3];
sx q[3];
rz(3.0454717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21005361) q[2];
sx q[2];
rz(-1.9320107) q[2];
sx q[2];
rz(-2.8741969) q[2];
rz(2.0343272) q[3];
sx q[3];
rz(-0.78247726) q[3];
sx q[3];
rz(-2.2926008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.7312819) q[0];
sx q[0];
rz(-1.5189497) q[0];
sx q[0];
rz(2.0547163) q[0];
rz(0.80401737) q[1];
sx q[1];
rz(-1.6537279) q[1];
sx q[1];
rz(-2.0860685) q[1];
rz(-2.0217996) q[2];
sx q[2];
rz(-2.9783772) q[2];
sx q[2];
rz(2.1676248) q[2];
rz(0.90292061) q[3];
sx q[3];
rz(-1.5048448) q[3];
sx q[3];
rz(2.0497143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
