OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(0.93044257) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247834) q[0];
sx q[0];
rz(-2.2342355) q[0];
sx q[0];
rz(-1.1532564) q[0];
x q[1];
rz(-0.47714699) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(-1.9805816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0336696) q[1];
sx q[1];
rz(-2.8406997) q[1];
sx q[1];
rz(1.184102) q[1];
rz(-pi) q[2];
rz(-1.1125621) q[3];
sx q[3];
rz(-2.4535865) q[3];
sx q[3];
rz(2.1306761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0779695) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(1.8135653) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-3.0096171) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(0.81726384) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5215223) q[0];
sx q[0];
rz(-1.9154473) q[0];
sx q[0];
rz(-2.9917813) q[0];
x q[1];
rz(1.9643289) q[2];
sx q[2];
rz(-2.6173008) q[2];
sx q[2];
rz(-0.92698586) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0751794) q[1];
sx q[1];
rz(-1.315409) q[1];
sx q[1];
rz(-2.4875719) q[1];
rz(1.9278139) q[3];
sx q[3];
rz(-0.5001874) q[3];
sx q[3];
rz(-0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2237504) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(0.98637995) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36689511) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(-1.6945217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7265978) q[0];
sx q[0];
rz(-1.5883755) q[0];
sx q[0];
rz(2.9265762) q[0];
x q[1];
rz(-1.8410728) q[2];
sx q[2];
rz(-2.8405361) q[2];
sx q[2];
rz(-1.9384055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9008873) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(-2.4918873) q[1];
rz(-pi) q[2];
rz(-1.9871739) q[3];
sx q[3];
rz(-1.3458816) q[3];
sx q[3];
rz(0.51605663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(-2.036371) q[2];
rz(-0.74622074) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(-2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063909) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(0.28494596) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(-2.7526061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8303191) q[0];
sx q[0];
rz(-1.3415601) q[0];
sx q[0];
rz(2.8515408) q[0];
x q[1];
rz(2.94611) q[2];
sx q[2];
rz(-2.0494378) q[2];
sx q[2];
rz(-0.10749707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3868689) q[1];
sx q[1];
rz(-1.9935973) q[1];
sx q[1];
rz(-0.60476426) q[1];
rz(2.3936478) q[3];
sx q[3];
rz(-1.4191295) q[3];
sx q[3];
rz(-2.9880854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(-2.8592498) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221955) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(-2.7815681) q[0];
rz(-0.64741627) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-0.49450758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7617103) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(2.7193927) q[0];
x q[1];
rz(-1.8525271) q[2];
sx q[2];
rz(-0.14094555) q[2];
sx q[2];
rz(2.3245387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8620017) q[1];
sx q[1];
rz(-1.630736) q[1];
sx q[1];
rz(2.5111591) q[1];
rz(-0.4440998) q[3];
sx q[3];
rz(-1.5478304) q[3];
sx q[3];
rz(0.75957739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(0.29850706) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(2.503094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9962149) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(0.021082489) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(1.235199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(1.8963277) q[0];
rz(-pi) q[1];
rz(-2.4371106) q[2];
sx q[2];
rz(-1.8330169) q[2];
sx q[2];
rz(0.30202497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1393309) q[1];
sx q[1];
rz(-1.9863223) q[1];
sx q[1];
rz(0.63654391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33552334) q[3];
sx q[3];
rz(-2.8918859) q[3];
sx q[3];
rz(1.6186796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49332508) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(2.7098999) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(-3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(0.11211638) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(2.0281866) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413234) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(-1.8510438) q[0];
rz(-pi) q[1];
rz(2.4593997) q[2];
sx q[2];
rz(-2.7436514) q[2];
sx q[2];
rz(-1.8856018) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50892936) q[1];
sx q[1];
rz(-0.76970657) q[1];
sx q[1];
rz(0.44429227) q[1];
rz(-pi) q[2];
rz(1.1464305) q[3];
sx q[3];
rz(-0.48454912) q[3];
sx q[3];
rz(0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(-0.12602885) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(-2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(-1.460176) q[0];
rz(1.8966282) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(-1.9326899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9019933) q[0];
sx q[0];
rz(-1.8327466) q[0];
sx q[0];
rz(0.082017935) q[0];
rz(-pi) q[1];
rz(-2.3102343) q[2];
sx q[2];
rz(-2.0908329) q[2];
sx q[2];
rz(1.3607197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9637451) q[1];
sx q[1];
rz(-0.81804619) q[1];
sx q[1];
rz(0.45958105) q[1];
rz(0.21945159) q[3];
sx q[3];
rz(-0.95153522) q[3];
sx q[3];
rz(2.0451562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(1.8219927) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.9086583) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(-2.7499054) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101333) q[0];
sx q[0];
rz(-2.0053232) q[0];
sx q[0];
rz(-1.5772485) q[0];
rz(-0.45231818) q[2];
sx q[2];
rz(-1.3903793) q[2];
sx q[2];
rz(-2.7197321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0605676) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(-0.4391567) q[1];
rz(-1.6438946) q[3];
sx q[3];
rz(-2.0953296) q[3];
sx q[3];
rz(0.91367164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.8048145) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-2.5706932) q[0];
rz(1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-0.16194078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8853332) q[0];
sx q[0];
rz(-2.4618039) q[0];
sx q[0];
rz(2.0392366) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5230721) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(1.7392841) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.098617741) q[1];
sx q[1];
rz(-2.9264755) q[1];
sx q[1];
rz(-2.2153562) q[1];
rz(-pi) q[2];
rz(3.083769) q[3];
sx q[3];
rz(-2.6177546) q[3];
sx q[3];
rz(-3.0081089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1841715) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.4282248) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.4409055) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(-1.3700925) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(-0.8926819) q[2];
sx q[2];
rz(-0.18845367) q[2];
sx q[2];
rz(-1.0343196) q[2];
rz(0.45392848) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];