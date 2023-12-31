OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(5.7603523) q[0];
sx q[0];
rz(8.8011959) q[0];
rz(3.4317598) q[1];
sx q[1];
rz(5.5640339) q[1];
sx q[1];
rz(13.066864) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1815163) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(2.5417679) q[0];
rz(-pi) q[1];
rz(2.1262769) q[2];
sx q[2];
rz(-0.17586389) q[2];
sx q[2];
rz(-1.442684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11195586) q[1];
sx q[1];
rz(-1.2458548) q[1];
sx q[1];
rz(0.8128266) q[1];
rz(-pi) q[2];
rz(-0.15668232) q[3];
sx q[3];
rz(-1.1374258) q[3];
sx q[3];
rz(-2.3436848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(1.2228489) q[2];
rz(1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37110776) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(-1.0789385) q[0];
rz(-1.3868015) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9294445) q[0];
sx q[0];
rz(-1.11709) q[0];
sx q[0];
rz(1.8871904) q[0];
x q[1];
rz(-1.7416679) q[2];
sx q[2];
rz(-1.2534281) q[2];
sx q[2];
rz(-2.3651809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14109719) q[1];
sx q[1];
rz(-0.54447237) q[1];
sx q[1];
rz(-2.673124) q[1];
x q[2];
rz(2.2802248) q[3];
sx q[3];
rz(-1.5449761) q[3];
sx q[3];
rz(-0.65755075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60454303) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(-1.4705307) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(-0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(-1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62832075) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(0.4011641) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7064352) q[2];
sx q[2];
rz(-1.7559768) q[2];
sx q[2];
rz(-2.8806825) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2565508) q[1];
sx q[1];
rz(-2.1532144) q[1];
sx q[1];
rz(-0.38862733) q[1];
rz(2.9754144) q[3];
sx q[3];
rz(-1.8718534) q[3];
sx q[3];
rz(-2.500246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.65163461) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(-2.4839694) q[2];
rz(-1.1714606) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(-1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.3477429) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.4452274) q[0];
rz(-1.6943278) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(0.34805527) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9075273) q[0];
sx q[0];
rz(-2.4287927) q[0];
sx q[0];
rz(-2.6839921) q[0];
rz(-0.15050998) q[2];
sx q[2];
rz(-1.9410656) q[2];
sx q[2];
rz(2.9589257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0192249) q[1];
sx q[1];
rz(-0.61673635) q[1];
sx q[1];
rz(1.7255746) q[1];
rz(-pi) q[2];
rz(-0.16163687) q[3];
sx q[3];
rz(-1.4378387) q[3];
sx q[3];
rz(-1.1383575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(-0.73741284) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(-0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(0.53043956) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.2984498) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2753678) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(-2.2989681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60646306) q[2];
sx q[2];
rz(-1.3242553) q[2];
sx q[2];
rz(-2.390887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7057719) q[1];
sx q[1];
rz(-2.768369) q[1];
sx q[1];
rz(-2.865764) q[1];
rz(-pi) q[2];
rz(-0.50000639) q[3];
sx q[3];
rz(-1.0370504) q[3];
sx q[3];
rz(1.223525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-0.47362622) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(-0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-2.6748437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022910718) q[0];
sx q[0];
rz(-2.4372299) q[0];
sx q[0];
rz(0.33606152) q[0];
rz(-pi) q[1];
rz(-0.42748638) q[2];
sx q[2];
rz(-1.7726521) q[2];
sx q[2];
rz(2.0919378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0007243) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
rz(3.0454037) q[3];
sx q[3];
rz(-1.4078762) q[3];
sx q[3];
rz(-0.067226203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5647282) q[2];
rz(-2.5148897) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(0.68157649) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8108114) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(2.7217641) q[0];
rz(-0.22142521) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2433462) q[0];
sx q[0];
rz(-1.4833437) q[0];
sx q[0];
rz(-1.5792219) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8553472) q[2];
sx q[2];
rz(-1.7562521) q[2];
sx q[2];
rz(2.4539349) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6429813) q[1];
sx q[1];
rz(-1.4652518) q[1];
sx q[1];
rz(1.5555744) q[1];
x q[2];
rz(-1.1057304) q[3];
sx q[3];
rz(-1.2601868) q[3];
sx q[3];
rz(3.0748933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(0.62492433) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4411366) q[0];
sx q[0];
rz(-1.6807798) q[0];
sx q[0];
rz(2.3814047) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74771379) q[2];
sx q[2];
rz(-1.3249825) q[2];
sx q[2];
rz(-0.73033787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18409477) q[1];
sx q[1];
rz(-3.0093319) q[1];
sx q[1];
rz(1.9930507) q[1];
rz(-0.30808361) q[3];
sx q[3];
rz(-0.83762533) q[3];
sx q[3];
rz(-1.777491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4079995) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(-2.5949809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1221736) q[0];
sx q[0];
rz(-1.1328567) q[0];
sx q[0];
rz(0.41284783) q[0];
x q[1];
rz(-2.390929) q[2];
sx q[2];
rz(-0.51807907) q[2];
sx q[2];
rz(2.2028365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4709028) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(-1.9914658) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2253753) q[3];
sx q[3];
rz(-1.1715874) q[3];
sx q[3];
rz(0.44772128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(1.6938422) q[2];
rz(2.0041806) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820213) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(-2.4243673) q[0];
rz(1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9428064) q[0];
sx q[0];
rz(-1.2278623) q[0];
sx q[0];
rz(-2.2862611) q[0];
rz(-pi) q[1];
rz(0.14214469) q[2];
sx q[2];
rz(-1.7927577) q[2];
sx q[2];
rz(0.97505002) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.087079436) q[1];
sx q[1];
rz(-2.2594249) q[1];
sx q[1];
rz(-0.41333945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35499771) q[3];
sx q[3];
rz(-0.98155752) q[3];
sx q[3];
rz(-2.0139351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83508867) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(1.1076526) q[2];
sx q[2];
rz(-1.9911498) q[2];
sx q[2];
rz(-0.1291612) q[2];
rz(-0.64745263) q[3];
sx q[3];
rz(-0.068131937) q[3];
sx q[3];
rz(1.6905231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
