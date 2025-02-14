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
rz(-2.2263865) q[0];
sx q[0];
rz(-2.6829166) q[0];
sx q[0];
rz(-0.31804481) q[0];
rz(-1.7755427) q[1];
sx q[1];
rz(-0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013115766) q[0];
sx q[0];
rz(-2.4139691) q[0];
sx q[0];
rz(-2.3250871) q[0];
rz(-pi) q[1];
rz(-0.89102052) q[2];
sx q[2];
rz(-1.2212379) q[2];
sx q[2];
rz(2.120979) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1789238) q[1];
sx q[1];
rz(-2.0589863) q[1];
sx q[1];
rz(-0.31655426) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99394057) q[3];
sx q[3];
rz(-1.9074884) q[3];
sx q[3];
rz(1.5722881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.723145) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-2.1212228) q[2];
rz(-2.4127035) q[3];
sx q[3];
rz(-2.708669) q[3];
sx q[3];
rz(-1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51434022) q[0];
sx q[0];
rz(-1.1662551) q[0];
sx q[0];
rz(3.1033206) q[0];
rz(-0.12380869) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(-1.570787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1403303) q[0];
sx q[0];
rz(-1.5902983) q[0];
sx q[0];
rz(0.01000761) q[0];
rz(-pi) q[1];
rz(-0.77124243) q[2];
sx q[2];
rz(-1.6420157) q[2];
sx q[2];
rz(-2.3598755) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1473737) q[1];
sx q[1];
rz(-0.0011073907) q[1];
sx q[1];
rz(-1.1785517) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6255076) q[3];
sx q[3];
rz(-1.4654963) q[3];
sx q[3];
rz(2.1886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52246419) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(-2.5627356) q[2];
rz(2.0959496) q[3];
sx q[3];
rz(-0.78543109) q[3];
sx q[3];
rz(2.3845909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-0.45397595) q[0];
rz(-2.9767735) q[1];
sx q[1];
rz(-0.77607981) q[1];
sx q[1];
rz(-1.4705315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1810581) q[0];
sx q[0];
rz(-1.9583324) q[0];
sx q[0];
rz(-0.18523943) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8118089) q[2];
sx q[2];
rz(-1.068914) q[2];
sx q[2];
rz(1.2041853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78782234) q[1];
sx q[1];
rz(-0.90836891) q[1];
sx q[1];
rz(-1.1819928) q[1];
x q[2];
rz(-0.64162125) q[3];
sx q[3];
rz(-2.372916) q[3];
sx q[3];
rz(0.25594433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37956023) q[2];
sx q[2];
rz(-1.4780059) q[2];
sx q[2];
rz(1.7851768) q[2];
rz(-0.49946579) q[3];
sx q[3];
rz(-1.6127337) q[3];
sx q[3];
rz(-2.0904026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9143963) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(1.1610441) q[1];
sx q[1];
rz(-0.5608905) q[1];
sx q[1];
rz(-0.11722359) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49699052) q[0];
sx q[0];
rz(-1.0287026) q[0];
sx q[0];
rz(-1.5250456) q[0];
x q[1];
rz(0.92553301) q[2];
sx q[2];
rz(-1.7286073) q[2];
sx q[2];
rz(0.35343808) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1438751) q[1];
sx q[1];
rz(-1.2923601) q[1];
sx q[1];
rz(-0.77634676) q[1];
rz(-2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(-1.6146631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3566572) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(-1.3680722) q[2];
rz(1.8562227) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(-0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042498978) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(-1.7294783) q[0];
rz(2.1389351) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(-1.5589421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5413645) q[0];
sx q[0];
rz(-1.1523655) q[0];
sx q[0];
rz(-0.61516841) q[0];
rz(0.095076128) q[2];
sx q[2];
rz(-0.72059599) q[2];
sx q[2];
rz(2.528185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4950558) q[1];
sx q[1];
rz(-0.55083129) q[1];
sx q[1];
rz(-1.6158094) q[1];
rz(-pi) q[2];
x q[2];
rz(3.047278) q[3];
sx q[3];
rz(-1.2504761) q[3];
sx q[3];
rz(1.3195147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5530508) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(-0.39349619) q[2];
rz(-0.95997512) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(-0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680962) q[0];
sx q[0];
rz(-3.0143026) q[0];
sx q[0];
rz(-0.18336503) q[0];
rz(2.6496437) q[1];
sx q[1];
rz(-0.79194561) q[1];
sx q[1];
rz(1.9221745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757965) q[0];
sx q[0];
rz(-1.4816493) q[0];
sx q[0];
rz(-1.6515948) q[0];
rz(-1.6492789) q[2];
sx q[2];
rz(-1.540479) q[2];
sx q[2];
rz(-0.91532367) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9522493) q[1];
sx q[1];
rz(-0.5533411) q[1];
sx q[1];
rz(-2.4207741) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7352571) q[3];
sx q[3];
rz(-2.6643848) q[3];
sx q[3];
rz(0.26613126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(-0.92665893) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(1.5026708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53967404) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(-0.71640054) q[0];
rz(0.40388233) q[1];
sx q[1];
rz(-1.5682181) q[1];
sx q[1];
rz(1.7049888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2879031) q[0];
sx q[0];
rz(-0.40659621) q[0];
sx q[0];
rz(-1.9415744) q[0];
rz(-pi) q[1];
rz(2.613343) q[2];
sx q[2];
rz(-1.6852323) q[2];
sx q[2];
rz(2.6334762) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0246504) q[1];
sx q[1];
rz(-2.7305909) q[1];
sx q[1];
rz(-2.890439) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3327504) q[3];
sx q[3];
rz(-1.0029965) q[3];
sx q[3];
rz(2.0140935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3639823) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(3.0762365) q[2];
rz(2.2743547) q[3];
sx q[3];
rz(-1.5012274) q[3];
sx q[3];
rz(1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.48677483) q[0];
sx q[0];
rz(-1.4210533) q[0];
sx q[0];
rz(2.4105893) q[0];
rz(-0.77516088) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(0.41466546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99299198) q[0];
sx q[0];
rz(-0.96382695) q[0];
sx q[0];
rz(-2.4803922) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3581267) q[2];
sx q[2];
rz(-1.050282) q[2];
sx q[2];
rz(-1.5368324) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9379256) q[1];
sx q[1];
rz(-1.3199184) q[1];
sx q[1];
rz(2.832042) q[1];
rz(-pi) q[2];
rz(3.0075668) q[3];
sx q[3];
rz(-1.7334337) q[3];
sx q[3];
rz(-2.7664879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9097462) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(-0.092378423) q[2];
rz(0.77629027) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6542776) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(3.1238632) q[0];
rz(-1.0461294) q[1];
sx q[1];
rz(-1.7283864) q[1];
sx q[1];
rz(-1.0606934) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720406) q[0];
sx q[0];
rz(-1.5020796) q[0];
sx q[0];
rz(-1.5848573) q[0];
rz(2.8834516) q[2];
sx q[2];
rz(-0.99810583) q[2];
sx q[2];
rz(-1.0024459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13421397) q[1];
sx q[1];
rz(-1.7430647) q[1];
sx q[1];
rz(0.80165792) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5117731) q[3];
sx q[3];
rz(-1.0477598) q[3];
sx q[3];
rz(-0.28429261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1549418) q[2];
sx q[2];
rz(-2.3255746) q[2];
sx q[2];
rz(-2.4972534) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(-2.2346066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3802721) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(-0.45183387) q[0];
rz(-1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(1.570328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1399978) q[0];
sx q[0];
rz(-0.17029914) q[0];
sx q[0];
rz(-3.1055121) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11650916) q[2];
sx q[2];
rz(-2.0429869) q[2];
sx q[2];
rz(1.0509059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45598509) q[1];
sx q[1];
rz(-2.7311014) q[1];
sx q[1];
rz(-1.2430771) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5785923) q[3];
sx q[3];
rz(-0.8484133) q[3];
sx q[3];
rz(2.9786199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62341225) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(-3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-0.51641974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8138206) q[0];
sx q[0];
rz(-1.0186503) q[0];
sx q[0];
rz(-3.0042197) q[0];
rz(0.2188006) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(-1.1458746) q[2];
sx q[2];
rz(-2.6617202) q[2];
sx q[2];
rz(-0.94750994) q[2];
rz(2.6811213) q[3];
sx q[3];
rz(-0.71131592) q[3];
sx q[3];
rz(0.55874353) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
