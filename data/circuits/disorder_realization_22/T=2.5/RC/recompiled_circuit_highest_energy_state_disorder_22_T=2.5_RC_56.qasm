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
rz(-2.9838188) q[0];
sx q[0];
rz(4.3133419) q[0];
sx q[0];
rz(11.316909) q[0];
rz(-0.14021048) q[1];
sx q[1];
rz(-1.6970716) q[1];
sx q[1];
rz(0.061847774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2551198) q[0];
sx q[0];
rz(-0.97772163) q[0];
sx q[0];
rz(-2.5307239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5704448) q[2];
sx q[2];
rz(-1.4803737) q[2];
sx q[2];
rz(0.11133678) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5420336) q[1];
sx q[1];
rz(-2.2890511) q[1];
sx q[1];
rz(-0.12964779) q[1];
rz(1.0855882) q[3];
sx q[3];
rz(-2.3923529) q[3];
sx q[3];
rz(3.0390714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42741117) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(2.8851435) q[2];
rz(-0.17476684) q[3];
sx q[3];
rz(-1.9758965) q[3];
sx q[3];
rz(-0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8137708) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-3.0178965) q[0];
rz(0.23873121) q[1];
sx q[1];
rz(-0.78014603) q[1];
sx q[1];
rz(0.10496584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86228131) q[0];
sx q[0];
rz(-2.263592) q[0];
sx q[0];
rz(0.19578085) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6506422) q[2];
sx q[2];
rz(-2.6770795) q[2];
sx q[2];
rz(2.2238942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2626896) q[1];
sx q[1];
rz(-0.57003747) q[1];
sx q[1];
rz(2.0712584) q[1];
rz(-1.9267004) q[3];
sx q[3];
rz(-2.2442563) q[3];
sx q[3];
rz(-1.9178903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3671942) q[2];
sx q[2];
rz(-1.4906733) q[2];
sx q[2];
rz(1.6801838) q[2];
rz(-1.8558308) q[3];
sx q[3];
rz(-1.4444084) q[3];
sx q[3];
rz(2.8694966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9937781) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(1.0149957) q[0];
rz(-0.89730942) q[1];
sx q[1];
rz(-1.4941447) q[1];
sx q[1];
rz(0.6764594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1508866) q[0];
sx q[0];
rz(-0.61752049) q[0];
sx q[0];
rz(-1.6157009) q[0];
x q[1];
rz(0.93134201) q[2];
sx q[2];
rz(-0.73340511) q[2];
sx q[2];
rz(2.1426107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3520917) q[1];
sx q[1];
rz(-1.3695696) q[1];
sx q[1];
rz(1.1835062) q[1];
rz(-pi) q[2];
rz(3.1234497) q[3];
sx q[3];
rz(-0.81644316) q[3];
sx q[3];
rz(1.3870365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2066388) q[2];
sx q[2];
rz(-0.97524869) q[2];
sx q[2];
rz(2.3812531) q[2];
rz(1.9350516) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(-2.2981203) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7059785) q[0];
sx q[0];
rz(-2.5137081) q[0];
sx q[0];
rz(1.5465558) q[0];
rz(-0.71890038) q[1];
sx q[1];
rz(-1.8214106) q[1];
sx q[1];
rz(2.9100606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084615413) q[0];
sx q[0];
rz(-1.1418793) q[0];
sx q[0];
rz(-2.5713628) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1377085) q[2];
sx q[2];
rz(-0.79266119) q[2];
sx q[2];
rz(-2.5156227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4930112) q[1];
sx q[1];
rz(-1.4775044) q[1];
sx q[1];
rz(-2.6562163) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1778963) q[3];
sx q[3];
rz(-1.0024973) q[3];
sx q[3];
rz(1.5400122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.71451688) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(2.5330949) q[2];
rz(0.19615873) q[3];
sx q[3];
rz(-1.720263) q[3];
sx q[3];
rz(0.99854809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8683559) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(-2.3625145) q[0];
rz(-0.92998663) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(-1.1735865) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0239651) q[0];
sx q[0];
rz(-1.5668586) q[0];
sx q[0];
rz(1.2327475) q[0];
x q[1];
rz(1.3814244) q[2];
sx q[2];
rz(-0.62386419) q[2];
sx q[2];
rz(2.3331235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8820937) q[1];
sx q[1];
rz(-2.8176687) q[1];
sx q[1];
rz(-1.7739576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7608123) q[3];
sx q[3];
rz(-1.1486067) q[3];
sx q[3];
rz(1.9334396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.32213) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(-2.4665311) q[2];
rz(0.86554646) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(-1.9338098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753767) q[0];
sx q[0];
rz(-3.1238811) q[0];
sx q[0];
rz(-2.2702763) q[0];
rz(-0.21698347) q[1];
sx q[1];
rz(-1.5419518) q[1];
sx q[1];
rz(1.8035696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965391) q[0];
sx q[0];
rz(-0.65700475) q[0];
sx q[0];
rz(-0.13583247) q[0];
rz(-pi) q[1];
rz(-0.14071847) q[2];
sx q[2];
rz(-0.85007668) q[2];
sx q[2];
rz(0.36397935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1068145) q[1];
sx q[1];
rz(-2.9902774) q[1];
sx q[1];
rz(-1.2447912) q[1];
rz(-pi) q[2];
rz(-1.7609672) q[3];
sx q[3];
rz(-2.0479879) q[3];
sx q[3];
rz(-2.793345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.013082144) q[2];
sx q[2];
rz(-1.6435813) q[2];
sx q[2];
rz(2.3415372) q[2];
rz(0.99177805) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(2.1984524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8868788) q[0];
sx q[0];
rz(-0.44097057) q[0];
sx q[0];
rz(0.15175858) q[0];
rz(-1.7249379) q[1];
sx q[1];
rz(-2.2817426) q[1];
sx q[1];
rz(2.3053665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0043512303) q[0];
sx q[0];
rz(-2.3286331) q[0];
sx q[0];
rz(1.0197958) q[0];
rz(3.128746) q[2];
sx q[2];
rz(-2.564866) q[2];
sx q[2];
rz(1.6614514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9735896) q[1];
sx q[1];
rz(-0.26339809) q[1];
sx q[1];
rz(-1.2348639) q[1];
rz(-pi) q[2];
rz(0.25612513) q[3];
sx q[3];
rz(-0.47893347) q[3];
sx q[3];
rz(-2.0253945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4862711) q[2];
sx q[2];
rz(-0.063491193) q[2];
sx q[2];
rz(-3.0625694) q[2];
rz(-2.4231353) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(0.92459905) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543168) q[0];
sx q[0];
rz(-0.59816718) q[0];
sx q[0];
rz(0.86791903) q[0];
rz(0.68880853) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(-0.93592962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.717632) q[0];
sx q[0];
rz(-2.0905295) q[0];
sx q[0];
rz(-1.3613767) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9804968) q[2];
sx q[2];
rz(-1.2563224) q[2];
sx q[2];
rz(1.0057698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1327428) q[1];
sx q[1];
rz(-2.8121083) q[1];
sx q[1];
rz(-1.4047755) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5722365) q[3];
sx q[3];
rz(-0.82732302) q[3];
sx q[3];
rz(1.2900521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1615289) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(2.1234546) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(0.97091466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-2.7629906) q[0];
sx q[0];
rz(-0.55736962) q[0];
sx q[0];
rz(-0.85365224) q[0];
rz(1.254982) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(-2.8010211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69215323) q[0];
sx q[0];
rz(-1.3763577) q[0];
sx q[0];
rz(-2.9783335) q[0];
x q[1];
rz(1.8179632) q[2];
sx q[2];
rz(-2.2296612) q[2];
sx q[2];
rz(-0.62550046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3601968) q[1];
sx q[1];
rz(-0.92222795) q[1];
sx q[1];
rz(2.1482723) q[1];
rz(-pi) q[2];
rz(2.0331618) q[3];
sx q[3];
rz(-1.3232854) q[3];
sx q[3];
rz(-0.68456291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7178932) q[2];
sx q[2];
rz(-1.7607949) q[2];
sx q[2];
rz(-2.763486) q[2];
rz(0.2374436) q[3];
sx q[3];
rz(-1.0708555) q[3];
sx q[3];
rz(0.28089359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0963999) q[0];
sx q[0];
rz(-2.0572331) q[0];
sx q[0];
rz(-0.64724809) q[0];
rz(1.9735533) q[1];
sx q[1];
rz(-0.85473514) q[1];
sx q[1];
rz(-2.7580269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76825324) q[0];
sx q[0];
rz(-3.1041234) q[0];
sx q[0];
rz(-2.5277407) q[0];
x q[1];
rz(0.52807669) q[2];
sx q[2];
rz(-2.5438692) q[2];
sx q[2];
rz(0.18250386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4870461) q[1];
sx q[1];
rz(-0.33815835) q[1];
sx q[1];
rz(0.66429331) q[1];
rz(-pi) q[2];
rz(2.0862574) q[3];
sx q[3];
rz(-1.664242) q[3];
sx q[3];
rz(2.6938113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26455227) q[2];
sx q[2];
rz(-0.88901797) q[2];
sx q[2];
rz(1.5569347) q[2];
rz(1.6282188) q[3];
sx q[3];
rz(-1.3686562) q[3];
sx q[3];
rz(-2.7148066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.3601111) q[0];
sx q[0];
rz(-1.4327015) q[0];
sx q[0];
rz(-2.844368) q[0];
rz(1.1813286) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(-2.0484106) q[2];
sx q[2];
rz(-1.3543159) q[2];
sx q[2];
rz(-1.3456624) q[2];
rz(-0.48458002) q[3];
sx q[3];
rz(-1.598834) q[3];
sx q[3];
rz(-1.9091189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
