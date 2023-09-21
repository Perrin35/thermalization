OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1113623) q[0];
sx q[0];
rz(-2.4863939) q[0];
sx q[0];
rz(2.1652048) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108609) q[0];
sx q[0];
rz(-2.2415677) q[0];
sx q[0];
rz(3.1109111) q[0];
x q[1];
rz(-2.5976074) q[2];
sx q[2];
rz(-1.7684801) q[2];
sx q[2];
rz(2.939784) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77289171) q[1];
sx q[1];
rz(-1.643631) q[1];
sx q[1];
rz(-2.021695) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3402432) q[3];
sx q[3];
rz(-1.0145463) q[3];
sx q[3];
rz(-2.706561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4684326) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(3.1288778) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(2.7761249) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8651543) q[0];
sx q[0];
rz(-1.8218453) q[0];
sx q[0];
rz(-2.6692078) q[0];
rz(-pi) q[1];
rz(0.89895504) q[2];
sx q[2];
rz(-2.0546753) q[2];
sx q[2];
rz(-1.3324141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47479113) q[1];
sx q[1];
rz(-0.7102237) q[1];
sx q[1];
rz(0.53479654) q[1];
rz(2.9404853) q[3];
sx q[3];
rz(-1.6893941) q[3];
sx q[3];
rz(0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.83313292) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(0.45062137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492309) q[0];
sx q[0];
rz(-1.5376687) q[0];
sx q[0];
rz(-1.4928198) q[0];
x q[1];
rz(-0.89183715) q[2];
sx q[2];
rz(-2.6747181) q[2];
sx q[2];
rz(-2.4613949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8038517) q[1];
sx q[1];
rz(-0.9461113) q[1];
sx q[1];
rz(1.9546263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2080473) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(-2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(-2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(1.2966688) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(2.8362714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353377) q[0];
sx q[0];
rz(-2.4843785) q[0];
sx q[0];
rz(-0.52657907) q[0];
x q[1];
rz(0.73037578) q[2];
sx q[2];
rz(-0.59130284) q[2];
sx q[2];
rz(-0.060746047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0381283) q[1];
sx q[1];
rz(-0.80059073) q[1];
sx q[1];
rz(2.994328) q[1];
rz(-pi) q[2];
rz(-2.3577945) q[3];
sx q[3];
rz(-2.4664306) q[3];
sx q[3];
rz(-2.6019707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2083464) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(-2.7944881) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-2.7358823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680506) q[0];
sx q[0];
rz(-1.2965856) q[0];
sx q[0];
rz(-0.6875086) q[0];
x q[1];
rz(1.97498) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(0.052978901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2629562) q[1];
sx q[1];
rz(-2.7107781) q[1];
sx q[1];
rz(2.2104435) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76463502) q[3];
sx q[3];
rz(-1.7406165) q[3];
sx q[3];
rz(3.0499383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.3173332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989527) q[0];
sx q[0];
rz(-1.5931804) q[0];
sx q[0];
rz(-1.0277261) q[0];
rz(-pi) q[1];
rz(-2.1520734) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(1.6473824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.085560872) q[1];
sx q[1];
rz(-2.3059658) q[1];
sx q[1];
rz(-0.0037770165) q[1];
rz(-pi) q[2];
rz(-0.96885724) q[3];
sx q[3];
rz(-0.67776206) q[3];
sx q[3];
rz(-2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(-1.9334531) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.5642501) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(1.1279001) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82176498) q[0];
sx q[0];
rz(-1.9263096) q[0];
sx q[0];
rz(2.6000644) q[0];
rz(-pi) q[1];
rz(0.77881323) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(-0.29701172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5269446) q[1];
sx q[1];
rz(-1.6741747) q[1];
sx q[1];
rz(-2.9109216) q[1];
rz(-pi) q[2];
rz(-0.96630649) q[3];
sx q[3];
rz(-1.5887727) q[3];
sx q[3];
rz(-3.0997961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82905519) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-2.5404239) q[2];
rz(-0.51856315) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475875) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(3.0684379) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97544599) q[0];
sx q[0];
rz(-1.5243634) q[0];
sx q[0];
rz(2.5509044) q[0];
rz(-1.0560889) q[2];
sx q[2];
rz(-2.0851496) q[2];
sx q[2];
rz(-0.27013847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29443024) q[1];
sx q[1];
rz(-1.8995598) q[1];
sx q[1];
rz(-0.86818236) q[1];
rz(2.153272) q[3];
sx q[3];
rz(-2.6675825) q[3];
sx q[3];
rz(-1.7085027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56117326) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(0.013307868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028037926) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-0.56217271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(1.1086585) q[0];
rz(-pi) q[1];
rz(1.4116889) q[2];
sx q[2];
rz(-2.6052887) q[2];
sx q[2];
rz(-1.2250587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0253898) q[1];
sx q[1];
rz(-2.2491978) q[1];
sx q[1];
rz(-2.5118026) q[1];
rz(-pi) q[2];
rz(-2.3638944) q[3];
sx q[3];
rz(-0.54293984) q[3];
sx q[3];
rz(-0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(-2.1600058) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(-2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-0.04709588) q[0];
rz(-1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809966) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(3.1050443) q[0];
rz(-pi) q[1];
rz(3.0907862) q[2];
sx q[2];
rz(-2.0127957) q[2];
sx q[2];
rz(-0.49651422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8216985) q[1];
sx q[1];
rz(-2.387429) q[1];
sx q[1];
rz(0.9174222) q[1];
rz(-0.40269561) q[3];
sx q[3];
rz(-1.6214569) q[3];
sx q[3];
rz(-2.5718873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(2.2052374) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3180278) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(1.5488831) q[2];
sx q[2];
rz(-2.1219325) q[2];
sx q[2];
rz(1.3757201) q[2];
rz(-1.9746642) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];