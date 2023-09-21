OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(-0.65519873) q[0];
sx q[0];
rz(0.97638786) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3800669) q[0];
sx q[0];
rz(-2.4702284) q[0];
sx q[0];
rz(1.5321561) q[0];
x q[1];
rz(0.54398529) q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-2.3789062) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(3.0607037) q[1];
rz(2.3402432) q[3];
sx q[3];
rz(-1.0145463) q[3];
sx q[3];
rz(0.43503161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(2.7761249) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7468837) q[0];
sx q[0];
rz(-0.53046662) q[0];
sx q[0];
rz(0.51325004) q[0];
x q[1];
rz(-0.89895504) q[2];
sx q[2];
rz(-2.0546753) q[2];
sx q[2];
rz(1.3324141) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9531627) q[1];
sx q[1];
rz(-0.97524446) q[1];
sx q[1];
rz(-1.9838536) q[1];
rz(-0.5378546) q[3];
sx q[3];
rz(-0.23306498) q[3];
sx q[3];
rz(-2.2190998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(0.16263738) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(0.45062137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-3.056884) q[0];
sx q[0];
rz(-1.9730294) q[0];
rz(-pi) q[1];
rz(1.1969823) q[2];
sx q[2];
rz(-1.8573559) q[2];
sx q[2];
rz(-2.8754004) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.140481) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(-0.66092296) q[1];
x q[2];
rz(-2.8353325) q[3];
sx q[3];
rz(-0.89950409) q[3];
sx q[3];
rz(1.1989532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11015636) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(1.1412096) q[0];
rz(1.2966688) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(-2.8362714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13320623) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(0.58831373) q[0];
rz(2.4112169) q[2];
sx q[2];
rz(-2.5502898) q[2];
sx q[2];
rz(3.0808466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2480337) q[1];
sx q[1];
rz(-0.78130022) q[1];
sx q[1];
rz(-1.4206738) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0852854) q[3];
sx q[3];
rz(-1.1122276) q[3];
sx q[3];
rz(1.6955299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(0.34710458) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(0.40571037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203128) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-0.41723199) q[0];
x q[1];
rz(-1.1666127) q[2];
sx q[2];
rz(-2.5378072) q[2];
sx q[2];
rz(-0.052978901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2629562) q[1];
sx q[1];
rz(-2.7107781) q[1];
sx q[1];
rz(-2.2104435) q[1];
rz(2.3769576) q[3];
sx q[3];
rz(-1.7406165) q[3];
sx q[3];
rz(0.09165435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(-1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(-0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(1.8242594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2989527) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(-1.0277261) q[0];
rz(-0.66770422) q[2];
sx q[2];
rz(-1.0945079) q[2];
sx q[2];
rz(-0.27031937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6588905) q[1];
sx q[1];
rz(-1.5679949) q[1];
sx q[1];
rz(2.3059694) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.156593) q[3];
sx q[3];
rz(-1.9337774) q[3];
sx q[3];
rz(-1.3604878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(-2.2463634) q[2];
rz(-1.9334531) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(-2.0136925) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989482) q[0];
sx q[0];
rz(-2.0751187) q[0];
sx q[0];
rz(1.9796611) q[0];
x q[1];
rz(-2.5008051) q[2];
sx q[2];
rz(-1.0377585) q[2];
sx q[2];
rz(-1.8719045) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5269446) q[1];
sx q[1];
rz(-1.6741747) q[1];
sx q[1];
rz(0.23067109) q[1];
x q[2];
rz(-0.021846847) q[3];
sx q[3];
rz(-2.1751746) q[3];
sx q[3];
rz(-1.5165839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82905519) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(2.5404239) q[2];
rz(2.6230295) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.4275838) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(2.1144497) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(-2.5491319) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6153529) q[0];
sx q[0];
rz(-2.549299) q[0];
sx q[0];
rz(-3.0583529) q[0];
rz(-pi) q[1];
rz(-2.0855037) q[2];
sx q[2];
rz(-2.0851496) q[2];
sx q[2];
rz(0.27013847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29443024) q[1];
sx q[1];
rz(-1.2420328) q[1];
sx q[1];
rz(0.86818236) q[1];
rz(1.1660277) q[3];
sx q[3];
rz(-1.824607) q[3];
sx q[3];
rz(-0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(-0.44211659) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.3046718) q[0];
rz(1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(2.5794199) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5931372) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(-2.0329342) q[0];
rz(0.093900605) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(1.7319861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1162029) q[1];
sx q[1];
rz(-0.89239489) q[1];
sx q[1];
rz(-0.62979001) q[1];
rz(-pi) q[2];
rz(2.7355248) q[3];
sx q[3];
rz(-1.9417524) q[3];
sx q[3];
rz(3.1379116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(-0.98158681) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(-0.04709588) q[0];
rz(-1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809966) q[0];
sx q[0];
rz(-0.91616917) q[0];
sx q[0];
rz(-3.1050443) q[0];
rz(1.677703) q[2];
sx q[2];
rz(-2.6968743) q[2];
sx q[2];
rz(2.5267548) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.88356) q[1];
sx q[1];
rz(-1.1415392) q[1];
sx q[1];
rz(-0.92991035) q[1];
rz(-pi) q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(0.88276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3180278) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(-1.5927096) q[2];
sx q[2];
rz(-2.1219325) q[2];
sx q[2];
rz(1.3757201) q[2];
rz(1.1669284) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];