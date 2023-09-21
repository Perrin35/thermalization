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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22098955) q[0];
sx q[0];
rz(-1.594829) q[0];
sx q[0];
rz(0.8997957) q[0];
x q[1];
rz(-2.5976074) q[2];
sx q[2];
rz(-1.7684801) q[2];
sx q[2];
rz(2.939784) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3687009) q[1];
sx q[1];
rz(-1.643631) q[1];
sx q[1];
rz(1.1198977) q[1];
rz(-pi) q[2];
rz(-2.4281265) q[3];
sx q[3];
rz(-2.2029075) q[3];
sx q[3];
rz(-1.6085094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67316002) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(2.6426278) q[2];
rz(0.4237825) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(-2.965062) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(0.36546779) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27643833) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(-0.47238484) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89895504) q[2];
sx q[2];
rz(-1.0869173) q[2];
sx q[2];
rz(-1.3324141) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6668015) q[1];
sx q[1];
rz(-2.431369) q[1];
sx q[1];
rz(-0.53479654) q[1];
x q[2];
rz(-0.20110735) q[3];
sx q[3];
rz(-1.4521986) q[3];
sx q[3];
rz(-0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-2.0642521) q[2];
rz(0.16263738) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(-0.61082947) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(0.45062137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.492309) q[0];
sx q[0];
rz(-1.5376687) q[0];
sx q[0];
rz(1.4928198) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1969823) q[2];
sx q[2];
rz(-1.2842368) q[2];
sx q[2];
rz(0.26619226) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0011117) q[1];
sx q[1];
rz(-1.8794267) q[1];
sx q[1];
rz(-0.66092296) q[1];
rz(-pi) q[2];
rz(-0.87617989) q[3];
sx q[3];
rz(-1.3324705) q[3];
sx q[3];
rz(0.5660457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(2.2971161) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-0.30532125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0083864) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(0.58831373) q[0];
rz(-1.1496454) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(2.2575833) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42974737) q[1];
sx q[1];
rz(-1.6763121) q[1];
sx q[1];
rz(-0.79515102) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3577945) q[3];
sx q[3];
rz(-2.4664306) q[3];
sx q[3];
rz(-2.6019707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(-0.34710458) q[2];
rz(-0.38337213) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(-2.7358823) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.273542) q[0];
sx q[0];
rz(-1.8450071) q[0];
sx q[0];
rz(-2.4540841) q[0];
x q[1];
rz(-0.26486764) q[2];
sx q[2];
rz(-1.0215853) q[2];
sx q[2];
rz(-0.4262281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5647445) q[1];
sx q[1];
rz(-1.2291359) q[1];
sx q[1];
rz(-2.8738352) q[1];
x q[2];
rz(2.898786) q[3];
sx q[3];
rz(-2.3620776) q[3];
sx q[3];
rz(-1.3047583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.4542788) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(1.3173332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.84264) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(-1.0277261) q[0];
rz(-0.66770422) q[2];
sx q[2];
rz(-1.0945079) q[2];
sx q[2];
rz(-0.27031937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4827022) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(2.3059694) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98499961) q[3];
sx q[3];
rz(-1.2078152) q[3];
sx q[3];
rz(-1.7811048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2736149) q[0];
sx q[0];
rz(-0.6379188) q[0];
sx q[0];
rz(2.5173553) q[0];
x q[1];
rz(0.77881323) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(-0.29701172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1612138) q[1];
sx q[1];
rz(-1.3413789) q[1];
sx q[1];
rz(1.4646261) q[1];
x q[2];
rz(-0.96630649) q[3];
sx q[3];
rz(-1.55282) q[3];
sx q[3];
rz(3.0997961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3125375) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475875) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(3.0684379) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-2.5491319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1661467) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(-0.59068824) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0855037) q[2];
sx q[2];
rz(-1.0564431) q[2];
sx q[2];
rz(2.8714542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2295099) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(2.0565226) q[1];
rz(-pi) q[2];
x q[2];
rz(2.153272) q[3];
sx q[3];
rz(-2.6675825) q[3];
sx q[3];
rz(1.43309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(3.1282848) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028037926) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.8369209) q[0];
rz(1.6944983) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87555128) q[0];
sx q[0];
rz(-2.3444359) q[0];
sx q[0];
rz(0.50811998) q[0];
rz(-2.1015342) q[2];
sx q[2];
rz(-1.6518403) q[2];
sx q[2];
rz(0.20866742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8852946) q[1];
sx q[1];
rz(-0.89026272) q[1];
sx q[1];
rz(0.93974944) q[1];
x q[2];
rz(-2.3638944) q[3];
sx q[3];
rz(-0.54293984) q[3];
sx q[3];
rz(-0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3146064) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(1.3866407) q[1];
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
rz(2.3536516) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(0.91584648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1282975) q[2];
sx q[2];
rz(-1.5248761) q[2];
sx q[2];
rz(-2.0890582) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25803265) q[1];
sx q[1];
rz(-1.1415392) q[1];
sx q[1];
rz(0.92991035) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(-2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64426595) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(-2.3406773) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(-0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3180278) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(0.50280747) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(-1.5927096) q[2];
sx q[2];
rz(-2.1219325) q[2];
sx q[2];
rz(1.3757201) q[2];
rz(-2.3021163) q[3];
sx q[3];
rz(-0.52121938) q[3];
sx q[3];
rz(-2.2073707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
