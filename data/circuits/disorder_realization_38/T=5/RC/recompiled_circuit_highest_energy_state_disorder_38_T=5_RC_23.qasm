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
rz(0.44904798) q[0];
sx q[0];
rz(4.2807978) q[0];
sx q[0];
rz(9.8674404) q[0];
rz(0.17065419) q[1];
sx q[1];
rz(3.9332665) q[1];
sx q[1];
rz(10.153833) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1478445) q[0];
sx q[0];
rz(-1.1708852) q[0];
sx q[0];
rz(0.73109674) q[0];
x q[1];
rz(-0.0021931113) q[2];
sx q[2];
rz(-0.38972112) q[2];
sx q[2];
rz(0.8726697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2393136) q[1];
sx q[1];
rz(-2.2444176) q[1];
sx q[1];
rz(0.31142733) q[1];
rz(0.40790848) q[3];
sx q[3];
rz(-0.98263697) q[3];
sx q[3];
rz(1.5666636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2810104) q[2];
sx q[2];
rz(-0.082110114) q[2];
sx q[2];
rz(-2.3561467) q[2];
rz(0.9663409) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53168374) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(0.63572836) q[0];
rz(2.5081778) q[1];
sx q[1];
rz(-0.76786357) q[1];
sx q[1];
rz(0.33087081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8757056) q[0];
sx q[0];
rz(-0.71015753) q[0];
sx q[0];
rz(-2.6214834) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.026658807) q[2];
sx q[2];
rz(-2.0377198) q[2];
sx q[2];
rz(-0.24031642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91869044) q[1];
sx q[1];
rz(-1.6044334) q[1];
sx q[1];
rz(-1.5123537) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5832354) q[3];
sx q[3];
rz(-1.6553989) q[3];
sx q[3];
rz(-1.0983262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.50515455) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(2.6804697) q[2];
rz(2.4165706) q[3];
sx q[3];
rz(-0.45261639) q[3];
sx q[3];
rz(-2.4729474) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.649699) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(-0.53832501) q[0];
rz(-0.0093983924) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(1.9422772) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56012404) q[0];
sx q[0];
rz(-2.631979) q[0];
sx q[0];
rz(1.6625924) q[0];
rz(0.17260562) q[2];
sx q[2];
rz(-2.4510018) q[2];
sx q[2];
rz(-2.0334854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6965683) q[1];
sx q[1];
rz(-2.4557132) q[1];
sx q[1];
rz(0.98937757) q[1];
rz(-pi) q[2];
rz(1.2869419) q[3];
sx q[3];
rz(-0.7470658) q[3];
sx q[3];
rz(-3.0892854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37375307) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(0.72244942) q[2];
rz(-0.59822285) q[3];
sx q[3];
rz(-3.1066419) q[3];
sx q[3];
rz(-1.2178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55906051) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(0.019056888) q[0];
rz(-1.6825153) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(-0.51024514) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048240926) q[0];
sx q[0];
rz(-0.66931319) q[0];
sx q[0];
rz(1.1360243) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8809659) q[2];
sx q[2];
rz(-2.2686743) q[2];
sx q[2];
rz(0.073402799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45106798) q[1];
sx q[1];
rz(-2.1803586) q[1];
sx q[1];
rz(-1.068348) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41401569) q[3];
sx q[3];
rz(-1.7018205) q[3];
sx q[3];
rz(1.6725292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14429188) q[2];
sx q[2];
rz(-1.68953) q[2];
sx q[2];
rz(-0.22225456) q[2];
rz(-0.079515919) q[3];
sx q[3];
rz(-0.54602081) q[3];
sx q[3];
rz(-2.5047746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6154489) q[0];
sx q[0];
rz(-2.0455102) q[0];
sx q[0];
rz(-1.8002321) q[0];
rz(-2.6306131) q[1];
sx q[1];
rz(-1.7781517) q[1];
sx q[1];
rz(-2.708639) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110033) q[0];
sx q[0];
rz(-1.827914) q[0];
sx q[0];
rz(-1.6297618) q[0];
rz(-0.21375653) q[2];
sx q[2];
rz(-1.9659489) q[2];
sx q[2];
rz(1.1590126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6420685) q[1];
sx q[1];
rz(-0.41710258) q[1];
sx q[1];
rz(-0.92315556) q[1];
rz(-pi) q[2];
rz(-0.28750395) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(2.7903729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0114835) q[2];
sx q[2];
rz(-0.92769647) q[2];
sx q[2];
rz(2.1437272) q[2];
rz(0.70895553) q[3];
sx q[3];
rz(-0.38781375) q[3];
sx q[3];
rz(-0.97878218) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82584941) q[0];
sx q[0];
rz(-0.5873) q[0];
sx q[0];
rz(-0.89383268) q[0];
rz(2.1122487) q[1];
sx q[1];
rz(-0.7411595) q[1];
sx q[1];
rz(-3.0267874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5328767) q[0];
sx q[0];
rz(-2.8136557) q[0];
sx q[0];
rz(-0.84285835) q[0];
rz(-pi) q[1];
rz(1.9283781) q[2];
sx q[2];
rz(-0.96065694) q[2];
sx q[2];
rz(-1.4550191) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56267589) q[1];
sx q[1];
rz(-1.8496184) q[1];
sx q[1];
rz(1.2827966) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5070262) q[3];
sx q[3];
rz(-1.8511417) q[3];
sx q[3];
rz(-0.75322039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7041695) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(2.9278921) q[2];
rz(2.5050488) q[3];
sx q[3];
rz(-2.611219) q[3];
sx q[3];
rz(2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8146424) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(0.11006926) q[0];
rz(-0.07668177) q[1];
sx q[1];
rz(-0.86014599) q[1];
sx q[1];
rz(-1.1383879) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58041328) q[0];
sx q[0];
rz(-1.7903768) q[0];
sx q[0];
rz(-2.7835569) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78779548) q[2];
sx q[2];
rz(-1.5681364) q[2];
sx q[2];
rz(3.015369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36409196) q[1];
sx q[1];
rz(-1.0887676) q[1];
sx q[1];
rz(-2.9338323) q[1];
x q[2];
rz(-1.7848694) q[3];
sx q[3];
rz(-1.3297289) q[3];
sx q[3];
rz(1.9578707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25973213) q[2];
sx q[2];
rz(-1.0819819) q[2];
sx q[2];
rz(-0.38121769) q[2];
rz(0.11738736) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(-0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748902) q[0];
sx q[0];
rz(-2.369304) q[0];
sx q[0];
rz(0.49852398) q[0];
rz(2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(-0.097749762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5437336) q[0];
sx q[0];
rz(-1.2654788) q[0];
sx q[0];
rz(1.1174563) q[0];
rz(-pi) q[1];
rz(-0.4011863) q[2];
sx q[2];
rz(-2.7747535) q[2];
sx q[2];
rz(2.7374379) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6548582) q[1];
sx q[1];
rz(-2.2796101) q[1];
sx q[1];
rz(-0.72388078) q[1];
rz(-3.1322543) q[3];
sx q[3];
rz(-1.2492325) q[3];
sx q[3];
rz(3.0967979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30663124) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(2.5508896) q[2];
rz(-3.0513884) q[3];
sx q[3];
rz(-0.43961757) q[3];
sx q[3];
rz(-2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10839323) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(0.75476187) q[0];
rz(-2.8056878) q[1];
sx q[1];
rz(-0.55481189) q[1];
sx q[1];
rz(2.1051443) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514318) q[0];
sx q[0];
rz(-2.2519172) q[0];
sx q[0];
rz(-1.4195042) q[0];
rz(-pi) q[1];
rz(1.5912362) q[2];
sx q[2];
rz(-1.0896434) q[2];
sx q[2];
rz(2.1060973) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17286138) q[1];
sx q[1];
rz(-1.6097428) q[1];
sx q[1];
rz(-0.51677468) q[1];
rz(-pi) q[2];
rz(-1.6290954) q[3];
sx q[3];
rz(-0.23007904) q[3];
sx q[3];
rz(2.0571902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-0.45647344) q[2];
rz(-2.5392695) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(-2.203234) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8896821) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(2.67814) q[0];
rz(0.015333029) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-0.32036805) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8863618) q[0];
sx q[0];
rz(-1.4275274) q[0];
sx q[0];
rz(1.4930288) q[0];
x q[1];
rz(-1.3213083) q[2];
sx q[2];
rz(-1.2659327) q[2];
sx q[2];
rz(2.9314205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.55904222) q[1];
sx q[1];
rz(-0.96331412) q[1];
sx q[1];
rz(2.5431125) q[1];
rz(-pi) q[2];
rz(-0.69971754) q[3];
sx q[3];
rz(-0.43439242) q[3];
sx q[3];
rz(-0.35067973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.028712332) q[2];
sx q[2];
rz(-0.72673231) q[2];
sx q[2];
rz(-2.1325364) q[2];
rz(2.3446337) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(2.8033281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0924031) q[0];
sx q[0];
rz(-1.6272463) q[0];
sx q[0];
rz(1.8820681) q[0];
rz(-3.0567723) q[1];
sx q[1];
rz(-1.6592574) q[1];
sx q[1];
rz(1.4316373) q[1];
rz(1.3176805) q[2];
sx q[2];
rz(-0.24948487) q[2];
sx q[2];
rz(-2.7560888) q[2];
rz(-2.8390213) q[3];
sx q[3];
rz(-1.745001) q[3];
sx q[3];
rz(0.59636084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
