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
rz(-1.3402101) q[0];
sx q[0];
rz(-2.8647381) q[0];
sx q[0];
rz(1.0387596) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(3.9781456) q[1];
sx q[1];
rz(9.8415924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2622803) q[0];
sx q[0];
rz(-1.4646155) q[0];
sx q[0];
rz(-2.4152629) q[0];
rz(-pi) q[1];
rz(2.7567467) q[2];
sx q[2];
rz(-0.59078465) q[2];
sx q[2];
rz(-2.7123775) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8326679) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(-0.69242386) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9454222) q[3];
sx q[3];
rz(-2.9769398) q[3];
sx q[3];
rz(-0.97033721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(-1.2539585) q[2];
rz(-1.5597255) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9480243) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(-0.76876202) q[0];
rz(-1.7747152) q[1];
sx q[1];
rz(-1.3221075) q[1];
sx q[1];
rz(2.1034525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4878142) q[0];
sx q[0];
rz(-0.013049203) q[0];
sx q[0];
rz(-2.292146) q[0];
x q[1];
rz(-1.9351472) q[2];
sx q[2];
rz(-2.2585105) q[2];
sx q[2];
rz(-1.5503413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1662054) q[1];
sx q[1];
rz(-2.176099) q[1];
sx q[1];
rz(2.0209842) q[1];
x q[2];
rz(-1.5350545) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(-1.8349748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64608964) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(2.3213279) q[2];
rz(-2.8881554) q[3];
sx q[3];
rz(-0.35496747) q[3];
sx q[3];
rz(2.7804815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.8609817) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(-0.88731998) q[0];
rz(2.6112828) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(1.1393772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.087505) q[0];
sx q[0];
rz(-1.2257135) q[0];
sx q[0];
rz(1.4201384) q[0];
x q[1];
rz(-1.0372889) q[2];
sx q[2];
rz(-0.59743687) q[2];
sx q[2];
rz(-0.00053279579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85998409) q[1];
sx q[1];
rz(-1.8164595) q[1];
sx q[1];
rz(1.5558234) q[1];
x q[2];
rz(-3.0501796) q[3];
sx q[3];
rz(-2.6210636) q[3];
sx q[3];
rz(0.5239858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.016971074) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-0.43928453) q[2];
rz(2.6708421) q[3];
sx q[3];
rz(-0.079340383) q[3];
sx q[3];
rz(1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.73862326) q[0];
sx q[0];
rz(-0.93337494) q[0];
sx q[0];
rz(-0.11548197) q[0];
rz(-0.11784095) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(1.3062564) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7457434) q[0];
sx q[0];
rz(-1.1538528) q[0];
sx q[0];
rz(-0.31503079) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4815062) q[2];
sx q[2];
rz(-1.266511) q[2];
sx q[2];
rz(3.0344998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0178284) q[1];
sx q[1];
rz(-2.1514847) q[1];
sx q[1];
rz(2.3821217) q[1];
rz(-pi) q[2];
x q[2];
rz(3.066411) q[3];
sx q[3];
rz(-2.1422221) q[3];
sx q[3];
rz(2.927185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8370886) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(-0.15466776) q[2];
rz(-1.539544) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(-2.535517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57347572) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(-0.83531761) q[0];
rz(2.4256445) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(3.0221525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7258541) q[0];
sx q[0];
rz(-1.5149759) q[0];
sx q[0];
rz(-3.124696) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7238462) q[2];
sx q[2];
rz(-1.1259534) q[2];
sx q[2];
rz(2.9314624) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99239527) q[1];
sx q[1];
rz(-2.1572504) q[1];
sx q[1];
rz(1.0066562) q[1];
rz(-2.4923513) q[3];
sx q[3];
rz(-2.5045145) q[3];
sx q[3];
rz(2.0811045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2401838) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(3.0933464) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(-0.72859305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7399087) q[0];
sx q[0];
rz(-0.70867276) q[0];
sx q[0];
rz(2.0106864) q[0];
rz(0.4785969) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(-1.4071677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78188183) q[0];
sx q[0];
rz(-0.053584307) q[0];
sx q[0];
rz(-0.028938541) q[0];
rz(-pi) q[1];
rz(2.8596609) q[2];
sx q[2];
rz(-2.273186) q[2];
sx q[2];
rz(-1.3739746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4987405) q[1];
sx q[1];
rz(-1.3428709) q[1];
sx q[1];
rz(0.1164611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7516486) q[3];
sx q[3];
rz(-1.3612124) q[3];
sx q[3];
rz(-1.4246032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5799134) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(1.178406) q[2];
rz(2.0922349) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2105763) q[0];
sx q[0];
rz(-1.1897621) q[0];
sx q[0];
rz(1.3354906) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(2.8335422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9481407) q[0];
sx q[0];
rz(-1.7823671) q[0];
sx q[0];
rz(-2.7443462) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1739866) q[2];
sx q[2];
rz(-2.0521297) q[2];
sx q[2];
rz(0.86908841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1738893) q[1];
sx q[1];
rz(-0.99882616) q[1];
sx q[1];
rz(2.553945) q[1];
x q[2];
rz(3.0940476) q[3];
sx q[3];
rz(-1.1851839) q[3];
sx q[3];
rz(-1.6507208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7293952) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(-3.0103053) q[2];
rz(0.73733759) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8312296) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(-0.53037733) q[0];
rz(-1.7165548) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(0.72149611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24151267) q[0];
sx q[0];
rz(-0.78966138) q[0];
sx q[0];
rz(-2.4898743) q[0];
rz(1.4230096) q[2];
sx q[2];
rz(-2.035454) q[2];
sx q[2];
rz(0.3275268) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3949273) q[1];
sx q[1];
rz(-0.29005602) q[1];
sx q[1];
rz(2.2941089) q[1];
rz(-pi) q[2];
rz(-0.78254487) q[3];
sx q[3];
rz(-1.492332) q[3];
sx q[3];
rz(2.9786547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7044652) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.4777769) q[2];
rz(1.1635228) q[3];
sx q[3];
rz(-1.082837) q[3];
sx q[3];
rz(-1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.32996938) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(0.60923088) q[0];
rz(-1.5513264) q[1];
sx q[1];
rz(-2.8158999) q[1];
sx q[1];
rz(-1.685166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44269366) q[0];
sx q[0];
rz(-0.44166587) q[0];
sx q[0];
rz(-2.1864258) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6563719) q[2];
sx q[2];
rz(-1.5169797) q[2];
sx q[2];
rz(-0.14419989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(1.2213329) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2595657) q[3];
sx q[3];
rz(-2.2891364) q[3];
sx q[3];
rz(0.73757987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7234601) q[2];
sx q[2];
rz(-0.89833608) q[2];
sx q[2];
rz(-2.2612803) q[2];
rz(2.9355925) q[3];
sx q[3];
rz(-2.7166631) q[3];
sx q[3];
rz(2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708165) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(0.01195512) q[0];
rz(-1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(2.5451122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76145455) q[0];
sx q[0];
rz(-1.5517457) q[0];
sx q[0];
rz(0.46940243) q[0];
x q[1];
rz(-1.311074) q[2];
sx q[2];
rz(-1.2190281) q[2];
sx q[2];
rz(-0.042366926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47190753) q[1];
sx q[1];
rz(-2.0457256) q[1];
sx q[1];
rz(-3.0529939) q[1];
x q[2];
rz(-1.7285198) q[3];
sx q[3];
rz(-1.8650569) q[3];
sx q[3];
rz(2.074948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(0.19700024) q[2];
rz(1.9893076) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(-2.6939189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780509) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(2.9327783) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(-0.094194407) q[2];
sx q[2];
rz(-0.8834367) q[2];
sx q[2];
rz(1.5987664) q[2];
rz(0.87416762) q[3];
sx q[3];
rz(-0.71577358) q[3];
sx q[3];
rz(-1.7795455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
