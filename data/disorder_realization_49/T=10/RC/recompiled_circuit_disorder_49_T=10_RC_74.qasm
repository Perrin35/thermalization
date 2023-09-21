OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(-0.0052069081) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(-1.985328) q[1];
sx q[1];
rz(-1.1896689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0931041) q[0];
sx q[0];
rz(-2.5079873) q[0];
sx q[0];
rz(0.60237576) q[0];
x q[1];
rz(2.4742545) q[2];
sx q[2];
rz(-0.22775209) q[2];
sx q[2];
rz(-1.8337133) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78566879) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(-2.0211401) q[1];
rz(-0.40640229) q[3];
sx q[3];
rz(-2.1592525) q[3];
sx q[3];
rz(1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(-2.544196) q[2];
rz(1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(0.11322583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7069488) q[0];
sx q[0];
rz(-1.5320677) q[0];
sx q[0];
rz(1.1452655) q[0];
rz(-0.073268854) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(-0.71436963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3791703) q[1];
sx q[1];
rz(-1.1602853) q[1];
sx q[1];
rz(-1.1344086) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96666386) q[3];
sx q[3];
rz(-1.3738487) q[3];
sx q[3];
rz(0.40991022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(1.0473898) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(0.12761322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022284431) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(0.22247252) q[0];
x q[1];
rz(-2.2682842) q[2];
sx q[2];
rz(-2.1224988) q[2];
sx q[2];
rz(1.1683977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4096654) q[1];
sx q[1];
rz(-1.8095784) q[1];
sx q[1];
rz(3.0569397) q[1];
x q[2];
rz(1.913194) q[3];
sx q[3];
rz(-2.4857593) q[3];
sx q[3];
rz(-2.0369903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(-0.37240949) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62896699) q[0];
sx q[0];
rz(-2.2722639) q[0];
sx q[0];
rz(0.15928282) q[0];
x q[1];
rz(-0.98948688) q[2];
sx q[2];
rz(-1.8539691) q[2];
sx q[2];
rz(-2.1894933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73788961) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(1.9562734) q[1];
rz(-0.59471547) q[3];
sx q[3];
rz(-0.82154951) q[3];
sx q[3];
rz(-2.6298414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-0.4549543) q[3];
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
rz(1.882778) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(1.1313261) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(2.4647443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35447435) q[0];
sx q[0];
rz(-1.1312383) q[0];
sx q[0];
rz(-1.6580824) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3555688) q[2];
sx q[2];
rz(-0.85126801) q[2];
sx q[2];
rz(2.1446251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0332898) q[1];
sx q[1];
rz(-1.3160333) q[1];
sx q[1];
rz(1.5974031) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058651794) q[3];
sx q[3];
rz(-0.39415112) q[3];
sx q[3];
rz(-2.4413721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(-2.3964264) q[1];
sx q[1];
rz(-2.7909653) q[1];
sx q[1];
rz(-2.2084592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720848) q[0];
sx q[0];
rz(-0.55969319) q[0];
sx q[0];
rz(-2.1493388) q[0];
rz(-pi) q[1];
rz(-2.6283693) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(-2.9483587) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3999403) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(1.2029349) q[1];
rz(1.798614) q[3];
sx q[3];
rz(-2.4711547) q[3];
sx q[3];
rz(-0.05712856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(0.28665001) q[2];
rz(0.24946985) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(-0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(-1.6749143) q[0];
rz(2.1215227) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2350378) q[0];
sx q[0];
rz(-0.28309238) q[0];
sx q[0];
rz(-1.4366158) q[0];
rz(-pi) q[1];
rz(-2.8011488) q[2];
sx q[2];
rz(-1.2203487) q[2];
sx q[2];
rz(-3.0710789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96879362) q[1];
sx q[1];
rz(-1.4152923) q[1];
sx q[1];
rz(-1.1771727) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7822958) q[3];
sx q[3];
rz(-2.2106417) q[3];
sx q[3];
rz(2.0169472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(0.83759585) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(2.7517095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1471257) q[0];
sx q[0];
rz(-1.9305221) q[0];
sx q[0];
rz(0.61867449) q[0];
rz(-pi) q[1];
rz(0.6673442) q[2];
sx q[2];
rz(-2.5255425) q[2];
sx q[2];
rz(-1.4916071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0583565) q[1];
sx q[1];
rz(-0.55587686) q[1];
sx q[1];
rz(0.91169375) q[1];
x q[2];
rz(0.26384683) q[3];
sx q[3];
rz(-0.90875328) q[3];
sx q[3];
rz(-2.2784233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3020246) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(-0.8977302) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(-0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927239) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(2.8651967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24647507) q[0];
sx q[0];
rz(-0.44024375) q[0];
sx q[0];
rz(-2.759139) q[0];
rz(-pi) q[1];
rz(-2.9133965) q[2];
sx q[2];
rz(-1.3820634) q[2];
sx q[2];
rz(-0.91918321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.285187) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(1.3473947) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8180088) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(-1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(0.35274831) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(2.1113077) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477662) q[0];
sx q[0];
rz(-0.86137912) q[0];
sx q[0];
rz(-1.8295656) q[0];
rz(1.9039541) q[2];
sx q[2];
rz(-0.91234708) q[2];
sx q[2];
rz(1.2308987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5926338) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(1.2195107) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61334765) q[3];
sx q[3];
rz(-1.8457335) q[3];
sx q[3];
rz(3.1230694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(0.34118787) q[2];
rz(0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-3.0497131) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13070233) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(2.4717992) q[1];
sx q[1];
rz(-2.5302946) q[1];
sx q[1];
rz(-1.5763462) q[1];
rz(1.7651991) q[2];
sx q[2];
rz(-1.9971078) q[2];
sx q[2];
rz(2.8304921) q[2];
rz(2.1445198) q[3];
sx q[3];
rz(-0.96521796) q[3];
sx q[3];
rz(-0.60346606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
