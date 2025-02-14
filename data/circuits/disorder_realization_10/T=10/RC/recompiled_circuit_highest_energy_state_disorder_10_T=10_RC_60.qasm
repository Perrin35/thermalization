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
rz(-2.7918474) q[0];
sx q[0];
rz(-2.1748769) q[0];
sx q[0];
rz(2.7803687) q[0];
rz(-0.53520441) q[1];
sx q[1];
rz(-1.1517795) q[1];
sx q[1];
rz(-1.3020383) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0143876) q[0];
sx q[0];
rz(-2.706658) q[0];
sx q[0];
rz(2.4282703) q[0];
rz(-pi) q[1];
rz(-1.6066685) q[2];
sx q[2];
rz(-2.2891993) q[2];
sx q[2];
rz(1.0234229) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8201789) q[1];
sx q[1];
rz(-0.4785236) q[1];
sx q[1];
rz(2.6027563) q[1];
rz(-pi) q[2];
rz(0.31631542) q[3];
sx q[3];
rz(-0.22238734) q[3];
sx q[3];
rz(2.5393644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5754622) q[2];
sx q[2];
rz(-2.8464816) q[2];
sx q[2];
rz(2.5556514) q[2];
rz(1.5015886) q[3];
sx q[3];
rz(-1.8668886) q[3];
sx q[3];
rz(1.2007825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457299) q[0];
sx q[0];
rz(-2.0450617) q[0];
sx q[0];
rz(1.1281661) q[0];
rz(2.3919171) q[1];
sx q[1];
rz(-1.3355037) q[1];
sx q[1];
rz(-1.4834167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7455341) q[0];
sx q[0];
rz(-1.1107971) q[0];
sx q[0];
rz(-2.3961623) q[0];
x q[1];
rz(-1.9374766) q[2];
sx q[2];
rz(-1.0016236) q[2];
sx q[2];
rz(1.8682084) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.0644599) q[1];
sx q[1];
rz(-1.049548) q[1];
sx q[1];
rz(1.8690994) q[1];
rz(0.34134369) q[3];
sx q[3];
rz(-2.0273547) q[3];
sx q[3];
rz(-1.2767291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7538606) q[2];
sx q[2];
rz(-0.78247491) q[2];
sx q[2];
rz(2.1369047) q[2];
rz(-0.0082155148) q[3];
sx q[3];
rz(-1.7693844) q[3];
sx q[3];
rz(2.7042702) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32366556) q[0];
sx q[0];
rz(-1.8765457) q[0];
sx q[0];
rz(1.0311968) q[0];
rz(0.19510706) q[1];
sx q[1];
rz(-0.61385265) q[1];
sx q[1];
rz(0.88417792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613286) q[0];
sx q[0];
rz(-1.5710281) q[0];
sx q[0];
rz(3.1383297) q[0];
rz(-pi) q[1];
rz(0.45978893) q[2];
sx q[2];
rz(-1.7722691) q[2];
sx q[2];
rz(0.10290621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9926874) q[1];
sx q[1];
rz(-2.7200448) q[1];
sx q[1];
rz(-1.804638) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8665119) q[3];
sx q[3];
rz(-2.1397192) q[3];
sx q[3];
rz(3.0771098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6241793) q[2];
sx q[2];
rz(-1.2531345) q[2];
sx q[2];
rz(-2.6873028) q[2];
rz(1.0034466) q[3];
sx q[3];
rz(-2.5277621) q[3];
sx q[3];
rz(1.7295037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91785947) q[0];
sx q[0];
rz(-2.7847325) q[0];
sx q[0];
rz(-0.80765635) q[0];
rz(1.3937021) q[1];
sx q[1];
rz(-1.7180157) q[1];
sx q[1];
rz(1.4180988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7515441) q[0];
sx q[0];
rz(-1.1466197) q[0];
sx q[0];
rz(-0.15985711) q[0];
x q[1];
rz(-1.3803064) q[2];
sx q[2];
rz(-0.87403008) q[2];
sx q[2];
rz(2.246169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6114072) q[1];
sx q[1];
rz(-0.3819335) q[1];
sx q[1];
rz(0.032738222) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7341081) q[3];
sx q[3];
rz(-2.0563246) q[3];
sx q[3];
rz(0.21903506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35004804) q[2];
sx q[2];
rz(-1.2309265) q[2];
sx q[2];
rz(-3.0080504) q[2];
rz(-2.7327025) q[3];
sx q[3];
rz(-3.0372527) q[3];
sx q[3];
rz(-2.1321808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23311663) q[0];
sx q[0];
rz(-2.3148843) q[0];
sx q[0];
rz(-2.6772461) q[0];
rz(2.6404479) q[1];
sx q[1];
rz(-0.25096384) q[1];
sx q[1];
rz(-0.73890013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66128816) q[0];
sx q[0];
rz(-1.6853852) q[0];
sx q[0];
rz(-0.25729723) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54292572) q[2];
sx q[2];
rz(-1.6871337) q[2];
sx q[2];
rz(0.1067208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6433324) q[1];
sx q[1];
rz(-2.6427489) q[1];
sx q[1];
rz(-2.6698334) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.471354) q[3];
sx q[3];
rz(-0.76858187) q[3];
sx q[3];
rz(-2.0266899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.240694) q[2];
sx q[2];
rz(-1.7186761) q[2];
sx q[2];
rz(0.23441976) q[2];
rz(-2.8068986) q[3];
sx q[3];
rz(-0.60939279) q[3];
sx q[3];
rz(1.7083683) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1471106) q[0];
sx q[0];
rz(-2.2587977) q[0];
sx q[0];
rz(0.91142803) q[0];
rz(1.8270739) q[1];
sx q[1];
rz(-0.85845566) q[1];
sx q[1];
rz(1.3988769) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11451498) q[0];
sx q[0];
rz(-2.0721779) q[0];
sx q[0];
rz(-1.03427) q[0];
rz(-pi) q[1];
rz(-3.0023742) q[2];
sx q[2];
rz(-1.5409383) q[2];
sx q[2];
rz(-0.6767737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8533852) q[1];
sx q[1];
rz(-2.7010272) q[1];
sx q[1];
rz(-1.8024995) q[1];
x q[2];
rz(-0.26246385) q[3];
sx q[3];
rz(-1.0499716) q[3];
sx q[3];
rz(2.9948522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92480245) q[2];
sx q[2];
rz(-0.81764597) q[2];
sx q[2];
rz(1.7693046) q[2];
rz(-1.4745845) q[3];
sx q[3];
rz(-1.0627012) q[3];
sx q[3];
rz(0.32349989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.40397662) q[0];
sx q[0];
rz(-2.1196899) q[0];
sx q[0];
rz(1.2603941) q[0];
rz(0.34126392) q[1];
sx q[1];
rz(-2.212846) q[1];
sx q[1];
rz(-2.0492679) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4173146) q[0];
sx q[0];
rz(-1.5515741) q[0];
sx q[0];
rz(2.6249983) q[0];
rz(-pi) q[1];
rz(0.34185648) q[2];
sx q[2];
rz(-1.3876788) q[2];
sx q[2];
rz(1.3951071) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2637529) q[1];
sx q[1];
rz(-1.8108551) q[1];
sx q[1];
rz(2.2113483) q[1];
x q[2];
rz(1.8856005) q[3];
sx q[3];
rz(-1.971759) q[3];
sx q[3];
rz(2.8481066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1616538) q[2];
sx q[2];
rz(-1.7540163) q[2];
sx q[2];
rz(0.8650583) q[2];
rz(-0.093552731) q[3];
sx q[3];
rz(-1.716194) q[3];
sx q[3];
rz(-3.0430999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.426067) q[0];
sx q[0];
rz(-0.97398296) q[0];
sx q[0];
rz(1.4564212) q[0];
rz(2.9086225) q[1];
sx q[1];
rz(-1.3825682) q[1];
sx q[1];
rz(1.2196541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15773179) q[0];
sx q[0];
rz(-1.2810575) q[0];
sx q[0];
rz(-0.74037551) q[0];
x q[1];
rz(-2.4039335) q[2];
sx q[2];
rz(-0.56471497) q[2];
sx q[2];
rz(1.6961404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0387242) q[1];
sx q[1];
rz(-1.7448493) q[1];
sx q[1];
rz(-2.7636488) q[1];
x q[2];
rz(2.574) q[3];
sx q[3];
rz(-2.360968) q[3];
sx q[3];
rz(-2.7989557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0435656) q[2];
sx q[2];
rz(-1.4412014) q[2];
sx q[2];
rz(-3.0230057) q[2];
rz(-2.3218527) q[3];
sx q[3];
rz(-0.44410646) q[3];
sx q[3];
rz(-2.8290101) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239301) q[0];
sx q[0];
rz(-0.95159641) q[0];
sx q[0];
rz(-2.1753525) q[0];
rz(0.36695925) q[1];
sx q[1];
rz(-0.96261135) q[1];
sx q[1];
rz(-0.56698322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0006994572) q[0];
sx q[0];
rz(-0.53017925) q[0];
sx q[0];
rz(2.7383366) q[0];
rz(-2.5091293) q[2];
sx q[2];
rz(-2.2603432) q[2];
sx q[2];
rz(2.0105711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7427269) q[1];
sx q[1];
rz(-2.4165771) q[1];
sx q[1];
rz(-2.8319977) q[1];
x q[2];
rz(1.1375327) q[3];
sx q[3];
rz(-0.71914395) q[3];
sx q[3];
rz(2.5844203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6734267) q[2];
sx q[2];
rz(-1.6750853) q[2];
sx q[2];
rz(1.0183938) q[2];
rz(-2.9112725) q[3];
sx q[3];
rz(-2.6789594) q[3];
sx q[3];
rz(-0.066430062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224943) q[0];
sx q[0];
rz(-1.007217) q[0];
sx q[0];
rz(-1.8100716) q[0];
rz(-1.9271756) q[1];
sx q[1];
rz(-1.4446832) q[1];
sx q[1];
rz(0.4164947) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40503026) q[0];
sx q[0];
rz(-0.89823898) q[0];
sx q[0];
rz(0.26215464) q[0];
rz(1.3168066) q[2];
sx q[2];
rz(-1.87623) q[2];
sx q[2];
rz(2.0485725) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.107376) q[1];
sx q[1];
rz(-1.9946846) q[1];
sx q[1];
rz(0.54383212) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0216398) q[3];
sx q[3];
rz(-1.4792253) q[3];
sx q[3];
rz(-1.6710323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.74849621) q[2];
sx q[2];
rz(-1.2747719) q[2];
sx q[2];
rz(1.3416802) q[2];
rz(-1.120535) q[3];
sx q[3];
rz(-1.7399961) q[3];
sx q[3];
rz(3.024658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8244767) q[0];
sx q[0];
rz(-1.308029) q[0];
sx q[0];
rz(-1.8608004) q[0];
rz(2.169213) q[1];
sx q[1];
rz(-2.0449816) q[1];
sx q[1];
rz(2.4349946) q[1];
rz(-1.8351716) q[2];
sx q[2];
rz(-0.96294667) q[2];
sx q[2];
rz(1.3542702) q[2];
rz(-1.6617358) q[3];
sx q[3];
rz(-1.6030238) q[3];
sx q[3];
rz(0.054790767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
