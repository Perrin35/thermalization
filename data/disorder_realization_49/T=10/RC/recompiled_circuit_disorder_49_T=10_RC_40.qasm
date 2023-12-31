OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(-1.9519238) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0484885) q[0];
sx q[0];
rz(-2.5079873) q[0];
sx q[0];
rz(0.60237576) q[0];
x q[1];
rz(-1.4283242) q[2];
sx q[2];
rz(-1.7490897) q[2];
sx q[2];
rz(-1.9879736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35808823) q[1];
sx q[1];
rz(-1.4231829) q[1];
sx q[1];
rz(-1.258177) q[1];
rz(-pi) q[2];
rz(1.0359997) q[3];
sx q[3];
rz(-2.4403604) q[3];
sx q[3];
rz(-1.086077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(1.3655837) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(-3.0283668) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43464381) q[0];
sx q[0];
rz(-1.5320677) q[0];
sx q[0];
rz(1.9963272) q[0];
rz(0.26248652) q[2];
sx q[2];
rz(-0.075857698) q[2];
sx q[2];
rz(-0.59466098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37564056) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(0.44771938) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9085957) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(-2.256957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-0.48015067) q[2];
sx q[2];
rz(-0.22228995) q[2];
rz(2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46626058) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(0.12761322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532967) q[0];
sx q[0];
rz(-1.7444381) q[0];
sx q[0];
rz(-2.2542473) q[0];
rz(0.80673809) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-0.15653175) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4096654) q[1];
sx q[1];
rz(-1.8095784) q[1];
sx q[1];
rz(3.0569397) q[1];
x q[2];
rz(2.8887799) q[3];
sx q[3];
rz(-2.1826934) q[3];
sx q[3];
rz(-1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(-2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(0.37240949) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83852864) q[0];
sx q[0];
rz(-1.4493363) q[0];
sx q[0];
rz(2.2785506) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98948688) q[2];
sx q[2];
rz(-1.8539691) q[2];
sx q[2];
rz(2.1894933) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1540893) q[1];
sx q[1];
rz(-1.9470864) q[1];
sx q[1];
rz(-2.980568) q[1];
x q[2];
rz(-2.1129205) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(-1.2937014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(0.27238971) q[2];
rz(2.2327936) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880496) q[0];
sx q[0];
rz(-1.4918259) q[0];
sx q[0];
rz(-2.700564) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73111515) q[2];
sx q[2];
rz(-1.7321246) q[2];
sx q[2];
rz(-0.71691712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0332898) q[1];
sx q[1];
rz(-1.8255594) q[1];
sx q[1];
rz(1.5441896) q[1];
rz(-pi) q[2];
rz(-0.39354126) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(-0.14906135) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60676735) q[0];
sx q[0];
rz(-1.2762428) q[0];
sx q[0];
rz(2.0539001) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51322333) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(2.9483587) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9698525) q[1];
sx q[1];
rz(-1.2029359) q[1];
sx q[1];
rz(3.1391341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9643781) q[3];
sx q[3];
rz(-2.2209077) q[3];
sx q[3];
rz(2.9110416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(0.28665001) q[2];
rz(-0.24946985) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(-2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0934802) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(-1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46465835) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(1.2901165) q[0];
x q[1];
rz(1.940735) q[2];
sx q[2];
rz(-1.2518034) q[2];
sx q[2];
rz(-1.5202886) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24504532) q[1];
sx q[1];
rz(-2.7198615) q[1];
sx q[1];
rz(1.1827724) q[1];
rz(-pi) q[2];
rz(-2.4909199) q[3];
sx q[3];
rz(-1.7400029) q[3];
sx q[3];
rz(0.31864877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-0.10350791) q[2];
rz(2.5497656) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(2.7517095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108111) q[0];
sx q[0];
rz(-0.99698742) q[0];
sx q[0];
rz(2.0033037) q[0];
rz(-pi) q[1];
rz(1.1578324) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(0.72380356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0695614) q[1];
sx q[1];
rz(-1.8998635) q[1];
sx q[1];
rz(1.1142932) q[1];
rz(-pi) q[2];
rz(-1.2478896) q[3];
sx q[3];
rz(-0.70526037) q[3];
sx q[3];
rz(-0.44912072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(2.4273382) q[2];
rz(0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-0.27639595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24647507) q[0];
sx q[0];
rz(-0.44024375) q[0];
sx q[0];
rz(-0.38245364) q[0];
rz(-pi) q[1];
rz(1.3771636) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(-0.60806882) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8918708) q[1];
sx q[1];
rz(-2.8686214) q[1];
sx q[1];
rz(2.194838) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8180088) q[3];
sx q[3];
rz(-0.97655481) q[3];
sx q[3];
rz(-1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(0.12750553) q[2];
rz(3.1048807) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(1.2344853) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(2.7888443) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.661783) q[0];
sx q[0];
rz(-2.3942238) q[0];
sx q[0];
rz(0.28967793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68600168) q[2];
sx q[2];
rz(-1.8324319) q[2];
sx q[2];
rz(0.54856448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0093617) q[1];
sx q[1];
rz(-1.8877601) q[1];
sx q[1];
rz(0.46225458) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.685931) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(2.4009005) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13070233) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(-2.4521811) q[3];
sx q[3];
rz(-2.0333615) q[3];
sx q[3];
rz(0.61483308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
