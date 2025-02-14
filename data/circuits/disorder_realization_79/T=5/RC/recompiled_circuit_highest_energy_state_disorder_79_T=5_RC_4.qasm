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
rz(1.5647178) q[0];
sx q[0];
rz(1.8160507) q[0];
sx q[0];
rz(9.9845822) q[0];
rz(-1.7436104) q[1];
sx q[1];
rz(-1.3275194) q[1];
sx q[1];
rz(1.7550069) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579642) q[0];
sx q[0];
rz(-1.0378519) q[0];
sx q[0];
rz(-0.29472875) q[0];
x q[1];
rz(1.6835911) q[2];
sx q[2];
rz(-1.8646282) q[2];
sx q[2];
rz(-0.71982671) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4475354) q[1];
sx q[1];
rz(-2.9172998) q[1];
sx q[1];
rz(2.7560225) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0069874115) q[3];
sx q[3];
rz(-1.8614359) q[3];
sx q[3];
rz(2.5537511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5379415) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(2.6998935) q[2];
rz(-0.80080992) q[3];
sx q[3];
rz(-1.2468612) q[3];
sx q[3];
rz(0.39113623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9073198) q[0];
sx q[0];
rz(-1.7664302) q[0];
sx q[0];
rz(-2.8148839) q[0];
rz(1.3455343) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(-0.84652841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3162982) q[0];
sx q[0];
rz(-0.40028737) q[0];
sx q[0];
rz(1.8497224) q[0];
rz(0.2791762) q[2];
sx q[2];
rz(-2.2344032) q[2];
sx q[2];
rz(0.070726591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7134756) q[1];
sx q[1];
rz(-1.5602116) q[1];
sx q[1];
rz(1.0416563) q[1];
rz(-pi) q[2];
rz(2.8277581) q[3];
sx q[3];
rz(-2.3646128) q[3];
sx q[3];
rz(1.8557567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.232406) q[2];
sx q[2];
rz(-1.0973944) q[2];
sx q[2];
rz(-0.48420134) q[2];
rz(-1.3062612) q[3];
sx q[3];
rz(-2.5808915) q[3];
sx q[3];
rz(-1.9561214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4803798) q[0];
sx q[0];
rz(-0.73037195) q[0];
sx q[0];
rz(0.54006201) q[0];
rz(-1.1506608) q[1];
sx q[1];
rz(-1.8692503) q[1];
sx q[1];
rz(-0.59404341) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4580989) q[0];
sx q[0];
rz(-2.9389295) q[0];
sx q[0];
rz(2.1432102) q[0];
x q[1];
rz(-2.4864628) q[2];
sx q[2];
rz(-1.3199886) q[2];
sx q[2];
rz(1.0162958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4773399) q[1];
sx q[1];
rz(-2.3788733) q[1];
sx q[1];
rz(-1.699317) q[1];
rz(-pi) q[2];
x q[2];
rz(1.698231) q[3];
sx q[3];
rz(-0.48526796) q[3];
sx q[3];
rz(-1.8125774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.322304) q[2];
sx q[2];
rz(-2.1257336) q[2];
sx q[2];
rz(-2.4508396) q[2];
rz(2.6639719) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(0.42948183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884884) q[0];
sx q[0];
rz(-2.7957343) q[0];
sx q[0];
rz(0.75230569) q[0];
rz(-0.90657702) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(0.22901542) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7389694) q[0];
sx q[0];
rz(-1.7949901) q[0];
sx q[0];
rz(-2.8408308) q[0];
rz(0.36831736) q[2];
sx q[2];
rz(-0.57007705) q[2];
sx q[2];
rz(-1.7961479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17002925) q[1];
sx q[1];
rz(-1.3336412) q[1];
sx q[1];
rz(-0.27598652) q[1];
rz(1.205577) q[3];
sx q[3];
rz(-2.1250279) q[3];
sx q[3];
rz(0.54868719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-3.0276073) q[2];
sx q[2];
rz(-2.1086741) q[2];
rz(-2.0756663) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(1.3037995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579987) q[0];
sx q[0];
rz(-0.086253919) q[0];
sx q[0];
rz(1.3294543) q[0];
rz(-2.6138002) q[1];
sx q[1];
rz(-1.3161491) q[1];
sx q[1];
rz(-2.7425308) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.746839) q[0];
sx q[0];
rz(-1.2195713) q[0];
sx q[0];
rz(0.082562692) q[0];
rz(0.73164263) q[2];
sx q[2];
rz(-1.2564617) q[2];
sx q[2];
rz(-0.86282496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1414772) q[1];
sx q[1];
rz(-2.2908604) q[1];
sx q[1];
rz(-1.9898371) q[1];
rz(-pi) q[2];
rz(-0.87598975) q[3];
sx q[3];
rz(-1.717698) q[3];
sx q[3];
rz(1.5796652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.993678) q[2];
sx q[2];
rz(-0.75239158) q[2];
sx q[2];
rz(1.0047151) q[2];
rz(0.29804722) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(2.0731488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816958) q[0];
sx q[0];
rz(-1.5157461) q[0];
sx q[0];
rz(1.1915278) q[0];
rz(0.53939348) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(0.76621145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3193885) q[0];
sx q[0];
rz(-2.4766716) q[0];
sx q[0];
rz(0.98916905) q[0];
x q[1];
rz(-0.79121235) q[2];
sx q[2];
rz(-2.3177766) q[2];
sx q[2];
rz(0.3283161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4529383) q[1];
sx q[1];
rz(-2.3804745) q[1];
sx q[1];
rz(-1.903694) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9412463) q[3];
sx q[3];
rz(-1.4806357) q[3];
sx q[3];
rz(2.3925635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3300276) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(-0.88328973) q[2];
rz(0.53089321) q[3];
sx q[3];
rz(-2.5575432) q[3];
sx q[3];
rz(-2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4009878) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(0.64939943) q[0];
rz(-2.7679288) q[1];
sx q[1];
rz(-2.3039736) q[1];
sx q[1];
rz(1.9992453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9050347) q[0];
sx q[0];
rz(-1.243374) q[0];
sx q[0];
rz(-2.9740646) q[0];
rz(-pi) q[1];
rz(-3.0685038) q[2];
sx q[2];
rz(-2.9311649) q[2];
sx q[2];
rz(3.0874012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2616407) q[1];
sx q[1];
rz(-1.2625719) q[1];
sx q[1];
rz(2.37642) q[1];
x q[2];
rz(2.029923) q[3];
sx q[3];
rz(-1.6111501) q[3];
sx q[3];
rz(-0.99276517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.00804) q[2];
sx q[2];
rz(-0.69591659) q[2];
sx q[2];
rz(2.7827061) q[2];
rz(-2.533203) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(2.2801094) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95978874) q[0];
sx q[0];
rz(-2.602808) q[0];
sx q[0];
rz(0.96610075) q[0];
rz(2.6995662) q[1];
sx q[1];
rz(-1.8531046) q[1];
sx q[1];
rz(2.7332222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20710564) q[0];
sx q[0];
rz(-1.4382595) q[0];
sx q[0];
rz(0.026959628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5277065) q[2];
sx q[2];
rz(-1.7183964) q[2];
sx q[2];
rz(1.4606295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8345771) q[1];
sx q[1];
rz(-2.53741) q[1];
sx q[1];
rz(1.8240364) q[1];
x q[2];
rz(0.73159742) q[3];
sx q[3];
rz(-1.3773736) q[3];
sx q[3];
rz(0.0091656589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8926706) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(-1.9425617) q[2];
rz(1.2976973) q[3];
sx q[3];
rz(-2.0632931) q[3];
sx q[3];
rz(0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2099828) q[0];
sx q[0];
rz(-1.7117806) q[0];
sx q[0];
rz(2.0344875) q[0];
rz(2.9529086) q[1];
sx q[1];
rz(-2.76177) q[1];
sx q[1];
rz(-1.3245378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134091) q[0];
sx q[0];
rz(-1.8556917) q[0];
sx q[0];
rz(0.5351184) q[0];
rz(-pi) q[1];
rz(0.021090551) q[2];
sx q[2];
rz(-0.97375768) q[2];
sx q[2];
rz(1.8151547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75211084) q[1];
sx q[1];
rz(-0.5403203) q[1];
sx q[1];
rz(-1.575927) q[1];
x q[2];
rz(-2.5741413) q[3];
sx q[3];
rz(-0.78560116) q[3];
sx q[3];
rz(-2.7170865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.79335) q[2];
sx q[2];
rz(-2.5620354) q[2];
sx q[2];
rz(-2.9938475) q[2];
rz(-0.82792264) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350236) q[0];
sx q[0];
rz(-0.35025418) q[0];
sx q[0];
rz(1.6178004) q[0];
rz(-1.2500457) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(-0.42125431) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68116248) q[0];
sx q[0];
rz(-1.958485) q[0];
sx q[0];
rz(-1.1739385) q[0];
rz(-pi) q[1];
rz(0.12261196) q[2];
sx q[2];
rz(-0.62527925) q[2];
sx q[2];
rz(1.774315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7819216) q[1];
sx q[1];
rz(-0.087347833) q[1];
sx q[1];
rz(-3.0328301) q[1];
rz(-pi) q[2];
rz(-2.6700745) q[3];
sx q[3];
rz(-2.6425411) q[3];
sx q[3];
rz(1.7744482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0201575) q[2];
sx q[2];
rz(-0.73548663) q[2];
sx q[2];
rz(-0.45041931) q[2];
rz(-1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(-2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12238518) q[0];
sx q[0];
rz(-1.6078124) q[0];
sx q[0];
rz(-1.5851198) q[0];
rz(0.71313329) q[1];
sx q[1];
rz(-1.5596371) q[1];
sx q[1];
rz(0.023991931) q[1];
rz(2.6789011) q[2];
sx q[2];
rz(-1.0924871) q[2];
sx q[2];
rz(0.15747216) q[2];
rz(-1.6433257) q[3];
sx q[3];
rz(-1.4909496) q[3];
sx q[3];
rz(-0.95076233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
