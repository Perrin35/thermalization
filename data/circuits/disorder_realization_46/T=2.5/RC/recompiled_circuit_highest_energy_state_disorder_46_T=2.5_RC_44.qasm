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
rz(0.47705874) q[0];
sx q[0];
rz(2.8395489) q[0];
sx q[0];
rz(10.680847) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(-1.3060275) q[1];
sx q[1];
rz(-1.1285055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94621979) q[0];
sx q[0];
rz(-1.2913307) q[0];
sx q[0];
rz(-2.165876) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41646012) q[2];
sx q[2];
rz(-2.4774733) q[2];
sx q[2];
rz(1.0937607) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30762954) q[1];
sx q[1];
rz(-2.0700196) q[1];
sx q[1];
rz(2.7021033) q[1];
rz(-pi) q[2];
rz(-1.1323186) q[3];
sx q[3];
rz(-2.1578159) q[3];
sx q[3];
rz(1.7917716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6065373) q[2];
sx q[2];
rz(-0.29416072) q[2];
sx q[2];
rz(-1.9425707) q[2];
rz(-0.27896518) q[3];
sx q[3];
rz(-1.4594892) q[3];
sx q[3];
rz(0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8798384) q[0];
sx q[0];
rz(-0.45670515) q[0];
sx q[0];
rz(-2.7754011) q[0];
rz(1.2884033) q[1];
sx q[1];
rz(-2.2346456) q[1];
sx q[1];
rz(-2.8790976) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2912882) q[0];
sx q[0];
rz(-1.2748963) q[0];
sx q[0];
rz(1.759191) q[0];
x q[1];
rz(-1.8196202) q[2];
sx q[2];
rz(-1.6273398) q[2];
sx q[2];
rz(-0.6565993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9001635) q[1];
sx q[1];
rz(-1.8594605) q[1];
sx q[1];
rz(2.7735387) q[1];
x q[2];
rz(3.1163773) q[3];
sx q[3];
rz(-1.8256911) q[3];
sx q[3];
rz(-2.1372634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3890106) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(-0.1788204) q[2];
rz(-0.016077476) q[3];
sx q[3];
rz(-0.5355081) q[3];
sx q[3];
rz(2.2107562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7054684) q[0];
sx q[0];
rz(-0.12691623) q[0];
sx q[0];
rz(-0.57417589) q[0];
rz(-0.15159675) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(-1.2281598) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77297091) q[0];
sx q[0];
rz(-0.54062343) q[0];
sx q[0];
rz(-1.2631046) q[0];
x q[1];
rz(2.1636859) q[2];
sx q[2];
rz(-1.6276858) q[2];
sx q[2];
rz(-1.4658734) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1054525) q[1];
sx q[1];
rz(-0.6195809) q[1];
sx q[1];
rz(-1.416227) q[1];
x q[2];
rz(-2.853573) q[3];
sx q[3];
rz(-1.1041118) q[3];
sx q[3];
rz(-0.53190255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0000618) q[2];
sx q[2];
rz(-0.63566339) q[2];
sx q[2];
rz(0.89111152) q[2];
rz(-0.38854232) q[3];
sx q[3];
rz(-1.5107692) q[3];
sx q[3];
rz(-0.33963206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5163088) q[0];
sx q[0];
rz(-0.044965222) q[0];
sx q[0];
rz(0.48211023) q[0];
rz(-0.66857839) q[1];
sx q[1];
rz(-1.5336288) q[1];
sx q[1];
rz(0.63749981) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26698819) q[0];
sx q[0];
rz(-0.98867304) q[0];
sx q[0];
rz(0.28018392) q[0];
rz(1.7155816) q[2];
sx q[2];
rz(-2.4630483) q[2];
sx q[2];
rz(-2.2172287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60739775) q[1];
sx q[1];
rz(-3.0574905) q[1];
sx q[1];
rz(-2.7264595) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96781625) q[3];
sx q[3];
rz(-0.42438904) q[3];
sx q[3];
rz(-2.5924204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.94347) q[2];
sx q[2];
rz(-1.024647) q[2];
sx q[2];
rz(-0.060297273) q[2];
rz(-0.30334011) q[3];
sx q[3];
rz(-0.077849418) q[3];
sx q[3];
rz(2.8764603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1497134) q[0];
sx q[0];
rz(-0.35357058) q[0];
sx q[0];
rz(0.41325945) q[0];
rz(-2.7464271) q[1];
sx q[1];
rz(-1.3576077) q[1];
sx q[1];
rz(-1.0023592) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226333) q[0];
sx q[0];
rz(-1.1341242) q[0];
sx q[0];
rz(0.60499259) q[0];
rz(-pi) q[1];
rz(1.7806899) q[2];
sx q[2];
rz(-1.5355459) q[2];
sx q[2];
rz(1.1473224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9262979) q[1];
sx q[1];
rz(-1.7855254) q[1];
sx q[1];
rz(-2.2533827) q[1];
rz(-1.2484364) q[3];
sx q[3];
rz(-1.768422) q[3];
sx q[3];
rz(1.884643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5137382) q[2];
sx q[2];
rz(-1.0971053) q[2];
sx q[2];
rz(1.6045137) q[2];
rz(-1.9419935) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(2.3698923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5143249) q[0];
sx q[0];
rz(-1.6384614) q[0];
sx q[0];
rz(-0.36870536) q[0];
rz(0.74327028) q[1];
sx q[1];
rz(-0.5629881) q[1];
sx q[1];
rz(-2.7875913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6021332) q[0];
sx q[0];
rz(-1.51855) q[0];
sx q[0];
rz(2.2118072) q[0];
x q[1];
rz(-1.2576302) q[2];
sx q[2];
rz(-1.0651565) q[2];
sx q[2];
rz(-0.87228105) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88601513) q[1];
sx q[1];
rz(-1.3929345) q[1];
sx q[1];
rz(-2.1484003) q[1];
x q[2];
rz(-1.9571113) q[3];
sx q[3];
rz(-2.1710178) q[3];
sx q[3];
rz(1.6627897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.02504286) q[2];
sx q[2];
rz(-1.9321238) q[2];
sx q[2];
rz(3.1361191) q[2];
rz(2.6239851) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(3.1143809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7356877) q[0];
sx q[0];
rz(-0.55957496) q[0];
sx q[0];
rz(2.9285808) q[0];
rz(-1.6464015) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(2.2921553) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54100688) q[0];
sx q[0];
rz(-1.6141506) q[0];
sx q[0];
rz(0.61886532) q[0];
rz(-pi) q[1];
x q[1];
rz(0.030730129) q[2];
sx q[2];
rz(-1.5275914) q[2];
sx q[2];
rz(1.876056) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7905459) q[1];
sx q[1];
rz(-2.4645269) q[1];
sx q[1];
rz(0.76404567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3091503) q[3];
sx q[3];
rz(-2.5980686) q[3];
sx q[3];
rz(1.8147008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9269632) q[2];
sx q[2];
rz(-1.233485) q[2];
sx q[2];
rz(0.54369533) q[2];
rz(-2.8804273) q[3];
sx q[3];
rz(-1.5615014) q[3];
sx q[3];
rz(2.9847667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9549304) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(0.18761158) q[0];
rz(-1.1232417) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(0.16256464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190348) q[0];
sx q[0];
rz(-1.8936689) q[0];
sx q[0];
rz(2.4486008) q[0];
rz(-0.17485042) q[2];
sx q[2];
rz(-1.0533337) q[2];
sx q[2];
rz(-0.78570056) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2073719) q[1];
sx q[1];
rz(-1.6160864) q[1];
sx q[1];
rz(-0.95691935) q[1];
rz(-pi) q[2];
rz(-0.97450337) q[3];
sx q[3];
rz(-1.664229) q[3];
sx q[3];
rz(0.0020041618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2742013) q[2];
sx q[2];
rz(-0.31495366) q[2];
sx q[2];
rz(-2.9401722) q[2];
rz(-2.6650186) q[3];
sx q[3];
rz(-1.3786992) q[3];
sx q[3];
rz(0.99982888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808712) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(2.5647105) q[0];
rz(0.30413973) q[1];
sx q[1];
rz(-2.0174593) q[1];
sx q[1];
rz(2.0374128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9228684) q[0];
sx q[0];
rz(-1.4822685) q[0];
sx q[0];
rz(1.1909816) q[0];
rz(-2.5726692) q[2];
sx q[2];
rz(-1.6732135) q[2];
sx q[2];
rz(-1.2889287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39958274) q[1];
sx q[1];
rz(-2.4192657) q[1];
sx q[1];
rz(-0.96256017) q[1];
rz(1.6100092) q[3];
sx q[3];
rz(-0.3914507) q[3];
sx q[3];
rz(1.1757562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96678174) q[2];
sx q[2];
rz(-0.2468144) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(-1.5481868) q[3];
sx q[3];
rz(-0.67845172) q[3];
sx q[3];
rz(2.9233176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9504647) q[0];
sx q[0];
rz(-0.8466962) q[0];
sx q[0];
rz(-2.6306613) q[0];
rz(-1.9966985) q[1];
sx q[1];
rz(-1.9547918) q[1];
sx q[1];
rz(1.6330382) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.78054) q[0];
sx q[0];
rz(-2.1117438) q[0];
sx q[0];
rz(-1.279328) q[0];
rz(2.050348) q[2];
sx q[2];
rz(-2.6065718) q[2];
sx q[2];
rz(-2.1370691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33433769) q[1];
sx q[1];
rz(-0.48408135) q[1];
sx q[1];
rz(2.7527807) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8415751) q[3];
sx q[3];
rz(-1.8138061) q[3];
sx q[3];
rz(-2.9307765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7213584) q[2];
sx q[2];
rz(-0.63174641) q[2];
sx q[2];
rz(1.7101804) q[2];
rz(-3.1259649) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90841993) q[0];
sx q[0];
rz(-1.7411727) q[0];
sx q[0];
rz(2.2199051) q[0];
rz(-1.635101) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(-1.9533515) q[2];
sx q[2];
rz(-1.3024223) q[2];
sx q[2];
rz(-0.060729751) q[2];
rz(1.9962068) q[3];
sx q[3];
rz(-3.0563995) q[3];
sx q[3];
rz(1.2264768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
