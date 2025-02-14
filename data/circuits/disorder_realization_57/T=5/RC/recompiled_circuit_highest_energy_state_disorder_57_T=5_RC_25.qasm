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
rz(1.8013826) q[0];
sx q[0];
rz(-0.27685452) q[0];
sx q[0];
rz(-1.0387596) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(-2.3050397) q[1];
sx q[1];
rz(0.41681448) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5688516) q[0];
sx q[0];
rz(-0.73264719) q[0];
sx q[0];
rz(2.9824663) q[0];
x q[1];
rz(2.7567467) q[2];
sx q[2];
rz(-2.550808) q[2];
sx q[2];
rz(2.7123775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8326679) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(-0.69242386) q[1];
x q[2];
rz(-2.9800426) q[3];
sx q[3];
rz(-1.6027502) q[3];
sx q[3];
rz(-2.3475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(1.2539585) q[2];
rz(1.5818671) q[3];
sx q[3];
rz(-2.4032148) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9480243) q[0];
sx q[0];
rz(-1.3688315) q[0];
sx q[0];
rz(0.76876202) q[0];
rz(-1.3668775) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(2.1034525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5032673) q[0];
sx q[0];
rz(-1.5794139) q[0];
sx q[0];
rz(1.5609972) q[0];
rz(1.9351472) q[2];
sx q[2];
rz(-0.88308217) q[2];
sx q[2];
rz(-1.5503413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8703192) q[1];
sx q[1];
rz(-0.73713494) q[1];
sx q[1];
rz(2.5802274) q[1];
rz(-pi) q[2];
rz(-1.5350545) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(-1.8349748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.495503) q[2];
sx q[2];
rz(-1.2868519) q[2];
sx q[2];
rz(2.3213279) q[2];
rz(0.25343728) q[3];
sx q[3];
rz(-0.35496747) q[3];
sx q[3];
rz(2.7804815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8609817) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(0.88731998) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-2.2128426) q[1];
sx q[1];
rz(2.0022154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47592012) q[0];
sx q[0];
rz(-2.7662656) q[0];
sx q[0];
rz(-0.39552839) q[0];
rz(1.0372889) q[2];
sx q[2];
rz(-2.5441558) q[2];
sx q[2];
rz(3.1410599) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2816086) q[1];
sx q[1];
rz(-1.3251332) q[1];
sx q[1];
rz(-1.5558234) q[1];
rz(-0.091413012) q[3];
sx q[3];
rz(-2.6210636) q[3];
sx q[3];
rz(-0.5239858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1246216) q[2];
sx q[2];
rz(-1.4554224) q[2];
sx q[2];
rz(-2.7023081) q[2];
rz(-2.6708421) q[3];
sx q[3];
rz(-0.079340383) q[3];
sx q[3];
rz(1.144217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(3.0261107) q[0];
rz(-0.11784095) q[1];
sx q[1];
rz(-1.203048) q[1];
sx q[1];
rz(1.8353362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7457434) q[0];
sx q[0];
rz(-1.1538528) q[0];
sx q[0];
rz(-2.8265619) q[0];
rz(1.9491862) q[2];
sx q[2];
rz(-0.94589627) q[2];
sx q[2];
rz(1.6922127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0178284) q[1];
sx q[1];
rz(-0.99010795) q[1];
sx q[1];
rz(-2.3821217) q[1];
x q[2];
rz(0.0751817) q[3];
sx q[3];
rz(-2.1422221) q[3];
sx q[3];
rz(-2.927185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.304504) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(-0.15466776) q[2];
rz(-1.6020487) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(2.535517) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57347572) q[0];
sx q[0];
rz(-1.0433759) q[0];
sx q[0];
rz(0.83531761) q[0];
rz(2.4256445) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(-0.11944019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12163945) q[0];
sx q[0];
rz(-3.0832735) q[0];
sx q[0];
rz(-1.8644237) q[0];
x q[1];
rz(1.7238462) q[2];
sx q[2];
rz(-2.0156392) q[2];
sx q[2];
rz(2.9314624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24158289) q[1];
sx q[1];
rz(-1.1092343) q[1];
sx q[1];
rz(-0.66630967) q[1];
x q[2];
rz(-2.4923513) q[3];
sx q[3];
rz(-0.63707817) q[3];
sx q[3];
rz(1.0604881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9014088) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(-3.0933464) q[3];
sx q[3];
rz(-1.9082853) q[3];
sx q[3];
rz(-0.72859305) q[3];
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
rz(-pi/2) q[0];
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
rz(2.401684) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(-1.1309062) q[0];
rz(-0.4785969) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(-1.4071677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3597108) q[0];
sx q[0];
rz(-3.0880083) q[0];
sx q[0];
rz(-3.1126541) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2931104) q[2];
sx q[2];
rz(-1.7847848) q[2];
sx q[2];
rz(3.1297562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9764912) q[1];
sx q[1];
rz(-0.25549421) q[1];
sx q[1];
rz(-2.0352023) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21295548) q[3];
sx q[3];
rz(-1.7476488) q[3];
sx q[3];
rz(2.9573754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56167928) q[2];
sx q[2];
rz(-2.734197) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(-2.0922349) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(-1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2105763) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(1.806102) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(0.30805045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086723344) q[0];
sx q[0];
rz(-0.44741524) q[0];
sx q[0];
rz(-0.50680508) q[0];
rz(-pi) q[1];
rz(2.5045583) q[2];
sx q[2];
rz(-2.5278628) q[2];
sx q[2];
rz(1.6047275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8561473) q[1];
sx q[1];
rz(-0.79570192) q[1];
sx q[1];
rz(2.281762) q[1];
rz(1.4542463) q[3];
sx q[3];
rz(-2.7532059) q[3];
sx q[3];
rz(-1.5248869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7293952) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(3.0103053) q[2];
rz(0.73733759) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(-2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103631) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(-0.53037733) q[0];
rz(-1.7165548) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(-0.72149611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58396858) q[0];
sx q[0];
rz(-0.97088214) q[0];
sx q[0];
rz(2.1198089) q[0];
rz(2.6725298) q[2];
sx q[2];
rz(-1.4387759) q[2];
sx q[2];
rz(-1.309883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87807314) q[1];
sx q[1];
rz(-1.7612447) q[1];
sx q[1];
rz(1.7909122) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.4777769) q[2];
rz(-1.9780698) q[3];
sx q[3];
rz(-1.082837) q[3];
sx q[3];
rz(1.6400953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32996938) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(-0.60923088) q[0];
rz(1.5513264) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.685166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44269366) q[0];
sx q[0];
rz(-0.44166587) q[0];
sx q[0];
rz(2.1864258) q[0];
x q[1];
rz(-3.0266043) q[2];
sx q[2];
rz(-2.6536332) q[2];
sx q[2];
rz(1.3249299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.506914) q[1];
sx q[1];
rz(-1.4089481) q[1];
sx q[1];
rz(1.1062201) q[1];
rz(-2.5129065) q[3];
sx q[3];
rz(-2.1910724) q[3];
sx q[3];
rz(-1.6339906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41813254) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(0.20600016) q[3];
sx q[3];
rz(-2.7166631) q[3];
sx q[3];
rz(0.4900842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708165) q[0];
sx q[0];
rz(-2.4801065) q[0];
sx q[0];
rz(-0.01195512) q[0];
rz(-1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(2.5451122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79967989) q[0];
sx q[0];
rz(-1.1014859) q[0];
sx q[0];
rz(-1.5494359) q[0];
rz(-pi) q[1];
rz(0.36294655) q[2];
sx q[2];
rz(-1.3273113) q[2];
sx q[2];
rz(-1.4371265) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4778127) q[1];
sx q[1];
rz(-0.48250178) q[1];
sx q[1];
rz(-1.4003808) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6634253) q[3];
sx q[3];
rz(-2.8088125) q[3];
sx q[3];
rz(-1.5733583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(-2.9445924) q[2];
rz(1.152285) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(2.6939189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780509) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(0.20881431) q[1];
sx q[1];
rz(-1.6046235) q[1];
sx q[1];
rz(-1.0135289) q[1];
rz(1.456719) q[2];
sx q[2];
rz(-2.4488505) q[2];
sx q[2];
rz(1.746576) q[2];
rz(0.98255929) q[3];
sx q[3];
rz(-1.1362024) q[3];
sx q[3];
rz(0.35432651) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
