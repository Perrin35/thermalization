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
rz(-0.62491971) q[0];
sx q[0];
rz(-1.3345557) q[0];
sx q[0];
rz(-0.053243756) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(6.1558131) q[1];
sx q[1];
rz(11.131412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.931728) q[0];
sx q[0];
rz(-1.4357982) q[0];
sx q[0];
rz(0.28688669) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8743319) q[2];
sx q[2];
rz(-1.4556985) q[2];
sx q[2];
rz(-2.5989344) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8813144) q[1];
sx q[1];
rz(-1.5173755) q[1];
sx q[1];
rz(-0.39398663) q[1];
x q[2];
rz(2.1987183) q[3];
sx q[3];
rz(-0.74340313) q[3];
sx q[3];
rz(-1.8441895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3367553) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(-2.8678144) q[2];
rz(-1.9233507) q[3];
sx q[3];
rz(-2.4680586) q[3];
sx q[3];
rz(0.28606733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022920595) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(-2.3619695) q[0];
rz(-0.66462213) q[1];
sx q[1];
rz(-1.2088935) q[1];
sx q[1];
rz(1.2453311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561403) q[0];
sx q[0];
rz(-0.13053556) q[0];
sx q[0];
rz(2.0640316) q[0];
rz(-pi) q[1];
x q[1];
rz(0.092463569) q[2];
sx q[2];
rz(-0.86148724) q[2];
sx q[2];
rz(0.48395983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1511953) q[1];
sx q[1];
rz(-1.6315236) q[1];
sx q[1];
rz(1.7833738) q[1];
rz(-pi) q[2];
rz(0.20540463) q[3];
sx q[3];
rz(-1.0213791) q[3];
sx q[3];
rz(-2.7618264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43210426) q[2];
sx q[2];
rz(-2.4435142) q[2];
sx q[2];
rz(-0.049987642) q[2];
rz(-1.6677808) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(-2.2059691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525928) q[0];
sx q[0];
rz(-2.279156) q[0];
sx q[0];
rz(-1.1784026) q[0];
rz(1.1475457) q[1];
sx q[1];
rz(-2.9541364) q[1];
sx q[1];
rz(2.5386834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5459203) q[0];
sx q[0];
rz(-1.4688558) q[0];
sx q[0];
rz(1.4096745) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88073894) q[2];
sx q[2];
rz(-1.1549486) q[2];
sx q[2];
rz(1.8772454) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5192831) q[1];
sx q[1];
rz(-1.9819248) q[1];
sx q[1];
rz(-2.4101188) q[1];
rz(2.1985377) q[3];
sx q[3];
rz(-0.78305093) q[3];
sx q[3];
rz(2.8179226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9511562) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(-0.3832761) q[2];
rz(-1.8303309) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(-0.12152984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7168147) q[0];
sx q[0];
rz(-2.1489693) q[0];
sx q[0];
rz(2.2130261) q[0];
rz(-0.83944744) q[1];
sx q[1];
rz(-1.8239559) q[1];
sx q[1];
rz(-0.80926698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0220714) q[0];
sx q[0];
rz(-2.4520527) q[0];
sx q[0];
rz(-0.9718231) q[0];
x q[1];
rz(1.7401314) q[2];
sx q[2];
rz(-1.9953097) q[2];
sx q[2];
rz(2.9237539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70185858) q[1];
sx q[1];
rz(-1.3254406) q[1];
sx q[1];
rz(-0.071003242) q[1];
rz(2.2146642) q[3];
sx q[3];
rz(-1.9044821) q[3];
sx q[3];
rz(-2.9179171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75260085) q[2];
sx q[2];
rz(-1.6099124) q[2];
sx q[2];
rz(-1.4468225) q[2];
rz(1.1270479) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(0.62895044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.22555722) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(-1.7115364) q[0];
rz(2.9043708) q[1];
sx q[1];
rz(-0.96313852) q[1];
sx q[1];
rz(2.846834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65564102) q[0];
sx q[0];
rz(-1.7136586) q[0];
sx q[0];
rz(0.44012196) q[0];
rz(-pi) q[1];
rz(-0.37173231) q[2];
sx q[2];
rz(-2.1991538) q[2];
sx q[2];
rz(-0.5231176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.393906) q[1];
sx q[1];
rz(-1.9853885) q[1];
sx q[1];
rz(1.0769597) q[1];
rz(-pi) q[2];
rz(-1.421991) q[3];
sx q[3];
rz(-1.0019046) q[3];
sx q[3];
rz(0.97967813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99990591) q[2];
sx q[2];
rz(-2.9408231) q[2];
sx q[2];
rz(3.0086009) q[2];
rz(0.76987949) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(-0.31908116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517476) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(1.4122562) q[0];
rz(3.0899561) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(0.19270611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.76336) q[0];
sx q[0];
rz(-0.51412778) q[0];
sx q[0];
rz(-3.1030416) q[0];
rz(-pi) q[1];
rz(-2.4102845) q[2];
sx q[2];
rz(-2.0530862) q[2];
sx q[2];
rz(2.363229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7900675) q[1];
sx q[1];
rz(-2.0242825) q[1];
sx q[1];
rz(-2.044911) q[1];
x q[2];
rz(0.64151836) q[3];
sx q[3];
rz(-1.2605485) q[3];
sx q[3];
rz(-1.1603242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.073033832) q[2];
sx q[2];
rz(-1.2898338) q[2];
sx q[2];
rz(1.3812836) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-0.88740715) q[3];
sx q[3];
rz(-1.7611586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16159049) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(-2.7951151) q[0];
rz(-1.0193635) q[1];
sx q[1];
rz(-0.66277021) q[1];
sx q[1];
rz(-1.4541218) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0270868) q[0];
sx q[0];
rz(-0.99106228) q[0];
sx q[0];
rz(2.113335) q[0];
x q[1];
rz(2.2242111) q[2];
sx q[2];
rz(-2.0163943) q[2];
sx q[2];
rz(0.080597045) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84216181) q[1];
sx q[1];
rz(-0.60363942) q[1];
sx q[1];
rz(-2.0718715) q[1];
rz(-2.1460974) q[3];
sx q[3];
rz(-0.18333215) q[3];
sx q[3];
rz(-0.14106476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.60160294) q[2];
sx q[2];
rz(-1.8113965) q[2];
sx q[2];
rz(-0.23923624) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-0.68238634) q[3];
sx q[3];
rz(0.95170963) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20188986) q[0];
sx q[0];
rz(-1.4677784) q[0];
sx q[0];
rz(1.7643167) q[0];
rz(-1.9210303) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(2.9753704) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9181053) q[0];
sx q[0];
rz(-1.1150196) q[0];
sx q[0];
rz(2.8546643) q[0];
rz(-pi) q[1];
rz(-0.34923415) q[2];
sx q[2];
rz(-2.376261) q[2];
sx q[2];
rz(-2.1527596) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4835755) q[1];
sx q[1];
rz(-1.3596351) q[1];
sx q[1];
rz(-1.6021452) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0395501) q[3];
sx q[3];
rz(-0.47787468) q[3];
sx q[3];
rz(0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.347747) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(-1.0469077) q[2];
rz(1.2547803) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5407402) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(-2.0284213) q[0];
rz(-0.81740776) q[1];
sx q[1];
rz(-1.6956885) q[1];
sx q[1];
rz(-0.36849749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2042646) q[0];
sx q[0];
rz(-2.6952792) q[0];
sx q[0];
rz(-2.242779) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0037874) q[2];
sx q[2];
rz(-0.74713444) q[2];
sx q[2];
rz(-1.2580308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1296325) q[1];
sx q[1];
rz(-2.3756177) q[1];
sx q[1];
rz(-0.38928826) q[1];
rz(2.542607) q[3];
sx q[3];
rz(-2.3040207) q[3];
sx q[3];
rz(-0.70928228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14785279) q[2];
sx q[2];
rz(-2.5233614) q[2];
sx q[2];
rz(1.3573793) q[2];
rz(2.3944858) q[3];
sx q[3];
rz(-2.1785469) q[3];
sx q[3];
rz(0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52637446) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(2.0988462) q[0];
rz(-2.3710947) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(2.3811293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2337991) q[0];
sx q[0];
rz(-1.9007508) q[0];
sx q[0];
rz(-0.25262649) q[0];
x q[1];
rz(1.7793525) q[2];
sx q[2];
rz(-2.957666) q[2];
sx q[2];
rz(2.6435801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.684243) q[1];
sx q[1];
rz(-1.4919229) q[1];
sx q[1];
rz(-2.7175987) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6728739) q[3];
sx q[3];
rz(-0.7112452) q[3];
sx q[3];
rz(-2.6576692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20327917) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(-2.7321613) q[2];
rz(-1.5583386) q[3];
sx q[3];
rz(-2.6810472) q[3];
sx q[3];
rz(2.515365) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4938477) q[0];
sx q[0];
rz(-1.0418325) q[0];
sx q[0];
rz(-2.7000725) q[0];
rz(-1.582675) q[1];
sx q[1];
rz(-1.9428923) q[1];
sx q[1];
rz(2.518242) q[1];
rz(1.5123488) q[2];
sx q[2];
rz(-1.4679906) q[2];
sx q[2];
rz(1.0124258) q[2];
rz(1.7033475) q[3];
sx q[3];
rz(-1.3302531) q[3];
sx q[3];
rz(2.1147685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
