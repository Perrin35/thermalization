OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(2.6846057) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39014434) q[0];
sx q[0];
rz(-1.4353308) q[0];
sx q[0];
rz(-0.22060237) q[0];
rz(0.692042) q[2];
sx q[2];
rz(-0.11728742) q[2];
sx q[2];
rz(-0.28809822) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0253804) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(1.8616574) q[1];
rz(-pi) q[2];
rz(1.0702707) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(-0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(-0.74222773) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(0.80274686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6532324) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(1.0351719) q[0];
rz(-pi) q[1];
rz(-1.0674092) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(-2.3561321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74233494) q[1];
sx q[1];
rz(-2.5207673) q[1];
sx q[1];
rz(-1.6595112) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0117399) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(-0.065404281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(-0.60570335) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-0.28513518) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(-1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0403255) q[0];
sx q[0];
rz(-0.87653941) q[0];
sx q[0];
rz(2.5800152) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8324098) q[2];
sx q[2];
rz(-0.90484607) q[2];
sx q[2];
rz(-3.0004629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(1.3819441) q[1];
rz(-pi) q[2];
rz(-0.47311802) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(0.28120041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.7335256) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4144856) q[0];
sx q[0];
rz(-1.3654728) q[0];
sx q[0];
rz(-1.7617102) q[0];
rz(0.71422691) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(-0.54745882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9208593) q[1];
sx q[1];
rz(-1.024106) q[1];
sx q[1];
rz(-0.85733719) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4227082) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-2.4812223) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.1594835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449074) q[0];
sx q[0];
rz(-2.3527745) q[0];
sx q[0];
rz(-0.44992723) q[0];
x q[1];
rz(-0.15578606) q[2];
sx q[2];
rz(-2.163379) q[2];
sx q[2];
rz(0.32183811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3573208) q[1];
sx q[1];
rz(-1.9958152) q[1];
sx q[1];
rz(2.1678863) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0347576) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(-2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5746675) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(3.040722) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8773008) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(0.63440462) q[0];
x q[1];
rz(2.6249044) q[2];
sx q[2];
rz(-0.54737216) q[2];
sx q[2];
rz(-1.9043497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74360352) q[1];
sx q[1];
rz(-1.7730224) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-pi) q[2];
rz(-0.50562596) q[3];
sx q[3];
rz(-2.6856832) q[3];
sx q[3];
rz(1.0491919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6763517) q[0];
sx q[0];
rz(-0.40333336) q[0];
sx q[0];
rz(-0.98531918) q[0];
x q[1];
rz(-0.22600941) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(1.7058536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2076599) q[1];
sx q[1];
rz(-1.7001517) q[1];
sx q[1];
rz(1.5473066) q[1];
x q[2];
rz(1.5457821) q[3];
sx q[3];
rz(-0.56846148) q[3];
sx q[3];
rz(0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-2.2040099) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(3.0623073) q[0];
rz(1.1212564) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.942873) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872902) q[0];
sx q[0];
rz(-2.0044921) q[0];
sx q[0];
rz(1.4951697) q[0];
rz(1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(1.5944634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.291154) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(0.50317851) q[1];
rz(0.46320398) q[3];
sx q[3];
rz(-1.3363046) q[3];
sx q[3];
rz(-0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(2.5532706) q[0];
rz(-0.026780216) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(-2.2081597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3073472) q[0];
sx q[0];
rz(-1.5114307) q[0];
sx q[0];
rz(0.10660118) q[0];
rz(-pi) q[1];
rz(0.8546631) q[2];
sx q[2];
rz(-2.6551506) q[2];
sx q[2];
rz(-1.9180627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48298207) q[1];
sx q[1];
rz(-2.254199) q[1];
sx q[1];
rz(-1.618209) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65255717) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(1.0126589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(-2.771647) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870435) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1856954) q[0];
sx q[0];
rz(-2.9736608) q[0];
sx q[0];
rz(2.3737565) q[0];
rz(2.3751971) q[2];
sx q[2];
rz(-0.99457127) q[2];
sx q[2];
rz(-0.11608427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.373133) q[1];
sx q[1];
rz(-2.2420954) q[1];
sx q[1];
rz(-2.9292604) q[1];
rz(-3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(-2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-0.20239057) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(2.6828962) q[2];
sx q[2];
rz(-1.9445322) q[2];
sx q[2];
rz(-1.4945488) q[2];
rz(0.17584569) q[3];
sx q[3];
rz(-2.1538215) q[3];
sx q[3];
rz(0.62721696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
