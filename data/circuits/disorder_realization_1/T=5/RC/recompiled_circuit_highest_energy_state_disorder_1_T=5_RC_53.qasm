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
rz(-0.048592351) q[0];
sx q[0];
rz(4.5669849) q[0];
sx q[0];
rz(9.1591747) q[0];
rz(-0.44759294) q[1];
sx q[1];
rz(4.6456479) q[1];
sx q[1];
rz(9.3001443) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7974313) q[0];
sx q[0];
rz(-2.0051885) q[0];
sx q[0];
rz(2.2940367) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9452478) q[2];
sx q[2];
rz(-1.4513029) q[2];
sx q[2];
rz(1.4422642) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9198299) q[1];
sx q[1];
rz(-1.3897093) q[1];
sx q[1];
rz(-2.9403375) q[1];
rz(0.83104317) q[3];
sx q[3];
rz(-2.9141015) q[3];
sx q[3];
rz(0.25615197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0894185) q[2];
sx q[2];
rz(-2.3201421) q[2];
sx q[2];
rz(-2.2649412) q[2];
rz(-0.18757251) q[3];
sx q[3];
rz(-0.98637527) q[3];
sx q[3];
rz(3.0228289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12980421) q[0];
sx q[0];
rz(-3.0846444) q[0];
sx q[0];
rz(-2.7563128) q[0];
rz(2.724559) q[1];
sx q[1];
rz(-2.4200491) q[1];
sx q[1];
rz(-1.0122274) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6979264) q[0];
sx q[0];
rz(-1.3529643) q[0];
sx q[0];
rz(2.1429254) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1181548) q[2];
sx q[2];
rz(-2.1760245) q[2];
sx q[2];
rz(-2.6626448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23823638) q[1];
sx q[1];
rz(-2.2975249) q[1];
sx q[1];
rz(1.958117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6777534) q[3];
sx q[3];
rz(-2.1968004) q[3];
sx q[3];
rz(-2.8371365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1483083) q[2];
sx q[2];
rz(-2.7316284) q[2];
sx q[2];
rz(0.5245463) q[2];
rz(2.2927393) q[3];
sx q[3];
rz(-0.69791228) q[3];
sx q[3];
rz(-2.2822288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16070212) q[0];
sx q[0];
rz(-0.41146678) q[0];
sx q[0];
rz(2.8149862) q[0];
rz(1.3702565) q[1];
sx q[1];
rz(-2.0854918) q[1];
sx q[1];
rz(0.036570963) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9243226) q[0];
sx q[0];
rz(-2.0170209) q[0];
sx q[0];
rz(-1.9666332) q[0];
rz(2.6447949) q[2];
sx q[2];
rz(-0.76947509) q[2];
sx q[2];
rz(2.0747607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83710872) q[1];
sx q[1];
rz(-1.3837985) q[1];
sx q[1];
rz(0.043882462) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0053592) q[3];
sx q[3];
rz(-0.32945368) q[3];
sx q[3];
rz(1.1245331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44564351) q[2];
sx q[2];
rz(-0.30569884) q[2];
sx q[2];
rz(1.1841527) q[2];
rz(0.46237692) q[3];
sx q[3];
rz(-1.0824883) q[3];
sx q[3];
rz(-2.6760127) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1831803) q[0];
sx q[0];
rz(-1.4875702) q[0];
sx q[0];
rz(0.91444772) q[0];
rz(0.86638266) q[1];
sx q[1];
rz(-1.2976846) q[1];
sx q[1];
rz(-0.31747174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4015071) q[0];
sx q[0];
rz(-2.751787) q[0];
sx q[0];
rz(0.29081776) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9725665) q[2];
sx q[2];
rz(-1.9155972) q[2];
sx q[2];
rz(1.2107236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7496192) q[1];
sx q[1];
rz(-1.8388868) q[1];
sx q[1];
rz(-0.15971649) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5741979) q[3];
sx q[3];
rz(-2.7614584) q[3];
sx q[3];
rz(-1.326732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0471961) q[2];
sx q[2];
rz(-2.6698343) q[2];
sx q[2];
rz(2.6867234) q[2];
rz(1.9589849) q[3];
sx q[3];
rz(-2.9136361) q[3];
sx q[3];
rz(-3.0003149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126548) q[0];
sx q[0];
rz(-1.8998572) q[0];
sx q[0];
rz(0.051359635) q[0];
rz(-0.32737577) q[1];
sx q[1];
rz(-0.86762571) q[1];
sx q[1];
rz(3.0816269) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7236008) q[0];
sx q[0];
rz(-1.2381885) q[0];
sx q[0];
rz(0.88653112) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4879901) q[2];
sx q[2];
rz(-1.8817668) q[2];
sx q[2];
rz(2.5429436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8096749) q[1];
sx q[1];
rz(-2.1109739) q[1];
sx q[1];
rz(-0.36560161) q[1];
rz(-pi) q[2];
rz(-0.19503212) q[3];
sx q[3];
rz(-0.86544207) q[3];
sx q[3];
rz(-2.5759199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3943693) q[2];
sx q[2];
rz(-0.9129492) q[2];
sx q[2];
rz(2.8933914) q[2];
rz(-0.40211755) q[3];
sx q[3];
rz(-2.6995903) q[3];
sx q[3];
rz(0.88002747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0070655951) q[0];
sx q[0];
rz(-1.2638673) q[0];
sx q[0];
rz(1.7279351) q[0];
rz(1.2029348) q[1];
sx q[1];
rz(-0.92034942) q[1];
sx q[1];
rz(-3.038182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0132842) q[0];
sx q[0];
rz(-3.04848) q[0];
sx q[0];
rz(0.72269209) q[0];
x q[1];
rz(-0.31423435) q[2];
sx q[2];
rz(-1.3948007) q[2];
sx q[2];
rz(0.90547784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6573665) q[1];
sx q[1];
rz(-1.5897337) q[1];
sx q[1];
rz(-3.0168946) q[1];
rz(-pi) q[2];
rz(-0.2528905) q[3];
sx q[3];
rz(-1.9441981) q[3];
sx q[3];
rz(1.4772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6695413) q[2];
sx q[2];
rz(-2.5504888) q[2];
sx q[2];
rz(0.21311398) q[2];
rz(1.6040246) q[3];
sx q[3];
rz(-3.0209916) q[3];
sx q[3];
rz(-0.11046256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14385001) q[0];
sx q[0];
rz(-0.85875964) q[0];
sx q[0];
rz(-0.48458883) q[0];
rz(-0.46802014) q[1];
sx q[1];
rz(-1.3222398) q[1];
sx q[1];
rz(3.0048634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28324652) q[0];
sx q[0];
rz(-1.8112881) q[0];
sx q[0];
rz(-1.1838311) q[0];
rz(0.71857326) q[2];
sx q[2];
rz(-2.6426275) q[2];
sx q[2];
rz(-2.8091407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51653105) q[1];
sx q[1];
rz(-1.704055) q[1];
sx q[1];
rz(-0.85162152) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4522631) q[3];
sx q[3];
rz(-1.9876936) q[3];
sx q[3];
rz(1.1030918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7318657) q[2];
sx q[2];
rz(-3.0146283) q[2];
sx q[2];
rz(2.6597888) q[2];
rz(1.2305772) q[3];
sx q[3];
rz(-1.9989719) q[3];
sx q[3];
rz(-3.0060911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9442673) q[0];
sx q[0];
rz(-2.8062286) q[0];
sx q[0];
rz(-2.770597) q[0];
rz(-0.15759298) q[1];
sx q[1];
rz(-1.3875562) q[1];
sx q[1];
rz(-2.2612459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1938095) q[0];
sx q[0];
rz(-1.0968465) q[0];
sx q[0];
rz(-3.0633846) q[0];
rz(-pi) q[1];
rz(-1.738638) q[2];
sx q[2];
rz(-0.8693822) q[2];
sx q[2];
rz(2.0244903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.62746) q[1];
sx q[1];
rz(-1.693502) q[1];
sx q[1];
rz(0.0038599188) q[1];
rz(-pi) q[2];
rz(1.861946) q[3];
sx q[3];
rz(-1.2978122) q[3];
sx q[3];
rz(-0.11364869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2788435) q[2];
sx q[2];
rz(-2.4767488) q[2];
sx q[2];
rz(-2.1319907) q[2];
rz(-0.75262117) q[3];
sx q[3];
rz(-1.8372583) q[3];
sx q[3];
rz(0.93839222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2272334) q[0];
sx q[0];
rz(-0.14150134) q[0];
sx q[0];
rz(-3.09521) q[0];
rz(1.6498238) q[1];
sx q[1];
rz(-1.3834508) q[1];
sx q[1];
rz(0.47328624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9289439) q[0];
sx q[0];
rz(-1.6814274) q[0];
sx q[0];
rz(-1.380019) q[0];
rz(2.1572596) q[2];
sx q[2];
rz(-0.55078673) q[2];
sx q[2];
rz(3.0760461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9732194) q[1];
sx q[1];
rz(-0.96046093) q[1];
sx q[1];
rz(1.8944505) q[1];
x q[2];
rz(1.5978053) q[3];
sx q[3];
rz(-0.031899769) q[3];
sx q[3];
rz(-0.75237319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16372323) q[2];
sx q[2];
rz(-2.551584) q[2];
sx q[2];
rz(-0.02656492) q[2];
rz(0.50655347) q[3];
sx q[3];
rz(-2.3285464) q[3];
sx q[3];
rz(-2.626239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16515054) q[0];
sx q[0];
rz(-0.9557752) q[0];
sx q[0];
rz(2.9859848) q[0];
rz(-0.43442976) q[1];
sx q[1];
rz(-0.54672086) q[1];
sx q[1];
rz(0.2880407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91950509) q[0];
sx q[0];
rz(-2.1489369) q[0];
sx q[0];
rz(2.2116106) q[0];
rz(-pi) q[1];
rz(-1.4819809) q[2];
sx q[2];
rz(-0.15379158) q[2];
sx q[2];
rz(-1.8078992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1120243) q[1];
sx q[1];
rz(-0.82371897) q[1];
sx q[1];
rz(-2.1805641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1959201) q[3];
sx q[3];
rz(-1.9145619) q[3];
sx q[3];
rz(1.6812066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7274999) q[2];
sx q[2];
rz(-0.45656559) q[2];
sx q[2];
rz(-0.045819316) q[2];
rz(3.0231061) q[3];
sx q[3];
rz(-0.89209569) q[3];
sx q[3];
rz(-0.74046016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4988149) q[0];
sx q[0];
rz(-1.6294263) q[0];
sx q[0];
rz(-1.908041) q[0];
rz(-1.9998101) q[1];
sx q[1];
rz(-0.70527609) q[1];
sx q[1];
rz(-1.3377778) q[1];
rz(-0.6836239) q[2];
sx q[2];
rz(-1.0928921) q[2];
sx q[2];
rz(3.0047807) q[2];
rz(-1.508504) q[3];
sx q[3];
rz(-1.6332165) q[3];
sx q[3];
rz(-1.6467057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
