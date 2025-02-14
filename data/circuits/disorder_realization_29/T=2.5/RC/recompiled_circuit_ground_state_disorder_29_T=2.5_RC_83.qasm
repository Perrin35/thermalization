OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(-0.23393272) q[0];
rz(-2.2995931) q[1];
sx q[1];
rz(-1.41398) q[1];
sx q[1];
rz(-1.6230621) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5748225) q[0];
sx q[0];
rz(-1.8810417) q[0];
sx q[0];
rz(-0.53939087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9286381) q[2];
sx q[2];
rz(-1.1446272) q[2];
sx q[2];
rz(2.1338685) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17659345) q[1];
sx q[1];
rz(-2.1273344) q[1];
sx q[1];
rz(2.2189369) q[1];
x q[2];
rz(0.99125864) q[3];
sx q[3];
rz(-2.391444) q[3];
sx q[3];
rz(1.0446435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5412377) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(-2.3360628) q[2];
rz(0.39189288) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(2.384757) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517752) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(2.7031194) q[0];
rz(1.8665727) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(-3.0335887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94919449) q[0];
sx q[0];
rz(-0.11365232) q[0];
sx q[0];
rz(-1.3495007) q[0];
x q[1];
rz(-0.89620583) q[2];
sx q[2];
rz(-0.87088481) q[2];
sx q[2];
rz(2.7090379) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17441544) q[1];
sx q[1];
rz(-1.6545873) q[1];
sx q[1];
rz(1.8248952) q[1];
rz(-pi) q[2];
rz(-1.2233954) q[3];
sx q[3];
rz(-0.51651556) q[3];
sx q[3];
rz(-1.7361189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3652304) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(-0.11621172) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-2.3381086) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6666343) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-2.3057002) q[0];
rz(-1.5180786) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(1.2380884) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1657175) q[0];
sx q[0];
rz(-1.200186) q[0];
sx q[0];
rz(1.0281282) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.798271) q[2];
sx q[2];
rz(-1.803071) q[2];
sx q[2];
rz(-2.6417123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7459384) q[1];
sx q[1];
rz(-1.6881583) q[1];
sx q[1];
rz(-3.0836041) q[1];
rz(-2.8594871) q[3];
sx q[3];
rz(-1.2938525) q[3];
sx q[3];
rz(1.1462584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8572924) q[2];
sx q[2];
rz(-0.14480536) q[2];
sx q[2];
rz(0.71883744) q[2];
rz(-2.5607064) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(2.412793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6179287) q[0];
sx q[0];
rz(-1.5939465) q[0];
sx q[0];
rz(-2.0148328) q[0];
rz(1.843938) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(2.0984971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9595056) q[0];
sx q[0];
rz(-2.8694177) q[0];
sx q[0];
rz(-1.4578117) q[0];
x q[1];
rz(1.0342233) q[2];
sx q[2];
rz(-1.3541647) q[2];
sx q[2];
rz(2.72675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.265643) q[1];
sx q[1];
rz(-1.7247206) q[1];
sx q[1];
rz(1.8791844) q[1];
rz(-pi) q[2];
rz(-0.65535069) q[3];
sx q[3];
rz(-1.5028302) q[3];
sx q[3];
rz(-2.7790536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2780693) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(-2.1045904) q[2];
rz(-2.1841124) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0015513) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(-0.087604372) q[0];
rz(-2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(1.3137438) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4632516) q[0];
sx q[0];
rz(-1.8049362) q[0];
sx q[0];
rz(1.0749504) q[0];
rz(-1.6694267) q[2];
sx q[2];
rz(-1.5336138) q[2];
sx q[2];
rz(-1.1929026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.524595) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(1.7434381) q[1];
rz(-pi) q[2];
x q[2];
rz(2.143712) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(-2.5469766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4801415) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(-0.56337774) q[2];
rz(1.2207458) q[3];
sx q[3];
rz(-1.3910339) q[3];
sx q[3];
rz(0.63108546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1221984) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(3.1296545) q[0];
rz(-2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(-2.8964002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48680533) q[0];
sx q[0];
rz(-1.1410895) q[0];
sx q[0];
rz(0.91581099) q[0];
x q[1];
rz(-1.7952314) q[2];
sx q[2];
rz(-1.3185474) q[2];
sx q[2];
rz(-0.66649619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4575544) q[1];
sx q[1];
rz(-2.0054617) q[1];
sx q[1];
rz(2.5358389) q[1];
rz(-pi) q[2];
rz(1.2782966) q[3];
sx q[3];
rz(-1.9579325) q[3];
sx q[3];
rz(-1.2237924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7616854) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(-2.6109429) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(-1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.948792) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(2.8444994) q[0];
rz(1.5104712) q[1];
sx q[1];
rz(-2.511697) q[1];
sx q[1];
rz(0.43103257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1265701) q[0];
sx q[0];
rz(-1.3608785) q[0];
sx q[0];
rz(0.2951584) q[0];
rz(-pi) q[1];
rz(-2.9210429) q[2];
sx q[2];
rz(-2.9052581) q[2];
sx q[2];
rz(-0.55921171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4342793) q[1];
sx q[1];
rz(-0.61798862) q[1];
sx q[1];
rz(-0.98843482) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2934612) q[3];
sx q[3];
rz(-2.0803422) q[3];
sx q[3];
rz(0.58125416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83583528) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(-2.8744899) q[2];
rz(0.032020656) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(-0.10572461) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(2.5015976) q[0];
rz(-0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(1.4124426) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873358) q[0];
sx q[0];
rz(-2.1001513) q[0];
sx q[0];
rz(-1.0933149) q[0];
rz(-pi) q[1];
rz(-1.475263) q[2];
sx q[2];
rz(-0.49411202) q[2];
sx q[2];
rz(1.4612559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.392726) q[1];
sx q[1];
rz(-2.7905474) q[1];
sx q[1];
rz(-0.24715329) q[1];
x q[2];
rz(0.85573393) q[3];
sx q[3];
rz(-2.6011356) q[3];
sx q[3];
rz(-3.1202701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60014805) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(-0.42204648) q[2];
rz(-0.47354928) q[3];
sx q[3];
rz(-0.62255064) q[3];
sx q[3];
rz(2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710829) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(-1.3516082) q[0];
rz(0.31556684) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(-1.1531166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4723268) q[0];
sx q[0];
rz(-1.6231704) q[0];
sx q[0];
rz(-3.0906424) q[0];
x q[1];
rz(1.6193689) q[2];
sx q[2];
rz(-1.0359541) q[2];
sx q[2];
rz(-2.4984307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0259195) q[1];
sx q[1];
rz(-0.5957091) q[1];
sx q[1];
rz(-2.9133948) q[1];
rz(2.3117468) q[3];
sx q[3];
rz(-2.6072558) q[3];
sx q[3];
rz(2.9754834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6844668) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(-2.7678164) q[2];
rz(-1.9155546) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(-1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1821197) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(-1.8208338) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(2.6722867) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824911) q[0];
sx q[0];
rz(-1.6880205) q[0];
sx q[0];
rz(3.0646695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9515418) q[2];
sx q[2];
rz(-2.1441438) q[2];
sx q[2];
rz(-0.1694451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0476956) q[1];
sx q[1];
rz(-2.9710431) q[1];
sx q[1];
rz(-2.6490414) q[1];
rz(-pi) q[2];
rz(2.8402249) q[3];
sx q[3];
rz(-1.8068988) q[3];
sx q[3];
rz(-0.88168854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70539537) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(-1.9311284) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(-2.1657522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217011) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(2.2484491) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(-2.2371815) q[2];
sx q[2];
rz(-2.7795962) q[2];
sx q[2];
rz(-2.4301651) q[2];
rz(-2.1410971) q[3];
sx q[3];
rz(-1.9420997) q[3];
sx q[3];
rz(1.6122769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
