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
rz(-1.6662958) q[0];
sx q[0];
rz(-1.8721606) q[0];
sx q[0];
rz(2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(0.7575922) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2494932) q[0];
sx q[0];
rz(-1.4612911) q[0];
sx q[0];
rz(3.0654869) q[0];
x q[1];
rz(1.2045317) q[2];
sx q[2];
rz(-2.0431113) q[2];
sx q[2];
rz(0.83845316) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1266844) q[1];
sx q[1];
rz(-1.1619301) q[1];
sx q[1];
rz(-1.8168628) q[1];
x q[2];
rz(-2.4663062) q[3];
sx q[3];
rz(-0.8199586) q[3];
sx q[3];
rz(0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0483094) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(-1.7993571) q[2];
rz(-1.3679158) q[3];
sx q[3];
rz(-1.0932837) q[3];
sx q[3];
rz(-1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9369478) q[0];
sx q[0];
rz(-2.1277675) q[0];
sx q[0];
rz(2.192705) q[0];
rz(-3.0401547) q[1];
sx q[1];
rz(-2.0815492) q[1];
sx q[1];
rz(2.2116275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5407911) q[0];
sx q[0];
rz(-1.7834429) q[0];
sx q[0];
rz(1.4669111) q[0];
x q[1];
rz(-3.1279966) q[2];
sx q[2];
rz(-0.50522035) q[2];
sx q[2];
rz(-1.9590953) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9252214) q[1];
sx q[1];
rz(-0.19268806) q[1];
sx q[1];
rz(0.57740258) q[1];
rz(-2.3306777) q[3];
sx q[3];
rz(-2.0457912) q[3];
sx q[3];
rz(2.7678633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(0.742221) q[2];
rz(0.46418515) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6716229) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(-1.6999014) q[0];
rz(2.9761159) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(-2.2742719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21737032) q[0];
sx q[0];
rz(-1.5708367) q[0];
sx q[0];
rz(-1.5702308) q[0];
rz(-pi) q[1];
rz(2.0763511) q[2];
sx q[2];
rz(-1.8095513) q[2];
sx q[2];
rz(-2.3468034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9359303) q[1];
sx q[1];
rz(-1.8962269) q[1];
sx q[1];
rz(-1.2742576) q[1];
rz(2.1380566) q[3];
sx q[3];
rz(-2.2972882) q[3];
sx q[3];
rz(-2.9851825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9637588) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(0.7589232) q[2];
rz(-0.80398503) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(-0.72833958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8300962) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-2.6575644) q[0];
rz(1.5361891) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(-1.4253634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7361476) q[0];
sx q[0];
rz(-0.63129866) q[0];
sx q[0];
rz(-0.59595705) q[0];
x q[1];
rz(2.8292848) q[2];
sx q[2];
rz(-1.7877525) q[2];
sx q[2];
rz(-1.7157451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5023574) q[1];
sx q[1];
rz(-2.4662152) q[1];
sx q[1];
rz(1.9278072) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68257733) q[3];
sx q[3];
rz(-2.5831476) q[3];
sx q[3];
rz(-0.48170939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.385685) q[2];
sx q[2];
rz(-1.0681095) q[2];
sx q[2];
rz(-0.29328406) q[2];
rz(0.048132345) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(2.8369246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3047979) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(-1.1037214) q[0];
rz(-0.91148218) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(0.10890659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2809362) q[0];
sx q[0];
rz(-1.0901648) q[0];
sx q[0];
rz(2.3970766) q[0];
rz(-pi) q[1];
rz(1.7519195) q[2];
sx q[2];
rz(-1.0430153) q[2];
sx q[2];
rz(0.15508791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5327235) q[1];
sx q[1];
rz(-0.91174284) q[1];
sx q[1];
rz(-0.20688914) q[1];
rz(-pi) q[2];
rz(2.6960228) q[3];
sx q[3];
rz(-0.62544367) q[3];
sx q[3];
rz(-1.3423962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(-2.8835127) q[2];
rz(2.7598925) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(-0.55735731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3968762) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(1.1599524) q[0];
rz(-2.5634735) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(2.8181308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2787298) q[0];
sx q[0];
rz(-2.5040309) q[0];
sx q[0];
rz(-0.80018534) q[0];
x q[1];
rz(2.696229) q[2];
sx q[2];
rz(-0.54071835) q[2];
sx q[2];
rz(3.0034844) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66993078) q[1];
sx q[1];
rz(-1.0012697) q[1];
sx q[1];
rz(-1.3013775) q[1];
rz(0.59185054) q[3];
sx q[3];
rz(-2.602792) q[3];
sx q[3];
rz(1.1065229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9385927) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(2.2933551) q[2];
rz(1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0271725) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(2.2681336) q[0];
rz(1.9050441) q[1];
sx q[1];
rz(-0.97573391) q[1];
sx q[1];
rz(0.083560856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3665584) q[0];
sx q[0];
rz(-1.8904443) q[0];
sx q[0];
rz(-1.699422) q[0];
rz(0.23891831) q[2];
sx q[2];
rz(-2.0878125) q[2];
sx q[2];
rz(-0.97415249) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8712743) q[1];
sx q[1];
rz(-1.1764929) q[1];
sx q[1];
rz(-2.32347) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6912937) q[3];
sx q[3];
rz(-0.45540998) q[3];
sx q[3];
rz(2.5728855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7053335) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(-1.3667038) q[2];
rz(-2.8691835) q[3];
sx q[3];
rz(-1.0206157) q[3];
sx q[3];
rz(-0.75850707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-2.2347539) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(2.2032264) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(2.3209007) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0291727) q[0];
sx q[0];
rz(-0.39390644) q[0];
sx q[0];
rz(-1.9256511) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9782045) q[2];
sx q[2];
rz(-2.738224) q[2];
sx q[2];
rz(2.5591116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8269315) q[1];
sx q[1];
rz(-2.1912327) q[1];
sx q[1];
rz(-0.13038306) q[1];
x q[2];
rz(1.197416) q[3];
sx q[3];
rz(-0.42694091) q[3];
sx q[3];
rz(-0.77746848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9390823) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(1.067777) q[2];
rz(1.1311401) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15602569) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(-0.40618968) q[0];
rz(1.4093026) q[1];
sx q[1];
rz(-2.1235762) q[1];
sx q[1];
rz(-1.0850151) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4403518) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(-1.3568272) q[0];
rz(-1.5753393) q[2];
sx q[2];
rz(-0.96049236) q[2];
sx q[2];
rz(1.6729599) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6500195) q[1];
sx q[1];
rz(-0.88188344) q[1];
sx q[1];
rz(-2.7116009) q[1];
rz(-0.75094838) q[3];
sx q[3];
rz(-0.72970684) q[3];
sx q[3];
rz(0.88353222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5857508) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(-2.2108868) q[2];
rz(-1.3961004) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(0.11955587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4569106) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(1.3071625) q[0];
rz(2.7556509) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-1.0345667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2292896) q[0];
sx q[0];
rz(-1.6737537) q[0];
sx q[0];
rz(1.3961755) q[0];
rz(-0.46073353) q[2];
sx q[2];
rz(-1.3855532) q[2];
sx q[2];
rz(-2.4590907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1423827) q[1];
sx q[1];
rz(-0.94005985) q[1];
sx q[1];
rz(-0.56908619) q[1];
rz(-pi) q[2];
rz(-0.12328445) q[3];
sx q[3];
rz(-2.3269231) q[3];
sx q[3];
rz(-0.88501634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6941541) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(-0.36699692) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3601396) q[0];
sx q[0];
rz(-1.1893138) q[0];
sx q[0];
rz(-1.0839533) q[0];
rz(-0.057859261) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(1.6690785) q[2];
sx q[2];
rz(-1.3986893) q[2];
sx q[2];
rz(-1.324011) q[2];
rz(0.41856159) q[3];
sx q[3];
rz(-2.3971315) q[3];
sx q[3];
rz(-2.286398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
