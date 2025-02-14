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
rz(-1.2217628) q[0];
sx q[0];
rz(-1.7696932) q[0];
sx q[0];
rz(2.0576117) q[0];
rz(-0.78217512) q[1];
sx q[1];
rz(3.5909619) q[1];
sx q[1];
rz(10.211791) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0660853) q[0];
sx q[0];
rz(-1.721816) q[0];
sx q[0];
rz(-1.1780894) q[0];
rz(1.5151565) q[2];
sx q[2];
rz(-2.4650827) q[2];
sx q[2];
rz(-1.8970053) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1842332) q[1];
sx q[1];
rz(-0.562698) q[1];
sx q[1];
rz(1.301728) q[1];
rz(0.88605864) q[3];
sx q[3];
rz(-1.3027667) q[3];
sx q[3];
rz(-2.4259329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3419753) q[2];
sx q[2];
rz(-1.819333) q[2];
sx q[2];
rz(0.70023099) q[2];
rz(1.7796984) q[3];
sx q[3];
rz(-1.9340065) q[3];
sx q[3];
rz(-0.12968682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6561061) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(-1.19278) q[0];
rz(-0.15395173) q[1];
sx q[1];
rz(-0.62869453) q[1];
sx q[1];
rz(-1.2492294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.10478) q[0];
sx q[0];
rz(-1.5082772) q[0];
sx q[0];
rz(0.17351242) q[0];
rz(1.4973499) q[2];
sx q[2];
rz(-2.5030862) q[2];
sx q[2];
rz(-0.28624619) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71603644) q[1];
sx q[1];
rz(-1.8752943) q[1];
sx q[1];
rz(1.4861121) q[1];
x q[2];
rz(0.90460299) q[3];
sx q[3];
rz(-1.8663916) q[3];
sx q[3];
rz(-1.4669661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2092756) q[2];
sx q[2];
rz(-1.2195769) q[2];
sx q[2];
rz(-0.90927124) q[2];
rz(0.13898177) q[3];
sx q[3];
rz(-0.2551955) q[3];
sx q[3];
rz(-2.2709258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22055498) q[0];
sx q[0];
rz(-1.9746566) q[0];
sx q[0];
rz(-1.9535109) q[0];
rz(2.2782169) q[1];
sx q[1];
rz(-1.2438351) q[1];
sx q[1];
rz(-0.9695425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351417) q[0];
sx q[0];
rz(-1.4124845) q[0];
sx q[0];
rz(0.16198762) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1968422) q[2];
sx q[2];
rz(-1.8965169) q[2];
sx q[2];
rz(1.0870799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2565793) q[1];
sx q[1];
rz(-1.6672214) q[1];
sx q[1];
rz(-2.9060591) q[1];
rz(-pi) q[2];
rz(-0.57200498) q[3];
sx q[3];
rz(-2.1915132) q[3];
sx q[3];
rz(2.2748924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36677507) q[2];
sx q[2];
rz(-2.884951) q[2];
sx q[2];
rz(-2.4859264) q[2];
rz(1.6171148) q[3];
sx q[3];
rz(-1.7888125) q[3];
sx q[3];
rz(0.85339439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26576385) q[0];
sx q[0];
rz(-2.0142856) q[0];
sx q[0];
rz(1.458459) q[0];
rz(-1.726285) q[1];
sx q[1];
rz(-1.7067319) q[1];
sx q[1];
rz(-1.2687792) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44132933) q[0];
sx q[0];
rz(-2.1887795) q[0];
sx q[0];
rz(2.000314) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2007723) q[2];
sx q[2];
rz(-2.0220304) q[2];
sx q[2];
rz(-0.68945092) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3038588) q[1];
sx q[1];
rz(-1.7540252) q[1];
sx q[1];
rz(-2.2732123) q[1];
rz(-pi) q[2];
rz(-1.2015518) q[3];
sx q[3];
rz(-2.8974468) q[3];
sx q[3];
rz(-2.5084054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65712523) q[2];
sx q[2];
rz(-1.7871658) q[2];
sx q[2];
rz(-1.8787059) q[2];
rz(2.9386988) q[3];
sx q[3];
rz(-2.1094567) q[3];
sx q[3];
rz(1.3304905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6577067) q[0];
sx q[0];
rz(-2.0575476) q[0];
sx q[0];
rz(2.6935691) q[0];
rz(1.1605284) q[1];
sx q[1];
rz(-1.0803761) q[1];
sx q[1];
rz(-1.0763268) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1370476) q[0];
sx q[0];
rz(-0.8595312) q[0];
sx q[0];
rz(-2.8595631) q[0];
x q[1];
rz(-0.33319039) q[2];
sx q[2];
rz(-1.5234064) q[2];
sx q[2];
rz(1.9006001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2158958) q[1];
sx q[1];
rz(-1.5091245) q[1];
sx q[1];
rz(-0.036176763) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.18783) q[3];
sx q[3];
rz(-1.5873529) q[3];
sx q[3];
rz(-1.7615033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5356323) q[2];
sx q[2];
rz(-1.2781906) q[2];
sx q[2];
rz(-0.34206259) q[2];
rz(-2.0096807) q[3];
sx q[3];
rz(-1.071238) q[3];
sx q[3];
rz(2.7560077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015942052) q[0];
sx q[0];
rz(-2.5378939) q[0];
sx q[0];
rz(-1.1018671) q[0];
rz(2.8562538) q[1];
sx q[1];
rz(-0.97831786) q[1];
sx q[1];
rz(0.91011059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.343276) q[0];
sx q[0];
rz(-1.7119223) q[0];
sx q[0];
rz(-0.1144764) q[0];
rz(-1.455972) q[2];
sx q[2];
rz(-1.4779203) q[2];
sx q[2];
rz(-2.6367413) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2294096) q[1];
sx q[1];
rz(-1.6278815) q[1];
sx q[1];
rz(-2.014671) q[1];
x q[2];
rz(-1.7496787) q[3];
sx q[3];
rz(-1.4174882) q[3];
sx q[3];
rz(-1.8453423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3204699) q[2];
sx q[2];
rz(-1.7197101) q[2];
sx q[2];
rz(2.087743) q[2];
rz(-2.7779135) q[3];
sx q[3];
rz(-1.4854919) q[3];
sx q[3];
rz(-2.6065839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66154552) q[0];
sx q[0];
rz(-1.0129901) q[0];
sx q[0];
rz(2.9841828) q[0];
rz(2.5120381) q[1];
sx q[1];
rz(-1.5903571) q[1];
sx q[1];
rz(-3.0622283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.019477) q[0];
sx q[0];
rz(-1.5003029) q[0];
sx q[0];
rz(-1.4792144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3936226) q[2];
sx q[2];
rz(-1.9594155) q[2];
sx q[2];
rz(1.3874229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0473898) q[1];
sx q[1];
rz(-1.6999287) q[1];
sx q[1];
rz(-2.6455621) q[1];
x q[2];
rz(1.5299876) q[3];
sx q[3];
rz(-0.58678484) q[3];
sx q[3];
rz(0.87848488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6307512) q[2];
sx q[2];
rz(-1.3151104) q[2];
sx q[2];
rz(-0.17244478) q[2];
rz(0.6423966) q[3];
sx q[3];
rz(-2.1746217) q[3];
sx q[3];
rz(-0.35014686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.0680577) q[0];
sx q[0];
rz(-1.302916) q[0];
sx q[0];
rz(0.56655836) q[0];
rz(-2.9134275) q[1];
sx q[1];
rz(-1.6888432) q[1];
sx q[1];
rz(-1.8295005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.909186) q[0];
sx q[0];
rz(-2.3433422) q[0];
sx q[0];
rz(-1.9203824) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0242537) q[2];
sx q[2];
rz(-0.84716958) q[2];
sx q[2];
rz(-0.1186419) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.560558) q[1];
sx q[1];
rz(-1.0315064) q[1];
sx q[1];
rz(-0.59625397) q[1];
x q[2];
rz(-1.9598403) q[3];
sx q[3];
rz(-1.5397289) q[3];
sx q[3];
rz(-0.16060747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18099004) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(0.20299882) q[2];
rz(2.3024043) q[3];
sx q[3];
rz(-2.8036717) q[3];
sx q[3];
rz(2.2904229) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508271) q[0];
sx q[0];
rz(-0.44742328) q[0];
sx q[0];
rz(-2.9946193) q[0];
rz(2.7087063) q[1];
sx q[1];
rz(-1.1914445) q[1];
sx q[1];
rz(-1.8870707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0085274335) q[0];
sx q[0];
rz(-1.0458032) q[0];
sx q[0];
rz(-1.9740482) q[0];
rz(-pi) q[1];
rz(-3.1077976) q[2];
sx q[2];
rz(-1.8189478) q[2];
sx q[2];
rz(-2.7128618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4476632) q[1];
sx q[1];
rz(-1.7433102) q[1];
sx q[1];
rz(1.2209784) q[1];
rz(2.9734475) q[3];
sx q[3];
rz(-1.5922308) q[3];
sx q[3];
rz(2.2529064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.36048421) q[2];
sx q[2];
rz(-1.0177178) q[2];
sx q[2];
rz(2.600889) q[2];
rz(-2.8402719) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(-1.8436684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28894579) q[0];
sx q[0];
rz(-1.1971373) q[0];
sx q[0];
rz(-2.4135015) q[0];
rz(-0.51746619) q[1];
sx q[1];
rz(-1.6518075) q[1];
sx q[1];
rz(2.1450086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1580799) q[0];
sx q[0];
rz(-1.5971113) q[0];
sx q[0];
rz(-1.9151494) q[0];
x q[1];
rz(1.4696737) q[2];
sx q[2];
rz(-0.99144672) q[2];
sx q[2];
rz(2.2183826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7000632) q[1];
sx q[1];
rz(-0.25694381) q[1];
sx q[1];
rz(-1.62613) q[1];
x q[2];
rz(0.49788614) q[3];
sx q[3];
rz(-2.9336946) q[3];
sx q[3];
rz(2.9836536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1835798) q[2];
sx q[2];
rz(-1.5778342) q[2];
sx q[2];
rz(-3.0713522) q[2];
rz(-0.0066512935) q[3];
sx q[3];
rz(-2.9700322) q[3];
sx q[3];
rz(0.10449617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2086647) q[0];
sx q[0];
rz(-1.7345971) q[0];
sx q[0];
rz(-1.4165184) q[0];
rz(2.3729462) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(1.6835536) q[2];
sx q[2];
rz(-1.5509477) q[2];
sx q[2];
rz(1.0890065) q[2];
rz(1.205659) q[3];
sx q[3];
rz(-1.8237999) q[3];
sx q[3];
rz(0.0013838125) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
