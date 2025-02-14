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
rz(2.8994695) q[0];
sx q[0];
rz(-2.4162633) q[0];
sx q[0];
rz(2.2348833) q[0];
rz(-0.46696219) q[1];
sx q[1];
rz(4.0443647) q[1];
sx q[1];
rz(9.180896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8236602) q[0];
sx q[0];
rz(-0.77046227) q[0];
sx q[0];
rz(2.5567358) q[0];
rz(2.2026625) q[2];
sx q[2];
rz(-2.1515232) q[2];
sx q[2];
rz(1.9730568) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4441057) q[1];
sx q[1];
rz(-0.19679697) q[1];
sx q[1];
rz(-0.37741952) q[1];
x q[2];
rz(-0.81772352) q[3];
sx q[3];
rz(-1.7264888) q[3];
sx q[3];
rz(0.45785357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68524086) q[2];
sx q[2];
rz(-0.46619236) q[2];
sx q[2];
rz(2.1166128) q[2];
rz(3.0023365) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(-2.9774184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76992947) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(1.2051693) q[0];
rz(1.6734164) q[1];
sx q[1];
rz(-1.877715) q[1];
sx q[1];
rz(1.7211627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98200765) q[0];
sx q[0];
rz(-2.3075342) q[0];
sx q[0];
rz(2.9247051) q[0];
x q[1];
rz(-1.3060644) q[2];
sx q[2];
rz(-1.6240623) q[2];
sx q[2];
rz(-0.85088986) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0542595) q[1];
sx q[1];
rz(-1.9078322) q[1];
sx q[1];
rz(-0.23566206) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51610701) q[3];
sx q[3];
rz(-0.81023216) q[3];
sx q[3];
rz(-0.50767553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9517407) q[2];
sx q[2];
rz(-2.9517089) q[2];
sx q[2];
rz(-2.3404549) q[2];
rz(2.7028911) q[3];
sx q[3];
rz(-1.2480241) q[3];
sx q[3];
rz(-1.1130822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7640215) q[0];
sx q[0];
rz(-2.8442597) q[0];
sx q[0];
rz(1.1391033) q[0];
rz(0.26690075) q[1];
sx q[1];
rz(-1.8903939) q[1];
sx q[1];
rz(-0.36531726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26317715) q[0];
sx q[0];
rz(-2.6179144) q[0];
sx q[0];
rz(-0.19864638) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2174805) q[2];
sx q[2];
rz(-2.036493) q[2];
sx q[2];
rz(0.52801029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68054861) q[1];
sx q[1];
rz(-2.7690384) q[1];
sx q[1];
rz(1.7303321) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86287873) q[3];
sx q[3];
rz(-0.85298698) q[3];
sx q[3];
rz(0.77117111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9109362) q[2];
sx q[2];
rz(-2.4012884) q[2];
sx q[2];
rz(3.018107) q[2];
rz(2.2423045) q[3];
sx q[3];
rz(-1.4495918) q[3];
sx q[3];
rz(1.9062769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3673636) q[0];
sx q[0];
rz(-1.6772567) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(-2.3956237) q[1];
sx q[1];
rz(-1.4624701) q[1];
sx q[1];
rz(-1.7568582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2413413) q[0];
sx q[0];
rz(-1.5129733) q[0];
sx q[0];
rz(0.034863254) q[0];
rz(-pi) q[1];
rz(1.8141995) q[2];
sx q[2];
rz(-2.1982773) q[2];
sx q[2];
rz(0.091077591) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31092584) q[1];
sx q[1];
rz(-1.7574903) q[1];
sx q[1];
rz(-2.5273538) q[1];
x q[2];
rz(-0.05352002) q[3];
sx q[3];
rz(-2.60703) q[3];
sx q[3];
rz(-0.10343119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7437462) q[2];
sx q[2];
rz(-2.3036239) q[2];
sx q[2];
rz(1.0154593) q[2];
rz(3.1189392) q[3];
sx q[3];
rz(-1.3709143) q[3];
sx q[3];
rz(-0.13815752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37462336) q[0];
sx q[0];
rz(-1.8648819) q[0];
sx q[0];
rz(-0.20481566) q[0];
rz(0.21518937) q[1];
sx q[1];
rz(-2.9209825) q[1];
sx q[1];
rz(-2.2664216) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272064) q[0];
sx q[0];
rz(-1.7191186) q[0];
sx q[0];
rz(-1.2640087) q[0];
x q[1];
rz(0.1540252) q[2];
sx q[2];
rz(-2.2999894) q[2];
sx q[2];
rz(-0.83524365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2332747) q[1];
sx q[1];
rz(-1.2728146) q[1];
sx q[1];
rz(-1.8012992) q[1];
x q[2];
rz(-0.14303) q[3];
sx q[3];
rz(-2.1969079) q[3];
sx q[3];
rz(-0.54007733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.508226) q[2];
sx q[2];
rz(-0.88819155) q[2];
sx q[2];
rz(1.4166895) q[2];
rz(-2.21777) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(-1.1317322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1205587) q[0];
sx q[0];
rz(-0.73796213) q[0];
sx q[0];
rz(2.293204) q[0];
rz(-1.5757164) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(-2.1956992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1740055) q[0];
sx q[0];
rz(-1.8552836) q[0];
sx q[0];
rz(-0.6993642) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20393238) q[2];
sx q[2];
rz(-2.3205767) q[2];
sx q[2];
rz(-1.7199788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1523244) q[1];
sx q[1];
rz(-0.75256189) q[1];
sx q[1];
rz(-2.4965733) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3729783) q[3];
sx q[3];
rz(-1.0879915) q[3];
sx q[3];
rz(1.609926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1762323) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(-1.7410834) q[2];
rz(1.6598353) q[3];
sx q[3];
rz(-1.3503431) q[3];
sx q[3];
rz(0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7401212) q[0];
sx q[0];
rz(-2.6417702) q[0];
sx q[0];
rz(1.2116145) q[0];
rz(-0.47677332) q[1];
sx q[1];
rz(-1.8649273) q[1];
sx q[1];
rz(1.635294) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7619373) q[0];
sx q[0];
rz(-2.5003644) q[0];
sx q[0];
rz(2.0564709) q[0];
rz(-pi) q[1];
rz(-1.2477307) q[2];
sx q[2];
rz(-1.808721) q[2];
sx q[2];
rz(-1.0987154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2989212) q[1];
sx q[1];
rz(-1.4475045) q[1];
sx q[1];
rz(2.4409632) q[1];
rz(-pi) q[2];
x q[2];
rz(2.659306) q[3];
sx q[3];
rz(-1.8332172) q[3];
sx q[3];
rz(-0.69604128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8022884) q[2];
sx q[2];
rz(-1.5172639) q[2];
sx q[2];
rz(2.9586672) q[2];
rz(-0.66169345) q[3];
sx q[3];
rz(-1.9293834) q[3];
sx q[3];
rz(-0.70484149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45109192) q[0];
sx q[0];
rz(-1.4649614) q[0];
sx q[0];
rz(2.9943384) q[0];
rz(2.9946949) q[1];
sx q[1];
rz(-0.92120996) q[1];
sx q[1];
rz(1.1796835) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3203293) q[0];
sx q[0];
rz(-2.593145) q[0];
sx q[0];
rz(-2.4571277) q[0];
rz(0.4381035) q[2];
sx q[2];
rz(-1.1313651) q[2];
sx q[2];
rz(2.8218486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1547566) q[1];
sx q[1];
rz(-0.97797003) q[1];
sx q[1];
rz(1.4378217) q[1];
rz(3.0249765) q[3];
sx q[3];
rz(-0.75748235) q[3];
sx q[3];
rz(-0.99441409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2008721) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(2.3737523) q[2];
rz(0.52669865) q[3];
sx q[3];
rz(-2.0898762) q[3];
sx q[3];
rz(-2.5832978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4621157) q[0];
sx q[0];
rz(-2.9556892) q[0];
sx q[0];
rz(1.0900849) q[0];
rz(1.3500301) q[1];
sx q[1];
rz(-1.7410024) q[1];
sx q[1];
rz(2.7166691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2239531) q[0];
sx q[0];
rz(-1.6275379) q[0];
sx q[0];
rz(-1.0992235) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.091179087) q[2];
sx q[2];
rz(-0.65185368) q[2];
sx q[2];
rz(0.21287795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7560727) q[1];
sx q[1];
rz(-1.913383) q[1];
sx q[1];
rz(2.6581061) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0941689) q[3];
sx q[3];
rz(-2.2096388) q[3];
sx q[3];
rz(2.1372319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2316042) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(-0.52499047) q[2];
rz(-1.7867583) q[3];
sx q[3];
rz(-2.2460008) q[3];
sx q[3];
rz(1.7790599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7056284) q[0];
sx q[0];
rz(-0.66052496) q[0];
sx q[0];
rz(3.0349162) q[0];
rz(2.0872929) q[1];
sx q[1];
rz(-1.7091457) q[1];
sx q[1];
rz(1.4073102) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808436) q[0];
sx q[0];
rz(-2.0404732) q[0];
sx q[0];
rz(0.90523714) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65347341) q[2];
sx q[2];
rz(-1.9536886) q[2];
sx q[2];
rz(-3.058397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52508721) q[1];
sx q[1];
rz(-0.27180755) q[1];
sx q[1];
rz(0.96801438) q[1];
rz(-0.65741567) q[3];
sx q[3];
rz(-2.2748526) q[3];
sx q[3];
rz(-2.2179066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6537031) q[2];
sx q[2];
rz(-2.6370878) q[2];
sx q[2];
rz(0.27296909) q[2];
rz(-0.87131396) q[3];
sx q[3];
rz(-1.6815691) q[3];
sx q[3];
rz(-0.13450809) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7731666) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(2.7586965) q[1];
sx q[1];
rz(-1.8538414) q[1];
sx q[1];
rz(-0.014898653) q[1];
rz(-2.0128925) q[2];
sx q[2];
rz(-2.0876035) q[2];
sx q[2];
rz(0.58135396) q[2];
rz(-3.0440108) q[3];
sx q[3];
rz(-0.80229901) q[3];
sx q[3];
rz(-2.7589519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
