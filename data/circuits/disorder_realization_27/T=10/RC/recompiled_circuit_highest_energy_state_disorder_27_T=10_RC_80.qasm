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
rz(2.6746305) q[1];
sx q[1];
rz(-0.90277201) q[1];
sx q[1];
rz(0.24388193) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31793247) q[0];
sx q[0];
rz(-0.77046227) q[0];
sx q[0];
rz(0.58485683) q[0];
rz(-2.4087231) q[2];
sx q[2];
rz(-2.3112765) q[2];
sx q[2];
rz(0.24093369) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0599614) q[1];
sx q[1];
rz(-1.3880127) q[1];
sx q[1];
rz(-1.4974529) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81772352) q[3];
sx q[3];
rz(-1.7264888) q[3];
sx q[3];
rz(-0.45785357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68524086) q[2];
sx q[2];
rz(-2.6754003) q[2];
sx q[2];
rz(1.0249798) q[2];
rz(0.13925615) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(2.9774184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.76992947) q[0];
sx q[0];
rz(-1.0924082) q[0];
sx q[0];
rz(1.9364233) q[0];
rz(1.4681762) q[1];
sx q[1];
rz(-1.877715) q[1];
sx q[1];
rz(1.42043) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6997864) q[0];
sx q[0];
rz(-1.4107293) q[0];
sx q[0];
rz(-2.3193441) q[0];
x q[1];
rz(1.3060644) q[2];
sx q[2];
rz(-1.6240623) q[2];
sx q[2];
rz(0.85088986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0542595) q[1];
sx q[1];
rz(-1.2337605) q[1];
sx q[1];
rz(-2.9059306) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0923493) q[3];
sx q[3];
rz(-2.25246) q[3];
sx q[3];
rz(0.18100706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9517407) q[2];
sx q[2];
rz(-0.18988374) q[2];
sx q[2];
rz(0.80113775) q[2];
rz(-0.43870157) q[3];
sx q[3];
rz(-1.8935685) q[3];
sx q[3];
rz(1.1130822) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7640215) q[0];
sx q[0];
rz(-0.29733297) q[0];
sx q[0];
rz(2.0024894) q[0];
rz(-0.26690075) q[1];
sx q[1];
rz(-1.8903939) q[1];
sx q[1];
rz(0.36531726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26317715) q[0];
sx q[0];
rz(-0.52367822) q[0];
sx q[0];
rz(-2.9429463) q[0];
x q[1];
rz(-2.2659698) q[2];
sx q[2];
rz(-0.77689028) q[2];
sx q[2];
rz(0.50621835) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.461044) q[1];
sx q[1];
rz(-2.7690384) q[1];
sx q[1];
rz(1.7303321) q[1];
rz(0.85478304) q[3];
sx q[3];
rz(-1.0589335) q[3];
sx q[3];
rz(-1.8292242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2306564) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(3.018107) q[2];
rz(-2.2423045) q[3];
sx q[3];
rz(-1.4495918) q[3];
sx q[3];
rz(1.2353157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-1.7742291) q[0];
sx q[0];
rz(-1.4643359) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(0.74596897) q[1];
sx q[1];
rz(-1.6791226) q[1];
sx q[1];
rz(-1.3847345) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2413413) q[0];
sx q[0];
rz(-1.5129733) q[0];
sx q[0];
rz(-0.034863254) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3273932) q[2];
sx q[2];
rz(-2.1982773) q[2];
sx q[2];
rz(0.091077591) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3900323) q[1];
sx q[1];
rz(-2.1728325) q[1];
sx q[1];
rz(1.3436418) q[1];
rz(-2.607658) q[3];
sx q[3];
rz(-1.5980532) q[3];
sx q[3];
rz(1.7202924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7437462) q[2];
sx q[2];
rz(-2.3036239) q[2];
sx q[2];
rz(1.0154593) q[2];
rz(-3.1189392) q[3];
sx q[3];
rz(-1.3709143) q[3];
sx q[3];
rz(0.13815752) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7669693) q[0];
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
rz(-0.78007215) q[0];
sx q[0];
rz(-2.8018537) q[0];
sx q[0];
rz(-1.1113412) q[0];
rz(-1.4007447) q[2];
sx q[2];
rz(-2.399246) q[2];
sx q[2];
rz(0.60630732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4103247) q[1];
sx q[1];
rz(-1.350623) q[1];
sx q[1];
rz(-2.8360044) q[1];
rz(-pi) q[2];
rz(2.9985626) q[3];
sx q[3];
rz(-0.9446848) q[3];
sx q[3];
rz(-2.6015153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.508226) q[2];
sx q[2];
rz(-0.88819155) q[2];
sx q[2];
rz(-1.4166895) q[2];
rz(-0.9238227) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(1.1317322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1205587) q[0];
sx q[0];
rz(-0.73796213) q[0];
sx q[0];
rz(-2.293204) q[0];
rz(1.5757164) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(-0.94589344) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96758715) q[0];
sx q[0];
rz(-1.8552836) q[0];
sx q[0];
rz(2.4422285) q[0];
rz(-1.7849493) q[2];
sx q[2];
rz(-2.3698167) q[2];
sx q[2];
rz(2.0146148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1523244) q[1];
sx q[1];
rz(-2.3890308) q[1];
sx q[1];
rz(-2.4965733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35870798) q[3];
sx q[3];
rz(-0.51878319) q[3];
sx q[3];
rz(-2.0174873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1762323) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(-1.7410834) q[2];
rz(-1.6598353) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4014715) q[0];
sx q[0];
rz(-0.49982247) q[0];
sx q[0];
rz(-1.9299782) q[0];
rz(0.47677332) q[1];
sx q[1];
rz(-1.2766653) q[1];
sx q[1];
rz(-1.5062987) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9325628) q[0];
sx q[0];
rz(-1.2878006) q[0];
sx q[0];
rz(0.98732938) q[0];
x q[1];
rz(-2.2231323) q[2];
sx q[2];
rz(-2.7428584) q[2];
sx q[2];
rz(0.14130558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.41667485) q[1];
sx q[1];
rz(-2.43201) q[1];
sx q[1];
rz(-2.9516944) q[1];
rz(-0.52495082) q[3];
sx q[3];
rz(-0.54407507) q[3];
sx q[3];
rz(-1.3348188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3393042) q[2];
sx q[2];
rz(-1.5172639) q[2];
sx q[2];
rz(-2.9586672) q[2];
rz(-2.4798992) q[3];
sx q[3];
rz(-1.9293834) q[3];
sx q[3];
rz(0.70484149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45109192) q[0];
sx q[0];
rz(-1.6766312) q[0];
sx q[0];
rz(-0.14725421) q[0];
rz(-2.9946949) q[1];
sx q[1];
rz(-2.2203827) q[1];
sx q[1];
rz(1.1796835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3203293) q[0];
sx q[0];
rz(-2.593145) q[0];
sx q[0];
rz(-2.4571277) q[0];
rz(-pi) q[1];
rz(0.83663655) q[2];
sx q[2];
rz(-0.61021611) q[2];
sx q[2];
rz(1.9883375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1547566) q[1];
sx q[1];
rz(-0.97797003) q[1];
sx q[1];
rz(1.4378217) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75408077) q[3];
sx q[3];
rz(-1.6508266) q[3];
sx q[3];
rz(2.4802993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2008721) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(-0.76784039) q[2];
rz(-0.52669865) q[3];
sx q[3];
rz(-2.0898762) q[3];
sx q[3];
rz(2.5832978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4621157) q[0];
sx q[0];
rz(-2.9556892) q[0];
sx q[0];
rz(2.0515077) q[0];
rz(-1.7915626) q[1];
sx q[1];
rz(-1.4005902) q[1];
sx q[1];
rz(0.42492351) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9176396) q[0];
sx q[0];
rz(-1.5140547) q[0];
sx q[0];
rz(1.0992235) q[0];
x q[1];
rz(-3.0504136) q[2];
sx q[2];
rz(-0.65185368) q[2];
sx q[2];
rz(2.9287147) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35985816) q[1];
sx q[1];
rz(-1.1175851) q[1];
sx q[1];
rz(-1.9537326) q[1];
x q[2];
rz(-2.0941689) q[3];
sx q[3];
rz(-0.93195385) q[3];
sx q[3];
rz(1.0043608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2316042) q[2];
sx q[2];
rz(-2.6308172) q[2];
sx q[2];
rz(-2.6166022) q[2];
rz(1.7867583) q[3];
sx q[3];
rz(-0.89559186) q[3];
sx q[3];
rz(-1.3625328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.43596426) q[0];
sx q[0];
rz(-2.4810677) q[0];
sx q[0];
rz(0.10667644) q[0];
rz(2.0872929) q[1];
sx q[1];
rz(-1.7091457) q[1];
sx q[1];
rz(-1.7342825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6808436) q[0];
sx q[0];
rz(-2.0404732) q[0];
sx q[0];
rz(-0.90523714) q[0];
rz(-pi) q[1];
rz(2.0402526) q[2];
sx q[2];
rz(-2.1699274) q[2];
sx q[2];
rz(1.208973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5104502) q[1];
sx q[1];
rz(-1.7235974) q[1];
sx q[1];
rz(-1.3451206) q[1];
rz(-pi) q[2];
rz(0.75026433) q[3];
sx q[3];
rz(-2.0553053) q[3];
sx q[3];
rz(1.110525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4878896) q[2];
sx q[2];
rz(-0.50450486) q[2];
sx q[2];
rz(-0.27296909) q[2];
rz(0.87131396) q[3];
sx q[3];
rz(-1.4600236) q[3];
sx q[3];
rz(-0.13450809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7731666) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(0.38289616) q[1];
sx q[1];
rz(-1.2877512) q[1];
sx q[1];
rz(3.126694) q[1];
rz(-2.5802774) q[2];
sx q[2];
rz(-1.9519162) q[2];
sx q[2];
rz(-1.219195) q[2];
rz(1.6712345) q[3];
sx q[3];
rz(-2.3681869) q[3];
sx q[3];
rz(-2.8988732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
