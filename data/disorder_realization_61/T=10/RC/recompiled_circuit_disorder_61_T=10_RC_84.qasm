OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(-1.5925621) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(-0.54418286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.005851) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(-2.6430921) q[0];
rz(-1.5266225) q[2];
sx q[2];
rz(-1.6720811) q[2];
sx q[2];
rz(-0.12461187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4678573) q[1];
sx q[1];
rz(-2.2506672) q[1];
sx q[1];
rz(0.52635898) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5892026) q[3];
sx q[3];
rz(-3.1079709) q[3];
sx q[3];
rz(2.606483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0698174) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(1.7791746) q[2];
rz(3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64489275) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(-1.171296) q[0];
rz(0.21121875) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(-1.404095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1275741) q[0];
sx q[0];
rz(-1.9919792) q[0];
sx q[0];
rz(0.76571) q[0];
rz(-0.48268433) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(-0.54373103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.029570015) q[1];
sx q[1];
rz(-2.8575172) q[1];
sx q[1];
rz(-3.0594538) q[1];
x q[2];
rz(-2.8849765) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(-1.4790505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8319548) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(2.1976166) q[2];
rz(2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29207644) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.7180432) q[0];
rz(-0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(-2.5779285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(1.4865727) q[0];
rz(-pi) q[1];
rz(1.3450422) q[2];
sx q[2];
rz(-2.1588237) q[2];
sx q[2];
rz(3.0622481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7970265) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(2.7558541) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54791252) q[3];
sx q[3];
rz(-2.4439545) q[3];
sx q[3];
rz(1.9426949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50239572) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(-2.0813265) q[2];
rz(-3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.32325) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(1.5441783) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(0.65778041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7155899) q[0];
sx q[0];
rz(-1.3960739) q[0];
sx q[0];
rz(1.5402373) q[0];
x q[1];
rz(-0.31099702) q[2];
sx q[2];
rz(-1.1592835) q[2];
sx q[2];
rz(0.81931537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0187877) q[1];
sx q[1];
rz(-1.0446761) q[1];
sx q[1];
rz(-0.29463525) q[1];
x q[2];
rz(2.8834881) q[3];
sx q[3];
rz(-1.8542284) q[3];
sx q[3];
rz(-0.62844814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0458935) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(1.0632769) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7970153) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(2.1176594) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.468058) q[0];
sx q[0];
rz(-1.5412016) q[0];
sx q[0];
rz(-3.1248321) q[0];
x q[1];
rz(-2.3847694) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(-0.28048453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.353646) q[1];
sx q[1];
rz(-1.2632571) q[1];
sx q[1];
rz(0.31378444) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0179126) q[3];
sx q[3];
rz(-1.9607753) q[3];
sx q[3];
rz(0.19815138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91288599) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(-1.0181001) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.6927032) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(1.1692283) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(-2.8170524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.823846) q[0];
sx q[0];
rz(-2.4577603) q[0];
sx q[0];
rz(2.2800287) q[0];
rz(-pi) q[1];
rz(2.4915699) q[2];
sx q[2];
rz(-0.92009244) q[2];
sx q[2];
rz(1.0879773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4536344) q[1];
sx q[1];
rz(-1.0369685) q[1];
sx q[1];
rz(-0.028793528) q[1];
rz(-1.8683897) q[3];
sx q[3];
rz(-1.1614359) q[3];
sx q[3];
rz(-2.5149432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67363182) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(-1.8661873) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642898) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(-1.4087079) q[0];
rz(0.45577058) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.2021525) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3287438) q[0];
sx q[0];
rz(-0.7757196) q[0];
sx q[0];
rz(-2.1562955) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9583148) q[2];
sx q[2];
rz(-2.2829208) q[2];
sx q[2];
rz(0.14856635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0899883) q[1];
sx q[1];
rz(-0.12433908) q[1];
sx q[1];
rz(0.73595562) q[1];
rz(-pi) q[2];
rz(2.5382302) q[3];
sx q[3];
rz(-2.8238736) q[3];
sx q[3];
rz(0.41894223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1414286) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705567) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-1.113168) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(1.5323458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198062) q[0];
sx q[0];
rz(-2.0397423) q[0];
sx q[0];
rz(1.297834) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1148909) q[2];
sx q[2];
rz(-2.0156983) q[2];
sx q[2];
rz(-2.0632495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0678565) q[1];
sx q[1];
rz(-1.590775) q[1];
sx q[1];
rz(-1.1205792) q[1];
rz(-pi) q[2];
rz(-2.9701482) q[3];
sx q[3];
rz(-2.532634) q[3];
sx q[3];
rz(-2.3122548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.9130075) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(-1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-2.6429208) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(-1.1876855) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.2845576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4391543) q[0];
sx q[0];
rz(-2.2674773) q[0];
sx q[0];
rz(-2.7657763) q[0];
rz(-pi) q[1];
rz(-0.2662439) q[2];
sx q[2];
rz(-1.2028482) q[2];
sx q[2];
rz(-0.14012303) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43209546) q[1];
sx q[1];
rz(-2.2120683) q[1];
sx q[1];
rz(-1.0048559) q[1];
rz(1.5981538) q[3];
sx q[3];
rz(-1.9606855) q[3];
sx q[3];
rz(1.1297806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5294042) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(0.55076304) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(-1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(-1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(2.3619161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44256193) q[0];
sx q[0];
rz(-2.583722) q[0];
sx q[0];
rz(-0.2691707) q[0];
x q[1];
rz(-2.6780307) q[2];
sx q[2];
rz(-0.87052411) q[2];
sx q[2];
rz(2.314032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.838678) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(-1.360449) q[1];
rz(-pi) q[2];
rz(-0.46882792) q[3];
sx q[3];
rz(-0.83823293) q[3];
sx q[3];
rz(1.9459141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(1.54281) q[2];
rz(2.4370082) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-0.012044756) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-0.96991878) q[2];
sx q[2];
rz(-2.6261332) q[2];
sx q[2];
rz(0.020608227) q[2];
rz(-2.5034954) q[3];
sx q[3];
rz(-1.4641855) q[3];
sx q[3];
rz(-2.2104213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];