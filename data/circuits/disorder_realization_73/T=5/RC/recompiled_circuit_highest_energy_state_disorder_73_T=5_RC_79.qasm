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
rz(-2.7849164) q[0];
sx q[0];
rz(-1.4465605) q[0];
sx q[0];
rz(0.90176982) q[0];
rz(1.3321441) q[1];
sx q[1];
rz(-2.6982215) q[1];
sx q[1];
rz(2.3763357) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0888521) q[0];
sx q[0];
rz(-2.46457) q[0];
sx q[0];
rz(0.96816008) q[0];
rz(-pi) q[1];
rz(2.8985326) q[2];
sx q[2];
rz(-1.2420734) q[2];
sx q[2];
rz(0.49561938) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3972348) q[1];
sx q[1];
rz(-1.9089948) q[1];
sx q[1];
rz(1.5282791) q[1];
rz(1.9992932) q[3];
sx q[3];
rz(-2.836578) q[3];
sx q[3];
rz(-1.3950789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8454933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(-1.6695401) q[2];
rz(1.7146401) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(-2.6078687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832837) q[0];
sx q[0];
rz(-1.6338209) q[0];
sx q[0];
rz(0.84877745) q[0];
rz(0.035004184) q[1];
sx q[1];
rz(-1.9947546) q[1];
sx q[1];
rz(-0.95357198) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7824421) q[0];
sx q[0];
rz(-1.1138565) q[0];
sx q[0];
rz(0.63553973) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30470253) q[2];
sx q[2];
rz(-1.1977473) q[2];
sx q[2];
rz(2.3359131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63605308) q[1];
sx q[1];
rz(-1.2314774) q[1];
sx q[1];
rz(-2.9341594) q[1];
rz(-pi) q[2];
rz(-2.4783432) q[3];
sx q[3];
rz(-2.2683072) q[3];
sx q[3];
rz(2.8467941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9780875) q[2];
sx q[2];
rz(-1.7814025) q[2];
sx q[2];
rz(-0.84954849) q[2];
rz(-2.9761525) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(2.3918242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(-2.7396696) q[0];
rz(2.9882714) q[1];
sx q[1];
rz(-0.17686495) q[1];
sx q[1];
rz(-2.0053999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72969681) q[0];
sx q[0];
rz(-1.9703016) q[0];
sx q[0];
rz(-3.1285498) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48353124) q[2];
sx q[2];
rz(-1.3267967) q[2];
sx q[2];
rz(0.86770832) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.043667533) q[1];
sx q[1];
rz(-0.33722028) q[1];
sx q[1];
rz(-0.80101669) q[1];
rz(2.4506684) q[3];
sx q[3];
rz(-0.57583416) q[3];
sx q[3];
rz(1.3802647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0574657) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-3.086536) q[2];
rz(1.347524) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(-1.0042892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.7771626) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(-1.5740016) q[0];
rz(0.69152999) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(2.9561668) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381993) q[0];
sx q[0];
rz(-1.1739587) q[0];
sx q[0];
rz(-1.3893886) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48443227) q[2];
sx q[2];
rz(-1.9524756) q[2];
sx q[2];
rz(-2.9228314) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5925827) q[1];
sx q[1];
rz(-1.8235778) q[1];
sx q[1];
rz(-1.2493253) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7774277) q[3];
sx q[3];
rz(-2.1071599) q[3];
sx q[3];
rz(0.11851507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7376248) q[2];
sx q[2];
rz(-1.4289958) q[2];
sx q[2];
rz(3.1363623) q[2];
rz(2.7541006) q[3];
sx q[3];
rz(-0.5046851) q[3];
sx q[3];
rz(-1.8385878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72652793) q[0];
sx q[0];
rz(-0.98353493) q[0];
sx q[0];
rz(2.5256185) q[0];
rz(1.9527324) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(-0.51330769) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5145288) q[0];
sx q[0];
rz(-1.5634057) q[0];
sx q[0];
rz(1.2891931) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31536021) q[2];
sx q[2];
rz(-0.7509481) q[2];
sx q[2];
rz(2.2436004) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3860046) q[1];
sx q[1];
rz(-1.6428448) q[1];
sx q[1];
rz(-2.5064038) q[1];
rz(-pi) q[2];
rz(3.0232459) q[3];
sx q[3];
rz(-1.1062628) q[3];
sx q[3];
rz(1.1778414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1335699) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(0.7424773) q[2];
rz(1.6978469) q[3];
sx q[3];
rz(-1.7323114) q[3];
sx q[3];
rz(-1.1013364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3105069) q[0];
sx q[0];
rz(-2.0948912) q[0];
sx q[0];
rz(2.781784) q[0];
rz(2.2911435) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(2.6757619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1547928) q[0];
sx q[0];
rz(-1.7868817) q[0];
sx q[0];
rz(-2.297202) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76325046) q[2];
sx q[2];
rz(-1.4723198) q[2];
sx q[2];
rz(-1.5710448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7342074) q[1];
sx q[1];
rz(-1.3300597) q[1];
sx q[1];
rz(-1.6233141) q[1];
x q[2];
rz(-0.96420607) q[3];
sx q[3];
rz(-2.0584985) q[3];
sx q[3];
rz(-1.694198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3992074) q[2];
sx q[2];
rz(-1.1804322) q[2];
sx q[2];
rz(-2.7772389) q[2];
rz(-0.24389167) q[3];
sx q[3];
rz(-0.049592169) q[3];
sx q[3];
rz(-1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5704982) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(1.1660227) q[0];
rz(-1.0150389) q[1];
sx q[1];
rz(-0.53703419) q[1];
sx q[1];
rz(2.7551415) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84595118) q[0];
sx q[0];
rz(-1.6014454) q[0];
sx q[0];
rz(1.3611193) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3101686) q[2];
sx q[2];
rz(-1.8156689) q[2];
sx q[2];
rz(2.819811) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1931175) q[1];
sx q[1];
rz(-1.5453234) q[1];
sx q[1];
rz(-1.6465882) q[1];
rz(-pi) q[2];
rz(1.3644045) q[3];
sx q[3];
rz(-2.5897621) q[3];
sx q[3];
rz(1.7655552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7361136) q[2];
sx q[2];
rz(-1.6454641) q[2];
sx q[2];
rz(-2.9138937) q[2];
rz(1.0730216) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(0.29357114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2597518) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(0.42801273) q[0];
rz(-2.1233066) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(-2.2705618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6414531) q[0];
sx q[0];
rz(-2.347751) q[0];
sx q[0];
rz(2.3026233) q[0];
x q[1];
rz(2.9350014) q[2];
sx q[2];
rz(-1.2204514) q[2];
sx q[2];
rz(1.3526582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4531252) q[1];
sx q[1];
rz(-0.18632132) q[1];
sx q[1];
rz(2.2800355) q[1];
rz(2.4610041) q[3];
sx q[3];
rz(-1.1201123) q[3];
sx q[3];
rz(0.87928538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.009306) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(1.2786678) q[2];
rz(0.57847413) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(1.9875897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4671675) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(-2.5417852) q[0];
rz(-1.8425875) q[1];
sx q[1];
rz(-1.6103585) q[1];
sx q[1];
rz(-0.64186796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079499809) q[0];
sx q[0];
rz(-0.60264093) q[0];
sx q[0];
rz(0.63061611) q[0];
rz(-pi) q[1];
rz(0.14417837) q[2];
sx q[2];
rz(-2.3308809) q[2];
sx q[2];
rz(-1.415103) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.228926) q[1];
sx q[1];
rz(-2.4335981) q[1];
sx q[1];
rz(-0.028381746) q[1];
rz(-0.75662002) q[3];
sx q[3];
rz(-1.9709658) q[3];
sx q[3];
rz(-2.4979046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.32144) q[2];
sx q[2];
rz(-1.7165311) q[2];
sx q[2];
rz(0.610262) q[2];
rz(-0.52014703) q[3];
sx q[3];
rz(-1.0545694) q[3];
sx q[3];
rz(0.49351969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(-0.34058061) q[0];
sx q[0];
rz(-1.2435253) q[0];
sx q[0];
rz(-0.93801671) q[0];
rz(0.34171379) q[1];
sx q[1];
rz(-1.382117) q[1];
sx q[1];
rz(-0.15728532) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324824) q[0];
sx q[0];
rz(-2.1868262) q[0];
sx q[0];
rz(1.8094713) q[0];
rz(0.089265169) q[2];
sx q[2];
rz(-2.3815739) q[2];
sx q[2];
rz(-1.6940728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57676901) q[1];
sx q[1];
rz(-1.8722459) q[1];
sx q[1];
rz(-0.5524854) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8870542) q[3];
sx q[3];
rz(-1.1738281) q[3];
sx q[3];
rz(-2.1074668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21891317) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(-1.867713) q[2];
rz(0.0056565469) q[3];
sx q[3];
rz(-1.7001245) q[3];
sx q[3];
rz(-0.98102942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2904749) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(0.46650096) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(-0.90032719) q[2];
sx q[2];
rz(-1.1569958) q[2];
sx q[2];
rz(-0.3668084) q[2];
rz(-2.9644184) q[3];
sx q[3];
rz(-2.7037947) q[3];
sx q[3];
rz(-1.6537742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
