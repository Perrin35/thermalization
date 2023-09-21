OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4607853) q[0];
sx q[0];
rz(-2.1587125) q[0];
sx q[0];
rz(2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(-1.2892105) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754817) q[0];
sx q[0];
rz(-2.5876343) q[0];
sx q[0];
rz(-2.1411116) q[0];
rz(3.1172328) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(0.18505219) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8569021) q[1];
sx q[1];
rz(-0.33387091) q[1];
sx q[1];
rz(0.88792172) q[1];
x q[2];
rz(-1.3207924) q[3];
sx q[3];
rz(-0.27640009) q[3];
sx q[3];
rz(1.1652511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(2.360789) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(1.82812) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.388893) q[0];
sx q[0];
rz(-1.3503195) q[0];
sx q[0];
rz(-1.7874784) q[0];
x q[1];
rz(-1.08554) q[2];
sx q[2];
rz(-1.8453571) q[2];
sx q[2];
rz(-0.95825125) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4310415) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(-0.1608506) q[1];
rz(-pi) q[2];
rz(1.4649269) q[3];
sx q[3];
rz(-1.4956116) q[3];
sx q[3];
rz(0.74010805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(0.95412811) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(-0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(2.9755039) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375219) q[0];
sx q[0];
rz(-1.4253758) q[0];
sx q[0];
rz(1.279633) q[0];
rz(-pi) q[1];
rz(2.4737349) q[2];
sx q[2];
rz(-2.4154426) q[2];
sx q[2];
rz(2.6454676) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7158311) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(-0.75158822) q[1];
rz(-0.2563128) q[3];
sx q[3];
rz(-2.267572) q[3];
sx q[3];
rz(-3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(-0.44149533) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-0.77484432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.144865) q[0];
sx q[0];
rz(-1.554368) q[0];
sx q[0];
rz(0.073547151) q[0];
rz(-pi) q[1];
rz(-0.36802937) q[2];
sx q[2];
rz(-1.6124915) q[2];
sx q[2];
rz(2.7421943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1244509) q[1];
sx q[1];
rz(-1.972773) q[1];
sx q[1];
rz(-3.0567398) q[1];
x q[2];
rz(-0.99153783) q[3];
sx q[3];
rz(-0.30617985) q[3];
sx q[3];
rz(1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(1.6397887) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-0.18138012) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-3.0338874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5126702) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(-0.098350071) q[0];
rz(-pi) q[1];
rz(-0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.2203072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6377423) q[1];
sx q[1];
rz(-0.78617326) q[1];
sx q[1];
rz(0.99650683) q[1];
x q[2];
rz(2.1343469) q[3];
sx q[3];
rz(-1.2754903) q[3];
sx q[3];
rz(-1.2122648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.084215) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(0.010820476) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(-1.5104729) q[0];
rz(-1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(2.1246134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4646343) q[0];
sx q[0];
rz(-2.5559253) q[0];
sx q[0];
rz(-2.2023489) q[0];
rz(-2.9404042) q[2];
sx q[2];
rz(-1.5474461) q[2];
sx q[2];
rz(0.77229283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.444126) q[1];
sx q[1];
rz(-2.0545469) q[1];
sx q[1];
rz(-0.53375785) q[1];
rz(1.059504) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(-0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-0.84386688) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49918136) q[0];
sx q[0];
rz(-1.6089604) q[0];
sx q[0];
rz(-0.092060815) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8104066) q[2];
sx q[2];
rz(-1.1761616) q[2];
sx q[2];
rz(-2.0926203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0724832) q[1];
sx q[1];
rz(-1.8169889) q[1];
sx q[1];
rz(2.2407131) q[1];
rz(-pi) q[2];
x q[2];
rz(0.089902417) q[3];
sx q[3];
rz(-1.1635167) q[3];
sx q[3];
rz(-2.390993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(2.7323639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054571) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(0.34988846) q[0];
rz(-pi) q[1];
rz(-2.6312469) q[2];
sx q[2];
rz(-2.5556457) q[2];
sx q[2];
rz(-2.3454391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48620263) q[1];
sx q[1];
rz(-1.4185925) q[1];
sx q[1];
rz(-0.57032077) q[1];
rz(-pi) q[2];
rz(1.8141361) q[3];
sx q[3];
rz(-2.0753324) q[3];
sx q[3];
rz(-0.53989053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(-0.0099442033) q[0];
rz(-2.0857281) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.508629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86635877) q[0];
sx q[0];
rz(-2.4320721) q[0];
sx q[0];
rz(-1.8190246) q[0];
x q[1];
rz(0.091392322) q[2];
sx q[2];
rz(-2.1835727) q[2];
sx q[2];
rz(2.6932655) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0548045) q[1];
sx q[1];
rz(-1.0541704) q[1];
sx q[1];
rz(0.36887849) q[1];
x q[2];
rz(-2.2739201) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(-1.5411351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(-0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.3814829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5436343) q[0];
sx q[0];
rz(-1.5289474) q[0];
sx q[0];
rz(2.6393487) q[0];
rz(1.4597458) q[2];
sx q[2];
rz(-2.0425218) q[2];
sx q[2];
rz(1.2088838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8379412) q[1];
sx q[1];
rz(-2.2397869) q[1];
sx q[1];
rz(2.6020357) q[1];
rz(-0.28678203) q[3];
sx q[3];
rz(-2.2638595) q[3];
sx q[3];
rz(1.9758488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6116466) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.8687517) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(0.045885432) q[2];
sx q[2];
rz(-1.3165717) q[2];
sx q[2];
rz(-1.7236621) q[2];
rz(2.6315401) q[3];
sx q[3];
rz(-1.6418213) q[3];
sx q[3];
rz(2.1591975) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];