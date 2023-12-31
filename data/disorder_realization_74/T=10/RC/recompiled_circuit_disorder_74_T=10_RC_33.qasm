OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(-0.98288012) q[0];
sx q[0];
rz(-2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(1.8523822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122983) q[0];
sx q[0];
rz(-1.1120783) q[0];
sx q[0];
rz(0.32231583) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1172328) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(-0.18505219) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9413486) q[1];
sx q[1];
rz(-1.779088) q[1];
sx q[1];
rz(1.8336481) q[1];
rz(-pi) q[2];
rz(1.3207924) q[3];
sx q[3];
rz(-0.27640009) q[3];
sx q[3];
rz(-1.1652511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(-0.1401976) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(-0.78080368) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.3134726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27539) q[0];
sx q[0];
rz(-1.35944) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-0.30828373) q[2];
sx q[2];
rz(-2.0364025) q[2];
sx q[2];
rz(0.75454933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12570757) q[1];
sx q[1];
rz(-1.7310377) q[1];
sx q[1];
rz(-1.4829774) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6766657) q[3];
sx q[3];
rz(-1.4956116) q[3];
sx q[3];
rz(2.4014846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(0.16608873) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9244708) q[0];
sx q[0];
rz(-2.8170601) q[0];
sx q[0];
rz(2.0425509) q[0];
rz(-pi) q[1];
rz(2.0735998) q[2];
sx q[2];
rz(-2.1192126) q[2];
sx q[2];
rz(1.8333679) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6925466) q[1];
sx q[1];
rz(-2.0448301) q[1];
sx q[1];
rz(-2.0708592) q[1];
rz(-2.8852799) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35963905) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(0.77484432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.144865) q[0];
sx q[0];
rz(-1.554368) q[0];
sx q[0];
rz(-0.073547151) q[0];
rz(-pi) q[1];
rz(-1.61548) q[2];
sx q[2];
rz(-1.9384906) q[2];
sx q[2];
rz(-1.1553264) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91037726) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(-1.3740205) q[1];
x q[2];
rz(-2.1500548) q[3];
sx q[3];
rz(-2.8354128) q[3];
sx q[3];
rz(-1.3907719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.6397887) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(0.60896215) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(0.10770527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5126702) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(-0.098350071) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5905215) q[2];
sx q[2];
rz(-2.1398009) q[2];
sx q[2];
rz(-2.0713194) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50385034) q[1];
sx q[1];
rz(-0.78617326) q[1];
sx q[1];
rz(-0.99650683) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0884573) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(-2.351458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.6311197) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(2.1246134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74422979) q[0];
sx q[0];
rz(-1.1083535) q[0];
sx q[0];
rz(2.7683393) q[0];
rz(-3.0252511) q[2];
sx q[2];
rz(-2.939072) q[2];
sx q[2];
rz(-2.4570738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6974666) q[1];
sx q[1];
rz(-1.0870458) q[1];
sx q[1];
rz(2.6078348) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0820886) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54414526) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(0.63465676) q[2];
rz(-2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(-0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25062659) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-0.84386688) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-0.73262501) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49918136) q[0];
sx q[0];
rz(-1.5326323) q[0];
sx q[0];
rz(3.0495318) q[0];
x q[1];
rz(0.90791038) q[2];
sx q[2];
rz(-0.50953509) q[2];
sx q[2];
rz(-2.8223035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19983229) q[1];
sx q[1];
rz(-0.70711771) q[1];
sx q[1];
rz(1.9553528) q[1];
rz(-0.089902417) q[3];
sx q[3];
rz(-1.9780759) q[3];
sx q[3];
rz(-2.390993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62961489) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(-0.037242446) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(-2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.1036296) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(-1.1720852) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-0.40922871) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17657875) q[0];
sx q[0];
rz(-1.8518378) q[0];
sx q[0];
rz(-0.9126419) q[0];
rz(-1.8843083) q[2];
sx q[2];
rz(-2.0743309) q[2];
sx q[2];
rz(-1.7538278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9600535) q[1];
sx q[1];
rz(-1.0078733) q[1];
sx q[1];
rz(1.7510508) q[1];
x q[2];
rz(-0.41142923) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-3.0761449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(0.55465737) q[2];
rz(-2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.6329637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471837) q[0];
sx q[0];
rz(-1.7315454) q[0];
sx q[0];
rz(-0.87662351) q[0];
x q[1];
rz(2.1855418) q[2];
sx q[2];
rz(-1.4960669) q[2];
sx q[2];
rz(1.1751307) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0867882) q[1];
sx q[1];
rz(-2.0874223) q[1];
sx q[1];
rz(-0.36887849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.056592654) q[3];
sx q[3];
rz(-2.2731299) q[3];
sx q[3];
rz(-3.1346723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14279723) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(1.3302594) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-0.56525266) q[0];
rz(0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.7601097) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488659) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(3.0548274) q[0];
rz(-pi) q[1];
rz(2.9276766) q[2];
sx q[2];
rz(-2.6579318) q[2];
sx q[2];
rz(2.1733401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8379412) q[1];
sx q[1];
rz(-2.2397869) q[1];
sx q[1];
rz(2.6020357) q[1];
x q[2];
rz(-1.8990717) q[3];
sx q[3];
rz(-2.4007113) q[3];
sx q[3];
rz(-1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(-0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6116466) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(1.8687517) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.8252774) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(-1.4894555) q[3];
sx q[3];
rz(-2.079439) q[3];
sx q[3];
rz(0.54872201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
