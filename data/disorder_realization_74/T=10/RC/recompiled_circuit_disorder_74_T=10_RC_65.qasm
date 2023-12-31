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
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(1.8523822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50531189) q[0];
sx q[0];
rz(-1.2828151) q[0];
sx q[0];
rz(2.050839) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1172328) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(-2.9565405) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2002441) q[1];
sx q[1];
rz(-1.3625047) q[1];
sx q[1];
rz(1.8336481) q[1];
rz(-pi) q[2];
rz(1.3025769) q[3];
sx q[3];
rz(-1.6383639) q[3];
sx q[3];
rz(-0.64642954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5916799) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(3.0060449) q[2];
rz(0.1401976) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.3134726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75269964) q[0];
sx q[0];
rz(-1.3503195) q[0];
sx q[0];
rz(-1.3541143) q[0];
rz(-pi) q[1];
rz(2.8333089) q[2];
sx q[2];
rz(-2.0364025) q[2];
sx q[2];
rz(0.75454933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(0.1608506) q[1];
x q[2];
rz(0.075606451) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(0.82270634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(-1.4387087) q[3];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-0.16608873) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040708) q[0];
sx q[0];
rz(-1.4253758) q[0];
sx q[0];
rz(1.8619596) q[0];
rz(-2.0735998) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(1.8333679) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3661763) q[1];
sx q[1];
rz(-1.1300547) q[1];
sx q[1];
rz(0.52904769) q[1];
x q[2];
rz(-1.8649678) q[3];
sx q[3];
rz(-0.73495451) q[3];
sx q[3];
rz(-0.24925772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(0.44149533) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(2.3667483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4962909) q[0];
sx q[0];
rz(-0.075356396) q[0];
sx q[0];
rz(0.2199748) q[0];
x q[1];
rz(-0.36802937) q[2];
sx q[2];
rz(-1.6124915) q[2];
sx q[2];
rz(2.7421943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1244509) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(0.084852858) q[1];
x q[2];
rz(-2.1500548) q[3];
sx q[3];
rz(-0.30617985) q[3];
sx q[3];
rz(-1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0094770771) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-2.9602125) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(3.0338874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62892249) q[0];
sx q[0];
rz(-1.218601) q[0];
sx q[0];
rz(-3.0432426) q[0];
rz(-pi) q[1];
rz(-0.92685076) q[2];
sx q[2];
rz(-2.0276208) q[2];
sx q[2];
rz(-2.9608375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50385034) q[1];
sx q[1];
rz(-0.78617326) q[1];
sx q[1];
rz(-2.1450858) q[1];
rz(-pi) q[2];
rz(-1.0072458) q[3];
sx q[3];
rz(-1.8661024) q[3];
sx q[3];
rz(-1.9293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(2.1246134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3973629) q[0];
sx q[0];
rz(-1.1083535) q[0];
sx q[0];
rz(0.37325333) q[0];
x q[1];
rz(0.11634155) q[2];
sx q[2];
rz(-0.2025207) q[2];
sx q[2];
rz(2.4570738) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14156995) q[1];
sx q[1];
rz(-1.1035898) q[1];
sx q[1];
rz(-1.0228553) q[1];
rz(2.3853112) q[3];
sx q[3];
rz(-0.68538266) q[3];
sx q[3];
rz(1.3553728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54414526) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909661) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(1.6183052) q[1];
sx q[1];
rz(-0.73262501) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0680925) q[0];
sx q[0];
rz(-1.4788027) q[0];
sx q[0];
rz(1.6091225) q[0];
rz(-2.2336823) q[2];
sx q[2];
rz(-0.50953509) q[2];
sx q[2];
rz(-2.8223035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4525758) q[1];
sx q[1];
rz(-0.92454443) q[1];
sx q[1];
rz(2.8313682) q[1];
rz(-1.9795498) q[3];
sx q[3];
rz(-1.6533274) q[3];
sx q[3];
rz(0.85588928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.62961489) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(-3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(-1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(0.40922871) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361355) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(0.34988846) q[0];
x q[1];
rz(0.52493237) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(2.8033825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.65539) q[1];
sx q[1];
rz(-1.4185925) q[1];
sx q[1];
rz(-0.57032077) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8141361) q[3];
sx q[3];
rz(-2.0753324) q[3];
sx q[3];
rz(2.6017021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2239969) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(1.508629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89440896) q[0];
sx q[0];
rz(-1.7315454) q[0];
sx q[0];
rz(-0.87662351) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6998859) q[2];
sx q[2];
rz(-0.61868762) q[2];
sx q[2];
rz(2.8512851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29533169) q[1];
sx q[1];
rz(-1.8896855) q[1];
sx q[1];
rz(2.1178513) q[1];
rz(-pi) q[2];
rz(-1.6375332) q[3];
sx q[3];
rz(-2.4373694) q[3];
sx q[3];
rz(0.080554068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(1.0288303) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(2.57634) q[0];
rz(-0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5436343) q[0];
sx q[0];
rz(-1.5289474) q[0];
sx q[0];
rz(-2.6393487) q[0];
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
rz(-1.6227464) q[1];
sx q[1];
rz(-1.1559876) q[1];
sx q[1];
rz(0.82621375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8990717) q[3];
sx q[3];
rz(-2.4007113) q[3];
sx q[3];
rz(-1.5981984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-0.59664574) q[2];
rz(-2.6560442) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116466) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(1.8252774) q[2];
sx q[2];
rz(-1.6152059) q[2];
sx q[2];
rz(-0.16441328) q[2];
rz(-2.9968895) q[3];
sx q[3];
rz(-0.51454138) q[3];
sx q[3];
rz(0.71458057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
