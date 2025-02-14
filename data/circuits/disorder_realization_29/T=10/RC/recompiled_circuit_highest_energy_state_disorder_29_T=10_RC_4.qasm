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
rz(0.88275498) q[0];
sx q[0];
rz(3.2864154) q[0];
sx q[0];
rz(10.176131) q[0];
rz(1.524628) q[1];
sx q[1];
rz(-2.1722062) q[1];
sx q[1];
rz(0.17925395) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0517901) q[0];
sx q[0];
rz(-1.8500916) q[0];
sx q[0];
rz(-1.0599049) q[0];
rz(-pi) q[1];
rz(-1.5380074) q[2];
sx q[2];
rz(-0.84814397) q[2];
sx q[2];
rz(-3.0422831) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36513953) q[1];
sx q[1];
rz(-1.498933) q[1];
sx q[1];
rz(2.8316569) q[1];
rz(1.6514014) q[3];
sx q[3];
rz(-0.63524073) q[3];
sx q[3];
rz(-2.1135534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2396607) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(-2.530976) q[2];
rz(-2.2527952) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(0.73386598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69773847) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(1.656172) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(-0.86404538) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5651503) q[0];
sx q[0];
rz(-0.2731384) q[0];
sx q[0];
rz(1.9594203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0338221) q[2];
sx q[2];
rz(-1.290375) q[2];
sx q[2];
rz(3.1299431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73693528) q[1];
sx q[1];
rz(-2.9184249) q[1];
sx q[1];
rz(0.88921247) q[1];
rz(0.78317376) q[3];
sx q[3];
rz(-2.0370954) q[3];
sx q[3];
rz(0.1649905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(0.19134276) q[2];
rz(0.30609104) q[3];
sx q[3];
rz(-1.4602665) q[3];
sx q[3];
rz(-2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3092344) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(-0.424463) q[0];
rz(-1.3407432) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(2.4901966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1196647) q[0];
sx q[0];
rz(-0.16379539) q[0];
sx q[0];
rz(-1.5271565) q[0];
x q[1];
rz(-2.7578951) q[2];
sx q[2];
rz(-1.4009152) q[2];
sx q[2];
rz(-0.023141247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0733903) q[1];
sx q[1];
rz(-1.7531839) q[1];
sx q[1];
rz(-0.50478151) q[1];
rz(-1.3440787) q[3];
sx q[3];
rz(-1.4170839) q[3];
sx q[3];
rz(0.24776974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1059619) q[2];
sx q[2];
rz(-1.1392081) q[2];
sx q[2];
rz(2.4448709) q[2];
rz(2.4524955) q[3];
sx q[3];
rz(-1.9870116) q[3];
sx q[3];
rz(-0.2564297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43831393) q[0];
sx q[0];
rz(-2.8499481) q[0];
sx q[0];
rz(2.1412204) q[0];
rz(1.6814303) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.2132852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1532261) q[0];
sx q[0];
rz(-1.8509764) q[0];
sx q[0];
rz(0.12616726) q[0];
rz(-0.72896616) q[2];
sx q[2];
rz(-1.9310968) q[2];
sx q[2];
rz(3.136022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69566947) q[1];
sx q[1];
rz(-1.6918105) q[1];
sx q[1];
rz(2.8112429) q[1];
rz(-pi) q[2];
rz(-0.090661006) q[3];
sx q[3];
rz(-1.034015) q[3];
sx q[3];
rz(1.8859552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0393684) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(0.064112045) q[2];
rz(-2.4168849) q[3];
sx q[3];
rz(-1.2084081) q[3];
sx q[3];
rz(2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940014) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-0.79175788) q[0];
rz(2.0296312) q[1];
sx q[1];
rz(-1.0589736) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9324786) q[0];
sx q[0];
rz(-2.1585585) q[0];
sx q[0];
rz(-2.5665119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9826775) q[2];
sx q[2];
rz(-2.5549485) q[2];
sx q[2];
rz(1.212709) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.89692749) q[1];
sx q[1];
rz(-1.7411971) q[1];
sx q[1];
rz(3.0651613) q[1];
rz(-1.8800635) q[3];
sx q[3];
rz(-0.59854186) q[3];
sx q[3];
rz(-0.33613294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2935334) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(1.1406356) q[2];
rz(-0.10661495) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(3.1383681) q[0];
rz(-3.0942753) q[1];
sx q[1];
rz(-1.4553757) q[1];
sx q[1];
rz(1.3501732) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605476) q[0];
sx q[0];
rz(-2.8821917) q[0];
sx q[0];
rz(1.3400643) q[0];
x q[1];
rz(2.5642942) q[2];
sx q[2];
rz(-0.98492981) q[2];
sx q[2];
rz(-2.455204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1803169) q[1];
sx q[1];
rz(-2.2928228) q[1];
sx q[1];
rz(-0.78603334) q[1];
rz(-pi) q[2];
rz(-0.030478625) q[3];
sx q[3];
rz(-1.7913831) q[3];
sx q[3];
rz(2.7762846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.99001592) q[2];
sx q[2];
rz(-2.9559957) q[2];
sx q[2];
rz(3.0604176) q[2];
rz(-0.89933991) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4261632) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(2.6370866) q[0];
rz(-2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(-0.02034932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74744862) q[0];
sx q[0];
rz(-1.6645303) q[0];
sx q[0];
rz(-1.9362437) q[0];
x q[1];
rz(-3.0145524) q[2];
sx q[2];
rz(-1.7774701) q[2];
sx q[2];
rz(-2.5673696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87096244) q[1];
sx q[1];
rz(-0.32901627) q[1];
sx q[1];
rz(-2.1966619) q[1];
x q[2];
rz(-0.048681569) q[3];
sx q[3];
rz(-1.657007) q[3];
sx q[3];
rz(1.8415368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(1.7669558) q[2];
rz(1.5291519) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(0.092718743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2328211) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(-0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(-0.94295162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97639192) q[0];
sx q[0];
rz(-1.8689253) q[0];
sx q[0];
rz(-2.0701029) q[0];
rz(-pi) q[1];
rz(0.23343691) q[2];
sx q[2];
rz(-2.1899411) q[2];
sx q[2];
rz(0.34282986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3642204) q[1];
sx q[1];
rz(-0.99507123) q[1];
sx q[1];
rz(-0.74933021) q[1];
rz(-pi) q[2];
rz(1.2528768) q[3];
sx q[3];
rz(-0.32264454) q[3];
sx q[3];
rz(1.8436197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1711787) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(1.1478434) q[2];
rz(0.36192274) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(-1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(-2.4107133) q[0];
sx q[0];
rz(-2.78237) q[0];
sx q[0];
rz(3.0391589) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(-0.817743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90392834) q[0];
sx q[0];
rz(-2.2127164) q[0];
sx q[0];
rz(-2.4556745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2417415) q[2];
sx q[2];
rz(-1.6616115) q[2];
sx q[2];
rz(-0.95041529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39108155) q[1];
sx q[1];
rz(-2.4121248) q[1];
sx q[1];
rz(-1.5424278) q[1];
rz(-pi) q[2];
rz(-0.93877403) q[3];
sx q[3];
rz(-1.9524487) q[3];
sx q[3];
rz(-1.8167239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7927336) q[2];
sx q[2];
rz(-2.0642955) q[2];
sx q[2];
rz(-0.31361541) q[2];
rz(2.6775728) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(-0.919842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2713276) q[0];
sx q[0];
rz(-0.79065228) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(-2.927921) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(-0.22458354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8471223) q[0];
sx q[0];
rz(-1.0884388) q[0];
sx q[0];
rz(-0.76089528) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2055254) q[2];
sx q[2];
rz(-2.1238616) q[2];
sx q[2];
rz(1.6286563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7763514) q[1];
sx q[1];
rz(-1.6851211) q[1];
sx q[1];
rz(-2.6120595) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8524695) q[3];
sx q[3];
rz(-1.895088) q[3];
sx q[3];
rz(1.3642023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1647722) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(-2.9898047) q[2];
rz(-0.031489059) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(-0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0781773) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(-2.0582485) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(-1.5189717) q[2];
sx q[2];
rz(-1.3695649) q[2];
sx q[2];
rz(2.8893378) q[2];
rz(-2.8121171) q[3];
sx q[3];
rz(-1.6678882) q[3];
sx q[3];
rz(-1.9141237) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
