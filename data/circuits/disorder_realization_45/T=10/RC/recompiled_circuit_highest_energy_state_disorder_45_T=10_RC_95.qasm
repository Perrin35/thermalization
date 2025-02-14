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
rz(1.9396012) q[0];
sx q[0];
rz(-0.48299462) q[0];
sx q[0];
rz(1.6311837) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(5.723602) q[1];
sx q[1];
rz(8.6948123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6522927) q[0];
sx q[0];
rz(-2.518939) q[0];
sx q[0];
rz(1.3176509) q[0];
rz(-2.5467186) q[2];
sx q[2];
rz(-1.383179) q[2];
sx q[2];
rz(1.7930195) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12355676) q[1];
sx q[1];
rz(-2.5606321) q[1];
sx q[1];
rz(3.0835129) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2830986) q[3];
sx q[3];
rz(-1.7904864) q[3];
sx q[3];
rz(-2.2030597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0240747) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(-0.81895858) q[2];
rz(3.135318) q[3];
sx q[3];
rz(-1.9082021) q[3];
sx q[3];
rz(2.1591469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5903951) q[0];
sx q[0];
rz(0.5994125) q[0];
rz(2.2741611) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(-0.74554602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4670938) q[0];
sx q[0];
rz(-1.8601396) q[0];
sx q[0];
rz(1.7036922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8023984) q[2];
sx q[2];
rz(-0.789398) q[2];
sx q[2];
rz(1.2123002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2170495) q[1];
sx q[1];
rz(-0.94247183) q[1];
sx q[1];
rz(1.373686) q[1];
rz(-pi) q[2];
rz(-2.1576513) q[3];
sx q[3];
rz(-0.7823669) q[3];
sx q[3];
rz(-0.38228211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.45592371) q[2];
sx q[2];
rz(-2.8440639) q[2];
sx q[2];
rz(1.4073184) q[2];
rz(-2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5832962) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(-2.1287647) q[0];
rz(0.60802513) q[1];
sx q[1];
rz(-1.5589747) q[1];
sx q[1];
rz(-1.8720522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758666) q[0];
sx q[0];
rz(-2.009543) q[0];
sx q[0];
rz(0.42515305) q[0];
x q[1];
rz(-1.6744587) q[2];
sx q[2];
rz(-2.4680063) q[2];
sx q[2];
rz(1.6885333) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3568282) q[1];
sx q[1];
rz(-2.2819789) q[1];
sx q[1];
rz(-1.8309092) q[1];
rz(0.57165159) q[3];
sx q[3];
rz(-0.56662512) q[3];
sx q[3];
rz(-1.0039312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030647) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(1.7648765) q[0];
rz(-1.6142913) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(2.6208904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0390065) q[0];
sx q[0];
rz(-1.9780567) q[0];
sx q[0];
rz(-2.7558221) q[0];
rz(-pi) q[1];
rz(-1.0215553) q[2];
sx q[2];
rz(-0.62473544) q[2];
sx q[2];
rz(2.3083589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4473089) q[1];
sx q[1];
rz(-1.0189302) q[1];
sx q[1];
rz(-2.0815316) q[1];
rz(2.2124452) q[3];
sx q[3];
rz(-0.7693253) q[3];
sx q[3];
rz(-1.3206467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6167831) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(-1.3516124) q[2];
rz(1.0194408) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.4627946) q[0];
sx q[0];
rz(-3.0426262) q[0];
rz(-0.70676604) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(2.6720572) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564559) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(-1.742247) q[0];
rz(2.1060313) q[2];
sx q[2];
rz(-1.240584) q[2];
sx q[2];
rz(2.2170628) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5818122) q[1];
sx q[1];
rz(-0.99039927) q[1];
sx q[1];
rz(-1.1066761) q[1];
rz(-0.23271493) q[3];
sx q[3];
rz(-0.56833) q[3];
sx q[3];
rz(-1.1956904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(-2.7169054) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(2.2556321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(1.8192044) q[0];
rz(-2.0139096) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(1.7599531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9499672) q[0];
sx q[0];
rz(-2.3268496) q[0];
sx q[0];
rz(3.0227468) q[0];
x q[1];
rz(2.6082615) q[2];
sx q[2];
rz(-2.0815649) q[2];
sx q[2];
rz(0.086670808) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0387011) q[1];
sx q[1];
rz(-2.3159317) q[1];
sx q[1];
rz(1.3172512) q[1];
rz(-2.9074677) q[3];
sx q[3];
rz(-0.30933274) q[3];
sx q[3];
rz(0.36946378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3809001) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(-0.5091269) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(-2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(-3.0522108) q[0];
rz(-1.0786169) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(-2.1481029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3099965) q[0];
sx q[0];
rz(-0.78959268) q[0];
sx q[0];
rz(-1.5168651) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87331302) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(-0.26604929) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9536994) q[1];
sx q[1];
rz(-0.39361289) q[1];
sx q[1];
rz(1.962349) q[1];
x q[2];
rz(-0.63103038) q[3];
sx q[3];
rz(-2.0027805) q[3];
sx q[3];
rz(-0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3420458) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(0.40204027) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(-1.3336257) q[0];
rz(-2.2857621) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(-0.91748253) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2846904) q[0];
sx q[0];
rz(-1.7698341) q[0];
sx q[0];
rz(-3.0375098) q[0];
rz(2.1945004) q[2];
sx q[2];
rz(-1.4741885) q[2];
sx q[2];
rz(-1.7283224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0304347) q[1];
sx q[1];
rz(-1.7923755) q[1];
sx q[1];
rz(1.9244003) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4948984) q[3];
sx q[3];
rz(-1.1521253) q[3];
sx q[3];
rz(-0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5069919) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(-1.290192) q[2];
rz(2.2318132) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057587) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(-1.9174891) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(1.3714429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200206) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(2.8369342) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54746898) q[2];
sx q[2];
rz(-1.1927422) q[2];
sx q[2];
rz(0.34502703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7595948) q[1];
sx q[1];
rz(-1.7773668) q[1];
sx q[1];
rz(0.91464154) q[1];
rz(-pi) q[2];
rz(0.84884642) q[3];
sx q[3];
rz(-1.7023785) q[3];
sx q[3];
rz(-2.2717486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.203043) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(-2.4533563) q[2];
rz(-0.31442434) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(-1.8800053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7670583) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(0.68069619) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(-2.4933955) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.307319) q[0];
sx q[0];
rz(-1.5569485) q[0];
sx q[0];
rz(-1.5939006) q[0];
rz(0.66301753) q[2];
sx q[2];
rz(-2.1923991) q[2];
sx q[2];
rz(-0.66689516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5773864) q[1];
sx q[1];
rz(-1.5930452) q[1];
sx q[1];
rz(-2.4425227) q[1];
rz(1.7675507) q[3];
sx q[3];
rz(-1.0621539) q[3];
sx q[3];
rz(3.0439049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8615243) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8661154) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(-0.38446174) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(2.7914417) q[2];
sx q[2];
rz(-1.8684917) q[2];
sx q[2];
rz(-2.6960052) q[2];
rz(2.3220358) q[3];
sx q[3];
rz(-0.71376505) q[3];
sx q[3];
rz(-2.5160088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
