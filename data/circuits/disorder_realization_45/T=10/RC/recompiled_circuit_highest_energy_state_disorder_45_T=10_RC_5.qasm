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
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(-2.411627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8529634) q[0];
sx q[0];
rz(-1.4242111) q[0];
sx q[0];
rz(-2.1781871) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59487409) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(1.7930195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3986721) q[1];
sx q[1];
rz(-1.6026596) q[1];
sx q[1];
rz(-0.58018654) q[1];
x q[2];
rz(-2.9128051) q[3];
sx q[3];
rz(-1.8513894) q[3];
sx q[3];
rz(-2.5737263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1175179) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(2.3226341) q[2];
rz(0.0062746127) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(2.1591469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(2.5421802) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(-2.3960466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065577995) q[0];
sx q[0];
rz(-1.6981372) q[0];
sx q[0];
rz(2.8498184) q[0];
rz(-pi) q[1];
rz(2.3815161) q[2];
sx q[2];
rz(-1.3323297) q[2];
sx q[2];
rz(-2.5395405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88953253) q[1];
sx q[1];
rz(-0.65450689) q[1];
sx q[1];
rz(2.8783074) q[1];
rz(-pi) q[2];
rz(-0.98394139) q[3];
sx q[3];
rz(-2.3592257) q[3];
sx q[3];
rz(2.7593105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-2.8440639) q[2];
sx q[2];
rz(1.7342742) q[2];
rz(-2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5832962) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(-2.1287647) q[0];
rz(0.60802513) q[1];
sx q[1];
rz(-1.5589747) q[1];
sx q[1];
rz(1.2695405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69508775) q[0];
sx q[0];
rz(-1.9534612) q[0];
sx q[0];
rz(1.0951359) q[0];
rz(-3.059194) q[2];
sx q[2];
rz(-0.90148704) q[2];
sx q[2];
rz(1.5853887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.041965466) q[1];
sx q[1];
rz(-1.3746975) q[1];
sx q[1];
rz(2.4134497) q[1];
x q[2];
rz(-0.49130398) q[3];
sx q[3];
rz(-1.2761371) q[3];
sx q[3];
rz(-0.06959411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.515392) q[2];
sx q[2];
rz(-2.0833368) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(3.1332704) q[3];
sx q[3];
rz(-2.1885927) q[3];
sx q[3];
rz(-2.9279809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-3.0030647) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(-1.7648765) q[0];
rz(-1.5273013) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(2.6208904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390065) q[0];
sx q[0];
rz(-1.163536) q[0];
sx q[0];
rz(-2.7558221) q[0];
x q[1];
rz(2.7815656) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(1.6619267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59078539) q[1];
sx q[1];
rz(-2.0001162) q[1];
sx q[1];
rz(-2.5270259) q[1];
x q[2];
rz(-2.2305829) q[3];
sx q[3];
rz(-1.1413594) q[3];
sx q[3];
rz(-2.8991606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(-1.3516124) q[2];
rz(-1.0194408) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(-1.9701689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.4627946) q[0];
sx q[0];
rz(-3.0426262) q[0];
rz(2.4348266) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(-2.6720572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5851368) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(1.742247) q[0];
rz(-pi) q[1];
rz(-2.1060313) q[2];
sx q[2];
rz(-1.240584) q[2];
sx q[2];
rz(-2.2170628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3217056) q[1];
sx q[1];
rz(-0.72607909) q[1];
sx q[1];
rz(-0.5989845) q[1];
x q[2];
rz(-1.7170226) q[3];
sx q[3];
rz(-2.1220045) q[3];
sx q[3];
rz(2.220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2609451) q[2];
sx q[2];
rz(-1.8069043) q[2];
sx q[2];
rz(-1.8355231) q[2];
rz(-2.7169054) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223618) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(-1.3223883) q[0];
rz(2.0139096) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(-1.3816396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19162543) q[0];
sx q[0];
rz(-0.81474308) q[0];
sx q[0];
rz(0.11884584) q[0];
rz(-pi) q[1];
rz(-2.1476949) q[2];
sx q[2];
rz(-1.111278) q[2];
sx q[2];
rz(-1.7651059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7378714) q[1];
sx q[1];
rz(-2.3624237) q[1];
sx q[1];
rz(2.8761151) q[1];
x q[2];
rz(2.9074677) q[3];
sx q[3];
rz(-0.30933274) q[3];
sx q[3];
rz(-0.36946378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3809001) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(-2.6603928) q[2];
rz(0.5091269) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(2.7595162) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(1.0786169) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(-0.99348974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29879323) q[0];
sx q[0];
rz(-1.6090819) q[0];
sx q[0];
rz(2.3596615) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69665945) q[2];
sx q[2];
rz(-0.91731812) q[2];
sx q[2];
rz(-0.6763263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76748874) q[1];
sx q[1];
rz(-1.2084157) q[1];
sx q[1];
rz(-0.15717536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.089608) q[3];
sx q[3];
rz(-1.0053653) q[3];
sx q[3];
rz(0.65549248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3420458) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(2.0461931) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(-2.5636165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(1.3336257) q[0];
rz(0.85583055) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(0.91748253) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30675754) q[0];
sx q[0];
rz(-1.4687755) q[0];
sx q[0];
rz(-1.3707042) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4063666) q[2];
sx q[2];
rz(-2.5114369) q[2];
sx q[2];
rz(-0.024261628) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92261295) q[1];
sx q[1];
rz(-2.726788) q[1];
sx q[1];
rz(-0.99402438) q[1];
rz(-pi) q[2];
rz(-0.6360953) q[3];
sx q[3];
rz(-0.75371742) q[3];
sx q[3];
rz(-1.6192644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5069919) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(-1.8514006) q[2];
rz(-2.2318132) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23583394) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(-1.9174891) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.3714429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200206) q[0];
sx q[0];
rz(-0.69584457) q[0];
sx q[0];
rz(-0.30465845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65176378) q[2];
sx q[2];
rz(-0.65417505) q[2];
sx q[2];
rz(2.4602565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4491357) q[1];
sx q[1];
rz(-2.4583011) q[1];
sx q[1];
rz(1.9016674) q[1];
rz(-1.7684494) q[3];
sx q[3];
rz(-2.4098793) q[3];
sx q[3];
rz(-0.55303516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.203043) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7670583) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(-0.57149291) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(2.4933955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4179736) q[0];
sx q[0];
rz(-0.026935808) q[0];
sx q[0];
rz(2.1108146) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66301753) q[2];
sx q[2];
rz(-0.94919357) q[2];
sx q[2];
rz(0.66689516) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1163016) q[1];
sx q[1];
rz(-0.87193438) q[1];
sx q[1];
rz(-1.5998597) q[1];
rz(-0.33721029) q[3];
sx q[3];
rz(-0.54223947) q[3];
sx q[3];
rz(2.6553939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28006831) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(-1.466922) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27547729) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(-2.7571309) q[1];
sx q[1];
rz(-1.5059595) q[1];
sx q[1];
rz(-0.62677871) q[1];
rz(-2.7914417) q[2];
sx q[2];
rz(-1.2731009) q[2];
sx q[2];
rz(0.44558744) q[2];
rz(-0.81955688) q[3];
sx q[3];
rz(-0.71376505) q[3];
sx q[3];
rz(-2.5160088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
