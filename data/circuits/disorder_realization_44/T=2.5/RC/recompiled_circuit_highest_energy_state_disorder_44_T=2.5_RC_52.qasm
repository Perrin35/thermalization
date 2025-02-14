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
rz(2.514578) q[0];
sx q[0];
rz(-0.63627565) q[0];
sx q[0];
rz(-1.0778435) q[0];
rz(-2.0650504) q[1];
sx q[1];
rz(-2.512518) q[1];
sx q[1];
rz(-2.10973) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92562308) q[0];
sx q[0];
rz(-1.5711938) q[0];
sx q[0];
rz(0.0041603869) q[0];
rz(0.59649845) q[2];
sx q[2];
rz(-0.41587999) q[2];
sx q[2];
rz(2.2592827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87776041) q[1];
sx q[1];
rz(-1.227637) q[1];
sx q[1];
rz(-2.7105217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21764619) q[3];
sx q[3];
rz(-2.7983694) q[3];
sx q[3];
rz(-1.1962593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91118139) q[2];
sx q[2];
rz(-1.1831256) q[2];
sx q[2];
rz(-0.46768701) q[2];
rz(-2.6905401) q[3];
sx q[3];
rz(-0.39877287) q[3];
sx q[3];
rz(2.8635645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69538799) q[0];
sx q[0];
rz(-2.9187293) q[0];
sx q[0];
rz(2.9872802) q[0];
rz(-2.800324) q[1];
sx q[1];
rz(-2.7879265) q[1];
sx q[1];
rz(0.42010677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25866461) q[0];
sx q[0];
rz(-0.9416343) q[0];
sx q[0];
rz(1.8101519) q[0];
rz(-pi) q[1];
rz(0.86067274) q[2];
sx q[2];
rz(-1.6063074) q[2];
sx q[2];
rz(-1.4461609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1359033) q[1];
sx q[1];
rz(-1.5694251) q[1];
sx q[1];
rz(0.47079177) q[1];
rz(-pi) q[2];
rz(-0.16815925) q[3];
sx q[3];
rz(-2.3871867) q[3];
sx q[3];
rz(1.8969632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61098617) q[2];
sx q[2];
rz(-2.2425118) q[2];
sx q[2];
rz(2.3685624) q[2];
rz(-0.84345877) q[3];
sx q[3];
rz(-1.5638899) q[3];
sx q[3];
rz(2.5696866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.968349) q[0];
sx q[0];
rz(-2.1143715) q[0];
sx q[0];
rz(0.20208836) q[0];
rz(-0.48187065) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(0.906382) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13440713) q[0];
sx q[0];
rz(-2.2762594) q[0];
sx q[0];
rz(2.5329068) q[0];
rz(2.0252805) q[2];
sx q[2];
rz(-1.8826199) q[2];
sx q[2];
rz(-0.21909595) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8123834) q[1];
sx q[1];
rz(-2.1076729) q[1];
sx q[1];
rz(1.7574134) q[1];
x q[2];
rz(-1.0624978) q[3];
sx q[3];
rz(-2.8064499) q[3];
sx q[3];
rz(0.80860521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0704982) q[2];
sx q[2];
rz(-1.1134032) q[2];
sx q[2];
rz(-2.5615198) q[2];
rz(1.5602559) q[3];
sx q[3];
rz(-1.7507078) q[3];
sx q[3];
rz(1.4238547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(2.4105014) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(0.047792338) q[0];
rz(1.2740678) q[1];
sx q[1];
rz(-1.85227) q[1];
sx q[1];
rz(-3.0963669) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52818283) q[0];
sx q[0];
rz(-0.65887132) q[0];
sx q[0];
rz(1.3263561) q[0];
rz(-2.2061646) q[2];
sx q[2];
rz(-1.3029939) q[2];
sx q[2];
rz(1.8647461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8020683) q[1];
sx q[1];
rz(-1.5652204) q[1];
sx q[1];
rz(0.0030203621) q[1];
rz(-pi) q[2];
rz(-0.43680059) q[3];
sx q[3];
rz(-1.068401) q[3];
sx q[3];
rz(-0.36239195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.73146003) q[2];
sx q[2];
rz(-1.2612017) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(-0.54640031) q[3];
sx q[3];
rz(-1.0011287) q[3];
sx q[3];
rz(0.31118292) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6379717) q[0];
sx q[0];
rz(-1.4542955) q[0];
sx q[0];
rz(-1.7244435) q[0];
rz(-1.1196989) q[1];
sx q[1];
rz(-1.3816625) q[1];
sx q[1];
rz(1.1823357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6851824) q[0];
sx q[0];
rz(-1.6401924) q[0];
sx q[0];
rz(-0.38763898) q[0];
x q[1];
rz(0.77967092) q[2];
sx q[2];
rz(-2.5354249) q[2];
sx q[2];
rz(-0.066628284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97302283) q[1];
sx q[1];
rz(-2.418879) q[1];
sx q[1];
rz(-3.0672468) q[1];
x q[2];
rz(-2.0145217) q[3];
sx q[3];
rz(-2.2699353) q[3];
sx q[3];
rz(-2.3026031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.917439) q[2];
sx q[2];
rz(-2.3771693) q[2];
sx q[2];
rz(0.42382851) q[2];
rz(2.0297) q[3];
sx q[3];
rz(-0.49911505) q[3];
sx q[3];
rz(2.4901938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424585) q[0];
sx q[0];
rz(-0.0722216) q[0];
sx q[0];
rz(2.1891731) q[0];
rz(-0.77596387) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(-1.5752569) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.562378) q[0];
sx q[0];
rz(-1.6088621) q[0];
sx q[0];
rz(-1.5505575) q[0];
rz(-pi) q[1];
rz(-3.1029019) q[2];
sx q[2];
rz(-0.96444791) q[2];
sx q[2];
rz(2.4232466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4500931) q[1];
sx q[1];
rz(-0.97292346) q[1];
sx q[1];
rz(1.3045425) q[1];
rz(-pi) q[2];
rz(1.906988) q[3];
sx q[3];
rz(-2.6264274) q[3];
sx q[3];
rz(-0.060843918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39885193) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(2.3952132) q[2];
rz(-2.0505203) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(-2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2090476) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(0.5433425) q[0];
rz(2.6047193) q[1];
sx q[1];
rz(-0.32506341) q[1];
sx q[1];
rz(-0.12319014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9730102) q[0];
sx q[0];
rz(-1.5401398) q[0];
sx q[0];
rz(2.2066096) q[0];
x q[1];
rz(0.51080475) q[2];
sx q[2];
rz(-1.1996562) q[2];
sx q[2];
rz(-2.5696511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.035288485) q[1];
sx q[1];
rz(-2.114944) q[1];
sx q[1];
rz(0.72334163) q[1];
rz(-pi) q[2];
rz(-2.3402392) q[3];
sx q[3];
rz(-0.94326245) q[3];
sx q[3];
rz(2.827044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3049551) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(-0.48509625) q[2];
rz(1.1525611) q[3];
sx q[3];
rz(-1.7771114) q[3];
sx q[3];
rz(-3.0936354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7570067) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(0.49569976) q[0];
rz(0.063830201) q[1];
sx q[1];
rz(-2.3921236) q[1];
sx q[1];
rz(-3.0705423) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71664229) q[0];
sx q[0];
rz(-2.0293689) q[0];
sx q[0];
rz(1.3032662) q[0];
rz(2.536473) q[2];
sx q[2];
rz(-1.4265043) q[2];
sx q[2];
rz(2.6168106) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9378375) q[1];
sx q[1];
rz(-1.8120288) q[1];
sx q[1];
rz(-0.15934847) q[1];
x q[2];
rz(-2.9661353) q[3];
sx q[3];
rz(-0.4532632) q[3];
sx q[3];
rz(0.82626617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12585982) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(0.86137548) q[2];
rz(-1.3043978) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(-1.5931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9651589) q[0];
sx q[0];
rz(-0.49089828) q[0];
sx q[0];
rz(1.7023671) q[0];
rz(0.08055117) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(1.390994) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8497365) q[0];
sx q[0];
rz(-0.69293368) q[0];
sx q[0];
rz(2.2676208) q[0];
x q[1];
rz(1.2660145) q[2];
sx q[2];
rz(-1.551318) q[2];
sx q[2];
rz(1.2661042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7841107) q[1];
sx q[1];
rz(-0.78178373) q[1];
sx q[1];
rz(-2.2428721) q[1];
rz(-1.9506878) q[3];
sx q[3];
rz(-1.9171439) q[3];
sx q[3];
rz(2.430918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3535658) q[2];
sx q[2];
rz(-1.8466419) q[2];
sx q[2];
rz(0.12727748) q[2];
rz(2.2150529) q[3];
sx q[3];
rz(-2.8997731) q[3];
sx q[3];
rz(-1.658879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38630286) q[0];
sx q[0];
rz(-0.64978623) q[0];
sx q[0];
rz(1.5007098) q[0];
rz(-0.81037784) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(-0.77218974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579073) q[0];
sx q[0];
rz(-2.2026688) q[0];
sx q[0];
rz(2.4993012) q[0];
rz(-pi) q[1];
rz(2.8462662) q[2];
sx q[2];
rz(-1.775303) q[2];
sx q[2];
rz(0.4819862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3898728) q[1];
sx q[1];
rz(-1.0214318) q[1];
sx q[1];
rz(0.013688727) q[1];
rz(-pi) q[2];
rz(-2.2507319) q[3];
sx q[3];
rz(-1.6107585) q[3];
sx q[3];
rz(-3.0639632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7605674) q[2];
sx q[2];
rz(-2.3367391) q[2];
sx q[2];
rz(-2.5123361) q[2];
rz(-0.73090807) q[3];
sx q[3];
rz(-2.2980502) q[3];
sx q[3];
rz(-1.6985016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44169852) q[0];
sx q[0];
rz(-2.6581673) q[0];
sx q[0];
rz(2.7375258) q[0];
rz(-0.40456698) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(0.28180939) q[2];
sx q[2];
rz(-1.6360313) q[2];
sx q[2];
rz(-2.4665063) q[2];
rz(-0.53530219) q[3];
sx q[3];
rz(-1.3131159) q[3];
sx q[3];
rz(2.5406607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
