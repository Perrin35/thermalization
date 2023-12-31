OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(3.9711877) q[0];
sx q[0];
rz(9.2708099) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(-2.8032803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1725537) q[0];
sx q[0];
rz(-1.8882897) q[0];
sx q[0];
rz(2.88455) q[0];
rz(-pi) q[1];
rz(-2.0182274) q[2];
sx q[2];
rz(-2.4123203) q[2];
sx q[2];
rz(2.4016618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6602064) q[1];
sx q[1];
rz(-1.2848789) q[1];
sx q[1];
rz(-1.6260765) q[1];
rz(1.0584352) q[3];
sx q[3];
rz(-1.5844987) q[3];
sx q[3];
rz(-2.0126359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.9677229) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8006111) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(-0.064963438) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.2423135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64105469) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(-0.48135249) q[0];
rz(-pi) q[1];
rz(0.1588891) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(-1.5140669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17774432) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(-2.6049155) q[1];
rz(0.70988016) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80766455) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-0.58369613) q[2];
rz(0.57404533) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(-0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-0.72845355) q[0];
rz(-1.6473673) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(1.0167936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66660488) q[0];
sx q[0];
rz(-3.0525065) q[0];
sx q[0];
rz(-0.41390093) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3967907) q[2];
sx q[2];
rz(-2.475127) q[2];
sx q[2];
rz(-1.1643861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46360717) q[1];
sx q[1];
rz(-1.7079759) q[1];
sx q[1];
rz(-0.82171085) q[1];
x q[2];
rz(-0.049116491) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(-1.0055055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8124354) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(-0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-2.8947815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033337489) q[0];
sx q[0];
rz(-0.70972432) q[0];
sx q[0];
rz(-1.8002585) q[0];
rz(0.47841448) q[2];
sx q[2];
rz(-1.3832428) q[2];
sx q[2];
rz(1.071655) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3226763) q[1];
sx q[1];
rz(-0.37441844) q[1];
sx q[1];
rz(2.9973292) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28624268) q[3];
sx q[3];
rz(-2.1677368) q[3];
sx q[3];
rz(2.3846574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(1.5412615) q[2];
rz(-2.4345496) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(-0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.33070579) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(-1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(-2.8932103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7800956) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(-1.6105152) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70154538) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(-2.0602351) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.76368139) q[1];
sx q[1];
rz(-1.7421725) q[1];
sx q[1];
rz(2.548449) q[1];
rz(-pi) q[2];
rz(2.8304843) q[3];
sx q[3];
rz(-2.4690383) q[3];
sx q[3];
rz(-0.5141408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(-0.38875368) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.7957934) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(3.016901) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42449441) q[0];
sx q[0];
rz(-1.6569123) q[0];
sx q[0];
rz(2.9647102) q[0];
x q[1];
rz(-0.6090392) q[2];
sx q[2];
rz(-1.3281203) q[2];
sx q[2];
rz(-2.2720624) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60285073) q[1];
sx q[1];
rz(-0.83487836) q[1];
sx q[1];
rz(-1.49453) q[1];
rz(2.9051404) q[3];
sx q[3];
rz(-1.7835622) q[3];
sx q[3];
rz(-0.80602431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51320118) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(-2.52264) q[2];
rz(1.0533054) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(-2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(1.0429617) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(1.0707062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8541504) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(0.48405148) q[0];
x q[1];
rz(-2.7942065) q[2];
sx q[2];
rz(-1.4457448) q[2];
sx q[2];
rz(1.320968) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0106525) q[1];
sx q[1];
rz(-0.7796692) q[1];
sx q[1];
rz(2.9995549) q[1];
rz(-0.058810874) q[3];
sx q[3];
rz(-2.8088514) q[3];
sx q[3];
rz(0.22418338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-1.0872844) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(-0.94747296) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763801) q[0];
sx q[0];
rz(-1.0351666) q[0];
sx q[0];
rz(3.0239848) q[0];
rz(2.410789) q[2];
sx q[2];
rz(-0.23554221) q[2];
sx q[2];
rz(-1.3256324) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.265043) q[1];
sx q[1];
rz(-1.3339086) q[1];
sx q[1];
rz(2.151728) q[1];
rz(-2.3341228) q[3];
sx q[3];
rz(-1.315457) q[3];
sx q[3];
rz(1.0367928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(2.6718111) q[2];
rz(1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(1.3289733) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(2.9387617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348939) q[0];
sx q[0];
rz(-1.5447504) q[0];
sx q[0];
rz(-3.1025725) q[0];
rz(-pi) q[1];
rz(2.9183396) q[2];
sx q[2];
rz(-1.7749783) q[2];
sx q[2];
rz(2.5399361) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.054246) q[1];
sx q[1];
rz(-3.0294703) q[1];
sx q[1];
rz(-2.0263158) q[1];
rz(2.8392302) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(-2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7245076) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(2.1450796) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(0.043047992) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(-0.28082401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107287) q[0];
sx q[0];
rz(-2.4952336) q[0];
sx q[0];
rz(-2.3848563) q[0];
rz(1.3781204) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(-3.0378621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3481969) q[1];
sx q[1];
rz(-1.6284202) q[1];
sx q[1];
rz(-3.1053931) q[1];
rz(-pi) q[2];
rz(2.093408) q[3];
sx q[3];
rz(-0.37915737) q[3];
sx q[3];
rz(2.2446333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(0.62310702) q[2];
rz(1.0021707) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-2.5429824) q[2];
sx q[2];
rz(-2.241588) q[2];
sx q[2];
rz(-0.66551756) q[2];
rz(-2.279083) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
