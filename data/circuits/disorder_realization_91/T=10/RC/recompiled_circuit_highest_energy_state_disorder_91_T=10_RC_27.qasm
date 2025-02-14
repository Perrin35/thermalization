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
rz(0.089601547) q[0];
sx q[0];
rz(2.2360585) q[0];
sx q[0];
rz(9.6146248) q[0];
rz(-0.39659652) q[1];
sx q[1];
rz(3.4822184) q[1];
sx q[1];
rz(10.034781) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34380117) q[0];
sx q[0];
rz(-2.7966948) q[0];
sx q[0];
rz(0.80356046) q[0];
rz(-2.0135569) q[2];
sx q[2];
rz(-1.4024078) q[2];
sx q[2];
rz(1.5284644) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.09874092) q[1];
sx q[1];
rz(-2.0804278) q[1];
sx q[1];
rz(-0.80511989) q[1];
rz(-pi) q[2];
rz(1.0759495) q[3];
sx q[3];
rz(-2.6013654) q[3];
sx q[3];
rz(-1.5264508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4492599) q[2];
sx q[2];
rz(-1.2168987) q[2];
sx q[2];
rz(0.24162351) q[2];
rz(-3.0786476) q[3];
sx q[3];
rz(-1.2017622) q[3];
sx q[3];
rz(-2.3285749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6908506) q[0];
sx q[0];
rz(-1.2797322) q[0];
sx q[0];
rz(-2.8513841) q[0];
rz(-2.8065575) q[1];
sx q[1];
rz(-0.93828833) q[1];
sx q[1];
rz(2.8163574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37595108) q[0];
sx q[0];
rz(-1.1055357) q[0];
sx q[0];
rz(-1.6164533) q[0];
x q[1];
rz(0.26518719) q[2];
sx q[2];
rz(-2.2329805) q[2];
sx q[2];
rz(-1.7005077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27694459) q[1];
sx q[1];
rz(-1.7200302) q[1];
sx q[1];
rz(-1.4008888) q[1];
rz(-pi) q[2];
rz(0.65167221) q[3];
sx q[3];
rz(-2.044994) q[3];
sx q[3];
rz(1.4693361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8364496) q[2];
sx q[2];
rz(-2.2025351) q[2];
sx q[2];
rz(-1.1391501) q[2];
rz(-0.75974733) q[3];
sx q[3];
rz(-3.0436438) q[3];
sx q[3];
rz(-2.7636512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515795) q[0];
sx q[0];
rz(-0.72999287) q[0];
sx q[0];
rz(2.3082025) q[0];
rz(-0.003412811) q[1];
sx q[1];
rz(-0.5178057) q[1];
sx q[1];
rz(-1.1606257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9728119) q[0];
sx q[0];
rz(-1.245975) q[0];
sx q[0];
rz(-0.030766597) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7951444) q[2];
sx q[2];
rz(-0.93502155) q[2];
sx q[2];
rz(-2.0047385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6033481) q[1];
sx q[1];
rz(-1.5711492) q[1];
sx q[1];
rz(-1.3360436) q[1];
rz(-pi) q[2];
rz(2.9594775) q[3];
sx q[3];
rz(-1.9312177) q[3];
sx q[3];
rz(-0.52626901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2680336) q[2];
sx q[2];
rz(-1.4159091) q[2];
sx q[2];
rz(-0.1669008) q[2];
rz(1.1901101) q[3];
sx q[3];
rz(-2.8968865) q[3];
sx q[3];
rz(-2.0580097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0077008) q[0];
sx q[0];
rz(-1.5950483) q[0];
sx q[0];
rz(0.85064763) q[0];
rz(-2.3386686) q[1];
sx q[1];
rz(-1.1576757) q[1];
sx q[1];
rz(1.1319152) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9029805) q[0];
sx q[0];
rz(-1.0762666) q[0];
sx q[0];
rz(0.20574768) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1430127) q[2];
sx q[2];
rz(-2.7957186) q[2];
sx q[2];
rz(-1.5550176) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4635361) q[1];
sx q[1];
rz(-2.4812323) q[1];
sx q[1];
rz(2.6956431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72992562) q[3];
sx q[3];
rz(-1.2582964) q[3];
sx q[3];
rz(2.4362628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6159281) q[2];
sx q[2];
rz(-1.5027081) q[2];
sx q[2];
rz(0.28044236) q[2];
rz(2.7403455) q[3];
sx q[3];
rz(-2.8420227) q[3];
sx q[3];
rz(1.7846599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2111557) q[0];
sx q[0];
rz(-1.3933975) q[0];
sx q[0];
rz(-0.18779553) q[0];
rz(-0.64741778) q[1];
sx q[1];
rz(-0.96962601) q[1];
sx q[1];
rz(2.5513249) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85344106) q[0];
sx q[0];
rz(-0.36124215) q[0];
sx q[0];
rz(-1.4807184) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3529604) q[2];
sx q[2];
rz(-1.8344387) q[2];
sx q[2];
rz(2.2396127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62377159) q[1];
sx q[1];
rz(-1.54269) q[1];
sx q[1];
rz(-0.080561056) q[1];
x q[2];
rz(-2.8742358) q[3];
sx q[3];
rz(-2.4594569) q[3];
sx q[3];
rz(2.8510044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4831627) q[2];
sx q[2];
rz(-1.9715318) q[2];
sx q[2];
rz(0.17407334) q[2];
rz(-2.6744794) q[3];
sx q[3];
rz(-0.73254782) q[3];
sx q[3];
rz(2.9113801) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9664522) q[0];
sx q[0];
rz(-1.3850965) q[0];
sx q[0];
rz(-2.0542282) q[0];
rz(-0.1611791) q[1];
sx q[1];
rz(-1.168707) q[1];
sx q[1];
rz(-2.2081614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3573489) q[0];
sx q[0];
rz(-1.1009365) q[0];
sx q[0];
rz(-0.35413262) q[0];
rz(1.8699588) q[2];
sx q[2];
rz(-0.38969061) q[2];
sx q[2];
rz(0.95012939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78466985) q[1];
sx q[1];
rz(-1.2402727) q[1];
sx q[1];
rz(0.4504942) q[1];
rz(-0.57755135) q[3];
sx q[3];
rz(-2.1114919) q[3];
sx q[3];
rz(-2.3531928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33438385) q[2];
sx q[2];
rz(-1.4391359) q[2];
sx q[2];
rz(1.8857694) q[2];
rz(0.0028336023) q[3];
sx q[3];
rz(-0.56003672) q[3];
sx q[3];
rz(0.66633666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6781219) q[0];
sx q[0];
rz(-0.016962873) q[0];
sx q[0];
rz(2.8609138) q[0];
rz(-0.30411389) q[1];
sx q[1];
rz(-0.76144832) q[1];
sx q[1];
rz(1.4842518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87697498) q[0];
sx q[0];
rz(-1.5533981) q[0];
sx q[0];
rz(-1.5511284) q[0];
rz(-0.16616352) q[2];
sx q[2];
rz(-0.78736178) q[2];
sx q[2];
rz(0.65262567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3918275) q[1];
sx q[1];
rz(-2.3778262) q[1];
sx q[1];
rz(2.3583724) q[1];
rz(1.2400886) q[3];
sx q[3];
rz(-2.7484012) q[3];
sx q[3];
rz(-1.1101983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7228399) q[2];
sx q[2];
rz(-2.3374228) q[2];
sx q[2];
rz(-1.0214825) q[2];
rz(2.569765) q[3];
sx q[3];
rz(-0.56838667) q[3];
sx q[3];
rz(-2.2233326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21917139) q[0];
sx q[0];
rz(-2.2497441) q[0];
sx q[0];
rz(-3.0944371) q[0];
rz(-1.4157408) q[1];
sx q[1];
rz(-1.9784617) q[1];
sx q[1];
rz(-0.39173752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83806709) q[0];
sx q[0];
rz(-2.6861992) q[0];
sx q[0];
rz(-0.097571578) q[0];
x q[1];
rz(1.2315627) q[2];
sx q[2];
rz(-1.7493141) q[2];
sx q[2];
rz(-2.0925131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4656745) q[1];
sx q[1];
rz(-1.5603702) q[1];
sx q[1];
rz(-1.584076) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7487031) q[3];
sx q[3];
rz(-1.3706012) q[3];
sx q[3];
rz(-1.4999215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91531104) q[2];
sx q[2];
rz(-0.87258029) q[2];
sx q[2];
rz(2.937781) q[2];
rz(-0.63621825) q[3];
sx q[3];
rz(-1.7812984) q[3];
sx q[3];
rz(3.0740331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5485452) q[0];
sx q[0];
rz(-0.021634463) q[0];
sx q[0];
rz(2.0245323) q[0];
rz(-1.2720269) q[1];
sx q[1];
rz(-2.2950324) q[1];
sx q[1];
rz(-2.5844432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37763835) q[0];
sx q[0];
rz(-1.391811) q[0];
sx q[0];
rz(2.0077487) q[0];
x q[1];
rz(-0.093634886) q[2];
sx q[2];
rz(-1.7358923) q[2];
sx q[2];
rz(-1.5993725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0275201) q[1];
sx q[1];
rz(-1.2344242) q[1];
sx q[1];
rz(0.26765021) q[1];
x q[2];
rz(2.689834) q[3];
sx q[3];
rz(-2.5521899) q[3];
sx q[3];
rz(1.2337934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40994) q[2];
sx q[2];
rz(-2.9059548) q[2];
sx q[2];
rz(-1.635599) q[2];
rz(-2.951494) q[3];
sx q[3];
rz(-2.5458769) q[3];
sx q[3];
rz(-0.070040919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6079123) q[0];
sx q[0];
rz(-0.30617014) q[0];
sx q[0];
rz(0.26468563) q[0];
rz(-0.39868042) q[1];
sx q[1];
rz(-2.61187) q[1];
sx q[1];
rz(0.20464373) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5335352) q[0];
sx q[0];
rz(-1.7711955) q[0];
sx q[0];
rz(-1.0193558) q[0];
x q[1];
rz(-2.8890633) q[2];
sx q[2];
rz(-2.1871532) q[2];
sx q[2];
rz(-1.7231307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4041022) q[1];
sx q[1];
rz(-0.82904094) q[1];
sx q[1];
rz(-3.0591687) q[1];
x q[2];
rz(0.47765215) q[3];
sx q[3];
rz(-1.700145) q[3];
sx q[3];
rz(-0.94633284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9201422) q[2];
sx q[2];
rz(-1.7998989) q[2];
sx q[2];
rz(1.3915001) q[2];
rz(-2.5642388) q[3];
sx q[3];
rz(-0.57830638) q[3];
sx q[3];
rz(-0.81380832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6089384) q[0];
sx q[0];
rz(-1.2464936) q[0];
sx q[0];
rz(-2.0410224) q[0];
rz(-0.34250034) q[1];
sx q[1];
rz(-2.1771912) q[1];
sx q[1];
rz(2.049581) q[1];
rz(0.17522801) q[2];
sx q[2];
rz(-1.2675076) q[2];
sx q[2];
rz(-2.8855973) q[2];
rz(2.5534775) q[3];
sx q[3];
rz(-0.93265231) q[3];
sx q[3];
rz(1.5476641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
