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
rz(2.7449961) q[1];
sx q[1];
rz(-0.34062579) q[1];
sx q[1];
rz(-0.61000282) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1413795) q[0];
sx q[0];
rz(-1.8166409) q[0];
sx q[0];
rz(2.8971998) q[0];
x q[1];
rz(2.0135569) q[2];
sx q[2];
rz(-1.4024078) q[2];
sx q[2];
rz(1.6131282) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0428517) q[1];
sx q[1];
rz(-1.0611649) q[1];
sx q[1];
rz(2.3364728) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8641256) q[3];
sx q[3];
rz(-2.0405117) q[3];
sx q[3];
rz(-1.0535002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4492599) q[2];
sx q[2];
rz(-1.2168987) q[2];
sx q[2];
rz(-0.24162351) q[2];
rz(3.0786476) q[3];
sx q[3];
rz(-1.9398305) q[3];
sx q[3];
rz(0.81301779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6908506) q[0];
sx q[0];
rz(-1.2797322) q[0];
sx q[0];
rz(-0.29020852) q[0];
rz(-2.8065575) q[1];
sx q[1];
rz(-0.93828833) q[1];
sx q[1];
rz(-0.32523528) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9262518) q[0];
sx q[0];
rz(-1.5299954) q[0];
sx q[0];
rz(2.6759139) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89128691) q[2];
sx q[2];
rz(-1.3625979) q[2];
sx q[2];
rz(0.03574275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5618111) q[1];
sx q[1];
rz(-2.9159286) q[1];
sx q[1];
rz(2.2975986) q[1];
rz(-2.4392862) q[3];
sx q[3];
rz(-0.78506535) q[3];
sx q[3];
rz(-2.5007574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8364496) q[2];
sx q[2];
rz(-0.93905753) q[2];
sx q[2];
rz(-1.1391501) q[2];
rz(0.75974733) q[3];
sx q[3];
rz(-3.0436438) q[3];
sx q[3];
rz(-0.37794149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29001319) q[0];
sx q[0];
rz(-0.72999287) q[0];
sx q[0];
rz(0.83339018) q[0];
rz(-0.003412811) q[1];
sx q[1];
rz(-2.623787) q[1];
sx q[1];
rz(1.1606257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0689499) q[0];
sx q[0];
rz(-2.8153689) q[0];
sx q[0];
rz(-1.6618927) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7951444) q[2];
sx q[2];
rz(-0.93502155) q[2];
sx q[2];
rz(1.1368542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6033481) q[1];
sx q[1];
rz(-1.5704435) q[1];
sx q[1];
rz(-1.3360436) q[1];
rz(-pi) q[2];
rz(-1.1228325) q[3];
sx q[3];
rz(-2.7395757) q[3];
sx q[3];
rz(-1.0075008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87355906) q[2];
sx q[2];
rz(-1.4159091) q[2];
sx q[2];
rz(-0.1669008) q[2];
rz(1.1901101) q[3];
sx q[3];
rz(-2.8968865) q[3];
sx q[3];
rz(1.083583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13389182) q[0];
sx q[0];
rz(-1.5465443) q[0];
sx q[0];
rz(-2.290945) q[0];
rz(0.8029241) q[1];
sx q[1];
rz(-1.9839169) q[1];
sx q[1];
rz(2.0096774) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7106774) q[0];
sx q[0];
rz(-1.7516023) q[0];
sx q[0];
rz(-2.0742832) q[0];
x q[1];
rz(-2.1430127) q[2];
sx q[2];
rz(-2.7957186) q[2];
sx q[2];
rz(1.5865751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4635361) q[1];
sx q[1];
rz(-0.66036036) q[1];
sx q[1];
rz(-0.44594958) q[1];
x q[2];
rz(1.1617125) q[3];
sx q[3];
rz(-2.2581824) q[3];
sx q[3];
rz(1.1339172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.930437) q[0];
sx q[0];
rz(-1.3933975) q[0];
sx q[0];
rz(-0.18779553) q[0];
rz(2.4941749) q[1];
sx q[1];
rz(-0.96962601) q[1];
sx q[1];
rz(2.5513249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3844073) q[0];
sx q[0];
rz(-1.2110855) q[0];
sx q[0];
rz(-0.033974302) q[0];
x q[1];
rz(-2.3529604) q[2];
sx q[2];
rz(-1.3071539) q[2];
sx q[2];
rz(0.90197998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62377159) q[1];
sx q[1];
rz(-1.54269) q[1];
sx q[1];
rz(-3.0610316) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8742358) q[3];
sx q[3];
rz(-2.4594569) q[3];
sx q[3];
rz(2.8510044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65842998) q[2];
sx q[2];
rz(-1.1700609) q[2];
sx q[2];
rz(0.17407334) q[2];
rz(-2.6744794) q[3];
sx q[3];
rz(-2.4090448) q[3];
sx q[3];
rz(-2.9113801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1751404) q[0];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0421222) q[0];
sx q[0];
rz(-2.5613031) q[0];
sx q[0];
rz(0.97162928) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0211391) q[2];
sx q[2];
rz(-1.9423123) q[2];
sx q[2];
rz(2.5132883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5111635) q[1];
sx q[1];
rz(-1.9952718) q[1];
sx q[1];
rz(1.9349348) q[1];
rz(-pi) q[2];
rz(0.83281886) q[3];
sx q[3];
rz(-0.76945451) q[3];
sx q[3];
rz(-0.11386816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33438385) q[2];
sx q[2];
rz(-1.7024567) q[2];
sx q[2];
rz(1.8857694) q[2];
rz(0.0028336023) q[3];
sx q[3];
rz(-2.5815559) q[3];
sx q[3];
rz(-0.66633666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46347076) q[0];
sx q[0];
rz(-0.016962873) q[0];
sx q[0];
rz(-2.8609138) q[0];
rz(-0.30411389) q[1];
sx q[1];
rz(-0.76144832) q[1];
sx q[1];
rz(-1.6573409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69347914) q[0];
sx q[0];
rz(-1.5904613) q[0];
sx q[0];
rz(3.124191) q[0];
x q[1];
rz(-1.4062469) q[2];
sx q[2];
rz(-0.79716792) q[2];
sx q[2];
rz(0.88594243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4485943) q[1];
sx q[1];
rz(-2.0830375) q[1];
sx q[1];
rz(-2.1650141) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2400886) q[3];
sx q[3];
rz(-0.39319143) q[3];
sx q[3];
rz(2.0313944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41875276) q[2];
sx q[2];
rz(-2.3374228) q[2];
sx q[2];
rz(-2.1201102) q[2];
rz(2.569765) q[3];
sx q[3];
rz(-0.56838667) q[3];
sx q[3];
rz(0.9182601) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21917139) q[0];
sx q[0];
rz(-0.89184856) q[0];
sx q[0];
rz(0.047155596) q[0];
rz(-1.4157408) q[1];
sx q[1];
rz(-1.163131) q[1];
sx q[1];
rz(0.39173752) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82041086) q[0];
sx q[0];
rz(-1.5279378) q[0];
sx q[0];
rz(-0.45351299) q[0];
x q[1];
rz(-2.0676631) q[2];
sx q[2];
rz(-0.38172445) q[2];
sx q[2];
rz(0.055502467) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89501666) q[1];
sx q[1];
rz(-1.5840753) q[1];
sx q[1];
rz(-0.010427019) q[1];
rz(-pi) q[2];
rz(-2.7487031) q[3];
sx q[3];
rz(-1.3706012) q[3];
sx q[3];
rz(-1.4999215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2262816) q[2];
sx q[2];
rz(-2.2690124) q[2];
sx q[2];
rz(0.20381168) q[2];
rz(-0.63621825) q[3];
sx q[3];
rz(-1.7812984) q[3];
sx q[3];
rz(-0.067559592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5485452) q[0];
sx q[0];
rz(-0.021634463) q[0];
sx q[0];
rz(-1.1170603) q[0];
rz(-1.2720269) q[1];
sx q[1];
rz(-0.8465603) q[1];
sx q[1];
rz(2.5844432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37763835) q[0];
sx q[0];
rz(-1.391811) q[0];
sx q[0];
rz(-2.0077487) q[0];
x q[1];
rz(2.0821758) q[2];
sx q[2];
rz(-0.18958986) q[2];
sx q[2];
rz(2.1185045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3342569) q[1];
sx q[1];
rz(-0.42667056) q[1];
sx q[1];
rz(-2.2183499) q[1];
rz(-pi) q[2];
rz(0.45175868) q[3];
sx q[3];
rz(-0.58940277) q[3];
sx q[3];
rz(1.2337934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7316526) q[2];
sx q[2];
rz(-2.9059548) q[2];
sx q[2];
rz(-1.635599) q[2];
rz(0.19009863) q[3];
sx q[3];
rz(-0.59571576) q[3];
sx q[3];
rz(0.070040919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53368038) q[0];
sx q[0];
rz(-0.30617014) q[0];
sx q[0];
rz(-2.876907) q[0];
rz(-0.39868042) q[1];
sx q[1];
rz(-0.52972263) q[1];
sx q[1];
rz(-0.20464373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7913403) q[0];
sx q[0];
rz(-0.58316427) q[0];
sx q[0];
rz(1.9406609) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25252931) q[2];
sx q[2];
rz(-2.1871532) q[2];
sx q[2];
rz(-1.4184619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73749048) q[1];
sx q[1];
rz(-2.3125517) q[1];
sx q[1];
rz(-0.082423969) q[1];
rz(-pi) q[2];
rz(0.47765215) q[3];
sx q[3];
rz(-1.4414476) q[3];
sx q[3];
rz(-2.1952598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9201422) q[2];
sx q[2];
rz(-1.7998989) q[2];
sx q[2];
rz(1.7500925) q[2];
rz(2.5642388) q[3];
sx q[3];
rz(-2.5632863) q[3];
sx q[3];
rz(-0.81380832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53265424) q[0];
sx q[0];
rz(-1.895099) q[0];
sx q[0];
rz(1.1005703) q[0];
rz(2.7990923) q[1];
sx q[1];
rz(-2.1771912) q[1];
sx q[1];
rz(2.049581) q[1];
rz(2.9663646) q[2];
sx q[2];
rz(-1.874085) q[2];
sx q[2];
rz(0.25599538) q[2];
rz(0.84273356) q[3];
sx q[3];
rz(-2.032654) q[3];
sx q[3];
rz(0.35498735) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
