OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(4.6012576) q[1];
sx q[1];
rz(9.6396946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41479933) q[0];
sx q[0];
rz(-0.3814632) q[0];
sx q[0];
rz(-2.3261855) q[0];
rz(-pi) q[1];
rz(2.6724042) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-2.5772622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54675198) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(-2.2307598) q[1];
rz(-2.4258852) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(-2.3472799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(1.365186) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7433886) q[0];
sx q[0];
rz(-0.87289116) q[0];
sx q[0];
rz(2.8845805) q[0];
x q[1];
rz(1.4977658) q[2];
sx q[2];
rz(-1.4023997) q[2];
sx q[2];
rz(-1.3161236) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.6078728) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(2.968722) q[1];
x q[2];
rz(-2.0131301) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(2.341552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74137694) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(2.7064586) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4774287) q[0];
sx q[0];
rz(-1.4633044) q[0];
sx q[0];
rz(2.7880923) q[0];
x q[1];
rz(1.7180195) q[2];
sx q[2];
rz(-1.3230811) q[2];
sx q[2];
rz(1.9927646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.16776925) q[1];
sx q[1];
rz(-0.5127738) q[1];
sx q[1];
rz(1.2330526) q[1];
x q[2];
rz(2.3782303) q[3];
sx q[3];
rz(-1.8244787) q[3];
sx q[3];
rz(-2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(0.91745013) q[0];
rz(0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53752758) q[0];
sx q[0];
rz(-0.91188216) q[0];
sx q[0];
rz(1.5327246) q[0];
rz(-pi) q[1];
rz(-2.8633966) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(2.1759335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3600196) q[1];
sx q[1];
rz(-0.85687602) q[1];
sx q[1];
rz(2.5931231) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78062765) q[3];
sx q[3];
rz(-2.160191) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.778487) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(2.0402133) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-2.662861) q[0];
rz(-2.1084673) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(-2.1889401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081776) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.4334701) q[0];
x q[1];
rz(-2.7258337) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(1.8005467) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6528656) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(-2.8847242) q[1];
x q[2];
rz(-1.6126552) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(-1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(-2.3099242) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.5166327) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(-2.9690202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5694002) q[0];
sx q[0];
rz(-1.5571556) q[0];
sx q[0];
rz(-1.1543247) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8820011) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(-1.358658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1291618) q[1];
sx q[1];
rz(-2.7134502) q[1];
sx q[1];
rz(-1.3151602) q[1];
rz(-0.63038007) q[3];
sx q[3];
rz(-1.3833589) q[3];
sx q[3];
rz(-2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(-2.0149569) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.3379898) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(-2.3849934) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.341757) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(-2.1205649) q[0];
x q[1];
rz(0.43811626) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(-2.9012836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25890299) q[1];
sx q[1];
rz(-1.8001302) q[1];
sx q[1];
rz(2.980466) q[1];
rz(-pi) q[2];
rz(-0.19161253) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(3.074926) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-2.1527122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8152155) q[0];
sx q[0];
rz(-1.2412374) q[0];
sx q[0];
rz(-1.9944847) q[0];
rz(0.093401508) q[2];
sx q[2];
rz(-1.5867234) q[2];
sx q[2];
rz(1.2707368) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7372072) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(1.7350446) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51560651) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(-3.072217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(-2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75893629) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(2.2682155) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72520032) q[0];
sx q[0];
rz(-1.7500449) q[0];
sx q[0];
rz(3.0730625) q[0];
rz(-1.4344425) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(-1.2863976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.020563515) q[1];
sx q[1];
rz(-0.68194807) q[1];
sx q[1];
rz(0.49973439) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7018413) q[3];
sx q[3];
rz(-1.0381178) q[3];
sx q[3];
rz(0.38910481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(2.7673289) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.3964765) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29165927) q[0];
sx q[0];
rz(-2.8072661) q[0];
sx q[0];
rz(0.26419421) q[0];
x q[1];
rz(-0.65812494) q[2];
sx q[2];
rz(-1.3326367) q[2];
sx q[2];
rz(2.4326774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4954088) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(2.3133548) q[1];
x q[2];
rz(-0.26225984) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044534279) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(-0.13327577) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(1.8272022) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(-2.6079569) q[3];
sx q[3];
rz(-1.613637) q[3];
sx q[3];
rz(-2.2231495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
