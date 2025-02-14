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
rz(-1.7173362) q[0];
sx q[0];
rz(-1.9219226) q[0];
sx q[0];
rz(-2.5665459) q[0];
rz(0.14248928) q[1];
sx q[1];
rz(-1.9228851) q[1];
sx q[1];
rz(-0.86427468) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48102114) q[0];
sx q[0];
rz(-2.8459274) q[0];
sx q[0];
rz(-2.2822102) q[0];
rz(-pi) q[1];
rz(-0.15344413) q[2];
sx q[2];
rz(-0.94509238) q[2];
sx q[2];
rz(-2.1476755) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.083347224) q[1];
sx q[1];
rz(-1.5687404) q[1];
sx q[1];
rz(2.6115692) q[1];
x q[2];
rz(-2.5745886) q[3];
sx q[3];
rz(-1.1409014) q[3];
sx q[3];
rz(-0.56741949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.064726) q[2];
sx q[2];
rz(-1.4996424) q[2];
sx q[2];
rz(2.1323252) q[2];
rz(2.3276954) q[3];
sx q[3];
rz(-2.8129357) q[3];
sx q[3];
rz(0.3391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086804248) q[0];
sx q[0];
rz(-1.7330994) q[0];
sx q[0];
rz(-0.24222294) q[0];
rz(1.2308661) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(0.79808527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4674007) q[0];
sx q[0];
rz(-0.62819203) q[0];
sx q[0];
rz(-0.19285658) q[0];
rz(-2.9954916) q[2];
sx q[2];
rz(-1.9105991) q[2];
sx q[2];
rz(2.3229342) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68951571) q[1];
sx q[1];
rz(-0.59714666) q[1];
sx q[1];
rz(-2.4680016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40753461) q[3];
sx q[3];
rz(-0.38713249) q[3];
sx q[3];
rz(-2.3705496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65608394) q[2];
sx q[2];
rz(-2.0723497) q[2];
sx q[2];
rz(2.3089224) q[2];
rz(-0.63140702) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(0.00028636534) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1239419) q[0];
sx q[0];
rz(-0.90921062) q[0];
sx q[0];
rz(-2.9538474) q[0];
rz(0.84367696) q[1];
sx q[1];
rz(-1.8673106) q[1];
sx q[1];
rz(2.5086596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.303514) q[0];
sx q[0];
rz(-1.2084586) q[0];
sx q[0];
rz(-3.0658609) q[0];
x q[1];
rz(1.9203414) q[2];
sx q[2];
rz(-2.2270906) q[2];
sx q[2];
rz(-1.2906769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6784981) q[1];
sx q[1];
rz(-2.5353574) q[1];
sx q[1];
rz(-2.5962968) q[1];
x q[2];
rz(-1.2237596) q[3];
sx q[3];
rz(-1.4489902) q[3];
sx q[3];
rz(-0.95637874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26744276) q[2];
sx q[2];
rz(-2.1953857) q[2];
sx q[2];
rz(-1.2713185) q[2];
rz(1.4960131) q[3];
sx q[3];
rz(-1.4964024) q[3];
sx q[3];
rz(0.97150826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520168) q[0];
sx q[0];
rz(-0.75279623) q[0];
sx q[0];
rz(2.4821607) q[0];
rz(-1.3176428) q[1];
sx q[1];
rz(-1.3392071) q[1];
sx q[1];
rz(2.3514294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1135318) q[0];
sx q[0];
rz(-1.4469997) q[0];
sx q[0];
rz(2.8040299) q[0];
rz(-pi) q[1];
rz(-1.538058) q[2];
sx q[2];
rz(-1.5509836) q[2];
sx q[2];
rz(2.0231501) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15406628) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(2.9996526) q[1];
rz(-pi) q[2];
rz(1.2983367) q[3];
sx q[3];
rz(-2.108223) q[3];
sx q[3];
rz(0.3934653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0822175) q[2];
sx q[2];
rz(-1.5267812) q[2];
sx q[2];
rz(2.8455632) q[2];
rz(-2.5791903) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(3.0276022) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.247308) q[0];
sx q[0];
rz(-1.1881275) q[0];
sx q[0];
rz(2.9344015) q[0];
rz(-1.0317135) q[1];
sx q[1];
rz(-1.9887911) q[1];
sx q[1];
rz(1.656104) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9821859) q[0];
sx q[0];
rz(-0.80538087) q[0];
sx q[0];
rz(2.5872487) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0096512) q[2];
sx q[2];
rz(-1.818294) q[2];
sx q[2];
rz(-2.0405318) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1904162) q[1];
sx q[1];
rz(-0.37280478) q[1];
sx q[1];
rz(-2.784381) q[1];
x q[2];
rz(1.6909825) q[3];
sx q[3];
rz(-0.73420364) q[3];
sx q[3];
rz(2.7199573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.038387211) q[2];
sx q[2];
rz(-0.81192553) q[2];
sx q[2];
rz(-2.0337598) q[2];
rz(-2.2191018) q[3];
sx q[3];
rz(-1.731571) q[3];
sx q[3];
rz(-0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6237727) q[0];
sx q[0];
rz(-2.6277442) q[0];
sx q[0];
rz(2.9129831) q[0];
rz(-0.27433431) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(-0.45023227) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73295702) q[0];
sx q[0];
rz(-1.8224323) q[0];
sx q[0];
rz(-2.3711414) q[0];
x q[1];
rz(3.1338056) q[2];
sx q[2];
rz(-1.8000364) q[2];
sx q[2];
rz(-2.0321329) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5094368) q[1];
sx q[1];
rz(-2.5773199) q[1];
sx q[1];
rz(0.48308259) q[1];
x q[2];
rz(-1.9852909) q[3];
sx q[3];
rz(-1.3236146) q[3];
sx q[3];
rz(2.7657764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.36394128) q[2];
sx q[2];
rz(-0.55299091) q[2];
sx q[2];
rz(-3.1138163) q[2];
rz(-2.3751496) q[3];
sx q[3];
rz(-1.1602297) q[3];
sx q[3];
rz(-0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9195093) q[0];
sx q[0];
rz(-1.5064025) q[0];
sx q[0];
rz(2.5191504) q[0];
rz(0.70603236) q[1];
sx q[1];
rz(-0.89203867) q[1];
sx q[1];
rz(-2.5546254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0737105) q[0];
sx q[0];
rz(-2.1176) q[0];
sx q[0];
rz(0.51954649) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3038425) q[2];
sx q[2];
rz(-1.0994229) q[2];
sx q[2];
rz(-1.4990569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0388958) q[1];
sx q[1];
rz(-1.6417112) q[1];
sx q[1];
rz(-0.50693477) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8587647) q[3];
sx q[3];
rz(-1.5786577) q[3];
sx q[3];
rz(-0.9031537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.065757699) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(-0.59448057) q[2];
rz(0.25932702) q[3];
sx q[3];
rz(-1.2649053) q[3];
sx q[3];
rz(-1.9243141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.038430564) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(2.5627947) q[0];
rz(1.831306) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(2.328918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8978764) q[0];
sx q[0];
rz(-1.0930632) q[0];
sx q[0];
rz(-1.7327139) q[0];
rz(-pi) q[1];
rz(1.9695639) q[2];
sx q[2];
rz(-2.5426455) q[2];
sx q[2];
rz(0.87397611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26749565) q[1];
sx q[1];
rz(-0.90522841) q[1];
sx q[1];
rz(2.0366621) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50847362) q[3];
sx q[3];
rz(-2.0805854) q[3];
sx q[3];
rz(-2.2126861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7909214) q[2];
sx q[2];
rz(-1.2988043) q[2];
sx q[2];
rz(-0.92362967) q[2];
rz(2.1212497) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(3.0237107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0807121) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(1.6545779) q[0];
rz(-0.05038536) q[1];
sx q[1];
rz(-2.3245508) q[1];
sx q[1];
rz(0.64687669) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23316194) q[0];
sx q[0];
rz(-1.9958423) q[0];
sx q[0];
rz(-2.2320497) q[0];
rz(-pi) q[1];
rz(2.0503309) q[2];
sx q[2];
rz(-0.96194907) q[2];
sx q[2];
rz(-0.11907585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3072435) q[1];
sx q[1];
rz(-1.780088) q[1];
sx q[1];
rz(2.7268975) q[1];
x q[2];
rz(-0.96109747) q[3];
sx q[3];
rz(-2.8831579) q[3];
sx q[3];
rz(2.1375382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50463027) q[2];
sx q[2];
rz(-1.5865822) q[2];
sx q[2];
rz(0.75616765) q[2];
rz(1.6715096) q[3];
sx q[3];
rz(-2.3556637) q[3];
sx q[3];
rz(-2.3362931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1038372) q[0];
sx q[0];
rz(-1.8917731) q[0];
sx q[0];
rz(1.1081498) q[0];
rz(-0.26421079) q[1];
sx q[1];
rz(-1.6725531) q[1];
sx q[1];
rz(0.60233751) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7450934) q[0];
sx q[0];
rz(-1.709793) q[0];
sx q[0];
rz(-1.2953912) q[0];
rz(0.71508821) q[2];
sx q[2];
rz(-2.9093766) q[2];
sx q[2];
rz(0.84399904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4385637) q[1];
sx q[1];
rz(-1.6191585) q[1];
sx q[1];
rz(2.2579184) q[1];
rz(-pi) q[2];
rz(-0.29903166) q[3];
sx q[3];
rz(-2.4686738) q[3];
sx q[3];
rz(-2.3367019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63310528) q[2];
sx q[2];
rz(-0.54163951) q[2];
sx q[2];
rz(-0.87149054) q[2];
rz(-1.2095215) q[3];
sx q[3];
rz(-2.8612374) q[3];
sx q[3];
rz(-2.5476294) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7021983) q[0];
sx q[0];
rz(-1.7560503) q[0];
sx q[0];
rz(-0.51881292) q[0];
rz(0.58204542) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(1.9459361) q[2];
sx q[2];
rz(-2.2715501) q[2];
sx q[2];
rz(3.0987433) q[2];
rz(-0.68861674) q[3];
sx q[3];
rz(-1.7027161) q[3];
sx q[3];
rz(-0.30218132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
