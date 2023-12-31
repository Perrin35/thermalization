OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87969765) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(2.4796955) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0628113) q[2];
sx q[2];
rz(-2.0487818) q[2];
sx q[2];
rz(0.21131549) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2204199) q[1];
sx q[1];
rz(-0.49282679) q[1];
sx q[1];
rz(-0.9763896) q[1];
rz(2.4926315) q[3];
sx q[3];
rz(-1.0421841) q[3];
sx q[3];
rz(1.1543857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0592247) q[0];
sx q[0];
rz(-1.1384283) q[0];
sx q[0];
rz(2.5655377) q[0];
x q[1];
rz(2.3697853) q[2];
sx q[2];
rz(-2.9559921) q[2];
sx q[2];
rz(1.1112569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37296346) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(1.782934) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50176974) q[3];
sx q[3];
rz(-1.0106777) q[3];
sx q[3];
rz(-3.113941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0057384) q[0];
sx q[0];
rz(-1.5691301) q[0];
sx q[0];
rz(1.1093596) q[0];
rz(-pi) q[1];
rz(2.2181182) q[2];
sx q[2];
rz(-1.4787276) q[2];
sx q[2];
rz(1.1084686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4003488) q[1];
sx q[1];
rz(-1.8596008) q[1];
sx q[1];
rz(-1.5653533) q[1];
rz(-2.9680786) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(0.80843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0125668) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(-0.16391779) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.4455459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3460798) q[0];
sx q[0];
rz(-0.68329408) q[0];
sx q[0];
rz(0.28738316) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38979001) q[2];
sx q[2];
rz(-2.9082657) q[2];
sx q[2];
rz(2.0248272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2913937) q[1];
sx q[1];
rz(-1.351601) q[1];
sx q[1];
rz(1.5047969) q[1];
x q[2];
rz(-2.0479855) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(-1.5968061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0754898) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(1.0243105) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(-2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-0.85420001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695278) q[0];
sx q[0];
rz(-2.2846691) q[0];
sx q[0];
rz(0.010849997) q[0];
rz(-pi) q[1];
rz(2.2703083) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(0.7427578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8728719) q[1];
sx q[1];
rz(-2.1495021) q[1];
sx q[1];
rz(2.242356) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79980005) q[3];
sx q[3];
rz(-0.721867) q[3];
sx q[3];
rz(-2.0197899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(-0.50950766) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-3.0153826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55805154) q[0];
sx q[0];
rz(-1.7833033) q[0];
sx q[0];
rz(1.6030747) q[0];
rz(-pi) q[1];
rz(-2.7563165) q[2];
sx q[2];
rz(-0.81309536) q[2];
sx q[2];
rz(-2.3713881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8773552) q[1];
sx q[1];
rz(-2.4009631) q[1];
sx q[1];
rz(1.8514368) q[1];
x q[2];
rz(1.8993127) q[3];
sx q[3];
rz(-1.3887172) q[3];
sx q[3];
rz(1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2386027) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(0.098408498) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0347621) q[0];
sx q[0];
rz(-0.97495279) q[0];
sx q[0];
rz(-1.725127) q[0];
rz(-pi) q[1];
rz(2.9497629) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(-1.9907469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4899788) q[1];
sx q[1];
rz(-1.7001121) q[1];
sx q[1];
rz(-0.08820769) q[1];
rz(0.26569326) q[3];
sx q[3];
rz(-2.4943647) q[3];
sx q[3];
rz(-2.2006187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(-0.87160814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5985142) q[0];
sx q[0];
rz(-0.075965479) q[0];
sx q[0];
rz(0.11673467) q[0];
rz(-1.7224738) q[2];
sx q[2];
rz(-1.1693839) q[2];
sx q[2];
rz(0.36827189) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(-0.70198595) q[1];
rz(-pi) q[2];
rz(-0.1400847) q[3];
sx q[3];
rz(-2.1685371) q[3];
sx q[3];
rz(2.7261581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-1.139572) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7302007) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(1.2325531) q[0];
x q[1];
rz(-0.47713251) q[2];
sx q[2];
rz(-1.7413365) q[2];
sx q[2];
rz(2.0810623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7414654) q[1];
sx q[1];
rz(-1.2418081) q[1];
sx q[1];
rz(0.16727438) q[1];
rz(-pi) q[2];
rz(1.8839621) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(-2.9477011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-2.8881853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2449269) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(3.0299597) q[0];
x q[1];
rz(2.9741653) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(1.6850922) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0903783) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(-1.8263837) q[1];
rz(-0.6110544) q[3];
sx q[3];
rz(-1.9285893) q[3];
sx q[3];
rz(-0.65294453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59166756) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(1.3118369) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(-2.3861804) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
