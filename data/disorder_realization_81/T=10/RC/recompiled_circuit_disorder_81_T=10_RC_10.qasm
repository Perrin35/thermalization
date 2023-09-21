OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(4.9217304) q[0];
sx q[0];
rz(11.187727) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(2.9773657) q[0];
rz(-pi) q[1];
rz(-0.11159201) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(0.0026207844) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3528459) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(2.7098141) q[1];
rz(2.2545635) q[3];
sx q[3];
rz(-1.6449882) q[3];
sx q[3];
rz(-0.85899734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(-2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3620136) q[0];
sx q[0];
rz(-1.5148666) q[0];
sx q[0];
rz(1.5520142) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19940168) q[2];
sx q[2];
rz(-1.6316895) q[2];
sx q[2];
rz(-0.25564889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9895049) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(1.3701887) q[1];
x q[2];
rz(-1.7632742) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(1.457085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.9821232) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(2.8306567) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(-2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504282) q[0];
sx q[0];
rz(-0.69520742) q[0];
sx q[0];
rz(1.7466963) q[0];
x q[1];
rz(-1.1782896) q[2];
sx q[2];
rz(-2.6143392) q[2];
sx q[2];
rz(2.1737827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0606544) q[1];
sx q[1];
rz(-1.9229691) q[1];
sx q[1];
rz(-1.1206131) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1060171) q[3];
sx q[3];
rz(-1.6764063) q[3];
sx q[3];
rz(1.5231903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8308668) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.8211676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62990084) q[0];
sx q[0];
rz(-1.1612079) q[0];
sx q[0];
rz(3.1303309) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45122066) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(-2.1860683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1530694) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(-2.4371229) q[1];
rz(-2.5370595) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(2.9648182) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-1.140973) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(1.8409761) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9601701) q[0];
sx q[0];
rz(-1.9248065) q[0];
sx q[0];
rz(1.5439073) q[0];
x q[1];
rz(-1.6426716) q[2];
sx q[2];
rz(-1.2051029) q[2];
sx q[2];
rz(0.14771151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0504426) q[1];
sx q[1];
rz(-1.0346518) q[1];
sx q[1];
rz(-1.9552783) q[1];
rz(0.50118581) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(2.3230769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(-0.19208433) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(0.011750301) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5884429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5396744) q[0];
sx q[0];
rz(-1.9214905) q[0];
sx q[0];
rz(2.7838216) q[0];
rz(-1.5137709) q[2];
sx q[2];
rz(-0.36643039) q[2];
sx q[2];
rz(2.1232405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99332422) q[1];
sx q[1];
rz(-2.4273708) q[1];
sx q[1];
rz(-3.0118914) q[1];
rz(-pi) q[2];
rz(-1.6586967) q[3];
sx q[3];
rz(-2.4368736) q[3];
sx q[3];
rz(-0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(1.2825512) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(-2.1645434) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86977406) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(3.0074044) q[0];
x q[1];
rz(1.5706967) q[2];
sx q[2];
rz(-1.4387555) q[2];
sx q[2];
rz(-3.0299203) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9338785) q[1];
sx q[1];
rz(-1.9095451) q[1];
sx q[1];
rz(-0.081688567) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2491751) q[3];
sx q[3];
rz(-1.2444032) q[3];
sx q[3];
rz(1.4332989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.2493856) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69985336) q[0];
sx q[0];
rz(-2.925736) q[0];
sx q[0];
rz(2.9237843) q[0];
rz(2.0935358) q[2];
sx q[2];
rz(-1.153423) q[2];
sx q[2];
rz(1.0440895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.318995) q[1];
sx q[1];
rz(-1.5486451) q[1];
sx q[1];
rz(1.5652565) q[1];
x q[2];
rz(1.9551679) q[3];
sx q[3];
rz(-2.4723408) q[3];
sx q[3];
rz(-2.3068908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.1784941) q[2];
rz(1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4998528) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(2.5440149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15411988) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(-1.2921635) q[0];
x q[1];
rz(2.0963247) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(1.8906821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55652009) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(-0.82533605) q[1];
rz(-0.58595539) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(-2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.2333599) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.6433158) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(-2.1854782) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(-2.4694209) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2490847) q[0];
sx q[0];
rz(-1.5810284) q[0];
sx q[0];
rz(0.0052878629) q[0];
x q[1];
rz(-0.51151885) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(-2.9715003) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8476665) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(2.6508209) q[1];
x q[2];
rz(-0.32398128) q[3];
sx q[3];
rz(-2.6880662) q[3];
sx q[3];
rz(-1.6739664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5621915) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-0.72369408) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(-1.9359246) q[2];
sx q[2];
rz(-0.42386133) q[2];
sx q[2];
rz(-0.78122666) q[2];
rz(2.1914992) q[3];
sx q[3];
rz(-2.6549669) q[3];
sx q[3];
rz(-2.0012729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];