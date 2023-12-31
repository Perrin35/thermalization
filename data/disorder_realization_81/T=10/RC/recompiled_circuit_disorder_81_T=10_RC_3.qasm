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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8060018) q[0];
sx q[0];
rz(-1.5561034) q[0];
sx q[0];
rz(3.0528085) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0300006) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(0.0026207844) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8494107) q[1];
sx q[1];
rz(-2.0690448) q[1];
sx q[1];
rz(1.317418) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88702918) q[3];
sx q[3];
rz(-1.6449882) q[3];
sx q[3];
rz(-2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-0.84428865) q[2];
rz(0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59250295) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.9553604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79016722) q[0];
sx q[0];
rz(-1.5520436) q[0];
sx q[0];
rz(0.055939527) q[0];
rz(1.6329174) q[2];
sx q[2];
rz(-1.7698235) q[2];
sx q[2];
rz(-1.8387427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1520878) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(1.771404) q[1];
rz(-pi) q[2];
rz(0.22963345) q[3];
sx q[3];
rz(-0.7080871) q[3];
sx q[3];
rz(-1.1585483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0120323) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(2.3348715) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(-2.321373) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0778753) q[0];
sx q[0];
rz(-0.88839196) q[0];
sx q[0];
rz(-0.14494411) q[0];
rz(-0.21913146) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-2.6205274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8169176) q[1];
sx q[1];
rz(-1.1500689) q[1];
sx q[1];
rz(-2.7540728) q[1];
rz(-pi) q[2];
rz(1.1060171) q[3];
sx q[3];
rz(-1.4651863) q[3];
sx q[3];
rz(1.5231903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(0.93079981) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-2.8500407) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(-2.2303175) q[0];
rz(-0.43831929) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.8211676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62990084) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(-3.1303309) q[0];
rz(0.40789149) q[2];
sx q[2];
rz(-2.6558999) q[2];
sx q[2];
rz(0.25024589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.79975407) q[1];
sx q[1];
rz(-2.2463887) q[1];
sx q[1];
rz(-1.231133) q[1];
rz(-0.60453316) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(1.7741007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8300366) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.8409761) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1039935) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(-3.0689737) q[0];
rz(-pi) q[1];
rz(0.18538961) q[2];
sx q[2];
rz(-0.37237793) q[2];
sx q[2];
rz(-2.7951954) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0504426) q[1];
sx q[1];
rz(-1.0346518) q[1];
sx q[1];
rz(-1.9552783) q[1];
x q[2];
rz(-0.26228321) q[3];
sx q[3];
rz(-1.4296921) q[3];
sx q[3];
rz(1.2348246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0098003) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(2.664393) q[2];
rz(-0.19208433) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3451097) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5531497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1588622) q[0];
sx q[0];
rz(-1.2356865) q[0];
sx q[0];
rz(-1.9431252) q[0];
x q[1];
rz(0.021868869) q[2];
sx q[2];
rz(-1.936603) q[2];
sx q[2];
rz(-1.0794229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3192056) q[1];
sx q[1];
rz(-0.86383312) q[1];
sx q[1];
rz(1.6824526) q[1];
x q[2];
rz(0.074514975) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(-0.45267347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(1.8590415) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5605374) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-0.67725956) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58986321) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(0.59318869) q[0];
rz(-pi) q[1];
rz(0.13204079) q[2];
sx q[2];
rz(-1.570895) q[2];
sx q[2];
rz(1.4591109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3358826) q[1];
sx q[1];
rz(-1.6478331) q[1];
sx q[1];
rz(-1.9105934) q[1];
rz(-pi) q[2];
rz(1.2491751) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(-1.7082937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7523505) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-2.1288669) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.8922071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218924) q[0];
sx q[0];
rz(-1.3601174) q[0];
sx q[0];
rz(1.6181437) q[0];
x q[1];
rz(2.2970389) q[2];
sx q[2];
rz(-2.484998) q[2];
sx q[2];
rz(3.055228) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82259761) q[1];
sx q[1];
rz(-1.5486451) q[1];
sx q[1];
rz(1.5763361) q[1];
x q[2];
rz(-1.1864248) q[3];
sx q[3];
rz(-2.4723408) q[3];
sx q[3];
rz(-2.3068908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(2.5440149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1406527) q[0];
sx q[0];
rz(-1.6093328) q[0];
sx q[0];
rz(-1.7059822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55677982) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(2.5533887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5850726) q[1];
sx q[1];
rz(-2.6954898) q[1];
sx q[1];
rz(-0.82533605) q[1];
rz(-2.5556373) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-0.023660252) q[0];
rz(-2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-0.67217174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2490847) q[0];
sx q[0];
rz(-1.5810284) q[0];
sx q[0];
rz(-3.1363048) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51151885) q[2];
sx q[2];
rz(-1.5626972) q[2];
sx q[2];
rz(2.9715003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.9085911) q[1];
sx q[1];
rz(-0.71725459) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4168596) q[3];
sx q[3];
rz(-1.9991176) q[3];
sx q[3];
rz(2.0314914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(-1.5036748) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(2.4178986) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(-1.9359246) q[2];
sx q[2];
rz(-0.42386133) q[2];
sx q[2];
rz(-0.78122666) q[2];
rz(1.9772114) q[3];
sx q[3];
rz(-1.8462528) q[3];
sx q[3];
rz(0.13312199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
