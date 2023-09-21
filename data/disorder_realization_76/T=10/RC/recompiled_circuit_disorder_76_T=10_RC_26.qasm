OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(3.0091118) q[0];
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(0.69256988) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4953295) q[0];
sx q[0];
rz(-1.4709934) q[0];
sx q[0];
rz(-2.0908337) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5508467) q[2];
sx q[2];
rz(-2.4060537) q[2];
sx q[2];
rz(2.9172446) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2107271) q[1];
sx q[1];
rz(-0.14666808) q[1];
sx q[1];
rz(-2.0861097) q[1];
rz(-pi) q[2];
rz(-2.8928738) q[3];
sx q[3];
rz(-2.0318444) q[3];
sx q[3];
rz(0.15795262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.6050603) q[2];
rz(-1.6202554) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3261616) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(1.2492299) q[0];
rz(2.5800887) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(-0.5805648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5843129) q[0];
sx q[0];
rz(-1.4627539) q[0];
sx q[0];
rz(0.31038196) q[0];
rz(-pi) q[1];
rz(2.850769) q[2];
sx q[2];
rz(-0.65867701) q[2];
sx q[2];
rz(-0.80292279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21743449) q[1];
sx q[1];
rz(-2.1719451) q[1];
sx q[1];
rz(-2.2648328) q[1];
rz(-pi) q[2];
rz(2.21653) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(-0.95228449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(2.2632329) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(2.3679249) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(-0.032827854) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531909) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(1.3283967) q[0];
rz(1.1187395) q[2];
sx q[2];
rz(-0.60048044) q[2];
sx q[2];
rz(-2.877419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5579538) q[1];
sx q[1];
rz(-1.4065521) q[1];
sx q[1];
rz(2.3418531) q[1];
x q[2];
rz(2.1081984) q[3];
sx q[3];
rz(-2.5279547) q[3];
sx q[3];
rz(0.19526853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(-0.39595655) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4375777) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(0.2335877) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(-2.3707726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17731006) q[0];
sx q[0];
rz(-1.5872123) q[0];
sx q[0];
rz(1.55127) q[0];
x q[1];
rz(-1.0806482) q[2];
sx q[2];
rz(-1.4785826) q[2];
sx q[2];
rz(-2.8871418) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96344906) q[1];
sx q[1];
rz(-1.2091067) q[1];
sx q[1];
rz(2.5835035) q[1];
rz(-pi) q[2];
rz(-0.1905258) q[3];
sx q[3];
rz(-1.3503721) q[3];
sx q[3];
rz(-2.4998553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(0.7129933) q[2];
rz(-2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(0.95389429) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.7011401) q[0];
rz(-3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(0.17000155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3225587) q[0];
sx q[0];
rz(-1.8043307) q[0];
sx q[0];
rz(2.3755431) q[0];
rz(-pi) q[1];
rz(-1.0988118) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(0.99265487) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0712142) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(2.1945164) q[1];
x q[2];
rz(0.65808987) q[3];
sx q[3];
rz(-2.1562025) q[3];
sx q[3];
rz(2.493849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99469441) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(1.7374932) q[2];
rz(1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(-1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(0.68516723) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72502575) q[0];
sx q[0];
rz(-1.5013114) q[0];
sx q[0];
rz(0.91992232) q[0];
rz(-0.14641996) q[2];
sx q[2];
rz(-0.78260566) q[2];
sx q[2];
rz(-0.27461068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31720686) q[1];
sx q[1];
rz(-0.82815352) q[1];
sx q[1];
rz(-0.87888996) q[1];
x q[2];
rz(1.54322) q[3];
sx q[3];
rz(-2.555489) q[3];
sx q[3];
rz(0.54857777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7489862) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(1.0990934) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(2.4122453) q[0];
rz(-2.8485281) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.1475295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345394) q[0];
sx q[0];
rz(-2.7211186) q[0];
sx q[0];
rz(0.18376952) q[0];
rz(2.2044264) q[2];
sx q[2];
rz(-1.2680149) q[2];
sx q[2];
rz(-0.62919754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9784669) q[1];
sx q[1];
rz(-1.6118057) q[1];
sx q[1];
rz(0.48467111) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.958193) q[3];
sx q[3];
rz(-2.0347188) q[3];
sx q[3];
rz(2.1277609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3532233) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-0.74907556) q[2];
rz(-2.4979112) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(-2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(2.3102982) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(0.39852279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24005323) q[0];
sx q[0];
rz(-2.0365289) q[0];
sx q[0];
rz(2.9630911) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6767098) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(2.523409) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9498082) q[1];
sx q[1];
rz(-2.5439918) q[1];
sx q[1];
rz(-1.530184) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83491915) q[3];
sx q[3];
rz(-2.4432126) q[3];
sx q[3];
rz(-1.8861119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1901671) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.297696) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(-2.4979533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289537) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(1.7522316) q[0];
rz(1.5147491) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(1.0983889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3745981) q[0];
sx q[0];
rz(-0.54531389) q[0];
sx q[0];
rz(2.0314412) q[0];
rz(2.056698) q[2];
sx q[2];
rz(-2.3294318) q[2];
sx q[2];
rz(0.19030262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92022773) q[1];
sx q[1];
rz(-2.5623294) q[1];
sx q[1];
rz(-3.1406162) q[1];
rz(-pi) q[2];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.8212574) q[3];
sx q[3];
rz(-2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(-0.51952726) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-2.675132) q[0];
rz(0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(2.5126273) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0590203) q[0];
sx q[0];
rz(-1.9452381) q[0];
sx q[0];
rz(3.0890205) q[0];
rz(-pi) q[1];
rz(2.2132657) q[2];
sx q[2];
rz(-1.7648342) q[2];
sx q[2];
rz(2.5460668) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2168857) q[1];
sx q[1];
rz(-0.44210426) q[1];
sx q[1];
rz(0.0048203992) q[1];
x q[2];
rz(-1.915515) q[3];
sx q[3];
rz(-2.260672) q[3];
sx q[3];
rz(-2.7310839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(2.0824599) q[2];
rz(-0.6774261) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(-0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8582936) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(-2.8876866) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(0.049747808) q[2];
sx q[2];
rz(-0.72409734) q[2];
sx q[2];
rz(-3.004336) q[2];
rz(-0.98388381) q[3];
sx q[3];
rz(-1.6279396) q[3];
sx q[3];
rz(2.770594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
