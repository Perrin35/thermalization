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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8938457) q[0];
sx q[0];
rz(-2.6129299) q[0];
sx q[0];
rz(1.7696487) q[0];
x q[1];
rz(2.0376671) q[2];
sx q[2];
rz(-0.97969998) q[2];
sx q[2];
rz(0.95962722) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8707667) q[1];
sx q[1];
rz(-1.4987136) q[1];
sx q[1];
rz(-1.4429528) q[1];
rz(-pi) q[2];
rz(1.1107221) q[3];
sx q[3];
rz(-0.51957909) q[3];
sx q[3];
rz(-0.67668623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1074368) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(1.6050603) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-0.03604123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(0.56150395) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(-0.5805648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626494) q[0];
sx q[0];
rz(-1.2622841) q[0];
sx q[0];
rz(1.4573775) q[0];
rz(-pi) q[1];
rz(2.850769) q[2];
sx q[2];
rz(-0.65867701) q[2];
sx q[2];
rz(-0.80292279) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3483352) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(-0.72850119) q[1];
x q[2];
rz(2.7892116) q[3];
sx q[3];
rz(-1.1412732) q[3];
sx q[3];
rz(-2.9126715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-2.2632329) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-2.3679249) q[0];
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
rz(2.2680227) q[0];
sx q[0];
rz(-2.7497254) q[0];
sx q[0];
rz(-0.64143945) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.122934) q[2];
sx q[2];
rz(-1.820192) q[2];
sx q[2];
rz(-2.2160335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82959158) q[1];
sx q[1];
rz(-0.81273505) q[1];
sx q[1];
rz(2.9144822) q[1];
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
rz(-1.025528) q[2];
sx q[2];
rz(0.27734217) q[2];
rz(-0.39595655) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(-0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70401496) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(-2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-0.77082005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642826) q[0];
sx q[0];
rz(-1.5543803) q[0];
sx q[0];
rz(-1.5903227) q[0];
x q[1];
rz(-1.3768251) q[2];
sx q[2];
rz(-2.6435404) q[2];
sx q[2];
rz(1.1454524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3168287) q[1];
sx q[1];
rz(-2.0889805) q[1];
sx q[1];
rz(1.1512685) q[1];
rz(-pi) q[2];
rz(2.2724857) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(-1.3645736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(-2.4285994) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.4404526) q[0];
rz(3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-0.17000155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029179137) q[0];
sx q[0];
rz(-0.83054435) q[0];
sx q[0];
rz(1.2519757) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0427809) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69668329) q[1];
sx q[1];
rz(-1.0001567) q[1];
sx q[1];
rz(-2.6615104) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2682297) q[3];
sx q[3];
rz(-1.0358827) q[3];
sx q[3];
rz(1.3267335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(-1.7374932) q[2];
rz(-1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(-0.45853841) q[0];
rz(-0.25587747) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(2.4564254) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2049094) q[0];
sx q[0];
rz(-0.65403599) q[0];
sx q[0];
rz(-1.685164) q[0];
rz(-pi) q[1];
rz(1.7148758) q[2];
sx q[2];
rz(-2.3428168) q[2];
sx q[2];
rz(0.069552334) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31720686) q[1];
sx q[1];
rz(-0.82815352) q[1];
sx q[1];
rz(2.2627027) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5983726) q[3];
sx q[3];
rz(-0.58610361) q[3];
sx q[3];
rz(0.54857777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3926065) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(-2.0424992) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-0.18923047) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012506164) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(2.4122453) q[0];
rz(2.8485281) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(1.1475295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345394) q[0];
sx q[0];
rz(-2.7211186) q[0];
sx q[0];
rz(-2.9578231) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36979923) q[2];
sx q[2];
rz(-0.97019201) q[2];
sx q[2];
rz(-2.415654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16312576) q[1];
sx q[1];
rz(-1.6118057) q[1];
sx q[1];
rz(2.6569215) q[1];
rz(-pi) q[2];
rz(-1.1833997) q[3];
sx q[3];
rz(-2.0347188) q[3];
sx q[3];
rz(-2.1277609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78836936) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(2.4979112) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(-2.3102982) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(0.39852279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14197037) q[0];
sx q[0];
rz(-0.49641434) q[0];
sx q[0];
rz(1.9103785) q[0];
rz(-2.987791) q[2];
sx q[2];
rz(-0.6066583) q[2];
sx q[2];
rz(-2.3369044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1917845) q[1];
sx q[1];
rz(-0.59760082) q[1];
sx q[1];
rz(1.530184) q[1];
rz(2.3066735) q[3];
sx q[3];
rz(-0.69838006) q[3];
sx q[3];
rz(-1.8861119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9514256) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(-1.8438967) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(-0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289537) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.7522316) q[0];
rz(-1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(-2.0432037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84868816) q[0];
sx q[0];
rz(-1.0875889) q[0];
sx q[0];
rz(-2.8781761) q[0];
rz(1.0848947) q[2];
sx q[2];
rz(-0.81216083) q[2];
sx q[2];
rz(0.19030262) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64975148) q[1];
sx q[1];
rz(-1.5713308) q[1];
sx q[1];
rz(-2.5623296) q[1];
rz(-2.5536355) q[3];
sx q[3];
rz(-2.8391317) q[3];
sx q[3];
rz(2.1289701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(-0.51952726) q[2];
rz(-1.8064921) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(-1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-2.675132) q[0];
rz(-2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-0.62896532) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6105729) q[0];
sx q[0];
rz(-1.5218698) q[0];
sx q[0];
rz(-1.945709) q[0];
rz(-pi) q[1];
rz(0.92832698) q[2];
sx q[2];
rz(-1.3767585) q[2];
sx q[2];
rz(-0.59552586) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91937376) q[1];
sx q[1];
rz(-2.0128951) q[1];
sx q[1];
rz(1.5730781) q[1];
rz(-pi) q[2];
rz(1.2260776) q[3];
sx q[3];
rz(-0.88092062) q[3];
sx q[3];
rz(-0.41050875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8978867) q[2];
sx q[2];
rz(-1.0906929) q[2];
sx q[2];
rz(1.0591327) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28329904) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(2.8876866) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(2.4181096) q[2];
sx q[2];
rz(-1.6037446) q[2];
sx q[2];
rz(-1.4708191) q[2];
rz(-1.6737291) q[3];
sx q[3];
rz(-0.58936215) q[3];
sx q[3];
rz(-2.0274558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
