OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(3.1199772) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0776805) q[0];
sx q[0];
rz(-0.67718107) q[0];
sx q[0];
rz(-1.022524) q[0];
rz(3.0388019) q[2];
sx q[2];
rz(-1.3445026) q[2];
sx q[2];
rz(0.13470995) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.057030765) q[1];
sx q[1];
rz(-1.7006526) q[1];
sx q[1];
rz(0.76743857) q[1];
rz(2.1943135) q[3];
sx q[3];
rz(-1.2926658) q[3];
sx q[3];
rz(-0.55752414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33064476) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(-2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9379399) q[0];
sx q[0];
rz(-1.5601888) q[0];
sx q[0];
rz(-1.6158576) q[0];
x q[1];
rz(-2.4992141) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(-0.013052879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.590608) q[1];
sx q[1];
rz(-1.3270757) q[1];
sx q[1];
rz(-0.79382146) q[1];
rz(-0.9886338) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(0.31103381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(0.70811159) q[2];
rz(2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-2.3143342) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(2.0522096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3141146) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(0.88721888) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1459648) q[2];
sx q[2];
rz(-2.1263188) q[2];
sx q[2];
rz(-0.59031634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5730126) q[1];
sx q[1];
rz(-1.9553361) q[1];
sx q[1];
rz(0.09210715) q[1];
rz(-pi) q[2];
rz(0.79077625) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(-0.63000597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(-1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(2.4147721) q[0];
rz(0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-0.35983905) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.115288) q[0];
sx q[0];
rz(-1.3867154) q[0];
sx q[0];
rz(0.91362761) q[0];
rz(-2.4904576) q[2];
sx q[2];
rz(-1.6614117) q[2];
sx q[2];
rz(3.0504984) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53107809) q[1];
sx q[1];
rz(-2.5602166) q[1];
sx q[1];
rz(-2.4060712) q[1];
x q[2];
rz(-1.9197649) q[3];
sx q[3];
rz(-1.2811321) q[3];
sx q[3];
rz(-1.7581913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-2.3763669) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9969479) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(-1.7255406) q[0];
rz(0.36711806) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9177592) q[0];
sx q[0];
rz(-0.98063722) q[0];
sx q[0];
rz(-2.607164) q[0];
x q[1];
rz(-2.6501422) q[2];
sx q[2];
rz(-1.5361668) q[2];
sx q[2];
rz(0.92726842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1977735) q[1];
sx q[1];
rz(-0.10995956) q[1];
sx q[1];
rz(0.52093671) q[1];
x q[2];
rz(0.27541311) q[3];
sx q[3];
rz(-1.7731035) q[3];
sx q[3];
rz(-1.3261258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96200213) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(1.996421) q[0];
rz(2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(2.9343658) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4743054) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(0.62367237) q[0];
rz(-pi) q[1];
rz(-2.5404262) q[2];
sx q[2];
rz(-1.9631557) q[2];
sx q[2];
rz(2.5208134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44240272) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(-1.8276023) q[1];
x q[2];
rz(-2.9793671) q[3];
sx q[3];
rz(-2.3174006) q[3];
sx q[3];
rz(1.970286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-0.72171372) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(0.26003626) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-0.2917372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7043982) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-2.722446) q[0];
rz(-pi) q[1];
rz(-1.8998434) q[2];
sx q[2];
rz(-2.8223158) q[2];
sx q[2];
rz(2.3925051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3953494) q[1];
sx q[1];
rz(-0.33848539) q[1];
sx q[1];
rz(2.6836718) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8802059) q[3];
sx q[3];
rz(-1.0687807) q[3];
sx q[3];
rz(1.7878143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-0.83731246) q[2];
rz(1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(-3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(-2.1829139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6129235) q[0];
sx q[0];
rz(-1.0597502) q[0];
sx q[0];
rz(-0.59707609) q[0];
x q[1];
rz(-1.1218698) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(1.6758855) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.465185) q[1];
sx q[1];
rz(-1.0470111) q[1];
sx q[1];
rz(-3.0182059) q[1];
rz(-pi) q[2];
rz(2.557425) q[3];
sx q[3];
rz(-1.468588) q[3];
sx q[3];
rz(-2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.1834043) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(-0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(2.6760496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532928) q[0];
sx q[0];
rz(-0.71209891) q[0];
sx q[0];
rz(3.0112991) q[0];
rz(1.1586458) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(-0.6790557) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0034345) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(0.65852965) q[1];
rz(-pi) q[2];
rz(-2.8903928) q[3];
sx q[3];
rz(-1.2217055) q[3];
sx q[3];
rz(-2.420345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(-1.643606) q[2];
rz(-2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964504) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(-0.75138584) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56349194) q[0];
sx q[0];
rz(-0.90257114) q[0];
sx q[0];
rz(-2.8816954) q[0];
rz(-pi) q[1];
rz(-1.8495314) q[2];
sx q[2];
rz(-0.50694743) q[2];
sx q[2];
rz(2.7915733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45422428) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(-0.99517676) q[1];
rz(2.9243745) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(-2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.8021884) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(1.8680686) q[2];
sx q[2];
rz(-2.3625629) q[2];
sx q[2];
rz(0.071803781) q[2];
rz(-1.0675666) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
