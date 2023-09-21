OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6591588) q[0];
sx q[0];
rz(-1.4855488) q[0];
sx q[0];
rz(1.8929338) q[0];
rz(-pi) q[1];
rz(2.2828322) q[2];
sx q[2];
rz(-2.3220064) q[2];
sx q[2];
rz(-0.31529266) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8391708) q[1];
sx q[1];
rz(-1.9170554) q[1];
sx q[1];
rz(-2.195921) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46408848) q[3];
sx q[3];
rz(-0.9872735) q[3];
sx q[3];
rz(0.97809631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(-2.5855529) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-2.1957943) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(0.10903407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5599247) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(1.3397564) q[0];
x q[1];
rz(1.3729587) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(-2.3193662) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2687159) q[1];
sx q[1];
rz(-1.1040338) q[1];
sx q[1];
rz(2.4750535) q[1];
rz(-pi) q[2];
rz(-0.52069943) q[3];
sx q[3];
rz(-1.7215014) q[3];
sx q[3];
rz(3.0050788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(-1.9272778) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(-0.24761565) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(0.48167357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96268481) q[0];
sx q[0];
rz(-2.509153) q[0];
sx q[0];
rz(-0.90895598) q[0];
x q[1];
rz(2.0273676) q[2];
sx q[2];
rz(-2.4764428) q[2];
sx q[2];
rz(0.84012023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.958365) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(-2.4135713) q[1];
x q[2];
rz(2.2948202) q[3];
sx q[3];
rz(-2.0623042) q[3];
sx q[3];
rz(-1.5566952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(-1.2447371) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-3.1367338) q[0];
rz(-pi) q[1];
rz(-1.3696026) q[2];
sx q[2];
rz(-1.4428713) q[2];
sx q[2];
rz(-0.97845562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4501805) q[1];
sx q[1];
rz(-1.0346864) q[1];
sx q[1];
rz(-0.24713534) q[1];
rz(-2.5707385) q[3];
sx q[3];
rz(-1.3802765) q[3];
sx q[3];
rz(-1.966147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(-1.6385471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016184729) q[0];
sx q[0];
rz(-2.5446919) q[0];
sx q[0];
rz(1.3422658) q[0];
rz(-pi) q[1];
rz(-1.4625174) q[2];
sx q[2];
rz(-1.6533972) q[2];
sx q[2];
rz(0.29512197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7170143) q[1];
sx q[1];
rz(-1.6728405) q[1];
sx q[1];
rz(1.2574408) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53346975) q[3];
sx q[3];
rz(-2.5662078) q[3];
sx q[3];
rz(-1.7076147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.903669) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(-0.26671985) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-2.3430603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7245367) q[0];
sx q[0];
rz(-0.31053156) q[0];
sx q[0];
rz(1.0945555) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49662923) q[2];
sx q[2];
rz(-1.2611258) q[2];
sx q[2];
rz(-1.3336099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3244464) q[1];
sx q[1];
rz(-1.518461) q[1];
sx q[1];
rz(-2.1423992) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.800662) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(-2.3199905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(2.6775449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(2.6130136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8115933) q[2];
sx q[2];
rz(-1.5903928) q[2];
sx q[2];
rz(-0.400825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6021767) q[1];
sx q[1];
rz(-2.5578941) q[1];
sx q[1];
rz(-1.0023487) q[1];
rz(-pi) q[2];
rz(-2.9052827) q[3];
sx q[3];
rz(-1.4910306) q[3];
sx q[3];
rz(-2.7304756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(3.1380838) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-2.7602957) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(-1.4415178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7521116) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.5556637) q[0];
x q[1];
rz(-0.054585329) q[2];
sx q[2];
rz(-0.61331257) q[2];
sx q[2];
rz(1.3837981) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.085233363) q[1];
sx q[1];
rz(-1.7006405) q[1];
sx q[1];
rz(-1.4443881) q[1];
rz(-0.12318792) q[3];
sx q[3];
rz(-2.8052969) q[3];
sx q[3];
rz(1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9986481) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(2.3118238) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(-0.98659602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164455) q[0];
sx q[0];
rz(-1.4129606) q[0];
sx q[0];
rz(-2.6437003) q[0];
rz(0.53194745) q[2];
sx q[2];
rz(-0.94594687) q[2];
sx q[2];
rz(-2.9192386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0341067) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-0.60955255) q[1];
rz(-pi) q[2];
rz(-1.9916233) q[3];
sx q[3];
rz(-1.7289366) q[3];
sx q[3];
rz(1.6588253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(2.6221361) q[3];
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
rz(2.7423994) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1499407) q[0];
sx q[0];
rz(-1.9016001) q[0];
sx q[0];
rz(1.8206235) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9344994) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(1.9050913) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.532383) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.7864368) q[1];
x q[2];
rz(1.2763001) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(-0.48081276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.538095) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(2.8812257) q[2];
sx q[2];
rz(-1.7454864) q[2];
sx q[2];
rz(-2.7228552) q[2];
rz(-2.5849871) q[3];
sx q[3];
rz(-1.6986596) q[3];
sx q[3];
rz(-1.2996775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
