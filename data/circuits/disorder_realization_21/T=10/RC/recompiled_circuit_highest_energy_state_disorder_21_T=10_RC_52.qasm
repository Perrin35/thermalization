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
rz(-2.8357808) q[0];
sx q[0];
rz(-0.6337136) q[0];
sx q[0];
rz(2.5330438) q[0];
rz(-0.67148709) q[1];
sx q[1];
rz(3.804764) q[1];
sx q[1];
rz(10.884203) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240294) q[0];
sx q[0];
rz(-1.8052973) q[0];
sx q[0];
rz(2.0962711) q[0];
rz(-pi) q[1];
rz(2.9179108) q[2];
sx q[2];
rz(-1.1304026) q[2];
sx q[2];
rz(-0.5612095) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8512177) q[1];
sx q[1];
rz(-0.77271116) q[1];
sx q[1];
rz(2.1593443) q[1];
rz(-pi) q[2];
rz(2.3816649) q[3];
sx q[3];
rz(-1.4455405) q[3];
sx q[3];
rz(1.6881936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3964316) q[2];
sx q[2];
rz(-2.4691984) q[2];
sx q[2];
rz(0.82610899) q[2];
rz(-2.7769026) q[3];
sx q[3];
rz(-1.5267905) q[3];
sx q[3];
rz(0.29432347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9593338) q[0];
sx q[0];
rz(-1.2428281) q[0];
sx q[0];
rz(0.44573319) q[0];
rz(-2.5511197) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(0.42764923) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3217309) q[0];
sx q[0];
rz(-0.45924265) q[0];
sx q[0];
rz(-1.6769354) q[0];
x q[1];
rz(2.1600464) q[2];
sx q[2];
rz(-0.84651154) q[2];
sx q[2];
rz(1.2638448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97274894) q[1];
sx q[1];
rz(-2.4764937) q[1];
sx q[1];
rz(2.0013627) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1366051) q[3];
sx q[3];
rz(-0.79646275) q[3];
sx q[3];
rz(0.69506587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4599956) q[2];
sx q[2];
rz(-2.5619016) q[2];
sx q[2];
rz(-1.0350636) q[2];
rz(1.5486108) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(0.88919324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11397938) q[0];
sx q[0];
rz(-3.1359105) q[0];
sx q[0];
rz(-3.0300544) q[0];
rz(1.5882209) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(-2.9389971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3521393) q[0];
sx q[0];
rz(-1.4367665) q[0];
sx q[0];
rz(-1.9699957) q[0];
rz(-pi) q[1];
rz(-2.8526788) q[2];
sx q[2];
rz(-1.2855069) q[2];
sx q[2];
rz(-2.4288175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9249768) q[1];
sx q[1];
rz(-1.0362715) q[1];
sx q[1];
rz(0.85304867) q[1];
rz(-pi) q[2];
rz(0.59059322) q[3];
sx q[3];
rz(-1.9265129) q[3];
sx q[3];
rz(-1.4416665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4488039) q[2];
sx q[2];
rz(-1.0705798) q[2];
sx q[2];
rz(0.98133522) q[2];
rz(2.8547309) q[3];
sx q[3];
rz(-0.57932866) q[3];
sx q[3];
rz(2.9441492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2364748) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(2.4531051) q[0];
rz(-0.14432898) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(-0.7483288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8546974) q[0];
sx q[0];
rz(-2.4208491) q[0];
sx q[0];
rz(0.067528226) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2950254) q[2];
sx q[2];
rz(-2.7477816) q[2];
sx q[2];
rz(1.9005601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0717582) q[1];
sx q[1];
rz(-0.35440517) q[1];
sx q[1];
rz(-0.56170424) q[1];
x q[2];
rz(2.2714628) q[3];
sx q[3];
rz(-1.3324454) q[3];
sx q[3];
rz(2.7173661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.051108483) q[2];
sx q[2];
rz(-1.313831) q[2];
sx q[2];
rz(-1.2468106) q[2];
rz(0.11940739) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(-2.8363805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7684105) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(0.61491948) q[0];
rz(0.83456314) q[1];
sx q[1];
rz(-1.395023) q[1];
sx q[1];
rz(2.3051197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040720721) q[0];
sx q[0];
rz(-0.03558579) q[0];
sx q[0];
rz(-1.8338704) q[0];
x q[1];
rz(0.70186285) q[2];
sx q[2];
rz(-2.3387675) q[2];
sx q[2];
rz(-0.35477625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53764105) q[1];
sx q[1];
rz(-1.0102235) q[1];
sx q[1];
rz(-3.0106697) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3616124) q[3];
sx q[3];
rz(-2.2298919) q[3];
sx q[3];
rz(-2.7798228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(0.92030805) q[2];
rz(0.80095428) q[3];
sx q[3];
rz(-2.6194173) q[3];
sx q[3];
rz(2.8029158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9698708) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(-0.26611662) q[0];
rz(-1.7099821) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(-2.5439579) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.584883) q[0];
sx q[0];
rz(-2.1765255) q[0];
sx q[0];
rz(-1.8023876) q[0];
rz(-0.18103894) q[2];
sx q[2];
rz(-1.5818051) q[2];
sx q[2];
rz(-0.77048466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3417334) q[1];
sx q[1];
rz(-1.6266154) q[1];
sx q[1];
rz(2.4444405) q[1];
rz(-pi) q[2];
rz(1.2500171) q[3];
sx q[3];
rz(-1.4238384) q[3];
sx q[3];
rz(-3.1023417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1282244) q[2];
sx q[2];
rz(-2.8743663) q[2];
sx q[2];
rz(2.5426148) q[2];
rz(0.61601764) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(1.0429355) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041158572) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(-0.11095412) q[1];
sx q[1];
rz(-1.3919421) q[1];
sx q[1];
rz(2.9999733) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34401152) q[0];
sx q[0];
rz(-2.251029) q[0];
sx q[0];
rz(1.6949095) q[0];
rz(-pi) q[1];
rz(-1.5779823) q[2];
sx q[2];
rz(-0.55139667) q[2];
sx q[2];
rz(1.0137272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1846022) q[1];
sx q[1];
rz(-0.58035589) q[1];
sx q[1];
rz(-0.34837153) q[1];
rz(1.8291437) q[3];
sx q[3];
rz(-2.606539) q[3];
sx q[3];
rz(-0.38219562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5854599) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(-0.11036135) q[2];
rz(2.6335671) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(-1.9917816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13462774) q[0];
sx q[0];
rz(-2.4899794) q[0];
sx q[0];
rz(1.2192669) q[0];
rz(0.02267516) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(0.55331826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0344201) q[0];
sx q[0];
rz(-1.4822472) q[0];
sx q[0];
rz(2.2293985) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99626096) q[2];
sx q[2];
rz(-1.8991514) q[2];
sx q[2];
rz(2.7789214) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4412429) q[1];
sx q[1];
rz(-1.2257595) q[1];
sx q[1];
rz(0.54919589) q[1];
rz(-pi) q[2];
rz(0.24711547) q[3];
sx q[3];
rz(-1.860642) q[3];
sx q[3];
rz(-0.71699504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40552178) q[2];
sx q[2];
rz(-0.61554474) q[2];
sx q[2];
rz(-1.7919398) q[2];
rz(-3.0810629) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(2.5133666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8796006) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(3.1349728) q[0];
rz(2.5026542) q[1];
sx q[1];
rz(-2.9220351) q[1];
sx q[1];
rz(2.7117597) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2868067) q[0];
sx q[0];
rz(-1.9298975) q[0];
sx q[0];
rz(-0.39639985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8877385) q[2];
sx q[2];
rz(-1.9905258) q[2];
sx q[2];
rz(1.3042637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2158969) q[1];
sx q[1];
rz(-1.5924504) q[1];
sx q[1];
rz(0.78288659) q[1];
rz(-3.0379574) q[3];
sx q[3];
rz(-2.5536827) q[3];
sx q[3];
rz(-2.8757446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(-3.1268934) q[2];
rz(-2.0059351) q[3];
sx q[3];
rz(-1.1123927) q[3];
sx q[3];
rz(-1.8358102) q[3];
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
rz(-3.084918) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(1.8490476) q[0];
rz(-0.1560642) q[1];
sx q[1];
rz(-2.3239457) q[1];
sx q[1];
rz(2.3110716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9915402) q[0];
sx q[0];
rz(-2.4084598) q[0];
sx q[0];
rz(1.4218279) q[0];
rz(1.8761329) q[2];
sx q[2];
rz(-2.3407901) q[2];
sx q[2];
rz(2.3126471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1546741) q[1];
sx q[1];
rz(-2.9103644) q[1];
sx q[1];
rz(-1.926941) q[1];
rz(2.7748608) q[3];
sx q[3];
rz(-1.368513) q[3];
sx q[3];
rz(2.3679581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0326651) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(2.6984974) q[2];
rz(2.9923934) q[3];
sx q[3];
rz(-2.7982893) q[3];
sx q[3];
rz(-0.33299115) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0014342) q[0];
sx q[0];
rz(-0.59068155) q[0];
sx q[0];
rz(-0.97472192) q[0];
rz(-2.1433266) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(-2.7985991) q[2];
sx q[2];
rz(-1.3806812) q[2];
sx q[2];
rz(1.5594202) q[2];
rz(-2.2343288) q[3];
sx q[3];
rz(-1.915692) q[3];
sx q[3];
rz(0.23585868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
