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
rz(0.30581185) q[0];
sx q[0];
rz(-2.5078791) q[0];
sx q[0];
rz(0.60854882) q[0];
rz(-0.67148709) q[1];
sx q[1];
rz(3.804764) q[1];
sx q[1];
rz(10.884203) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5175633) q[0];
sx q[0];
rz(-1.3362954) q[0];
sx q[0];
rz(2.0962711) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9179108) q[2];
sx q[2];
rz(-1.1304026) q[2];
sx q[2];
rz(0.5612095) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2903749) q[1];
sx q[1];
rz(-0.77271116) q[1];
sx q[1];
rz(2.1593443) q[1];
x q[2];
rz(1.3988092) q[3];
sx q[3];
rz(-2.3233052) q[3];
sx q[3];
rz(0.23556544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7451611) q[2];
sx q[2];
rz(-2.4691984) q[2];
sx q[2];
rz(0.82610899) q[2];
rz(-2.7769026) q[3];
sx q[3];
rz(-1.6148022) q[3];
sx q[3];
rz(-0.29432347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9593338) q[0];
sx q[0];
rz(-1.8987645) q[0];
sx q[0];
rz(2.6958595) q[0];
rz(2.5511197) q[1];
sx q[1];
rz(-2.1069374) q[1];
sx q[1];
rz(-2.7139434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70156258) q[0];
sx q[0];
rz(-1.1143346) q[0];
sx q[0];
rz(3.0892526) q[0];
x q[1];
rz(0.98154624) q[2];
sx q[2];
rz(-0.84651154) q[2];
sx q[2];
rz(1.8777478) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44438293) q[1];
sx q[1];
rz(-0.97543179) q[1];
sx q[1];
rz(-0.31636379) q[1];
x q[2];
rz(-2.3451361) q[3];
sx q[3];
rz(-1.5743619) q[3];
sx q[3];
rz(-2.2623747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6815971) q[2];
sx q[2];
rz(-0.57969105) q[2];
sx q[2];
rz(-2.106529) q[2];
rz(-1.5486108) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(-0.88919324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11397938) q[0];
sx q[0];
rz(-0.005682156) q[0];
sx q[0];
rz(-0.11153829) q[0];
rz(-1.5533718) q[1];
sx q[1];
rz(-2.7370743) q[1];
sx q[1];
rz(2.9389971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78945335) q[0];
sx q[0];
rz(-1.4367665) q[0];
sx q[0];
rz(-1.9699957) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3417046) q[2];
sx q[2];
rz(-2.7383907) q[2];
sx q[2];
rz(0.099911913) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2059171) q[1];
sx q[1];
rz(-2.1723803) q[1];
sx q[1];
rz(0.66608534) q[1];
rz(-pi) q[2];
rz(0.58835619) q[3];
sx q[3];
rz(-2.4632471) q[3];
sx q[3];
rz(2.7916186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4488039) q[2];
sx q[2];
rz(-1.0705798) q[2];
sx q[2];
rz(-0.98133522) q[2];
rz(2.8547309) q[3];
sx q[3];
rz(-2.562264) q[3];
sx q[3];
rz(0.19744344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9051179) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(2.4531051) q[0];
rz(0.14432898) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(-2.3932638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37667021) q[0];
sx q[0];
rz(-0.85205305) q[0];
sx q[0];
rz(1.5115949) q[0];
x q[1];
rz(0.84656723) q[2];
sx q[2];
rz(-0.39381105) q[2];
sx q[2];
rz(1.9005601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0717582) q[1];
sx q[1];
rz(-0.35440517) q[1];
sx q[1];
rz(-2.5798884) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87012989) q[3];
sx q[3];
rz(-1.3324454) q[3];
sx q[3];
rz(-2.7173661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0904842) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(1.8947821) q[2];
rz(3.0221853) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(2.8363805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.3731821) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(2.5266732) q[0];
rz(-2.3070295) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(0.83647299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2671605) q[0];
sx q[0];
rz(-1.5800486) q[0];
sx q[0];
rz(1.6051588) q[0];
rz(-0.70186285) q[2];
sx q[2];
rz(-2.3387675) q[2];
sx q[2];
rz(-2.7868164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.96326085) q[1];
sx q[1];
rz(-1.6815917) q[1];
sx q[1];
rz(-2.1352467) q[1];
x q[2];
rz(0.74252601) q[3];
sx q[3];
rz(-2.1603185) q[3];
sx q[3];
rz(1.7537011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5323083) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(-2.2212846) q[2];
rz(-0.80095428) q[3];
sx q[3];
rz(-0.52217537) q[3];
sx q[3];
rz(-0.33867684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9698708) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(-0.26611662) q[0];
rz(-1.4316106) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(-0.59763479) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.584883) q[0];
sx q[0];
rz(-0.96506715) q[0];
sx q[0];
rz(-1.8023876) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18103894) q[2];
sx q[2];
rz(-1.5818051) q[2];
sx q[2];
rz(-2.371108) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1823765) q[1];
sx q[1];
rz(-0.87494779) q[1];
sx q[1];
rz(1.6435502) q[1];
rz(-1.8915755) q[3];
sx q[3];
rz(-1.4238384) q[3];
sx q[3];
rz(0.039250936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0133682) q[2];
sx q[2];
rz(-2.8743663) q[2];
sx q[2];
rz(-2.5426148) q[2];
rz(0.61601764) q[3];
sx q[3];
rz(-2.7003761) q[3];
sx q[3];
rz(2.0986572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.041158572) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(-0.11095412) q[1];
sx q[1];
rz(-1.3919421) q[1];
sx q[1];
rz(-0.14161938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6017766) q[0];
sx q[0];
rz(-2.4519073) q[0];
sx q[0];
rz(-0.15183555) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5779823) q[2];
sx q[2];
rz(-2.590196) q[2];
sx q[2];
rz(-1.0137272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31894058) q[1];
sx q[1];
rz(-1.7590862) q[1];
sx q[1];
rz(-2.5892808) q[1];
x q[2];
rz(-2.0911679) q[3];
sx q[3];
rz(-1.7014352) q[3];
sx q[3];
rz(-1.4121233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5854599) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(0.50802556) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(-1.149811) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13462774) q[0];
sx q[0];
rz(-0.65161324) q[0];
sx q[0];
rz(1.9223258) q[0];
rz(3.1189175) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(2.5882744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39529453) q[0];
sx q[0];
rz(-2.2263701) q[0];
sx q[0];
rz(0.11179278) q[0];
rz(0.3855386) q[2];
sx q[2];
rz(-2.1111167) q[2];
sx q[2];
rz(-1.413942) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.63444126) q[1];
sx q[1];
rz(-2.5025625) q[1];
sx q[1];
rz(-2.5386058) q[1];
x q[2];
rz(-2.257709) q[3];
sx q[3];
rz(-2.7629768) q[3];
sx q[3];
rz(0.0061638262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40552178) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(1.7919398) q[2];
rz(-3.0810629) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(-0.6282261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8796006) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(0.0066198786) q[0];
rz(0.63893843) q[1];
sx q[1];
rz(-0.21955755) q[1];
sx q[1];
rz(2.7117597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.00365) q[0];
sx q[0];
rz(-1.2009504) q[0];
sx q[0];
rz(1.957264) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8877385) q[2];
sx q[2];
rz(-1.9905258) q[2];
sx q[2];
rz(-1.3042637) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4749516) q[1];
sx q[1];
rz(-2.3534497) q[1];
sx q[1];
rz(1.6013384) q[1];
rz(-pi) q[2];
rz(1.5019526) q[3];
sx q[3];
rz(-0.9864583) q[3];
sx q[3];
rz(-0.39019728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85475737) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(-3.1268934) q[2];
rz(2.0059351) q[3];
sx q[3];
rz(-1.1123927) q[3];
sx q[3];
rz(1.8358102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056674615) q[0];
sx q[0];
rz(-2.882353) q[0];
sx q[0];
rz(-1.8490476) q[0];
rz(0.1560642) q[1];
sx q[1];
rz(-2.3239457) q[1];
sx q[1];
rz(-2.3110716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6097858) q[0];
sx q[0];
rz(-1.6702819) q[0];
sx q[0];
rz(-2.2983944) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2654598) q[2];
sx q[2];
rz(-2.3407901) q[2];
sx q[2];
rz(-0.82894553) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78962425) q[1];
sx q[1];
rz(-1.3543173) q[1];
sx q[1];
rz(-0.08190544) q[1];
x q[2];
rz(-1.3545348) q[3];
sx q[3];
rz(-1.929708) q[3];
sx q[3];
rz(-0.87417904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0326651) q[2];
sx q[2];
rz(-2.445161) q[2];
sx q[2];
rz(0.44309524) q[2];
rz(-2.9923934) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(2.8086015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14015848) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(-2.1433266) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(-2.6218565) q[2];
sx q[2];
rz(-0.39032857) q[2];
sx q[2];
rz(-2.6662568) q[2];
rz(2.0988437) q[3];
sx q[3];
rz(-0.73560148) q[3];
sx q[3];
rz(-0.92675496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
