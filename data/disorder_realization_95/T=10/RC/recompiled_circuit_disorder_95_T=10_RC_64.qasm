OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(0.89755091) q[0];
rz(2.907213) q[1];
sx q[1];
rz(-2.8657764) q[1];
sx q[1];
rz(-2.0770567) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25047725) q[0];
sx q[0];
rz(-0.63397898) q[0];
sx q[0];
rz(1.4527713) q[0];
rz(-pi) q[1];
rz(2.0374374) q[2];
sx q[2];
rz(-2.4541514) q[2];
sx q[2];
rz(-1.825037) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41649109) q[1];
sx q[1];
rz(-2.2033268) q[1];
sx q[1];
rz(-1.7822669) q[1];
x q[2];
rz(0.75818054) q[3];
sx q[3];
rz(-1.4361793) q[3];
sx q[3];
rz(1.4338223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(-1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.9663548) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-2.2251341) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(-0.8786456) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7322313) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(2.2344927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36388134) q[2];
sx q[2];
rz(-2.4907787) q[2];
sx q[2];
rz(1.3712856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6138184) q[1];
sx q[1];
rz(-1.4973745) q[1];
sx q[1];
rz(-2.7223177) q[1];
x q[2];
rz(-1.3602123) q[3];
sx q[3];
rz(-2.8254291) q[3];
sx q[3];
rz(0.27526835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8615222) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(-1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(2.0203967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7911581) q[0];
sx q[0];
rz(-0.8236304) q[0];
sx q[0];
rz(-0.50823786) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.063935117) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(2.8901697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0908302) q[1];
sx q[1];
rz(-1.8659667) q[1];
sx q[1];
rz(2.7961618) q[1];
rz(-pi) q[2];
rz(-1.6882012) q[3];
sx q[3];
rz(-1.0798287) q[3];
sx q[3];
rz(-1.9672058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(-2.0846941) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43198904) q[0];
sx q[0];
rz(-0.63392144) q[0];
sx q[0];
rz(-0.043179913) q[0];
x q[1];
rz(-2.1999947) q[2];
sx q[2];
rz(-0.59855748) q[2];
sx q[2];
rz(1.7184005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1310281) q[1];
sx q[1];
rz(-2.4782135) q[1];
sx q[1];
rz(2.4478711) q[1];
rz(-pi) q[2];
rz(0.096480358) q[3];
sx q[3];
rz(-0.49177846) q[3];
sx q[3];
rz(-1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.1486357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0557077) q[0];
sx q[0];
rz(-0.62725337) q[0];
sx q[0];
rz(1.6968615) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6905641) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(0.23652467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8685638) q[1];
sx q[1];
rz(-1.3820573) q[1];
sx q[1];
rz(-2.7815458) q[1];
rz(-pi) q[2];
rz(0.43907459) q[3];
sx q[3];
rz(-2.3665161) q[3];
sx q[3];
rz(0.72417688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(-1.032069) q[2];
rz(0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(3.1177974) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(0.39598879) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9128742) q[0];
sx q[0];
rz(-1.657907) q[0];
sx q[0];
rz(0.74761439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80870734) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(-1.2869814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8099765) q[1];
sx q[1];
rz(-0.14252256) q[1];
sx q[1];
rz(-0.82377164) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19595512) q[3];
sx q[3];
rz(-1.2489508) q[3];
sx q[3];
rz(2.4623507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6017194) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71935463) q[0];
sx q[0];
rz(-2.1033759) q[0];
sx q[0];
rz(-0.60441916) q[0];
x q[1];
rz(-2.6518875) q[2];
sx q[2];
rz(-1.4486794) q[2];
sx q[2];
rz(-2.7011938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8163029) q[1];
sx q[1];
rz(-1.5029229) q[1];
sx q[1];
rz(-1.7273278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1738449) q[3];
sx q[3];
rz(-1.839404) q[3];
sx q[3];
rz(-2.8318162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(1.582675) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1655501) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(-2.3483777) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17692716) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(2.4307876) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78302661) q[1];
sx q[1];
rz(-1.4604124) q[1];
sx q[1];
rz(-0.13659887) q[1];
x q[2];
rz(1.2113167) q[3];
sx q[3];
rz(-2.3593138) q[3];
sx q[3];
rz(2.4855011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28875479) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(3.0158214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.754697) q[0];
sx q[0];
rz(-1.7594975) q[0];
sx q[0];
rz(2.9625091) q[0];
x q[1];
rz(-1.3734829) q[2];
sx q[2];
rz(-1.8981877) q[2];
sx q[2];
rz(1.2158074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8217433) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(2.7601932) q[1];
rz(-0.29571663) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(0.28702345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(2.9848849) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-1.1402003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894293) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(1.2055231) q[0];
x q[1];
rz(1.0802202) q[2];
sx q[2];
rz(-0.51418257) q[2];
sx q[2];
rz(1.0487923) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0298437) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(-2.7010121) q[1];
x q[2];
rz(-1.8885918) q[3];
sx q[3];
rz(-0.80298775) q[3];
sx q[3];
rz(0.12249891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(2.2429788) q[2];
rz(-0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-0.5207516) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(2.3788135) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(0.58581523) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
