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
rz(-2.2440417) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9165186) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(2.201447) q[0];
rz(-pi) q[1];
rz(-2.2034982) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(3.0245568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(0.27889241) q[1];
rz(-pi) q[2];
rz(-1.3863871) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(2.409639) q[2];
rz(-2.1814363) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-2.2251341) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(-2.2629471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6966349) q[0];
sx q[0];
rz(-0.96224552) q[0];
sx q[0];
rz(2.6702325) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36388134) q[2];
sx q[2];
rz(-2.4907787) q[2];
sx q[2];
rz(1.770307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6138184) q[1];
sx q[1];
rz(-1.4973745) q[1];
sx q[1];
rz(-2.7223177) q[1];
rz(-pi) q[2];
x q[2];
rz(0.068275498) q[3];
sx q[3];
rz(-1.2618511) q[3];
sx q[3];
rz(-2.6451049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8615222) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6140401) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7911581) q[0];
sx q[0];
rz(-2.3179623) q[0];
sx q[0];
rz(-2.6333548) q[0];
x q[1];
rz(1.7240702) q[2];
sx q[2];
rz(-2.7445265) q[2];
sx q[2];
rz(-3.0561471) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9417291) q[1];
sx q[1];
rz(-2.6911096) q[1];
sx q[1];
rz(2.4099039) q[1];
rz(-pi) q[2];
rz(0.21568732) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(1.4195201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(2.1901954) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(0.11638164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560293) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(-1.602519) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0668683) q[2];
sx q[2];
rz(-1.2328086) q[2];
sx q[2];
rz(2.4525814) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9434005) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(1.1073768) q[1];
rz(-pi) q[2];
rz(-1.5192401) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(-1.6980905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.1486357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7589353) q[0];
sx q[0];
rz(-1.6446582) q[0];
sx q[0];
rz(2.1942684) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1158475) q[2];
sx q[2];
rz(-0.75136853) q[2];
sx q[2];
rz(-0.92912208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83541218) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(2.6447891) q[1];
rz(-pi) q[2];
rz(-2.4162021) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(0.52291742) q[3];
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
rz(2.1095236) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(0.39598879) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9128742) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(2.3939783) q[0];
rz(2.2588737) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(0.27872745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1604662) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(1.6756945) q[1];
rz(1.8984853) q[3];
sx q[3];
rz(-1.7565691) q[3];
sx q[3];
rz(-0.9542619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(-2.5653429) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.2111838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241087) q[0];
sx q[0];
rz(-0.78290126) q[0];
sx q[0];
rz(0.80362513) q[0];
rz(-pi) q[1];
rz(2.6518875) q[2];
sx q[2];
rz(-1.6929132) q[2];
sx q[2];
rz(-2.7011938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8853828) q[1];
sx q[1];
rz(-1.4146283) q[1];
sx q[1];
rz(-3.0728818) q[1];
rz(-pi) q[2];
rz(-0.32254036) q[3];
sx q[3];
rz(-2.1493559) q[3];
sx q[3];
rz(-1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3367735) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-1.0205644) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.221009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2894665) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(-0.37129398) q[0];
x q[1];
rz(-1.7336573) q[2];
sx q[2];
rz(-1.7454299) q[2];
sx q[2];
rz(2.3101431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80291058) q[1];
sx q[1];
rz(-1.435034) q[1];
sx q[1];
rz(-1.6822097) q[1];
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
rz(0.28875479) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(-0.6048454) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-0.12577122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.754697) q[0];
sx q[0];
rz(-1.7594975) q[0];
sx q[0];
rz(-0.17908355) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8081886) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(-2.8508027) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3198493) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-2.7601932) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3689234) q[3];
sx q[3];
rz(-1.2807506) q[3];
sx q[3];
rz(-1.9162852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(0.62018958) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(-2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-2.0013924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70779078) q[0];
sx q[0];
rz(-1.3085438) q[0];
sx q[0];
rz(-3.0387525) q[0];
x q[1];
rz(-2.033038) q[2];
sx q[2];
rz(-1.8046364) q[2];
sx q[2];
rz(-3.0548981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89123594) q[1];
sx q[1];
rz(-1.48238) q[1];
sx q[1];
rz(2.9524515) q[1];
rz(-pi) q[2];
rz(-0.79348989) q[3];
sx q[3];
rz(-1.7975382) q[3];
sx q[3];
rz(-1.2236809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7662979) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(-0.0071772655) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(2.6208411) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(2.3788135) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(-2.5614212) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
