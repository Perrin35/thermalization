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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39660397) q[0];
sx q[0];
rz(-2.19967) q[0];
sx q[0];
rz(-0.086358503) q[0];
rz(-0.35383309) q[2];
sx q[2];
rz(-2.1733123) q[2];
sx q[2];
rz(-1.2474071) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.861046) q[1];
sx q[1];
rz(-1.7409054) q[1];
sx q[1];
rz(0.64330805) q[1];
rz(-pi) q[2];
rz(2.3834121) q[3];
sx q[3];
rz(-1.4361793) q[3];
sx q[3];
rz(1.7077703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.477318) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-2.409639) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.7094973) q[3];
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
rz(-2.9663548) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-2.2251341) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(2.2629471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40936138) q[0];
sx q[0];
rz(-1.1890113) q[0];
sx q[0];
rz(0.9071) q[0];
rz(2.7777113) q[2];
sx q[2];
rz(-0.65081396) q[2];
sx q[2];
rz(1.770307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0658873) q[1];
sx q[1];
rz(-1.1527219) q[1];
sx q[1];
rz(1.4904406) q[1];
rz(-pi) q[2];
x q[2];
rz(0.068275498) q[3];
sx q[3];
rz(-1.2618511) q[3];
sx q[3];
rz(0.49648778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(-1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(2.0203967) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7911581) q[0];
sx q[0];
rz(-0.8236304) q[0];
sx q[0];
rz(-0.50823786) q[0];
rz(3.0776575) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(-0.25142297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9417291) q[1];
sx q[1];
rz(-2.6911096) q[1];
sx q[1];
rz(0.73168879) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9259053) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(-1.7220725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6275416) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(0.78805584) q[0];
rz(-2.0846941) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560293) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(1.5390736) q[0];
rz(-pi) q[1];
rz(0.38168455) q[2];
sx q[2];
rz(-1.0978062) q[2];
sx q[2];
rz(-0.7009398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1981922) q[1];
sx q[1];
rz(-1.0775837) q[1];
sx q[1];
rz(-1.1073768) q[1];
rz(-0.4898407) q[3];
sx q[3];
rz(-1.6162989) q[3];
sx q[3];
rz(-0.10304606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(-1.9929569) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589353) q[0];
sx q[0];
rz(-1.6446582) q[0];
sx q[0];
rz(-2.1942684) q[0];
x q[1];
rz(0.89678905) q[2];
sx q[2];
rz(-1.2090346) q[2];
sx q[2];
rz(-2.0828431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27302882) q[1];
sx q[1];
rz(-1.7595353) q[1];
sx q[1];
rz(-2.7815458) q[1];
x q[2];
rz(-2.7025181) q[3];
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
rz(-2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(0.023795279) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500568) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43564046) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(3.0138426) q[0];
x q[1];
rz(-0.66531078) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(-0.88924185) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3316162) q[1];
sx q[1];
rz(-0.14252256) q[1];
sx q[1];
rz(0.82377164) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0422802) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(-0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-2.3997967) q[1];
sx q[1];
rz(-1.2111838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21748397) q[0];
sx q[0];
rz(-0.78290126) q[0];
sx q[0];
rz(2.3379675) q[0];
rz(0.48970512) q[2];
sx q[2];
rz(-1.6929132) q[2];
sx q[2];
rz(-0.44039886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2562099) q[1];
sx q[1];
rz(-1.4146283) q[1];
sx q[1];
rz(0.068710879) q[1];
x q[2];
rz(-2.8190523) q[3];
sx q[3];
rz(-0.99223677) q[3];
sx q[3];
rz(1.4417779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3367735) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(-0.15469805) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120136) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(1.194186) q[0];
rz(-pi) q[1];
rz(-1.4079354) q[2];
sx q[2];
rz(-1.7454299) q[2];
sx q[2];
rz(-2.3101431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3386821) q[1];
sx q[1];
rz(-1.435034) q[1];
sx q[1];
rz(-1.6822097) q[1];
rz(0.82151316) q[3];
sx q[3];
rz(-1.320208) q[3];
sx q[3];
rz(-0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9916358) q[0];
sx q[0];
rz(-1.3949252) q[0];
sx q[0];
rz(1.3791023) q[0];
rz(-pi) q[1];
rz(-2.618082) q[2];
sx q[2];
rz(-0.38041174) q[2];
sx q[2];
rz(-1.3695804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51024918) q[1];
sx q[1];
rz(-1.2878839) q[1];
sx q[1];
rz(-2.3307073) q[1];
rz(0.29571663) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(2.8545692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.003309) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-1.1402003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4338019) q[0];
sx q[0];
rz(-1.3085438) q[0];
sx q[0];
rz(-3.0387525) q[0];
rz(2.0613725) q[2];
sx q[2];
rz(-2.6274101) q[2];
sx q[2];
rz(-2.0928004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89123594) q[1];
sx q[1];
rz(-1.6592126) q[1];
sx q[1];
rz(-0.18914117) q[1];
rz(-1.8885918) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(0.89861384) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.7858508) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-0.58017147) q[3];
sx q[3];
rz(-2.4709354) q[3];
sx q[3];
rz(1.2448268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];