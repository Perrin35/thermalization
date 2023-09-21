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
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(-1.0645359) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9165186) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(0.9401456) q[0];
rz(-pi) q[1];
rz(1.1041553) q[2];
sx q[2];
rz(-2.4541514) q[2];
sx q[2];
rz(1.825037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(2.8627002) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1944794) q[3];
sx q[3];
rz(-0.76768657) q[3];
sx q[3];
rz(0.0038113468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(2.409639) q[2];
rz(-2.1814363) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(-0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(2.2629471) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71683305) q[0];
sx q[0];
rz(-2.3905907) q[0];
sx q[0];
rz(-2.1483833) q[0];
rz(0.61848817) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(2.6478812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2616927) q[1];
sx q[1];
rz(-2.7163134) q[1];
sx q[1];
rz(-2.9628423) q[1];
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
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(-2.98996) q[3];
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
rz(2.6140401) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(-1.4100769) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58223984) q[0];
sx q[0];
rz(-1.2057349) q[0];
sx q[0];
rz(-0.75612005) q[0];
rz(-pi) q[1];
rz(-3.0776575) q[2];
sx q[2];
rz(-1.9629515) q[2];
sx q[2];
rz(2.8901697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0908302) q[1];
sx q[1];
rz(-1.275626) q[1];
sx q[1];
rz(-2.7961618) q[1];
rz(0.21568732) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(1.4195201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48556337) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(1.602519) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0747244) q[2];
sx q[2];
rz(-1.908784) q[2];
sx q[2];
rz(-0.68901125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0013106) q[1];
sx q[1];
rz(-1.1661068) q[1];
sx q[1];
rz(2.6005122) q[1];
rz(-pi) q[2];
rz(3.0451123) q[3];
sx q[3];
rz(-0.49177846) q[3];
sx q[3];
rz(1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(0.25203618) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(-1.9929569) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38265739) q[0];
sx q[0];
rz(-1.4969345) q[0];
sx q[0];
rz(0.9473243) q[0];
rz(-2.6905641) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(0.23652467) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3061805) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(2.6447891) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1762189) q[3];
sx q[3];
rz(-2.2567281) q[3];
sx q[3];
rz(1.8358313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(-1.032069) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(0.2063624) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8800031) q[0];
sx q[0];
rz(-2.3149009) q[0];
sx q[0];
rz(-1.6893301) q[0];
rz(-pi) q[1];
rz(2.2588737) q[2];
sx q[2];
rz(-1.025891) q[2];
sx q[2];
rz(2.8628652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57950912) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(0.097192055) q[1];
rz(0.19595512) q[3];
sx q[3];
rz(-1.2489508) q[3];
sx q[3];
rz(-2.4623507) q[3];
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
rz(-0.74832908) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-1.07871) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.2111838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71935463) q[0];
sx q[0];
rz(-2.1033759) q[0];
sx q[0];
rz(0.60441916) q[0];
rz(-1.4326101) q[2];
sx q[2];
rz(-1.0850564) q[2];
sx q[2];
rz(1.1952343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3019575) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(1.1595999) q[1];
rz(-pi) q[2];
rz(-2.022642) q[3];
sx q[3];
rz(-0.65330905) q[3];
sx q[3];
rz(0.89316955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3367735) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.221009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3295791) q[0];
sx q[0];
rz(-0.81389577) q[0];
sx q[0];
rz(1.9474067) q[0];
rz(0.17692716) q[2];
sx q[2];
rz(-1.4104341) q[2];
sx q[2];
rz(-0.71080506) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.358566) q[1];
sx q[1];
rz(-1.6811803) q[1];
sx q[1];
rz(3.0049938) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2113167) q[3];
sx q[3];
rz(-2.3593138) q[3];
sx q[3];
rz(0.6560916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(2.5382036) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-0.12577122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1544979) q[0];
sx q[0];
rz(-0.25941601) q[0];
sx q[0];
rz(0.82018606) q[0];
x q[1];
rz(1.7681098) q[2];
sx q[2];
rz(-1.8981877) q[2];
sx q[2];
rz(-1.9257853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77482624) q[1];
sx q[1];
rz(-2.3407196) q[1];
sx q[1];
rz(-1.9701387) q[1];
rz(-2.5500507) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(1.8464309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(0.62018958) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70779078) q[0];
sx q[0];
rz(-1.8330488) q[0];
sx q[0];
rz(0.10284013) q[0];
x q[1];
rz(-1.0802202) q[2];
sx q[2];
rz(-2.6274101) q[2];
sx q[2];
rz(1.0487923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(1.4807832) q[1];
rz(-pi) q[2];
rz(-1.8885918) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89467775) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-2.649879) q[2];
sx q[2];
rz(-1.9946873) q[2];
sx q[2];
rz(-0.36110525) q[2];
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