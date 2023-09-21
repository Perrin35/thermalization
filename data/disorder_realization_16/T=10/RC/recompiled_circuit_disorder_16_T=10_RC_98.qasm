OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(-2.9876246) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51988039) q[0];
sx q[0];
rz(-1.8147239) q[0];
sx q[0];
rz(1.8983311) q[0];
rz(-pi) q[1];
rz(0.89262427) q[2];
sx q[2];
rz(-1.2783588) q[2];
sx q[2];
rz(-0.48722789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4813862) q[1];
sx q[1];
rz(-1.2848789) q[1];
sx q[1];
rz(1.5155161) q[1];
rz(-pi) q[2];
rz(-1.0584352) q[3];
sx q[3];
rz(-1.5844987) q[3];
sx q[3];
rz(-1.1289568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.9677229) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(-3.0766292) q[0];
rz(0.57463542) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.8992791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984098) q[0];
sx q[0];
rz(-1.3230723) q[0];
sx q[0];
rz(1.7032911) q[0];
x q[1];
rz(2.9827036) q[2];
sx q[2];
rz(-2.1974568) q[2];
sx q[2];
rz(1.6275258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2989267) q[1];
sx q[1];
rz(-1.0416404) q[1];
sx q[1];
rz(1.3859205) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70988016) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-2.3185454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80766455) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(2.5578965) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-2.1247991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25123888) q[0];
sx q[0];
rz(-1.6523424) q[0];
sx q[0];
rz(-1.606705) q[0];
x q[1];
rz(2.0605893) q[2];
sx q[2];
rz(-2.042633) q[2];
sx q[2];
rz(0.29957289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2336725) q[1];
sx q[1];
rz(-2.3111812) q[1];
sx q[1];
rz(-2.9552712) q[1];
x q[2];
rz(0.049116491) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(-2.1360872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(-2.0920848) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-0.24681117) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8762159) q[0];
sx q[0];
rz(-2.2582158) q[0];
sx q[0];
rz(2.9486297) q[0];
rz(-pi) q[1];
rz(-0.47841448) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(-2.0699376) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4775131) q[1];
sx q[1];
rz(-1.2004566) q[1];
sx q[1];
rz(1.5143637) q[1];
x q[2];
rz(2.187192) q[3];
sx q[3];
rz(-1.3351001) q[3];
sx q[3];
rz(-2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.3274308) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(2.8932103) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7800956) q[0];
sx q[0];
rz(-1.0307612) q[0];
sx q[0];
rz(1.6105152) q[0];
rz(-pi) q[1];
rz(2.9551198) q[2];
sx q[2];
rz(-1.4154134) q[2];
sx q[2];
rz(1.1764256) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69265428) q[1];
sx q[1];
rz(-0.9874953) q[1];
sx q[1];
rz(1.7765691) q[1];
rz(0.6487209) q[3];
sx q[3];
rz(-1.7626764) q[3];
sx q[3];
rz(-0.81024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7053232) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(0.1246917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.161675) q[0];
sx q[0];
rz(-1.3945762) q[0];
sx q[0];
rz(-1.4833223) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6090392) q[2];
sx q[2];
rz(-1.3281203) q[2];
sx q[2];
rz(2.2720624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4253937) q[1];
sx q[1];
rz(-0.73912207) q[1];
sx q[1];
rz(3.0576586) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3966339) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(-1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(2.0708864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844855) q[0];
sx q[0];
rz(-2.6575408) q[0];
sx q[0];
rz(-0.0012782106) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4378909) q[2];
sx q[2];
rz(-1.226236) q[2];
sx q[2];
rz(2.8466356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0106525) q[1];
sx q[1];
rz(-2.3619235) q[1];
sx q[1];
rz(2.9995549) q[1];
rz(-pi) q[2];
rz(2.8093852) q[3];
sx q[3];
rz(-1.5899961) q[3];
sx q[3];
rz(-1.4022049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(2.8179742) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(2.0543082) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-0.94747296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763801) q[0];
sx q[0];
rz(-2.1064261) q[0];
sx q[0];
rz(-3.0239848) q[0];
x q[1];
rz(1.4119554) q[2];
sx q[2];
rz(-1.3961332) q[2];
sx q[2];
rz(-2.5607002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2944813) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(2.860445) q[1];
x q[2];
rz(0.8074699) q[3];
sx q[3];
rz(-1.315457) q[3];
sx q[3];
rz(1.0367928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(2.6718111) q[2];
rz(1.8404768) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(1.3289733) q[0];
rz(0.7912311) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-0.20283094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77712599) q[0];
sx q[0];
rz(-1.6098032) q[0];
sx q[0];
rz(-1.5447306) q[0];
x q[1];
rz(-1.780026) q[2];
sx q[2];
rz(-1.3522569) q[2];
sx q[2];
rz(-2.2184559) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.087346615) q[1];
sx q[1];
rz(-0.11212238) q[1];
sx q[1];
rz(-2.0263158) q[1];
rz(-pi) q[2];
rz(-0.30236249) q[3];
sx q[3];
rz(-0.48067579) q[3];
sx q[3];
rz(-0.16929786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(2.3341808) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080169454) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(-0.28082401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107287) q[0];
sx q[0];
rz(-2.4952336) q[0];
sx q[0];
rz(-2.3848563) q[0];
rz(0.73239399) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(-0.64901272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2320497) q[1];
sx q[1];
rz(-0.068040158) q[1];
sx q[1];
rz(-1.0104936) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2383575) q[3];
sx q[3];
rz(-1.3849764) q[3];
sx q[3];
rz(-0.18248724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-2.5785057) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6476718) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(2.5429824) q[2];
sx q[2];
rz(-0.90000464) q[2];
sx q[2];
rz(2.4760751) q[2];
rz(0.35691805) q[3];
sx q[3];
rz(-1.9580943) q[3];
sx q[3];
rz(3.0860268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];