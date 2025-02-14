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
rz(-1.7595093) q[0];
sx q[0];
rz(-2.3926662) q[0];
sx q[0];
rz(1.393526) q[0];
rz(0.8028318) q[1];
sx q[1];
rz(-0.79610151) q[1];
sx q[1];
rz(-2.610745) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40924588) q[0];
sx q[0];
rz(-1.0826001) q[0];
sx q[0];
rz(2.7335579) q[0];
rz(-1.5535001) q[2];
sx q[2];
rz(-1.7375542) q[2];
sx q[2];
rz(1.1823428) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72361354) q[1];
sx q[1];
rz(-2.6043677) q[1];
sx q[1];
rz(-1.4239271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.051732) q[3];
sx q[3];
rz(-1.6593907) q[3];
sx q[3];
rz(-2.5930694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.406245) q[2];
sx q[2];
rz(-0.9333868) q[2];
sx q[2];
rz(-0.99785844) q[2];
rz(2.2570299) q[3];
sx q[3];
rz(-1.1793143) q[3];
sx q[3];
rz(-0.71200371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3957921) q[0];
sx q[0];
rz(-0.53676787) q[0];
sx q[0];
rz(-2.0641548) q[0];
rz(-2.5224345) q[1];
sx q[1];
rz(-1.8749323) q[1];
sx q[1];
rz(1.9177297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.773945) q[0];
sx q[0];
rz(-1.2804693) q[0];
sx q[0];
rz(0.18240697) q[0];
x q[1];
rz(1.5101132) q[2];
sx q[2];
rz(-3.1113495) q[2];
sx q[2];
rz(-1.9918421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.465816) q[1];
sx q[1];
rz(-0.73938939) q[1];
sx q[1];
rz(0.63111102) q[1];
rz(-2.2961462) q[3];
sx q[3];
rz(-0.92638141) q[3];
sx q[3];
rz(-3.0385142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18664843) q[2];
sx q[2];
rz(-0.84543219) q[2];
sx q[2];
rz(-1.0364944) q[2];
rz(-1.9519818) q[3];
sx q[3];
rz(-1.2396953) q[3];
sx q[3];
rz(3.0724683) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67967296) q[0];
sx q[0];
rz(-2.9000977) q[0];
sx q[0];
rz(1.4659721) q[0];
rz(-1.2566603) q[1];
sx q[1];
rz(-0.64777056) q[1];
sx q[1];
rz(0.57674903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38274469) q[0];
sx q[0];
rz(-0.85475105) q[0];
sx q[0];
rz(0.76356739) q[0];
rz(-pi) q[1];
rz(3.0664394) q[2];
sx q[2];
rz(-0.45991969) q[2];
sx q[2];
rz(-3.0004139) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27611342) q[1];
sx q[1];
rz(-0.41655585) q[1];
sx q[1];
rz(-1.6989442) q[1];
rz(-2.0641608) q[3];
sx q[3];
rz(-2.9558792) q[3];
sx q[3];
rz(0.41994219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74167788) q[2];
sx q[2];
rz(-2.8027813) q[2];
sx q[2];
rz(-2.3728288) q[2];
rz(-0.1712884) q[3];
sx q[3];
rz(-1.7007217) q[3];
sx q[3];
rz(-0.35458529) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6866368) q[0];
sx q[0];
rz(-0.88382116) q[0];
sx q[0];
rz(-2.6052642) q[0];
rz(2.3388011) q[1];
sx q[1];
rz(-1.2840459) q[1];
sx q[1];
rz(-3.0435496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5216828) q[0];
sx q[0];
rz(-2.8044716) q[0];
sx q[0];
rz(-1.328892) q[0];
rz(2.1613902) q[2];
sx q[2];
rz(-0.88104311) q[2];
sx q[2];
rz(-0.99926567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9598011) q[1];
sx q[1];
rz(-1.0611294) q[1];
sx q[1];
rz(2.9087429) q[1];
x q[2];
rz(-0.40880568) q[3];
sx q[3];
rz(-1.173436) q[3];
sx q[3];
rz(1.8026343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98196491) q[2];
sx q[2];
rz(-2.1106796) q[2];
sx q[2];
rz(0.31943303) q[2];
rz(-2.5751513) q[3];
sx q[3];
rz(-0.74685493) q[3];
sx q[3];
rz(1.9802861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.40120688) q[0];
sx q[0];
rz(-1.8219319) q[0];
sx q[0];
rz(-2.5886986) q[0];
rz(-0.7792019) q[1];
sx q[1];
rz(-0.82256493) q[1];
sx q[1];
rz(2.353277) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6471651) q[0];
sx q[0];
rz(-2.8719423) q[0];
sx q[0];
rz(-0.2616051) q[0];
x q[1];
rz(0.11190986) q[2];
sx q[2];
rz(-0.21006489) q[2];
sx q[2];
rz(0.0029759759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32994575) q[1];
sx q[1];
rz(-2.0036942) q[1];
sx q[1];
rz(1.2654773) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3230927) q[3];
sx q[3];
rz(-1.2799529) q[3];
sx q[3];
rz(-1.820562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7439482) q[2];
sx q[2];
rz(-2.0136191) q[2];
sx q[2];
rz(-0.20565847) q[2];
rz(-1.3501984) q[3];
sx q[3];
rz(-0.74068991) q[3];
sx q[3];
rz(-0.010644309) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81269294) q[0];
sx q[0];
rz(-0.44726547) q[0];
sx q[0];
rz(1.5146259) q[0];
rz(0.80500066) q[1];
sx q[1];
rz(-1.4470419) q[1];
sx q[1];
rz(2.8256493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63271967) q[0];
sx q[0];
rz(-0.53180662) q[0];
sx q[0];
rz(-0.72152941) q[0];
rz(0.50675372) q[2];
sx q[2];
rz(-0.67595081) q[2];
sx q[2];
rz(-1.6981704) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22623569) q[1];
sx q[1];
rz(-1.2958044) q[1];
sx q[1];
rz(2.8746241) q[1];
x q[2];
rz(2.6942433) q[3];
sx q[3];
rz(-2.7998689) q[3];
sx q[3];
rz(2.0096092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67779764) q[2];
sx q[2];
rz(-0.69910502) q[2];
sx q[2];
rz(-0.15764906) q[2];
rz(-2.6080103) q[3];
sx q[3];
rz(-1.491051) q[3];
sx q[3];
rz(-1.5177479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5428298) q[0];
sx q[0];
rz(-2.6451126) q[0];
sx q[0];
rz(-0.98094034) q[0];
rz(-2.6988103) q[1];
sx q[1];
rz(-1.1863703) q[1];
sx q[1];
rz(-0.47169366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2344246) q[0];
sx q[0];
rz(-1.8629854) q[0];
sx q[0];
rz(-2.0227595) q[0];
rz(-pi) q[1];
rz(1.4002789) q[2];
sx q[2];
rz(-1.2954752) q[2];
sx q[2];
rz(-1.4845276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.737094) q[1];
sx q[1];
rz(-1.5570159) q[1];
sx q[1];
rz(3.0062129) q[1];
rz(-pi) q[2];
rz(1.7552492) q[3];
sx q[3];
rz(-1.9058936) q[3];
sx q[3];
rz(1.9012251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20785759) q[2];
sx q[2];
rz(-2.2722878) q[2];
sx q[2];
rz(3.0233439) q[2];
rz(2.5737428) q[3];
sx q[3];
rz(-1.7852781) q[3];
sx q[3];
rz(2.3387199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0589013) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(-0.68728224) q[0];
rz(1.8742689) q[1];
sx q[1];
rz(-1.9272389) q[1];
sx q[1];
rz(2.5158688) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3104048) q[0];
sx q[0];
rz(-0.021569546) q[0];
sx q[0];
rz(-1.789272) q[0];
x q[1];
rz(2.5753046) q[2];
sx q[2];
rz(-1.0121965) q[2];
sx q[2];
rz(-1.2142912) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81258147) q[1];
sx q[1];
rz(-2.1454403) q[1];
sx q[1];
rz(-3.126815) q[1];
x q[2];
rz(2.1504907) q[3];
sx q[3];
rz(-1.8324346) q[3];
sx q[3];
rz(-1.2708479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9132793) q[2];
sx q[2];
rz(-1.2484756) q[2];
sx q[2];
rz(2.0466364) q[2];
rz(0.45371184) q[3];
sx q[3];
rz(-0.79115051) q[3];
sx q[3];
rz(-0.18610826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5943282) q[0];
sx q[0];
rz(-2.6481977) q[0];
sx q[0];
rz(3.0739947) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.2815579) q[1];
sx q[1];
rz(3.0060815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9477959) q[0];
sx q[0];
rz(-2.3289177) q[0];
sx q[0];
rz(-0.77085797) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4082528) q[2];
sx q[2];
rz(-0.83627273) q[2];
sx q[2];
rz(-2.5275309) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77733946) q[1];
sx q[1];
rz(-0.65636623) q[1];
sx q[1];
rz(-0.47231625) q[1];
x q[2];
rz(-0.28279808) q[3];
sx q[3];
rz(-1.334189) q[3];
sx q[3];
rz(-0.11219926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4542666) q[2];
sx q[2];
rz(-2.8886075) q[2];
sx q[2];
rz(-0.48736408) q[2];
rz(-2.7519915) q[3];
sx q[3];
rz(-2.1890958) q[3];
sx q[3];
rz(-0.065936955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49366632) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(-1.3847463) q[0];
rz(1.0549649) q[1];
sx q[1];
rz(-2.2672548) q[1];
sx q[1];
rz(-2.8816282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1459634) q[0];
sx q[0];
rz(-0.71486799) q[0];
sx q[0];
rz(1.9674106) q[0];
rz(-0.91935888) q[2];
sx q[2];
rz(-0.16460379) q[2];
sx q[2];
rz(2.2155264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.00947) q[1];
sx q[1];
rz(-2.1984221) q[1];
sx q[1];
rz(0.043379003) q[1];
x q[2];
rz(-2.412739) q[3];
sx q[3];
rz(-2.8480004) q[3];
sx q[3];
rz(2.7135682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7197623) q[2];
sx q[2];
rz(-2.7037342) q[2];
sx q[2];
rz(-1.7913294) q[2];
rz(-1.0849902) q[3];
sx q[3];
rz(-1.2144054) q[3];
sx q[3];
rz(2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0706901) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(1.4821953) q[1];
sx q[1];
rz(-1.2980325) q[1];
sx q[1];
rz(1.4549805) q[1];
rz(-2.0878592) q[2];
sx q[2];
rz(-2.0894433) q[2];
sx q[2];
rz(3.1002239) q[2];
rz(-1.0709892) q[3];
sx q[3];
rz(-1.3060112) q[3];
sx q[3];
rz(-1.3597091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
