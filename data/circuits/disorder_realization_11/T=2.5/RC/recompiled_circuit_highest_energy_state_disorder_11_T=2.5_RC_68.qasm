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
rz(0.053057916) q[0];
sx q[0];
rz(0.36202708) q[0];
sx q[0];
rz(10.590635) q[0];
rz(2.2447658) q[1];
sx q[1];
rz(4.5936102) q[1];
sx q[1];
rz(11.138154) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0071636) q[0];
sx q[0];
rz(-1.0593793) q[0];
sx q[0];
rz(1.5017139) q[0];
rz(-pi) q[1];
rz(-0.27730242) q[2];
sx q[2];
rz(-2.8805974) q[2];
sx q[2];
rz(2.9948611) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62318351) q[1];
sx q[1];
rz(-1.759638) q[1];
sx q[1];
rz(3.0075226) q[1];
rz(-2.0588097) q[3];
sx q[3];
rz(-1.8086156) q[3];
sx q[3];
rz(0.89624119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4521788) q[2];
sx q[2];
rz(-3.1248326) q[2];
sx q[2];
rz(-2.9407799) q[2];
rz(0.14874841) q[3];
sx q[3];
rz(-3.1368308) q[3];
sx q[3];
rz(2.8301921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.1397322) q[0];
sx q[0];
rz(-2.5480324) q[0];
sx q[0];
rz(-2.1019905) q[0];
rz(-0.014558583) q[1];
sx q[1];
rz(-1.2332375) q[1];
sx q[1];
rz(-1.5537517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4400875) q[0];
sx q[0];
rz(-1.0453512) q[0];
sx q[0];
rz(-2.3129025) q[0];
rz(1.7660242) q[2];
sx q[2];
rz(-0.075610925) q[2];
sx q[2];
rz(-1.4812673) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88513597) q[1];
sx q[1];
rz(-0.26139046) q[1];
sx q[1];
rz(1.4381203) q[1];
rz(-pi) q[2];
rz(-0.33287853) q[3];
sx q[3];
rz(-1.3886222) q[3];
sx q[3];
rz(1.5150439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62850922) q[2];
sx q[2];
rz(-1.6078948) q[2];
sx q[2];
rz(1.3879363) q[2];
rz(1.7657109) q[3];
sx q[3];
rz(-1.0466156) q[3];
sx q[3];
rz(-2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3917711) q[0];
sx q[0];
rz(-2.9048558) q[0];
sx q[0];
rz(-2.5340875) q[0];
rz(1.5433743) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(-2.1764596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34244363) q[0];
sx q[0];
rz(-1.5502872) q[0];
sx q[0];
rz(3.1375225) q[0];
rz(-pi) q[1];
rz(-1.1095382) q[2];
sx q[2];
rz(-1.8806305) q[2];
sx q[2];
rz(0.43333915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9671792) q[1];
sx q[1];
rz(-1.7153746) q[1];
sx q[1];
rz(-3.105393) q[1];
rz(-pi) q[2];
rz(0.85609658) q[3];
sx q[3];
rz(-0.12663791) q[3];
sx q[3];
rz(2.6117725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8049916) q[2];
sx q[2];
rz(-0.6707297) q[2];
sx q[2];
rz(0.86656183) q[2];
rz(1.1066412) q[3];
sx q[3];
rz(-1.5525147) q[3];
sx q[3];
rz(1.6710056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57257819) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(-1.5255852) q[0];
rz(-0.010604803) q[1];
sx q[1];
rz(-0.0037825982) q[1];
sx q[1];
rz(-0.7503646) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.01975) q[0];
sx q[0];
rz(-1.6656309) q[0];
sx q[0];
rz(-0.66480831) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60549824) q[2];
sx q[2];
rz(-1.7586305) q[2];
sx q[2];
rz(2.2772636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1018095) q[1];
sx q[1];
rz(-1.0652115) q[1];
sx q[1];
rz(1.4391446) q[1];
rz(-2.8298805) q[3];
sx q[3];
rz(-1.4501713) q[3];
sx q[3];
rz(-1.8256622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7287207) q[2];
sx q[2];
rz(-1.0888638) q[2];
sx q[2];
rz(-1.9000165) q[2];
rz(-0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(-0.84398213) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4960957) q[0];
sx q[0];
rz(-3.0912919) q[0];
sx q[0];
rz(2.2182933) q[0];
rz(-2.3410489) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(-2.961535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6430017) q[0];
sx q[0];
rz(-2.8992436) q[0];
sx q[0];
rz(-1.7447628) q[0];
rz(-pi) q[1];
rz(1.6512958) q[2];
sx q[2];
rz(-1.4729807) q[2];
sx q[2];
rz(-0.43392935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5023401) q[1];
sx q[1];
rz(-2.709734) q[1];
sx q[1];
rz(-0.60245241) q[1];
rz(-pi) q[2];
rz(-1.2881785) q[3];
sx q[3];
rz(-0.37943951) q[3];
sx q[3];
rz(-1.5490378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30183145) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(-1.7054455) q[2];
rz(-1.885421) q[3];
sx q[3];
rz(-1.6585645) q[3];
sx q[3];
rz(-0.095631599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.489478) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(0.44334626) q[0];
rz(-1.1706785) q[1];
sx q[1];
rz(-0.0010633855) q[1];
sx q[1];
rz(-2.7098157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3928089) q[0];
sx q[0];
rz(-1.8816299) q[0];
sx q[0];
rz(2.3883467) q[0];
x q[1];
rz(1.3987006) q[2];
sx q[2];
rz(-1.4884713) q[2];
sx q[2];
rz(2.4351127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0178294) q[1];
sx q[1];
rz(-2.4430954) q[1];
sx q[1];
rz(-0.61459728) q[1];
rz(-pi) q[2];
rz(2.5434142) q[3];
sx q[3];
rz(-2.7477187) q[3];
sx q[3];
rz(-1.0306213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40054896) q[2];
sx q[2];
rz(-0.86047518) q[2];
sx q[2];
rz(1.384037) q[2];
rz(2.4436229) q[3];
sx q[3];
rz(-2.3845086) q[3];
sx q[3];
rz(2.82011) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31239241) q[0];
sx q[0];
rz(-1.3425403) q[0];
sx q[0];
rz(0.37305748) q[0];
rz(-2.864605) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(2.3854947) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38967237) q[0];
sx q[0];
rz(-0.7451267) q[0];
sx q[0];
rz(-2.5650399) q[0];
rz(-1.937708) q[2];
sx q[2];
rz(-1.9607301) q[2];
sx q[2];
rz(-2.6904047) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3792808) q[1];
sx q[1];
rz(-1.168615) q[1];
sx q[1];
rz(0.97907514) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9916261) q[3];
sx q[3];
rz(-1.2927353) q[3];
sx q[3];
rz(-0.39390491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9296391) q[2];
sx q[2];
rz(-0.55277199) q[2];
sx q[2];
rz(0.92145222) q[2];
rz(0.067342162) q[3];
sx q[3];
rz(-1.9332956) q[3];
sx q[3];
rz(-1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15143722) q[0];
sx q[0];
rz(-2.86148) q[0];
sx q[0];
rz(0.13023278) q[0];
rz(-2.3479334) q[1];
sx q[1];
rz(-0.001611324) q[1];
sx q[1];
rz(2.8614955) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1156688) q[0];
sx q[0];
rz(-1.5172795) q[0];
sx q[0];
rz(-1.6497067) q[0];
rz(-pi) q[1];
rz(-0.21799223) q[2];
sx q[2];
rz(-1.6577066) q[2];
sx q[2];
rz(2.4119968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62781266) q[1];
sx q[1];
rz(-2.4901721) q[1];
sx q[1];
rz(0.86077229) q[1];
x q[2];
rz(-2.6032366) q[3];
sx q[3];
rz(-2.3545485) q[3];
sx q[3];
rz(-0.53883906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.085999504) q[2];
sx q[2];
rz(-1.4693825) q[2];
sx q[2];
rz(-0.80297339) q[2];
rz(1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(-0.83139658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11237535) q[0];
sx q[0];
rz(-0.0028828415) q[0];
sx q[0];
rz(-0.10920864) q[0];
rz(2.7681007) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(2.5901897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5057988) q[0];
sx q[0];
rz(-0.83573816) q[0];
sx q[0];
rz(0.73979017) q[0];
rz(-pi) q[1];
rz(-1.4457747) q[2];
sx q[2];
rz(-1.3879965) q[2];
sx q[2];
rz(-1.3318544) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99433339) q[1];
sx q[1];
rz(-0.87169391) q[1];
sx q[1];
rz(2.2963013) q[1];
rz(-pi) q[2];
rz(0.40355669) q[3];
sx q[3];
rz(-1.3440895) q[3];
sx q[3];
rz(1.3317396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6138844) q[2];
sx q[2];
rz(-1.3091427) q[2];
sx q[2];
rz(-1.8180397) q[2];
rz(1.2644794) q[3];
sx q[3];
rz(-1.8529961) q[3];
sx q[3];
rz(0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5453813) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(0.7363466) q[0];
rz(-2.9296854) q[1];
sx q[1];
rz(-1.0286464) q[1];
sx q[1];
rz(-1.5440936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1101226) q[0];
sx q[0];
rz(-1.8499062) q[0];
sx q[0];
rz(1.5520679) q[0];
rz(-0.030363516) q[2];
sx q[2];
rz(-1.5734451) q[2];
sx q[2];
rz(2.9752258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4344132) q[1];
sx q[1];
rz(-0.79977555) q[1];
sx q[1];
rz(2.4753285) q[1];
rz(-pi) q[2];
rz(-2.5251901) q[3];
sx q[3];
rz(-0.52866259) q[3];
sx q[3];
rz(-1.5746436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.836901) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(1.2738073) q[2];
rz(-1.4443719) q[3];
sx q[3];
rz(-0.074051753) q[3];
sx q[3];
rz(-1.6465638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.973751) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(-1.5372859) q[1];
sx q[1];
rz(-2.2289386) q[1];
sx q[1];
rz(-2.9569721) q[1];
rz(-1.5617076) q[2];
sx q[2];
rz(-1.5307003) q[2];
sx q[2];
rz(-1.8673473) q[2];
rz(-3.01916) q[3];
sx q[3];
rz(-2.1748016) q[3];
sx q[3];
rz(0.55007269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
