OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(3.0732529) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(3.9897444) q[1];
sx q[1];
rz(6.2453649) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4662287) q[0];
sx q[0];
rz(-2.9328406) q[0];
sx q[0];
rz(1.1842313) q[0];
rz(-1.4191364) q[2];
sx q[2];
rz(-2.4745387) q[2];
sx q[2];
rz(-0.59282138) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1483128) q[1];
sx q[1];
rz(-0.82511745) q[1];
sx q[1];
rz(1.6402022) q[1];
rz(2.9904198) q[3];
sx q[3];
rz(-1.2281979) q[3];
sx q[3];
rz(1.3959988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35090703) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(3.0154199) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(-0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.44777563) q[0];
sx q[0];
rz(-2.2651146) q[0];
sx q[0];
rz(-2.5449261) q[0];
rz(-1.5860575) q[1];
sx q[1];
rz(-1.7998453) q[1];
sx q[1];
rz(-1.3635051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3421905) q[0];
sx q[0];
rz(-1.0418833) q[0];
sx q[0];
rz(-2.0687813) q[0];
x q[1];
rz(3.0599669) q[2];
sx q[2];
rz(-0.55742369) q[2];
sx q[2];
rz(-2.4068085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9047782) q[1];
sx q[1];
rz(-1.8493376) q[1];
sx q[1];
rz(2.4836471) q[1];
rz(-pi) q[2];
rz(-2.3724243) q[3];
sx q[3];
rz(-2.2393919) q[3];
sx q[3];
rz(-1.9070966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8186701) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(-2.4070516) q[2];
rz(-0.62720403) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4194141) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(2.869379) q[0];
rz(2.294337) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(2.8289657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011742) q[0];
sx q[0];
rz(-0.22589382) q[0];
sx q[0];
rz(2.9001539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1001415) q[2];
sx q[2];
rz(-2.0690284) q[2];
sx q[2];
rz(2.2217896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12007512) q[1];
sx q[1];
rz(-2.8719963) q[1];
sx q[1];
rz(1.8345548) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39194312) q[3];
sx q[3];
rz(-0.80229811) q[3];
sx q[3];
rz(-0.18303579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12038885) q[2];
sx q[2];
rz(-0.52336064) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(0.96873823) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(0.43280861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.3880436) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(0.62336212) q[0];
rz(0.81758824) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-3.0923016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.277963) q[0];
sx q[0];
rz(-0.81252126) q[0];
sx q[0];
rz(0.066594007) q[0];
rz(-1.2055779) q[2];
sx q[2];
rz(-2.3962002) q[2];
sx q[2];
rz(3.0129907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.019776736) q[1];
sx q[1];
rz(-2.7275118) q[1];
sx q[1];
rz(2.2112234) q[1];
rz(-pi) q[2];
rz(-1.0443346) q[3];
sx q[3];
rz(-1.9658078) q[3];
sx q[3];
rz(-2.928424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7307044) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-0.98497406) q[2];
rz(-2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(0.27339098) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34489283) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(-0.087619089) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(0.87096754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6773416) q[0];
sx q[0];
rz(-2.8528677) q[0];
sx q[0];
rz(0.10084734) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9268553) q[2];
sx q[2];
rz(-0.41741727) q[2];
sx q[2];
rz(0.88127121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4082143) q[1];
sx q[1];
rz(-0.92792643) q[1];
sx q[1];
rz(-1.2182359) q[1];
rz(-0.80645873) q[3];
sx q[3];
rz(-1.6975003) q[3];
sx q[3];
rz(2.5106631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.093420371) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(-2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46309328) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(3.034806) q[0];
rz(-1.1865901) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(-2.1170763) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2183285) q[0];
sx q[0];
rz(-0.54696333) q[0];
sx q[0];
rz(-0.57759072) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1149939) q[2];
sx q[2];
rz(-2.1264646) q[2];
sx q[2];
rz(1.232604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5441106) q[1];
sx q[1];
rz(-2.1253617) q[1];
sx q[1];
rz(-2.5764562) q[1];
x q[2];
rz(-1.939219) q[3];
sx q[3];
rz(-1.107186) q[3];
sx q[3];
rz(-1.1045477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.391905) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(2.270703) q[2];
rz(1.8093367) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9695327) q[0];
sx q[0];
rz(-1.7344069) q[0];
sx q[0];
rz(0.55091888) q[0];
rz(-0.46539601) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-0.23682061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65695545) q[0];
sx q[0];
rz(-1.4716363) q[0];
sx q[0];
rz(-0.037845503) q[0];
rz(-2.4004164) q[2];
sx q[2];
rz(-0.20222649) q[2];
sx q[2];
rz(1.2349595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19141087) q[1];
sx q[1];
rz(-1.4958053) q[1];
sx q[1];
rz(1.8412776) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3053815) q[3];
sx q[3];
rz(-0.85994342) q[3];
sx q[3];
rz(1.4399432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44234309) q[2];
sx q[2];
rz(-0.68168679) q[2];
sx q[2];
rz(-2.701475) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(1.346689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8286164) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(3.112088) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(0.11880076) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1664365) q[0];
sx q[0];
rz(-1.527589) q[0];
sx q[0];
rz(2.8309566) q[0];
rz(2.253184) q[2];
sx q[2];
rz(-2.4783122) q[2];
sx q[2];
rz(2.1036489) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58662215) q[1];
sx q[1];
rz(-1.8645617) q[1];
sx q[1];
rz(-2.657386) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9769067) q[3];
sx q[3];
rz(-2.0593004) q[3];
sx q[3];
rz(-2.1736341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7416731) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(1.5681533) q[2];
rz(-1.9741156) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(1.2742111) q[3];
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
rz(0.90299273) q[0];
sx q[0];
rz(-0.82505834) q[0];
sx q[0];
rz(-2.9883244) q[0];
rz(-2.080147) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.9877888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3183347) q[0];
sx q[0];
rz(-1.2055921) q[0];
sx q[0];
rz(-1.2814786) q[0];
rz(2.597528) q[2];
sx q[2];
rz(-0.96368921) q[2];
sx q[2];
rz(-2.0008759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8037387) q[1];
sx q[1];
rz(-1.5152463) q[1];
sx q[1];
rz(1.7800203) q[1];
rz(-pi) q[2];
rz(3.0204569) q[3];
sx q[3];
rz(-0.88705685) q[3];
sx q[3];
rz(-0.63431206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.29685059) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(-1.6142169) q[2];
rz(-2.1447003) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8940354) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(-2.9515008) q[0];
rz(-2.4709573) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(-1.488283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9207536) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(-3.0160745) q[0];
rz(-pi) q[1];
rz(-1.8411631) q[2];
sx q[2];
rz(-1.1425945) q[2];
sx q[2];
rz(-2.2141475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6675122) q[1];
sx q[1];
rz(-2.0733359) q[1];
sx q[1];
rz(2.4221592) q[1];
rz(-0.58191396) q[3];
sx q[3];
rz(-2.1776608) q[3];
sx q[3];
rz(-2.5556263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0037447475) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(2.2424662) q[2];
rz(2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(-2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.0638194) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(-0.28941119) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(0.32439705) q[2];
sx q[2];
rz(-0.73642052) q[2];
sx q[2];
rz(-3.0036075) q[2];
rz(1.1424941) q[3];
sx q[3];
rz(-0.29853587) q[3];
sx q[3];
rz(1.5034911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
