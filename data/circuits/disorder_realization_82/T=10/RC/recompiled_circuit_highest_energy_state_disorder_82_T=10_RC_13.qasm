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
rz(-3.1373625) q[0];
sx q[0];
rz(-2.5047996) q[0];
sx q[0];
rz(-1.9544344) q[0];
rz(-2.1597593) q[1];
sx q[1];
rz(-2.667006) q[1];
sx q[1];
rz(-2.2003953) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6146916) q[0];
sx q[0];
rz(-2.6462835) q[0];
sx q[0];
rz(-0.0035347819) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12856483) q[2];
sx q[2];
rz(-1.3701653) q[2];
sx q[2];
rz(-3.0561997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2814215) q[1];
sx q[1];
rz(-2.1263391) q[1];
sx q[1];
rz(1.9940328) q[1];
rz(-pi) q[2];
rz(-2.0273846) q[3];
sx q[3];
rz(-1.5043882) q[3];
sx q[3];
rz(-2.4650064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5473951) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(2.1217864) q[2];
rz(2.8504573) q[3];
sx q[3];
rz(-0.19164339) q[3];
sx q[3];
rz(1.8118793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52784598) q[0];
sx q[0];
rz(-2.3044523) q[0];
sx q[0];
rz(1.6380731) q[0];
rz(1.232052) q[1];
sx q[1];
rz(-2.3161395) q[1];
sx q[1];
rz(-1.253461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952714) q[0];
sx q[0];
rz(-1.1635096) q[0];
sx q[0];
rz(2.8952959) q[0];
rz(1.2792781) q[2];
sx q[2];
rz(-2.0644858) q[2];
sx q[2];
rz(-1.8367298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4644891) q[1];
sx q[1];
rz(-1.3687412) q[1];
sx q[1];
rz(1.4239763) q[1];
x q[2];
rz(0.74442975) q[3];
sx q[3];
rz(-1.9081313) q[3];
sx q[3];
rz(2.6312089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2848009) q[2];
sx q[2];
rz(-2.5992726) q[2];
sx q[2];
rz(2.8112603) q[2];
rz(3.0387943) q[3];
sx q[3];
rz(-1.2448575) q[3];
sx q[3];
rz(2.5225294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5264346) q[0];
sx q[0];
rz(-2.0005083) q[0];
sx q[0];
rz(2.0779628) q[0];
rz(-3.0377153) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(1.1054543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1310467) q[0];
sx q[0];
rz(-1.5233114) q[0];
sx q[0];
rz(3.1094883) q[0];
rz(-1.8806522) q[2];
sx q[2];
rz(-0.83854691) q[2];
sx q[2];
rz(1.0266048) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7522831) q[1];
sx q[1];
rz(-1.0573309) q[1];
sx q[1];
rz(-3.0016243) q[1];
rz(-pi) q[2];
rz(0.01171578) q[3];
sx q[3];
rz(-1.5289617) q[3];
sx q[3];
rz(-1.5277629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1198472) q[2];
sx q[2];
rz(-2.3357119) q[2];
sx q[2];
rz(2.9160807) q[2];
rz(1.3269904) q[3];
sx q[3];
rz(-2.3046389) q[3];
sx q[3];
rz(2.5007611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.91318146) q[0];
sx q[0];
rz(-1.5629733) q[0];
sx q[0];
rz(-2.5033503) q[0];
rz(1.111521) q[1];
sx q[1];
rz(-2.1941954) q[1];
sx q[1];
rz(-0.18724719) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3950662) q[0];
sx q[0];
rz(-1.8573666) q[0];
sx q[0];
rz(-0.92978433) q[0];
rz(-pi) q[1];
rz(0.78088607) q[2];
sx q[2];
rz(-1.5263288) q[2];
sx q[2];
rz(1.5574297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3447002) q[1];
sx q[1];
rz(-0.53965599) q[1];
sx q[1];
rz(-0.15426387) q[1];
rz(-pi) q[2];
rz(0.10113206) q[3];
sx q[3];
rz(-2.4256878) q[3];
sx q[3];
rz(0.16815312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.036666544) q[2];
sx q[2];
rz(-1.5035524) q[2];
sx q[2];
rz(-2.6724913) q[2];
rz(0.30096287) q[3];
sx q[3];
rz(-0.260869) q[3];
sx q[3];
rz(-1.4891967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.8982573) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(-0.54018706) q[0];
rz(-0.46145269) q[1];
sx q[1];
rz(-1.5216454) q[1];
sx q[1];
rz(2.8823421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.777939) q[0];
sx q[0];
rz(-2.724132) q[0];
sx q[0];
rz(2.4221942) q[0];
rz(-0.86452338) q[2];
sx q[2];
rz(-2.5062423) q[2];
sx q[2];
rz(-0.29624789) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4673678) q[1];
sx q[1];
rz(-1.4150054) q[1];
sx q[1];
rz(-2.1464661) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9018039) q[3];
sx q[3];
rz(-1.3223303) q[3];
sx q[3];
rz(1.0120165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3429012) q[2];
sx q[2];
rz(-2.4553802) q[2];
sx q[2];
rz(-2.2314824) q[2];
rz(2.5765007) q[3];
sx q[3];
rz(-1.7730954) q[3];
sx q[3];
rz(-0.71649396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703576) q[0];
sx q[0];
rz(-2.7483181) q[0];
sx q[0];
rz(1.2275335) q[0];
rz(0.50499376) q[1];
sx q[1];
rz(-1.0161437) q[1];
sx q[1];
rz(-1.203677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01194309) q[0];
sx q[0];
rz(-0.27515689) q[0];
sx q[0];
rz(-2.6889474) q[0];
rz(-pi) q[1];
rz(-2.1799548) q[2];
sx q[2];
rz(-0.74660245) q[2];
sx q[2];
rz(-0.37140977) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5281637) q[1];
sx q[1];
rz(-2.2417896) q[1];
sx q[1];
rz(0.74929418) q[1];
rz(1.0685445) q[3];
sx q[3];
rz(-2.6274791) q[3];
sx q[3];
rz(-0.83465605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1167072) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(2.1825979) q[2];
rz(2.8420908) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(1.6925156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2285948) q[0];
sx q[0];
rz(-1.2499502) q[0];
sx q[0];
rz(2.6477497) q[0];
rz(2.2834868) q[1];
sx q[1];
rz(-1.162642) q[1];
sx q[1];
rz(1.7869305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0070902) q[0];
sx q[0];
rz(-1.0225774) q[0];
sx q[0];
rz(-0.28315084) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0787215) q[2];
sx q[2];
rz(-1.4495696) q[2];
sx q[2];
rz(1.4722784) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87397447) q[1];
sx q[1];
rz(-0.81470352) q[1];
sx q[1];
rz(-1.4741298) q[1];
x q[2];
rz(0.51952298) q[3];
sx q[3];
rz(-0.96599993) q[3];
sx q[3];
rz(-2.7407569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4570423) q[2];
sx q[2];
rz(-2.5164618) q[2];
sx q[2];
rz(-0.54614145) q[2];
rz(2.9438733) q[3];
sx q[3];
rz(-2.1106014) q[3];
sx q[3];
rz(-0.24699591) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207183) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(2.0350631) q[0];
rz(2.5571892) q[1];
sx q[1];
rz(-1.1738651) q[1];
sx q[1];
rz(1.0027142) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2276777) q[0];
sx q[0];
rz(-1.4105689) q[0];
sx q[0];
rz(1.1315617) q[0];
rz(0.8876351) q[2];
sx q[2];
rz(-0.37345593) q[2];
sx q[2];
rz(-1.4591726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20368491) q[1];
sx q[1];
rz(-1.6069429) q[1];
sx q[1];
rz(-0.51236492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5224696) q[3];
sx q[3];
rz(-1.2803708) q[3];
sx q[3];
rz(-2.8324236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5661085) q[2];
sx q[2];
rz(-0.9245975) q[2];
sx q[2];
rz(-1.3893348) q[2];
rz(0.53650457) q[3];
sx q[3];
rz(-1.5080695) q[3];
sx q[3];
rz(-2.2947521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0480807) q[0];
sx q[0];
rz(-2.7219613) q[0];
sx q[0];
rz(2.804948) q[0];
rz(-3.0632784) q[1];
sx q[1];
rz(-1.9322194) q[1];
sx q[1];
rz(1.9858817) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2172359) q[0];
sx q[0];
rz(-1.2343533) q[0];
sx q[0];
rz(-0.84678679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1572273) q[2];
sx q[2];
rz(-1.6317657) q[2];
sx q[2];
rz(-0.79198906) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1286238) q[1];
sx q[1];
rz(-1.4358724) q[1];
sx q[1];
rz(1.0781481) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8248057) q[3];
sx q[3];
rz(-1.5419349) q[3];
sx q[3];
rz(-0.76735088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0597421) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(2.3183863) q[2];
rz(-2.1246223) q[3];
sx q[3];
rz(-0.40139324) q[3];
sx q[3];
rz(-3.1034191) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80242771) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(-1.7145351) q[0];
rz(1.6629705) q[1];
sx q[1];
rz(-2.069811) q[1];
sx q[1];
rz(-2.2907385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9768391) q[0];
sx q[0];
rz(-1.098159) q[0];
sx q[0];
rz(-0.12906277) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1673604) q[2];
sx q[2];
rz(-0.15359226) q[2];
sx q[2];
rz(3.0418015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6809594) q[1];
sx q[1];
rz(-0.52648989) q[1];
sx q[1];
rz(0.10641702) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5062648) q[3];
sx q[3];
rz(-1.7108265) q[3];
sx q[3];
rz(-2.9794995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6829546) q[2];
sx q[2];
rz(-1.7913603) q[2];
sx q[2];
rz(0.82009912) q[2];
rz(1.594515) q[3];
sx q[3];
rz(-1.4328512) q[3];
sx q[3];
rz(-1.1944176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9783258) q[0];
sx q[0];
rz(-1.6303202) q[0];
sx q[0];
rz(-1.6250961) q[0];
rz(0.77240472) q[1];
sx q[1];
rz(-0.62320566) q[1];
sx q[1];
rz(-2.1355245) q[1];
rz(1.3428691) q[2];
sx q[2];
rz(-0.89222903) q[2];
sx q[2];
rz(-1.0747838) q[2];
rz(2.0396654) q[3];
sx q[3];
rz(-2.0004604) q[3];
sx q[3];
rz(-0.018044005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
