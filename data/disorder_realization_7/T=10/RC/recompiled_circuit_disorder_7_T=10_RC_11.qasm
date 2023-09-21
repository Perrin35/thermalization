OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(2.9404844) q[0];
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(2.5147658) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97323595) q[0];
sx q[0];
rz(-1.7587979) q[0];
sx q[0];
rz(-0.12113916) q[0];
x q[1];
rz(-0.35194273) q[2];
sx q[2];
rz(-1.8004187) q[2];
sx q[2];
rz(-2.4156092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1605692) q[1];
sx q[1];
rz(-1.261928) q[1];
sx q[1];
rz(2.1916051) q[1];
x q[2];
rz(0.17625531) q[3];
sx q[3];
rz(-1.3000543) q[3];
sx q[3];
rz(0.83812974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8132919) q[0];
sx q[0];
rz(-2.11587) q[0];
sx q[0];
rz(-2.6585447) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7254843) q[2];
sx q[2];
rz(-0.62972087) q[2];
sx q[2];
rz(-0.94905084) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0774035) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(2.2662152) q[1];
rz(-pi) q[2];
rz(0.38856216) q[3];
sx q[3];
rz(-1.7266759) q[3];
sx q[3];
rz(0.673783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(0.57998747) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-1.0823762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28856453) q[0];
sx q[0];
rz(-1.0596501) q[0];
sx q[0];
rz(-1.9938064) q[0];
rz(-pi) q[1];
rz(1.9252831) q[2];
sx q[2];
rz(-1.4282994) q[2];
sx q[2];
rz(0.49152495) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0116926) q[1];
sx q[1];
rz(-1.3492279) q[1];
sx q[1];
rz(-1.3264015) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.247418) q[3];
sx q[3];
rz(-1.7824899) q[3];
sx q[3];
rz(2.4993103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(1.6484377) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64788139) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(-0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(2.6534973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(3.0547633) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5627665) q[2];
sx q[2];
rz(-2.0818713) q[2];
sx q[2];
rz(1.5103112) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51302233) q[1];
sx q[1];
rz(-1.4161308) q[1];
sx q[1];
rz(3.0111074) q[1];
x q[2];
rz(-2.5591764) q[3];
sx q[3];
rz(-1.1629259) q[3];
sx q[3];
rz(1.2553314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(-0.50746894) q[2];
rz(-1.1821702) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5421211) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-2.7713293) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(2.945074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614583) q[0];
sx q[0];
rz(-0.41813865) q[0];
sx q[0];
rz(2.3925376) q[0];
x q[1];
rz(-3.1095805) q[2];
sx q[2];
rz(-2.0316761) q[2];
sx q[2];
rz(-2.8061158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49911753) q[1];
sx q[1];
rz(-2.827008) q[1];
sx q[1];
rz(-2.5423074) q[1];
x q[2];
rz(-0.27627857) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(-2.7460263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(-0.20203461) q[0];
rz(-1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0575384) q[0];
sx q[0];
rz(-1.6469427) q[0];
sx q[0];
rz(2.8094859) q[0];
x q[1];
rz(-3.0300573) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(2.6967449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(2.1761314) q[1];
rz(-1.1145669) q[3];
sx q[3];
rz(-0.8884512) q[3];
sx q[3];
rz(-1.8519459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(-2.1202309) q[0];
rz(2.0102603) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(2.9313415) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1649095) q[0];
sx q[0];
rz(-1.7047791) q[0];
sx q[0];
rz(3.0747736) q[0];
x q[1];
rz(-1.636593) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(-2.3492299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0381895) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(1.3693387) q[1];
rz(-0.93123575) q[3];
sx q[3];
rz(-0.82914549) q[3];
sx q[3];
rz(0.88013807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(-0.58724171) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(0.97672021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92102345) q[0];
sx q[0];
rz(-2.0221634) q[0];
sx q[0];
rz(1.736182) q[0];
x q[1];
rz(0.080539695) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(-1.5921519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50423065) q[1];
sx q[1];
rz(-2.0081875) q[1];
sx q[1];
rz(-0.38423844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55240734) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(-1.6747024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(0.89909536) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2313401) q[0];
sx q[0];
rz(-1.9950584) q[0];
sx q[0];
rz(-2.2146241) q[0];
rz(-pi) q[1];
rz(-2.1295206) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(0.65288359) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5324459) q[1];
sx q[1];
rz(-1.6858613) q[1];
sx q[1];
rz(-1.902641) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9409688) q[3];
sx q[3];
rz(-1.5705974) q[3];
sx q[3];
rz(1.1699333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1961394) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2755651) q[0];
sx q[0];
rz(-2.4336928) q[0];
sx q[0];
rz(0.78535725) q[0];
x q[1];
rz(-2.8357387) q[2];
sx q[2];
rz(-1.9216188) q[2];
sx q[2];
rz(2.2985814) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6798903) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(-3.0479792) q[1];
rz(-pi) q[2];
rz(-2.7361717) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(-2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(3.0449384) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(1.2697521) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-2.4772714) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(-1.312064) q[3];
sx q[3];
rz(-1.6862292) q[3];
sx q[3];
rz(-0.31173691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];