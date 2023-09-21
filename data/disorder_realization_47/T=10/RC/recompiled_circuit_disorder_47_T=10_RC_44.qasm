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
rz(-2.2934409) q[1];
sx q[1];
rz(3.1037722) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8671994) q[0];
sx q[0];
rz(-1.6490071) q[0];
sx q[0];
rz(-1.7645416) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0231789) q[2];
sx q[2];
rz(-0.9127494) q[2];
sx q[2];
rz(-0.40068914) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5303858) q[1];
sx q[1];
rz(-1.621765) q[1];
sx q[1];
rz(-0.74688046) q[1];
rz(-pi) q[2];
x q[2];
rz(1.171265) q[3];
sx q[3];
rz(-2.7683308) q[3];
sx q[3];
rz(-1.319862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7906856) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(-0.12617271) q[3];
sx q[3];
rz(-1.7254555) q[3];
sx q[3];
rz(-3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(-1.3417473) q[1];
sx q[1];
rz(-1.7780875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7994021) q[0];
sx q[0];
rz(-1.0418833) q[0];
sx q[0];
rz(2.0687813) q[0];
x q[1];
rz(-0.5559276) q[2];
sx q[2];
rz(-1.5276507) q[2];
sx q[2];
rz(2.3748929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.99216695) q[1];
sx q[1];
rz(-2.4352695) q[1];
sx q[1];
rz(-2.7041433) q[1];
rz(-pi) q[2];
rz(2.2927106) q[3];
sx q[3];
rz(-2.1697681) q[3];
sx q[3];
rz(-2.2357383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3229225) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(2.4070516) q[2];
rz(-2.5143886) q[3];
sx q[3];
rz(-1.6278798) q[3];
sx q[3];
rz(0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(0.27221361) q[0];
rz(-2.294337) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(-2.8289657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011742) q[0];
sx q[0];
rz(-0.22589382) q[0];
sx q[0];
rz(-0.24143879) q[0];
rz(-pi) q[1];
rz(-2.4467144) q[2];
sx q[2];
rz(-2.4701397) q[2];
sx q[2];
rz(-1.7365255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4362267) q[1];
sx q[1];
rz(-1.6402906) q[1];
sx q[1];
rz(1.8314929) q[1];
x q[2];
rz(-1.1945046) q[3];
sx q[3];
rz(-0.84405758) q[3];
sx q[3];
rz(-2.4220667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12038885) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(2.1728544) q[3];
sx q[3];
rz(-2.0961943) q[3];
sx q[3];
rz(0.43280861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75354904) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(-2.5182305) q[0];
rz(2.3240044) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(-3.0923016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.277963) q[0];
sx q[0];
rz(-2.3290714) q[0];
sx q[0];
rz(-0.066594007) q[0];
x q[1];
rz(0.85929112) q[2];
sx q[2];
rz(-1.3261194) q[2];
sx q[2];
rz(1.7161075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4781487) q[1];
sx q[1];
rz(-1.2423007) q[1];
sx q[1];
rz(-0.25681396) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0443346) q[3];
sx q[3];
rz(-1.1757848) q[3];
sx q[3];
rz(-2.928424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-2.1566186) q[2];
rz(0.90304053) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(0.27339098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7966998) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(0.087619089) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(-2.2706251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5721711) q[0];
sx q[0];
rz(-1.8580125) q[0];
sx q[0];
rz(1.6006908) q[0];
rz(-pi) q[1];
rz(-1.9268553) q[2];
sx q[2];
rz(-2.7241754) q[2];
sx q[2];
rz(0.88127121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4082143) q[1];
sx q[1];
rz(-0.92792643) q[1];
sx q[1];
rz(-1.9233568) q[1];
x q[2];
rz(1.3887651) q[3];
sx q[3];
rz(-0.77266274) q[3];
sx q[3];
rz(0.80881892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0481723) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(-2.7887153) q[2];
rz(-2.2757163) q[3];
sx q[3];
rz(-0.70169774) q[3];
sx q[3];
rz(-2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46309328) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(-0.10678664) q[0];
rz(-1.9550025) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(2.1170763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56652503) q[0];
sx q[0];
rz(-1.1199513) q[0];
sx q[0];
rz(-1.8917811) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5365731) q[2];
sx q[2];
rz(-1.9540678) q[2];
sx q[2];
rz(-2.5503416) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34895996) q[1];
sx q[1];
rz(-1.0981202) q[1];
sx q[1];
rz(-0.93797586) q[1];
rz(-pi) q[2];
rz(-2.5173353) q[3];
sx q[3];
rz(-0.58371021) q[3];
sx q[3];
rz(-0.39238413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.391905) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(-0.8708896) q[2];
rz(-1.332256) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9695327) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-0.55091888) q[0];
rz(-2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-2.904772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2315002) q[0];
sx q[0];
rz(-1.5331368) q[0];
sx q[0];
rz(1.670027) q[0];
rz(-pi) q[1];
rz(1.7083488) q[2];
sx q[2];
rz(-1.4220861) q[2];
sx q[2];
rz(2.6580722) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6432861) q[1];
sx q[1];
rz(-2.8611538) q[1];
sx q[1];
rz(-1.8449057) q[1];
rz(0.83644609) q[3];
sx q[3];
rz(-1.340938) q[3];
sx q[3];
rz(-0.072007192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6992496) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(-2.701475) q[2];
rz(2.0464499) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(1.346689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(-0.029504689) q[0];
rz(-0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(0.11880076) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58177452) q[0];
sx q[0];
rz(-1.2604598) q[0];
sx q[0];
rz(1.6161726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45784874) q[2];
sx q[2];
rz(-2.0689788) q[2];
sx q[2];
rz(2.904441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58662215) q[1];
sx q[1];
rz(-1.8645617) q[1];
sx q[1];
rz(-0.48420669) q[1];
rz(1.2715862) q[3];
sx q[3];
rz(-0.51338235) q[3];
sx q[3];
rz(-1.3083003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7416731) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(-1.5681533) q[2];
rz(-1.167477) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(1.2742111) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2385999) q[0];
sx q[0];
rz(-0.82505834) q[0];
sx q[0];
rz(-0.15326823) q[0];
rz(-1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.1538039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.518084) q[0];
sx q[0];
rz(-0.46184807) q[0];
sx q[0];
rz(-0.64103809) q[0];
rz(2.2527163) q[2];
sx q[2];
rz(-2.0098915) q[2];
sx q[2];
rz(0.097749226) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8968618) q[1];
sx q[1];
rz(-1.3618999) q[1];
sx q[1];
rz(-0.056785866) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0204569) q[3];
sx q[3];
rz(-0.88705685) q[3];
sx q[3];
rz(2.5072806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29685059) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(-1.6142169) q[2];
rz(-2.1447003) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(2.2629288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24755724) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(-0.19009185) q[0];
rz(0.67063531) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(1.6533096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22083902) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(-0.12551813) q[0];
rz(-pi) q[1];
rz(-0.52942099) q[2];
sx q[2];
rz(-2.6396857) q[2];
sx q[2];
rz(-1.6255962) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59939811) q[1];
sx q[1];
rz(-2.290526) q[1];
sx q[1];
rz(-0.69516121) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87741239) q[3];
sx q[3];
rz(-2.0392232) q[3];
sx q[3];
rz(1.7978158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1378479) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(0.89912644) q[2];
rz(2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(-2.9639444) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-2.4317447) q[2];
sx q[2];
rz(-1.3550497) q[2];
sx q[2];
rz(1.4646127) q[2];
rz(1.9990986) q[3];
sx q[3];
rz(-2.8430568) q[3];
sx q[3];
rz(-1.6381016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];