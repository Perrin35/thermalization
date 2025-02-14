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
rz(2.6707668) q[0];
sx q[0];
rz(-0.45007053) q[0];
sx q[0];
rz(-0.39029628) q[0];
rz(-2.9974239) q[1];
sx q[1];
rz(-1.6427957) q[1];
sx q[1];
rz(-2.1013451) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9269104) q[0];
sx q[0];
rz(-1.2217064) q[0];
sx q[0];
rz(-0.30544282) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39019806) q[2];
sx q[2];
rz(-2.7510561) q[2];
sx q[2];
rz(-0.85987505) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8017218) q[1];
sx q[1];
rz(-2.2328908) q[1];
sx q[1];
rz(2.8784196) q[1];
x q[2];
rz(-1.5055429) q[3];
sx q[3];
rz(-2.3480519) q[3];
sx q[3];
rz(-0.82543594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1990004) q[2];
sx q[2];
rz(-0.60574836) q[2];
sx q[2];
rz(-1.595363) q[2];
rz(2.6664901) q[3];
sx q[3];
rz(-2.4964156) q[3];
sx q[3];
rz(1.1554385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5763181) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(0.43462547) q[0];
rz(-2.077153) q[1];
sx q[1];
rz(-0.39355215) q[1];
sx q[1];
rz(-2.2501066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28434174) q[0];
sx q[0];
rz(-1.9280757) q[0];
sx q[0];
rz(2.7140359) q[0];
rz(-2.5847748) q[2];
sx q[2];
rz(-0.80792431) q[2];
sx q[2];
rz(-2.9064532) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.621755) q[1];
sx q[1];
rz(-0.72417604) q[1];
sx q[1];
rz(-2.4228079) q[1];
rz(-pi) q[2];
rz(-2.2300612) q[3];
sx q[3];
rz(-1.5943822) q[3];
sx q[3];
rz(3.0167836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6330304) q[2];
sx q[2];
rz(-2.5789301) q[2];
sx q[2];
rz(2.0821345) q[2];
rz(-1.9153473) q[3];
sx q[3];
rz(-1.6112593) q[3];
sx q[3];
rz(-3.0219769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0934963) q[0];
sx q[0];
rz(-2.7432848) q[0];
sx q[0];
rz(-0.68921971) q[0];
rz(2.7293909) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(-1.162792) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6434518) q[0];
sx q[0];
rz(-1.2664794) q[0];
sx q[0];
rz(0.16022072) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1751129) q[2];
sx q[2];
rz(-0.89077026) q[2];
sx q[2];
rz(0.71875188) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4913669) q[1];
sx q[1];
rz(-2.0511481) q[1];
sx q[1];
rz(-0.65746376) q[1];
rz(-pi) q[2];
rz(-0.27966313) q[3];
sx q[3];
rz(-1.2919909) q[3];
sx q[3];
rz(1.6300843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73843655) q[2];
sx q[2];
rz(-1.4102035) q[2];
sx q[2];
rz(-0.23846826) q[2];
rz(0.23369914) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(-0.15271798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9018263) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(2.8357586) q[0];
rz(-0.26890525) q[1];
sx q[1];
rz(-0.79335672) q[1];
sx q[1];
rz(0.3522402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21111317) q[0];
sx q[0];
rz(-1.5310627) q[0];
sx q[0];
rz(-1.4369523) q[0];
x q[1];
rz(2.1969123) q[2];
sx q[2];
rz(-1.4123823) q[2];
sx q[2];
rz(-2.2974654) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8246987) q[1];
sx q[1];
rz(-0.76008633) q[1];
sx q[1];
rz(2.555067) q[1];
rz(-2.8958578) q[3];
sx q[3];
rz(-0.80348368) q[3];
sx q[3];
rz(-1.6197698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29117808) q[2];
sx q[2];
rz(-1.6213657) q[2];
sx q[2];
rz(-0.16695437) q[2];
rz(-3.0252365) q[3];
sx q[3];
rz(-2.6226624) q[3];
sx q[3];
rz(1.2701344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6312234) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(-1.2315617) q[1];
sx q[1];
rz(-1.4304588) q[1];
sx q[1];
rz(-0.11428741) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60576263) q[0];
sx q[0];
rz(-1.0530942) q[0];
sx q[0];
rz(2.1735031) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8631526) q[2];
sx q[2];
rz(-1.2781218) q[2];
sx q[2];
rz(-3.1029683) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3975349) q[1];
sx q[1];
rz(-1.1611132) q[1];
sx q[1];
rz(2.4593705) q[1];
rz(3.1350713) q[3];
sx q[3];
rz(-2.4355585) q[3];
sx q[3];
rz(-1.8268585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0564702) q[2];
sx q[2];
rz(-1.2135999) q[2];
sx q[2];
rz(-1.4465793) q[2];
rz(0.96212402) q[3];
sx q[3];
rz(-2.646793) q[3];
sx q[3];
rz(-1.3303293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98244786) q[0];
sx q[0];
rz(-1.5626937) q[0];
sx q[0];
rz(-2.6084117) q[0];
rz(1.5348596) q[1];
sx q[1];
rz(-2.3695562) q[1];
sx q[1];
rz(2.3981222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3822298) q[0];
sx q[0];
rz(-2.8385525) q[0];
sx q[0];
rz(-2.8735301) q[0];
x q[1];
rz(1.942988) q[2];
sx q[2];
rz(-2.5649338) q[2];
sx q[2];
rz(1.3194336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20821321) q[1];
sx q[1];
rz(-2.4543896) q[1];
sx q[1];
rz(-0.40172462) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0251158) q[3];
sx q[3];
rz(-2.4719596) q[3];
sx q[3];
rz(-0.39519924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47034904) q[2];
sx q[2];
rz(-0.84731421) q[2];
sx q[2];
rz(-0.49873763) q[2];
rz(2.4399452) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(-2.4832895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43211234) q[0];
sx q[0];
rz(-1.9146336) q[0];
sx q[0];
rz(3.0416601) q[0];
rz(-1.2808895) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(2.4023712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4866132) q[0];
sx q[0];
rz(-0.57679048) q[0];
sx q[0];
rz(2.3705215) q[0];
rz(1.5042502) q[2];
sx q[2];
rz(-0.4302667) q[2];
sx q[2];
rz(-0.86003785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2249965) q[1];
sx q[1];
rz(-2.1526622) q[1];
sx q[1];
rz(0.47486979) q[1];
x q[2];
rz(1.5600863) q[3];
sx q[3];
rz(-2.3751343) q[3];
sx q[3];
rz(-0.86957726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40333906) q[2];
sx q[2];
rz(-2.5103939) q[2];
sx q[2];
rz(1.6048019) q[2];
rz(-1.4295476) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(-0.79609377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0159863) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(-2.0269537) q[0];
rz(0.58427018) q[1];
sx q[1];
rz(-2.22157) q[1];
sx q[1];
rz(3.027473) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3875562) q[0];
sx q[0];
rz(-0.68337959) q[0];
sx q[0];
rz(2.6793733) q[0];
rz(-pi) q[1];
rz(2.1340964) q[2];
sx q[2];
rz(-0.85342583) q[2];
sx q[2];
rz(-0.12441758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89592592) q[1];
sx q[1];
rz(-1.9787496) q[1];
sx q[1];
rz(2.308564) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8382254) q[3];
sx q[3];
rz(-1.4576313) q[3];
sx q[3];
rz(-2.6306977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4383661) q[2];
sx q[2];
rz(-2.4195636) q[2];
sx q[2];
rz(-0.88877338) q[2];
rz(-1.599954) q[3];
sx q[3];
rz(-1.2069353) q[3];
sx q[3];
rz(2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090076598) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(-0.28710452) q[0];
rz(-1.9376532) q[1];
sx q[1];
rz(-2.2724889) q[1];
sx q[1];
rz(1.513011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0368691) q[0];
sx q[0];
rz(-1.6469886) q[0];
sx q[0];
rz(-2.8073187) q[0];
rz(-pi) q[1];
rz(1.0059297) q[2];
sx q[2];
rz(-0.73137368) q[2];
sx q[2];
rz(1.5648901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87344681) q[1];
sx q[1];
rz(-1.6153187) q[1];
sx q[1];
rz(2.4067307) q[1];
rz(1.318526) q[3];
sx q[3];
rz(-3.1384094) q[3];
sx q[3];
rz(1.1109655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.63605753) q[2];
sx q[2];
rz(-1.9745741) q[2];
sx q[2];
rz(0.30910811) q[2];
rz(1.2067893) q[3];
sx q[3];
rz(-2.759582) q[3];
sx q[3];
rz(2.8822854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1768271) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(2.4859909) q[0];
rz(-0.62878311) q[1];
sx q[1];
rz(-1.7318334) q[1];
sx q[1];
rz(1.7579196) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1882598) q[0];
sx q[0];
rz(-1.5273222) q[0];
sx q[0];
rz(-3.093958) q[0];
x q[1];
rz(-0.83413275) q[2];
sx q[2];
rz(-1.9185432) q[2];
sx q[2];
rz(-1.8054838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.261499) q[1];
sx q[1];
rz(-1.4215018) q[1];
sx q[1];
rz(-1.5890676) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8877533) q[3];
sx q[3];
rz(-0.51341265) q[3];
sx q[3];
rz(-1.9771345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3722966) q[2];
sx q[2];
rz(-1.6326222) q[2];
sx q[2];
rz(-0.24202913) q[2];
rz(0.58487839) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(-0.52403319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723421) q[0];
sx q[0];
rz(-0.62234288) q[0];
sx q[0];
rz(-0.55707669) q[0];
rz(-0.88049018) q[1];
sx q[1];
rz(-1.4028032) q[1];
sx q[1];
rz(1.0407851) q[1];
rz(-3.1351719) q[2];
sx q[2];
rz(-1.3646135) q[2];
sx q[2];
rz(-2.7529181) q[2];
rz(-1.7193033) q[3];
sx q[3];
rz(-2.5703493) q[3];
sx q[3];
rz(-0.66019365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
