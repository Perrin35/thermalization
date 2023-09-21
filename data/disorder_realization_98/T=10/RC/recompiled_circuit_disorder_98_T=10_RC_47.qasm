OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(-0.87061849) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(2.9918848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(2.2358405) q[0];
x q[1];
rz(0.82586536) q[2];
sx q[2];
rz(-2.4561433) q[2];
sx q[2];
rz(-0.65537383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61923164) q[1];
sx q[1];
rz(-1.3904966) q[1];
sx q[1];
rz(-0.038832263) q[1];
x q[2];
rz(-2.9505694) q[3];
sx q[3];
rz(-2.1130307) q[3];
sx q[3];
rz(-2.1477826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(0.40701436) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23806439) q[2];
sx q[2];
rz(-1.2859584) q[2];
sx q[2];
rz(-0.1711947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0374239) q[1];
sx q[1];
rz(-1.8579322) q[1];
sx q[1];
rz(-1.8038521) q[1];
rz(-pi) q[2];
rz(-0.84516256) q[3];
sx q[3];
rz(-2.3361428) q[3];
sx q[3];
rz(-2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(-1.7920378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91857598) q[0];
sx q[0];
rz(-0.71764676) q[0];
sx q[0];
rz(-0.36803228) q[0];
rz(0.037319855) q[2];
sx q[2];
rz(-1.5048426) q[2];
sx q[2];
rz(-2.0685591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6408491) q[1];
sx q[1];
rz(-1.8060246) q[1];
sx q[1];
rz(2.6231223) q[1];
rz(-pi) q[2];
rz(1.5924256) q[3];
sx q[3];
rz(-2.7154185) q[3];
sx q[3];
rz(-2.222995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-0.90399495) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(-0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-0.036380336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43544337) q[0];
sx q[0];
rz(-0.94067803) q[0];
sx q[0];
rz(-0.26421996) q[0];
rz(-2.5600299) q[2];
sx q[2];
rz(-1.2956603) q[2];
sx q[2];
rz(-3.0951701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42602793) q[1];
sx q[1];
rz(-1.702311) q[1];
sx q[1];
rz(0.87042602) q[1];
x q[2];
rz(0.39156885) q[3];
sx q[3];
rz(-2.8885926) q[3];
sx q[3];
rz(0.60245017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(-2.4711117) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(2.9072445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097965) q[0];
sx q[0];
rz(-2.095982) q[0];
sx q[0];
rz(-1.964142) q[0];
x q[1];
rz(-1.0767897) q[2];
sx q[2];
rz(-1.6550118) q[2];
sx q[2];
rz(0.46352026) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4489331) q[1];
sx q[1];
rz(-1.7198791) q[1];
sx q[1];
rz(2.3922608) q[1];
rz(0.078951051) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-0.22053545) q[2];
rz(0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.1436499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7276579) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(-1.6137705) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4684858) q[2];
sx q[2];
rz(-1.3654725) q[2];
sx q[2];
rz(0.38773195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5864582) q[1];
sx q[1];
rz(-1.1291593) q[1];
sx q[1];
rz(-3.0416136) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.850071) q[3];
sx q[3];
rz(-1.0696628) q[3];
sx q[3];
rz(2.6078893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75366655) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(-1.2333262) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9702643) q[0];
sx q[0];
rz(-1.5818705) q[0];
sx q[0];
rz(-1.1867255) q[0];
rz(-pi) q[1];
rz(2.5160518) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(0.90460888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26112939) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(0.1111828) q[1];
x q[2];
rz(-0.60621467) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-1.0160149) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163517) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.4455147) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381666) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(2.0057136) q[0];
x q[1];
rz(0.9332946) q[2];
sx q[2];
rz(-2.1142695) q[2];
sx q[2];
rz(-0.18804929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82842365) q[1];
sx q[1];
rz(-2.1143882) q[1];
sx q[1];
rz(1.6924752) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85556742) q[3];
sx q[3];
rz(-1.1547525) q[3];
sx q[3];
rz(-2.1895529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(-0.02773157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8160307) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(0.82657878) q[0];
rz(-0.8717732) q[2];
sx q[2];
rz(-1.36424) q[2];
sx q[2];
rz(2.2733462) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30222826) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(2.539413) q[1];
x q[2];
rz(1.5042217) q[3];
sx q[3];
rz(-1.5449617) q[3];
sx q[3];
rz(1.0915826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1256844) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.998418) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89850241) q[0];
sx q[0];
rz(-1.9542964) q[0];
sx q[0];
rz(2.6570508) q[0];
rz(2.5970039) q[2];
sx q[2];
rz(-1.1368183) q[2];
sx q[2];
rz(-1.4101654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.232302) q[1];
sx q[1];
rz(-1.5850987) q[1];
sx q[1];
rz(2.4453352) q[1];
rz(-pi) q[2];
rz(-2.8994843) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(-0.22656245) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(-1.9526341) q[2];
sx q[2];
rz(-1.1321862) q[2];
sx q[2];
rz(1.95375) q[2];
rz(-0.26279454) q[3];
sx q[3];
rz(-1.1017208) q[3];
sx q[3];
rz(2.8337939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
