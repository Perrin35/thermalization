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
rz(-0.57617968) q[0];
sx q[0];
rz(-1.4644858) q[0];
sx q[0];
rz(-0.65281868) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(5.6607487) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.884812) q[0];
sx q[0];
rz(-0.33397618) q[0];
sx q[0];
rz(1.2080753) q[0];
rz(0.60127778) q[2];
sx q[2];
rz(-1.2922505) q[2];
sx q[2];
rz(2.8097866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6316526) q[1];
sx q[1];
rz(-1.7544244) q[1];
sx q[1];
rz(1.0650403) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9838324) q[3];
sx q[3];
rz(-0.96786849) q[3];
sx q[3];
rz(-1.9229398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5114674) q[2];
sx q[2];
rz(-1.3518535) q[2];
sx q[2];
rz(0.54568616) q[2];
rz(0.68583471) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(-0.95664501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.086455258) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(2.0461244) q[0];
rz(2.144004) q[1];
sx q[1];
rz(-1.2350524) q[1];
sx q[1];
rz(1.197061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.469968) q[0];
sx q[0];
rz(-1.8238245) q[0];
sx q[0];
rz(-0.96970153) q[0];
x q[1];
rz(-0.95112339) q[2];
sx q[2];
rz(-2.24555) q[2];
sx q[2];
rz(-0.63979606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61651308) q[1];
sx q[1];
rz(-0.54821205) q[1];
sx q[1];
rz(-2.1384778) q[1];
x q[2];
rz(2.6359699) q[3];
sx q[3];
rz(-2.7693336) q[3];
sx q[3];
rz(2.3501808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7257488) q[2];
sx q[2];
rz(-1.3943322) q[2];
sx q[2];
rz(1.7572629) q[2];
rz(-1.3437126) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(1.869092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62190732) q[0];
sx q[0];
rz(-2.3540731) q[0];
sx q[0];
rz(0.4271048) q[0];
rz(1.2318132) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(2.5274091) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6525567) q[0];
sx q[0];
rz(-1.2430191) q[0];
sx q[0];
rz(1.2851738) q[0];
rz(2.3947507) q[2];
sx q[2];
rz(-1.7201506) q[2];
sx q[2];
rz(0.17490444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78313561) q[1];
sx q[1];
rz(-1.453521) q[1];
sx q[1];
rz(-3.0774119) q[1];
x q[2];
rz(-1.8133468) q[3];
sx q[3];
rz(-1.0772675) q[3];
sx q[3];
rz(-2.7427615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68982879) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(0.27099398) q[2];
rz(-0.48218918) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86376205) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(2.7226287) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.9326262) q[1];
sx q[1];
rz(3.0640501) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6200752) q[0];
sx q[0];
rz(-1.8635611) q[0];
sx q[0];
rz(-2.1174681) q[0];
rz(2.783804) q[2];
sx q[2];
rz(-2.3364106) q[2];
sx q[2];
rz(1.5186269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3540707) q[1];
sx q[1];
rz(-1.7131299) q[1];
sx q[1];
rz(2.700129) q[1];
rz(-pi) q[2];
rz(-0.61309149) q[3];
sx q[3];
rz(-0.35260751) q[3];
sx q[3];
rz(-2.4168454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0858687) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(-0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.5643352) q[0];
sx q[0];
rz(-0.99522796) q[0];
sx q[0];
rz(1.7876392) q[0];
rz(-2.1063781) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.530102) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28883067) q[0];
sx q[0];
rz(-1.0438598) q[0];
sx q[0];
rz(2.2599392) q[0];
x q[1];
rz(-2.8132854) q[2];
sx q[2];
rz(-1.7264778) q[2];
sx q[2];
rz(-2.1552483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6616231) q[1];
sx q[1];
rz(-1.1692059) q[1];
sx q[1];
rz(1.5127403) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0102481) q[3];
sx q[3];
rz(-2.7274611) q[3];
sx q[3];
rz(-1.9589748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80663854) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(0.063892603) q[2];
rz(2.4334) q[3];
sx q[3];
rz(-1.4539366) q[3];
sx q[3];
rz(-1.3341058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42234364) q[0];
sx q[0];
rz(-0.023119211) q[0];
sx q[0];
rz(-2.0656021) q[0];
rz(0.05323449) q[1];
sx q[1];
rz(-0.85317555) q[1];
sx q[1];
rz(-2.7630973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814893) q[0];
sx q[0];
rz(-1.3162587) q[0];
sx q[0];
rz(-1.7053563) q[0];
rz(1.1302275) q[2];
sx q[2];
rz(-0.097245596) q[2];
sx q[2];
rz(0.85376109) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15684756) q[1];
sx q[1];
rz(-1.3917006) q[1];
sx q[1];
rz(-0.26247929) q[1];
x q[2];
rz(0.58121292) q[3];
sx q[3];
rz(-0.85678116) q[3];
sx q[3];
rz(-2.9670144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9516912) q[2];
sx q[2];
rz(-1.6804164) q[2];
sx q[2];
rz(2.107673) q[2];
rz(2.2377491) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(3.097539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251625) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(0.59301162) q[0];
rz(-3.016839) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(-0.4932901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6007227) q[0];
sx q[0];
rz(-2.713795) q[0];
sx q[0];
rz(-0.85036253) q[0];
rz(-pi) q[1];
rz(-2.9134995) q[2];
sx q[2];
rz(-1.5465294) q[2];
sx q[2];
rz(-2.9920141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3354937) q[1];
sx q[1];
rz(-0.68365288) q[1];
sx q[1];
rz(0.70913507) q[1];
rz(-2.7760963) q[3];
sx q[3];
rz(-2.8474244) q[3];
sx q[3];
rz(-2.1902167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3844246) q[2];
sx q[2];
rz(-1.9809456) q[2];
sx q[2];
rz(1.9942795) q[2];
rz(2.3186963) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61854521) q[0];
sx q[0];
rz(-1.2724027) q[0];
sx q[0];
rz(0.4185032) q[0];
rz(-3.0283527) q[1];
sx q[1];
rz(-1.4694045) q[1];
sx q[1];
rz(-1.07771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200054) q[0];
sx q[0];
rz(-1.5527871) q[0];
sx q[0];
rz(-1.8593345) q[0];
rz(-pi) q[1];
rz(-0.32025614) q[2];
sx q[2];
rz(-2.3204707) q[2];
sx q[2];
rz(2.1964354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.019494836) q[1];
sx q[1];
rz(-1.6198525) q[1];
sx q[1];
rz(-8*pi/13) q[1];
rz(3.0642068) q[3];
sx q[3];
rz(-0.65801453) q[3];
sx q[3];
rz(-2.8946427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8354127) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(1.6861247) q[2];
rz(-2.9742187) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(1.0719871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-0.48285943) q[0];
sx q[0];
rz(-1.3128244) q[0];
sx q[0];
rz(-2.7400548) q[0];
rz(2.7032779) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(-1.4567136) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4493713) q[0];
sx q[0];
rz(-1.6047932) q[0];
sx q[0];
rz(3.1319989) q[0];
rz(-pi) q[1];
rz(2.8598665) q[2];
sx q[2];
rz(-0.67517074) q[2];
sx q[2];
rz(-3.0670241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46692013) q[1];
sx q[1];
rz(-0.90306696) q[1];
sx q[1];
rz(2.2064462) q[1];
x q[2];
rz(2.4026186) q[3];
sx q[3];
rz(-0.65118507) q[3];
sx q[3];
rz(0.84924752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7598286) q[2];
sx q[2];
rz(-2.9324053) q[2];
sx q[2];
rz(0.71167243) q[2];
rz(0.034218637) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(-1.155352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782368) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(2.3427298) q[0];
rz(1.0460188) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(0.23652133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4170161) q[0];
sx q[0];
rz(-1.5484705) q[0];
sx q[0];
rz(1.6315559) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4948214) q[2];
sx q[2];
rz(-1.0318709) q[2];
sx q[2];
rz(2.981271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1139086) q[1];
sx q[1];
rz(-0.56237537) q[1];
sx q[1];
rz(2.4706868) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13670758) q[3];
sx q[3];
rz(-2.4259704) q[3];
sx q[3];
rz(-0.8807883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.11655) q[2];
sx q[2];
rz(-0.5566842) q[2];
sx q[2];
rz(-0.64209783) q[2];
rz(-0.56946483) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(3.0552982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4702598) q[0];
sx q[0];
rz(-1.8076121) q[0];
sx q[0];
rz(-2.1847834) q[0];
rz(-1.9500465) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(-0.79193476) q[2];
sx q[2];
rz(-1.6583233) q[2];
sx q[2];
rz(0.59817373) q[2];
rz(-1.996205) q[3];
sx q[3];
rz(-1.26401) q[3];
sx q[3];
rz(-1.1372529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
