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
rz(0.12198099) q[0];
sx q[0];
rz(3.0966336) q[0];
sx q[0];
rz(9.3293204) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(3.0756693) q[1];
sx q[1];
rz(8.8311721) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5636162) q[0];
sx q[0];
rz(-1.2145743) q[0];
sx q[0];
rz(-2.6249159) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4727137) q[2];
sx q[2];
rz(-0.19193412) q[2];
sx q[2];
rz(2.8366983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2954457) q[1];
sx q[1];
rz(-0.6653924) q[1];
sx q[1];
rz(-2.800368) q[1];
x q[2];
rz(-0.27563654) q[3];
sx q[3];
rz(-2.079981) q[3];
sx q[3];
rz(-2.07723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.060071271) q[2];
sx q[2];
rz(-1.8895443) q[2];
sx q[2];
rz(2.579465) q[2];
rz(-2.8020322) q[3];
sx q[3];
rz(-1.5226676) q[3];
sx q[3];
rz(-1.08574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.958441) q[0];
sx q[0];
rz(-0.14587942) q[0];
sx q[0];
rz(1.9802144) q[0];
rz(-2.9005652) q[1];
sx q[1];
rz(-0.82254219) q[1];
sx q[1];
rz(2.8384812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68837092) q[0];
sx q[0];
rz(-1.3728598) q[0];
sx q[0];
rz(-3.0730547) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6952405) q[2];
sx q[2];
rz(-1.7156938) q[2];
sx q[2];
rz(3.0595487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54359879) q[1];
sx q[1];
rz(-2.139787) q[1];
sx q[1];
rz(1.8813716) q[1];
x q[2];
rz(-2.6637023) q[3];
sx q[3];
rz(-1.9955561) q[3];
sx q[3];
rz(2.6997944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3031528) q[2];
sx q[2];
rz(-1.8586321) q[2];
sx q[2];
rz(-0.34070936) q[2];
rz(-2.7029612) q[3];
sx q[3];
rz(-0.80629587) q[3];
sx q[3];
rz(3.0830234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4230147) q[0];
sx q[0];
rz(-2.5231762) q[0];
sx q[0];
rz(2.3037236) q[0];
rz(2.2184929) q[1];
sx q[1];
rz(-1.2078614) q[1];
sx q[1];
rz(-0.34742483) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9073931) q[0];
sx q[0];
rz(-0.1009909) q[0];
sx q[0];
rz(1.1924465) q[0];
rz(1.3143888) q[2];
sx q[2];
rz(-1.6056345) q[2];
sx q[2];
rz(-2.9098791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4975884) q[1];
sx q[1];
rz(-2.1198744) q[1];
sx q[1];
rz(-2.1357029) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33073552) q[3];
sx q[3];
rz(-1.3130672) q[3];
sx q[3];
rz(-2.7455519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.122494) q[2];
sx q[2];
rz(-0.118003) q[2];
sx q[2];
rz(0.72354358) q[2];
rz(-1.2875693) q[3];
sx q[3];
rz(-0.55793327) q[3];
sx q[3];
rz(-3.0961228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7751854) q[0];
sx q[0];
rz(-2.7025096) q[0];
sx q[0];
rz(1.084569) q[0];
rz(1.1868125) q[1];
sx q[1];
rz(-1.2326406) q[1];
sx q[1];
rz(-1.4694227) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5795008) q[0];
sx q[0];
rz(-3.1070624) q[0];
sx q[0];
rz(0.51106913) q[0];
rz(-pi) q[1];
rz(-0.78408636) q[2];
sx q[2];
rz(-2.6166281) q[2];
sx q[2];
rz(-0.80277473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96257229) q[1];
sx q[1];
rz(-2.031759) q[1];
sx q[1];
rz(2.712972) q[1];
rz(-pi) q[2];
rz(1.2354492) q[3];
sx q[3];
rz(-1.2307388) q[3];
sx q[3];
rz(0.53676134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12545086) q[2];
sx q[2];
rz(-0.90105385) q[2];
sx q[2];
rz(-0.93552843) q[2];
rz(-2.9901796) q[3];
sx q[3];
rz(-0.97984034) q[3];
sx q[3];
rz(-0.51101959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147515) q[0];
sx q[0];
rz(-1.5837357) q[0];
sx q[0];
rz(2.3874808) q[0];
rz(0.61112815) q[1];
sx q[1];
rz(-1.1437623) q[1];
sx q[1];
rz(0.67620826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0128764) q[0];
sx q[0];
rz(-0.36660796) q[0];
sx q[0];
rz(2.0373197) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1235672) q[2];
sx q[2];
rz(-1.8844731) q[2];
sx q[2];
rz(2.3540135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83482997) q[1];
sx q[1];
rz(-1.3993629) q[1];
sx q[1];
rz(-1.7333561) q[1];
rz(2.5557983) q[3];
sx q[3];
rz(-1.5021281) q[3];
sx q[3];
rz(-0.60528994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0040697441) q[2];
sx q[2];
rz(-0.71228945) q[2];
sx q[2];
rz(-2.3884657) q[2];
rz(-0.5641886) q[3];
sx q[3];
rz(-2.0337532) q[3];
sx q[3];
rz(2.3851725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.18993987) q[0];
sx q[0];
rz(-2.6615182) q[0];
sx q[0];
rz(-0.22115627) q[0];
rz(-3.0033374) q[1];
sx q[1];
rz(-1.4555376) q[1];
sx q[1];
rz(2.5855605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1336687) q[0];
sx q[0];
rz(-1.8147665) q[0];
sx q[0];
rz(1.3707815) q[0];
x q[1];
rz(-0.37555917) q[2];
sx q[2];
rz(-1.9347768) q[2];
sx q[2];
rz(-3.08422) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8670571) q[1];
sx q[1];
rz(-2.9867771) q[1];
sx q[1];
rz(-3.0068924) q[1];
x q[2];
rz(-3.0228244) q[3];
sx q[3];
rz(-1.1221544) q[3];
sx q[3];
rz(1.4172163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14744645) q[2];
sx q[2];
rz(-2.3799956) q[2];
sx q[2];
rz(2.2046294) q[2];
rz(-2.7312036) q[3];
sx q[3];
rz(-0.97987163) q[3];
sx q[3];
rz(0.65569896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74528247) q[0];
sx q[0];
rz(-0.87438011) q[0];
sx q[0];
rz(-0.96258798) q[0];
rz(0.67086041) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(-0.76505351) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017185171) q[0];
sx q[0];
rz(-1.9694424) q[0];
sx q[0];
rz(0.17332698) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1704438) q[2];
sx q[2];
rz(-3.0191688) q[2];
sx q[2];
rz(2.3442307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9817762) q[1];
sx q[1];
rz(-2.4764531) q[1];
sx q[1];
rz(0.90250166) q[1];
rz(1.9387127) q[3];
sx q[3];
rz(-1.4548317) q[3];
sx q[3];
rz(-0.74527568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44136167) q[2];
sx q[2];
rz(-0.92741489) q[2];
sx q[2];
rz(-2.3564763) q[2];
rz(2.650812) q[3];
sx q[3];
rz(-1.1622585) q[3];
sx q[3];
rz(2.6665915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15071507) q[0];
sx q[0];
rz(-3.0325723) q[0];
sx q[0];
rz(2.9955067) q[0];
rz(2.8699005) q[1];
sx q[1];
rz(-0.79333317) q[1];
sx q[1];
rz(0.21154107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1871227) q[0];
sx q[0];
rz(-2.3516708) q[0];
sx q[0];
rz(0.23440897) q[0];
x q[1];
rz(3.0258133) q[2];
sx q[2];
rz(-1.2010788) q[2];
sx q[2];
rz(-3.0860975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.817343) q[1];
sx q[1];
rz(-1.6403664) q[1];
sx q[1];
rz(2.9051498) q[1];
rz(-1.7320427) q[3];
sx q[3];
rz(-0.60617709) q[3];
sx q[3];
rz(2.2796749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81296706) q[2];
sx q[2];
rz(-2.2058637) q[2];
sx q[2];
rz(2.3198371) q[2];
rz(-1.3385319) q[3];
sx q[3];
rz(-2.0526363) q[3];
sx q[3];
rz(-0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61447918) q[0];
sx q[0];
rz(-0.69632691) q[0];
sx q[0];
rz(-1.3592199) q[0];
rz(-0.95787734) q[1];
sx q[1];
rz(-0.33708894) q[1];
sx q[1];
rz(-2.8639796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2569297) q[0];
sx q[0];
rz(-2.1799934) q[0];
sx q[0];
rz(2.4692187) q[0];
rz(-0.77582716) q[2];
sx q[2];
rz(-1.3444573) q[2];
sx q[2];
rz(-3.0507342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23944868) q[1];
sx q[1];
rz(-2.7056597) q[1];
sx q[1];
rz(0.3050584) q[1];
rz(-pi) q[2];
rz(-2.7294159) q[3];
sx q[3];
rz(-2.7320903) q[3];
sx q[3];
rz(3.1016027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9684888) q[2];
sx q[2];
rz(-2.5471881) q[2];
sx q[2];
rz(1.28432) q[2];
rz(-2.3641018) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(-0.29977453) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023329968) q[0];
sx q[0];
rz(-1.4813923) q[0];
sx q[0];
rz(0.0059286038) q[0];
rz(-1.3395576) q[1];
sx q[1];
rz(-2.1052994) q[1];
sx q[1];
rz(-0.48066995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67424947) q[0];
sx q[0];
rz(-3.1200069) q[0];
sx q[0];
rz(-2.4186446) q[0];
rz(-2.1176012) q[2];
sx q[2];
rz(-2.4158187) q[2];
sx q[2];
rz(0.99601907) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0152032) q[1];
sx q[1];
rz(-1.9815832) q[1];
sx q[1];
rz(-2.4900764) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3044429) q[3];
sx q[3];
rz(-1.7031809) q[3];
sx q[3];
rz(1.3416106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9504451) q[2];
sx q[2];
rz(-2.5823249) q[2];
sx q[2];
rz(0.17678235) q[2];
rz(-0.92987972) q[3];
sx q[3];
rz(-1.1888489) q[3];
sx q[3];
rz(1.0298347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69242351) q[0];
sx q[0];
rz(-1.7208736) q[0];
sx q[0];
rz(2.0509913) q[0];
rz(1.0724267) q[1];
sx q[1];
rz(-1.4785531) q[1];
sx q[1];
rz(-1.4428152) q[1];
rz(1.5771796) q[2];
sx q[2];
rz(-2.0626358) q[2];
sx q[2];
rz(-0.067968702) q[2];
rz(1.4789875) q[3];
sx q[3];
rz(-2.1740365) q[3];
sx q[3];
rz(-2.4400644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
