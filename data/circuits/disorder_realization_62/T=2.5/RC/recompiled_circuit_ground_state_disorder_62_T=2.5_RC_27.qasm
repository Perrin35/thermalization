OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.59392053) q[0];
sx q[0];
rz(-2.5424818) q[0];
sx q[0];
rz(0.17679086) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(-2.63201) q[1];
sx q[1];
rz(-0.086960763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158543) q[0];
sx q[0];
rz(-0.70924067) q[0];
sx q[0];
rz(-2.5410557) q[0];
rz(-1.7638788) q[2];
sx q[2];
rz(-2.1336485) q[2];
sx q[2];
rz(2.4342997) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1198169) q[1];
sx q[1];
rz(-1.4939185) q[1];
sx q[1];
rz(-1.6989811) q[1];
rz(2.7417408) q[3];
sx q[3];
rz(-2.4824871) q[3];
sx q[3];
rz(2.6100998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1842492) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-0.59172612) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-0.44132909) q[3];
sx q[3];
rz(0.87619585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1205207) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(1.0907809) q[0];
rz(-2.4202994) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(-0.24478197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0140959) q[0];
sx q[0];
rz(-1.2359706) q[0];
sx q[0];
rz(-1.1017407) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82142395) q[2];
sx q[2];
rz(-1.9000179) q[2];
sx q[2];
rz(3.1067348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.822387) q[1];
sx q[1];
rz(-2.165395) q[1];
sx q[1];
rz(-2.2317504) q[1];
rz(-pi) q[2];
rz(-2.8448288) q[3];
sx q[3];
rz(-1.7113842) q[3];
sx q[3];
rz(2.1641017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.327534) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(-0.84612334) q[2];
rz(-0.36492473) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(-0.97682166) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1035128) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(0.14347759) q[0];
rz(1.6733276) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(-0.066468261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9003948) q[0];
sx q[0];
rz(-1.4816435) q[0];
sx q[0];
rz(2.4597079) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9793545) q[2];
sx q[2];
rz(-0.58149946) q[2];
sx q[2];
rz(-1.7844019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4926949) q[1];
sx q[1];
rz(-1.6872511) q[1];
sx q[1];
rz(1.6687188) q[1];
rz(1.7371561) q[3];
sx q[3];
rz(-0.8373973) q[3];
sx q[3];
rz(0.93033965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1116144) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(0.31271333) q[2];
rz(-2.8800268) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(-2.3436782) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12544352) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(3.0545767) q[0];
rz(-1.7866987) q[1];
sx q[1];
rz(-1.5785297) q[1];
sx q[1];
rz(-0.048390128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78260471) q[0];
sx q[0];
rz(-1.6986956) q[0];
sx q[0];
rz(-0.77632287) q[0];
x q[1];
rz(2.9913285) q[2];
sx q[2];
rz(-1.74904) q[2];
sx q[2];
rz(1.5706289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.684206) q[1];
sx q[1];
rz(-2.4438639) q[1];
sx q[1];
rz(-3.0569424) q[1];
x q[2];
rz(0.88416962) q[3];
sx q[3];
rz(-2.3034952) q[3];
sx q[3];
rz(-2.9626486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71063572) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(0.44130138) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(1.8080447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271673) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(-0.72702485) q[0];
rz(-1.6983039) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(-1.4220994) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2724394) q[0];
sx q[0];
rz(-2.5792988) q[0];
sx q[0];
rz(2.8773688) q[0];
rz(-pi) q[1];
rz(0.24025323) q[2];
sx q[2];
rz(-1.7723871) q[2];
sx q[2];
rz(0.83490419) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.369097) q[1];
sx q[1];
rz(-2.9926692) q[1];
sx q[1];
rz(-2.5293674) q[1];
x q[2];
rz(0.87559132) q[3];
sx q[3];
rz(-0.79447047) q[3];
sx q[3];
rz(-0.21847413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9299341) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(1.5809853) q[2];
rz(1.3382781) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(-2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018589858) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(0.71989584) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(-2.8954411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1937597) q[0];
sx q[0];
rz(-1.594055) q[0];
sx q[0];
rz(-2.860753) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75743875) q[2];
sx q[2];
rz(-1.2691109) q[2];
sx q[2];
rz(-1.2011004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0689905) q[1];
sx q[1];
rz(-1.700987) q[1];
sx q[1];
rz(-2.4058002) q[1];
x q[2];
rz(2.6443374) q[3];
sx q[3];
rz(-0.75102931) q[3];
sx q[3];
rz(0.3270783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64421946) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(-0.79968828) q[2];
rz(-0.52404809) q[3];
sx q[3];
rz(-2.7553813) q[3];
sx q[3];
rz(0.016949765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2972357) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(2.2204087) q[0];
rz(1.720403) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(2.1465819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5000347) q[0];
sx q[0];
rz(-1.3480268) q[0];
sx q[0];
rz(1.1831229) q[0];
rz(-1.7172377) q[2];
sx q[2];
rz(-0.42420039) q[2];
sx q[2];
rz(-0.37728024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.095671818) q[1];
sx q[1];
rz(-2.0883882) q[1];
sx q[1];
rz(-1.0817097) q[1];
rz(-2.7217676) q[3];
sx q[3];
rz(-1.1683162) q[3];
sx q[3];
rz(2.4703006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2034188) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(-2.7040238) q[2];
rz(-2.3214052) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1831128) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(2.130765) q[0];
rz(-0.084835947) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(-2.5929677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2261467) q[0];
sx q[0];
rz(-1.6024764) q[0];
sx q[0];
rz(1.394954) q[0];
rz(-2.7090576) q[2];
sx q[2];
rz(-1.8709595) q[2];
sx q[2];
rz(-1.6520239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1774166) q[1];
sx q[1];
rz(-0.89277041) q[1];
sx q[1];
rz(-1.7551588) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88415481) q[3];
sx q[3];
rz(-3.0100689) q[3];
sx q[3];
rz(0.92708528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34714547) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(-2.6084206) q[2];
rz(-2.063607) q[3];
sx q[3];
rz(-2.2205133) q[3];
sx q[3];
rz(-2.6326411) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0549523) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(0.78224283) q[0];
rz(-0.070146322) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(2.9152962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0840112) q[0];
sx q[0];
rz(-1.6414483) q[0];
sx q[0];
rz(-1.5018626) q[0];
rz(-pi) q[1];
rz(1.802472) q[2];
sx q[2];
rz(-2.4501928) q[2];
sx q[2];
rz(1.1065799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4676241) q[1];
sx q[1];
rz(-1.4611164) q[1];
sx q[1];
rz(-0.53916559) q[1];
x q[2];
rz(-1.4833916) q[3];
sx q[3];
rz(-0.25185302) q[3];
sx q[3];
rz(-0.15628584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(0.26257026) q[3];
sx q[3];
rz(-1.5841443) q[3];
sx q[3];
rz(-0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(-2.9578399) q[0];
rz(-2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(2.5180838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79270169) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(-1.9135273) q[0];
rz(0.97098668) q[2];
sx q[2];
rz(-2.1040593) q[2];
sx q[2];
rz(0.87109921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7694313) q[1];
sx q[1];
rz(-1.4648533) q[1];
sx q[1];
rz(1.5936567) q[1];
rz(-pi) q[2];
rz(1.1961837) q[3];
sx q[3];
rz(-1.5996029) q[3];
sx q[3];
rz(-1.0956941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(3.1039216) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(-0.17879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2035718) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(-2.6899295) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(1.2757318) q[2];
sx q[2];
rz(-2.0841334) q[2];
sx q[2];
rz(1.1708553) q[2];
rz(-1.7601109) q[3];
sx q[3];
rz(-1.2997205) q[3];
sx q[3];
rz(-1.442853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
