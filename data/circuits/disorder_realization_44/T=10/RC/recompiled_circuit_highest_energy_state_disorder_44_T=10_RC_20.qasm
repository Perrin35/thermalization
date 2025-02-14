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
rz(1.1062082) q[0];
sx q[0];
rz(-2.455403) q[0];
sx q[0];
rz(-1.9880779) q[0];
rz(1.310362) q[1];
sx q[1];
rz(-0.49340931) q[1];
sx q[1];
rz(-0.44841132) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759244) q[0];
sx q[0];
rz(-2.2184555) q[0];
sx q[0];
rz(-1.5041385) q[0];
rz(-0.37031074) q[2];
sx q[2];
rz(-2.8680621) q[2];
sx q[2];
rz(-0.47765484) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7668575) q[1];
sx q[1];
rz(-2.285706) q[1];
sx q[1];
rz(2.292407) q[1];
rz(-pi) q[2];
rz(1.5228062) q[3];
sx q[3];
rz(-1.0700883) q[3];
sx q[3];
rz(-1.9013311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51556921) q[2];
sx q[2];
rz(-1.7531351) q[2];
sx q[2];
rz(2.5173729) q[2];
rz(0.37912399) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(2.8930801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9666331) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(-1.0354743) q[0];
rz(1.8677208) q[1];
sx q[1];
rz(-2.404411) q[1];
sx q[1];
rz(2.6587291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55601701) q[0];
sx q[0];
rz(-0.10679467) q[0];
sx q[0];
rz(0.44321816) q[0];
rz(-pi) q[1];
rz(1.4522885) q[2];
sx q[2];
rz(-1.7359455) q[2];
sx q[2];
rz(1.120795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.00512) q[1];
sx q[1];
rz(-1.1389975) q[1];
sx q[1];
rz(0.010331945) q[1];
rz(-pi) q[2];
rz(1.6522406) q[3];
sx q[3];
rz(-1.8907428) q[3];
sx q[3];
rz(0.68405747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1172993) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(-1.8710322) q[2];
rz(2.471586) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4033177) q[0];
sx q[0];
rz(-2.6955695) q[0];
sx q[0];
rz(2.4025412) q[0];
rz(-0.96145472) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(-2.3846073) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805646) q[0];
sx q[0];
rz(-1.3098728) q[0];
sx q[0];
rz(-2.3105614) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1951094) q[2];
sx q[2];
rz(-0.91478148) q[2];
sx q[2];
rz(0.13740787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18338246) q[1];
sx q[1];
rz(-0.64309769) q[1];
sx q[1];
rz(0.22497589) q[1];
rz(1.8784896) q[3];
sx q[3];
rz(-1.5041927) q[3];
sx q[3];
rz(2.8860725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(-0.70510954) q[2];
rz(-2.8201568) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(-0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0099156378) q[0];
sx q[0];
rz(-0.83808815) q[0];
sx q[0];
rz(1.2257082) q[0];
rz(-1.2340087) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(-0.066224901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4915337) q[0];
sx q[0];
rz(-2.8797132) q[0];
sx q[0];
rz(-1.9240379) q[0];
rz(-pi) q[1];
rz(-2.8429864) q[2];
sx q[2];
rz(-0.82393194) q[2];
sx q[2];
rz(-2.2366984) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.855763) q[1];
sx q[1];
rz(-1.0914088) q[1];
sx q[1];
rz(-2.781032) q[1];
x q[2];
rz(-1.6965167) q[3];
sx q[3];
rz(-1.9341521) q[3];
sx q[3];
rz(1.1896127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0464728) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(-0.44060102) q[2];
rz(-2.1353841) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42136583) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(-1.6857326) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855749) q[0];
sx q[0];
rz(-0.7472207) q[0];
sx q[0];
rz(-1.78563) q[0];
rz(-pi) q[1];
rz(-0.22730374) q[2];
sx q[2];
rz(-1.6844201) q[2];
sx q[2];
rz(-0.34662468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.887874) q[1];
sx q[1];
rz(-1.591133) q[1];
sx q[1];
rz(0.65171297) q[1];
rz(-pi) q[2];
rz(-1.7248575) q[3];
sx q[3];
rz(-0.96821456) q[3];
sx q[3];
rz(-1.3101206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8684034) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(-0.91147649) q[2];
rz(-0.8738001) q[3];
sx q[3];
rz(-2.3995212) q[3];
sx q[3];
rz(0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28354302) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-3.0112596) q[0];
rz(0.48769543) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(2.3977051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7672507) q[0];
sx q[0];
rz(-3.0453186) q[0];
sx q[0];
rz(-1.0262579) q[0];
x q[1];
rz(-2.1962847) q[2];
sx q[2];
rz(-1.4577366) q[2];
sx q[2];
rz(0.94719749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0186613) q[1];
sx q[1];
rz(-2.5187188) q[1];
sx q[1];
rz(-2.0707002) q[1];
rz(-pi) q[2];
rz(-1.7508932) q[3];
sx q[3];
rz(-1.0922179) q[3];
sx q[3];
rz(-0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1193715) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(-2.7458701) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(-0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(-1.047026) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.2774757) q[1];
sx q[1];
rz(-2.7488757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34907863) q[0];
sx q[0];
rz(-0.70988467) q[0];
sx q[0];
rz(1.841808) q[0];
rz(-pi) q[1];
rz(0.057633295) q[2];
sx q[2];
rz(-2.1890292) q[2];
sx q[2];
rz(-2.1098441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9153292) q[1];
sx q[1];
rz(-1.2955503) q[1];
sx q[1];
rz(1.5987426) q[1];
rz(-pi) q[2];
rz(0.32870558) q[3];
sx q[3];
rz(-1.8227326) q[3];
sx q[3];
rz(3.1303034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0543694) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(-1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8488309) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(-2.8344179) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(0.90075341) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7194289) q[0];
sx q[0];
rz(-1.3086645) q[0];
sx q[0];
rz(1.5533226) q[0];
x q[1];
rz(1.0020761) q[2];
sx q[2];
rz(-2.4042712) q[2];
sx q[2];
rz(-1.3739283) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5985721) q[1];
sx q[1];
rz(-1.2391587) q[1];
sx q[1];
rz(-0.33473067) q[1];
rz(2.9489904) q[3];
sx q[3];
rz(-0.64036548) q[3];
sx q[3];
rz(2.7964724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8073392) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(-2.5267498) q[2];
rz(2.0452512) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(-0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480302) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(0.50931859) q[0];
rz(-0.18109426) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(0.96955713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2354202) q[0];
sx q[0];
rz(-0.83011043) q[0];
sx q[0];
rz(-1.3614015) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0779586) q[2];
sx q[2];
rz(-1.615287) q[2];
sx q[2];
rz(0.26598334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4050715) q[1];
sx q[1];
rz(-0.35686359) q[1];
sx q[1];
rz(-0.31330152) q[1];
rz(-pi) q[2];
rz(1.9593616) q[3];
sx q[3];
rz(-2.3747184) q[3];
sx q[3];
rz(0.23393133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8934882) q[2];
sx q[2];
rz(-2.9073145) q[2];
sx q[2];
rz(2.060037) q[2];
rz(2.9099416) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.3807826) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(-0.64714062) q[0];
rz(2.6864247) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(0.761935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021558048) q[0];
sx q[0];
rz(-1.2823449) q[0];
sx q[0];
rz(-0.7201654) q[0];
rz(-0.89770395) q[2];
sx q[2];
rz(-2.2676149) q[2];
sx q[2];
rz(0.76329939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0207723) q[1];
sx q[1];
rz(-1.6050395) q[1];
sx q[1];
rz(2.4832151) q[1];
rz(-1.1033789) q[3];
sx q[3];
rz(-1.5776792) q[3];
sx q[3];
rz(-0.21871834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.053293856) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(2.2376412) q[2];
rz(1.5229185) q[3];
sx q[3];
rz(-1.0205597) q[3];
sx q[3];
rz(1.9797549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91697964) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-3.0147973) q[1];
sx q[1];
rz(-2.3190111) q[1];
sx q[1];
rz(0.93228985) q[1];
rz(-0.2520855) q[2];
sx q[2];
rz(-1.132195) q[2];
sx q[2];
rz(1.8762527) q[2];
rz(-0.42999646) q[3];
sx q[3];
rz(-1.5734133) q[3];
sx q[3];
rz(-1.6071241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
