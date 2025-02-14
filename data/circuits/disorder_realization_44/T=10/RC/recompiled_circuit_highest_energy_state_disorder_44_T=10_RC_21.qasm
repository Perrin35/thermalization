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
rz(-2.0353844) q[0];
sx q[0];
rz(2.455403) q[0];
sx q[0];
rz(10.578293) q[0];
rz(1.310362) q[1];
sx q[1];
rz(-0.49340931) q[1];
sx q[1];
rz(-0.44841132) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656682) q[0];
sx q[0];
rz(-2.2184555) q[0];
sx q[0];
rz(1.6374542) q[0];
rz(1.6719867) q[2];
sx q[2];
rz(-1.3162344) q[2];
sx q[2];
rz(-2.2806666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.37473512) q[1];
sx q[1];
rz(-0.85588662) q[1];
sx q[1];
rz(-0.84918569) q[1];
rz(-1.6187864) q[3];
sx q[3];
rz(-1.0700883) q[3];
sx q[3];
rz(1.2402616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-0.62421978) q[2];
rz(-2.7624687) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(2.8930801) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17495951) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(1.0354743) q[0];
rz(1.2738719) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(2.6587291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55601701) q[0];
sx q[0];
rz(-0.10679467) q[0];
sx q[0];
rz(2.6983745) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16629433) q[2];
sx q[2];
rz(-1.6876843) q[2];
sx q[2];
rz(0.4695732) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.711593) q[1];
sx q[1];
rz(-1.5801799) q[1];
sx q[1];
rz(2.0026155) q[1];
x q[2];
rz(1.489352) q[3];
sx q[3];
rz(-1.8907428) q[3];
sx q[3];
rz(2.4575352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.024293385) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(1.8710322) q[2];
rz(-2.471586) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73827493) q[0];
sx q[0];
rz(-2.6955695) q[0];
sx q[0];
rz(2.4025412) q[0];
rz(2.1801379) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(2.3846073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805646) q[0];
sx q[0];
rz(-1.8317199) q[0];
sx q[0];
rz(0.83103128) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6967825) q[2];
sx q[2];
rz(-2.399657) q[2];
sx q[2];
rz(-2.7049989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5685421) q[1];
sx q[1];
rz(-1.4366163) q[1];
sx q[1];
rz(-2.5106984) q[1];
x q[2];
rz(1.7875761) q[3];
sx q[3];
rz(-0.31459537) q[3];
sx q[3];
rz(1.5217239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(-0.70510954) q[2];
rz(-2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0099156378) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(1.9158844) q[0];
rz(1.2340087) q[1];
sx q[1];
rz(-1.0154513) q[1];
sx q[1];
rz(-0.066224901) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2628675) q[0];
sx q[0];
rz(-1.4811133) q[0];
sx q[0];
rz(-1.8171726) q[0];
rz(-pi) q[1];
rz(2.340256) q[2];
sx q[2];
rz(-1.3532012) q[2];
sx q[2];
rz(-2.269553) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5410616) q[1];
sx q[1];
rz(-0.59128535) q[1];
sx q[1];
rz(-0.97452428) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8230571) q[3];
sx q[3];
rz(-0.38357601) q[3];
sx q[3];
rz(-0.84795241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0464728) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(-0.44060102) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(-0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42136583) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(-1.6857326) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855749) q[0];
sx q[0];
rz(-2.394372) q[0];
sx q[0];
rz(1.78563) q[0];
x q[1];
rz(2.6728476) q[2];
sx q[2];
rz(-0.25368099) q[2];
sx q[2];
rz(-1.6799334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29042915) q[1];
sx q[1];
rz(-2.4896087) q[1];
sx q[1];
rz(-3.1080721) q[1];
rz(1.4167352) q[3];
sx q[3];
rz(-2.1733781) q[3];
sx q[3];
rz(-1.8314721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8684034) q[2];
sx q[2];
rz(-1.5634147) q[2];
sx q[2];
rz(0.91147649) q[2];
rz(-2.2677926) q[3];
sx q[3];
rz(-2.3995212) q[3];
sx q[3];
rz(3.1349283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.8580496) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-3.0112596) q[0];
rz(0.48769543) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(-0.74388751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9209393) q[0];
sx q[0];
rz(-1.6531117) q[0];
sx q[0];
rz(-0.049985188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7623474) q[2];
sx q[2];
rz(-2.5073176) q[2];
sx q[2];
rz(0.77855643) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6106657) q[1];
sx q[1];
rz(-2.1082889) q[1];
sx q[1];
rz(2.8100138) q[1];
rz(1.3906995) q[3];
sx q[3];
rz(-2.0493747) q[3];
sx q[3];
rz(0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0222212) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(2.7458701) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(-3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
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
rz(0.39271694) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.792514) q[0];
sx q[0];
rz(-2.431708) q[0];
sx q[0];
rz(1.2997846) q[0];
x q[1];
rz(-1.4899859) q[2];
sx q[2];
rz(-0.62056345) q[2];
sx q[2];
rz(-2.0106237) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32875529) q[1];
sx q[1];
rz(-2.8649674) q[1];
sx q[1];
rz(-0.098621086) q[1];
rz(-0.32870558) q[3];
sx q[3];
rz(-1.3188601) q[3];
sx q[3];
rz(3.1303034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.087223209) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(1.3537815) q[3];
sx q[3];
rz(-1.9446707) q[3];
sx q[3];
rz(-0.27031171) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29276174) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(2.8344179) q[0];
rz(-1.81709) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(0.90075341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7194289) q[0];
sx q[0];
rz(-1.8329282) q[0];
sx q[0];
rz(-1.5882701) q[0];
rz(-0.91761748) q[2];
sx q[2];
rz(-1.9412882) q[2];
sx q[2];
rz(-2.8965184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4169654) q[1];
sx q[1];
rz(-2.6748383) q[1];
sx q[1];
rz(0.80893597) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7124655) q[3];
sx q[3];
rz(-0.9441388) q[3];
sx q[3];
rz(0.10654813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8073392) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(0.61484289) q[2];
rz(-2.0452512) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(-2.443327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7480302) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(0.50931859) q[0];
rz(2.9604984) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(0.96955713) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345393) q[0];
sx q[0];
rz(-1.7248123) q[0];
sx q[0];
rz(0.75169433) q[0];
rz(-pi) q[1];
rz(-1.0636341) q[2];
sx q[2];
rz(-1.615287) q[2];
sx q[2];
rz(0.26598334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0694094) q[1];
sx q[1];
rz(-1.2320215) q[1];
sx q[1];
rz(1.4563926) q[1];
rz(-pi) q[2];
rz(0.84253715) q[3];
sx q[3];
rz(-1.304783) q[3];
sx q[3];
rz(-1.0501978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8934882) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(2.9099416) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(-0.20802465) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3807826) q[0];
sx q[0];
rz(-0.2247227) q[0];
sx q[0];
rz(-2.494452) q[0];
rz(2.6864247) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(0.761935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2357016) q[0];
sx q[0];
rz(-2.3755223) q[0];
sx q[0];
rz(2.7188042) q[0];
rz(0.89770395) q[2];
sx q[2];
rz(-2.2676149) q[2];
sx q[2];
rz(-0.76329939) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1208204) q[1];
sx q[1];
rz(-1.5365531) q[1];
sx q[1];
rz(-0.65837752) q[1];
rz(-pi) q[2];
rz(-1.5555218) q[3];
sx q[3];
rz(-0.46746436) q[3];
sx q[3];
rz(-1.338442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0882988) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(2.2376412) q[2];
rz(1.6186742) q[3];
sx q[3];
rz(-1.0205597) q[3];
sx q[3];
rz(1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91697964) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-0.12679535) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(2.0218718) q[2];
sx q[2];
rz(-1.3430165) q[2];
sx q[2];
rz(-2.7271885) q[2];
rz(1.5679172) q[3];
sx q[3];
rz(-1.1408014) q[3];
sx q[3];
rz(-0.037527966) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
