OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.51337564) q[0];
sx q[0];
rz(-1.6802508) q[0];
sx q[0];
rz(2.2267377) q[0];
rz(-0.86923161) q[1];
sx q[1];
rz(-2.4183122) q[1];
sx q[1];
rz(1.7887315) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2268838) q[0];
sx q[0];
rz(-1.0585203) q[0];
sx q[0];
rz(-0.73988503) q[0];
rz(-pi) q[1];
rz(1.1961785) q[2];
sx q[2];
rz(-2.2785419) q[2];
sx q[2];
rz(-1.5951235) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1654403) q[1];
sx q[1];
rz(-0.70164499) q[1];
sx q[1];
rz(-1.8681504) q[1];
x q[2];
rz(0.21463359) q[3];
sx q[3];
rz(-2.1691536) q[3];
sx q[3];
rz(-0.90902381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73504084) q[2];
sx q[2];
rz(-0.81321365) q[2];
sx q[2];
rz(2.1347894) q[2];
rz(-1.462228) q[3];
sx q[3];
rz(-1.6961325) q[3];
sx q[3];
rz(-1.6026976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891069) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(-0.12595969) q[0];
rz(2.0055298) q[1];
sx q[1];
rz(-1.845153) q[1];
sx q[1];
rz(2.0819285) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39108178) q[0];
sx q[0];
rz(-1.568887) q[0];
sx q[0];
rz(-1.5715412) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66920377) q[2];
sx q[2];
rz(-2.0811709) q[2];
sx q[2];
rz(-0.83590311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4821533) q[1];
sx q[1];
rz(-1.1350613) q[1];
sx q[1];
rz(-2.3547291) q[1];
rz(-0.13215315) q[3];
sx q[3];
rz(-0.36214585) q[3];
sx q[3];
rz(-1.3782448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5560567) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(2.5118828) q[2];
rz(-1.2398047) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(-2.0201717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.020551) q[0];
sx q[0];
rz(-2.910055) q[0];
sx q[0];
rz(1.2303906) q[0];
rz(-0.70368272) q[1];
sx q[1];
rz(-2.2129702) q[1];
sx q[1];
rz(-1.5603265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9890922) q[0];
sx q[0];
rz(-1.9643503) q[0];
sx q[0];
rz(3.1253184) q[0];
rz(-pi) q[1];
rz(1.376344) q[2];
sx q[2];
rz(-2.299435) q[2];
sx q[2];
rz(0.096708111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3250503) q[1];
sx q[1];
rz(-2.5245164) q[1];
sx q[1];
rz(0.80775921) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3558667) q[3];
sx q[3];
rz(-1.8453987) q[3];
sx q[3];
rz(-0.532632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33502093) q[2];
sx q[2];
rz(-1.0806012) q[2];
sx q[2];
rz(-3.0873155) q[2];
rz(-1.0862167) q[3];
sx q[3];
rz(-2.351206) q[3];
sx q[3];
rz(-0.49200341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7739968) q[0];
sx q[0];
rz(-3.1258899) q[0];
sx q[0];
rz(0.98454654) q[0];
rz(-1.4060075) q[1];
sx q[1];
rz(-1.6008629) q[1];
sx q[1];
rz(-2.9071009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4658452) q[0];
sx q[0];
rz(-2.391531) q[0];
sx q[0];
rz(1.1613599) q[0];
x q[1];
rz(-0.245204) q[2];
sx q[2];
rz(-1.8498932) q[2];
sx q[2];
rz(-2.9989105) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4602743) q[1];
sx q[1];
rz(-1.360962) q[1];
sx q[1];
rz(-0.55822452) q[1];
x q[2];
rz(0.99458875) q[3];
sx q[3];
rz(-1.595701) q[3];
sx q[3];
rz(-1.8751609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24388127) q[2];
sx q[2];
rz(-0.66456866) q[2];
sx q[2];
rz(-2.393874) q[2];
rz(-2.054935) q[3];
sx q[3];
rz(-0.99435157) q[3];
sx q[3];
rz(-0.35191107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.20027593) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(1.548832) q[0];
rz(-3.0039655) q[1];
sx q[1];
rz(-2.177114) q[1];
sx q[1];
rz(2.4639938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9842868) q[0];
sx q[0];
rz(-1.5781227) q[0];
sx q[0];
rz(-0.013989438) q[0];
x q[1];
rz(-0.65998544) q[2];
sx q[2];
rz(-2.8980245) q[2];
sx q[2];
rz(-3.0483831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69285821) q[1];
sx q[1];
rz(-3.0046607) q[1];
sx q[1];
rz(-1.6137692) q[1];
x q[2];
rz(-2.4047071) q[3];
sx q[3];
rz(-1.3892655) q[3];
sx q[3];
rz(-1.5893861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4981093) q[2];
sx q[2];
rz(-2.6124239) q[2];
sx q[2];
rz(0.37437487) q[2];
rz(0.74786413) q[3];
sx q[3];
rz(-1.4948083) q[3];
sx q[3];
rz(1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0102274) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(0.22450547) q[0];
rz(-0.35047105) q[1];
sx q[1];
rz(-1.9572565) q[1];
sx q[1];
rz(1.1856461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.810038) q[0];
sx q[0];
rz(-1.4228978) q[0];
sx q[0];
rz(0.055813544) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6780532) q[2];
sx q[2];
rz(-1.4847418) q[2];
sx q[2];
rz(1.254509) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.396186) q[1];
sx q[1];
rz(-0.90382677) q[1];
sx q[1];
rz(1.0140258) q[1];
rz(-pi) q[2];
rz(0.91573651) q[3];
sx q[3];
rz(-2.1300507) q[3];
sx q[3];
rz(2.9489234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.4993569) q[2];
sx q[2];
rz(-2.2580052) q[2];
sx q[2];
rz(-0.16481608) q[2];
rz(2.3757101) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(-0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931452) q[0];
sx q[0];
rz(-3.1146545) q[0];
sx q[0];
rz(2.7222166) q[0];
rz(-1.558149) q[1];
sx q[1];
rz(-2.7132468) q[1];
sx q[1];
rz(-1.3126119) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2079577) q[0];
sx q[0];
rz(-2.5267753) q[0];
sx q[0];
rz(1.2864248) q[0];
rz(0.75757005) q[2];
sx q[2];
rz(-0.84572863) q[2];
sx q[2];
rz(0.49882327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.928682) q[1];
sx q[1];
rz(-2.8281459) q[1];
sx q[1];
rz(-0.9153961) q[1];
x q[2];
rz(1.4410331) q[3];
sx q[3];
rz(-1.8365873) q[3];
sx q[3];
rz(2.5523976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4075564) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(2.7835795) q[2];
rz(-2.9426306) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(-3.0437886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1426549) q[0];
sx q[0];
rz(-1.0336579) q[0];
sx q[0];
rz(2.348483) q[0];
rz(2.2456887) q[1];
sx q[1];
rz(-1.221012) q[1];
sx q[1];
rz(0.47826794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2390564) q[0];
sx q[0];
rz(-2.0180114) q[0];
sx q[0];
rz(-1.4330401) q[0];
rz(1.9321279) q[2];
sx q[2];
rz(-1.981297) q[2];
sx q[2];
rz(2.362006) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9517802) q[1];
sx q[1];
rz(-2.5972022) q[1];
sx q[1];
rz(1.0448827) q[1];
x q[2];
rz(1.0555154) q[3];
sx q[3];
rz(-1.6664423) q[3];
sx q[3];
rz(-3.1351922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6200977) q[2];
sx q[2];
rz(-1.564881) q[2];
sx q[2];
rz(-0.0077956789) q[2];
rz(-1.1826285) q[3];
sx q[3];
rz(-1.7595485) q[3];
sx q[3];
rz(1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1496534) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(-2.7329408) q[0];
rz(1.0952134) q[1];
sx q[1];
rz(-1.6588255) q[1];
sx q[1];
rz(-2.2832787) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3214071) q[0];
sx q[0];
rz(-1.5044065) q[0];
sx q[0];
rz(1.4231349) q[0];
rz(-pi) q[1];
rz(-1.8377531) q[2];
sx q[2];
rz(-2.6745195) q[2];
sx q[2];
rz(2.2427223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9954736) q[1];
sx q[1];
rz(-1.3750556) q[1];
sx q[1];
rz(-0.11222307) q[1];
rz(0.59767234) q[3];
sx q[3];
rz(-1.3857171) q[3];
sx q[3];
rz(1.5782034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3506763) q[2];
sx q[2];
rz(-2.3886267) q[2];
sx q[2];
rz(-0.6616627) q[2];
rz(2.3412797) q[3];
sx q[3];
rz(-1.9789663) q[3];
sx q[3];
rz(0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8372339) q[0];
sx q[0];
rz(-2.9999314) q[0];
sx q[0];
rz(2.4270571) q[0];
rz(2.6658658) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(0.48019662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59360945) q[0];
sx q[0];
rz(-2.3385323) q[0];
sx q[0];
rz(-0.70269967) q[0];
rz(-pi) q[1];
rz(-2.8101361) q[2];
sx q[2];
rz(-2.1961947) q[2];
sx q[2];
rz(-1.3783704) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.678315) q[1];
sx q[1];
rz(-2.2147605) q[1];
sx q[1];
rz(-1.9309082) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28292029) q[3];
sx q[3];
rz(-1.991254) q[3];
sx q[3];
rz(-0.69129163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51018888) q[2];
sx q[2];
rz(-1.1978585) q[2];
sx q[2];
rz(-0.27604827) q[2];
rz(-0.43802842) q[3];
sx q[3];
rz(-1.5166401) q[3];
sx q[3];
rz(-0.62694222) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.940687) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(-3.1055462) q[1];
sx q[1];
rz(-1.5279952) q[1];
sx q[1];
rz(-1.5855018) q[1];
rz(-0.27388688) q[2];
sx q[2];
rz(-1.0928434) q[2];
sx q[2];
rz(2.6437987) q[2];
rz(-2.8874019) q[3];
sx q[3];
rz(-1.7028451) q[3];
sx q[3];
rz(-2.9447671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
