OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(2.8154362) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155976) q[0];
sx q[0];
rz(-1.7863569) q[0];
sx q[0];
rz(-2.9247345) q[0];
rz(-0.60249451) q[2];
sx q[2];
rz(-1.7817111) q[2];
sx q[2];
rz(0.22533016) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1118288) q[1];
sx q[1];
rz(-1.3297237) q[1];
sx q[1];
rz(-1.8087216) q[1];
x q[2];
rz(-3.1278586) q[3];
sx q[3];
rz(-0.78679689) q[3];
sx q[3];
rz(2.9247583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8618384) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(-0.88511434) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(-2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-2.2639993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7769593) q[0];
sx q[0];
rz(-1.3898802) q[0];
sx q[0];
rz(0.69676708) q[0];
rz(-pi) q[1];
rz(-2.8009731) q[2];
sx q[2];
rz(-2.7567731) q[2];
sx q[2];
rz(-2.91586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6192012) q[1];
sx q[1];
rz(-0.11208216) q[1];
sx q[1];
rz(1.2680608) q[1];
rz(-0.33438501) q[3];
sx q[3];
rz(-1.2063724) q[3];
sx q[3];
rz(0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7426804) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14625064) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(-2.8515942) q[0];
rz(2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(-0.07382948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66729546) q[0];
sx q[0];
rz(-1.4814261) q[0];
sx q[0];
rz(1.6325634) q[0];
rz(0.87186558) q[2];
sx q[2];
rz(-2.7060894) q[2];
sx q[2];
rz(-1.0812024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0970043) q[1];
sx q[1];
rz(-1.6267596) q[1];
sx q[1];
rz(-1.7486649) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31432398) q[3];
sx q[3];
rz(-0.36800185) q[3];
sx q[3];
rz(-0.60955334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(-2.0329287) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17764238) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(-1.1886764) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68543363) q[0];
sx q[0];
rz(-1.478727) q[0];
sx q[0];
rz(2.409056) q[0];
rz(-pi) q[1];
rz(-0.75886274) q[2];
sx q[2];
rz(-2.7042411) q[2];
sx q[2];
rz(-2.8749089) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77364591) q[1];
sx q[1];
rz(-1.918643) q[1];
sx q[1];
rz(-2.0639973) q[1];
rz(-1.4401682) q[3];
sx q[3];
rz(-1.4454953) q[3];
sx q[3];
rz(-3.0076722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58549762) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(2.0193224) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(-2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68200237) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(-2.7686152) q[0];
rz(0.22398082) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7861159) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(-1.0771846) q[0];
rz(-pi) q[1];
rz(-2.3773642) q[2];
sx q[2];
rz(-0.35596213) q[2];
sx q[2];
rz(-2.228235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74354467) q[1];
sx q[1];
rz(-0.053357031) q[1];
sx q[1];
rz(1.3392901) q[1];
rz(-0.76977323) q[3];
sx q[3];
rz(-2.4379745) q[3];
sx q[3];
rz(2.9911656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-2.3357847) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(-2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39078113) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(0.85995752) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(-1.8213173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4754007) q[0];
sx q[0];
rz(-1.1328508) q[0];
sx q[0];
rz(-1.9802718) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8684902) q[2];
sx q[2];
rz(-0.58758508) q[2];
sx q[2];
rz(0.36537974) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7399866) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(2.2120038) q[1];
rz(0.28989132) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(2.1260726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0423353) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(2.5406204) q[2];
rz(-0.48505923) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(-1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(-0.55554187) q[0];
rz(-3.1069966) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.7506036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7000274) q[0];
sx q[0];
rz(-1.0374984) q[0];
sx q[0];
rz(-1.1839068) q[0];
rz(-pi) q[1];
rz(-0.58629845) q[2];
sx q[2];
rz(-0.43160298) q[2];
sx q[2];
rz(1.7238215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.890306) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(0.89819737) q[1];
rz(3.0656747) q[3];
sx q[3];
rz(-2.3644991) q[3];
sx q[3];
rz(-1.0958375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.81327072) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(2.7461046) q[2];
rz(-1.288712) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(-0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(-1.4920374) q[0];
rz(2.18816) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(-1.7038201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6839562) q[0];
sx q[0];
rz(-0.55011612) q[0];
sx q[0];
rz(1.2162627) q[0];
rz(-pi) q[1];
rz(-1.4346801) q[2];
sx q[2];
rz(-1.4021177) q[2];
sx q[2];
rz(2.8240734) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86086035) q[1];
sx q[1];
rz(-2.3128465) q[1];
sx q[1];
rz(-2.4651205) q[1];
rz(-pi) q[2];
rz(-1.1233166) q[3];
sx q[3];
rz(-1.5888702) q[3];
sx q[3];
rz(1.3656838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(2.9628741) q[2];
rz(0.86137613) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(-0.3716968) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9361967) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(-0.73668346) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(2.0827983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878384) q[0];
sx q[0];
rz(-0.61843473) q[0];
sx q[0];
rz(2.0536325) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86973127) q[2];
sx q[2];
rz(-0.5898925) q[2];
sx q[2];
rz(1.9260977) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4150548) q[1];
sx q[1];
rz(-1.448436) q[1];
sx q[1];
rz(-1.9499669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1684253) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(-2.6228867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(-1.6607364) q[2];
rz(-2.7311834) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720298) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-2.5323396) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(-1.7094918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8490484) q[0];
sx q[0];
rz(-1.8869072) q[0];
sx q[0];
rz(2.6228117) q[0];
rz(-2.6851995) q[2];
sx q[2];
rz(-2.6420339) q[2];
sx q[2];
rz(-1.511614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4869625) q[1];
sx q[1];
rz(-1.6143867) q[1];
sx q[1];
rz(-0.55200465) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7581975) q[3];
sx q[3];
rz(-2.7029576) q[3];
sx q[3];
rz(2.8965829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6802406) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(1.6067778) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(-2.6721568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(2.7728511) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(0.22696162) q[2];
sx q[2];
rz(-1.771275) q[2];
sx q[2];
rz(2.7834409) q[2];
rz(0.41714824) q[3];
sx q[3];
rz(-2.683831) q[3];
sx q[3];
rz(2.8575069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
