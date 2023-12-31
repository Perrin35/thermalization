OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51973242) q[0];
sx q[0];
rz(-0.73017263) q[0];
sx q[0];
rz(0.61868389) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47317998) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(-0.48130408) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4788471) q[1];
sx q[1];
rz(-1.4884243) q[1];
sx q[1];
rz(-1.808951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21858469) q[3];
sx q[3];
rz(-2.1544475) q[3];
sx q[3];
rz(1.0047131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(0.69357187) q[0];
rz(2.0454848) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(-2.9512761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7305671) q[0];
sx q[0];
rz(-1.9733631) q[0];
sx q[0];
rz(0.21582027) q[0];
rz(-pi) q[1];
rz(-0.53493494) q[2];
sx q[2];
rz(-1.9133647) q[2];
sx q[2];
rz(1.554622) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7479334) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(-0.02971239) q[1];
rz(-2.042949) q[3];
sx q[3];
rz(-0.90001366) q[3];
sx q[3];
rz(2.9428279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(-1.2197536) q[2];
rz(2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-2.232961) q[3];
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
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1093729) q[0];
sx q[0];
rz(-1.6633699) q[0];
sx q[0];
rz(0.60034445) q[0];
rz(-0.2072316) q[2];
sx q[2];
rz(-1.5079632) q[2];
sx q[2];
rz(-2.732423) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6229912) q[1];
sx q[1];
rz(-0.90497436) q[1];
sx q[1];
rz(1.4856505) q[1];
rz(-pi) q[2];
rz(-1.9565342) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(0.08337534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(-2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-2.9555292) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(1.8444555) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870677) q[0];
sx q[0];
rz(-1.5274807) q[0];
sx q[0];
rz(-1.4403507) q[0];
rz(-pi) q[1];
rz(-2.0879891) q[2];
sx q[2];
rz(-2.811811) q[2];
sx q[2];
rz(1.5151378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.101673) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(0.69570978) q[1];
x q[2];
rz(-1.6468871) q[3];
sx q[3];
rz(-2.9615059) q[3];
sx q[3];
rz(2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-0.91919351) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-2.2633973) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(1.426288) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8910687) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(0.10352452) q[0];
rz(-pi) q[1];
rz(0.23767383) q[2];
sx q[2];
rz(-1.2665247) q[2];
sx q[2];
rz(1.301847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6225699) q[1];
sx q[1];
rz(-1.6791324) q[1];
sx q[1];
rz(-0.83801724) q[1];
x q[2];
rz(1.8910847) q[3];
sx q[3];
rz(-2.3158145) q[3];
sx q[3];
rz(-0.68009963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(3.0997979) q[2];
rz(-3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(2.5700263) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.556658) q[0];
sx q[0];
rz(-0.99781636) q[0];
sx q[0];
rz(0.95227382) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2366441) q[2];
sx q[2];
rz(-1.4706503) q[2];
sx q[2];
rz(-3.1040994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7476269) q[1];
sx q[1];
rz(-0.52226258) q[1];
sx q[1];
rz(2.0678492) q[1];
x q[2];
rz(-2.3267641) q[3];
sx q[3];
rz(-1.0199254) q[3];
sx q[3];
rz(1.2448454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-0.56419939) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(0.66551048) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367046) q[0];
sx q[0];
rz(-0.24222736) q[0];
sx q[0];
rz(1.5750984) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7475142) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(1.9213898) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4588786) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(-2.0290124) q[1];
x q[2];
rz(2.9969278) q[3];
sx q[3];
rz(-1.8105227) q[3];
sx q[3];
rz(0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15726382) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(1.7187913) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(0.33871067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074293) q[0];
sx q[0];
rz(-1.9142262) q[0];
sx q[0];
rz(-0.15983454) q[0];
rz(-pi) q[1];
rz(-2.68967) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(0.6349596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.72570669) q[1];
sx q[1];
rz(-2.6542414) q[1];
sx q[1];
rz(-1.8094256) q[1];
rz(-pi) q[2];
rz(1.8975774) q[3];
sx q[3];
rz(-1.8339001) q[3];
sx q[3];
rz(0.14740482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.609628) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7394373) q[0];
sx q[0];
rz(-1.5097152) q[0];
sx q[0];
rz(-0.020045965) q[0];
rz(-0.058768674) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(1.8975443) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82397205) q[1];
sx q[1];
rz(-0.52597731) q[1];
sx q[1];
rz(-3.0082263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6406052) q[3];
sx q[3];
rz(-0.6456635) q[3];
sx q[3];
rz(-1.071969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6167986) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(3.0537135) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164924) q[0];
sx q[0];
rz(-0.52545588) q[0];
sx q[0];
rz(1.0111965) q[0];
rz(-pi) q[1];
rz(2.9160203) q[2];
sx q[2];
rz(-1.8403056) q[2];
sx q[2];
rz(0.05664209) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2828196) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(1.3057083) q[1];
x q[2];
rz(2.2757169) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(2.3072484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2075656) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(2.9539625) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(-1.4245695) q[2];
sx q[2];
rz(-1.8478827) q[2];
sx q[2];
rz(-2.293496) q[2];
rz(0.68708146) q[3];
sx q[3];
rz(-2.1344746) q[3];
sx q[3];
rz(1.8070756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
