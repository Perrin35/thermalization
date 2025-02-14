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
rz(0.053057916) q[0];
sx q[0];
rz(-2.7795656) q[0];
sx q[0];
rz(1.9757353) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(-1.7133763) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7390763) q[0];
sx q[0];
rz(-1.5105643) q[0];
sx q[0];
rz(0.51243685) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25144724) q[2];
sx q[2];
rz(-1.5000952) q[2];
sx q[2];
rz(1.9858907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1686656) q[1];
sx q[1];
rz(-1.7024689) q[1];
sx q[1];
rz(1.3802856) q[1];
rz(2.0479508) q[3];
sx q[3];
rz(-0.53864281) q[3];
sx q[3];
rz(-0.25687309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6894138) q[2];
sx q[2];
rz(-0.016760085) q[2];
sx q[2];
rz(0.20081271) q[2];
rz(-2.9928442) q[3];
sx q[3];
rz(-3.1368308) q[3];
sx q[3];
rz(2.8301921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1397322) q[0];
sx q[0];
rz(-0.59356028) q[0];
sx q[0];
rz(-2.1019905) q[0];
rz(3.1270341) q[1];
sx q[1];
rz(-1.9083551) q[1];
sx q[1];
rz(-1.5878409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36984617) q[0];
sx q[0];
rz(-2.2622007) q[0];
sx q[0];
rz(-2.2798674) q[0];
rz(-pi) q[1];
rz(1.3755685) q[2];
sx q[2];
rz(-3.0659817) q[2];
sx q[2];
rz(1.6603254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0224198) q[1];
sx q[1];
rz(-1.8298365) q[1];
sx q[1];
rz(-3.1062192) q[1];
rz(-pi) q[2];
rz(1.3782937) q[3];
sx q[3];
rz(-1.2436335) q[3];
sx q[3];
rz(-0.006803676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62850922) q[2];
sx q[2];
rz(-1.6078948) q[2];
sx q[2];
rz(-1.3879363) q[2];
rz(1.3758818) q[3];
sx q[3];
rz(-2.0949771) q[3];
sx q[3];
rz(-2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917711) q[0];
sx q[0];
rz(-0.23673683) q[0];
sx q[0];
rz(-2.5340875) q[0];
rz(-1.5433743) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(2.1764596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34244363) q[0];
sx q[0];
rz(-1.5502872) q[0];
sx q[0];
rz(0.0040702013) q[0];
rz(2.1943614) q[2];
sx q[2];
rz(-0.5493702) q[2];
sx q[2];
rz(1.6877162) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9671792) q[1];
sx q[1];
rz(-1.426218) q[1];
sx q[1];
rz(-0.036199613) q[1];
x q[2];
rz(-0.083250982) q[3];
sx q[3];
rz(-1.4752582) q[3];
sx q[3];
rz(-0.18886177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8049916) q[2];
sx q[2];
rz(-2.470863) q[2];
sx q[2];
rz(-2.2750308) q[2];
rz(2.0349515) q[3];
sx q[3];
rz(-1.589078) q[3];
sx q[3];
rz(1.6710056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5690145) q[0];
sx q[0];
rz(-2.4622279) q[0];
sx q[0];
rz(1.6160075) q[0];
rz(3.1309879) q[1];
sx q[1];
rz(-3.1378101) q[1];
sx q[1];
rz(-2.3912281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5724154) q[0];
sx q[0];
rz(-0.67052013) q[0];
sx q[0];
rz(2.9886093) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7979787) q[2];
sx q[2];
rz(-2.1641693) q[2];
sx q[2];
rz(2.5636755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039783104) q[1];
sx q[1];
rz(-1.0652115) q[1];
sx q[1];
rz(1.4391446) q[1];
rz(-pi) q[2];
rz(1.4441277) q[3];
sx q[3];
rz(-1.2614246) q[3];
sx q[3];
rz(-0.29361967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41287199) q[2];
sx q[2];
rz(-1.0888638) q[2];
sx q[2];
rz(-1.2415761) q[2];
rz(3.1359172) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(2.2976105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64549696) q[0];
sx q[0];
rz(-0.050300751) q[0];
sx q[0];
rz(0.92329931) q[0];
rz(-2.3410489) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(-2.961535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411856) q[0];
sx q[0];
rz(-1.5292455) q[0];
sx q[0];
rz(-1.8096258) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4902968) q[2];
sx q[2];
rz(-1.4729807) q[2];
sx q[2];
rz(2.7076633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48985607) q[1];
sx q[1];
rz(-1.3313313) q[1];
sx q[1];
rz(-0.36291562) q[1];
rz(1.2050797) q[3];
sx q[3];
rz(-1.4673181) q[3];
sx q[3];
rz(0.24170719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30183145) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(-1.4361471) q[2];
rz(-1.2561717) q[3];
sx q[3];
rz(-1.4830282) q[3];
sx q[3];
rz(-0.095631599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489478) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(-2.6982464) q[0];
rz(-1.1706785) q[1];
sx q[1];
rz(-0.0010633855) q[1];
sx q[1];
rz(-2.7098157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1373119) q[0];
sx q[0];
rz(-2.3385424) q[0];
sx q[0];
rz(2.7025168) q[0];
x q[1];
rz(3.058039) q[2];
sx q[2];
rz(-1.3992893) q[2];
sx q[2];
rz(-0.87860859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37913191) q[1];
sx q[1];
rz(-2.1239696) q[1];
sx q[1];
rz(1.1198612) q[1];
x q[2];
rz(-0.59817846) q[3];
sx q[3];
rz(-2.7477187) q[3];
sx q[3];
rz(2.1109714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7410437) q[2];
sx q[2];
rz(-2.2811175) q[2];
sx q[2];
rz(-1.384037) q[2];
rz(-0.69796973) q[3];
sx q[3];
rz(-2.3845086) q[3];
sx q[3];
rz(-0.32148263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8292002) q[0];
sx q[0];
rz(-1.3425403) q[0];
sx q[0];
rz(-0.37305748) q[0];
rz(-0.27698764) q[1];
sx q[1];
rz(-0.00033683446) q[1];
sx q[1];
rz(2.3854947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5146268) q[0];
sx q[0];
rz(-1.9494162) q[0];
sx q[0];
rz(-2.4832804) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.937708) q[2];
sx q[2];
rz(-1.1808625) q[2];
sx q[2];
rz(-0.45118794) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4226729) q[1];
sx q[1];
rz(-0.70168272) q[1];
sx q[1];
rz(0.91928457) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30310615) q[3];
sx q[3];
rz(-1.167093) q[3];
sx q[3];
rz(1.2991326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9296391) q[2];
sx q[2];
rz(-2.5888207) q[2];
sx q[2];
rz(-2.2201404) q[2];
rz(-3.0742505) q[3];
sx q[3];
rz(-1.208297) q[3];
sx q[3];
rz(-1.4837846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15143722) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(-3.0113599) q[0];
rz(-2.3479334) q[1];
sx q[1];
rz(-0.001611324) q[1];
sx q[1];
rz(-0.28009716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6009499) q[0];
sx q[0];
rz(-1.6495935) q[0];
sx q[0];
rz(-0.053683563) q[0];
x q[1];
rz(2.9236004) q[2];
sx q[2];
rz(-1.6577066) q[2];
sx q[2];
rz(-0.72959585) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62781266) q[1];
sx q[1];
rz(-0.65142056) q[1];
sx q[1];
rz(0.86077229) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6032366) q[3];
sx q[3];
rz(-2.3545485) q[3];
sx q[3];
rz(0.53883906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.085999504) q[2];
sx q[2];
rz(-1.6722101) q[2];
sx q[2];
rz(-0.80297339) q[2];
rz(-1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(0.83139658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0292173) q[0];
sx q[0];
rz(-0.0028828415) q[0];
sx q[0];
rz(3.032384) q[0];
rz(0.373492) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(0.55140299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38574252) q[0];
sx q[0];
rz(-1.047121) q[0];
sx q[0];
rz(-0.68501212) q[0];
rz(-pi) q[1];
rz(-2.5481651) q[2];
sx q[2];
rz(-2.9205236) q[2];
sx q[2];
rz(1.204837) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0950737) q[1];
sx q[1];
rz(-1.0380901) q[1];
sx q[1];
rz(0.84360509) q[1];
rz(1.3250454) q[3];
sx q[3];
rz(-1.9634523) q[3];
sx q[3];
rz(0.33473864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6138844) q[2];
sx q[2];
rz(-1.3091427) q[2];
sx q[2];
rz(-1.323553) q[2];
rz(1.8771133) q[3];
sx q[3];
rz(-1.2885965) q[3];
sx q[3];
rz(0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5453813) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(2.4052461) q[0];
rz(0.21190724) q[1];
sx q[1];
rz(-2.1129463) q[1];
sx q[1];
rz(-1.597499) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6074267) q[0];
sx q[0];
rz(-1.5527927) q[0];
sx q[0];
rz(-0.27915633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5681463) q[2];
sx q[2];
rz(-1.6011597) q[2];
sx q[2];
rz(1.7372436) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4344132) q[1];
sx q[1];
rz(-2.3418171) q[1];
sx q[1];
rz(2.4753285) q[1];
rz(-0.44477099) q[3];
sx q[3];
rz(-1.8666779) q[3];
sx q[3];
rz(2.5887161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3046917) q[2];
sx q[2];
rz(-2.3062314) q[2];
sx q[2];
rz(1.2738073) q[2];
rz(-1.4443719) q[3];
sx q[3];
rz(-0.074051753) q[3];
sx q[3];
rz(1.4950289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.973751) q[0];
sx q[0];
rz(-1.5580307) q[0];
sx q[0];
rz(1.8488484) q[0];
rz(-1.6043067) q[1];
sx q[1];
rz(-0.91265408) q[1];
sx q[1];
rz(0.18462054) q[1];
rz(2.9188048) q[2];
sx q[2];
rz(-0.041112715) q[2];
sx q[2];
rz(1.0512752) q[2];
rz(1.3956232) q[3];
sx q[3];
rz(-2.5268231) q[3];
sx q[3];
rz(-2.3781621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
