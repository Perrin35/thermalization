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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(-0.21633202) q[0];
rz(-2.1842015) q[1];
sx q[1];
rz(-0.67617813) q[1];
sx q[1];
rz(1.6300936) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51115655) q[0];
sx q[0];
rz(-1.5611708) q[0];
sx q[0];
rz(-2.5853392) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6681603) q[2];
sx q[2];
rz(-1.8634999) q[2];
sx q[2];
rz(2.6930489) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2701157) q[1];
sx q[1];
rz(-1.9028712) q[1];
sx q[1];
rz(-2.6144652) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2735748) q[3];
sx q[3];
rz(-1.8773762) q[3];
sx q[3];
rz(1.0667232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1067074) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(1.0484288) q[2];
rz(0.18167051) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(-2.488193) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6492017) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(-2.6948068) q[0];
rz(-2.2611639) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(2.3588038) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0271757) q[0];
sx q[0];
rz(-1.3265298) q[0];
sx q[0];
rz(0.063112325) q[0];
x q[1];
rz(2.5704284) q[2];
sx q[2];
rz(-0.70020247) q[2];
sx q[2];
rz(-2.458606) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70322733) q[1];
sx q[1];
rz(-1.8706873) q[1];
sx q[1];
rz(-0.48734003) q[1];
rz(-pi) q[2];
rz(-0.47329013) q[3];
sx q[3];
rz(-0.9447228) q[3];
sx q[3];
rz(1.9260709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.030152628) q[2];
sx q[2];
rz(-1.6609265) q[2];
sx q[2];
rz(2.1622369) q[2];
rz(-2.7566946) q[3];
sx q[3];
rz(-1.220547) q[3];
sx q[3];
rz(2.9773007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591831) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(2.3517877) q[0];
rz(-0.18665953) q[1];
sx q[1];
rz(-1.4013545) q[1];
sx q[1];
rz(2.1479215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2659205) q[0];
sx q[0];
rz(-1.8199931) q[0];
sx q[0];
rz(0.56206352) q[0];
rz(-pi) q[1];
rz(2.5285401) q[2];
sx q[2];
rz(-1.7568577) q[2];
sx q[2];
rz(1.2836054) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6911654) q[1];
sx q[1];
rz(-1.9948927) q[1];
sx q[1];
rz(2.0753839) q[1];
x q[2];
rz(1.8422115) q[3];
sx q[3];
rz(-0.86835734) q[3];
sx q[3];
rz(-2.0347119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8197202) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(2.7703088) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-2.0472417) q[3];
sx q[3];
rz(2.546052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6868941) q[0];
sx q[0];
rz(-1.0201539) q[0];
sx q[0];
rz(1.4042847) q[0];
rz(-2.4644201) q[1];
sx q[1];
rz(-1.9858457) q[1];
sx q[1];
rz(1.6489395) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9335404) q[0];
sx q[0];
rz(-1.3553936) q[0];
sx q[0];
rz(-0.24788863) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1588232) q[2];
sx q[2];
rz(-2.6361536) q[2];
sx q[2];
rz(2.7014521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6552393) q[1];
sx q[1];
rz(-2.7905373) q[1];
sx q[1];
rz(0.69610657) q[1];
rz(-pi) q[2];
rz(-1.9203606) q[3];
sx q[3];
rz(-0.75569433) q[3];
sx q[3];
rz(0.44103482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45381418) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(2.7933534) q[2];
rz(-1.6866775) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(-1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9652047) q[0];
sx q[0];
rz(-0.56469733) q[0];
sx q[0];
rz(0.89214605) q[0];
rz(-2.6773894) q[1];
sx q[1];
rz(-1.8966388) q[1];
sx q[1];
rz(-1.8468599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034445914) q[0];
sx q[0];
rz(-2.4305516) q[0];
sx q[0];
rz(-0.76816316) q[0];
rz(0.53769298) q[2];
sx q[2];
rz(-1.6811996) q[2];
sx q[2];
rz(0.7796208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66440873) q[1];
sx q[1];
rz(-2.5818965) q[1];
sx q[1];
rz(-2.2137292) q[1];
x q[2];
rz(0.89035676) q[3];
sx q[3];
rz(-2.0687752) q[3];
sx q[3];
rz(-0.80163664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(2.5757705) q[2];
rz(3.0602509) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(-2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.6767947) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(-1.5860522) q[0];
rz(1.0606891) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(-2.5097844) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.466296) q[0];
sx q[0];
rz(-1.6501556) q[0];
sx q[0];
rz(2.9064889) q[0];
rz(-pi) q[1];
x q[1];
rz(1.00497) q[2];
sx q[2];
rz(-2.3477051) q[2];
sx q[2];
rz(-0.53295202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7673847) q[1];
sx q[1];
rz(-2.4904618) q[1];
sx q[1];
rz(-0.99094772) q[1];
rz(-pi) q[2];
rz(1.7177714) q[3];
sx q[3];
rz(-2.2490361) q[3];
sx q[3];
rz(-0.036546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1598728) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(2.6427606) q[2];
rz(1.3129129) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(1.7074728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(2.4420807) q[0];
rz(-2.6761159) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(2.2043998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584211) q[0];
sx q[0];
rz(-1.4025549) q[0];
sx q[0];
rz(2.9267163) q[0];
rz(-pi) q[1];
x q[1];
rz(1.466339) q[2];
sx q[2];
rz(-1.3323931) q[2];
sx q[2];
rz(1.6676211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.921007) q[1];
sx q[1];
rz(-1.7821508) q[1];
sx q[1];
rz(1.2799954) q[1];
rz(-pi) q[2];
rz(1.2770416) q[3];
sx q[3];
rz(-1.2806727) q[3];
sx q[3];
rz(-1.4690831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.297544) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(-0.37332264) q[2];
rz(1.1897872) q[3];
sx q[3];
rz(-0.50896421) q[3];
sx q[3];
rz(-2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.3968286) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(2.3936791) q[0];
rz(-2.3751936) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(3.1386197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52067962) q[0];
sx q[0];
rz(-1.3786043) q[0];
sx q[0];
rz(2.2635133) q[0];
rz(-pi) q[1];
rz(-1.8393821) q[2];
sx q[2];
rz(-3.0680484) q[2];
sx q[2];
rz(3.1036249) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6312478) q[1];
sx q[1];
rz(-1.5273792) q[1];
sx q[1];
rz(-1.7249291) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19488867) q[3];
sx q[3];
rz(-1.3759383) q[3];
sx q[3];
rz(3.0611567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11091867) q[2];
sx q[2];
rz(-1.3838394) q[2];
sx q[2];
rz(-2.0745011) q[2];
rz(-3.057632) q[3];
sx q[3];
rz(-2.6559918) q[3];
sx q[3];
rz(0.77897227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(-3.0352266) q[0];
rz(-2.1616409) q[1];
sx q[1];
rz(-1.4915024) q[1];
sx q[1];
rz(-0.76593691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.879012) q[0];
sx q[0];
rz(-1.1404788) q[0];
sx q[0];
rz(0.58027123) q[0];
x q[1];
rz(-2.9418403) q[2];
sx q[2];
rz(-1.6589266) q[2];
sx q[2];
rz(2.5823808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3043266) q[1];
sx q[1];
rz(-2.2464753) q[1];
sx q[1];
rz(1.5262239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25772734) q[3];
sx q[3];
rz(-3.1270087) q[3];
sx q[3];
rz(2.7590883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6672259) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-2.4412947) q[2];
rz(2.4750366) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(-1.0880281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67548442) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(1.3960557) q[0];
rz(1.667977) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(-0.6667164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4182435) q[0];
sx q[0];
rz(-0.94941345) q[0];
sx q[0];
rz(0.80051144) q[0];
rz(1.2856977) q[2];
sx q[2];
rz(-0.93109967) q[2];
sx q[2];
rz(-0.49298795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.085070193) q[1];
sx q[1];
rz(-1.4398972) q[1];
sx q[1];
rz(0.31899778) q[1];
rz(-pi) q[2];
rz(1.0239086) q[3];
sx q[3];
rz(-1.0677665) q[3];
sx q[3];
rz(1.5707113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0746158) q[2];
sx q[2];
rz(-2.4846027) q[2];
sx q[2];
rz(-0.89078772) q[2];
rz(0.43205076) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(-2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73410949) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(-1.8941849) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(1.9966077) q[2];
sx q[2];
rz(-0.37526423) q[2];
sx q[2];
rz(3.0107977) q[2];
rz(-1.6040989) q[3];
sx q[3];
rz(-1.0819482) q[3];
sx q[3];
rz(-0.29500189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
